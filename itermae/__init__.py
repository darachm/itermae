#!/usr/bin/env python3

# Importing packages for programming, odds and ends
import time
import statistics
import sys
import gzip
import string
import argparse
import re
import itertools
import copy
import io

# Importing packages for the heart of it, fuzzy regex and SeqIO classes
import regex
from Bio import SeqIO
from Bio import Seq, SeqRecord

class MatchScores:
    """
    This just makes an object to hold these three where they're easy to type
    (as attributes not keyed dict). Well, and a flatten function for printing.
    """
    def __init__(self, substitutions, insertions, deletions):
        self.substitutions = substitutions
        self.insertions = insertions
        self.deletions = deletions
    def flatten(self):
        return str(self.substitutions)+"_"+\
            str(self.insertions)+"_"+\
            str(self.deletions)


class GroupStats:
    """
    This just makes an object to hold these three where they're easy to type
    (as attributes not keyed dict). Well, and a flatten function for printing.
    """
    def __init__(self, start, end, quality):
        self.start = start 
        self.end = end 
        self.length = self.end - self.start
        self.quality = quality
    def flatten(self):
        return str(self.start)+"_"+str(self.end)+"_"+str(self.length)


class SeqHolder: 
    """
    This is the main holder of sequences, and does the matching and stuff.
    I figured a Class might make it a bit tidier.
    """
    def __init__(self, input_record, verbosity=4):
        # So the .seqs holds the sequences accessed by the matching, and there's
        # a dummyspacer in there just for making outputs where you want that
        # for later partitioning. Input is input.
        self.seqs = {
            'dummyspacer': SeqRecord.SeqRecord(Seq.Seq("X"),id="dummyspacer"),
            'input': input_record }
        self.seqs['dummyspacer'].letter_annotations['phred_quality'] = [40]
        self.verbosity = verbosity
        # These two dicts hold the scores for each match operation (in order),
        # and the start end length statistics for each matched group.
        self.match_scores = {}
        self.group_stats = {}

    def apply_operation(self, match_id, input_group, regex):
        """
        This applies the matches, saves how it did, and saves extracted groups.
        Details commented below.
        """

        # Try to find the input, if it ain't here then just return
        try: 
            self.seqs[input_group]
        except:
            self.match_scores[match_id] = MatchScores(None,None,None)
            return self

        if self.verbosity >= 3:
            print("\n["+str(time.time())+"] : attempting to match : "+
                str(regex)+" against "+self.seqs[input_group].seq,
                file=sys.stderr)

        # Here we execute the actual meat of the business.
        # Note that the input is made uppercase.
        fuzzy_match = regex.search( str(self.seqs[input_group].seq).upper() )

        if self.verbosity >= 3:
            print("\n["+str(time.time())+"] : match is : "+str(fuzzy_match),
                file=sys.stderr)

        try:
            # This is making and storing an object for just accessing these
            # numbers nicely in the arguments for forming outputs and filtering.
            self.match_scores[match_id] = MatchScores(*fuzzy_match.fuzzy_counts)

            # Then for each of the groups matched by the regex
            for match_name in fuzzy_match.groupdict():
    
                # We stick into the holder a slice of the input seq, that is 
                # the matched # span of this matching group. So, extract.
                self.seqs[match_name] = \
                    self.seqs[input_group][slice(*fuzzy_match.span(match_name))]

                self.seqs[match_name].description = "" 
                # This is to fix a bug where the ID is stuck into the 
                # description and gets unpacked on forming outputs

                # Then we record the start, end, and length of the matched span
                self.group_stats[match_name] = \
                    GroupStats(*fuzzy_match.span(match_name),
                        quality=self.seqs[match_name].letter_annotations['phred_quality']
                        )

        except:
            self.match_scores[match_id] = MatchScores(None,None,None)

    def build_context(self):
        """
        This just unpacks group match stats/scores into an environment that
        the filter can then use to ... well ... filter. 
        """

        # This is context for the filters, so is operating more as values,
        # as opposed to the context_seq which is operating with SeqRecords
        self.context_filter = { **self.group_stats , **self.match_scores }
        for i in self.seqs:
            if i in self.context_filter.keys():
                self.context_filter[i].seq = self.seqs[i].seq
                    # Also adding on the actual sequence, so it's accessible

        # Then unpack the sequences as a context for building the output 
        # sequences, this is different so that the qualities get stuck with
        # the bases of the groups
        self.context_seq = { **self.seqs }

    def evaluate_filter_of_output(self,output_dict):
        """
        This tests a defined filter on the 'seq_holder' object
        """

        try:
            return eval(output_dict['filter'],globals(),self.context_filter)
        except:
            if self.verbosity >= 2:
                print("\n["+str(time.time())+"] : This read "+
                    self.seqs['input'].id+" failed to evaluate the filter "+
                    str(output_dict['filter']),file=sys.stderr)
            return False

    def build_output(self,output_dict):
        """
        This builds the output
        """

        try:
            out_seq = SeqRecord.SeqRecord(Seq.Seq(""))
            out_seq = eval(output_dict['seq'],globals(),self.context_seq)
            out_seq.id = str(eval(output_dict['id'],globals(),self.context_seq))
            return out_seq
        except:
            if self.verbosity >= 2:
                print("\n["+str(time.time())+"] : This read "+
                    self.seqs['input'].id+" failed to build the output "+
                    "id: "+str(output_dict['id'])+
                    "seq: "+str(output_dict['seq']) ,file=sys.stderr)
            return None

    def format_report(self,label,output_seq,evaluated_filters):
        """
        This is for formatting a standard report line for the reporting function
        """
        
        return ( "\""+label+"\",\""+
            str(self.seqs['input'].id)+"\",\""+
            str(self.seqs['input'].seq)+"\",\""+
            str(output_seq)+"\",\""+
            str(evaluated_filters)+"\",\""+
            "-".join([ i+"_"+self.group_stats[i].flatten() 
                        for i in self.group_stats ] )+
            "\"" )


def format_sam_record(record_id, sequence, qualities, tags,
        flag='0', reference_name='*', 
        mapping_position='0', mapping_quality='255', cigar_string='*',
        reference_name_of_mate='=', position_of_mate='0', template_length='0'
    ):
    return "\t".join([
            record_id,
            flag,
            reference_name,
            mapping_position,
            mapping_quality,
            cigar_string,
            reference_name_of_mate,
            position_of_mate,
            template_length,
            sequence,
            qualities,
            tags
        ])


def read_sam_file(fh):
    """
    This is a minimal reader, just for getting the fields I like and emiting
    SeqRecord objects, sort of like SeqIO. Putting SAM tags in description.
    """
    for i in fh.readlines():
        fields = i.rstrip('\n').split('\t')
        yield SeqRecord.SeqRecord(
            Seq.Seq(fields[9]),
            id=fields[0],
            letter_annotations={'phred_quality':[ord(i)-33 for i in fields[10]]},
            description=fields[11]
            )


def read_txt_file(fh):
    """
    This just treats one sequence per line as a SeqRecord.
    """
    for i in fh.readlines():
        seq = i.rstrip()
        yield SeqRecord.SeqRecord( Seq.Seq(seq), id=seq )


def open_appropriate_input_format(in_fh, format_name):
    if   format_name == 'fastq':
        return SeqIO.parse(in_fh, format_name)
    elif format_name == 'sam':
        return iter(read_sam_file(in_fh))
    elif format_name == 'fasta':
        return SeqIO.parse(in_fh, format_name)
    elif format_name == 'txt':
        return iter(read_txt_file(in_fh))
    else:
        print("I don't know that input file format name '"+format_name+
            "'. I will try and use the provided format name in BioPython "+
            "SeqIO, and we will find out together if that works.",
            file=sys.stderr) 
        return SeqIO.parse(in_fh, format_name)


def open_input_fh(file_string,gzipped=False):
    if file_string.upper() == 'STDIN':
        if gzipped:
            print("I can't handle gzipped inputs on STDIN ! "+
                "You shouldn't see this error, it shoulda been caught in "+
                "the launcher script.",file=sys.stderr) 
            exit(1)
        else:
            return open(sys.stdin,'rt')
    else:
        if gzipped:
            return gzip.open(file_string,'rt',encoding='ascii')
        else:
            return open(file_string,'rt')


def open_output_fh(file_string):
    if file_string.upper() == 'STDOUT':
        return sys.stdout
    elif file_string.upper() == 'STDERR':
        return sys.stderr
    else:
        return open(file_string,'a')


def reader(configuration):
    """
    This reads inputs, calls the `chop` function on each one, and sorts it
    off to outputs. So this is called by the main function, and is mostly about
    handling the I/O. 
    """

    ### Open up file handles

    # Input
    input_seqs = open_appropriate_input_format(
        open_input_fh(configuration['input'],configuration['input_gzipped']),
        configuration['input_format'])

    # Outputs - passed records, failed records, report file
    output_fh = open_output_fh(configuration['output'])
    try:
        failed_fh = open_output_fh(configuration['failed']),
    except:
        failed_fh = None
    try:
        report_fh = open_output_fh(configuration['report'])
    except:
        report_fh = None

    # Do the chop-ing...
    for each_seq in input_seqs:
            # Each sequence, one by one...
        chop(
            seq_holder=SeqHolder(each_seq,verbosity=configuration['verbosity']),  
            operations_array=configuration['matches'],
            outputs_array=configuration['output_groups'],
            out_format=configuration['output_format'],
            input_format=configuration['input_format'],
            output_fh=output_fh, failed_fh=failed_fh, report_fh=report_fh,
            verbosity=configuration['verbosity']
            )

    for i in [ input_seqs, output_fh, failed_fh, report_fh] :
        try:
            i.close()
        except:
            pass

    return(0)


def write_out_seq(seq,fh,format,which):
    if format == "sam":
        print( format_sam_record( seq.id, str(seq.seq),
                ''.join([chr(i+33) for i in 
                        seq.letter_annotations['phred_quality']]),
                "IE:Z:"+str(which) ),file=fh)
    elif format == "txt":
        print( str(seq.seq), file=fh)
    else:
        try:
            SeqIO.write(seq, fh, format) 
        except:
            print("I don't know '"+format+"' format, and it doesn't seem to "+
                "work with BioPython.", file=sys.stderr) 
            exit(1)


def chop(
    seq_holder,
    operations_array, outputs_array, 
    out_format, input_format,
    output_fh, failed_fh, report_fh,
    verbosity
    ):
    """
    This one takes each record, applies the operations, evaluates the filters,
    generates outputs, and writes them to output handles as appropriate.
    It's a bit messy, so I've tried to make it clear with comments to break it
    up into sections.
    """

    # If qualities are missing, add them as just 40
    if 'phred_quality' not in seq_holder.seqs['input'].letter_annotations.keys():
        seq_holder.seqs['input'].letter_annotations['phred_quality'] = [40]*len(seq_holder.seqs['input'])
        if verbosity >= 2:
            print("\n["+str(time.time())+"] : adding missing qualities of 40 "+
                "to sequence.",
                file=sys.stderr)

    # For chop grained verbosity, report
    if verbosity >= 2:
        print("\n["+str(time.time())+"] : starting to process : "+
            seq_holder.seqs['input'].id+"\n  "+seq_holder.seqs['input'].seq+"\n  "+ 
            str(seq_holder.seqs['input'].letter_annotations['phred_quality']),
            file=sys.stderr)

    # This should fail if you didn't specify anything taking from input stream!
    assert operations_array[0]['input'] == "input", (
        "can't find the sequence named `input`, rather we see `"+
        operations_array[0]['input']+"` in the holder, so breaking. You should "+
        "have the first operation start with `input` as a source." )

    # ITERATING THROUGH THE MATCHING

    # First, apply each operation !
    for operation_number, operation in enumerate(operations_array):

        seq_holder.apply_operation( 'match_'+str(operation_number),
                operation['input'], operation['regex'] )

    # Now seq_holder should have a lot of goodies, match scores and group stats
    # and matched sequences groups.
    # All these values allow us to apply filters :

    ### Filtering and generating outputs

    # First unpacking matches and scores into an internal environment for
    # the filter 'eval's
    seq_holder.build_context()

    # Then we eval the filters and build outputs, for each output
    output_records = []
    for each_output in outputs_array:
        output_records.append(
            ( seq_holder.evaluate_filter_of_output(each_output), 
                seq_holder.build_output(each_output) )
        )

    passed_filters = not any([ i == False for i,j in output_records ])

    # So if we should write this per-record report
    if report_fh != None:
        if passed_filters:
            print( seq_holder.format_report("PassedFilters",
                    seq_holder.seqs['input'].seq,
                    str([i is not None for i in output_records]) )
                ,file=report_fh)
        else:
            print( seq_holder.format_report("FailedAtLeastOutput",
                    seq_holder.seqs['input'].seq,
                    str([i is not None for i in output_records]) )
                ,file=report_fh)

    if failed_fh != None:
        if not passed_filters:
            write_out_seq(seq_holder.seqs['input'], failed_fh, input_format,0)

    for which, output_record in enumerate(output_records):
        if output_record is not None:
            write_out_seq(output_record, output_fh, out_format, which)



