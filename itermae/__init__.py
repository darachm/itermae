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

    def apply_filters(self, filters):
        """
        This is for applying written filters to the results, so you can fail
        things that don't look right by position of the matches, or the
        statistics of each match. 
        First we unpack all the group and match stats/scores, so you can
        access them in defining filters easily.
        Then we're just straight eval'ing in that context, because I'm not
        thinking about security at all.
        """

        env_thing = { **self.group_stats , **self.match_scores }
        for i in self.seqs:
            if i in env_thing.keys():
                env_thing[i].seq = self.seqs[i].seq

        return_object = []
        try:
            for i in range(len(filters)):
                return_object.append(False)
                # Here we evaluate them but using that dictionary as the
                # global dictionary, because done is better than dogma.
                try:
                    if eval(filters[i],globals(),env_thing):
                        return_object[i] = True
                except:
                    return_object[i] = False
        except:
            return_object.append(False)

        return return_object

    def build_output(self,output_id_def,output_seq_def):
        """
        Similar thing as above, but just making it flat of all the seqs
        so you can build what you want in the outputs. First we make the output
        seq object, then the ID (which can have seqs in it, as part of the ID, 
        so like extracted UMIs or sample-indicies).
        """

        env_thing = { **self.seqs }

        out_seq = SeqRecord.SeqRecord(Seq.Seq(""))
        out_seq = eval(output_seq_def,globals(),env_thing)
        out_seq.id = str(eval(output_id_def,globals(),env_thing))

        return out_seq

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


def reader(
        input_file, is_gzipped, 
        operations_array, filters , outputs_array,
        in_format, out_format, output_file, failed_file, report_file,
        verbosity
        ):
    """
    This reads inputs, calls the `chop` function on each one, and sorts it
    off to outputs. So this is called by the main function, and is all about
    the IO. 
    """

    #
    # Open up file handles
    #

    # Input
    if input_file == "STDIN": # STDIN is default
        if is_gzipped: # sys.stdin is opened as txt, not bytes
            print("I can't handle gzipped inputs on STDIN ! Un-gzip for me. "+
                "Or write to file, and point me that-a-way. \n\n"+
                "You shouldn't see this error, it shoulda been caught in the "+
                "launcher script.",file=sys.stderr) 
            exit(1)
        input_text = sys.stdin
    else:
        if is_gzipped:
            input_text = gzip.open(input_file,'rt',encoding='ascii')
        else:
            input_text = open(input_file,"rt") 

    input_seqs = open_appropriate_input_format(input_text, in_format)

#    # Input
#    if input_file == "STDIN": # STDIN is default
#        if is_gzipped: # sys.stdin is opened as txt, not bytes
#            print("I can't handle gzipped inputs on STDIN ! Un-gzip for me. "+
#                "Or write to file, and point me that-a-way. \n\n"+
#                "You shouldn't see this error, it shoulda been caught in the "+
#                "launcher script.",file=sys.stderr) 
#            exit(1)
#        input_fh = sys.stdin
#        input_seqs = open_appropriate_input_format(input_fh, in_format)
#    else:
#        input_fh = open(input_file,"rb") 
#        if is_gzipped:
#            input_fh_gz = gzip.open(input_fh,'rt',encoding='ascii')
#            input_seqs = open_appropriate_input_format(input_fh_gz, in_format)
#        else:
#            input_seqs = open_appropriate_input_format(
#                    io.TextIOWrapper(input_fh),
#                    in_format)

    # Outputs - passed records, failed records, report file
    if output_file == "STDOUT":
        output_fh = sys.stdout
    else: # If you've specified a filepath, then that's here
        output_fh = open(output_file,"a")
    # If no failed file specified, then we're just ignoring it
    if failed_file is None:
        failed_fh = None
    # But if specified, then it gets written
    else:
        failed_fh = open(failed_file,"a")
    # Same for optional report
    if report_file == None:
        report_fh = None
    else:
        report_fh = open(report_file,"a")

    # Do the chop-ing...
    for each_seq in input_seqs:
        # Each sequence, one by one...

        chop(
            seq_holder=SeqHolder(each_seq,verbosity=verbosity),  
            operations_array=operations_array, filters=filters, 
            outputs_array=outputs_array, 
            in_format=in_format,
            out_format=out_format, 
            output_fh=output_fh, failed_fh=failed_fh, report_fh=report_fh,
            verbosity=verbosity
            )

    input_text.close()

    return(0)


def chop(
    seq_holder,
    operations_array, filters, outputs_array, 
    in_format, out_format,
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
    assert operations_array[0][0] == "input", (
        "can't find the sequence named `input`, rather we see `"+
        operations_array[0][0]+"` in the holder, so breaking. You should "+
        "have the first operation start with `input` as a source." )

    #
    #
    # ITERATING THROUGH THE MATCHING
    #
    #

    # First, apply each operation !
    for operation_number, operation in enumerate(operations_array):

        if operation_number > 26:
            print("Listen, here's the deal. I did not anticipate anyone would "+
                "be trying more than a few operations, so the IDs are just "+
                "one letter. So, use <=26 operations, or rewrite it "+
                "yourself around where it calls `enumerate(operations_array)`.",
                file=sys.stderr)
            exit(1)

        seq_holder.apply_operation( string.ascii_lowercase[operation_number],
                operation[0],operation[1] )

    # Now seq_holder should have a lot of goodies, match scores and group stats
    # and matched sequences groups.
    # All these values allow us to apply filters :

    #
    #
    # APPLYING FILTERS
    #
    #

    evaluated_filters = seq_holder.apply_filters(filters) 

    # This evaluated_filters should be boolean list. So did we pass all filters?
    # If not then do this
    if not all(evaluated_filters):

        if verbosity >= 2:
            print("\n["+str(time.time())+"] : match is : evaluated the "+
                "filters as : "+str(evaluated_filters)+" and so failed.", 
                file=sys.stderr)

        # So if we should write this per-record report
        if report_fh != None:
            print( seq_holder.format_report("FailedFilter",
                    seq_holder.seqs['input'].seq, evaluated_filters)
                ,file=report_fh)

        if failed_fh != None:
            SeqIO.write(seq_holder.seqs['input'], failed_fh, "fastq")

        return 0

    # Or if we passed all filters, then we try to form the outputs
    else:

        #
        #
        # FORMING OUTPUTS
        #
        #
        try:

            # We attempt to form the correct output records
            output_records = [ seq_holder.build_output(i, j) for i, j in outputs_array ]
            # So this will fail us out of the 'try' if it doesn't form.
            # i is the output_arrays ID spec, and j is sequence spec.

            # Format the output record as appropriate
            for which, output_record in enumerate(output_records):
                if out_format == "sam":
                    print(
                        format_sam_record(
                            output_record.id, str(output_record.seq),
                            ''.join([chr(i+33) for i in 
                                    output_record.letter_annotations['phred_quality']]
                                    ),
                            "IE:Z:"+str(which)
                        ),
                        file=output_fh
                    )
                elif out_format == "txt":
                    print(
                        str(output_record.seq),
                        file=output_fh
                    )
                else:
                    try:
                        SeqIO.write(output_record, output_fh, out_format) 
                    except:
                        print("I don't know '"+out_format+"' format, "+
                            "I know sam, txt, fastq, and fasta.",
                            file=sys.stderr) 
                        exit(1)

                # If we want to write the report, we make it
                if report_fh != None:
                    print( seq_holder.format_report("Passed",
                            output_record.seq, evaluated_filters)
                        ,file=report_fh)

            if verbosity >= 2:
                print("\n["+str(time.time())+"] : evaluated the "+
                    "filters as : "+str(evaluated_filters)+" and so passed.", 
                    file=sys.stderr)

            return 0

        #
        #
        # FAILED FORMATTING FAILS
        #
        #
        except:

            if verbosity >= 2:
                print("\n["+str(time.time())+"] : failed upon forming the "+
                    "output.", file=sys.stderr)

            # If we want to write the report, we make it
            if report_fh != None:
                print( 
                    seq_holder.format_report("FailedDirectivesToMakeOutputSeq",
                        seq_holder.seqs['input'].seq, evaluated_filters)
                    ,file=report_fh)

            if failed_fh != None:
                SeqIO.write(seq_holder.seqs['input'], failed_fh, "fastq")

            return 0



