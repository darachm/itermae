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

        env_thing_stats = { **self.group_stats , **self.match_scores }
        env_thing_seqs = { **self.seqs }

        return_object = []
        try:
            for i in range(len(filters)):
                # Here we evaluate them but using that dictionary as the
                # global dictionary, because done is better than dogma.
                try:
                    if eval(filters[i],globals(),env_thing_stats):
                        return_object[i] = True
                    else:
                        return([False])
                except:
                    pass
                try:
                    if eval(filters[i],globals(),env_thing_seqs):
                        return_object[i] = True
                    else:
                        return([False])
                except:
                    pass
        except:
            return([False])

        return return_object

    def build_output(self,output_id_def,output_seq_def):
        """
        Similar thing as above, but just making it flat of all the seqs
        so you can build what you want in the outputs. First we make the output
        seq object, then the ID (which can have seqs in it, as part of the ID, 
        so like extracted UMIs or sample-indicies).
        """

        env_thing = { **self.seqs }
        out_seq = eval(output_seq_def,globals(),env_thing)
        out_seq.id = str(eval(output_id_def,globals(),env_thing))

        return copy.deepcopy(out_seq)

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


def reader(
        input_file, is_gzipped, 
        operations_array, filters , outputs_array,
        out_format, output_file, failed_file, report_file,
        verbosity
        ):
    """
    This reads inputs, calls the `chop` function on each one, and sorts it
    off to outputs. So this is called by the main function, and is all about
    the IO. 
    """

    ### Open up file handles

    # If that's STDIN, which is default, then we're taking sequences by STDIN
    if input_file is "STDIN":
        if is_gzipped:
            with gzip.open(sys.stdin,"rt") as input_file_gz:
                input_seqs = SeqIO.parse(input_file_gz,"fastq")
                # No idea if this works, is at all sensible
        else:
            input_seqs = SeqIO.parse(sys.stdin,"fastq")
    else:
        # Or if it's gzipped then it's from a gzipped file (but no gzipped
        # STDIN plz, just zcat it
        if is_gzipped:
            with gzip.open(input_file,"rt") as input_file_gz:
                input_seqs = SeqIO.parse(input_file_gz,"fastq")
        # Otherwise is a flat file I assume
        else:
            input_seqs = SeqIO.parse(input_file,"fastq")

    # Opening up output file handles, will hand them off to each chop 
    if output_file is "STDOUT":
        output_fh = sys.stdout
    # If you've specified a file, then that's here
    else:
        output_fh = open(output_file,"a")
    # If no failed file specified, then we're just ignoring it
    if failed_file is None:
        failed_fh = None
    # But if specified, then it gets written
    else:
        failed_fh = open(failed_file,"a")
    # Same for optional report
    if report_file is None:
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
            out_format=out_format, 
            output_fh=output_fh, failed_fh=failed_fh, report_fh=report_fh,
            verbosity=verbosity
            )

    return(0)


def chop(
    seq_holder,
    operations_array, filters, outputs_array, 
    out_format,
    output_fh, failed_fh, report_fh,
    verbosity
    ):
    """
    This one takes each record, applies the operations, evaluates the filters,
    generates outputs, and writes them to output handles as appropriate.
    It's a bit messy, so I've tried to make it clear with comments to break it
    up into sections.
    """

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
        if report_fh is not None:
            print( seq_holder.format_report("FailedFilter",
                    seq_holder.seqs['input'].seq, evaluated_filters)
                ,file=report_fh)

        if failed_fh is not None:
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
                if out_format == "SAM":
                    print(
                        "\t".join([
                            output_record.id,
                            "0", "*", "0", "255", "*", "=", "0", "0", 
                            str(output_record.seq),
                            ''.join([chr(i+33) for i in 
                                    output_record.letter_annotations['phred_quality']]
                                    ),
                            "XI:"+str(which)
                            ])
                        ,file=output_fh)
                elif out_format == "FASTQ":
                    SeqIO.write(output_record, output_fh, "fastq") 
                elif out_format == "FASTA":
                    SeqIO.write(output_record, output_fh, "fasta") 
                else:
                    print("I don't know '"+out_format+"' format, "+
                        "exiting over that. I know SAM, FASTQ, and FASTA.") 
                    exit(1)

                # If we want to write the report, we make it
                if report_fh is not None:
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
            if report_fh is not None:
                print( 
                    seq_holder.format_report("FailedDirectivesToMakeOutputSeq",
                        seq_holder.seqs['input'].seq, evaluated_filters)
                    ,file=report_fh)

            if failed_fh is not None:
                SeqIO.write(seq_holder.seqs['input'], failed_fh, "fastq")

            return 0



