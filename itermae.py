#!/usr/bin/env python3

# Importing packages for programming, odds and ends
import argparse
import re
import itertools
import json
import time
import statistics
import sys
import gzip
import string

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
    def __init__(self, start, end):
        self.start = start 
        self.end = end 
        self.length = self.end - self.start
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
    
                # Then we record the start, end, and length of the matched span
                self.group_stats[match_name] = GroupStats(*fuzzy_match.span(match_name))

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

        return_object = []
        try:
            for each_filter in filters:
                # Here we evaluate them but using that dictionary as the
                # global dictionary, because done is better than dogma.
                if eval(each_filter,globals(),env_thing):
                    return_object.append(True)
                else:
                    return([False])
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





def reader(input_file, is_gzipped, 
        output_file, failed_file, report_file,
        operations_array, filters , outputs_array,
        out_format,
        verbosity
        ):

    """
    This reads inputs, calls the `chop` function on each one, and sorts it
    off to outputs.
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
        output = sys.stdout
    # If you've specified a file, then that's here
    else:
        output = open(output_file,"a")
    # If no failed file specified, then we're just ignoring it
    if failed_file is None:
        failed = None
    # But if specified, then it gets written
    else:
        failed = open(failed_file,"a")
    # Same for optional report
    if report_file is None:
        report = None
    else:
        report = open(report_file,"a")

    # Making a dummyspacer thing and seq_holder, this gets copied for each
    # chop function call, so it's where the internals can hold matches for
    # subsequent group searches. Dummy spacer is for odd formatting outputs.
    # Not very necessary now that I put it out as SAM...

    ### Do the chop-ing

    for each_seq in input_seqs:
        # Each sequence, one by one...

        chop(
            SeqHolder(each_seq,verbosity=verbosity),  
            operations_array, filters, outputs_array, # Things todo!
            output, failed, report, # File handles
            out_format, # Er... what format for output?
            verbosity
            )

    return(0)


def chop(
    seq_holder,
    operations_array, filters, outputs_array, 
    output, failed, report,
    out_format,
    verbosity
    ):
    """
    This one takes each record, applies the operations, evaluates the filters,
    generates outputs, and writes them to output handles as appropriate.
    """

    # Chop grained verbosity
    if verbosity >= 2:
        print("\n["+str(time.time())+"] : starting to process : "+
            seq_holder.seqs['input'].id+"\n  "+seq_holder.seqs['input'].seq+"\n  "+ 
            str(seq_holder.seqs['input'].letter_annotations['phred_quality']),
            file=sys.stderr)

    # This should fail if you didn't specify anything taking 
    # from input stream!
    assert operations_array[0][0] == "input", (
        "can't find the sequence named `input`, rather we see `"+
        operations_array[0][0]+"` in the holder, so breaking. You should "+
        "have the first operation start with `input` as a source." )

    for operation_number, operation in enumerate(operations_array):

        if operation_number > 26:
            print("Listen, here's the deal. I did not anticipate anyone would "+
                "be trying more than a few operations, so the IDs are just "+
                "one letter. So, use fewer operations, or rewrite it "+
                "yourself around when it calls `enumerate(operations_array)`.",
                file=sys.stderr)

        seq_holder.apply_operation( string.ascii_lowercase[operation_number],
                operation[0],operation[1] )

    # All these values allow use to apply filters, using this
    # function

    evaluated_filters = seq_holder.apply_filters(filters) 

    # This evaluated_filters should be logical list. So did we pass all filters?
    if not all(evaluated_filters):

        if verbosity >= 2:
            print("\n["+str(time.time())+"] : match is : evaluated the "+
                "filters as : "+str(evaluated_filters)+" and so failed.", 
                file=sys.stderr)

        # So if we should write this per-record report
        if report is not None:
            print( seq_holder.format_report("FailedFilter",
                    seq_holder.seqs['input'].seq, evaluated_filters)
                ,file=report)

        if failed is not None:
            SeqIO.write(seq_holder.seqs['input'], failed, "fastq")

        return 0

    else:


        try:
            # We attempt to form the correct output records

            output_records = [ seq_holder.build_output(i, j) for i, j in outputs_array ]
            # So this will fail us out of the 'try' if it doesn't form

            # Otherwise, just the record and if it passed
            for which, output_record in enumerate(output_records):
                if out_format == "sam":
                    print(
                        "\t".join([
                            output_record.id,
                            "0", "*", "0", "255", "*", "=", "0", "0", 
                            str(output_record.seq),
                            ''.join([chr(i+33) for i in output_record.letter_annotations['phred_quality']]),
                            "XI:"+str(which)
                            ])
                        ,file=output)
                elif out_format == "fastq":
                    SeqIO.write(output_record, output, "fastq") 
                else:
                    print("I don't "+out_format+" format, looping over here") 

                # If we want to write the report, we make it
                if report is not None:
                    print( seq_holder.format_report("Passed",
                            output_record.seq, evaluated_filters)
                        ,file=report)

            if verbosity >= 2:
                print("\n["+str(time.time())+"] : evaluated the "+
                    "filters as : "+str(evaluated_filters)+" and so passed.", 
                    file=sys.stderr)

            return 0

        except:

            if verbosity >= 2:
                print("\n["+str(time.time())+"] : failed upon forming the "+
                    "output.", file=sys.stderr)

            # If we want to write the report, we make it
            if report is not None:
                print( 
                    seq_holder.format_report("FailedDirectivesToMakeOutputSeq",
                        seq_holder.seqs['input'].seq, evaluated_filters)
                    ,file=report)

            if failed is not None:
                SeqIO.write(input_record, failed, "fastq")

            return 0



if __name__ == '__main__':
    # This `if` statement only runs if it's called alone, I believe. So I think 
    # you can import these functions for use in other scripts, but I'm not 
    # intending for that. This is supposed to be standalone.

    # Using argparse module to define the arguments

    # Name and description of this program
    parser = argparse.ArgumentParser(description=""+
        "itermae - Tool for iteratively chopping up each FASTQ read using "+
        "fuzzy regular expressions, for many many reads."+
        "\n\n"+
        "For best performance, we recommend you use this with GNU parallel, "+
        "see examples in the README.")

    # Input determination 
    parser.add_argument("--input",default="STDIN",
        help="We expect this flag is normally NOT USED, and so I will expect "+
            "FASTQ(Z) input on the STDIN. However, you can specify a path "+
            "using this flag to read from a file. This is not parallelized, "+
            "so it's slower. But it's there. "+
            "Importantly, ALL SEQUENCE is coerced to uppercase letters before "+
            "matching, for speed. So keep it uppercase in the regex.")
    # Is it gzip'd ?
    parser.add_argument("-z","--gzipped",action="store_true",
        help="This flag indicates that input, be it STDIN or a file path, "+
            "is GZIPPED FASTQ. Note that this implies that any STDIN usage is "+
            "not parallel, so we do not expect this to be used in combination "+
            "with the STDIN input. But it's there.")

    # Output determination
    parser.add_argument("--output",default="STDOUT",
        help="we expect this flag is normally NOT USED, and so output would "+
            "be send to STDOUT. However, it's here incase you'd like to write "+
            "it to a file. Note I do not add a suffix for you.")
    parser.add_argument("--output-format",default='sam',
        help="The output format specification. Default is an unmapped sam "+
            "('sam') so that the tabular nature can permit joining with other "+
            "parsed output from the reverse read, but you can also specify "+
            "'fastq' to get a FASTQ back.")
    parser.add_argument("--failed",default=None,
        help="Optional name of output file for failed reads, for debugging. "+
            "If you say 'STDOUT' then it'll go there, so you could use that "+
            "to collect failed and passed reads in the same 'sam' file for "+
            "example.")
    parser.add_argument("--report",default=None,
        help="Add this flag with a filename to print a report of per-read "+
            "statistics. That's a lot of disk writes btw, but "+
            "would be good in combination with a '--limit' argument "+
            "so that you can spec out the kinds of noise you got "+
            "in your data and debug the running.")

    # verbosity
    parser.add_argument("-v","--verbose",action="count",default=0,
        help="How much debugging issues should I pipe out to STDERR?"+
            "None is nothing, "+
            "-v is setup messages and start-stop messsages, "+
            "-v -v is worker-level details, "+
            "-v -v -v is chop-level details, "+
            "-v -v -v -v is each operation level details. "+
            "Keep in mind these are going to STDERR.")

    ### CLI operation and output specification
    # Operations
    parser.add_argument("--operation","-o",action="append",
        help="The pattern to match and extract. This has a specific and sort "+
            "of complicated syntax. Copy and tweak examples, and/or refer to "+
            "the documentation via the README.md file, or to the original "+
            "regex module documentation. Each operation is done in the order "+
            "you specify, so later operations can use previous matching "+
            "groups as inputs. But, you can't name a group 'input' or "+
            "'dummyspacer', I'm using those names ! "+
            "Importantly, WRITE ALL DNA BASES AS UPPERCASE. I won't coerce "+
            "it in the regex because that messes up group names.")
    # Filter specification
    parser.add_argument("--filter","-f",action="append",
        help="Filters for eliminating reads that don't pass some criteria. "+
            "You write these based on the groups captured in the operations "+
            "above. So a group named 'group' would have a start at "+
            "'group_start', same for the 'group_end' and 'group_length'. "+
            "You're welcome. Also, `statistics` package is loaded, so you "+
            "can use those expressions for means, medians, etc.")
    # Outputs
    parser.add_argument("--output-id",action="append",
        help="A list of output ID definitions, in the same order as for "+
            "output-seq (see that one, below probably). This is evaluated for "+
            "reads that pass filter. You can access 'input.id' to get the "+
            "original FASTQ read ID, and 'group.seq' to get the sequence of a "+
            "particular group, so for example you could do "+
            "'input.id+\"_\"+index.seq' to append the index sequence to the "+
            "FASTQ ID. Unlike output-seq, this a string not a Sequence object.")
    parser.add_argument("--output-seq",action="append",
        help="A list of output seq definitions, in the same order as for "+
            "output-ids. This is evaluated on the groups being Biopython "+
            "SeqRecords, so you you can paste them together like "+
            "'sample+barcode' or 'sample+dummyspacer+barcode+dummyspacer+umi' "+
            "if you'd like a spacer X character to delimit them specifically. "+
            "You can't access .id or _length properties or the like, just "+
            "groups. Unlike output-id, this a Sequence object not a string.")

    ### Parse the arguments, check them.

    args = parser.parse_args()

    # Operations, outputs are read as an array of dicts to keep it ordered.
    operations_array = []
    outputs_array = []
    # Here we read on through to add those on, and complain loudly if someone
    # tries to use our reserved names.
    try:
        for each in args.operation:
            if each.find("<dummyspacer>") > 0:
                print("Hey, you can't name a capture group "+
                    "'dummyspacer', I'm using that! Pick a different name."
                    ,file=sys.stderr)
                exit(1)

            if each.find("<input>") > 0:
                print("Hey, you can't name a capture group "+
                    "'input', I'm using that! Pick a different name."
                    ,file=sys.stderr)
                exit(1)
            # Here, we use the ` > ` to specify the flow of input to
            # the regex
            (input_string, regex_string) = re.split("\s>\s",each.strip())
            compiled_regex = regex.compile(
                regex_string.strip(), # We use this regex
                regex.BESTMATCH # And we use the BESTMATCH strategy, I think
                )
            # append that to the operations array
            operations_array.append( [input_string.strip(), compiled_regex] )
    except:
        # Failure likely from lack of operations to do
        print("Wait a second, I don't understand the operations to be done! "+
            "Are there any? Maybe there's small part I'm choking on? Maybe "+
            "try adding steps in one at a time in an interactive context with "+
            "'--limit' set, to debug easier. Exiting...",file=sys.stderr)
        exit(1)

    # Next we build the array of outputs, by combining the named output IDs
    # and seq specifications.
    try:
        for each_id, each_seq in zip(args.output_id, args.output_seq):
            # append that to the outputs array
            outputs_array.append( [each_id, each_seq] )
    except:
        # Failure likely from lack of operations to do
        print("Wait a second, I don't understand the outputs to be done! "+
            "Are there any? Maybe there's small part I'm choking on? Maybe "+
            "try adding steps in one at a time in an interactive context with "+
            "'--limit' set, to debug easier. Exiting...",file=sys.stderr)
        exit(1)

    if args.verbose >= 1:
        print("\n["+str(time.time())+"] : "+
            "I'm reading in something, applying these operations of "+
            "alignment:\n",file=sys.stderr)
        for each in operations_array:
            print("  - from : "+each[0]+"\n"+
                "    extract groups with regex : '"+str(each[1])
                ,file=sys.stderr)

    if args.verbose >= 1:
        print("\n["+str(time.time())+"] : ...and with these filters:\n",
            file=sys.stderr)
        try:
            for i in args.filter:
                print("  - "+i,file=sys.stderr)
        except:
            print("  ( no filters defined )",file=sys.stderr)

    if args.verbose >= 1:
        print("\n["+str(time.time())+"] : "+
            "Then I'm going to construct outputs that look like:\n",
            file=sys.stderr)
        for each in outputs_array:
            print("  - With ID of : "+each[0]+"\n"+
                "    and the sequence is the group(s) : "+each[1],file=sys.stderr)

    # If it's omitted, we believe that means no filter, and we make it True
    # because it gets `eval`'d in the function. 
    if args.filter is None:
        args.filter = ["True == True"]

    if args.verbose >= 1:
        print("\n["+str(time.time())+"] : Then, I'm going to write out a "+
            args.output_format+" format file to "+
            args.output+"",file=sys.stderr)
        if args.report is not None:
            print("\n["+str(time.time())+"] : and a report to '"+
                vars(args)["report"]+".",file=sys.stderr)

# I should implement the below ... but really?
#    # checking file existance for outputs, zipping together the 
#    # output base with each of the three. 
#    exit_flag = 0
#    for each in zip( [vars(args)["output-base"]]*20,            \
#                    ["_fail.fastq", "_pass.fastq",              \
#                        "_report.fastq", "_report.csv" ] ):
#        # At this stage, the tuple is joined to make the filename
#        this_path = ''.join(each)
#        # If the write-report flag is off and the path is the report,
#        # then this won't trip True for that path existing
#import os.path
#        if os.path.isfile(this_path) and              \
#                not( not(args.write_report) and       \
#                    (this_path.find("_report")>=0) ):
#            print("\n"+"["+str(time.time())+"]"+" : "+"File "+this_path+
#                " exits, so I'm quitting before you ask me to do "+
#                "something you might regret.")
#            exit_flag = 1
#    if exit_flag == 1:
#        exit(1)

    # We begin
    if args.verbose >= 1:
        print("\n["+str(time.time())+"] : BEGIN RUNNING",file=sys.stderr)
    
    reader(
        vars(args)["input"],
        args.gzipped,
        vars(args)["output"],
        vars(args)["failed"],
        vars(args)["report"],
        operations_array, 
        args.filter ,
        outputs_array,
        args.output_format,
        args.verbose
        )

    print("\n"+"["+str(time.time())+"]"+" : "+
        "All worked 'till the work is done --- or some fatal error.",file=sys.stderr)

    exit(0)
