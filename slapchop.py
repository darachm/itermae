#!/usr/bin/env python3

# 
# description to be written ...
#

import re
import json
import regex
import os.path
import argparse
import itertools
import multiprocessing 
from Bio import Seq, SeqRecord

#### The reader function, the thing paralleled
    # It takes a bunch of arguments
def reader(
    input_line_queue, input_fastq,
    pass_lock,   pass_fastq,
    fail_lock,   fail_fastq,
    report_lock, report_csv,
    bite_size,   limit_size, if_write_report, verbosity,
    operations_array, filters ,
    output_seq_spec, output_id_spec
    ):

    # Always looping
    while True:

        # Waits for the input line marker from the queue
        current_pos = input_line_queue.get()

        # If it's a exitpill, declare the end and return.
        # It's not a poisonpill, because some people have lost
        # actual people to suicide.
        if current_pos == "exitpill":
            if verbosity > 1:
                print("\n"+multiprocessing.current_process().name+
                    " with pid "+
                    str(multiprocessing.current_process().pid)+
                    " found a exit pill and is exiting"
                    )
            input_line_queue.put("exitpill")
            return(0)

        # Next we try to use the limit_size to see if we should stop
        if limit_size is not None:
            if current_pos > limit_size:
                if verbosity > 1:
                    print("\n"+multiprocessing.current_process().name+
                        " with pid "+
                        str(multiprocessing.current_process().pid)+
                        " exceeded the limit and is exiting"
                        )
                input_line_queue.put("exitpill")
                return(0)

        # Tell us what you're doing?
        if verbosity > 1:
            print("\n"+multiprocessing.current_process().name+
                " trying to read a chunk of size "+
                str(bite_size)+" fastq records."
                )

        # If we have the current_pos from above, then no-one else
        # is reading the file, so let's open it
        ifqp = open(input_fastq,"r")
        # Go to that position
        ifqp.seek(current_pos)
        # And read a chunk. We use an iterator because we have to be
        # line-oriented but the position marker is byte-oriented.
        chunk = []
        for i in range(bite_size*4):
            chunk.append(ifqp.readline())
        # We read back the current position in bytes
        current_pos = ifqp.tell()
        # And detect if we're at the end of the file, if so, exitpill
        if ifqp.readline() == "":
            if verbosity > 1:
                print("\n"+multiprocessing.current_process().name+
                    " with pid "+
                    str(multiprocessing.current_process().pid)+
                    " found the end and is done"
                    )
            input_line_queue.put("exitpill")
            return(0)
        # Otherwise, then just put the current position back in the
        # queue
        else:
            input_line_queue.put(current_pos)

        # Let's make some arrays to hold three kinds of records.
        # First is those that pass filters, then those that fail.
        # The last is the report if you've made that option
        pass_records = []
        fail_records = []
        # We make report_records empty even if not using it
        report_records = []

        # We use itertools to slice 4 line chunks (this is fastq)
        for slice_base in itertools.islice(range(len(chunk)),0,None,4):

            # if it's empty, we must have butted against the end
            # of the file, but we already inserted the exitpill
            if chunk[slice_base] == "":
                break
    
            # Otherwise, call chop and save the returns 
            (passed, output_record, report_object) = \
                chop(
                    # This is the actual fastq record
                    chunk[slice(slice_base,(slice_base+4))],
                    # This is the operations to do, filters, 
                    # then if to save report records,
                    # and how verbose to be
                    operations_array, filters,
                    output_seq_spec, output_id_spec,
                    if_write_report, verbosity
                    )

            # This is a per-record test
            if passed:
                pass_records.append(output_record)
            else:
                fail_records.append(output_record)
            if if_write_report:
                report_records.append(report_object)
        
        # Then we write out the records with the appropriate locks
        with pass_lock:
            with open(pass_fastq,"a") as f:
                for i in pass_records:
                    print(str(i.id)+"\n"+str(i.seq)+"\n"+"+"+"\n"+
                            i.letter_annotations['phred_quality']+"\n",
                        file=f)

        with fail_lock:
            with open(fail_fastq,"a") as f:
                for i in fail_records:
                    print(str(i.id)+"\n"+str(i.seq)+"\n"+"+"+"\n"+
                            i.letter_annotations['phred_quality']+"\n",
                        file=f)

        with report_lock:
            with open(report_csv,"a") as f:
                for i in report_records:
                    print(i,file=f)



def chop(
    record, operations_array, filters,
    output_seq_spec, output_id_spec,
    if_write_report = False, verbosity=0
    ):

    # Making the input record from the raw strings
    input_record = SeqRecord.SeqRecord(Seq.Seq(record[1].rstrip()),
        id = record[0].rstrip().split(" ")[0])
    input_record.letter_annotations['phred_quality'] = \
        record[3][0:len(record[1].rstrip())]

    # We make some holders for these operations
    scores_holder = dict()
    seq_holder = dict()
    seq_holder['input'] = input_record 

#rewrite as a class ????

    # Chop grained verbosity
    if verbosity > 2:
        print("\n"+multiprocessing.current_process().name+
            " starting to process :\n  "+
            input_record.id+"\n  "+
            input_record.seq+"\n  "+
            input_record.letter_annotations['phred_quality']
            )

    for each_operation in operations_array:

        # The first element is the name, the next two are used later
        operation_name = each_operation[0]
        operation = each_operation[1:]

        # This should fail if you didn't specify anything taking 
        # from input stream!
        if operation[0] not in seq_holder.keys():
            if verbosity > 3:
                print("\n"+multiprocessing.current_process().name+
                " can't find the sequence named `"+
                operation[0]+"` in the holder, so continuing."
                )
            continue

        if verbosity > 3:
            print("\n"+multiprocessing.current_process().name+
                " attempting to match : "+operation[1]+
                " against "+seq_holder[operation[0]].seq
                )

        # Here we execute the actual meat of the business
        fuzzy_match = regex.search(
            # We use this regex
            operation[1], 
            # to search on this sequence
            str(seq_holder[operation[0]].seq),
            # And we use the BESTMATCH strategy, I think
            regex.BESTMATCH )

        if verbosity > 3:
            print("\n"+multiprocessing.current_process().name+
                " match is : "+
                str(fuzzy_match)
                )

        # This is fine, just means the pattern couldn't match at all
        if fuzzy_match is None:
            continue
        # If we did match, then we store them in places
        else:
            # We use tuples to store all the details of the kinds
            # of errors that allowed the match
            (scores_holder[operation_name+'_substitutions'],
                scores_holder[operation_name+'_insertions'],
                scores_holder[operation_name+'_deletions']
                ) = fuzzy_match.fuzzy_counts
            # Then for each of the groups matched by the regex
            for match_name in fuzzy_match.groupdict():
                # We stick into the holder
                # a slice of the input seq, that is the matched
                # span of this matching group
                seq_holder[match_name] = \
                    seq_holder[operation[0]]\
                    [slice(*fuzzy_match.span(match_name))]
                # Then we record the start, end, and length of the
                # matched span
                (scores_holder[match_name+'_start'],
                    scores_holder[match_name+'_end']
                    ) = fuzzy_match.span(match_name)
                scores_holder[match_name+'_length'] = \
                    (scores_holder[match_name+'_end'] - 
                        scores_holder[match_name+'_start'])

    # All these values allow use to apply filters, using this
    # function
    evaluated_filters = evaluate_filters(filters, scores_holder)

    # This evaluated_filters should be logical list
    if not all(evaluated_filters):

        if verbosity > 2:
            print("\n"+multiprocessing.current_process().name+
                " evaluated the filters as : "+
                str(evaluated_filters)+
                " and so failed."
                )

        # So if we should write this per-record report
        if if_write_report:
            return((False,input_record,
                "\"FailedFilterOnThisInput\","+
                "\""+input_record.id+"\",\""+
                    input_record.seq+"\",\""+
                    re.sub("\"","\\\"",re.sub(",","\,",
                        json.dumps(scores_holder)))+"\""
                ))
        # If this json dump is empty, it might be because it didn't
        # ever match the first operation, so then just died without
        # building that object

    else:

        try:
            # We attempt to form the correct output record based on
            # the arguments given
            output_record = evaluate_output_directives(
                output_seq_spec,output_id_spec,seq_holder) 

            if verbosity > 2:
                print("\n"+multiprocessing.current_process().name+
                    " evaluated the filters as : "+
                    str(evaluated_filters)+
                    " and so passed!"
                    )

            # If we want to write the report, we make it
            if if_write_report:
                return((True,output_record,
                    "\"Passed\","+
                    "\""+output_record.id+"\",\""+
                        output_record.seq+"\",\""+
                        re.sub("\"","\\\"",re.sub(",","\,",
                            json.dumps(scores_holder)))+"\""
                    ))
            # Otherwise, just the record and if it passed
            else:
                return((True,output_record,""))

        except:

            if verbosity > 2:
                print("\n"+multiprocessing.current_process().name+
                    " failed upon forming the output."
                    )

            if if_write_report:
                return((False,input_record,
                    "\"FailedDirectivesToMakeOutputSeq\","+
                    "\""+input_record.id+"\",\""+
                        input_record.seq+"\",\""+
                        re.sub("\"","\\\"",re.sub(",","\,",
                            json.dumps(scores_holder)))+"\""
                    ))
            else:
                return((False,input_record,""))

def evaluate_output_directives(output_seq, output_id, seq_holder):
    # Here we evaluate them but using that dictionary as the global
    # dictionary, because done is better than dogma.
    return_record = eval(output_seq,{},seq_holder)
    return_record.id = eval(output_id,{},seq_holder)
    return(return_record)


def evaluate_filters(filters,scores_holder):
    return_object = []
    try:
        for each_filter in filters:
            # Here we evaluate them but using that dictionary as the
            # global dictionary, because done is better than dogma.
            if eval(each_filter,{},scores_holder):
                return_object.append(True)
            else:
                return([False])
    except:
        return([False])
    return(return_object)


if __name__ == '__main__':

#### defining arguments with argparse module

    # Name and description
    parser = argparse.ArgumentParser(description=""+
        "slapchop.py")
    parser.add_argument("input-fastq",
        help="The FASTQ formatted file to process.")

    # Resources details
    parser.add_argument("--processes",default=1)
    parser.add_argument("--bite-size",default=1000,
        help="The size of bites to chomp off of the input file for "+
            "multi-process, also the max size of the disk write "+
            "caching."
        )
    parser.add_argument("--limit",default=None,
        help="The limit of reads to process, useful for just "+
            "proofing that your operations are actually working, "+
            "and for collecting stats in the report for setting "+
            "filters."
        )

    # verbosity
    parser.add_argument("-v","--verbose",action="count",default=0,
        help=" 0 is nothing, "+
            "1 is setup messages and start-stop messsages, "+
            "2 is worker-level details, "+
            "3 is chop-level details, "+
            "4 is each operation level details."+
            "All to standard out, so be ready for it."
            )

    # Operations
    parser.add_argument("--operation","-o",action="append",
        help="The pattern to match and extract. This has a "+
            "specific and sort of complicated syntax. Refer to "+
            "the documentation via the README.md file."+
            "They are chained in the order you specify.")

    # Filter specification
    parser.add_argument("--filter","-f",action="append",
        help="A filter for eliminating reads that don't pass some "+
            "alignment based cutoff. This is specified per "+
            "operation, so remember the name from above. Syntax: ''")

    # Output stream
    parser.add_argument("--output-id",
        help="format for the output file id, per read that passes "+
            "filter",
        default="input.id")
    parser.add_argument("--output-seq",
        help="format for the output file seq, per read that "+
            "passes filter",
        default="input")

    # Output files
    parser.add_argument("output-base",
        help="Base name for the output, will be used to make, "+
            "for example: basename.fastq, basename.report")
    parser.add_argument("--write-report",action='store_true',
        help="Add this flag to print a report of per-read "+
            "statistics, that's a lot of disk writes btw, but "+
            "would be good in combination with a --limit argument "+
            "so that you can spec out the kinds of noise you got "+
            "in your data.")
    parser.add_argument("--maxQueueSize",help="in gigs",default=10)

#### Parse, clean up, and possibly report arguments

    args = parser.parse_args()

    # Convert to integer so I don't have to later
    args.bite_size = int(args.bite_size)
    try:
        args.limit = int(args.limit)
    except:
        pass

    # Operations are read as an array of dicts to keep it ordered
    operations_array = []

    try:
        for each in args.operation:
            # Split on the colon to specify the name on the left of it
            (name, instruction) = each.split(":")
            # similarly use the ` > ` to specify the flow of input to
            # the regex
            (input_string, regex_string) = instruction.strip().split(" > ")
            # append that to the operations array
            operations_array.append( [name, input_string, regex_string] )
    except:
        # Failure likely from lack of operations to do
        print("\n"+"Wait a second, there's no operations to be done! "+
            "What am I doing here? What is my purpose? "+
            "As they said in Darkstar: 'Let there be light.' "+
            "Exiting...")
        exit(1)

    if args.verbose > 0:
        print("\n"+"I'm reading in '"+vars(args)["input-fastq"]+"', "+
            "applying these operations of alignment :\n")
        for each in operations_array:
            print("- "+each[0]+" :\n"+
                "  from : "+each[1]+"\n"+
                "  extract groups with regex : '"+each[2])

    if args.verbose > 0:
        print("\n"+"...and with these filters:")
        try:
            for i in args.filter:
                print("  "+i)
        except:
            print("  ( no filters defined )")

    if args.verbose > 0:
        print("\n"+"Then, I'm going to write out a FASTQ file to '"+
            vars(args)["output-base"]+".fastq'",end="")
        if args.write_report:
            print(" and a report to '"+
                vars(args)["output-base"]+"_report.csv'",end="")
        print(".")

    if args.verbose > 0:
        print("\n"+"I will proceed to process the file with "+
            str(args.processes)+" processes operating in chunks of "+
            str(args.bite_size)+" records.")

        if args.limit is not None:
            print("But, I am limiting myself to the first "+
                str(args.limit)+" records for this run.")
        
    # checking file existance for outputs, zipping together the 
    # output base with each of the three. 
    exit_flag = 0
    for each in zip( [vars(args)["output-base"]]*20,            \
                    ["_fail.fastq", "_pass.fastq",              \
                        "_report.fastq", "_report.csv" ] ):
        # At this stage, the tuple is joined to make the filename
        this_path = ''.join(each)
        # If the write-report flag is off and the path is the report,
        # then this won't trip True for that path existing
        if os.path.isfile(this_path) and              \
                not( not(args.write_report) and       \
                    (this_path.find("_report")>=0) ):
            print("\n"+"File "+this_path+
                " exits, so I'm quitting before you ask me to do "+
                "something you might regret.")
            exit_flag = 1
    if exit_flag == 1:
        exit(1)

    # We begin
    if args.verbose > 0:
        print("\nBEGIN\n")
    
#### Running the actual functions using multiprocessing

    # Make a manager object
    manager  = multiprocessing.Manager()
    # Make a cue for holding the place in the input file
    input_line_queue = multiprocessing.Queue()
    # We start at 0
    input_line_queue.put(0)
    # We make a few locks, as a tuple
    (pass_lock, fail_lock, report_lock) = (
        manager.Lock(), manager.Lock(), manager.Lock() )

    # We start some jobs, one per process

    jobs = []

    for i in range(1,int(args.processes)+1):

        jobs.append(
            multiprocessing.Process(
                target=reader,
                args=(
                    input_line_queue, vars(args)["input-fastq"],
                    pass_lock, vars(args)["output-base"]+"_pass.fastq",
                    fail_lock, vars(args)["output-base"]+"_fail.fastq",
                    report_lock, vars(args)["output-base"]+"_report.csv",
                    args.bite_size, args.limit, 
                    args.write_report, args.verbose,
                    operations_array, args.filter ,
                    args.output_seq, args.output_id
                    ),
                name="Comrade"+str(i)
                )
            )
        if args.verbose > 0:
            print("\n"+"Comrade"+str(i)+" starting up...")

    for i in jobs:
        i.start()

    for i in jobs:
        if i.is_alive():
            i.join()

    print("\n"+
        "All worked 'till the work is done --- or some fatal error.")
    exit(0)
