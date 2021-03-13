#!/usr/bin/env python3

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
import yaml

import regex
from Bio import SeqIO
from Bio import Seq, SeqRecord

##### Input configuration handling utilities

# IUPAC dictionary for translating codes to regex.
# Note the inclusion of * and + for repeats.
# from http://www.bioinformatics.org/sms/iupac.html
iupac_codes = { # only used for the configuration file input!
    'A':'A', 'C':'C', 'T':'T', 'G':'G',
    'R':'[AG]', 'Y':'[CT]', 'S':'[GC]', 'W':'[AT]',
    'K':'[GT]', 'M':'[AC]',
    'B':'[CGT]', 'D':'[AGT]', 'H':'[ACT]', 'V':'[ACG]',
    'N':'[ATCGN]', '*':'.*', '+':'.+' }


def config_from_file(file_path):
    """
    Tries to parse a configuration YAML file, and form a dictionary to pass into
    the main itermae reader function.
    """

    configuration = {}
    # Verbosely attempt to read it, and I want a hard exit if it's not parsed
    try:
        with open(file_path,'r') as f:
            config = yaml.load(f,Loader=yaml.SafeLoader)
    except:
        print('I failed to parse the supplied YAML file path name.',
            file=sys.stderr)
        raise
    # Looking for verbosity instruction global, if not global, then in 'outputs'
    try:
        try:
            verbosity = config['verbosity']
        except:
            verbosity = config['output']['verbosity']
    except:
        verbosity = 0 # else, just keep it bottled it up and tell no one 0_0
    configuration['verbosity'] = verbosity
    # Immediately use that verbostiy
    if verbosity >= 1:
        print("Reading and processing the configuration file '"+
            str(file_path)+"'.",file=sys.stderr)

    # Building array of matches objects, so input and compiled regex
    matches_array = []
    if verbosity >= 1:
        print("Processing each match:",file=sys.stderr)
    for each in config['matches']:
        try:
            each['use']
        except:
            each['use'] = 'input'
        if verbosity >= 1:
            print("    Taking '"+each['use']+"'. \n", end="",file=sys.stderr)
        if len(re.sub(r'(.)\1+',r'\1',each['marking'])) > len(set(each['marking'])):
            print("Error in reading yaml config! "+
                "It looks like you've repeated a group marking "+
                "character to match in multiple places. I do not support "+
                "that, use a different character.",file=sys.stderr)
            raise
        if len(each['pattern']) != len(each['marking']):
            print("Error in reading yaml config! "+
                "The pattern and marking you've defined are of "+
                "different lengths. I need them to be the same length.",
                file=sys.stderr)
            raise
        regex_groups = dict()
        group_order = list() # This is to keep track of the order in which
            # the groups are being defined in the paired lines
        for character, mark in zip(each['pattern'],each['marking']):
            if mark not in group_order:
                group_order.append(mark)
            try:
                regex_groups[mark] += iupac_codes[character.upper()]
                    # This is adding on the pattern for a certain marked
                    # matching group, as zipped above, and we're using
                    # IUPAC codes to turn ambiguity codes into ranges
                    # Note that it is converted to upper case!
            except:
                regex_groups[mark] = iupac_codes[character.upper()]
        regex_string = '' # building this now
        i = 0 # this is for keeping track of the untitled group numbering
        for mark in group_order:
            if 'name' in each['marked_groups'][mark].keys():
                check_reserved_name(each['marked_groups'][mark]['name'])
            else:
                each['marked_groups'][mark]['name'] = "untitled_group"+str(i)
                i += 1
            if verbosity >= 1:
                print("        Found group '"+mark+"' with pattern '"+
                    regex_groups[mark]+"'",end="",file=sys.stderr)
            try: # trying to build a repeat range, if supplied
                if 'repeat_min' not in each['marked_groups'][mark].keys():
                    each['marked_groups'][mark]['repeat_min'] = \
                        each['marked_groups'][mark]['repeat']
                if 'repeat_max' not in each['marked_groups'][mark].keys():
                    each['marked_groups'][mark]['repeat_max'] = \
                        each['marked_groups'][mark]['repeat']
                regex_groups[mark] = ('('+regex_groups[mark]+')'+
                    '{'+str(each['marked_groups'][mark]['repeat_min'])+','+
                        str(each['marked_groups'][mark]['repeat_max'])+'}'
                    )
                if verbosity >= 1:
                    print(", repeated between "+
                        str(each['marked_groups'][mark]['repeat_min'])+
                        " and "+
                        str(each['marked_groups'][mark]['repeat_max'])+
                        " times",end="",file=sys.stderr)
            except:
                pass
            error_array = [] # Then building the error tolerance spec
            try: 
                error_array.append(
                    "e<="+str(each['marked_groups'][mark]['allowed_errors']) )
            except:
                pass # This part takes up so much room because of try excepts...
            try: 
                error_array.append(
                    "i<="+str(each['marked_groups'][mark]['allowed_insertions']) )
            except:
                pass
            try: 
                error_array.append(
                    "d<="+str(each['marked_groups'][mark]['allowed_deletions']) )
            except:
                pass
            try: 
                error_array.append(
                    "s<="+str(each['marked_groups'][mark]['allowed_substitutions']) )
            except:
                pass
            if len(error_array):
                error_string = "{"+','.join(error_array)+"}"
            else:
                error_string = ""
            if verbosity >= 1:
                print(".\n",end="",file=sys.stderr)
            regex_string += ( "(?<"+each['marked_groups'][mark]['name']+
                ">"+regex_groups[mark]+")"+error_string )
        # Okay, then use the built up regex_string to compile it
        compiled_regex = regex.compile( regex_string, regex.BESTMATCH )
        # And save it with the input source used, in array
        matches_array.append( {'input':each['use'], 'regex':compiled_regex} )
    configuration['matches'] = matches_array

    if verbosity >= 1:
        print("Processing output specifications.",file=sys.stderr)
    output_list = config['output']['list'] # I do need some outputs, or fail
    outputs_array = [] 
    i = 0 # this is for naming untitled outputs sequentially
    for each in output_list:
        try:
            each['id']
        except:
            each['id'] = 'input.id' # default, the input.id
        try:
            each['name']
        except:
            each['name'] = 'untitled_output_'+str(i)
            i += 1
        try:
            each['filter']
        except:
            each['filter'] = 'True' # so will pass if not provided
        if verbosity >= 1:
            print("    Parsing output specification of '"+each['name']+"', "+
                "ID is '"+each['id']+"', filter outputs it if '"+
                each['filter']+"', with sequence derived of '"+
                each['seq']+"'.",file=sys.stderr)
        outputs_array.append( {
                'name':each['name'],
                'filter':[ each['filter'],
                    compile(each['filter'],'<string>','eval',optimize=2) ],
                'id':[ each['id'], 
                    compile(each['id'],'<string>','eval',optimize=2) ],
                'seq':[ each['seq'],
                    compile(each['seq'],'<string>','eval',optimize=2) ]
            })
    configuration['output_groups'] = outputs_array

    # Passing through rest or setting defaults
    try:
        configuration['input'] = config['input']['from']
    except:
        configuration['input'] = 'STDIN'
    try:
        configuration['input_format'] = config['input']['format']
    except:
        configuration['input_format'] = 'fastq'
    try:
        configuration['input_gzipped'] = config['input']['gzipped']
    except:
        configuration['input_gzipped'] = False
    try:
        configuration['output'] = config['output']['to']
    except:
        configuration['output'] = 'STDOUT'
    try:
        configuration['output_format'] = config['output']['format']
    except:
        configuration['output_format'] = 'sam'
    try:
        configuration['failed'] = config['output']['failed']
    except:
        configuration['failed'] = None
    try:
        configuration['report'] = config['output']['report']
    except:
        configuration['report'] = None

    return configuration


def config_from_args(args_copy):
    """
    Make configuration object from arguments provided. Should be the same as 
    the config_from_yaml output, if supplied the same.
    """

    configuration = {}
    verbosity = configuration['verbosity'] = args_copy.verbose

    # Make matches array
    try:
        matches_array = []
        for each in args_copy.match:
            for capture_name in re.findall('<(.*?)>',each):
                check_reserved_name(capture_name)
            try:
                (input_string, regex_string) = re.split("\s>\s",each.strip())
            except:
                input_string = 'input' # default to just use raw input
                regex_string = each.strip()
            compiled_regex = regex.compile(
                regex_string.strip(), # We use this regex
                regex.BESTMATCH # And we use the BESTMATCH strategy, I think
                )
            matches_array.append( {'input':input_string.strip(), 'regex':compiled_regex} )
    except:
        print("I failed to build matches array from the arguments supplied.",
            file=sys.stderr)
        raise
    configuration['matches'] = matches_array

    # Adding in defaults for outputs, may be redundant with argparse settings...
    if args_copy.output_id is None:
        args_copy.output_id = ['input.id']
    if args_copy.output_filter is None:
        args_copy.output_filter = ['True']

    # Normalizing all singletons to same length
    maximum_number_of_outputs = max( [len(args_copy.output_id), 
        len(args_copy.output_seq), len(args_copy.output_filter)] )
    if len(args_copy.output_id) == 1:
        args_copy.output_id = args_copy.output_id * maximum_number_of_outputs
    if len(args_copy.output_seq) == 1:
        args_copy.output_seq = args_copy.output_seq * maximum_number_of_outputs
    if len(args_copy.output_filter) == 1:
        args_copy.output_filter = args_copy.output_filter * maximum_number_of_outputs
    if not len(args_copy.output_id) == len(args_copy.output_seq) == len(args_copy.output_filter):
        print("The output IDs, seqs, and filters are of unequal sizes. "+
            "Make them equal, or only define one (and it will be reused "+
            "across all).",file=sys.stderr)
        raise

    try:
        i = 0
        outputs_array = [] 
        for idz, seqz, filterz in zip(args_copy.output_id,args_copy.output_seq,args_copy.output_filter):
            this_name = 'output_'+str(i)
            i += 1
            outputs_array.append( {   
                    'name': this_name,
                    'filter': [ filterz,
                        compile(filterz,'<string>','eval',optimize=2) ],
                    'id': [ idz, compile(idz,'<string>','eval',optimize=2) ],
                    'seq': [ seqz, compile(seqz,'<string>','eval',optimize=2) ] 
                })
    except:
        print("I failed to build outputs array from the arguments supplied.",
            file=sys.stderr)
        raise
    configuration['output_groups'] = outputs_array
    
    # Passing through the rest, defaults should be set in argparse defs
    configuration['input'] = args_copy.input
    configuration['input_gzipped'] = args_copy.gzipped
    configuration['input_format'] = args_copy.input_format
    configuration['output'] = args_copy.output
    configuration['output_format'] = args_copy.output_format
    configuration['failed'] = args_copy.failed
    configuration['report'] = args_copy.report
 
    return configuration


def check_reserved_name(name,reserved_names=['dummyspacer','input']):
    """
    This is just to check that the name is not one of these, if so, error out.
    """
    if name in reserved_names:
        print("Hey, you can't name a capture group "+
            (" or ".join(reserved_names[ [(i == name) for i in reserved_names]]))+
            ", I'm using that/those! Pick a different name.",
            file=sys.stderr)
        raise


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
        return str(self.substitutions)+"_"+str(self.insertions)+"_"+\
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
        # Note that the input is made uppercase!
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
            return eval(output_dict['filter'][1],globals(),self.context_filter)
        except:
            if self.verbosity >= 3:
                print("\n["+str(time.time())+"] : This read "+
                    self.seqs['input'].id+" failed to evaluate the filter "+
                    str(output_dict['filter'][0]),file=sys.stderr)
            return False

    def build_output(self,output_dict):
        """
        This builds the output
        """

        try:
            out_seq = SeqRecord.SeqRecord(Seq.Seq(""))
            out_seq = eval(output_dict['seq'][1],globals(),self.context_seq)
            out_seq.id = str(eval(output_dict['id'][1],globals(),self.context_seq))
            return out_seq
        except:
            if self.verbosity >= 3:
                print("\n["+str(time.time())+"] : This read "+
                    self.seqs['input'].id+" failed to build the output "+
                    "id: "+str(output_dict['id'][0])+
                    "seq: "+str(output_dict['seq'][0]) ,file=sys.stderr)
            return None

    def format_report(self,label,output_seq):
        """
        This is for formatting a standard report line for the reporting function
        """

        if output_seq is None:
            output_seq = SeqRecord.SeqRecord('X',
                id='ERROR',
                letter_annotations={'phred_quality':[0]})

        try:
            output_string = ( str(output_seq.id)+"\",\""+
                str(output_seq.seq)+"\",\""+
                phred_number_array_to_joined_string(
                    output_seq.letter_annotations['phred_quality']) )
        except:
            output_string = "*,*,*"

        return ( "\""+label+"\",\""+
            str(self.seqs['input'].id)+"\",\""+
            str(self.seqs['input'].seq)+"\",\""+
            phred_number_array_to_joined_string(self.seqs['input'].letter_annotations['phred_quality'])+"\",\""+
            output_string+"\",\""+
            "-".join([ i+"_"+self.group_stats[i].flatten() 
                        for i in self.group_stats ] )+
            "\"" ) # See group_stats method for what these are (start stop len)


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


def phred_letter_to_number(x):
    return ord(x)-33


def phred_number_to_letter(x):
    return chr(x+33)


def phred_number_array_to_joined_string(x):
    return str("".join([ phred_number_to_letter(i) for i in x]))


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
            letter_annotations={'phred_quality':
                [phred_letter_to_number(i) for i in fields[10]]},
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
            raise
        else:
            return sys.stdin
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
        failed_fh = open_output_fh(configuration['failed'])
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
                phred_number_array_to_joined_string(seq.letter_annotations['phred_quality']),
                "IE:Z:"+str(which) ),file=fh)
    elif format == "txt":
        print( str(seq.seq), file=fh)
    else:
        SeqIO.write(seq, fh, format) 


def chop( seq_holder, operations_array, outputs_array, out_format, input_format,
    output_fh, failed_fh, report_fh, verbosity ):
    """
    This one takes each record, applies the operations, evaluates the filters,
    generates outputs, and writes them to output handles as appropriate.
    """

    # If qualities are missing, add them as just 40
    if 'phred_quality' not in seq_holder.seqs['input'].letter_annotations.keys():
        seq_holder.seqs['input'].letter_annotations['phred_quality'] = [40]*len(seq_holder.seqs['input'])

        if verbosity >= 2:
            print("\n["+str(time.time())+"] : adding missing qualities of 40 "+
                "to sequence.", file=sys.stderr)

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

    # Next, iterate through the matches, applying each one
    for operation_number, operation in enumerate(operations_array):

        seq_holder.apply_operation( 'match_'+str(operation_number),
                operation['input'], operation['regex'] )

    # Now seq_holder should have a lot of matches, match scores and group stats,
    # and matched sequences groups. All these values allow us to apply filters
    # We unpack matches and scores into an internal environment for the filters
    seq_holder.build_context()

    # Then we eval the filters and build outputs, for each output
    output_records = []
    for each_output in outputs_array:
        output_records.append( { 
                'name': each_output['name'],
                'filter_result': seq_holder.evaluate_filter_of_output(each_output), 
                'output': seq_holder.build_output(each_output) 
            } )

    # This is just if we pass all the filters provided
    passed_filters = not any( 
            [ i['filter_result'] == False for i in output_records ] )

    # Then we can make the report CSV if asked for (mainly for debugging/tuning)
    if report_fh != None:
        for output_record in output_records:
            if output_record['filter_result']:
                print( seq_holder.format_report( 
                        "PassedFilterFor_"+output_record['name'], 
                        output_record['output'] ) ,file=report_fh)
            else:
                print( seq_holder.format_report( 
                        "FailedFilterFor_"+output_record['name'], 
                        output_record['output'] ) ,file=report_fh)

    # Finally, write all the outputs, to main stream if passed, otherwise to
    # the failed output (if provided)
    for output_record in output_records:
        if output_record['filter_result'] and output_record['output'] is not None:
            write_out_seq(output_record['output'], output_fh, out_format, 
                output_record['name'])
        elif failed_fh != None:
            write_out_seq(seq_holder.seqs['input'], failed_fh, input_format, 
                output_record['name'])

