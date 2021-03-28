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
        yield SeqRecord.SeqRecord( Seq.Seq(seq), id=seq, description="")


def fix_desc(seq_record):
    """
    According to https://github.com/biopython/biopython/issues/398 ,
    BioPython mimics an old weird behavior by outputting the ID in the
    description field. There's a fix for the FASTA writer, but not the
    FASTQ ... so here we munge that by removing the ID from the description.

    ... except that I don't know where to put this.
    """
    seq_record.description = re.sub(str(seq_record.id),"",seq_record.description)
    return seq_record


def write_out_seq(seq,fh,format,which):
    if format == "sam":
        print( format_sam_record( seq.id, str(seq.seq),
                phred_number_array_to_joined_string(seq.letter_annotations['phred_quality']),
                "IE:Z:"+str(which) ),file=fh)
        # We ignore printing the description anywhere - if you need it, concat
        # it onto the ID
    elif format == "txt":
        print( str(seq.seq), file=fh)
    else:
        SeqIO.write(seq, fh, format) 





class Configuration:
    """
    This is for configuring itermae, from YAML or CLI arguments.
    """

    def __init__(self):
        self.verbosity = 0
        self.matches_array = []
        self.outputs_array = []
        self.untitled_group_number = 0
        self.untitled_output_number = 0
        self.input = 'STDIN'
        self.input_format = 'fastq'
        self.gzipped = False
        self.to = 'STOUT'
        self.output_format = 'sam'
        self.failed = None
        self.report = None
        self.matches_array = []
        self.outputs_array = []
        self.output_fh = None
        self.failed_fh = None
        self.report_fh = None

        # IUPAC dictionary for translating codes to regex.
        # from http://www.bioinformatics.org/sms/iupac.html
        # Note the inclusion of * and + for repeats.
        self.iupac_codes = { # only used for the configuration file input!
            'A':'A', 'C':'C', 'T':'T', 'G':'G',
            'R':'[AG]', 'Y':'[CT]', 'S':'[GC]', 'W':'[AT]',
            'K':'[GT]', 'M':'[AC]',
            'B':'[CGT]', 'D':'[AGT]', 'H':'[ACT]', 'V':'[ACG]',
            'N':'[ATCGN]', '*':'.*', '+':'.+' }

    def open_input_fh(self):
        if self.input.upper() == 'STDIN':
            if self.gzipped:
                print("I can't handle gzipped inputs on STDIN ! "+
                    "You shouldn't see this error, it shoulda been caught in "+
                    "the launcher script.",file=sys.stderr) 
                raise
            else:
                self.input_fh = sys.stdin
        else:
            if self.gzipped:
                self.input_fh = gzip.open(self.input,'rt',encoding='ascii')
            else:
                self.input_fh = open(self.input,'rt')

    def open_appropriate_input_format(self):
        if   self.input_format == 'fastq':
            self.input_seqs = SeqIO.parse(self.input_fh, self.input_format)
        elif self.input_format == 'sam':
            self.input_seqs = iter(read_sam_file(self.input_fh))
        elif self.input_format == 'fasta':
            self.input_seqs = SeqIO.parse(self.input_fh, self.input_format)
        elif self.input_format == 'txt':
            self.input_seqs = iter(read_txt_file(self.input_fh))
        else:
            print("I don't know that input file format name '"+self.input_format+
                "'. I will try and use the provided format name in BioPython "+
                "SeqIO, and we will find out together if that works.",
                file=sys.stderr) 
            self.input_seqs = SeqIO.parse(self.input_fh, self.input_format)

    def get_input_seqs(self):
        """
        This calls 
        to set the `input_seqs` attribute to be an iterator of BioPython
        sequence records.
        """
        self.open_input_fh()
        self.open_appropriate_input_format()

    def open_output_fh(self,file_string):
        if file_string is None:
            return None
        if file_string.upper() == 'STDOUT':
            return sys.stdout
        elif file_string.upper() == 'STDERR':
            return sys.stderr
        else:
            return open(file_string,'a')

    def open_output(self):
        self.output_fh = self.open_output_fh(self.output)

    def open_report(self):
        self.report_fh = self.open_output_fh(self.report)

    def open_failed(self):
        self.failed_fh = self.open_output_fh(self.failed)

    def close_fhs(self):
        for i in [ self.input_seqs, self.output_fh, self.failed_fh, self.report_fh] :
            try:
                i.close()
            except:
                pass

    def check_reserved_name(self,name,
            reserved_names=['dummyspacer','input','id','description'] ):
        """
        This is just to check that the name is not one of these, if so, error out.
        - 'dummyspacer' is so you can pop an X into your sequence as a separator
            delimiter for later processing
        - 'input' is the input group, the original one
        - 'id' is the input ID, here just as 'id' so it's easy to find
        - 'description' is for mapping over the FASTQ description
        """
        if name in reserved_names:
            print("Hey, you can't name a capture group "+
                (" or ".join(reserved_names[ [(i == name) for i in reserved_names]]))+
                ", I'm using that/those! Pick a different name.",
                file=sys.stderr)
            raise

    def config_from_file(self,file_path):
        """
        Tries to parse a configuration YAML file to update this configuration
        object. Recommend you run this config first, then config_from_args.
        Pass in the file path as an argument.
        """

        if file_path == False:
            return
    
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
                self.verbosity = config['verbosity']
            except:
                self.verbosity = config['output']['verbosity']
        except:
            pass

        # Immediately use that verbostiy
        if self.verbosity >= 1:
            print("Reading and processing the configuration file '"+
                str(file_path)+"'.",file=sys.stderr)
    
        # Building array of matches objects, so input and compiled regex
        if self.verbosity >= 1:
            print("Processing each match:",file=sys.stderr)
        for each in config['matches']:
            try:
                each['use']
            except:
                each['use'] = 'input'
            if self.verbosity >= 1:
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
            pattern_groups = dict()
            group_order = list() # This is to keep track of the order in which
                # the groups are being defined in the paired lines
            for character, mark in zip(each['pattern'],each['marking']):
                if mark not in group_order:
                    group_order.append(mark)
                try:
                    pattern_groups[mark] += character.upper()
                except:
                    pattern_groups[mark] = character.upper()

            regex_string = '' # building this now
            for mark in group_order:

                if 'name' in each['marked_groups'][mark].keys():
                    self.check_reserved_name(each['marked_groups'][mark]['name'])
                else:
                    each['marked_groups'][mark]['name'] = "untitled_group"+\
                        str(self.untitled_group_number)
                    self.untitled_group_number += 1

                pattern_string = ""
                if len(set(pattern_groups[mark])) == 1:
                    pattern_string = self.iupac_codes[pattern_groups[mark][0].upper()]
                else:
                    for character in pattern_groups[mark]:
                        pattern_string += self.iupac_codes[character.upper()]
                        # This is adding on the pattern for a certain marked
                        # matching group, as zipped above, and we're using
                        # IUPAC codes to turn ambiguity codes into ranges
                        # Note that it is converted to upper case!

                if self.verbosity >= 1:
                    print("        Found group '"+mark+"' with pattern '"+
                        pattern_string+"'",end="",file=sys.stderr)

                try: # trying to build a repeat range, if supplied
                    if 'repeat_min' not in each['marked_groups'][mark].keys():
                        each['marked_groups'][mark]['repeat_min'] = \
                            each['marked_groups'][mark]['repeat']
                    if 'repeat_max' not in each['marked_groups'][mark].keys():
                        each['marked_groups'][mark]['repeat_max'] = \
                            each['marked_groups'][mark]['repeat']
                    pattern_string = ('('+pattern_string+')'+
                        '{'+str(each['marked_groups'][mark]['repeat_min'])+','+
                            str(each['marked_groups'][mark]['repeat_max'])+'}'
                        )
                    if self.verbosity >= 1:
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
                if self.verbosity >= 1:
                    print(".\n",end="",file=sys.stderr)

                regex_string += ( "(?<"+each['marked_groups'][mark]['name']+
                    ">"+pattern_string+")"+error_string )

            # Okay, then use the built up regex_string to compile it
            compiled_regex = regex.compile( regex_string, regex.BESTMATCH )
            # And save it with the input source used, in array
            self.matches_array.append( {'input':each['use'], 'regex':compiled_regex} )
    
        if self.verbosity >= 1:
            print("Processing output specifications.",file=sys.stderr)

        output_list = config['output']['list'] # I do need some outputs, or fail
        for each in output_list:

            try:
                each['id']
            except:
                each['id'] = 'id' # default, the id
            try:
                each['description']
            except:
                each['description'] = 'description' # default pass through from in
            try:
                each['name']
            except:
                each['name'] = 'untitled_output_'+str(self.untitled_output_number)
                self.untitled_output_number += 1
            try:
                each['filter']
            except:
                each['filter'] = 'True' # so will pass if not provided

            if self.verbosity >= 1:
                print("    Parsing output specification of '"+each['name']+"', "+
                    "ID is '"+each['id']+"' (input ID is 'id'), filter outputs "+
                    "to accept only if '"+each['filter']+"' is True, with "+
                    "sequence derived from '"+each['seq']+"', and a description "+
                    "of '"+each['description']+"' ('description' is input "+
                    "description').",file=sys.stderr)

            self.outputs_array.append( {
                    'name':each['name'],
                    'filter':[ each['filter'],
                        compile(each['filter'],'<string>','eval',optimize=2) ],
                    'id':[ each['id'], 
                        compile(each['id'],'<string>','eval',optimize=2) ],
                    'seq':[ each['seq'],
                        compile(each['seq'],'<string>','eval',optimize=2) ],
                    'description':[ each['description'],
                        compile(each['description'],'<string>','eval',optimize=2) ]
                })

    def config_from_args(self,args_copy):
        """
        Make configuration object from arguments provided. Should be the same as 
        the config_from_yaml output, if supplied the same.
        """
    
        if args_copy.verbose:
            self.verbosity = args_copy.verbose

        for each in args_copy.match:
            try:
                for capture_name in re.findall('<(.*?)>',each):
                    self.check_reserved_name(capture_name)
                try:
                    (input_string, regex_string) = re.split("\s>\s",each.strip())
                except:
                    input_string = 'input' # default to just use raw input
                    regex_string = each.strip()
                compiled_regex = regex.compile(
                    regex_string.strip(), # We use this regex
                    regex.BESTMATCH # And we use the BESTMATCH strategy, I think
                    )
                self.matches_array.append( {'input':input_string.strip(), 'regex':compiled_regex} )
            except:
                print("I failed to build matches array from the arguments supplied.",
                    file=sys.stderr)
                raise

        # Adding in defaults for outputs. Can't do that with argparse, I think,
        # because this needs to be appending
        maximum_number_of_outputs = max( [len(args_copy.output_id), 
            len(args_copy.output_seq), len(args_copy.output_filter),
            len(args_copy.output_description)] )

        if maximum_number_of_outputs:
            if args_copy.output_id is []:
                args_copy.output_id = ['id']
            if args_copy.output_filter is []:
                args_copy.output_filter = ['True']
            if args_copy.output_description is []:
                args_copy.output_description = ['description']

            # Normalizing all singletons to same length
            if len(args_copy.output_id) == 1:
                args_copy.output_id = args_copy.output_id * maximum_number_of_outputs
            if len(args_copy.output_seq) == 1:
                args_copy.output_seq = args_copy.output_seq * maximum_number_of_outputs
            if len(args_copy.output_filter) == 1:
                args_copy.output_filter = args_copy.output_filter * maximum_number_of_outputs
            if len(args_copy.output_description) == 1:
                args_copy.output_description = args_copy.output_description * maximum_number_of_outputs
            if not ( len(args_copy.output_id) == len(args_copy.output_seq) == 
                    len(args_copy.output_filter) == len(args_copy.output_description) ):
                print("The output IDs, seqs, descriptions, and filters are of unequal "+
                    "sizes. Make them equal, or only define one each and it will be "+
                    "reused across all.",file=sys.stderr)
                raise
    
            try:
                i = 0
                for idz, seqz, filterz, description in zip(args_copy.output_id, args_copy.output_seq, args_copy.output_filter, args_copy.output_description) :
                    this_name = 'output_'+str(i)
                    i += 1
                    self.outputs_array.append( {   
                            'name': this_name,
                            'filter': [ filterz,
                                compile(filterz,'<string>','eval',optimize=2) ],
                            'id': [ idz, compile(idz,'<string>','eval',optimize=2) ],
                            'seq': [ seqz, compile(seqz,'<string>','eval',optimize=2) ] ,
                            'description':[ description, compile(description,'<string>','eval',optimize=2) ]
                        })
            except:
                print("I failed to build outputs array from the arguments supplied.",
                    file=sys.stderr)
                raise
        
        # Passing through the rest, defaults should be set in argparse defs
        if args_copy.input is not None:
            self.input = args_copy.input
        if args_copy.input_format is not None:
            self.input_format = args_copy.input_format
        if args_copy.gzipped is not None:
            self.gzipped = args_copy.gzipped
        if args_copy.output is not None:
            self.output = args_copy.output
        if args_copy.output_format is not None:
            self.output_format = args_copy.output_format
        if args_copy.failed is not None:
            self.failed = args_copy.failed
        if args_copy.report is not None:
            self.report = args_copy.report

    def summary(self):
        return_string = ('Configured as:'+
            '\n    input from: '+self.input+
            '\n    input format: '+self.input_format+
            '\n    is it gzipped?: '+str(self.gzipped)+
            '\n    output APPENDING to: '+self.output+
            '\n    output format is: '+self.output_format+
            '\n    failed being APPENDED to file: '+str(self.failed)+
            '\n    report being APPENDED to file: '+str(self.report)+
            '\n    with verbosity set at: '+str(self.verbosity)+
            '\n    doing these matches:')
        for each in self.matches_array:
            return_string += '\n        input: '+each['input']
            return_string += '\n        regex: '+str(each['regex'])
        return_string += '\n    writing these outputs:'
        for each in self.outputs_array:
            return_string += '\n        id: '+str(each['id'][0])
            return_string += '\n        description: '+str(each['description'][0])
            return_string += '\n        seq: '+str(each['seq'][0])
            return_string += '\n        filter: '+str(each['filter'][0])
        return return_string

    def reader(self):
        """
        This reads inputs, calls the `chop` function on each one, and sorts it
        off to outputs. So this is called by the main function, and is mostly about
        handling the I/O. 
        """
    
        ### Open up file handles
    
        # Input
        self.get_input_seqs()
    
        # Outputs - passed records, failed records, report file
        self.open_output()
        self.open_failed()
        self.open_report()
    
        # Do the chop-ing...
        for each_seq in self.input_seqs:
    
            # CAUTION
            # The below is a munge. 
            # According to https://github.com/biopython/biopython/issues/398 ,
            # BioPython mimics an old tool's weird behavior by outputting the 
            # ID in the description field. The fix for it relies on a comparing
            # a white-space 'split' to remove the ID if it's in the description.
            # So that doesn't work if you modify the ID or so, so I remove right
            # after parsing.
            each_seq.description = re.sub(str(each_seq.id),"",
                each_seq.description).lstrip()

            seq_holder = SeqHolder(each_seq,configuration=self)
            seq_holder.chop()
    
        self.close_fhs()
    
        return(0)



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
    def __init__(self, input_record, configuration):
        # So the .seqs holds the sequences accessed by the matching, and there's
        # a dummyspacer in there just for making outputs where you want that
        # for later partitioning. Input is input.
        self.seqs = {
            'dummyspacer': SeqRecord.SeqRecord(Seq.Seq("X"),id="dummyspacer"),
            'input': input_record }
        self.seqs['dummyspacer'].letter_annotations['phred_quality'] = [40]
        self.configuration = configuration
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

        if self.configuration.verbosity >= 3:
            print("\n["+str(time.time())+"] : attempting to match : "+
                str(regex)+" against "+self.seqs[input_group].seq,
                file=sys.stderr)

        # Here we execute the actual meat of the business.
        # Note that the input is made uppercase!
        fuzzy_match = regex.search( str(self.seqs[input_group].seq).upper() )

        if self.configuration.verbosity >= 3:
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

                #self.seqs[match_name].description = "" 
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

        # Then one for the IDs, so we're setting the input ID as 'id', and then
        # each group name just refers to the sequence. I assume folks are not
        # wanting to be putting seq qualities in the ID. We do make 
        # 'description' available if that's important
        self.context_id = { 
            'id': self.seqs['input'].id , 
            'description': self.seqs['input'].description , 
            **{ i: str(self.seqs[i].seq) for i in self.seqs } }

    def evaluate_filter_of_output(self,output_dict):
        """
        This tests a defined filter on the 'seq_holder' object
        """

        try:
            return eval(output_dict['filter'][1],globals(),self.context_filter)
        except:
            if self.configuration.verbosity >= 3:
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
            out_seq.id = str(eval(output_dict['id'][1],globals(),self.context_id))
            out_seq.description = str(eval(output_dict['description'][1],globals(),self.context_id))
            return out_seq
        except:
            if self.configuration.verbosity >= 3:
                print("\n["+str(time.time())+"] : This read "+
                    self.seqs['input'].id+" failed to build the output of "+
                    "id: '"+str(output_dict['id'][0])+"', and "+
                    "seq: '"+str(output_dict['seq'][0])+"'." ,file=sys.stderr)
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

    def chop(self):
        """
        This one takes each record, applies the operations, evaluates the filters,
        generates outputs, and writes them to output handles as appropriate.
        """
    
        # If qualities are missing, add them as just 40
        if 'phred_quality' not in self.seqs['input'].letter_annotations.keys():
            self.seqs['input'].letter_annotations['phred_quality'] = [40]*len(self.seqs['input'])
    
            if self.configuration.verbosity >= 2:
                print("\n["+str(time.time())+"] : adding missing qualities of 40 "+
                    "to sequence.", file=sys.stderr)
    
        # For chop grained self.configuration.verbosity, report
        if self.configuration.verbosity >= 2:
            print("\n["+str(time.time())+"] : starting to process : "+
                self.seqs['input'].id+"\n  "+self.seqs['input'].seq+"\n  "+ 
                str(self.seqs['input'].letter_annotations['phred_quality']),
                file=sys.stderr)
    
        # This should fail if you didn't specify anything taking from input stream!
        assert self.configuration.matches_array[0]['input'] == "input", (
            "can't find the sequence named `input`, rather we see `"+
            self.configuration.matches_array[0]['input']+"` in the holder, so breaking. You should "+
            "have the first operation start with `input` as a source." )
    
        # Next, iterate through the matches, applying each one
        for operation_number, operation in enumerate(self.configuration.matches_array):
    
            self.apply_operation( 'match_'+str(operation_number),
                    operation['input'], operation['regex'] )
    
        # Now self should have a lot of matches, match scores and group stats,
        # and matched sequences groups. All these values allow us to apply filters
        # We unpack matches and scores into an internal environment for the filters
        self.build_context()
    
        # Then we eval the filters and build outputs, for each output
        output_records = []
        for each_output in self.configuration.outputs_array:
            output_records.append( { 
                    'name': each_output['name'],
                    'filter_result': self.evaluate_filter_of_output(each_output), 
                    'output': self.build_output(each_output) 
                } )
    
        # This is just if we pass all the filters provided
        passed_filters = not any( 
                [ i['filter_result'] == False for i in output_records ] )
    
        # Then we can make the report CSV if asked for (mainly for debugging/tuning)
        if self.configuration.report_fh != None:
            for output_record in output_records:
                if output_record['filter_result']:
                    print( self.format_report( 
                            "PassedFilterFor_"+output_record['name'], 
                            output_record['output'] ) ,file=self.configuration.report_fh)
                else:
                    print( self.format_report( 
                            "FailedFilterFor_"+output_record['name'], 
                            output_record['output'] ) ,file=self.configuration.report_fh)
    
        # Finally, write all the outputs, to main stream if passed, otherwise to
        # the failed output (if provided)
        for output_record in output_records:
            if output_record['filter_result'] and output_record['output'] is not None:
                write_out_seq(output_record['output'], self.configuration.output_fh, self.configuration.output_format, 
                    output_record['name'])
            elif self.configuration.failed_fh != None:
                write_out_seq(self.seqs['input'], self.configuration.failed_fh, self.configuration.input_format, 
                    output_record['name'])




