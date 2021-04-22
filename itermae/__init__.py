#!/usr/bin/env python3

import time
import statistics
import sys
import gzip
import string
import argparse
import re
import itertools

import yaml
import regex
from Bio import SeqIO
from Bio import Seq, SeqRecord

# TODO pass description to flags field, but this requires lots of warnings
# and caveats to the users that they will have to preface the right SAM tag
# headers and such!
def format_sam_record(record_id, sequence, qualities, tags,
        flag='0', reference_name='*', mapping_position='0', 
        mapping_quality='255', cigar_string='*', reference_name_of_mate='=', 
        position_of_mate='0', template_length='0' ):
    """This formats arguments into a string for outputting as a SAM format
    record. This is missing the ability to handle descriptions or tags other
    than the one indicating which group this is.

    :param record_id: the ID of the read
    :type record_id: str
    :param sequence: the nucleotide sequence to output
    :type sequence: str
    :param qualities: per-base qualities, encoded as letters (ASCII I think)
    :type qualities: str
    :param tags: any tags to add to the tags field, see SAM specification for
        the proper format of these, defaults to blank
    :type tags: str, optional
    :param flag: the bit-flag, defaults to '0'
    :type flag: str, optional
    :param reference_name: name of reference template, defaults to '*'
    :type reference_name: str, optional
    :param mapping_position: mapping position, defaults to '0'
    :type mapping_position: str, optional
    :param mapping_quality: mapping quality, defaults to '255'
    :type mapping_quality: str, optional 
    :param cigar_string: CIGAR string of mutations relative to reference, 
        defaults to '*'
    :type cigar_string: str, optional
    :param reference_name_of_mate: reference name of mate, defaults to '='
    :type reference_name_of_mate: str, optional
    :param position_of_mate: position of mate, defaults to '0'
    :type position_of_mate: str, optional
    :param template_length: length of template, defaults to '0' 
    :type template_length: str, optional

    :return: returns string of the fields, tab-separated for output as a
        SAM record
    :rtype: str
    """

    return "\t".join([ record_id, flag, reference_name, mapping_position, 
            mapping_quality, cigar_string, reference_name_of_mate, 
            position_of_mate, template_length, sequence, qualities, tags ])


def phred_letter_to_number(letter):
    """Simple function to turn a PHRED score from letter to number. That's it.

    :param letter: PHRED score letter
    :type letter: str
    :return: Returns just PHRED score (Illumina 1.8+, I believe) 
        score corresponding to the letter score
    :rtype: int I think
    """
    return ord(letter)-33


def phred_number_to_letter(score):
    """Simple function to turn a PHRED score from number to letter. That's it.

    :param score: PHRED score number
    :type score: int
    :return: Returns just PHRED score (Illumina 1.8+, I believe) 
        letter corresponding to the numeric score
    :rtype: str
    """
    return chr(score+33)


def phred_number_array_to_joined_string(score_array):
    """Turn a list of PHRED score numbers to a letter string. That's it.

    :param score_array: PHRED score array
    :type score_array: list of int
    :return: Returns a string of the PHRED scores (Illumina 1.8+, I believe) 
        converted from a numeric list to a letter string.
    :rtype: str
    """
    return str("".join([ phred_number_to_letter(i) for i in score_array]))


def read_sam_file(fh):
    """This is a minimal SAM reader, just for getting the fields I like and 
    yielding SeqRecord objects, sort of like BioPython SeqIO. Here, we are
    putting SAM tags into the description field so it should be possible to 
    pass those through, but that's not well designed yet.
    
    :param fh: file handle to read
    :type fh: file handle opened by
    :return: yields SeqRecords
    :rtype: Bio.SeqRecord.SeqRecord
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
    """Reads a text file, and yields SeqRecords where the string in the line
    is the sequence and the ID of the record.
    
    :param fh: file handle opened by 
    :type fh: file handle opened by 
    :return: yields SeqRecords
    :rtype: Bio.SeqRecord.SeqRecord
    """
    for i in fh.readlines():
        seq = i.rstrip()
        yield SeqRecord.SeqRecord( Seq.Seq(seq), id=seq, description="")


# TODO consider moving the 'which' bit to something specified in the 
# build_context, sort of like 'id' and 'description'
def write_out_seq(seq,fh,format,which):
    """This little utility just handles which of the four formats to print out,
    and for SAM appends a tag with which match this is, using the IE tag.
    
    :param seq: The SeqRecord to write
    :type seq: Bio.SeqRecord
    :param fh: file handle
    :type fh: file handle returned by
    :param format: which format to output, one of 'sam', 'txt', or something
        that Bio.SeqIO will recognize
    :type format: str
    :param which: which output this is, so for SAM this appened to a tag, but
        is ignored for the other formats
    :type which: str
    :return: nothing, it writes to a file
    :rtype: None
    """
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
    """This class is for configuring itermae, from YAML or CLI arguments.
    No arguments for initializing, it will set default values.
    Then you use the configuration methods.
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
        self.output = 'STDOUT'
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
        """Opens file-handle based on the configuration.
        Requires `input` to be set.

        :raises ValueError: Can't handle gzipped inputs on STDIN.
        """
        if self.input.upper() == 'STDIN':
            if self.gzipped:
                raise ValueError("I can't handle gzipped inputs on STDIN ! "
                    "You shouldn't see this error, it shoulda been caught in "
                    "the launcher script.") 
            else:
                self.input_fh = sys.stdin
        else:
            if self.gzipped:
                self.input_fh = gzip.open(self.input,'rt',encoding='ascii')
            else:
                self.input_fh = open(self.input,'rt')

    def open_appropriate_input_format(self):
        """Uses `input_format` and `input_fh` to set iterators
        of SeqRecords from the appropriate inputs, in `input_seqs`.
        Tries to handle all formats known, but will try with SeqIO
        in case there's one I didn't think about.
        """
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
        """This calls `open_input_fh()` to set the `input_fh` attribute,
        then calls `open_appropriate_input_format` to use this and the 
        `input_format` attribute to save an iterator of SeqRecords
        into `input_seqs`.

        Note this is inconsistent with design of the output, will pick one or
        the other ... later.
        """
        self.open_input_fh()
        self.open_appropriate_input_format()

    def open_output_fh(self,file_string):
        """Opens output file handle, which can then be written to later with 
        a format specification.

        Note this is inconsistent with design of the input, will pick one or
        the other ... later.

        :param file_string: file to wrote to, or STDOUT or STDERR
        :type file_string: str
        :return: file string for appending output
        :rtype: file handle returned by `open()`
        """
        if file_string is None:
            return None
        if file_string.upper() == 'STDOUT':
            return sys.stdout
        elif file_string.upper() == 'STDERR':
            return sys.stderr
        else:
            return open(file_string,'a')

    def close_fhs(self):
        """This is for cleaning up, and tries to close file handles at
        `input_seqs`, `ouput_fh`, `failed_fh`, `report_fh`.
        """
        for i in [ self.input_seqs, self.output_fh, self.failed_fh, self.report_fh] :
            try:
                i.close()
            except:
                pass

    def check_reserved_name(self,name,
            reserved_names=['dummyspacer','input','id','description'] ):
        """This checks if the name is one of a reserved list, and raises error
        if so. These names are reserved for these reasons:

        - `dummyspacer` is so you can pop an X into your sequence as a separator
        delimiter for later processing
        - `input` is the input group, the original one
        - `id` is the input ID, here just as `id` so it`s easy to find
        - `description` is for mapping over the FASTQ description

        :param name: name of group
        :type name: str
        :raises ValueError: raised if you're using one of the reserved names...
        """
        if name in reserved_names:
            raise ValueError("Hey, you can't name a capture group "+
                (" or ".join(reserved_names[ [(i == name) for i in reserved_names]]))+
                ", I'm using that/those! Pick a different name.")

    def config_from_file(self,file_path):
        """Tries to parse a configuration YAML file to update this configuration
        object. Pass in the file path as an argument.
        Recommend you run this config first, then config_from_args, as done in
        `bin/itermae`.
        
        :param file_path: file path to configure from, expecting it to point to
            an appropriately formatted YAML file
        :type file_path: str
        :raises ValueError: Failure to parse the supplied YAML
        :raise KeyError: You need to define a group called `pattern:` 
            inside each of the list inside of `matches:`
        :raise ValueError: Error in yaml config, you`ve repeated a group 
            marking character to match in multiple places
        :raise ValueError: Error in yaml config, the pattern and marking 
            you`ve defined are of different lengths
        :raise ValueError: Error in yaml config
        :raise KeyError: Marked roup in `marking:` field does not have
            corresponding entry in `marked_groups:`. 
        :raise ValueError: Either the supplied `filter`, `id`, `seq`, or 
            `description` expression for a match group does not look like a 
            python expression
        """

        if file_path is None:
            return
    
        try:
            with open(file_path,'r') as f:
                config = yaml.load(f,Loader=yaml.SafeLoader)
        except Exception as error:
            raise ValueError(repr(error)+" : "
                "I failed to parse the supplied YAML file at that path.")

        # Looking for verbosity instruction global, if not global, then in 'outputs'
        try:
            self.verbosity = config['verbosity']
        except:
            pass

        try:
            self.input = config['input_from']
        except:
            pass
        try:
            self.input_format = config['input_format']
        except:
            pass
        try:
            self.gzipped = config['input_gzipped']
        except:
            pass
        try:
            self.output = config['output_to']
        except:
            pass
        try:
            self.output_format = config['output_format']
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
            try:
                assert each['pattern']
            except Exception as error:
                raise KeyError("You need to define a group called 'pattern:' "
                    "inside each of the list (denoted by '-'s) inside of "
                    "'matches:' - what is the sequence pattern to match?")
            if self.verbosity >= 1:
                print("    Taking '"+each['use']+"'. \n", end="",file=sys.stderr)
            if len(re.sub(r'(.)\1+',r'\1',each['marking'])) > len(set(each['marking'])):
                raise ValueError("Error in reading yaml config! It looks like "
                    "you've repeated a group marking character to match in "
                    "multiple places. I do not support that, "
                    "use a different character.")
            if len(each['pattern']) != len(each['marking']):
                raise ValueError("Error in reading yaml config! "
                    "The pattern and marking you've defined are of "
                    "different lengths. I need them to be the same length.")
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

                try:
                    assert each['marked_groups'][mark]
                except Exception as error:
                    raise KeyError("You've marked a group in the 'marking:' "
                        "field but have not supplied a corresponding entry in "
                        "'marked_groups:', hence "+repr(error)+".")

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

        output_list = config['output_list'] # I do need some outputs, or fail
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

            try:
                self.outputs_array.append( {
                        'name':each['name'], # These are on oneline for error
                            # readability about which one is the problem
                        'filter':[ each['filter'], compile(each['filter'],'<string>','eval',optimize=2) ],
                        'id':[ each['id'], compile(each['id'],'<string>','eval',optimize=2) ],
                        'seq':[ each['seq'], compile(each['seq'],'<string>','eval',optimize=2) ],
                        'description':[ each['description'], compile(each['description'],'<string>','eval',optimize=2) ]
                    })
            except Exception as error:
                raise ValueError(repr(error)+" : "
                    "Either the supplied 'filter', 'id', 'seq', "
                    "or 'description' expression for a match group does "
                    "not look like a python expression - are all "
                    "non-group-name parts in quotes? Are group-names and "
                    "other parts connected with + signs?")


        try:
            self.failed = config['output_failed']
        except:
            pass
        try:
            self.report = config['output_report']
        except:
            pass


    def config_from_args(self,args_copy):
        """Make configuration object from arguments provided. Should be the 
        same as the config_from_yaml output, if supplied the same.

        :param args_copy: pass in the `argparse` args object after collecting
            the startup command line arguments
        :type args_copy: argparse object, I think
        :raise ValueError: I failed to build the regular expression for a match
        :raise ValueError: The output IDs, seqs, descriptions, and filters are 
            of unequal sizes, make them equal or only define one of each 
        :raise ValueError: Either the supplied `filter`, `id`, `seq`, or 
            `description` expression for a match group does not look like a 
            python expression
        """
    
        if args_copy.verbose:
            self.verbosity = args_copy.verbose

        for each in args_copy.match:
            try:
                for capture_name in re.findall('<(.*?)>',each):
                    self.check_reserved_name(capture_name)
                try:
                    (input_string, regex_string) = re.split(r"\s>\s",each.strip())
                except:
                    input_string = 'input' # default to just use raw input
                    regex_string = each.strip()
                compiled_regex = regex.compile(
                    regex_string.strip(), # We use this regex
                    regex.BESTMATCH # And we use the BESTMATCH strategy, I think
                    )
                self.matches_array.append( {'input':input_string.strip(), 'regex':compiled_regex} )
            except Exception as error:
                raise ValueError(repr(error)+" : "
                    "I failed to build the regular expression from the "
                    "command-line argument supplied.")

        # Adding in defaults for outputs. Can't do that with argparse, I think,
        # because this needs to be appending. First add in defaults, but
        # absolutely first need an output_seq to be defined for it to try this:
        if args_copy.output_seq:
            if not args_copy.output_id:
                args_copy.output_id = ['id']
            if not args_copy.output_filter:
                args_copy.output_filter = ['True']
            if not args_copy.output_description:
                args_copy.output_description = ['description']
            # Then normalize the length 1 to max length
            maximum_number_of_outputs = max( [len(args_copy.output_id), 
                len(args_copy.output_seq), len(args_copy.output_filter),
                len(args_copy.output_description)] )
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
                raise ValueError("The output IDs, seqs, descriptions, and "
                    "filters are of unequal sizes. Make them equal, or only "
                    "define one each and it will be reused across all."+
                    repr(( len(args_copy.output_id), len(args_copy.output_seq),
                        len(args_copy.output_filter), len(args_copy.output_description) )) )
    
            i = 0
            for idz, seqz, filterz, description in zip(args_copy.output_id, args_copy.output_seq, args_copy.output_filter, args_copy.output_description) :
                this_name = 'untitled_output_'+str(i)
                i += 1
                try:
                    self.outputs_array.append( {   
                            'name': this_name,
                            'filter': [ filterz, compile(filterz,'<string>','eval',optimize=2) ],
                            'id': [ idz, compile(idz,'<string>','eval',optimize=2) ],
                            'seq': [ seqz, compile(seqz,'<string>','eval',optimize=2) ] ,
                            'description':[ description, compile(description,'<string>','eval',optimize=2) ]
                        })
                except Exception as error:
                    raise ValueError(repr(error)+" : "
                        "Either the supplied 'filter', 'id', 'seq', "
                        "or 'description' expression for a match group does "
                        "not look like a python expression - are all "
                        "non-group-name parts in quotes? Are group-names and "
                        "other parts connected with + signs?")
        
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
            return_string += '\n        - input: '+each['input']
            return_string += '\n          regex: '+str(each['regex'])
        return_string += '\n    writing these outputs:'
        for each in self.outputs_array:
            return_string += '\n        - id: '+str(each['id'][0])
            return_string += '\n          description: '+str(each['description'][0])
            return_string += '\n          seq: '+str(each['seq'][0])
            return_string += '\n          filter: '+str(each['filter'][0])
        return return_string

    def reader(self):
        """This reads inputs, calls the `chop` method on each one, and sorts 
        it off to outputs. So this is called by the main function, and is 
        mostly about handling the I/O and handing it to the `chop` function.
        Thus, this depends on the `Configuration` class being properly 
        configured with all the appropriate values.
        """
    
        # Input
        self.get_input_seqs()
    
        # Outputs - passed records, failed records, report file
        self.output_fh = self.open_output_fh(self.output)
        self.report_fh = self.open_output_fh(self.report)
        self.failed_fh = self.open_output_fh(self.failed)
    
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
    

class MatchScores:
    """This is a little class just to hold the three scores under attributes,
    such that they're easier to type for writing filters. Also, it flattens
    them for debug report printing.

    :param substitions: number to store under `.substitions` attribute
    :type substitions: int
    :param insertions: number to store under `.insertions` attribute
    :type insertions: int
    :param deletions: number to store under `.deletions` attribute
    :type deletions: int
    """
    def __init__(self, substitutions, insertions, deletions):
        self.substitutions = substitutions
        self.insertions = insertions
        self.deletions = deletions
    def flatten(self):
        """Flatten this object for printing debug reports.
    
        :return: string in form substitutions_insertions_deletions
        :rtype: str
        """
        return str(self.substitutions)+"_"+str(self.insertions)+"_"+\
            str(self.deletions)


class GroupStats:
    """Object for conveniently holding parameters from the match, so that
    they're easier to type for filters/output specification, and to flatten
    for debug printing.

    :param start: number to store under `.start` attribute
    :type start: int
    :param end: number to store under `.end` attribute
    :type end: int
    :param length: number to store under `.length` attribute
    :type length: int
    :param quality: list of numbers to store under `.quality` attribute
    :type quality: list of int
    :param quality_string: string of the quality array under PHRED encodings
    :type quality_string: string
    """

    def __init__(self, start, end, seq, quality):
        self.start = start 
        self.end = end 
        self.length = self.end - self.start
        self.seq = seq
        self.quality = quality
        self.quality_string = phred_number_array_to_joined_string(quality)

    def flatten(self):
        """Flatten this object for printing debug reports, but just for
        the start, end, length attributes. Not quality.
    
        :return: string in form start_end_length
        :rtype: str
        """
        return str(self.start)+"_"+str(self.end)+"_"+str(self.length)

    def __eq__(self,other):
        """Attention! This is a hack to allow for using the group's name
        (ie 'barcode') instead of accessing the '.seq' method.
        """
        return str(self.seq.seq) == other


class SeqHolder: 
    """This is the main holder of sequences, and has methods for doing matching,
    building contexts, filtering, etcetra. Basically there is one of these
    initialized per input, then each operation is done with this object, then
    it generates the appropriate outputs and `chop` actually writes them.
    Used in `chop`.

    The `.seqs` attribute holds the sequences accessed by the matching,
    initialized with the `input_record` SeqRecord and a `dummyspacer` for
    output formatting with a separator.
    
    :param input_record: an input SeqRecord object
    :type input_record: Bio.SeqRecord.SeqRecord
    :param configuration: the whole program's Configuration object, with
        appropriate file-handles opened up and defaults set
    :type configuration: itermae.Configuration
#    :raises [ErrorType]: [ErrorDescription]
#    :return: [ReturnDescription]
#    :rtype: [ReturnType]
    """
    def __init__(self, input_record, configuration):
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
        """This applies the given match to the `SeqHolder` object, and saves 
        how it did internally.

        :param match_id: what name should we call this match? This is useful
            for debugging reports and filtering only.
        :type match_id: str
        :param input_group: which input group to use, by name of the group
        :type input_group: str
        :param regex: the regular expression to apply, complete with named 
            groups to save for subsequent match operations
        :type regex: regex compiled regular expression object
        :return: self, this is just done so it can exit early if no valid input
        :rtype: itermae.SeqHolder
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
                        seq=self.seqs[match_name],
                        quality=self.seqs[match_name].letter_annotations['phred_quality']
                        )

        except:
            self.match_scores[match_id] = MatchScores(None,None,None)

    def build_context(self):
        """This unpacks group match stats/scores into an environment that
        the filter can then use to ... well ... filter. 
        """

        # This is context for the filters, so is operating more as values,
        # as opposed to the context_seq which is operating with SeqRecords
        self.context_filter = { **self.group_stats , **self.match_scores }

        # Then unpack the sequences as a context for building the output 
        # sequences, this is different so that the qualities get stuck with
        # the bases of the groups
        self.context_seq = { **self.seqs }

        # Then one for the IDs, so we're setting the input ID as 'id', and then
        # each group name just refers to the sequence. And I finally put seq 
        # qualities in the ID. We do make 'description' available if needed
        self.context_id = { 
            'id': self.seqs['input'].id , 
            'description': self.seqs['input'].description , 
            **{ i: str(self.seqs[i].seq) for i in self.seqs } ,
            **{ i+'_quality': self.group_stats[i].quality_string 
                    for i in self.group_stats } }

    def evaluate_filter_of_output(self,output_dict):
        """This tests a user-defined filter on the 'seq_holder' object.
        This has already been `compile`'d, and here we just attempt to
        evaluate these to `True`, where `True` is passing the filter.
        Exceptions are blocked by using `try`/`except` so that it can fail on 
        a single match and move onto the next match/read.

        :param output_dict: a dictionary of outputs to form, as generated from
            the configuration initialization
        :type output_dict: dict
        :return: `True` if the filter passed and the output should be generated
        :rtype: bool
        """

        try:
            filter_result = eval(output_dict['filter'][1],globals(),self.context_filter)
            if self.configuration.verbosity >= 3:
                print("\n["+str(time.time())+"] : This read "+
                    self.seqs['input'].id+" successfully evaluated the filter "+
                    str(output_dict['filter'][0])+" as "+str(filter_result),
                    file=sys.stderr)
            return filter_result
        except Exception as error:
            if self.configuration.verbosity >= 3:
                print("\n["+str(time.time())+"] : This read "+
                    self.seqs['input'].id+" failed to evaluate the filter "+
                    str(output_dict['filter'][0])+". "+
                    repr(error),file=sys.stderr)
            return False

    def build_output(self,output_dict):
        """Builds the output from the `SeqHolder` object according to the 
        outputs in `output_dict`.

        :param output_dict: a dictionary of outputs to form, as generated from
            the configuration initialization
        :type output_dict: dict
        :return: the successfully built SeqRecord, or None if it fails
        :rtype: Bio.SeqRecord.SeqRecord or None
        """

        try:
            output_seq = eval(output_dict['seq'][1],globals(),self.context_seq)
            out_seq = SeqRecord.SeqRecord(
                seq = Seq.Seq(str(output_seq.seq)) ,
                id = str(eval(output_dict['id'][1],globals(),self.context_id)) ,
                description = str(eval(output_dict['description'][1],globals(),self.context_id)) ,
                letter_annotations = {'phred_quality':output_seq.letter_annotations['phred_quality']}
            )
            if self.configuration.verbosity >= 3:
                print("\n["+str(time.time())+"] : This read "+
                    self.seqs['input'].id+" successfully built the output of "+
                    "id: '"+str(output_dict['id'][0])+"', and "+
                    "description: '"+str(output_dict['description'][0])+"', and "+
                    "seq: '"+str(output_dict['seq'][0])+"'." ,file=sys.stderr)
            return out_seq
        except Exception as error:
            if self.configuration.verbosity >= 3:
                print("\n["+str(time.time())+"] : This read "+
                    self.seqs['input'].id+" failed to build the output of "+
                    "id: '"+str(output_dict['id'][0])+"', and "+
                    "description: '"+str(output_dict['description'][0])+"', and "+
                    "seq: '"+str(output_dict['seq'][0])+"'. "+
                    repr(error) ,file=sys.stderr)
            return None

    def format_report(self,label,output_seq):
        """Formats a standard report line for the debug reporting function.

        :param label: what type of report line this is, so a string describing
            how it went - passed? Failed?
        :type label: str
        :param label: the attempt at generating an output SeqRecord, so either
            one that was formed or None
        :type label: Bio.SeqRecord.SeqRecord or None
        :return: the string for the report
        :rtype: str
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
        """This executes the intended purpose of the `SeqRecord` object, and is
        called once. It uses the configured object to apply each match
        operation as best it can with the sequences it is given or can generate,
        then writes the outputs in the specified formats to specified places
        as configured.
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
                phred_number_array_to_joined_string(self.seqs['input'].letter_annotations['phred_quality']),
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
                if self.configuration.verbosity >= 3:
                    print("\n["+str(time.time())+"] : wrote out output '"+
                        output_record['name']+"' for this input",
                        file=sys.stderr)
            elif self.configuration.failed_fh != None:
                write_out_seq(self.seqs['input'], self.configuration.failed_fh, self.configuration.input_format, 
                    output_record['name'])
                if self.configuration.verbosity >= 3:
                    print("\n["+str(time.time())+"] : output "+
                        output_record['name']+" failed, written to fail file\n",
                        file=sys.stderr)
