#!/usr/bin/env python3

import pytest

import itermae

# Required for testing apply_operation, as regex is generated in the command
# line program
import regex
# Required for testing SeqHolder etc. without the file IO of main script
from Bio import SeqIO
# For doing combinations of the full-file parameters and formats
import itertools
# Required for full-file input/output testing
import subprocess
# For checking
import re

#### Ye Tests

# Same test for coming from a YML file
#@pytest.fixture
def configuration_yaml():
    configuration = itermae.Configuration()
    configuration.config_from_file("itermae/data/example_schema.yml")
    return configuration

#@pytest.fixture
def benchmark_config():
    return {
        'verbosity': 1,
        'input_from': 'STDIN',
        'input_format': 'fastq',
        'input_gzipped': False,
        'matches': [
            {   'use': 'input',
                'pattern': 'NGTCCACGAGGTCTCTNCGTACGCTG',
                'marking': 'ABBBBBBBBBBBBBBBCDDDDDDDDD',
                'marked_groups': {
                    'A': { 'name': 'sampleIndex',
                        'repeat': 5 },
                    'B': { 'name': 'prefix',
                        'allowed_errors': 2 },
                    'C': { 'name': 'barcode',
                        'repeat_min': 18,
                        'repeat_max': 22 },
                    'D': { 'allowed_insertions': 1,
                        'allowed_deletions': 2,
                        'allowed_substititions': 2 }
                }
            },
            {   'use': 'barcode',
                'pattern': 'N',
                'marking': 'z',
                'marked_groups': {
                    'z': { 'name': 'first_five_barcode',
                        'repeat': 5 },
                }
            }
        ],
        'output_list': [
            {   'name': 'short_barcodes',
                'seq': 'barcode',
                'filter': 'True',
                'description': 'description' },
            {   'name': 'sampleIndex',
                'seq': 'sampleIndex',
                'filter': 'sampleIndex.length == 5',
                'description': 'description+" this is just the sampleIndex"' },
            {   'name': 'usualbarcode',
                'seq': 'barcode',
                'id': 'id',
                'description': 'description+" sample="+sampleIndex' },
            {   'name': 'other_barcodes',
                'seq': 'barcode',
                'filter': 'True',
                'description': 'description+" other_barcodes"' }
        ],
        'output_to': 'STDOUT',
        'output_format': 'fasta',
        'output_failed': 'failed.fastq',
        'output_report': 'report.csv',
    }


#def ztest_configuration_args(configuration_yaml):
def ztest_configuration_args(conf,benchmark):

    assert benchmark['verbosity']     == getattr(conf,'verbosity')
    assert benchmark['input_from']    == getattr(conf,'input')
    assert benchmark['input_format']  == getattr(conf,'input_format')
    assert benchmark['input_gzipped'] == getattr(conf,'gzipped')
    assert benchmark['output_to']     == getattr(conf,'output')
    assert benchmark['output_format'] == getattr(conf,'output_format')
    assert benchmark['output_failed'] == getattr(conf,'failed')
    assert benchmark['output_report'] == getattr(conf,'report')

    for i in range(len(benchmark['output_list'])):
        try: # Names should be equal
            assert( benchmark['output_list'][i]['name'] == 
                getattr(conf,'outputs_array')[i]['name'] )
        except: # or they should be set in conf to something
            assert getattr(conf,'outputs_array')[i]['name'] 
        # Output sequence is required, we make no attempt to check 'compile'
        assert ( benchmark['output_list'][i]['seq'] == 
                getattr(conf,'outputs_array')[i]['seq'][0] )
        try: # The filter specification should be equal,
            assert ( benchmark['output_list'][i]['filter'] == 
                getattr(conf,'outputs_array')[i]['filter'][0] )
        except: # or should be not False. True == does not work ... mysterious
            assert ( False != getattr(conf,'outputs_array')[i]['filter'][0] )
        try: # The description specification should be equal,
            assert ( benchmark['output_list'][i]['description'] == 
                getattr(conf,'outputs_array')[i]['description'][0] )
        except: # or set to 'description'
            assert ( 'description' == getattr(conf,'outputs_array')[i]['description'][0] )

    for i in range(len(benchmark['matches'])):
        # 'use' is turned into 'input'
        assert( benchmark['matches'][i]['use'] == 
                getattr(conf,'matches_array')[i]['input'] )
        for j in benchmark['matches'][i]['marked_groups']:
            try: # for all matches, what's the benchmark group name?
                bench_name = benchmark['matches'][i]['marked_groups'][j]['name']
            except: # or is generic untitled group name
                bench_name = 'untitled_group\d+'
            group_pattern = re.compile("\(\?<"+bench_name+">")
            # do we find that group name?
            assert( group_pattern.search(str(getattr(conf,'matches_array')[i]['regex'] )))
            try: # next we move onto pulling in any repeat specs
                repeat_min = str(benchmark['matches'][i]['marked_groups'][j]['repeat_min'])
            except: # ... laboriously to avoid catching errors
                repeat_min = None
            try:
                repeat_max = str(benchmark['matches'][i]['marked_groups'][j]['repeat_max'])
            except:
                repeat_max = None
            try:
                repeat = str(benchmark['matches'][i]['marked_groups'][j]['repeat'])
            except:
                repeat = None
            if repeat_min and repeat_max: # then depending on what we have available, we search
                assert ( re.search("\{"+repeat_min+","+repeat_max+"}",
                        str(getattr(conf,'matches_array')[i]['regex'] )))
            if repeat_min and repeat and ( repeat_max == None ): 
                assert ( re.search("\{"+repeat_min+","+repeat+"}",
                        str(getattr(conf,'matches_array')[i]['regex'] )))
            if ( repeat_min == None ) and repeat and repeat_max: 
                assert ( re.search("\{"+repeat+","+repeat_max+"}",
                        str(getattr(conf,'matches_array')[i]['regex'] )))
            if ( repeat_min == None ) and repeat and ( repeat_max == None ): 
                assert ( re.search("\{"+repeat+","+repeat+"}",
                        str(getattr(conf,'matches_array')[i]['regex'] )))

ztest_configuration_args(configuration_yaml(),benchmark_config())

def ztest_configuration_args(benchmark):
    class A: pass # Here I just need a burner object that accepts new attributes
    args_test = A()
    args_test.input = 'STDIN'
    args_test.input_format = 'fastq'
    args_test.gzipped = False
    args_test.match = ["input > (?P<sampleIndex>[ATCGN]{5,5})(?P<upPrime>GTCCTCGAGGTCTCT){e<=1}(?P<barcode>[ATCGN]{18,22})(?P<downPrime>CGTACGCTG){e<=1}"]
    args_test.output_seq = ["barcode"]
    args_test.output_id = ["id"]
    args_test.output_description = None
    args_test.output_filter = None
    args_test.verbose = 1
    args_test.output = "STDOUT"
    args_test.output_format = "fasta"
    args_test.failed = 'failed.fastq'
    args_test.report = 'report.csv'
    conf = itermae.Configuration()
    conf.config_from_args(args_copy=args_test)

    assert benchmark['verbosity']     == getattr(conf,'verbosity')
    assert benchmark['input_from']    == getattr(conf,'input')
    assert benchmark['input_from']    == getattr(conf,'input')
    assert benchmark['input_format']  == getattr(conf,'input_format')
    assert benchmark['input_gzipped'] == getattr(conf,'gzipped')
    assert benchmark['output_to']     == getattr(conf,'output')
    assert benchmark['output_format'] == getattr(conf,'output_format')
    assert benchmark['output_format'] == getattr(conf,'output_format')
    assert benchmark['output_failed'] == getattr(conf,'failed')
    assert benchmark['output_report'] == getattr(conf,'report')

ztest_configuration_args(benchmark_config())


# Test MatchScores class
#@pytest.fixture
def matchscore():
    return itermae.MatchScores(1,2,3)
def test_matchscore_subs(matchscore):
    assert matchscore.substitutions == 1
def test_matchscore_ins(matchscore):
    assert matchscore.insertions == 2
def test_matchscore_dels(matchscore):
    assert matchscore.deletions == 3
def test_matchscore_flatten(matchscore):
    assert matchscore.flatten() == "1_2_3"

matchscore()

# Test GroupStats class
#@pytest.fixture
def groupstats():
    return itermae.GroupStats(5,15,[36]*10)
def test_groupstats_start(groupstats):
    assert groupstats.start == 5
def test_groupstats_end(groupstats):
    assert groupstats.end == 15
def test_groupstats_length(groupstats):
    assert groupstats.length == 10
def test_groupstats_quality(groupstats):
    assert groupstats.quality == [36]*10
def test_groupstats_flatten(groupstats):
    assert groupstats.flatten() == "5_15_10"

groupstats()

# Setup inputs
# I'm not testing BioPython SeqIO, I assume that's good.
# This instead uses that to return a list of the records in the file.
#@pytest.fixture
def fastqfile():
    return SeqIO.parse("itermae/data/barseq.fastq","fastq")

## SeqHolder Tests


# Test SeqHolder dummyspacer is right
def test_seqholder_dummy(fastqfile):
    for i in fastqfile:
        seqholder = itermae.SeqHolder(i,verbosity=3)
        # Is the sequence X?
        assert seqholder.seqs['dummyspacer'].seq == 'X'
        # Is the number we just put there 40?
        assert seqholder.seqs['dummyspacer'].letter_annotations['phred_quality'] == [40] 

ztest_configuration_args(configuration_yaml(),benchmark_config())

# Test that SeqHolder can apply_operation, then since we're there testing
# that it finds the right groups for each seq, and passes or fails filters 
# appropriately.
def test_seqholder_match_filter(fastqfile):
    for seq, pos_pass, qual_pass, seq_pass, sequences_found, \
        seq_targets, report_targets \
        in zip(fastqfile,
            [ i == 1 for i in [1,1,0,1,0,0,1,1,1] ],
            [ i == 1 for i in [0,0,0,1,0,0,1,1,0] ],
            [ i == 1 for i in [1,0,0,0,0,0,0,1,0] ],
            [   set(['dummyspacer','input','sample','fixed1','rest','tag','strain']),
                set(['dummyspacer','input','sample','fixed1','rest','tag','strain','fixed2','UMItail']),
                set(['dummyspacer','input']),
                set(['dummyspacer','input','sample','fixed1','rest','tag','strain','fixed2','UMItail']),
                set(['dummyspacer','input']),
                set(['dummyspacer','input']),
                set(['dummyspacer','input','sample','fixed1','rest','tag','strain','fixed2','UMItail']),
                set(['dummyspacer','input','sample','fixed1','rest','tag','strain','fixed2','UMItail']),
                set(['dummyspacer','input','sample','fixed1','rest','tag','strain','fixed2','UMItail']),
            ],
            [   ('ExampleToyReads:1:exampleFlowCell:1:10000:10000:10000_TTCAC','TCAGTCGTAGCAGTTCGATG'),
                ('ExampleToyReads:1:exampleFlowCell:1:10000:10001:10001_GCTTC', 'TGGCAGACACACGCTACA'),
                (None,None),
                ('ExampleToyReads:1:exampleFlowCell:1:10000:10003:10003_CTACT', 'GATGCACTGCGTTCCATGTT'),
                (None,None),
                (None,None),
                ('ExampleToyReads:1:exampleFlowCell:1:10000:10006:10006_TCGGC', 'ATTCTGAGCGGTGCCATAGT'),
                ('ExampleToyReads:1:exampleFlowCell:1:10000:10007:10007_AGGAG', 'ATAAGTTAGACAGGTCAGC'),
                ('ExampleToyReads:1:exampleFlowCell:1:10000:10008:10008_ACGTA', 'CACACGCACGAATTTGCATA')
                ],
            [   '"pass","ExampleToyReads:1:exampleFlowCell:1:10000:10000:10000","TTCACGTCCTCGAGGTCTCTTCAGTCGTAGCAGTTCGATGCGTACGCTACAGGTCGACGGTAAGAGAGGGATGTG","TCAGTCGTAGCAGTTCGATG","FilterResults","sample_0_5_5-fixed1_5_17_12-rest_17_75_58-tag_0_3_3-strain_3_23_20"',
                '"pass","ExampleToyReads:1:exampleFlowCell:1:10000:10001:10001","GCTTCGTCCTCGAGGTCTATTGGCAGACACACGCTACACGTACGCTGCAGGTCGAGGGCACGCGAGAGATGTGTG","TGGCAGACACACGCTACA","FilterResults","sample_0_5_5-fixed1_5_17_12-rest_17_75_58-tag_0_3_3-strain_3_21_18-fixed2_21_36_15-UMItail_36_53_17"',
                '"fail_to_form","ExampleToyReads:1:exampleFlowCell:1:10000:10002:10002","CCCGGCGTTCGGGGAAGGACGTCAATAGTCACACAGTCCTTGACGGTATAATAACCACCATCATGGCGACCATCC","TGGCAGACACACGCTACA","FilterResults",""',
                '"pass","ExampleToyReads:1:exampleFlowCell:1:10000:10003:10003","CTACTGTCCACGAGGTCTCTGATGCACTGCGTTCCATGTTCGTACGCTGCAGGTCGACGGAAGGAGCGCGATGTG","GATGCACTGCGTTCCATGTT","FilterResults","sample_0_5_5-fixed1_5_17_12-rest_17_75_58-tag_0_3_3-strain_3_23_20-fixed2_23_38_15-UMItail_38_55_17"',
                '"fail_to_form","ExampleToyReads:1:exampleFlowCell:1:10000:10004:10004","AAATTAGGGTCAACGCTACCTGTAGGAAGTGTCCGCATAAAGTGCACCGCATGGAAATGAAGACGGCCATTAGCT","GATGCACTGCGTTCCATGTT","FilterResults",""',
                '"fail_to_form","ExampleToyReads:1:exampleFlowCell:1:10000:10005:10005","CCCGGCGTTCGGGGAAGGACGTCAATAGTCACACAGTCCTTGACGGTATAATAACCACCATCATGGCGACCATCC","GATGCACTGCGTTCCATGTT","FilterResults",""',
                '"pass","ExampleToyReads:1:exampleFlowCell:1:10000:10006:10006","TCGGCGTCCTCGAGGTCTCTATTCTGAGCGGTGCCATAGTCGTACGCTGCAGGTCGACCGAAGGTGGGAGATGTG","ATTCTGAGCGGTGCCATAGT","FilterResults","sample_0_5_5-fixed1_5_17_12-rest_17_75_58-tag_0_3_3-strain_3_23_20-fixed2_23_38_15-UMItail_38_55_17"',
                '"pass","ExampleToyReads:1:exampleFlowCell:1:10000:10007:10007","AGGAGGTCCTCGAGGTCTCTATAAGTTAGACAGGTCAGCCGTACGCTGCAGGTCGACAGCTGGCGCGCGATGTGA","ATAAGTTAGACAGGTCAGC","FilterResults","sample_0_5_5-fixed1_5_17_12-rest_17_75_58-tag_0_3_3-strain_3_22_19-fixed2_22_37_15-UMItail_37_54_17"',
                '"pass","ExampleToyReads:1:exampleFlowCell:1:10000:10008:10008","ACGTAGTCCACGAGGTCTCTCACACGCACGAATTTGCATACGTACGCTGCAGGTCGACTGGAAGGGCGGGATGTG","CACACGCACGAATTTGCATA","FilterResults","sample_0_5_5-fixed1_5_17_12-rest_17_75_58-tag_0_3_3-strain_3_23_20-fixed2_23_38_15-UMItail_38_55_17"'
            ]
            ):
        # Read in the sequence to the holder
        seqholder = itermae.SeqHolder(seq,verbosity=0)
        # Apply operations
        seqholder.apply_operation('a','input',
            regex.compile("(?P<sample>[ATCG]{5})(?P<fixed1>GTCCACGAGGTC){e<=2}(?P<rest>TCT.*){e<=1}",
                regex.BESTMATCH) )
        seqholder.apply_operation('b','rest',
            regex.compile("(?P<tag>TCT){e<=1}(?P<strain>[ATCG]{10,26})(CGTACGCTGC){e<=2}",
                regex.BESTMATCH) )
        seqholder.apply_operation('c','rest',
            regex.compile("(?P<fixed2>CGTACGCTGCAGGTC)(?<UMItail>GAC[ATCG]G[ATCG]A[ATCG]G[ATCG]G[ATCG]G[ATCG]GAT){s<=2}",
                regex.BESTMATCH) )
        # Are the right sequences found/matched for each read?
        assert set(seqholder.seqs.keys()) == sequences_found 
        # Does it pass a position filter?
        seqholder.build_context()
        first_filter = 'sample.length == 5 and rest.start >= 15'
        first_id = "id+'_'+sample"
        first_seq = "strain"
        first_desc = "description"
        filter_result = seqholder.evaluate_filter_of_output(
                {   'name':'test',
                    'filter': [ first_filter,
                        compile(first_filter,'<string>','eval',optimize=2) ] ,
                    'id': [ first_id, 
                        compile(first_id,'<string>','eval',optimize=2) ] ,
                    'seq': [ first_seq,
                        compile(first_seq,'<string>','eval',optimize=2) ],
                    'description': [ first_desc,
                        compile(first_desc,'<string>','eval',optimize=2) ]
                }) 
        assert pos_pass == filter_result
        # Does it pass a quality filter, with statistics?
        second_filter = 'statistics.mean(fixed1.quality) >= 33.5'
        filter_result = seqholder.evaluate_filter_of_output(
                {   'name':'test',
                    'filter': [ second_filter, 
                        compile(second_filter,'<string>','eval',optimize=2) ],
                    'id': [ first_id, 
                        compile(first_id,'<string>','eval',optimize=2) ] ,
                    'seq': [ first_seq,
                        compile(first_seq,'<string>','eval',optimize=2) ],
                    'description': [ first_desc,
                        compile(first_desc,'<string>','eval',optimize=2) ]
                }) 
        assert qual_pass == filter_result
        # Does it pass a specific sequence filter?
        third_filter = 'sample.seq == "TTCAC" or sample.seq == "AGGAG"'
        filter_result = seqholder.evaluate_filter_of_output(
                {   'name':'test',
                    'filter': [ third_filter,
                        compile(third_filter,'<string>','eval',optimize=2) ] ,
                    'id': [ first_id, 
                        compile(first_id,'<string>','eval',optimize=2) ] ,
                    'seq': [ first_seq,
                        compile(first_seq,'<string>','eval',optimize=2) ],
                    'description': [ first_desc,
                        compile(first_desc,'<string>','eval',optimize=2) ]
                }) 
        assert seq_pass == filter_result
        # Then test outputs
        built_output = seqholder.build_output(
                {   'name':'test',
                    'filter': [ 'True', compile('True','<string>','eval',optimize=2) ] ,
                    'id': [ first_id, 
                        compile(first_id,'<string>','eval',optimize=2) ] ,
                    'seq': [ first_seq,
                        compile(first_seq,'<string>','eval',optimize=2) ],
                    'description': [ first_desc,
                        compile(first_desc,'<string>','eval',optimize=2) ]
                }) 
        # Are the right outputs constructed?
        if built_output is None:
            assert seq_targets == ( None, None)
        else:
            assert seq_targets == ( built_output.id, built_output.seq ) 

## Full Tests

# Using the YML configuration file, but only with one format of output tested
def test_full_combinations_yml():
    results = subprocess.run(
        'itermae --config itermae/data/test_schema.yml',
        shell=True,capture_output=True,encoding='utf-8')
    filename = 'itermae/data/test_outputs/barseq_ymltest.sam'
#    with open(filename,'w') as f:
#        f.write(results.stdout)
    with open(filename,'r') as f:
        expected_file = f.readlines()
    for i,j in zip(results.stdout.split('\n'),expected_file):
        assert str(i) == str(j.rstrip('\n'))

# Operations for argument-specified options
one_operation_string = (
    '-m "input > (?P<sampleIndex>[ATCGN]{5,5})(?P<upPrime>GTCCTCGAGGTCTCT){e<=1}(?P<barcode>[ATCGN]{18,22})(?P<downPrime>CGTACGCTG){e<=1}" '+
    '-os "barcode" -oi "id+\\"_\\"+sampleIndex" '
    )
two_operation_string = (
    '-m "input > (?P<sampleIndex>[ATCGN]{5,5})(?P<rest>(GTCCTCGAGGTCTCT){e<=1}[ATCGN]*)" '+
    '-m "rest  > (?P<upPrime>GTCCTCGAGGTCTCT){e<=1}(?P<barcode>[ATCGN]{18,22})(?P<downPrime>CGTACGCTG){e<=1}" '+
    '-os "barcode" -oi "id+\\"_\\"+sampleIndex" '+
    '-os "upPrime+barcode+downPrime" -oi "id+\\"_withFixedFlanking_\\"+sampleIndex" '
    )

# Each operation applied to shortread FASTQ, non-gzipped
def test_full_combinations():
    operations_list = [one_operation_string, two_operation_string]
    for input_format, which_ops, output_format in itertools.product(
            ['fastq','fasta','sam','txt'],
            [0,1],
            ['fastq','fasta','sam','txt'],
        ):
        results = subprocess.run(
            'cat itermae/data/barseq.'+input_format+' | '+
            'itermae '+
            '--input-format '+input_format+' '+
            operations_list[which_ops]+
            '--verbose --output-format '+output_format ,
            shell=True,capture_output=True,encoding='utf-8')
        filename = 'itermae/data/test_outputs/barseq_combinations_'
        filename += 'input_'+input_format+'_'
        filename += str(which_ops+1)+'_ops'
        filename += '.'+output_format
#        with open(filename,'w') as f:
#            f.write(results.stdout)
        with open(filename,'r') as f:
            expected_file = f.readlines()
        for i,j in zip(results.stdout.split('\n'),expected_file):
            assert str(i) == str(j.rstrip('\n'))


# Each operation applied to shortread FASTQ, gzipped
def test_full_combinations_gzipped():
    operations_list = [one_operation_string, two_operation_string]
    for input_format, which_ops, output_format in itertools.product(
            ['fastq','fasta','sam','txt'],
            [0,1],
            ['fastq','fasta','sam','txt'],
        ):
        results = subprocess.run(
            'itermae '+
            '--input itermae/data/barseq.'+input_format+'.gz --gzipped '+
            '--input-format '+input_format+' '+
            operations_list[which_ops]+
            '--verbose --output-format '+output_format ,
            shell=True,capture_output=True,encoding='utf-8')
        filename = 'itermae/data/test_outputs/barseq_combinations_'
        filename += 'input_'+input_format+'_'
        filename += str(which_ops+1)+'_ops'
        filename += '.'+output_format
#        with open(filename,'w') as f:
#            f.write(results.stdout)
        with open(filename,'r') as f:
            expected_file = f.readlines()
        for i,j in zip(results.stdout.split('\n'),expected_file):
            assert str(i) == str(j.rstrip('\n'))
