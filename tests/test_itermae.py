#!/usr/bin/env python3

import pytest

import itermae

# Required for testing apply_operation, as regex is generated in the command
# line program
import regex
# Required for testing SeqHolder etc. without the file IO of main script
from Bio import SeqIO, Seq, SeqRecord
# For doing combinations of the full-file parameters and formats
import itertools
# Required for full-file input/output testing
import subprocess
# For checking
import re

#### Ye Tests

# Same test for coming from a YML file
@pytest.fixture
def configuration_yaml():
    configuration = itermae.Configuration()
    configuration.config_from_file("itermae/data/example_schema.yml")
    return configuration

@pytest.fixture
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


def test_configuration_yaml(configuration_yaml,benchmark_config):

    assert benchmark_config['verbosity']     == getattr(configuration_yaml,'verbosity')
    assert benchmark_config['input_from']    == getattr(configuration_yaml,'input')
    assert benchmark_config['input_format']  == getattr(configuration_yaml,'input_format')
    assert benchmark_config['input_gzipped'] == getattr(configuration_yaml,'gzipped')
    assert benchmark_config['output_to']     == getattr(configuration_yaml,'output')
    assert benchmark_config['output_format'] == getattr(configuration_yaml,'output_format')
    assert benchmark_config['output_failed'] == getattr(configuration_yaml,'failed')
    assert benchmark_config['output_report'] == getattr(configuration_yaml,'report')

    for i in range(len(benchmark_config['output_list'])):
        try: # Names should be equal
            assert( benchmark_config['output_list'][i]['name'] == 
                getattr(configuration_yaml,'outputs_array')[i]['name'] )
        except: # or they should be set in conf to something
            assert getattr(configuration_yaml,'outputs_array')[i]['name'] 
        # Output sequence is required, we make no attempt to check 'compile'
        assert ( benchmark_config['output_list'][i]['seq'] == 
                getattr(configuration_yaml,'outputs_array')[i]['seq'][0] )
        try: # The filter specification should be equal,
            assert ( benchmark_config['output_list'][i]['filter'] == 
                getattr(configuration_yaml,'outputs_array')[i]['filter'][0] )
        except: # or should be not False. True == does not work ... mysterious
            assert ( False != getattr(configuration_yaml,'outputs_array')[i]['filter'][0] )
        try: # The description specification should be equal,
            assert ( benchmark_config['output_list'][i]['description'] == 
                getattr(configuration_yaml,'outputs_array')[i]['description'][0] )
        except: # or set to 'description'
            assert ( 'description' == getattr(configuration_yaml,'outputs_array')[i]['description'][0] )

    for i in range(len(benchmark_config['matches'])):
        # 'use' is turned into 'input'
        assert( benchmark_config['matches'][i]['use'] == 
                getattr(configuration_yaml,'matches_array')[i]['input'] )
        for j in benchmark_config['matches'][i]['marked_groups']:
            try: # for all matches, what's the benchmark_config group name?
                bench_name = benchmark_config['matches'][i]['marked_groups'][j]['name']
            except: # or is generic untitled group name
                bench_name = r'untitled_group\d+'
            group_pattern = re.compile(r"\(\?<"+bench_name+">")
            # do we find that group name?
            assert( group_pattern.search(str(getattr(configuration_yaml,'matches_array')[i]['regex'] )))
            try: # next we move onto pulling in any repeat specs
                repeat_min = str(benchmark_config['matches'][i]['marked_groups'][j]['repeat_min'])
            except: # ... laboriously to avoid catching errors
                repeat_min = None
            try:
                repeat_max = str(benchmark_config['matches'][i]['marked_groups'][j]['repeat_max'])
            except:
                repeat_max = None
            try:
                repeat = str(benchmark_config['matches'][i]['marked_groups'][j]['repeat'])
            except:
                repeat = None
            if repeat_min and repeat_max: # then depending on what we have available, we search
                assert ( re.search(r"\{"+repeat_min+","+repeat_max+"}",
                        str(getattr(configuration_yaml,'matches_array')[i]['regex'] )))
            if repeat_min and repeat and ( repeat_max == None ): 
                assert ( re.search(r"\{"+repeat_min+","+repeat+"}",
                        str(getattr(configuration_yaml,'matches_array')[i]['regex'] )))
            if ( repeat_min == None ) and repeat and repeat_max: 
                assert ( re.search(r"\{"+repeat+","+repeat_max+"}",
                        str(getattr(configuration_yaml,'matches_array')[i]['regex'] )))
            if ( repeat_min == None ) and repeat and ( repeat_max == None ): 
                assert ( re.search(r"\{"+repeat+","+repeat+"}",
                        str(getattr(configuration_yaml,'matches_array')[i]['regex'] )))

def test_configuration_args(benchmark_config):
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

    assert benchmark_config['verbosity']     == getattr(conf,'verbosity')
    assert benchmark_config['input_from']    == getattr(conf,'input')
    assert benchmark_config['input_from']    == getattr(conf,'input')
    assert benchmark_config['input_format']  == getattr(conf,'input_format')
    assert benchmark_config['input_gzipped'] == getattr(conf,'gzipped')
    assert benchmark_config['output_to']     == getattr(conf,'output')
    assert benchmark_config['output_format'] == getattr(conf,'output_format')
    assert benchmark_config['output_format'] == getattr(conf,'output_format')
    assert benchmark_config['output_failed'] == getattr(conf,'failed')
    assert benchmark_config['output_report'] == getattr(conf,'report')


# Test MatchScores class
@pytest.fixture
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

# Test GroupStats class
@pytest.fixture
def groupstats():
    return itermae.GroupStats(5,15,
        SeqRecord.SeqRecord(Seq.Seq('ATCGATCGAT')),[36]*10)

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

def test_groupstats_repr(groupstats):
    assert type(repr(groupstats)) == type(str())


# Setup inputs
# I'm not testing BioPython SeqIO, I assume that's good.
# This instead uses that to return a list of the records in the file.
@pytest.fixture
def fastqfile():
    return SeqIO.parse("itermae/data/barseq.fastq","fastq")

## SeqHolder Tests

# Test that SeqHolder can apply_operation, then since we're there testing
# that it finds the right groups for each seq, and passes or fails filters 
# appropriately.
def test_seqholder_match_filter(fastqfile,configuration_yaml):
    for seq, pos_pass, qual_pass, seq_pass, sequences_found, seq_targets \
        in zip(fastqfile,
            [ i == 1 for i in [1,1,1,1,1,1,1,1,1] ],
            [ i == 1 for i in [0,0,1,1,0,0,0,1,1] ],
            [ i == 1 for i in [1,0,0,0,0,0,0,0,0] ],
            [   set(['dummyspacer','input','sample','fixed1','rest','tag','strain']),
                set(['dummyspacer','input','sample','fixed1','rest','tag','strain','fixed2','UMItail']),
                set(['dummyspacer','input','sample','fixed1','rest','tag','strain','fixed2','UMItail']),
                set(['dummyspacer','input','sample','fixed1','rest','tag','strain','fixed2','UMItail']),
                set(['dummyspacer','input','sample','fixed1','rest','tag','strain']),
                set(['dummyspacer','input','sample','fixed1','rest','tag','strain','fixed2','UMItail']),
                set(['dummyspacer','input','sample','fixed1','rest','tag','strain']),
                set(['dummyspacer','input','sample','fixed1','rest','tag','strain','fixed2','UMItail']),
                set(['dummyspacer','input','sample','fixed1','rest','tag','strain','fixed2','UMItail']),
            ],
            [   ('NB501157:100:H5J5LBGX2:1:11101:10000:10043_TTCAC', 'TCAGTCGTAGCAGTTCGATG'),
                ('NB501157:100:H5J5LBGX2:1:11101:10000:10138_GCTTC', 'TGGGCAGACACAACGCTACA'),
                ('NB501157:100:H5J5LBGX2:1:11101:10000:16613_GCTTC','GACAGACTGATAACCCTTGC'),
                ('NB501157:100:H5J5LBGX2:1:11101:10000:19701_CTACT', 'GATGCACTGCGTTCCATGTT'),
                ('NB501157:100:H5J5LBGX2:1:11101:10000:5096_TAAGT','AGGGCTCGTCGATTCGTCTT'),
                ('NB501157:100:H5J5LBGX2:1:11101:10000:6068_CTACT','GCAGATAATACACTGTCACC'),
                ('NB501157:100:H5J5LBGX2:1:11101:10000:8488_CATAA','TCGAGGGGTTACATACG'),
                ('NB501157:100:H5J5LBGX2:1:11101:10001:10798_TCTAG','GAGGCTACGGTACGTTCCTT'),
                ('NB501157:100:H5J5LBGX2:1:11101:10001:11700_CGCAA','TGCGCCACATAGTATAAAT'),
                ]
            ):
        # Read in the sequence to the holder
        seqholder = itermae.SeqHolder(seq,configuration=configuration_yaml)
        # Is the dummy X?
        assert seqholder.seqs['dummyspacer'].seq == 'X'
        # Is the number we just put there 40?
        assert seqholder.seqs['dummyspacer'].letter_annotations['phred_quality'] == [40] 
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
        third_filter = 'sample == "TTCAC" or sample == "AGGAG"'
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

## Operations for argument-specified options
#one_operation_string = (
#    '-m "input > (?P<sampleIndex>[ATCGN]{5,5})(?P<upPrime>GTCCTCGAGGTCTCT){e<=1}(?P<barcode>[ATCGN]{18,22})(?P<downPrime>CGTACGCTG){e<=1}" '+
#    '-os "barcode" -oi "id+\\"_\\"+sampleIndex" '
#    )
#two_operation_string = (
#    '-m "input > (?P<sampleIndex>[ATCGN]{5,5})(?P<rest>(GTCCTCGAGGTCTCT){e<=1}[ATCGN]*)" '+
#    '-m "rest  > (?P<upPrime>GTCCTCGAGGTCTCT){e<=1}(?P<barcode>[ATCGN]{18,22})(?P<downPrime>CGTACGCTG){e<=1}" '+
#    '-os "barcode" -oi "id+\\"_\\"+sampleIndex" '+
#    '-os "upPrime+barcode+downPrime" -oi "id+\\"_withFixedFlanking_\\"+sampleIndex" '
#    )

## Each operation applied to shortread FASTQ, non-gzipped
#def test_full_combinations():
#    operations_list = [one_operation_string, two_operation_string]
#    for input_format, which_ops, output_format in itertools.product(
#            ['fastq','fasta','sam','txt'],
#            [0,1],
#            ['fastq','fasta','sam','txt'],
#        ):
#        results = subprocess.run(
#            'cat itermae/data/barseq.'+input_format+' | '+
#            'itermae '+
#            '--input-format '+input_format+' '+
#            operations_list[which_ops]+
#            '--verbose --output-format '+output_format ,
#            shell=True,capture_output=True,encoding='utf-8')
#        filename = 'itermae/data/test_outputs/barseq_combinations_'
#        filename += 'input_'+input_format+'_'
#        filename += str(which_ops+1)+'_ops'
#        filename += '.'+output_format
##        with open(filename,'w') as f:
##            f.write(results.stdout)
#        with open(filename,'r') as f:
#            expected_file = f.readlines()
#        for i,j in zip(results.stdout.split('\n'),expected_file):
#            assert str(i) == str(j.rstrip('\n'))

## Each operation applied to shortread FASTQ, gzipped
#def test_full_combinations_gzipped():
#    operations_list = [one_operation_string, two_operation_string]
#    for input_format, which_ops, output_format in itertools.product(
#            ['fastq','fasta','sam','txt'],
#            [0,1],
#            ['fastq','fasta','sam','txt'],
#        ):
#        results = subprocess.run(
#            'itermae '+
#            '--input itermae/data/barseq.'+input_format+'.gz --gzipped '+
#            '--input-format '+input_format+' '+
#            operations_list[which_ops]+
#            '--verbose --output-format '+output_format ,
#            shell=True,capture_output=True,encoding='utf-8')
#        filename = 'itermae/data/test_outputs/barseq_combinations_'
#        filename += 'input_'+input_format+'_'
#        filename += str(which_ops+1)+'_ops'
#        filename += '.'+output_format
##        with open(filename,'w') as f:
##            f.write(results.stdout)
#        with open(filename,'r') as f:
#            expected_file = f.readlines()
#        for i,j in zip(results.stdout.split('\n'),expected_file):
#            assert str(i) == str(j.rstrip('\n'))

# Buncha full tests:
# Syntax of names is test_full_, then the test number, 
# then the type of matches (simple one op, or complex UMI bit),
# then the target type

input_dicts = [
    {   'input_from': 'itermae/data/tests/test_inputs/barseq.sam.gz',
        'input_format': 'sam', 'input_gzipped': 'true', 
        'has_desc':False, 'seq_as_id':False},
    {   'input_from': 'itermae/data/tests/test_inputs/barseq.fastq.gz',
        'input_format': 'fastq', 'input_gzipped': 'true',
        'has_desc':True, 'seq_as_id':False},
    {   'input_from': 'itermae/data/tests/test_inputs/barseq.fasta.gz',
        'input_format': 'fasta', 'input_gzipped': 'true',
        'has_desc':True, 'seq_as_id':False},
    {   'input_from': 'itermae/data/tests/test_inputs/barseq.txt.gz',
        'input_format': 'fastq', 'input_gzipped': 'true', 
        'has_desc':False, 'seq_as_id':True},
]

match_yaml_blocks = [
"""matches:
    -   use: 'input'
        pattern: 'NGTCCTCGAGGTCTCT'
        marking: 'ABBBBBBBBBBBBBBB'
        marked_groups:
            A:
                name: sampleIndex
                repeat: 5
            B:
                name: rest"""
,
"""matches:
    -   use: 'input'
        pattern: 'NGTCCTCGAGGTCTCT+'
        marking: 'ABBBBBBBBBBBBBBBB'
        marked_groups:
            A:
                name: sampleIndex
                repeat: 5
            B:
                name: rest
                allowed_errors: 1 
    -   use: rest 
        pattern: 'GTCCTCGAGGTCTCTNCGTACGCTG+'
        marking: 'AAAAAAAAAAAAAAABCCCCCCCCCD'
        marked_groups:
            A:
                name: upPrime
                allowed_errors: 1
            B:
                name: barcode
                repeat_min: 18
                repeat_max: 22
            C:
                name: downPrime
                allowed_errors: 1
            D:
                name: downstream
    -   use: downstream
        pattern: 'CGTACGCTGCAGGTCGACNGNANGNGNGNGAT'
        marking: 'AAAAAAAAAAAAAAAAAABBBBBBBBBBBCCC'
        marked_groups:
            A:
                name: fixed_pre_umi
                allowed_errors: 2
            B:
                name: interspersed_umi
                allowed_errors: 1
            C:
                name: tail
                allowed_errors: 1
"""
]

output_dicts = [
    {   'output_to': 'STDOUT', 'output_format': 'sam' },
    {   'output_to': 'STDOUT', 'output_format': 'fastq' },
    {   'output_to': 'STDOUT', 'output_format': 'fasta' },
    {   'output_to': 'STDOUT', 'output_format': 'txt' }
]

output_yaml_blocks = [        
"""output_list: 
    -   name: filtered_by_sampleIndex
        filter: \'sampleIndex == "GCTTC"\' 
        seq: 'input' 
"""
,
"""output_list: 
    -   name: barcode 
        filter: 'statistics.mean(barcode.quality) >= 25' 
        id: 'id+"_"+sampleIndex'
        seq: 'barcode' 
        description: 'description'
    -   name: sampleIndex 
        filter: 'sampleIndex.length >= 3'
        id: 'id+"_withFixedFlanking_"+sampleIndex'
        seq: 'upPrime+barcode+downPrime' 
"""
]

def making_a_full_test(config_file_path, 
        which_input, which_output, which_matches, which_outputs ):
    this_input_dict = input_dicts[which_input]
    this_output_dict = output_dicts[which_output]
    this_match_yaml_block = match_yaml_blocks[which_matches]
    this_output_yaml_block = output_yaml_blocks[which_outputs]
    config_file = config_file_path / "config.yml"
    config_file.write_text(
        'input_from: '+this_input_dict['input_from']+"\n"+
        'input_format: '+this_input_dict['input_format']+"\n"+
        'input_gzipped: '+this_input_dict['input_gzipped']+"\n"+
        this_match_yaml_block+"\n"+
        'output_to: '+this_output_dict['output_to']+"\n"+
        'output_format: '+this_output_dict['output_format']+"\n"+
        this_output_yaml_block
    )
    results = subprocess.run(
        'itermae -v --config '+str(config_file),
        shell=True,capture_output=True,encoding='utf-8')
    filename = ('itermae/data/tests/test_outputs/'+
        'matches-'+str(which_matches)+
        '_outputs-'+str(which_outputs)+
        '_hasID-'+  str(this_input_dict['seq_as_id'])+
        '_hasDesc-'+str(this_input_dict['has_desc'])+
        '.'+this_output_dict['output_format'])
    with open(filename,'w') as f:
        f.write(results.stdout)
    with open(filename,'r') as f:
        expected_file = f.readlines()
    for i,j in zip(results.stdout.split('\n'),expected_file):
        assert str(i) == str(j.rstrip('\n'))

def test_full_1(tmp_path):
    making_a_full_test(tmp_path,0,0,0,0)
        
