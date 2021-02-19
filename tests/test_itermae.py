#!/usr/bin/env python3

import pytest

import itermae

# Required for testing apply_operation, as regex is generated in the command
# line program
import regex
# Required for testing SeqHolder etc. without the file IO of main script
from Bio import SeqIO
# Required for full-file input/output testing
import subprocess

#
#
# Ye Tests
#
#

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

#
# SeqHolder Tests
#

# Test SeqHolder verbosity
def test_seqholder_verbosity(fastqfile):
    for i in fastqfile:
        seqholder = itermae.SeqHolder(i,verbosity=3)
        # Is the verbosity right? (this is a gimme)
        assert seqholder.verbosity == 3

# Test SeqHolder dummyspacer is right
def test_seqholder_dummy(fastqfile):
    for i in fastqfile:
        seqholder = itermae.SeqHolder(i,verbosity=3)
        # Is the sequence X?
        assert seqholder.seqs['dummyspacer'].seq == 'X'
        # Is the number we just put there 40?
        assert seqholder.seqs['dummyspacer'].letter_annotations['phred_quality'] == [40] 

# Setup inputs
# I'm not testing BioPython SeqIO, I assume that's good.
# This instead uses that to return a list of the records in the file.
@pytest.fixture
def fastqfile():
    return SeqIO.parse("itermae/data/toy.fastq","fastq")

# Test that SeqHolder can apply_operation, then since we're there testing
# that it finds the right groups for each seq, and passes or fails filters 
# appropriately.
# Shortread FASTQ input
def test_seqholder_match_filter(fastqfile):
    for seq, pos_pass, qual_pass, seq_pass, sequences_found, \
        seq_targets, report_targets \
        in zip(fastqfile,
            [ [i == 1] for i in [1,1,0,1,0,0,1,1,1] ],
            [ [i == 1] for i in [0,0,0,1,0,0,1,1,0] ],
            [ [i == 1] for i in [1,0,0,0,0,0,0,1,0] ],
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
        assert pos_pass == seqholder.apply_filters(['sample.length == 5 and rest.start >= 15']) 
        # Does it pass a quality filter, with statistics?
        assert qual_pass == seqholder.apply_filters(['statistics.mean(fixed1.quality) >= 33.5'])
        # Does it pass a specific sequence filter?
        assert seq_pass == seqholder.apply_filters(['sample.seq == "TTCAC" or sample.seq == "AGGAG"'])
        # Then build outputs
        (this_output_id, this_output_seq) = (None,None)
        try:
            this_output = seqholder.build_output(
                    "input.id+'_'+sample.seq",
                    "strain" ) 
            (this_output_id, this_output_seq) = (this_output.id,this_output.seq)
            this_report = seqholder.format_report("pass",this_output.seq,"FilterResults")
        except:
            this_report = seqholder.format_report("fail_to_form",this_output.seq,"FilterResults")
        # Are the right outputs constructed?
        assert seq_targets == ( this_output_id, this_output_seq ) 
        # Are the right filter reports constructed?
        assert report_targets == this_report

#
#
# Full Tests
#
#

#
# shortread FASTQ, one operation
#

# As it says on the tin
def test_full_shortread_FASTQ_one_operation_to_fasta():
    results = subprocess.run(
        'cat itermae/data/barseq.fastq | '+
            'itermae '+
            '-o "input > (?P<sampleIndex>[ATCGN]{5,5})(?P<upPrime>GTCCTCGAGGTCTCT){e<=1}(?P<barcode>[ATCGN]{18,22})(?P<downPrime>CGTACGCTG){e<=1}" '+
            '-oseq "barcode" -oid "input.id+\\"_\\"+sampleIndex.seq" '+
            '--verbose -of "fasta"',
        shell=True,capture_output=True,encoding='utf-8')
#    with open('itermae/data/barseq_one_operation.fasta','w') as f:
#        f.write(results.stdout)
    with open('itermae/data/barseq_one_operation.fasta','r') as f:
        expected_file = f.readlines()
    for i,j in zip(results.stdout.split('\n'),expected_file):
        assert str(i) == str(j.rstrip('\n'))

# As it says on the tin
def test_full_shortread_FASTQ_one_operation_to_fastq():
    results = subprocess.run(
        'cat itermae/data/barseq.fastq | '+
            'itermae '+
            '-o "input > (?P<sampleIndex>[ATCGN]{5,5})(?P<upPrime>GTCCTCGAGGTCTCT){e<=1}(?P<barcode>[ATCGN]{18,22})(?P<downPrime>CGTACGCTG){e<=1}" '+
            '-oseq "barcode" -oid "input.id+\\"_\\"+sampleIndex.seq" '+
            '--verbose -of "fastq"',
        shell=True,capture_output=True,encoding='utf-8')
#    with open('itermae/data/barseq_one_operation.fastq','w') as f:
#        f.write(results.stdout)
    with open('itermae/data/barseq_one_operation.fastq','r') as f:
        expected_file = f.readlines()
    for i,j in zip(results.stdout.split('\n'),expected_file):
        assert str(i) == str(j.rstrip('\n'))

# As it says on the tin
def test_full_shortread_FASTQ_one_operation_to_sam():
    results = subprocess.run(
        'cat itermae/data/barseq.fastq | '+
            'itermae '+
            '-o "input > (?P<sampleIndex>[ATCGN]{5,5})(?P<upPrime>GTCCTCGAGGTCTCT){e<=1}(?P<barcode>[ATCGN]{18,22})(?P<downPrime>CGTACGCTG){e<=1}" '+
            '-oseq "barcode" -oid "input.id+\\"_\\"+sampleIndex.seq" '+
            '--verbose -of "sam"',
        shell=True,capture_output=True,encoding='utf-8')
#    with open('itermae/data/barseq_one_operation.sam','w') as f:
#        f.write(results.stdout)
    with open('itermae/data/barseq_one_operation.sam','r') as f:
        expected_file = f.readlines()
    for i,j in zip(results.stdout.split('\n'),expected_file):
        assert str(i) == str(j.rstrip('\n'))

#
# shortread FASTQ, two operation
#

# As it says on the tin
def test_full_shortread_FASTQ_two_operations_to_fasta():
    results = subprocess.run(
        'cat itermae/data/barseq.fastq | '+
            'itermae '+
            '-o "input > (?P<sampleIndex>[ATCGN]{5,5})(?P<rest>(GTCCTCGAGGTCTCT){e<=1}[ATCGN]*)" '+
            '-o "rest  > (?P<upPrime>GTCCTCGAGGTCTCT){e<=1}(?P<barcode>[ATCGN]{18,22})(?P<downPrime>CGTACGCTG){e<=1}" '+
            '-oseq "barcode" -oid "input.id+\\"_\\"+sampleIndex.seq" '+
            '-oseq "upPrime+barcode+downPrime" -oid "input.id+\\"_withFixedFlanking_\\"+sampleIndex.seq" '+
            '--verbose -of "fasta"',
        shell=True,capture_output=True,encoding='utf-8')
#    with open('itermae/data/barseq_two_operations.fasta','w') as f:
#        f.write(results.stdout)
    with open('itermae/data/barseq_two_operations.fasta','r') as f:
        expected_file = f.readlines()
    for i,j in zip(results.stdout.split('\n'),expected_file):
        assert str(i) == str(j.rstrip('\n'))

# As it says on the tin
def test_full_shortread_FASTQ_two_operations_to_fastq():
    results = subprocess.run(
        'cat itermae/data/barseq.fastq | '+
            'itermae '+
            '-o "input > (?P<sampleIndex>[ATCGN]{5,5})(?P<rest>(GTCCTCGAGGTCTCT){e<=1}[ATCGN]*)" '+
            '-o "rest  > (?P<upPrime>GTCCTCGAGGTCTCT){e<=1}(?P<barcode>[ATCGN]{18,22})(?P<downPrime>CGTACGCTG){e<=1}" '+
            '-oseq "barcode" -oid "input.id+\\"_\\"+sampleIndex.seq" '+
            '-oseq "upPrime+barcode+downPrime" -oid "input.id+\\"_withFixedFlanking_\\"+sampleIndex.seq" '+
            '--verbose -of "fastq"',
        shell=True,capture_output=True,encoding='utf-8')
#    with open('itermae/data/barseq_two_operations.fastq','w') as f:
#        f.write(results.stdout)
    with open('itermae/data/barseq_two_operations.fastq','r') as f:
        expected_file = f.readlines()
    for i,j in zip(results.stdout.split('\n'),expected_file):
        assert str(i) == str(j.rstrip('\n'))

# As it says on the tin
def test_full_shortread_FASTQ_two_operations_to_sam():
    results = subprocess.run(
        'cat itermae/data/barseq.fastq | '+
            'itermae '+
            '-o "input > (?P<sampleIndex>[ATCGN]{5,5})(?P<rest>(GTCCTCGAGGTCTCT){e<=1}[ATCGN]*)" '+
            '-o "rest  > (?P<upPrime>GTCCTCGAGGTCTCT){e<=1}(?P<barcode>[ATCGN]{18,22})(?P<downPrime>CGTACGCTG){e<=1}" '+
            '-oseq "barcode" -oid "input.id+\\"_\\"+sampleIndex.seq" '+
            '-oseq "upPrime+barcode+downPrime" -oid "input.id+\\"_withFixedFlanking_\\"+sampleIndex.seq" '+
            '--verbose -of "sam"',
        shell=True,capture_output=True,encoding='utf-8')
#    with open('itermae/data/barseq_two_operations.sam','w') as f:
#        f.write(results.stdout)
    with open('itermae/data/barseq_two_operations.sam','r') as f:
        expected_file = f.readlines()
    for i,j in zip(results.stdout.split('\n'),expected_file):
        assert str(i) == str(j.rstrip('\n'))

#
# shortread SAM, two operation
#

# As it says on the tin
def test_full_shortread_SAM_two_operations_to_fasta():
    results = subprocess.run(
        'cat itermae/data/barseq.sam | '+
            'itermae '+
            '-if sam '+
            '-o "input > (?P<sampleIndex>[ATCGN]{5,5})(?P<rest>(GTCCTCGAGGTCTCT){e<=1}[ATCGN]*)" '+
            '-o "rest  > (?P<upPrime>GTCCTCGAGGTCTCT){e<=1}(?P<barcode>[ATCGN]{18,22})(?P<downPrime>CGTACGCTG){e<=1}" '+
            '-oseq "barcode" -oid "input.id+\\"_\\"+sampleIndex.seq" '+
            '-oseq "upPrime+barcode+downPrime" -oid "input.id+\\"_withFixedFlanking_\\"+sampleIndex.seq" '+
            '--verbose -of "fasta"',
        shell=True,capture_output=True,encoding='utf-8')
#    with open('itermae/data/barseq_sam_two_operations.fasta','w') as f:
#        f.write(results.stdout)
    with open('itermae/data/barseq_two_operations.fasta','r') as f:
        expected_file = f.readlines()
    for i,j in zip(results.stdout.split('\n'),expected_file):
        assert str(i) == str(j.rstrip('\n'))

# As it says on the tin
def test_full_shortread_SAM_two_operations_to_fastq():
    results = subprocess.run(
        'cat itermae/data/barseq.sam | '+
            'itermae '+
            '-if sam '+
            '-o "input > (?P<sampleIndex>[ATCGN]{5,5})(?P<rest>(GTCCTCGAGGTCTCT){e<=1}[ATCGN]*)" '+
            '-o "rest  > (?P<upPrime>GTCCTCGAGGTCTCT){e<=1}(?P<barcode>[ATCGN]{18,22})(?P<downPrime>CGTACGCTG){e<=1}" '+
            '-oseq "barcode" -oid "input.id+\\"_\\"+sampleIndex.seq" '+
            '-oseq "upPrime+barcode+downPrime" -oid "input.id+\\"_withFixedFlanking_\\"+sampleIndex.seq" '+
            '--verbose -of "fastq"',
        shell=True,capture_output=True,encoding='utf-8')
#    with open('itermae/data/barseq_sam_two_operations.fastq','w') as f:
#        f.write(results.stdout)
    with open('itermae/data/barseq_two_operations.fastq','r') as f:
        expected_file = f.readlines()
    for i,j in zip(results.stdout.split('\n'),expected_file):
        assert str(i) == str(j.rstrip('\n'))

# As it says on the tin
def test_full_shortread_SAM_two_operations_to_sam():
    results = subprocess.run(
        'cat itermae/data/barseq.sam | '+
            'itermae '+
            '-if sam '+
            '-o "input > (?P<sampleIndex>[ATCGN]{5,5})(?P<rest>(GTCCTCGAGGTCTCT){e<=1}[ATCGN]*)" '+
            '-o "rest  > (?P<upPrime>GTCCTCGAGGTCTCT){e<=1}(?P<barcode>[ATCGN]{18,22})(?P<downPrime>CGTACGCTG){e<=1}" '+
            '-oseq "barcode" -oid "input.id+\\"_\\"+sampleIndex.seq" '+
            '-oseq "upPrime+barcode+downPrime" -oid "input.id+\\"_withFixedFlanking_\\"+sampleIndex.seq" '+
            '--verbose -of "sam"',
        shell=True,capture_output=True,encoding='utf-8')
#    with open('itermae/data/barseq_sam_two_operations.sam','w') as f:
#        f.write(results.stdout)
    with open('itermae/data/barseq_two_operations.sam','r') as f:
        expected_file = f.readlines()
    for i,j in zip(results.stdout.split('\n'),expected_file):
        assert str(i) == str(j.rstrip('\n'))

#
# shortread FASTA, two operation
#

# As it says on the tin
def test_full_shortread_fasta_two_operations_to_fasta():
    results = subprocess.run(
        'cat itermae/data/barseq.fasta | '+
            'itermae '+
            '-if fasta '+
            '-o "input > (?P<sampleIndex>[ATCGN]{5,5})(?P<rest>(GTCCTCGAGGTCTCT){e<=1}[ATCGN]*)" '+
            '-o "rest  > (?P<upPrime>GTCCTCGAGGTCTCT){e<=1}(?P<barcode>[ATCGN]{18,22})(?P<downPrime>CGTACGCTG){e<=1}" '+
            '-oseq "barcode" -oid "input.id+\\"_\\"+sampleIndex.seq" '+
            '-oseq "upPrime+barcode+downPrime" -oid "input.id+\\"_withFixedFlanking_\\"+sampleIndex.seq" '+
            '--verbose -of "fasta"',
        shell=True,capture_output=True,encoding='utf-8')
    with open('itermae/data/barseq_two_operations.fasta','r') as f:
        expected_file = f.readlines()
    for i,j in zip(results.stdout.split('\n'),expected_file):
        assert str(i) == str(j.rstrip('\n'))

# As it says on the tin
def test_full_shortread_fasta_two_operations_to_fastq():
    results = subprocess.run(
        'cat itermae/data/barseq.fasta | '+
            'itermae '+
            '-if fasta '+
            '-o "input > (?P<sampleIndex>[ATCGN]{5,5})(?P<rest>(GTCCTCGAGGTCTCT){e<=1}[ATCGN]*)" '+
            '-o "rest  > (?P<upPrime>GTCCTCGAGGTCTCT){e<=1}(?P<barcode>[ATCGN]{18,22})(?P<downPrime>CGTACGCTG){e<=1}" '+
            '-oseq "barcode" -oid "input.id+\\"_\\"+sampleIndex.seq" '+
            '-oseq "upPrime+barcode+downPrime" -oid "input.id+\\"_withFixedFlanking_\\"+sampleIndex.seq" '+
            '--verbose -of "fastq"',
        shell=True,capture_output=True,encoding='utf-8')
#    with open('itermae/data/barseq_max_quality_two_operations.fastq','w') as f:
#        f.write(results.stdout)
    with open('itermae/data/barseq_max_quality_two_operations.fastq','r') as f:
        expected_file = f.readlines()
    for i,j in zip(results.stdout.split('\n'),expected_file):
        assert str(i) == str(j.rstrip('\n'))

# As it says on the tin
def test_full_shortread_fasta_two_operations_to_fasta():
    results = subprocess.run(
        'cat itermae/data/barseq.fasta | '+
            'itermae '+
            '-if fasta '+
            '-o "input > (?P<sampleIndex>[ATCGN]{5,5})(?P<rest>(GTCCTCGAGGTCTCT){e<=1}[ATCGN]*)" '+
            '-o "rest  > (?P<upPrime>GTCCTCGAGGTCTCT){e<=1}(?P<barcode>[ATCGN]{18,22})(?P<downPrime>CGTACGCTG){e<=1}" '+
            '-oseq "barcode" -oid "input.id+\\"_\\"+sampleIndex.seq" '+
            '-oseq "upPrime+barcode+downPrime" -oid "input.id+\\"_withFixedFlanking_\\"+sampleIndex.seq" '+
            '--verbose -of "sam"',
        shell=True,capture_output=True,encoding='utf-8')
#    with open('itermae/data/barseq_max_quality_two_operations.sam','w') as f:
#        f.write(results.stdout)
    with open('itermae/data/barseq_max_quality_two_operations.sam','r') as f:
        expected_file = f.readlines()
    for i,j in zip(results.stdout.split('\n'),expected_file):
        assert str(i) == str(j.rstrip('\n'))


#
# shortread raw txt, two operation
#

# As it says on the tin
def test_full_shortread_txt_two_operations_to_fasta():
    results = subprocess.run(
        'cat itermae/data/barseq.txt | '+
            'itermae '+
            '-if txt '+
            '-o "input > (?P<sampleIndex>[ATCGN]{5,5})(?P<rest>(GTCCTCGAGGTCTCT){e<=1}[ATCGN]*)" '+
            '-o "rest  > (?P<upPrime>GTCCTCGAGGTCTCT){e<=1}(?P<barcode>[ATCGN]{18,22})(?P<downPrime>CGTACGCTG){e<=1}" '+
            '-oseq "barcode" -oid "input.id+\\"_\\"+sampleIndex.seq" '+
            '-oseq "upPrime+barcode+downPrime" -oid "input.id+\\"_withFixedFlanking_\\"+sampleIndex.seq" '+
            '--verbose -of "fasta"',
        shell=True,capture_output=True,encoding='utf-8')
#    with open('itermae/data/barseq_max_quality_seq_is_id_two_operations.fasta','w') as f:
#        f.write(results.stdout)
    with open('itermae/data/barseq_max_quality_seq_is_id_two_operations.fasta','r') as f:
        expected_file = f.readlines()
    for i,j in zip(results.stdout.split('\n'),expected_file):
        assert str(i) == str(j.rstrip('\n'))

# As it says on the tin
def test_full_shortread_txt_two_operations_to_fastq():
    results = subprocess.run(
        'cat itermae/data/barseq.txt | '+
            'itermae '+
            '-if txt '+
            '-o "input > (?P<sampleIndex>[ATCGN]{5,5})(?P<rest>(GTCCTCGAGGTCTCT){e<=1}[ATCGN]*)" '+
            '-o "rest  > (?P<upPrime>GTCCTCGAGGTCTCT){e<=1}(?P<barcode>[ATCGN]{18,22})(?P<downPrime>CGTACGCTG){e<=1}" '+
            '-oseq "barcode" -oid "input.id+\\"_\\"+sampleIndex.seq" '+
            '-oseq "upPrime+barcode+downPrime" -oid "input.id+\\"_withFixedFlanking_\\"+sampleIndex.seq" '+
            '--verbose -of "fastq"',
        shell=True,capture_output=True,encoding='utf-8')
#    with open('itermae/data/barseq_max_quality_seq_is_id_two_operations.fastq','w') as f:
#        f.write(results.stdout)
    with open('itermae/data/barseq_max_quality_seq_is_id_two_operations.fastq','r') as f:
        expected_file = f.readlines()
    for i,j in zip(results.stdout.split('\n'),expected_file):
        assert str(i) == str(j.rstrip('\n'))

# As it says on the tin
def test_full_shortread_txt_two_operations_to_sam():
    results = subprocess.run(
        'cat itermae/data/barseq.txt | '+
            'itermae '+
            '-if txt '+
            '-o "input > (?P<sampleIndex>[ATCGN]{5,5})(?P<rest>(GTCCTCGAGGTCTCT){e<=1}[ATCGN]*)" '+
            '-o "rest  > (?P<upPrime>GTCCTCGAGGTCTCT){e<=1}(?P<barcode>[ATCGN]{18,22})(?P<downPrime>CGTACGCTG){e<=1}" '+
            '-oseq "barcode" -oid "input.id+\\"_\\"+sampleIndex.seq" '+
            '-oseq "upPrime+barcode+downPrime" -oid "input.id+\\"_withFixedFlanking_\\"+sampleIndex.seq" '+
            '--verbose -of "sam"',
        shell=True,capture_output=True,encoding='utf-8')
#    with open('itermae/data/barseq_max_quality_seq_is_id_two_operations.sam','w') as f:
#        f.write(results.stdout)
    with open('itermae/data/barseq_max_quality_seq_is_id_two_operations.sam','r') as f:
        expected_file = f.readlines()
    for i,j in zip(results.stdout.split('\n'),expected_file):
        assert str(i) == str(j.rstrip('\n'))






##### Template for full functional tests
##### As it says on the tin
####def test_full_shortread_FASTQ_two_operations_to_sam():
####    results = subprocess.run(
####        'cat itermae/data/barseq.fastq | '+
####            'itermae '+
####            '-o "input > (?P<sampleIndex>[ATCGN]{5,5})(?P<rest>(GTCCTCGAGGTCTCT){e<=1}[ATCGN]*)" '+
####            '-o "rest  > (?P<upPrime>GTCCTCGAGGTCTCT){e<=1}(?P<barcode>[ATCGN]{18,22})(?P<downPrime>CGTACGCTG){e<=1}" '+
####            '-oseq "barcode" -oid "input.id+\\"_\\"+sampleIndex.seq" '+
####            '-oseq "upPrime+barcode+downPrime" -oid "input.id+\\"_withFixedFlanking_\\"+sampleIndex.seq" '+
####            '--verbose -of "sam"',
####        shell=True,capture_output=True,encoding='utf-8')
####    with open('itermae/data/barseq_two_operations.fastq','w') as f:
####        f.write(results.stdout)
#####    with open('itermae/data/barseq_two_operations.fastq','r') as f:
#####        expected_file = f.readlines()
#####    for i,j in zip(results.stdout.split('\n'),expected_file):
#####        assert str(i) == str(j.rstrip('\n'))
####


