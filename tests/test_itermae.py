#!/usr/bin/env python3

import pytest

import itermae
# Requires that you've installed it, so make the dist-files and install locally
# The Makefile rule `test` does this

import regex
from Bio import SeqIO

#
#
# Class Tests
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
#
# Full Tests
#
#

#
# Shortread FASTQ input
#

# Setup input

@pytest.fixture
def fastqfile():
    return SeqIO.parse("itermae/data/toy.fastq","fastq")
def test_seqholder_verbosity(fastqfile):
    for i in fastqfile:
        seqholder = itermae.SeqHolder(i,verbosity=3)
        assert seqholder.verbosity == 3
def test_seqholder_dummy(fastqfile):
    for i in fastqfile:
        seqholder = itermae.SeqHolder(i,verbosity=3)
        assert seqholder.seqs['dummyspacer'].seq == 'X'

def test_seqholder_match_filter(fastqfile):
    for seq, pos_pass, qual_pass, seq_pass in zip(fastqfile,
            [ [i == 1] for i in [1,1,0,1,0,0,1,1,1] ],
            [ [i == 1] for i in [0,0,0,1,0,0,1,1,0] ],
            [ [i == 1] for i in [1,0,0,0,0,0,0,1,0] ]
            ):
        seqholder = itermae.SeqHolder(seq,verbosity=3)
        seqholder.apply_operation('a','input',
            regex.compile("(?P<sample>[ATCG]{5})(?P<fixed1>GTCCACGAGGTC){e<=2}(?P<rest>TCT.*){e<=1}",
                regex.BESTMATCH) )
        seqholder.apply_operation('b','rest',
            regex.compile("(?P<tag>TCT){e<=1}(?P<strain>[ATCG]{10,26})(CGTACGCTGC){e<=2}",
                regex.BESTMATCH) )
        seqholder.apply_operation('c','rest',
            regex.compile("(?P<fixed2>CGTACGCTGCAGGTC)(?<UMItail>GAC[ATCG]G[ATCG]A[ATCG]G[ATCG]G[ATCG]G[ATCG]GAT){s<=2}",
                regex.BESTMATCH) )
        assert pos_pass == seqholder.apply_filters(['sample.length == 5 and rest.start >= 15']) 
        assert qual_pass == seqholder.apply_filters(['statistics.mean(fixed1.quality) >= 33.5'])
        assert seq_pass == seqholder.apply_filters(['sample.seq == "TTCAC" or sample.seq == "AGGAG"'])

def test_seqholder_outputs(fastqfile):
    for seq, seq_targets, filter_targets in zip(fastqfile,
            [ 
                ('ExampleToyReads:1:exampleFlowCell:1:10000:10000:10000_TTCAC','TCAGTCGTAGCAGTTCGATG'),
                ('ExampleToyReads:1:exampleFlowCell:1:10000:10001:10001_GCTTC', 'TGGCAGACACACGCTACA'),
                (None,None),
                ('ExampleToyReads:1:exampleFlowCell:1:10000:10003:10003_CTACT', 'GATGCACTGCGTTCCATGTT'),
                (None,None),
                (None,None),
                ('ExampleToyReads:1:exampleFlowCell:1:10000:10006:10006_TCGGC', 'ATTCTGAGCGGTGCCATAGT'),
                ('ExampleToyReads:1:exampleFlowCell:1:10000:10007:10007_AGGAG', 'ATAAGTTAGACAGGTCAGC'),
                ('ExampleToyReads:1:exampleFlowCell:1:10000:10008:10008_ACGTA', 'CACACGCACGAATTTGCATA')
                ],
                            [
                '"pass","ExampleToyReads:1:exampleFlowCell:1:10000:10000:10000","TTCACGTCCTCGAGGTCTCTTCAGTCGTAGCAGTTCGATGCGTACGCTACAGGTCGACGGTAAGAGAGGGATGTG","TCAGTCGTAGCAGTTCGATG","FilterResults","sample_0_5_5-fixed1_5_17_12-rest_17_75_58-tag_0_3_3-strain_3_23_20"',
                '"pass","ExampleToyReads:1:exampleFlowCell:1:10000:10001:10001","GCTTCGTCCTCGAGGTCTATTGGCAGACACACGCTACACGTACGCTGCAGGTCGAGGGCACGCGAGAGATGTGTG","TGGCAGACACACGCTACA","FilterResults","sample_0_5_5-fixed1_5_17_12-rest_17_75_58-tag_0_3_3-strain_3_21_18-fixed2_21_36_15-UMItail_36_53_17"',
                '"fail_to_form","ExampleToyReads:1:exampleFlowCell:1:10000:10002:10002","CCCGGCGTTCGGGGAAGGACGTCAATAGTCACACAGTCCTTGACGGTATAATAACCACCATCATGGCGACCATCC","TGGCAGACACACGCTACA","FilterResults",""',
                '"pass","ExampleToyReads:1:exampleFlowCell:1:10000:10003:10003","CTACTGTCCACGAGGTCTCTGATGCACTGCGTTCCATGTTCGTACGCTGCAGGTCGACGGAAGGAGCGCGATGTG","GATGCACTGCGTTCCATGTT","FilterResults","sample_0_5_5-fixed1_5_17_12-rest_17_75_58-tag_0_3_3-strain_3_23_20-fixed2_23_38_15-UMItail_38_55_17"',
                '"fail_to_form","ExampleToyReads:1:exampleFlowCell:1:10000:10004:10004","AAATTAGGGTCAACGCTACCTGTAGGAAGTGTCCGCATAAAGTGCACCGCATGGAAATGAAGACGGCCATTAGCT","GATGCACTGCGTTCCATGTT","FilterResults",""',
                '"fail_to_form","ExampleToyReads:1:exampleFlowCell:1:10000:10005:10005","CCCGGCGTTCGGGGAAGGACGTCAATAGTCACACAGTCCTTGACGGTATAATAACCACCATCATGGCGACCATCC","GATGCACTGCGTTCCATGTT","FilterResults",""',
                '"pass","ExampleToyReads:1:exampleFlowCell:1:10000:10006:10006","TCGGCGTCCTCGAGGTCTCTATTCTGAGCGGTGCCATAGTCGTACGCTGCAGGTCGACCGAAGGTGGGAGATGTG","ATTCTGAGCGGTGCCATAGT","FilterResults","sample_0_5_5-fixed1_5_17_12-rest_17_75_58-tag_0_3_3-strain_3_23_20-fixed2_23_38_15-UMItail_38_55_17"',
                '"pass","ExampleToyReads:1:exampleFlowCell:1:10000:10007:10007","AGGAGGTCCTCGAGGTCTCTATAAGTTAGACAGGTCAGCCGTACGCTGCAGGTCGACAGCTGGCGCGCGATGTGA","ATAAGTTAGACAGGTCAGC","FilterResults","sample_0_5_5-fixed1_5_17_12-rest_17_75_58-tag_0_3_3-strain_3_22_19-fixed2_22_37_15-UMItail_37_54_17"',
                '"pass","ExampleToyReads:1:exampleFlowCell:1:10000:10008:10008","ACGTAGTCCACGAGGTCTCTCACACGCACGAATTTGCATACGTACGCTGCAGGTCGACTGGAAGGGCGGGATGTG","CACACGCACGAATTTGCATA","FilterResults","sample_0_5_5-fixed1_5_17_12-rest_17_75_58-tag_0_3_3-strain_3_23_20-fixed2_23_38_15-UMItail_38_55_17"'
                ],
            ):
        seqholder = itermae.SeqHolder(seq,verbosity=3)
        seqholder.apply_operation('a','input',
            regex.compile("(?P<sample>[ATCG]{5})(?P<fixed1>GTCCACGAGGTC){e<=2}(?P<rest>TCT.*){e<=1}",
                regex.BESTMATCH) )
        seqholder.apply_operation('b','rest',
            regex.compile("(?P<tag>TCT){e<=1}(?P<strain>[ATCG]{10,26})(CGTACGCTGC){e<=2}",
                regex.BESTMATCH) )
        seqholder.apply_operation('c','rest',
            regex.compile("(?P<fixed2>CGTACGCTGCAGGTC)(?<UMItail>GAC[ATCG]G[ATCG]A[ATCG]G[ATCG]G[ATCG]G[ATCG]GAT){s<=2}",
                regex.BESTMATCH) )
        (this_output_id, this_output_seq) = (None,None)
        try:
            this_output = seqholder.build_output(
                    "input.id+'_'+sample.seq",
                    "strain" ) 
            (this_output_id, this_output_seq) = (this_output.id,this_output.seq)
            this_result = seqholder.format_report("pass",this_output.seq,"FilterResults")
        except:
            this_result = seqholder.format_report("fail_to_form",this_output.seq,"FilterResults")
        assert seq_targets == ( this_output_id, this_output_seq ) 
        assert filter_targets == this_result
