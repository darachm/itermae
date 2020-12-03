#!/usr/bin/env python3

import pytest

import itermae

import regex
from Bio import SeqIO

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

@pytest.fixture
def fastqfile():
    return SeqIO.parse("example-data/toy.fastq","fastq")
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
        try:
            seqholder.apply_operation('a','input',
                regex.compile("(?P<sample>[ATCG]{5})(?P<fixed1>GTCCACGAGGTC){e<=2}(?P<rest>TCT.*){e<=1}",
                    regex.BESTMATCH) )
            seqholder.apply_operation('b','rest',
                regex.compile("(?P<tag>TCT){e<=1}(?P<strain>[ATCG]{10,26})(CGTACGCTGC){e<=2}",
                    regex.BESTMATCH) )
            seqholder.apply_operation('c','rest',
                regex.compile("(?P<fixed2>CGTACGCTGCAGGTC)(?<UMItail>GAC[ATCG]G[ATCG]A[ATCG]G[ATCG]G[ATCG]G[ATCG]GAT){s<=2}",
                    regex.BESTMATCH) )
        except:
            pass
        assert pos_pass == seqholder.apply_filters(['sample.length == 5 and rest.start >= 15']) 
        assert qual_pass == seqholder.apply_filters(['statistics.mean(fixed1.quality) >= 33.5'])
        assert seq_pass == seqholder.apply_filters(['sample.seq == "TTCAC" or sample.seq == "AGGAG"'])


#            try:
#                this_output = this.build_output(
#                        "input.id+'_'+sample.seq+'_'+umi1.seq+umi2.seq+umi3.seq",
#                        "strain" ) 
#                outputs.append( ( this_output.id, this_output.seq ) )
#                reports.append(this.format_report("pass",this_output.seq,
#                        filter_results))
#            except:
#                outputs.append(None)
#                reports.append(this.format_report("fail",this.seqs['input'],
#                        filter_results))
#        self.assertListEqual(filter_results,[[True],[True],[False],[True]])
#        self.assertListEqual(outputs,[None,None,None,
#                ('NB501157:100:H5J5LBGX2:1:11101:10000:19701_CTACT_GAG',
#                    'GATGCACTGCGTTCCATGTT')
#                ])
#
#if __name__ == '__main__':
#    unittest.main()
