#!/usr/bin/env python3

import unittest

from itermae import *
from Bio import SeqIO, Seq, SeqRecord
import gzip

class TestItermae(unittest.TestCase):
    def setUp(self):
        self.input_fastq  = "example-data/test.fastq"
        self.input_fastqz = "example-data/test.fastqz"

    def test_MatchScores(self):
        this = MatchScores(1,2,3)
        self.assertEqual(this.substitutions,1)
        self.assertEqual(this.insertions,2)
        self.assertEqual(this.deletions,3)
        self.assertEqual(this.flatten(),"1_2_3")

    def test_GroupStats(self):
        this = GroupStats(5,15)
        self.assertEqual(this.start,5)
        self.assertEqual(this.end,15)
        self.assertEqual(this.length,10)
        self.assertEqual(this.flatten(),"5_15_10")

    def test_SeqHolder_initial(self):
        input_seqs = SeqIO.parse(self.input_fastq,"fastq")
        for i in input_seqs:
            this = SeqHolder(i,verbosity=3)
            self.assertEqual(this.verbosity,3)
            self.assertDictEqual(this.match_scores,{})
            self.assertDictEqual(this.group_stats,{})
            self.assertEqual(this.seqs['dummyspacer'].seq,'X')
            self.assertEqual(this.seqs['dummyspacer'].id,'dummyspacer')

    def test_SeqHolder_matchAndFilter(self):
        input_seqs = SeqIO.parse(self.input_fastq,"fastq")
        filter_results = []
        for i in input_seqs:
            this = SeqHolder(i,verbosity=0)
            this.apply_operation('a','input',
                regex.compile("(?P<sample>[ATCG]{5})(?P<fixed1>GTCCACGAGGTC){e<=1}(?P<rest>TCT.*){e<=1}",
                    regex.BESTMATCH) )
            this.apply_operation('b','rest',
                regex.compile("(?P<tag>TCT){e<=1}(?P<strain>[ATCG]{10,26})CGTACGCTGCAGGTCGAC",
                    regex.BESTMATCH) )
            this.apply_operation('c','rest',
                regex.compile("(?P<fixed2>CGTACGCTGCAGGTC)(?<UMItail>GAC[ATCG]G[ATCG]A[ATCG]G[ATCG]G[ATCG]G[ATCG]GAT){s<=2}",
                    regex.BESTMATCH) )
            this.apply_operation('a','UMItail',
                regex.compile("(GAC(?P<umi1>[ATCG])G(?<umi2>[ATCG])A(?<umi3>[ATCG])G(?<umi4>[ATCG])G(?<umi5>[ATCG])G(?<umi6>[ATCG])G){e<=2}",
                    regex.BESTMATCH) )
            filter_results.append(
                this.apply_filters(['sample.length == 5 and rest.start >= 15'])
                )
        self.assertListEqual(filter_results,[[True],[True],[False],[True]])

#        --output-seq "strain" \
#        --output-id "input.id+'_'+sample.seq" --report report \
#        --filter "sample.length == 5 and rest.start >= 16"  -v --output-format sam
            #self.assertEqual(this.verbosity,)
    

        
    
