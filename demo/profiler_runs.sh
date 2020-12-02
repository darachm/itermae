#!/usr/bin/env bash

# Demos with profiling
# Single core of course

# First, simple one operation one output
seq 50 | xargs -I{} cat example-data/barseq.fastq \
    | python3 -m cProfile -o demo/simple_opOne_outOne.profile \
        bin/itermae \
            -o "input > (?P<sampleIndex>[ATCGN]{5,5})(?P<rest>(GTCCACGAG){e<=1}[ATCGN]*)" \
            -oseq "rest" -oid "input.id+\"_\"+sampleIndex.seq" \
            -of "fasta" \
    > demo/simple_opOne_outOne.fasta

# Long complex operation, one output
seq 50 | xargs -I{} cat example-data/barseq.fastq \
    | python3 -m cProfile -o demo/long_opOne_outOne.profile \
        bin/itermae \
            -o "input > (?P<sampleIndex>[ATCGN]{5,5})(?P<upPrime>GTCCTCGAGGTCTCT){e<=1}(?P<barcode>[ATCGN]{18,22})(?P<downPrime>CGTACGCTG){e<=1}" \
            -oseq "barcode" -oid "input.id+\"_\"+sampleIndex.seq" \
            -of "fasta" \
    > demo/long_opOne_outOne.fasta

# Two simpler operations, one output
seq 50 | xargs -I{} cat example-data/barseq.fastq \
    | python3 -m cProfile -o demo/opTwo_outOne.profile \
        bin/itermae \
            -o "input > (?P<sampleIndex>[ATCGN]{5,5})(?P<rest>(GTCCTCGAGGTCTCT){e<=1}[ATCGN]*)" \
            -o "rest  > (?P<upPrime>GTCCTCGAGGTCTCT){e<=1}(?P<barcode>[ATCGN]{18,22})(?P<downPrime>CGTACGCTG){e<=1}" \
            -oseq "barcode" -oid "input.id+\"_\"+sampleIndex.seq" \
            -of "fasta" \
    > demo/opTwo_outOne.fasta

# Two simpler operations, two outputs
seq 50 | xargs -I{} cat example-data/barseq.fastq \
    | python3 -m cProfile -o demo/opTwo_outTwo.profile \
        bin/itermae \
            -o "input > (?P<sampleIndex>[ATCGN]{5,5})(?P<rest>(GTCCTCGAGGTCTCT){e<=1}[ATCGN]*)" \
            -o "rest  > (?P<upPrime>GTCCTCGAGGTCTCT){e<=1}(?P<barcode>[ATCGN]{18,22})(?P<downPrime>CGTACGCTG){e<=1}" \
            -oseq "barcode" -oid "input.id+\"_\"+sampleIndex.seq" \
            -oseq "upPrime+barcode+downPrime" -oid "input.id+\"_withFixedFlanking_\"+sampleIndex.seq" \
            -of "fasta" \
    > demo/opTwo_outTwo.fasta

