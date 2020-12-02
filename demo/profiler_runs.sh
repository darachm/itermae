#!/usr/bin/env bash

# Demos with profiling
# Single core of course

# First, simple one operation one output
seq 10 | xargs -I{} cat ../example-data/barseq.fastq \
    | python3 -m cProfile -o simple_opOne_outOne.profile \
        ../bin/itermae \
            -o "input > (?P<sampleIndex>[ATCGN]{5,5})(?P<rest>(GTCCACGAG){e<=1}[ATCGN]*)" \
            -oseq "rest" -oid "input.id+\"_\"+sampleIndex.seq" \
            -of "fasta" \
    > /dev/null

# Long complex operation, one output
seq 10 | xargs -I{} cat ../example-data/barseq.fastq \
    | python3 -m cProfile -o long_opOne_outOne.profile \
        ../bin/itermae \
            -o "input > (?P<sampleIndex>[ATCGN]{5,5})[ATCGN]*(?P<upPrime>GTCCACGAGGTCTCT){e<=2}(?P<barcode>[ATCGN]{18,22})(?P<downPrime>CGTACGCTG){e<=1}" \
            -oseq "barcode" -oid "input.id+\"_\"+sampleIndex.seq" \
            -of "fasta" \
    > /dev/null

# Two simpler operations, one output
seq 10 | xargs -I{} cat ../example-data/barseq.fastq \
    | python3 -m cProfile -o opTwo_outOne.profile \
        ../bin/itermae \
            -o "input > (?P<sampleIndex>[ATCGN]{5,5})(?P<rest>(GTCCACGAG){e<=1}[ATCGN]*)" \
            -o "rest  > (?P<upPrime>GTCCACGAGGTCTCT){e<=2}(?P<barcode>[ATCGN]{18,22})(?P<downPrime>CGTACGCTG){e<=1}" \
            -oseq "barcode" -oid "input.id+\"_\"+sampleIndex.seq" \
            -of "fasta" \
    > /dev/null

# Two simpler operations, two outputs
seq 10 | xargs -I{} cat ../example-data/barseq.fastq \
    | python3 -m cProfile -o opTwo_outTwo.profile \
        ../bin/itermae \
            -o "input > (?P<sampleIndex>[ATCGN]{5,5})(?P<rest>(GTCCACGAG){e<=1}[ATCGN]*)" \
            -o "rest  > (?P<upPrime>GTCCACGAGGTCTCT){e<=2}(?P<barcode>[ATCGN]{18,22})(?P<downPrime>CGTACGCTG){e<=1}" \
            -oseq "barcode" -oid "input.id+\"_\"+sampleIndex.seq" \
            -oseq "upPrime+barcode+downPrime" -oid "input.id+\"_withFixedFlanking_\"+sampleIndex.seq" \
            -of "fasta" \
    > /dev/null

