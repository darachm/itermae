#!/usr/bin/env bash

# Demos with profiling.
# Single core of course.
# 'seq' and 'xargs' is just used to repeat the input file 50 times

# First, simple one operation one output
seq 50 | xargs -I{} cat itermae/data/barseq.fastq \
    | python3 -m cProfile -o profiling_tests/simple_opOne_outOne.profile \
        bin/itermae \
            -m "input > (?P<sampleIndex>[ATCGN]{5,5})(?P<rest>(GTCCACGAG){e<=1}[ATCGN]*)" \
            -os "rest" -oi "input.id+\"_\"+sampleIndex.seq" \
            --output-format "fasta" \
    > profiling_tests/simple_opOne_outOne.fasta

# Long complex operation, one output
seq 50 | xargs -I{} cat itermae/data/barseq.fastq \
    | python3 -m cProfile -o profiling_tests/long_opOne_outOne.profile \
        bin/itermae \
            -m "input > (?P<sampleIndex>[ATCGN]{5,5})(?P<upPrime>GTCCTCGAGGTCTCT){e<=1}(?P<barcode>[ATCGN]{18,22})(?P<downPrime>CGTACGCTG){e<=1}" \
            -os "barcode" -oi "input.id+\"_\"+sampleIndex.seq" \
            --output-format "fasta" \
    > profiling_tests/long_opOne_outOne.fasta

# Two simpler operations, one output
seq 50 | xargs -I{} cat itermae/data/barseq.fastq \
    | python3 -m cProfile -o profiling_tests/opTwo_outOne.profile \
        bin/itermae \
            -m "input > (?P<sampleIndex>[ATCGN]{5,5})(?P<rest>(GTCCTCGAGGTCTCT){e<=1}[ATCGN]*)" \
            -m "rest  > (?P<upPrime>GTCCTCGAGGTCTCT){e<=1}(?P<barcode>[ATCGN]{18,22})(?P<downPrime>CGTACGCTG){e<=1}" \
            -os "barcode" -oi "input.id+\"_\"+sampleIndex.seq" \
            --output-format "fasta" \
    > profiling_tests/opTwo_outOne.fasta

# Two simpler operations, two outputs
seq 50 | xargs -I{} cat itermae/data/barseq.fastq \
    | python3 -m cProfile -o profiling_tests/opTwo_outTwo.profile \
        bin/itermae \
            -m "input > (?P<sampleIndex>[ATCGN]{5,5})(?P<rest>(GTCCTCGAGGTCTCT){e<=1}[ATCGN]*)" \
            -m "rest  > (?P<upPrime>GTCCTCGAGGTCTCT){e<=1}(?P<barcode>[ATCGN]{18,22})(?P<downPrime>CGTACGCTG){e<=1}" \
            -os "barcode" -oi "input.id+\"_\"+sampleIndex.seq" \
            -os "upPrime+barcode+downPrime" -oi "input.id+\"_withFixedFlanking_\"+sampleIndex.seq" \
            --output-format "fasta" \
    > profiling_tests/opTwo_outTwo.fasta

