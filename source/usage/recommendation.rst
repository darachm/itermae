Recommended procedure
===================================

This tool is a bit complex to use, so we make recommendations here for how to
effectively test and iteratively tune the tool for your application.
These are recommendations, of course.

We assume you are probably trying to parse sequence patterns (like barcodes,
sample indicies, tags, UMIs, variant loci) out of many Illumina short-reads,
but a similar procedure should work for other technologies (so long as it's in
FASTQ, SAM, FASTA, or txt file format).

1. Generate a representative subset of your sequences to work on. 
I'd recommend no fewer than 100, so that you see common errors. 
You should avoid taking the very first sequences, as those tend to be a bit
erroneous.

.. jupyter-execute::
    :stderr:
    :raises:

    head -n 1000 itermae/data/barseq.fastq \
        | tail -n 400 > test_set.fastq

Note that I like to use the ``\`` to escape/ignore the newline and continue
with a ``|`` character on the next line. 
Just a personal style, I find it to be clearer.

.. jupyter-execute::
    :stderr:
    :raises:

    wc -l test_set.fastq

2. Copy the example YAML configuration file and modify it to make a minimal
matching pattern. This is located at ``itermae/data/example_schema.yaml``.


