.. _examples:

Examples
========

The examples are run using using 
`jupyter-sphinx <https://jupyter-sphinx.readthedocs.io/en/latest/>`_, 
so should be reproducible for you if you install ``itermae`` 
and run each one or all of them (with ``make docs``) 
from inside the package repo.

All the below examples are using the 
:ref:`command-line configuration <cli-config>`.
This requires some familiarity with regular expressions, and the ``regex``
module syntax for error tolerance (fuzziness).
In the future, I hope to also write a 
:ref:`YAML config file <yaml-config>` 
for each, but for now the :doc:`tutorial` focuses on that.

If you encounter any errors, turn on verbose modes with ``-v`` up to ``-vvv``.

.. jupyter-kernel:: bash
    :id: bashy

Boring minimal example
-------------------------

This just reads in the first two records (8 lines in a FASTQ) 
and outputs them again.

.. jupyter-execute::
    :stderr:
    :raises:
 
    head -n 8 itermae/data/tests/test_inputs/barseq.fastq \
        | itermae -m "input > ." -os "input"

We could use this as a format converter, by using ``--output-format``.

.. jupyter-execute::
    :stderr:
    :raises:
 
    head -n 8 itermae/data/tests/test_inputs/barseq.fastq \
        | itermae -m "input > ." -os "input" --output-format fasta

And we can read different input formats by specifying ``--input-format``.

.. jupyter-execute::
    :stderr:
    :raises:
 
    head -n 4 itermae/data/tests/test_inputs/barseq.fasta \
        | itermae -m "input > ." -os "input" \
            --input-format fasta --output-format fastq

Parsing an sample index and UMI from the beginning of short-reads
------------------------------------------------------------------

For our first real example, we'll parse 
`iSeq2.0 <https://doi.org/10.1016/j.cels.2019.03.005>`_ reads.
In this design, here on the forward read, the first eight bases are degenerate
primer-incorporated UMIs. The next six bases are index, then we have a fixed
priming sequence of ``TTAATATGGACTAAAGGAGGCTTTT``.

Here we'll sort through ten reads from that file, pick out the UMI and sample
index code, and print out the priming site as the main sequence:

.. jupyter-execute::
    :stderr:
    :raises:
 
    cat itermae/data/iseq2_example.fastq \
        | itermae \
            -m "input > ^(?P<umi>[ATCGN]{6,9})(?P<sample>[ATCGN]{6,6})(?<fixed>TTAATATGGACTAAAGGAGGCTTTT){e<=3}" \
            -os "fixed" \
            -od "description+' umi='+umi+' sample='+sample" \
            --output-format fasta

That's not very useful, so to extract the barcode here I'm going to split this
up into two matches. First, we use the first fixed sequence to split the
read into UMI, sample code, and the rest. Then we used the sequence just
outside the barcode to find the barcode. And we also delete the illumina
multiplexing tags from the description, to clean that up a bit:

.. jupyter-execute::
    :stderr:
    :raises:
 
    cat itermae/data/iseq2_example.fastq \
        | itermae \
            -m "input > ^(?P<umi>[ATCGN]{6,9})(?P<sample>[ATCGN]{6,6})(?<fixed1>(TTAATATGGACTAAAGGAGGCTTTT){e<=3})(?<rest>[ATCGN]*$)" \
            -m "rest > (?P<fixed_up>TATCGGTACC){e<=1}(?<barcode>[ATCGN]{0,40})(?P<fixed_down>GATAACTTCG){e<=1}" \
            -os "barcode" \
            -od "'umi='+umi+' sample='+sample" \
            --output-format fasta

And that's ready for downstream processing (clustering).

Extracting a complex UMI from a variable position 
-------------------------------------------------------------------------

``itermae`` was written originally (as ``SLAPCHOP``, using sequence alignments)
for this use case. 
This is barseq of the 
`yeast deletion collection tags <https://doi.org/10.1126/science.285.5429.901>`_,
but incorporating a sample-multiplexing tag in the first five bases and a UMI
in the reverse primer.
The collection, as any biological entity, is mutable and it has been
`re-annotated <https://doi.org/10.1101/gr.093955.109>`_
to characterize that the 20-mer barcodes are actually now of variable length
(10-26) and at least one has very similar sequence to the reverse priming site
(potentially a deletion of the barcode?).
Additionally, the use of a UMI in a low-input sample inspired the use of
UMIs with semi-fixed positions (to prevent priming off of similar random UMIs).

.. image:: /img/parse_diagram_1.svg

The below command parses that:

.. jupyter-execute::
    :stderr:
    :raises:
 
    head -n 40 itermae/data/tests/test_inputs/barseq.fastq \
        | itermae \
            -m "input > (?P<sample>[ATCG]{5})(?P<fixed1>GTCCACGAGGTC){e<=2}(?P<rest>TCT.*){e<=1}" \
            -m "rest > (?P<tag>TCT){e<=1}(?P<strain>[ATCG]{10,26})(CGTACGCTGC){e<=2}" \
            -m "rest > (?P<fixed2>CGTACGCTGCAGGTC)(?<UMItail>GAC[ATCG]G[ATCG]A[ATCG]G[ATCG]G[ATCG]G[ATCG]GAT){s<=2}" \
            -os "strain" \
            -od "'sample='+sample+' umi='+UMItail" \
            --output-format fasta

Splitting barcodes from one read into different records
------------------------------------------------------------------

Different barcodes in an amplicon design likely require different 
post-processing, such as clustering-based error correction.
Parameters for this are tuned for barcodes of different library and 
barcode-space complexity, and so ``itermae`` takes advantage of SAM tags to
mark different outputs to permit forking one parsed file into multiple
downstream clustering workflows.

For example:

.. jupyter-execute::
    :stderr:
    :raises:
 
    cat itermae/data/iseq2_example.fastq \
        | itermae \
            -m "input > ^(?P<umi>[ATCGN]{6,9})(?P<sample>[ATCGN]{6,6})(?<fixed1>(TTAATATGGACTAAAGGAGGCTTTT){e<=3})(?<rest>[ATCGN]*$)" \
            -m "rest > (?P<fixed_up>TATCGGTACC){e<=1}(?<barcode>[ATCGN]{0,40})(?P<fixed_down>GATAACTTCG){e<=1}" \
            -os "sample" \
            -os "barcode" \
            --output-format sam

Note that there are alternating lines of the sample barcode (~6-base)
and strain barcode (~26-base), with tags of ``IE:Z:untitled_output_0``
and ``IE:Z:untitled_output_1``. This would permit splitting these with 
something like ``... | grep "IE:Z:untitled_output_0" | ...`` in the pipeline.
( Note: the YAML API, as detailed in :doc:`tutorial`, permits naming outputs. )


Parallelization - with parallel
--------------------------------------

``itermae`` originally attempted to launch and manage multiple-processes within
one launcher program. Inspired by a memory-leak (from escalation of variables
to global from within the ``regex`` module), I decided to focus ``itermae``
as a pipe-in pipe-out do-one-thing-well command-line tool.

Instead, I now let `GNU parallel <https://www.gnu.org/software/parallel/>`_ 
do the hard work. This can be a little strange to write, but is made much
easier with the `YAML config <yaml-config>`_ interface.
One could also use other strategies, like splitting files with ``split``, but
I have found this one to be stable, well-supported, and performant.

For example, here I parallelize the above. It feeds one chunk at a time 
(``-N 1``), where chunks are 4 lines (``-l 4``), pipes it in to itermae 
(``--pipe``), keeps the order of input-output (``--keep-order``), and
uses ``--quote`` to protect all the funny regex characters:

.. jupyter-execute::
    :stderr:
    :raises:
 
    cat itermae/data/iseq2_example.fastq \
        | parallel --pipe -l 4 --keep-order -N 1 --quote \
            itermae \
            -m "input > ^(?P<umi>[ATCGN]{6,9})(?P<sample>[ATCGN]{6,6})(?<fixed1>(TTAATATGGACTAAAGGAGGCTTTT){e<=3})(?<rest>[ATCGN]*$)" \
            -m "rest > (?P<fixed_up>TATCGGTACC){e<=1}(?<barcode>[ATCGN]{0,40})(?P<fixed_down>GATAACTTCG){e<=1}" \
            -os "sample" \
            -os "barcode" \
            --output-format sam

For actual large-runs, I recommend setting ``-N 100000``, such that good sized
chunks are run per the overhead of each ``itermae`` configuration and setup 
stage. The number of jobs run defaults to run one job per CPU, but can be 
regulated with a ``-j 4`` option.

GNU parallel has `extensive documentation and tutorials <https://www.gnu.org/software/parallel/index.html#Tutorial>`_.

.. Should add examples with special considerations about really big/long sam 
   files, ie PacBio data, but I'm still working that out!
