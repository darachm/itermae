Tutorial, demo
===================================

This tool is a bit complex to use, so we make recommendations here for how to
effectively test and iteratively tune the tool for your application.
We assume you are probably trying to parse sequence patterns (like barcodes,
sample indicies, tags, UMIs, variant loci) out of many Illumina short-reads,
but a similar procedure should work for other goals and file formats.

The below code snippets are being run from inside the ``itermae`` repo 
directory, as bash commands.
Typological note, I like to use the ``\`` 
to escape/ignore newlines so I can continue with a ``|`` 
character on the next line. 
This style is slightly unconventional, but I find it to be clearer.

1. Install it
-------------------

:doc:`Install the package </install>`. The below is written as if it is a
``pip`` installation - for container-based usage replace ``itermae`` with::

    singularity exec shub://darachm/itermae:latest itermae

in the below examples.

2. Work on a representative subset
--------------------------------------

Generate a representative subset of your sequences to work on. 
I'd recommend no fewer than 100, so that you see common errors. 
You should avoid taking the very first sequences, as those tend to be a bit
erroneous.

.. jupyter-execute::
    :stderr:
    :raises:

    head -n 1000 itermae/data/barseq.fastq \
        | tail -n 400 > test_set.fastq


.. jupyter-execute::
    :stderr:
    :raises:

    head -n 4 test_set.fastq

3. Begin editing a minimal configuration file
-------------------------------------------------

There is an example YAML configuration file available at 
``itermae/data/example_schema.yml`` that you can copy::

    cp itermae/data/example_schema.yml test_config.yml

Or, just write your own.
There's three sections to write, as detailed in 
:ref:`the YAML config docs <yaml-config>`.

- ``input`` - describe where the reads are coming from and what format they're
  in.
- ``matches`` - specify what patterns to match and what groups to save.
- ``output`` - specify what groups to write, where to write it, and in what
  format.

Note that ``input`` is optional, and if omitted we will assume you are 
pipe-ing (with ``|``) a FASTQ file in, so something like ::

    zcat yourFile.fastqz | itermae --config test_config.yml > output.sam

If you do specify an ``input`` block, you can specify what format the reads
are in and if they are gzipped - like so::

    input:
        format: some_gzipped_fastq_file.fastqz
        gzipped: true

For ``matches``, each match is indented with a ``-``. 
You specify what sequence you're using with the ``use:`` key, which is by
default ``input`` (the input sequence). So start with::

    matches:
        - use: input

I recommend you then paste in the sequence of what you expect
your library to look like, then write under it what group each sequence is
part of for your first match. Put in front of these keys 
``pattern:`` and ``marking:`` (picked so that it's the same length as 
``pattern``). For example::

    matches:
        - use: input
          pattern: NNNNNGTCCTCGAGGTCTCTNNNNNNNNNNNNNNNNNNNNCGTACGCTGCAGGTC
          marking: aaaaaBBBBBBBBBBBBBBBccccccccccccccccccccDDDDDDDDDDDDDDD

Here,
``a`` is the first five bases (sample index), 
``B`` is fixed primer sequence,
``c`` is a ~20 base barcode,
and ``D`` is other fixed sequence.

Then below that you want to specify what each of these groups is called,
and any rules about the matching. You can specify the ``name:`` and how long
the letter ``repeat:``'s for or (``repeat_min:`` and ``repeat_max:``).
You can specify error tolerance by specifying how many of any kind of errors 
are allowed (``allowed_errors:``) or particular types of errors 
(``allowed_insertions:``, ``allowed_substitutions:``, ``allowed_deletions:``).
Like so::

    matches:
        - use: input
          pattern: NNNNNGTCCTCGAGGTCTCTNNNNNNNNNNNNNNNNNNNNCGTACGCTGCAGGTC
          marking: aaaaaBBBBBBBBBBBBBBBccccccccccccccccccccDDDDDDDDDDDDDDD
          marked_groups:  
              a:          
                  name: sampleIndex 
                  repeat: 5            
              B:                       
                  name: prefix
                  allowed_errors: 2 
              c:                    
                  name: barcode
                  repeat_min: 18 
                  repeat_max: 22
              D:  
                  allowed_insertions: 1 
                  allowed_deletions: 2
                  allowed_substititions: 2

Note that if you use one of these repeat parameters with a group that is all 
one letter (like a pattern of ``NNNNN``), 
it will collapse that into one character repeated
for as long as you specify (using ``repeat:`` and the like). 
If it's multiple characters (like ``GN``), it will
repeat the whole pattern (like ``GNGNGNGNGN`` if ``repeat: 5``).

For ``output:``, you can specify where to go with ``to:`` and what format with
``format:``. Default ``to:`` is standard output.
Here we'll write it to some FASTQ file::

    output:
        to: some_output_file.fastQ
        format: FASTQ

Below we'll output FASTA to standard output, to demo.

Then we specify a ``list:`` of the different outputs to generate. Here we will
write a first record that is named 'barcode'. It will use the same 'id' field
as the input record, put the sample index sequence in the 'description' field,
and the sequence will just be the 'barcode' matched above::

    output: 
        format: FASTA
        list: 
            -   name: 'barcode'
                id: 'input'
                description: '\"barcode=\"+barcode'
                seq: 'barcode' 

Note that for modifying the ``id:``, ``description:``, or ``seq:``, you've got
to put any plain text in quotes (``"barcode="`` above) and append (``+``) it
to the group sequences you want to append (like ``+barcode``).

Finally, you can set ``verbosity:`` to one of several levels. You can also
set verbosity with the command line argument ``-v``. Command-line directives
are added in after the YAML configuration file is read.

We'll save the total configuration to a file ``test_config.yml``.

.. jupyter-execute::
    :stderr:
    :raises:
    :hide-code:
    :hide-output:

    echo "matches:
        - use: input
          pattern: NNNNNGTCCTCGAGGTCTCTNNNNNNNNNNNNNNNNNNNNCGTACGCTGCAGGTC
          marking: aaaaaBBBBBBBBBBBBBBBccccccccccccccccccccDDDDDDDDDDDDDDD
          marked_groups:
              a:
                  name: sampleIndex
                  repeat: 5
              B:
                  name: prefix
                  allowed_errors: 2
              c:
                  name: barcode
                  repeat_min: 18
                  repeat_max: 22
              D:
                  allowed_insertions: 1
                  allowed_deletions: 2
                  allowed_substititions: 2
    output:
        format: fasta
        list:
            -   name: 'barcode'
                id: 'id'
                description: '\"sample=\"+sampleIndex'
                seq: 'barcode'
    " > test_config.yml

.. jupyter-execute::
    :stderr:
    :raises:
    :hide-code:

    cat test_config.yml

Now we run that, without any verbosity (no ``-v`` on the command-line or 
``verbosity: 1`` in the YAML). 
What do we get? Here we just look at the head of the results.

.. jupyter-execute::
    :stderr:
    :raises:

    cat test_set.fastq | itermae --config test_config.yml | head


4. Seeing more about errors and troubleshooting these
----------------------------------------------------------

Well this is nice that it works, but it would sure be more useful if I showed
you some errors. Here, I'm going to put some errors in the YAML config, and
show you how to see and fix these.

If I know about the error, I have tried to raise a descriptive exception that
explains what to do. If that doesn't make sense then 
`raise an issue at the GitLab repo <https://gitlab.com/darachm/itermae/-/issues>`_.

I am running the below with verbosity set on three by putting ``-vvv`` at the
end. I'm going to limit the inputs to 2 records so that it doesn't output much.

Error in the YAML keys
^^^^^^^^^^^^^^^^^^^^^^^^^^
.. jupyter-execute::
    :stderr:
    :raises:
    :hide-code:
    :hide-output:
 
    echo "matches:
        - use: input
          paddern: NNNNNGTCCTCGAGGTCTCTNNNNNNNNNNNNNNNNNNNNCGTACGCTGCAGGTC
          marking: aaaaaBBBBBBBBBBBBBBBccccccccccccccccccccDDDDDDDDDDDDDDD
          marked_groups:
              a:
                  name: sampleIndex
                  repeat: 5
              B:
                  name: prefix
                  allowed_errors: 2
              #c:
              #    name: barcode
              #    repeat_min: 18
              #    repeat_max: 22
              #D:
              #    allowed_insertions: 1
              #    allowed_deletions: 2
              #    allowed_substititions: 2
    output:
        format: fastaz
        list:
            -   name: 'barcode'
                id: 'input'
                description: 'sample=sampleIndex'
                seq: 'barcode sample'
    " > test_config.yml
 
.. jupyter-execute::
    :stderr:
    :raises:
 
    head -n 8 test_set.fastq | itermae --config test_config.yml -vvv

Ah, I mis-spelled ``pattern:`` as ``paddern:``. This is a silly error, but
that's what it will look like.

Recycling markings
^^^^^^^^^^^^^^^^^^^^^^^^^^
.. jupyter-execute::
    :stderr:
    :raises:
    :hide-code:
    :hide-output:
 
    echo "matches:
        - use: input
          pattern: NNNNNGTCCTCGAGGTCTCTNNNNNNNNNNNNNNNNNNNNCGTACGCTGCAGGTC
          marking: aaaaaBBBBBBBBBBBBBBBccccccccccccccccccccBBBBBBBBBBBBBBB
          marked_groups:
              a:
                  name: sampleIndex
                  repeat: 5
              B:
                  name: prefix
                  allowed_errors: 2
              #c:
              #    name: barcode
              #    repeat_min: 18
              #    repeat_max: 22
              #D:
              #    allowed_insertions: 1
              #    allowed_deletions: 2
              #    allowed_substititions: 2
    output:
        format: fastaz
        list:
            -   name: 'barcode'
                id: 'input'
                description: 'sample=sampleIndex'
                seq: 'barcode sample'
    " > test_config.yml
 
.. jupyter-execute::
    :stderr:
    :raises:
 
    head -n 8 test_set.fastq | itermae --config test_config.yml -vvv

Ah! There is an error in the YAML config::

    marking: aaaaaBBBBBBBBBBBBBBBccccccccccccccccccccBBBBBBBBBBBBBBB

should be::

    marking: aaaaaBBBBBBBBBBBBBBBccccccccccccccccccccDDDDDDDDDDDDDDD

If you want to capture multiple parts as one group, capture them as multiple 
groups and paste them together later.

Missing the ``marked_groups:`` entry for a group
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. jupyter-execute::
    :stderr:
    :raises:
    :hide-code:
    :hide-output:
 
    echo "matches:
        - use: input
          pattern: NNNNNGTCCTCGAGGTCTCTNNNNNNNNNNNNNNNNNNNNCGTACGCTGCAGGTC
          marking: aaaaaBBBBBBBBBBBBBBBccccccccccccccccccccDDDDDDDDDDDDDDD
          marked_groups:
              a:
                  name: sampleIndex
                  repeat: 5
              B:
                  name: prefix
                  allowed_errors: 2
              #c:
              #    name: barcode
              #    repeat_min: 18
              #    repeat_max: 22
              #D:
              #    allowed_insertions: 1
              #    allowed_deletions: 2
              #    allowed_substititions: 2
    output:
        format: fastaz
        list:
            -   name: 'barcode'
                id: 'input'
                description: 'sample=sampleIndex'
                seq: 'barcode sample'
    " > test_config.yml
 
.. jupyter-execute::
    :stderr:
    :raises:
 
    head -n 8 test_set.fastq | itermae --config test_config.yml -vvv

See last line.

Error in syntax of defining output description
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. jupyter-execute::
    :stderr:
    :raises:
    :hide-code:
    :hide-output:
 
    echo "matches:
        - use: input
          pattern: NNNNNGTCCTCGAGGTCTCTNNNNNNNNNNNNNNNNNNNNCGTACGCTGCAGGTC
          marking: aaaaaBBBBBBBBBBBBBBBccccccccccccccccccccDDDDDDDDDDDDDDD
          marked_groups:
              a:
                  name: sampleIndex
                  repeat: 5
              B:
                  name: prefix
                  allowed_errors: 2
              c:
                  name: barcode
                  repeat_min: 18
                  repeat_max: 22
              D:
                  allowed_insertions: 1
                  allowed_deletions: 2
                  allowed_substititions: 2
    output:
        format: fastaz
        list:
            -   name: 'barcode'
                id: 'input'
                description: 'sample=sampleIndex'
                seq: 'barcode_sample'
    " > test_config.yml
 
.. jupyter-execute::
    :stderr:
    :raises:
 
    head -n 8 test_set.fastq | itermae --config test_config.yml -vvv

This line::

    description: 'sample=sampleIndex'

should be::

    description: '"sample="+sampleIndex'

Because the ``sample=`` part is just text pasted on front of the 'sampleIndex'
matched group. Used ``+`` to paste groups and/or text together!

Error in syntax of defining output description
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. jupyter-execute::
    :stderr:
    :raises:
    :hide-code:
    :hide-output:
 
    echo "matches:
        - use: input
          pattern: NNNNNGTCCTCGAGGTCTCTNNNNNNNNNNNNNNNNNNNNCGTACGCTGCAGGTC
          marking: aaaaaBBBBBBBBBBBBBBBccccccccccccccccccccDDDDDDDDDDDDDDD
          marked_groups:
              a:
                  name: sampleIndex
                  repeat: 5
              B:
                  name: prefix
                  allowed_errors: 2
              c:
                  name: barcode
                  repeat_min: 18
                  repeat_max: 22
              D:
                  allowed_insertions: 1
                  allowed_deletions: 2
                  allowed_substititions: 2
    output:
        format: fastaz
        list:
            -   name: 'barcode'
                id: 'input'
                description: '\"sample=\"+sampleIndex'
                seq: 'barcode_sample'
    " > test_config.yml
 
.. jupyter-execute::
    :stderr:
    :raises:
 
    head -n 8 test_set.fastq | itermae --config test_config.yml -vvv

Looks great? Nope! Note that no sequence is output, this is just verbose output.
We see that we start to process each read and attempt to match. The first
read fails to find a match, which is fine because there's not a good match.
But we find a match on the second, but then have "failed to build the output". 
What's wrong?

First problem - the ``seq:`` is set to 'barcode_sample'. Note that we match
a group called 'barcode' and a group called 'sampleIndex', but not 
'barcode_sample'. Instead, let's try ``seq: 'barcode+sampleIndex'`` to paste
them together.

On this one I will hide the verbosity to show the output:

.. jupyter-execute::
    :stderr:
    :raises:
    :hide-code:
    :hide-output:
 
    echo "matches:
        - use: input
          pattern: NNNNNGTCCTCGAGGTCTCTNNNNNNNNNNNNNNNNNNNNCGTACGCTGCAGGTC
          marking: aaaaaBBBBBBBBBBBBBBBccccccccccccccccccccDDDDDDDDDDDDDDD
          marked_groups:
              a:
                  name: sampleIndex
                  repeat: 5
              B:
                  name: prefix
                  allowed_errors: 2
              c:
                  name: barcode
                  repeat_min: 18
                  repeat_max: 22
              D:
                  allowed_insertions: 1
                  allowed_deletions: 2
                  allowed_substititions: 2
    output:
        format: fasta
        list:
            -   name: 'barcode'
                id: 'input'
                description: '\"sample=\"+sampleIndex'
                seq: 'barcode+sampleIndex'
    " > test_config.yml
 
.. jupyter-execute::
    :stderr:
    :raises:
 
    head -n 8 test_set.fastq | itermae --config test_config.yml 

Well the sequence looks like it has the 'sampleIndex' at the end, but...
huh? The ID is the sequence of the input file! That's because I specified::

    id: 'input'

which sets the ID as the 'input' sequence group - the input sequence.
Instead, we can use this field like the 'description' field - this is 
especially useful for passing metadata through formats like SAM.
Here we stick the 'sampleIndex' onto the ID in SAM.

.. jupyter-execute::
    :stderr:
    :raises:
    :hide-code:
    :hide-output:
 
    echo "matches:
        - use: input
          pattern: NNNNNGTCCTCGAGGTCTCTNNNNNNNNNNNNNNNNNNNNCGTACGCTGCAGGTC
          marking: aaaaaBBBBBBBBBBBBBBBccccccccccccccccccccDDDDDDDDDDDDDDD
          marked_groups:
              a:
                  name: sampleIndex
                  repeat: 5
              B:
                  name: prefix
                  allowed_errors: 2
              c:
                  name: barcode
                  repeat_min: 18
                  repeat_max: 22
              D:
                  allowed_insertions: 1
                  allowed_deletions: 2
                  allowed_substititions: 2
    output:
        format: sam
        list:
            -   name: 'barcode'
                id: 'id+\"_sample=\"+sampleIndex'
                description: '\"sample=\"+sampleIndex'
                seq: 'barcode'
    " > test_config.yml
 
.. jupyter-execute::
    :stderr:
    :raises:
 
    head -n 8 test_set.fastq | itermae --config test_config.yml 

That is set by this config


.. jupyter-execute::
    :stderr:
    :raises:
    :hide-code:
 
    cat test_config.yml 

So: 

* generally error messages should be informative, and if not please submit
  a GitLab issue so that I can know and fix that 
* if nothing is being output, there is likely an error in the pattern not
  matching or output not forming. Run with maximal verbosity ``-vvv``, then
  you can compare the outputs and see how far it gets, and where it fails.

5. Parallel-ing with parallel
-------------------------------

I recommend you start by debugging on a small file interactively with a single
command, but most of the time you'll want to actually be eventually be running
this in parallel.

``itermae`` is designed to avoid the complexity and issues of internal
multi-processing, and instead rely on Linux/Unix pipes to have the user do
the parallelization. Chiefly, this is intended to be used with GNU ``parallel``,
since that's a stable tool that performs well and is readily available.

Basically, ``parallel`` will launch multiple ``itermae`` instances, feed chunks
of the file into each job, and collate the output back into one stream
(for each of STDOUT and for STDERR).

I recommend installing GNU ``parallel`` from source, since your package manager
may be out of date or change what program is in the ``parallel`` package 
between releases (I'm looking at you ubuntu). Here's the installation
instructions from the 20210322 release::

    = Full installation =

    Full installation of GNU Parallel is as simple as:

        wget https://ftpmirror.gnu.org/parallel/parallel-20210322.tar.bz2
        wget https://ftpmirror.gnu.org/parallel/parallel-20210322.tar.bz2.sig
        gpg parallel-20210322.tar.bz2.sig
        bzip2 -dc parallel-20210322.tar.bz2 | tar xvf -
        cd parallel-20210322
        ./configure && make && sudo make install

Then start by just getting familiar with running your sequences into 
``parallel``. I recommend using these settings:

* ``--pipe`` pipes the input into each process as STDIN
* ``-l 4`` denotes each record is 4 lines - change this for FASTA, SAM, etc
* ``--keep-order`` maintains the order of input/output 
* ``-N 10000`` denotes how many sequence records to handle per launch of
  ``itermae``, 10000 works - much much larger will consume more RAM and too
  small consumes too much overhead, but it's flexible

Try it out with ``cat``-ing one record out per job. 

.. jupyter-execute::
    :stderr:
    :raises:
 
    echo "=========="
    echo "Compare"
    echo "=========="
    head -n 16 test_set.fastq 
    echo "=========="
    echo "To"
    echo "=========="
    head -n 16 test_set.fastq | parallel --pipe -l 4 --keep-order -N 1 cat

If you run it multiple times without ``--keep-order``, the order should change.
But for bioinformatics you may not want that.

And now we can slot in the ``itermae`` call, here using that config from before.

.. jupyter-execute::
    :stderr:
    :raises:
 
    cat test_set.fastq | parallel --pipe -l 4 --keep-order \
            -N 100 itermae --config test_config.yml \
        | head -n 10

6. Run on the whole file
-------------------------------

Then, just run it on your entire file and save the results.
Here, ``-v`` is useful for run-level configuration and messages:

.. jupyter-execute::
    :stderr:
    :raises:

    cat itermae/data/barseq.fastq \
        | parallel --pipe -l 4 --keep-order -N 1000 \
            itermae --config test_config.yml -v \
        > chopped_outputs.sam
    wc -l chopped_outputs.sam

Of the 1000 input records, in 724 we find out matches and can form the 
desired outputs.
