Recommended procedure
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
are allowed (``allowed_errors``) or particular types of errors 
(``allowed_insertions``, ``allowed_substitutions``, ``allowed_deletions``).
::

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
for as long as you specify. If it's multiple characters (like ``GN``), it will
repeat the whole pattern (like ``GNGNGNGNGN`` if ``repeat: 5``).

For ``output:``, you can specify where to go with ``to:`` and what format with
``format:``. Default ``to:`` is standard output, we'll change the format to
FASTA::

    ouptut:
        format: FASTA

Then we specify a ``list:`` of the different outputs to generate. Here we will
write a first record that is named 'barcode'. It will use the same 'id' field
as the input record, put the sample index sequence in the 'description' field,
and the sequence will just be the 'barcode' matched above::

    output: 
        format: FASTA
        list: 
            -   name: barcode 
                id: input 
                description: '\"barcode=\"+barcode'
                seq: barcode 

We'll save the above to a file ``test_config.yml``.

.. jupyter-execute::
    :stderr:
    :raises:
    :hide-code:
    :hide-output:

    echo "verbosity: 1 
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
    output: 
        format: fasta
        list:
            -   name: barcode 
                id: input 
                description: 'barcode'
                seq: barcode 
        report: report
    #            -   name: barcode 
    #                filter: 'barcode.length >= 3' 
    #                id: input 
    #                description: '\"barcode of \"+barcode'
    #                seq: barcode 
    #            -   name: sampleIndex 
    #                filter: 'sampleIndex.length >= 3'
    #                seq: sampleIndex
    #                description: 'description+\" is the input description\"'
    #            -   name: demo 
    #                id:  'id+\"_\"+sampleIndex'
    #                seq: 'sampleIndex+dummyspacer+first_five_barcode+dummyspacer+barcode'
    " > test_config.yml

.. jupyter-execute::
    :stderr:
    :raises:
    :hide-code:

    cat test_config.yml


Now we run that, finally

.. jupyter-execute::
    :stderr:
    :raises:

    cat test_set.fastq | itermae --config test_config.yml 




er::

        -   use: barcode 
            pattern: N   
            marking: A
            marked_groups:
                A:
                    name: first_five_bases_of_the_barcode 
                    repeat: 5

GACNGNANGNGNGNGATGTG
