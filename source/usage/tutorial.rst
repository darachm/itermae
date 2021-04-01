Tutorial / demos
===================================

Here we make recommendations here for how to
effectively test and iteratively tune the tool for your application.
We assume you are probably trying to parse sequence patterns (like barcodes,
sample indicies, tags, UMIs, variant loci) out of many Illumina short-reads,
but a similar procedure should work for other goals and file formats.

The below focuses on using the YAML config file.
For quick examples of using command-line arguments, 
see :doc:`examples`.

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

in the below examples. This is especially useful when using ``parallel`` below,
as that is included in the ``itermae`` container.

2. Work on a representative subset
--------------------------------------

Generate a representative subset of your sequences to work on. 
I'd recommend no fewer than 100, so that you see common errors. 
You should avoid taking the very first sequences of a real Illumina FASTQ file, 
as those tend to be a bit erroneous.

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

    cat yourFile.fastq | itermae --config test_config.yml > output.sam

or for a gzip'd FASTQ (FASTQZ) file::

    zcat yourFile.fastqz | itermae --config test_config.yml > output.sam

If you do specify an ``input`` block in the YAML config file, 
you can specify what format the reads are in and if they are gzipped - like so::

    input:
        format: some_gzipped_fastq_file.fastqz
        gzipped: true

Then for the ``matches`` section, each match is indented with a ``-``
to denote that it is an item in a list. 
You specify what sequence you're using with the ``use:`` key, which is by
default ``input`` (the input sequence). So start with::

    matches:
        - use: input

I recommend you then paste in the sequence of what you expect
your library to look like, then write under it what group each sequence is
part of for your first match. Put in front of these keys 
``pattern:`` and ``marking:`` (picked so that it's the same length as 
``pattern``). Note that they need to be indented from the ``-`` a little.
For example::

    matches:
        - use: input
          pattern: NNNNNGTCCTCGAGGTCTCTNNNNNNNNNNNNNNNNNNNNCGTACGCTGCAGGTC
          marking: aaaaaBBBBBBBBBBBBBBBccccccccccccccccccccDDDDDDDDDDDDDDD

Here,
``a`` is the first five bases (sample index), 
``B`` is fixed primer sequence,
``c`` is a ~20 base barcode,
and ``D`` is other fixed sequence.

Then below that you want to specify in a section called ``marked_groups:``
what each of these groups is called,
and any rules about the matching. 
Note that this is indented in from the previous section.
You can specify the ``name:`` and how long
the pattern ``repeat:``'s for or (``repeat_min:`` and ``repeat_max:``).
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

For ``output:``, you can specify a file path to write to with ``to:`` 
(or just leave it on the default of standard output) and what format with
``format:``.  Here we'll write it to some FASTQ file path::

    output:
        to: some_output_file.fastq
        format: FASTQ

Below we'll output FASTA to standard output, to demonstrate.

Then we specify a ``list:`` of the different outputs to generate. Here we will
write a first record that is named 'barcode'. It will use the same 'id' field
as the input record, put the sample index sequence in the 'description' field,
and the sequence will just be the 'barcode' matched above::

    output: 
        format: FASTA
        list: 
            -   name: 'barcode'
                id: 'input'
                description: 'description+" sample="+sampleIndex'
                seq: 'barcode' 

Note that for modifying the ``id:``, ``description:``, or ``seq:``, you've got
to put any plain text in quotes (``" sample="`` above) and append (``+``) it
to the group sequences you want to append (like ``+sampleIndex``). 
``description`` contains the original description, so the above is appending 
``" sample="+sampleIndex``
on to that.

Finally, you can set ``verbosity:`` to one of several levels. 
Rather, I would recommend that you use the command line argument ``-v``,
as this is more readable and changeable in debugging. 
Command-line directives
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
                description: 'description+\" sample=\"+sampleIndex'
                seq: 'barcode'
    " > test_config.yml

.. jupyter-execute::
    :stderr:
    :raises:
    :hide-code:

    cat test_config.yml

Now we run that, without any verbosity (no ``-v`` on the command-line or 
``verbosity:`` in the YAML). 
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

If I know about the type of error, 
I have tried to raise a descriptive exception that
explains what to do. If that doesn't make sense then please
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
                description: 'description+\" sample=\"+sampleIndex'
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
                description: 'description+\" sample=\"+sampleIndex'
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
                description: 'description+\" sample=\"+sampleIndex'
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
                description: 'description+sample=+sampleIndex'
                seq: 'barcode_sample'
    " > test_config.yml
 
.. jupyter-execute::
    :stderr:
    :raises:
 
    head -n 8 test_set.fastq | itermae --config test_config.yml -vvv

This line::

    description: 'description+sample=sampleIndex'

should be::

    description: 'description+" sample="+sampleIndex'

Because the ``sample=`` part is just text pasted inbetween the 'description'
and the 'sampleIndex' matched group. 
Use ``+`` to paste groups and/or quoted text together!

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
                description: 'description+\" sample=\"+sampleIndex'
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
                description: 'description+\" sample=\"+sampleIndex'
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
                description: 'description+\" sample=\"+sampleIndex'
                seq: 'barcode'
    " > test_config.yml
 
.. jupyter-execute::
    :stderr:
    :raises:
 
    head -n 8 test_set.fastq | itermae --config test_config.yml 

So: 

* generally error messages should be informative, and if not please submit
  a GitLab issue so that I can know and fix that 
* if nothing is being output, there is likely an error in the pattern not
  matching or output not forming. Run with maximal verbosity ``-vvv``, then
  you can compare the outputs and see how far it gets, and where it fails.


5. Filtering
-------------------------------

You may very well want to filter the reads on a variety of properties.
You do this by adding a ``filter:`` to an output in the output list that
will output that sequence if it passes the filter. Such as::

    output:
        format: fasta
        list:
            -   name: 'barcode'
                id: 'input'
                description: 'description+\" sample=\"+sampleIndex'
                seq: 'barcode+sampleIndex'
                filter: 'barcode.length >= 20'

The filter is evaluated as python, so you can use things like ``>=`` or
``and`` or ``or`` to combine multiple statements in the filter.
There are several of such properties available internally per matched group, 
such as:

* ``some_group.start`` - specifies where in the read ``some_group`` starts
* ``some_group.end`` - specifies where in the read ``some_group`` ends
* ``some_group.length`` - specifies the length of ``some_group``
* ``some_group.quality`` - stores a numeric array of the PHRED qualities 
    associated with the sequence in ``some_group``

For that last one especially, the module ``statistics`` is loaded so that you
may make use of expressions such as 
``statistics.median(some_group.quality) >= 30``. 
See the `statistics <https://docs.python.org/3/library/statistics.html>`_ module
for more functions, like 
``statistics.mean(some_group.quality) >= 30`` or 
``statistics.geometric_mean(some_group.quality) >= 30``.

There are match-level properties too. Each match is named ``match_0`` or
``match_1`` etc in the order that it is specified (in YAML or command line),
so these properties can also be used in a filter:

* ``some_match.substitutions`` - stores how many substitutions were necessary 
  for a match
* ``some_match.insertions`` - stores how many insertions were necessary 
  for a match
* ``some_match.deletions`` - stores how many deletions were necessary 
  for a match

I do not anticipate these to be readily useful, but they are available in case
you envision some useful edge case, like 
``( statistics.mean(some_group.quality) >= 30 and match_0.substitions == 0 ) or ( statistics.mean(some_group.quality) <= 30 and match_0.substitions >= 0) )`` 
... it's there if you need it.

Adding that into the previous configuration, here the file we have built up:

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
                description: 'description+\" sample=\"+sampleIndex'
                seq: 'barcode'
                filter: 'statistics.mean(barcode.quality) >= 30'
    " > test_config.yml
 
.. jupyter-execute::
    :stderr:
    :raises:
    :hide-code:
 
    cat test_config.yml 


6. Parallel-ing with parallel
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

* ``--quote`` protects any funny regex characters from being interpreted as BASH
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
    head -n 16 test_set.fastq \
        | parallel --quote --pipe -l 4 --keep-order -N 1 cat

If you run it multiple times without ``--keep-order``, the order should change.
But for bioinformatics you may not want that.

And now we can slot in the ``itermae`` call, here using that config from before,
and restricting it to just the first ten outputs:

.. jupyter-execute::
    :stderr:
    :raises:
 
    cat test_set.fastq | parallel --quote --pipe -l 4 --keep-order \
            -N 100 itermae --config test_config.yml \
        | head -n 10

7. Run on the whole file
-------------------------------

Then, just run it on your entire file and save the results.
Here, ``-v`` is useful for run-level configuration and messages:

.. jupyter-execute::
    :stderr:
    :raises:
 
    cat itermae/data/barseq.fastq \
        | parallel --quote --pipe -l 4 --keep-order -N 1000 \
            itermae --config test_config.yml -v \
        > chopped_outputs.sam


Of the 1000 input records, we find this many that have matches and can form
the desired outputs:

.. jupyter-execute::
    :stderr:
    :raises:
 
    wc -l chopped_outputs.sam
