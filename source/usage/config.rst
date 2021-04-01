
.. _configuration:

Configuration
====================================

There's two ways of configuring this tool.
Internally, ``itermae`` first reads the config file you give it
(with ``--config`` argument), then it will apply arguments you supply on
the command-line.
The :doc:`tutorial` is focused on building the YAML config,
and the :doc:`examples` is focused on the command-line interface.
Below is a summary of each:

.. _yaml-config:

Option 1 - a separate YAML file
---------------------------------

I recommend this option, but will maintain the other option because it offers
some direct-to-regex flexibility.

YAML is what
^^^^^^^^^^^^^^^^^

`YAML <https://yaml.org/>`_ is a "human-friendly data serialization standard",
namely a way of structuring a data file that's more readable/editable by us.
For example, a list looks like this::

    - element 1
    - element 2
    - element 3

and dictionaries look like this ::

    key1: 'value 1'
    key2: value 2

and a dictionary in a list in a dictionary looks like::

    top-level dictionary:
        - some list item
        - inner dictionary key 1: value 1
          inner dictionary key 2: value 2
          inner dictionary key 3: value 3
        - another list item

where quotes aren't necessary if it's unambiguous. 
You can test out YAML code 
`in this parser here <https://nodeca.github.io/js-yaml/>`_.

Sections in the config
^^^^^^^^^^^^^^^^^^^^^^^^^^^

``iterame`` expects three top-level keys, like so::

    input:
    matches:
    output:

In ``input:`` you can specify:

* ``from:`` what filepath, if not STDIN, the reads are coming from
* ``format:`` what format, if not FASTQ, the reads are in
* ``gzipped:`` if the input file is gzip compressed

In ``matches:`` you can specify a list (meaning each list item is indented 
and the first line of the list item it's preceeded with a ``-``) 
that specifies:

* ``use:`` which sequence group to use (or the default of 'input')
* ``pattern:`` the sequence pattern to match (in A, T, C, G and IUPAC 
  wildcards such as N - also * means zero or more of anything and + means one
  or more of anything)
* ``marking:`` the corresponding group that each position in the pattern belongs
  to, used in ``marked_groups:`` below
* ``marked_groups:`` for each of the markings used above, several properties
  about the group marked

    * ``name:`` name of the marked group, like 'index' or 'barcode' - this
      allows you to use it as input in later matches, filters, or in 
      output-building
    * ``repeat:`` the number of times to repeat the group - if the group is a
      single character (like ``N``) then the group is just that single character
      repeated that many times, but if it is multiple characters then the entire
      pattern is repeated
    * ``repeat_min:`` to specify a variable range, the minimum
    * ``repeat_max:`` to specify a variable range, the maximum
    * ``allowed_substitutions:`` the maximum number of substitutions allowed 
    * ``allowed_insertions:`` the maximum number of insertions allowed 
    * ``allowed_deletions:`` the maximum number of deletions allowed 
    * ``allowed_errors:`` the maximum number of all error types allowed 

In ``output:`` you can specify:

* ``to:`` what filepath, if not STDOUT, the reads are to be written to
  (keep in mind all messages and debug info goes to standard error)
* ``format:`` what format, if not SAM, to write to (SAM is default because it
  plays nicely with tab-delimited table filtering/splitting/joining)
* ``report:`` an optional filepath, if provided then ``itermae`` will generate
  a per-output report of what was written out, if it was filtered, if it failed,
  and some statistics about the matches 
  (see ``itermae.SeqHolder.format_report``)
* ``failed:`` an optional filepath, if provided then all input reads that fail
  the matches and/or filters will just be printed to this, by default they are
  just forgotten
* ``list:`` a list of outputs to form, for each you can specify

    * ``name:`` name of the output
    * ``filter:`` filter condition to determine if output should be made,
      see :doc:`tutorial` for a more detailed description of how to use this
    * ``id:`` what to put in the ID field of the output sequence record,
      by default this is the input ID ('id'), but you can add on groups or
      text using ``+`` (for example ``id: id+'_someUMIgroup='+matchedUMIgroup``
    * ``description:`` similarly to the ``id:`` field, this is the output field
      of 'description' and is by default the input ('description') 
      but can be added to using ``+`` (for example: 
      ``description: description+' someUMIgroup='+matchedUMIgroup``)
    * ``seq:`` the actual sequence record (and associated quality scores) to
      output, created from sequence groups matched and potentially concatenated
      together with the ``+`` operator (for example: 
      ``seq: sampleIndex+barcode``)
 
Here's an example of a config file::

    input:
        from: STDIN
        format: fastq
        gzipped: false # Default
    matches:
        -   use: input 
            pattern: NNNGTCCACGAGGTCTCTNNNCGTACGCTG 
            marking: AAABBBBBBBBBBBBBBBCCCDDDDDDDDD 
            marked_groups:
                A:
                    name: sampleIndex 
                    repeat: 5            
                B:                       
                    name: prefix
                    allowed_errors: 2 
                C:                    
                    name: barcode
                    repeat_min: 18 
                    repeat_max: 22
                D:  
                    allowed_insertions: 1 
                    allowed_deletions: 2
                    allowed_substititions: 2
        -   use: barcode 
            pattern: N   
            marking: A
            marked_groups:
                A:
                    name: first_five_barcode 
                    repeat: 5
    output: 
        to: STDOUT 
        format: sam 
        list: 
            -   name: barcode 
                filter: 'barcode.length >= 3' 
                id: input 
                description: '"barcode of "+barcode'
                seq: barcode 
                filter: 'statistics.mean(barcode.quality) >= 30'
            -   name: sampleIndex 
                filter: 'sampleIndex.length >= 3'
                seq: sampleIndex
                description: 'description+" is the input description"'
            -   name: demo 
                id:  'id+"_"+sampleIndex'
                seq: 'sampleIndex+dummyspacer+first_five_barcode+dummyspacer+barcode'
        failed: failed 
        report: report.csv 

For clarification, please 
`submit an issue on the GitLab repo <https://gitlab.com/darachm/itermae/-/issues>`_.
I would appreciate the feedback and your help in pointing out problems.

.. _cli-config:

Option 2 - command-line arguments
---------------------------------

You can also just configure ``itermae`` using command-line arugments.
This has the major difference that you have to write full regular-expressions 
in this mode, as opposed to the pattern/marking-groups interface in the
YAML config file.
For additional help with that, see the 
`regex module <https://pypi.org/project/regex/>`_.

Arguments are documented in the command help (see below), but I should make
clear the model of handling matches and outputs.
For both of these modes, I am expected one or more match or output instructions,
so you specify that argument multiple times. As a schematic example::

    itermae -i input_file -m firstMatchRegex -m secondMatchRegex \
        -os firstOutputSequence -os secondOutputSequence

This would apply the ``firstMatchRegex``, then the ``secondMatchRegex``, then
try to generate an output with the ``firstOutputSequence`` and then the
``secondOutputSequence``. 

For output there are multiple options per output, so if you just specify one
``--output-filter`` and multiple ``--output-seq`` options, I will recycle the
same ``--output-filter`` for all of the ``--output-seq``'s.

.. jupyter-kernel:: bash
    :id: bashy

.. jupyter-execute::
    :stderr:
    :raises:

    itermae 

