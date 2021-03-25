Concept
=======

This tool is designed to robustly parse DNA amplicon-sequencing from designs 
that use complex DNA barcoding architectures. 

Often, positional information is used to extract these. 
This can frequently fail in extracting the wrong sequence due to random or 
systematic deviation from the expectation. 
DNA sequence (especially degenerately synthesized/cloned sequence) is mutable 
both due to technical errors (in PCR and sequencing) or biological mutagenesis.
Both variable lengths of barcodes or technical insertion/deletion errors can 
throw an entire read out of frame, hence regular-expression based methods 
can more robustly extract patterns of sequence from reads. 
However, these are brittle to substitution errors in the "fixed" sequences 
if error-tolerant (ie 'fuzzy') regular expressions are not used. 
Building regular expressions to parse complete amplicon architectures with
complex look-back/forward expressions may be difficult for non-expert users 
to write, debug, and update. 

Therefore a tool that combines chains of multiple simple matches
with flexible filtering expressions can filter out artefactual matches to
efficiently identify expected constructs and yield desired parsed sequence.
The aim of ``itermae`` is to allow a user to flexibly apply 
multiple fuzzy regular-expressions to a sequence-record, in series on the 
same read, and then flexibly filter and build outputs based on what matches. 

For example, here is one FASTQ read with relevant patterns extracted and 
re-arranged for downstream use, all using fuzzy regular expressions:

.. image:: /img/parse_diagram_1.svg

There are three main groups of options to configure:

* Input/Output : what formats to read and write, and where from/to
* Matches : what matches to apply, to what groups to save, and how error-tolerant to be
* Outputs : what criteria on which to filter and what outputs to construct from matched groups

These can be :doc:`configured <usage/config>` with 
:ref:`command-line arguments <cli-config>`
or with a
:ref:`YAML config file <yaml-config>`,
refer to those pages for details.

To begin using, we've outlined a
:doc:`recommended procedure <usage/recommendation>`
for testing configurations.
See :doc:`examples <usage/examples>` for more ideas.

:doc:`Contribution/development <development>` information is available with
current plans and ways to contribute code, ideas, criticism, and bugs.
Also see the :doc:`API <package>` for function-level debugging/development.

