Examples
========


For this page the examples are run using using `jupyter-execute`_.
To do this,

.. _jupyter-execute: a link

.. jupyter-kernel:: bash
    :id: bashy

.. jupyter-execute::
    :stderr:
    :raises:

    itermae 



Parsing an sample index and UMI from the beginning of short-reads
------------------------------------------------------------------


iseq barcodes
input of FASTQ


Extracting a complex UMI from a variable position 
-------------------------------------------------------------------------

sobaseq
complicated yaml grouping


Splitting barcodes from one read into different records
------------------------------------------------------------------

iseq barcodes, but output as different records
multiple outputs

Parallelization - with parallel
--------------------------------------

``itermae`` no longer attempts to run as multiple-processes within one program.
Instead, we keep it as CLI and let parallel do the hard work

examples thereof 

Processing PacBio long reads
--------------------------------------

special considerations about really big/long sam files
