==============================
API
==============================

``itermae`` is intended to be used as a command line utility, but here's 
list of the internal functions and classes to orient y'all for debugging (and 
contributing to development?).

-----------------------
Main functions
-----------------------
.. autofunction:: itermae.reader
.. autofunction:: itermae.chop

-----------------------
Classes
-----------------------
.. autofunction:: itermae.MatchScores
.. autofunction:: itermae.GroupStats
.. autofunction:: itermae.SeqHolder

-----------------------
Configuration utilities
-----------------------

.. autofunction:: itermae.config_from_file
.. autofunction:: itermae.config_from_args
.. autofunction:: itermae.check_reserved_name

-----------------------
File format utilities
-----------------------
.. autofunction:: itermae.open_appropriate_input_format
.. autofunction:: itermae.read_sam_file
.. autofunction:: itermae.read_txt_file
.. autofunction:: itermae.write_out_seq
.. autofunction:: itermae.format_sam_record

----------------------
Misc utilities
----------------------

.. autofunction:: itermae.open_input_fh
.. autofunction:: itermae.open_output_fh
.. autofunction:: itermae.fix_desc
.. autofunction:: itermae.phred_letter_to_number
.. autofunction:: itermae.phred_number_to_letter
.. autofunction:: itermae.phred_number_array_to_joined_string



