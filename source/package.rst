==============================
API - class/function docs
==============================

``itermae`` is intended to be used as a command line utility, but here's 
list of the internal functions and classes to orient y'all for debugging (and 
contributing to development?).

Essentially, the ``bin/itermae`` launcher CLI script reads arguments in,
creates a ``Configuration`` class, tells it to configure with certain arguments,
then tells it to start reading with method ``.reader()``.
Internally, it creates a ``SeqHolder`` object for each input read, 
which handles all the sequence intermediates and outputing.
There are a few other little utility functions/modules.

The below is automatically generated from the function-level docstrings:

.. automodule:: itermae
    :members:
