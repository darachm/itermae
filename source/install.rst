Install
=======

Containers
----------------

We recommend that you don't install ``itermae``.
Seriously, you can pull the Singularity container 
from `here <https://singularity-hub.org/collections/4537>`_
by running::

    singularity exec shub://darachm/itermae:latest itermae

That container contains ``itermae`` along with 
`GNU parallel <https://www.gnu.org/software/parallel/>`_ 
and a few other utilites.

A docker image is forthcoming, it's on the todo list.
The plan is to make an authoritative Docker container, and have the Singularity
container just be built from the Docker.

Actual installation
----------------------

If you'd like to install the package/utility fully, 
you can do that by using pip::

    python3 -m pip install itermae

Remember that you can specify a particular version::

    python3 -m pip install itermae==v0.4.2
