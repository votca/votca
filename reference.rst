Reference
=========

Programs
--------

Mapping file
------------

The root node always has to be cg\_molecule. It can contain the
following keywords:

.. include:: ../csg/mapping.rst

Topology file
-------------

The XMLtopology file

.. include:: ../csg/topol.rst

Settings file
-------------

All options for the iterative script are stored in an xml file.
[sec:ref\_options]

.. include:: ../csg/csg_defaults.rst

Scripts
-------

Scripts are used by and . The script table commonly used (compare
``csg_call –list``):

Script calls can be overwritten by adding a line with the 3rd column
changed to ``csg_table`` in directory.
