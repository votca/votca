Reference
=========

Programs
--------

.. include:: ../csg/csg_boltzmann.rst
.. include:: ../csg/csg_call.rst
.. include:: ../csg/csg_defaults.rst
.. include:: ../csg/csg_density.rst
.. include:: ../csg/csg_dlptopol.rst
.. include:: ../csg/csg_dump.rst
.. include:: ../csg/csg_fluctuations.rst
.. include:: ../csg/csg_fmatch.rst
.. include:: ../csg/csg_gmxtopol.rst
.. include:: ../csg/csg_imc_solve.rst
.. include:: ../csg/csg_inverse.rst
.. include:: ../csg/csg_map.rst
.. include:: ../csg/csg_orientcorr.rst
.. include:: ../csg/csg_part_dist.rst
.. include:: ../csg/csg_partial_rdf.rst
.. include:: ../csg/csg_property.rst
.. include:: ../csg/csg_radii.rst
.. include:: ../csg/csg_resample.rst
.. include:: ../csg/csg_reupdate.rst
.. include:: ../csg/csg_sphericalorder.rst
.. include:: ../csg/csg_stat.rst
.. include:: ../csg/csg_traj_force.rst

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
