Reference
=========

Programs
--------

csg_boltzmann
^^^^^^^^^^^^^

.. include:: csg_boltzmann.rst

csg_call
^^^^^^^^

.. include:: csg_call.rst

csg_density
^^^^^^^^^^^
.. include:: csg_density.rst

csg_dlptopol
^^^^^^^^^^^^
.. include:: csg_dlptopol.rst

csg_dump
^^^^^^^^
.. include:: csg_dump.rst

csg_fmatch
^^^^^^^^^^
.. include:: csg_fmatch.rst

csg_gmxtopol
^^^^^^^^^^^^
.. include:: csg_gmxtopol.rst

csg_imc_solve
^^^^^^^^^^^^^

.. include:: csg_imc_solve.rst

csg_inverse
^^^^^^^^^^^

.. include:: csg_inverse.rst

csg_map
^^^^^^^

.. include:: csg_map.rst

csg_property
^^^^^^^^^^^^
.. include:: csg_property.rst

csg_resample
^^^^^^^^^^^^

.. include:: csg_resample.rst

csg_reupdate
^^^^^^^^^^^^

.. include:: csg_reupdate.rst

csg_stat
^^^^^^^^

.. include:: csg_stat.rst


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

Script calls can be overwritten by adding a similar line with the 3rd column
changed from the original ``csg_table`` in directory list in the
``<scriptdir>`` property of the settings xml file.

List of the default scripts:

.. include:: csg_table.rst

csgapps
-------

Extra user-contibuted applications based on the csg application class, that aren't needed for the normal workflows. Installable with ``-DBUILD_CSGAPPS=ON``.

csg_fluctuations
^^^^^^^^^^^^^^^^

.. include:: csg_fluctuations.rst

csg_orientcorr
^^^^^^^^^^^^^^

.. include:: csg_orientcorr.rst

csg_part_dist
^^^^^^^^^^^^^

.. include:: csg_part_dist.rst

csg_partial_rdf
^^^^^^^^^^^^^^^

.. include:: csg_partial_rdf.rst

csg_radii
^^^^^^^^^

.. include:: csg_radii.rst

csg_sphericalorder
^^^^^^^^^^^^^^^^^^

.. include:: csg_sphericalorder.rst

csg_traj_force
^^^^^^^^^^^^^^

.. include:: csg_traj_force.rst
