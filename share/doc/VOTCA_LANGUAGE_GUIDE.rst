VOTCA Internal Contributor Language Guide
=========================================

This language guide has been created to establish rules to keep VOTCA's
code consistent between repositories. In the past, there has been
difficulty in translating functionality between repositories and within
the same repositories because different properties have been used to
describe the same object attribute. For general programming guidelines 
look at the `Developers Guide <share/doc/DEVELOPERS_GUIDE.rst>`__

Types and Ids
-------------

As an example, consider the csg bead object which had at one time
contained the name, type and id attribute. The name of a bead is ill
defined and could be unique but was not guaranteed to be so.

If a bead were named C5 it was not clear if this was an arbitrary bead
name, or if it was the 5th carbon atom in the system. In any case the
name attribute is not needed because if a unique id is needed the id of
the bead could be used whereas if the type of the bead was needed the
type attribute could be used. As such, the name method and attribute has
been removed from the object.

As a general rule, objects should not have a name method or attribute
rather, any attribute that is not unique to an object should be
indicated with a type method and attribute.

::

    std::string getBeadType();
    std::string getResidueType();

To indicate a unique attribute an id should be used.

::

    int getBeadId();
    int getMoleculeId();

Units in VOTCA
--------------

VOTCA tried as much as possible to standarize units across both CSG and
XTP. Externally, we parse in the units of the respective file format,
e.g. ``.xyz`` ``Angstrom``, ``.gro`` ``nm``. Internally, we convert all
parsed units to:

-  CSG: length ``nm``, energy ``kJ/mol`` and time ``ps``
-  XTP: length ``bohr``, energy ``Hartree`` and time ``ps``

Indexing in VOTCA
-----------------

All indeces in VOTCA start at ``0``. This is useful, because C++ arrays
start at index 0.

Apart from special cases all indices and integers in votca should be
``votca::Index`` which is a typedef for ``long int``. ``.size()``
methods of std::containers return an ``unsigned long int`` and should be
cast to ``votca::Index``. i.e: ``votca::Index(vector.size())``
