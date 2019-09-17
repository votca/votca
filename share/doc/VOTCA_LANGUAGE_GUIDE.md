# VOTCA Internal Contributor Language Guide

This language guide has been created to establish rules to keep VOTCAs code
consistent between repositories. In the past, there has been difficulty in
translating functionality between repositories and within the same repositories
because different properties have been used to describe the same object
attribute.

## Types and Ids

As an example consider the csg bead object which had at one time contained the
name, type and id attribute. The name of a bead is ill defined and could be
unique but was not guaranteed to be so. 

If a bead were named C5 it was not clear if this was an arbitrary bead name, or
if it was the 5th carbon atom in the system. Any any case the name attribute is
not needed because if a unique id is needed the id of the bead could be used
whereas if the type of the bead was needed the type attribute could be used.
As such, the name method and attribute has been removed from the object. 

As a general rule, objects should not have a name method or attribute rather,
any attribute that is not unique to an object should be indicated with a type
method and attribute. 

    std::string getBeadType();
    std::string getResidueType();

To indicate a unique attribute an id should be used. 

    int getBeadId();
    int getMoleculeId();
