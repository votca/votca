#ifndef MOL_PAIR
#define MOL_PAIR

#include <string>
#include "qm_molecule.h"
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <string.h>

class mol_pair {
private:

    mol_and_orb *_morb1, *_morb2;

public:

    mol_pair(mol_and_orb* morb1, mol_and_orb* morb2);

    ~mol_pair();

    void write_conf(const char *file_type);
};

// int strcmp ( const char * str1, const char * str2 );

inline mol_pair::mol_pair(mol_and_orb* morb1, mol_and_orb* morb2)
{
    _morb1 = morb1;
    _morb2 = morb2;
}

inline mol_pair::~mol_pair()
{
    delete _morb1; // do we need this?
    delete _morb2;
}

#endif
