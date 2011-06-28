#ifndef _DUMP_ATOMS_BJ_H
#define	_DUMP_ATOMS_BJ_H

#include "qmpair.h"
#include "qmcalculator.h"

/**
    \brief Writes xyz coordinates [Angstroem] for all atoms with a given molecule in the origin.

Callname: dumpatoms

 For creating input files used in the TINKER program to compute the self-consistent polarization xyz coordinates [Angstroem] of all atoms are written to files xyz_N, where N (starting at 0) is the number of the molecule of interest whose site energy is to be computed with TINKER. This molecule will be the first entry in xyz_N. Nearest image convention is used after molecule N has been placed in the center of the box. Additionally a cutoff(given in nm with default 50nm) based on centers of mass of molecules can be employed to reduce the number of atoms. The structure of the file xyz_N is: molecule number N(the one to be charged), molecule number K(to which the coordinates xyz of the atoms belong), coordinates X Y Z in Angstreom (TINKER needs Angstroem).

*/
class DumpAtomsBJ : public QMCalculator
{
public:

    DumpAtomsBJ() {};
    ~DumpAtomsBJ() {};
    const char *Description() { return "Write three files per pair for integrals in presence of surrouding partial charges"; }

    void Initialize(QMTopology *top, Property *options);
    bool EvaluateFrame(QMTopology *top);
    void WriteAtoms(Topology *atop, Molecule *mol);
   
private:
    double _dump_cutoff;
    Property * _options;
};

#endif	/* _DUMP_ATOMS_BJ_H */

