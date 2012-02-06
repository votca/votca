/* 
 * File:   EMultipole.cc
 * Author: poelking
 * 
 * Created on December 27, 2011, 4:10 PM
 */

#include "emultipole.h"
#include <vector>

namespace votca { namespace ctp {

/**
 * \brief Read in parameters for Thole model and convergence loop
 * @param top
 * @param opt
 */
void EMultipole::Initialize(QMTopology *top, Property *opt) {

    string key;
    string xmlfile;

    /* ---- OPTIONS.XML Structure ----
     * <emultipole>
     *
     *      <multipoles></multipoles>
     *
     *      <tholeparam>
     *          <cutoff></cutoff>
     *          <expdamp></expdamp>
     *          <scaling></scaling>
     *      </tholeparam>
     *
     *      <convparam>
     *          <omegSOR></omegSOR>
     *          <maxiter></maxiter>
     *          <tolerance></tolerance>
     *      </convparam>
     */

    key = "options.emultipole.multipoles";

        if ( opt->exists(key) ) {
            xmlfile = opt->get(key).as< string >();
        }

    key = "options.emultipole.tholeparam";

        if ( opt->exists(key+".cutoff") ) {
            _cutoff = opt->get(key+".cutoff").as< double >();
            if (_cutoff) { _useCutoff = true; }
        }
        if ( opt->exists(key+".expdamp") ) {
            _aDamp = opt->get(key+".expdamp").as< double >();
            if (_aDamp) { _useExp = true; }
        }
         if ( opt->exists(key+".scaling") ) {
            _scale1 = opt->get(key+".scaling").as< vector<double> >();
            if (0 < _scale1.size() && _scale1.size() < 4) {
                _useScaling = true; }
            else {
                _useScaling = false;
                cout << "WARNING: 1-N SCALING SWITCHED OFF" << endl; }
        }

    key = "options.emultipole.convparam";

        if ( opt->exists(key+".omegSOR") ) {
            _omegSOR = opt->get(key+".omegSOR").as< float >();
        }
        else { _omegSOR = 0.75; }

        if ( opt->exists(key+".maxiter") ) {
            _maxIter = opt->get(key+".maxiter").as< int >();
        }
        else { _maxIter = 512; }

        if ( opt->exists(key+".tolerance") ) {
            _epsTol = opt->get(key+".tolerance").as< double >();
        }
        else { _epsTol = 0.001; } 

    // Import charges from XML source
    this->GetMPoles(xmlfile);

    // Equip topology with charges
    this->EquipTop(top);

    // Restate input data to screen
    this->PrintInfo("Thole"); this->PrintInfo("MPoles");

}

/**
 * \brief Read in multipole data for atom types from XML file.
 * @param xmlfile
 */
void EMultipole::GetMPoles(string &xmlfile) {
     cout << "Import multipoles from " << xmlfile << endl;

     // TODO Replace rsdno with rsdname; also affects ::EquipTop
     // TODO and ::ChargeMol; teach topology about  residue names

     // map to store polarizability tensors
     //     string -> atom identifier
     //     matrix -> polarizability tensor
     map < string, matrix >           polTs;
     map < string, matrix >::iterator polit;

     // map to store charges [ anion, neutral, cation ]
     //     string -> atom identifier
     //     int    -> charge state (-1, 1, 1)
     //     double -> charge
     map < string, map<int, double> >           chrgs;
     map < string, map<int, double> >::iterator chrgit;


     // state sfor which data is available;
     // neutral assumed as precondition
     bool anion;
     bool cation;

    /* ---- MULTIPOLES.XML Structure ----
     *  <molecule>
     *      <molname></molname>
     * 
     *      <residue>   
     *          <rsdno></rsdno>
     *   
     *          <atom>
     *              <atmname></atmname>
     *              <chrgA></chrgA>
     *              <chrgN></chrgN>
     *              <chrgC></chrgC>
     *              <polTensor></polTensor>
     *          </atom>  
     *     
     *      </residue>     
     *  </molecule>
     */

    // Load options from specified xml file
    Property opt;
    load_property_from_xml( opt, xmlfile.c_str() );

    // Loop through molecules
    list< Property * > mols = opt.Select("molecules.molecule");
    list< Property * >::iterator mit;
    for ( mit = mols.begin();
          mit != mols.end();
          mit++ ) {

          string molname = (*mit)->get("molname").as< string >();

          anion = true;
          cation = true;

          // Loop through residues
          list< Property * > rsds = (*mit)->Select("residue");
          list< Property * >::iterator rit;
          for ( rit = rsds.begin();
                rit != rsds.end();
                rit++ ) {

                string rsdno = (*rit)->get("rsdno").as< string >();

                // Loop through atoms
                list< Property * > atms = (*rit)->Select("atom");
                list< Property * >::iterator ait;
                for ( ait = atms.begin();
                      ait != atms.end();
                      ait++ ) {

                      string atmname = (*ait)->get("atmname").as< string >();

                      // Putting atmname first should speed up look-up process
                      // later on in ::Induce, ::...
                      string namekey = atmname+"_"+rsdno+"_"+molname;

                      // Check whether data already stored. If so, continue.
                      polit = polTs.find(namekey);
                      chrgit = chrgs.find(namekey);
                      if ( polit != polTs.end() &&
                           chrgit != chrgs.end() ) { continue; }

                      // Retrieve charges (anion, neutral, cation)                      
                      map< int, double > chrg;
                      double chrgA;
                      double chrgN;
                      double chrgC;
                      
                      try { chrgA = (*ait)->get("chrgA").as< double >();}
                      catch ( std::runtime_error ) { anion = false; }

                      chrgN = (*ait)->get("chrgN").as< double >();

                      try { chrgC = (*ait)->get("chrgC").as< double >();}
                      catch ( std::runtime_error ) { cation = false; }

                      if (anion) { chrg[-1] = chrgA; }
                      chrg[0] = chrgN;
                      if (cation) { chrg[1] = chrgC; }

                      // Retrieve dipole polarizability tensor
                      vector<double> polV =
                         (*ait)->get("polTensor").as< vector<double> >();

                      if ( polV.size() != 9 ) {
                          cout << "Invalid input for polarizability tensor of "
                               << "atom " << atmname
                               << "in residue " << rsdno
                               << "in molecule " << molname
                               << "." << endl;
                          cout << "Format is: xx yx zx xy yy zy xz yz zz"
                               << "." << endl;
                          throw std::runtime_error("Redo tensor input.");
                      }

                      matrix polT = matrix( vec(polV[0], polV[3], polV[6]),
                                            vec(polV[1], polV[4], polV[7]),
                                            vec(polV[2], polV[5], polV[8]) );
                      
                      // Add data to tensor and charge maps
                      polTs[ namekey ] = polT;
                      chrgs[ namekey ] = chrg;

                } /* exit loop over atoms */
          } /* exit loop over residues */

          cout << "Molecule " << molname << ": "
               << " Found data for";
          if (anion)  { cout <<  " negative"; }
                        cout << ", neutral";
          if (cation) { cout << ", positive"; }
          cout << " state. " << endl;

    } /* exit loop over molecules */

    if ( polTs.size() != chrgs.size() ) {
         cout << "Inconsistent data set size "
              << "for charges and polarizabilities"
              << "." << endl;
         throw std::runtime_error("Check " + xmlfile + ".");
    }

    cout << "Found multipole data for "
         << polTs.size() << " atom types." << endl;

    // TODO Move out to atom class
    _polTs  = polTs;
    _polit  = _polTs.begin();
    _chrgs  = chrgs;
    _chrgit = _chrgs.begin();
}

/**
 * \brief Equip topology with atomistic charges
 * @param top
 * @return
 */
void EMultipole::EquipTop(QMTopology *top) {
    cout << "EMultipole::EquipTop *** RETURN ***" << endl;
    return;
    // NOTE Starting here, substitute <Molecule>
    // NOTE with QMCrgunit
    // NOTE top->CrgUnits() instead of top->Molecules()
    // NOTE Also iterate over atoms instead of beads
    vector < Bead * >::iterator atmit;

    for ( atmit = top->Beads().begin();
          atmit != top->Beads().end();
          atmit++ ) {

          string keyname =  (*atmit)->getName()+ "_" +
                            boost::lexical_cast< string, int > 
                                   ((*atmit)->getResnr()) + "_" +
                            (*atmit)->getMolecule()->getName();

          // UPDATE HERE
          (*atmit)->setQ( _chrgs[keyname][-1] ); //>>>
          (*atmit)->setQ( _chrgs[keyname][1] );  //>>>
          (*atmit)->setQ( _chrgs[keyname][0] );  //>>>
          //<<< (*atmit)->setQs( _chrgs[keyname] );
          //<<< (*atmit)->setPolT( _poliTs[keyname] );
      }
}

/**
 * \brief Compute electrostatic contribution to site energies
 * @param top
 * @return
 */
bool EMultipole::EvaluateFrame(QMTopology *top) {
    cout << "EMultipole::EvalFrame *** RETURN ***" << endl;
    return 0;
    // Site energies are attributes of charge units (= segments);
    // hence, first charge each segment appropriately,
    // then compute energy of that configuration
    // and compare to energy of neutral state.

    // Set electrostatic prefactor: All energies calculated via multipole
    // interaction tensors have to be multiplied by this factor => yields
    // energy in eV.
    //       m           10⁹   nm
    //    J --- e² = eV ----- ---- e²
    //       C²           e    C²
    double _r4PiEps0 = 1.602176487e-19 * 1.000e9 * 1/(4*M_PI*8.854187817e-12);

    // If polarizabilites are given in A°**3, this already includes 1/ 4PiEps0;
    // Hence multiply all energies calculated via multipole interaction tensors
    // by this factor => yields ind. dpl. moments in C*nm
    double _A3toNM3 = 1.000e-3;
    
    // Initialize vector to store induced dipole moments
    if ( !_muInd.size() ) {
         vector < Bead *>::iterator ait;
         for ( ait = top->Beads().begin();
               ait != top->Beads().end();
               ait++ ) {
               _muInd.push_back( vec(0,0,0) );
               _pField.push_back( vec(0,0,0) );
               _iField.push_back( vec(0,0,0) );
        }
    }

    // TODO Starting here, substitute <Molecule> with QMCrgunit;
    // TODO top->CrgUnits() instead of top->Molecules()
    // TODO Also change from bead to atom class.

    vector < Molecule * >::iterator chuit;
    for ( chuit = top->Molecules().begin();
          chuit != top->Molecules().end();
          chuit++ ) {

          // Zero out energy contributions for all states
          _pInter[-1] = _pInter[0] = _pInter[1] = 0.0;
          _iInter[-1] = _iInter[0] = _iInter[1] = 0.0;
          _pIntra[-1] = _pIntra[0] = _pIntra[1] = 0.0;
          _iIntra[-1] = _iIntra[0] = _iIntra[1] = 0.0;

          int state;
          // state -1 <=> anion
          // ...

          // TODD Add function here to determine which charge states
          // TODO to investigate for this molecule

          if (_anion) { 
              state = -1;
              ChargeMol(*chuit, state);
              Induce(top);
              CalcIntEnergy(top, state);
              Depolarize(top);
          }          
          if (_cation) {
              state = 1;
              ChargeMol(*chuit, state);
              Induce(top);
              CalcIntEnergy(top, state);
              Depolarize(top);
          }          
          if (true) {
              state = 0;
              ChargeMol(*chuit, state);
              Induce(top);
              CalcIntEnergy(top, state);
              Depolarize(top);
          }
          
          // Output results
          cout << "Segment " << (*chuit)->getName()
               << (*chuit)->getId() << endl;
          
          if (_cation) {
              double eNC = (_pInter[+1] + _iInter[+1])
                         - (_pInter[0] + _iInter[0]);
              
              cout << "E(\"Cation\") - E(\"Neutral\"): " << eNC << endl;
                   
              // TODO Call members of charge unit to set energies
          }
          
          if (_anion) {
              double eNA = (_pInter[-1] + _iInter[-1])
                         - (_pInter[0] + _iInter[0]); 
              
              cout << "E(\"Anion\")  - E(\"Neutral\"): " << eNA << endl;
                   
              // TODO Call members of charge unit to set energies
          }

    } /* exit loop over segments */

    return 1;
}


/**
 * \brief Charge atoms in molecule appropriately
 * @param top
 * @return
 */
void EMultipole::ChargeMol(Molecule *mol, int state) {

    // TODO Move this to Molecule.h such that
    // TODO mol->setChrgStat(state) suffices in ::EvaluateFrame

    for ( int i = 0; i < mol->BeadCount(); i++ ) {
        Bead *atm = mol->getBead(i);
        string key = atm->getName()+ "_" +                  //>>>
                     boost::lexical_cast< string, int >     //>>>
                            (atm->getResnr()) + "_" +       //>>>
                     mol->getName();                        //>>>
        mol->getBead(i)->setQ(_chrgs[key][state]);          //>>>
        //<<< atm->setQStat(state);

    }
}


/**
 * \brief Calculate induced dipole moments on the basis of Thole's dipole
 *        interaction model
 * @param top
 */
void EMultipole::Induce(QMTopology *top) {

    // Calculate field generated by charges, perm. multipoles

    vector < Bead *>::iterator ait;
    vector < Bead *>::iterator bit;
    for ( ait = top->Beads().begin();
          ait != top->Beads().end();
          ait++ ) {
          Bead *atm = *ait;
          
          for ( bit = ait+1;
                bit != top->Beads().end();
                bit++ ) {
                Bead *btm = *bit;

                if (btm->getMolecule() == atm->getMolecule()) { continue; }

                // Add field seen by a due to b (and vice versa) to field vector
                // Pay attention to directionality: R_ab != R_ba
                vec R = top->BCShortestConnection(
                        atm->getPos(), btm->getPos());

                _pField[atm->getId()] += - T1(R)*btm->getQ();
                _pField[btm->getId()] +=   T1(R)*atm->getQ();

                // TODO Apply 1-n scaling to intramolecular interactions.
          }
    }


  // Calculate induced moments in self-consistent manner

  bool converged = false;
  int  iter = 0;

  while (!converged && iter <= _maxIter) {

    iter++;

    // Reset field due to induced moments
    vector <vec>::iterator vit;
    for ( vit = _iField.begin();
          vit != _iField.end();
          vit++ ) {
          *vit = vec(0,0,0);
    }

    // Update induction field
    for ( ait = top->Beads().begin();
          ait != top->Beads().end();
          ait++ ) {

          Bead *atm = *ait;

          for ( bit = ait+1;
                bit != top->Beads().end();
                bit++ ) {

                Bead *btm = *bit;

                if (btm == atm) { continue; }

                int a = atm->getId();
                int b = btm->getId();

                vec R = top->BCShortestConnection(
                        atm->getPos(), btm->getPos());
                        // => R_b - R_a

                double u = abs(R) / pow(
                    _dpolT[atm->getId()].get(0,0) *
                    _dpolT[btm->getId()].get(0,0) , 1/6);

                _iField[a] = T2(R,u) * _muInd[b];
                _iField[b] = T2(R,u) * _muInd[a];

                bool intra = false;
                if (atm->getMolecule() == btm->getMolecule()) {intra = true;}
                // TODO Apply 1-n scaling to intramolecular interactions.

         } /* exit atom a loop */

    } /* exit atom b loop */

    // Update induced moments using induction field
    converged = true;

    for ( ait = top->Beads().begin();
          ait != top->Beads().end();
          ait++ ) {

          Bead *atm = *ait;
          int a = atm->getId();

          vec muIndNew = _dpolT[a] * (_pField[a] + _iField[a]);

          // Test for convergence
          double deltaMu = abs(muIndNew - _muInd[a])/abs(_muInd[a]);
          if ( deltaMu > _epsTol ) { converged = false; }

          _muInd[a] = muIndNew;
    }

  } /* exit converge loop */

  if (iter > _maxIter ) {
      cout << "Warning: Induced dipole moments may not be converged." << endl;
  }

}


/**
 * \brief Zero out induced dipoles moment
 * @param top
 */
void EMultipole::Depolarize(QMTopology *top) {

    vector < Bead *>::iterator ait;
    for ( ait = top->Beads().begin();
          ait != top->Beads().end();
          ait++ ) {
          _muInd.push_back( vec(0,0,0) );
    }
}


/**
 * \brief Calculate energy for specific charge state of system
 * @param top
 * @param state
 */
void EMultipole::CalcIntEnergy(QMTopology *top, int &state) {

    // pick bead from molecule, pick bead from topology
    // check: intra- or inter- ?
    // if intra: apply 1-n scaling
    vector < Bead *>::iterator ait;
    for ( ait = top->Beads().begin();
          ait != top->Beads().end();
          ait++ ) {

          Bead *atm = *ait;

          vector < Bead *>::iterator bit;
          for ( bit = ait+1;
                bit != top->Beads().end();
                bit++ ) {

                Bead *btm = *bit;

                // In principal, intramolecular interactions are already
                // included in QM energies. Hence, focus on intermolecular
                // energy and keep track of intramolecular contributions just
                // to document. Accordingly:

                // "intra-action" or interaction ?
                bool intra = false;
                if ( atm->getMolecule() == btm->getMolecule() ) {intra = true;}
                
                // It matters whether R = R_ba or = R_ab;
                // to work with equations further below, R = R_b - R_a.
                vec R = top->BCShortestConnection(atm->getPos(), btm->getPos());

                double u = abs(R) / pow(
                    _dpolT[atm->getId()].get(0,0) *
                    _dpolT[btm->getId()].get(0,0) , 1/6);

                double ePerm = btm->getQ() * this->T0(R) * atm->getQ();

                // Note the factor 1/2, which accounts for the work necessary
                // to induce multipoles.
                double eInd = 0.5*(
                    _muInd[btm->getId()] *   this->T1(R,u) * atm->getQ() +
                    _muInd[atm->getId()] * - this->T1(R,u) * btm->getQ() );

                // TODO Apply 1-n scaling to intramolecular interactions.
                // TODO This scaling, however, is in this context only
                // TODO important for the induction process in ::Induce.

                if (intra) {
                    _pIntra[state] += ePerm;
                    _iIntra[state] += eInd;
                }
                else {
                    _pInter[state] += ePerm;
                    _iInter[state] += eInd;
                }
          } /* exit loop over atom b */
    } /* exit loop over atom a */
}

/**
  * \brief Print info to screen
  */

void EMultipole::PrintInfo(const string &key) {

  if (key == "Thole") {
	cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++ " << endl;
    cout << endl;
    cout << "Thole Model Parameters ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ " << endl;
    cout << "Use cutoff?        " << _useCutoff << endl;
    if (_useCutoff) {
        cout << "Cutoff [nm]        " << _cutoff << endl;
    }
    cout << "Use exp. damping?  " << _useExp << endl;
    if (_useExp) {
        cout << "Exp. damp. factor  " << _aDamp << endl;
    }
    cout << "Use 1-n scaling?   " << _useScaling << endl;
    if (_useScaling) {
        cout << "1-n scaling        ";
        for (int i = 2; i < 5; i++) {
            cout << "1-" << i << ": " << _scale1[i-2] << " | "; }
        cout << endl;
    }
    cout << "Convergence Loop Parameters ~~~~~~~~~~~~~~~~~~~~~~~~~~ " << endl;
    cout << "Omega SOR          " << _omegSOR << endl;
    cout << "Tolerance          " << _epsTol << endl;
    cout << "Max. iteration     " << _maxIter << endl << endl;
    
  }

  else if (key == "MPoles") {
    cout << "Charges ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ " << endl;
    for ( _chrgit = _chrgs.begin(); _chrgit != _chrgs.end(); _chrgit++ ) {
          cout << _chrgit->first << " ";
          map <int, double> tmp = _chrgit->second;
          cout << " -1:" << tmp[-1];
          cout <<  " | 0:" << tmp[0];
          cout << " | +1:" << tmp[1] << endl;
    }
    cout << "Polarizability Tensors ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ " << endl;
    for ( _polit = _polTs.begin(); _polit != _polTs.end(); _polit ++ ) {
          cout << _polit->first << endl;
          cout << _polit->second;
    }
    cout << endl;
  }
}

}} /* exit namespaces votca, ctp */

