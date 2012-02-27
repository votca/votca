/* 
 * File:   EMultipole.cc
 * Author: poelking
 * 
 * Created on December 27, 2011, 4:10 PM
 */

#include "emultipole.h"
#include <votca/tools/globals.h>

namespace votca { namespace ctp {

/**
 * \brief Read in parameters for Thole model and convergence loop
 * @param top
 * @param opt
 */
void EMultipole::Initialize(Topology *top, Property *opt) {

    cout << endl <<  "... ... Parametrizing Thole model";

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

    // Restate input data to screen
    if (true) { // OVERRIDE verbose
        this->PrintInfo("Thole"); this->PrintInfo("MPoles");
    }
}

/**
 * \brief Read in multipole data for atom types from XML file.
 * @param xmlfile
 */
void EMultipole::GetMPoles(string &xmlfile) {

     cout << endl << "... ... Import multipoles from " << xmlfile;

     // map to store polarizability tensors
     //     string -> atom identifier
     //     matrix -> polarizability tensor
     map < string, matrix >           ptensors;
     map < string, matrix >::iterator ptit;

     // map to store charges [ anion, neutral, cation ]
     //     string -> atom identifier
     //     int    -> charge state (-1, 0, 1)
     //     double -> charge
     map < string, map<int, double> >           chrgs;
     map < string, map<int, double> >::iterator chrgit;

     // map to keep track of which segments allow for which charge states
     //    string -> segment identifier
     //    int    -> charge state
     //    bool   -> has charge state?
     map < string, map<int, bool> > hasChrg;


    /* ---- MULTIPOLES.XML Structure ----
     *  <segment>
     *      <name></name>
     * 
     *      <fragment>
     *          <name></name>
     *   
     *          <atom>
     *              <name></name>
     *              <chrgA></chrgA>
     *              <chrgN></chrgN>
     *              <chrgC></chrgC>
     *              <ptensor></ptensor>
     *          </atom>  
     *     
     *      </fragment>
     *  </segment>
     */

    // Load options from specified xml file
    Property opt;
    load_property_from_xml(opt, xmlfile.c_str());

    // Loop through molecules
    list< Property * > segs = opt.Select("segments.segment");
    list< Property * >::iterator sit;
    for ( sit = segs.begin();
          sit != segs.end();
          sit++ ) {

          string segname = (*sit)->get("name").as< string >();

          bool anion = true;
          bool cation = true;
          bool neutral = true;

          // Loop through residues
          list< Property * > frags = (*sit)->Select("fragment");
          list< Property * >::iterator fit;
          for ( fit = frags.begin();
                fit != frags.end();
                fit++ ) {

                string fragname = (*fit)->get("name").as< string >();

                // Loop through atoms
                list< Property * > atms = (*fit)->Select("atom");
                list< Property * >::iterator ait;
                for ( ait = atms.begin();
                      ait != atms.end();
                      ait++ ) {

                      string atmname = (*ait)->get("name").as< string >();

                      // Putting atmname first should speed up look-up process
                      // later on in ::EquipTop, ::...
                      string namekey = atmname+"_"+fragname+"_"+segname;

                      // Check whether data already stored. If so, continue.
                      ptit = ptensors.find(namekey);
                      chrgit = chrgs.find(namekey);
                      if ( ptit != ptensors.end() &&
                           chrgit != chrgs.end() ) {
                          cout << "... ... WARNING: multiple charge def.s for "
                                  "(ATOM_FRAG_SEG) " << namekey << endl;
                          continue;
                      }

                      // Retrieve charges (anion, neutral, cation)                      
                      map< int, double > chrg;
                      double chrgA;
                      double chrgN;
                      double chrgC;
                      
                      try { chrgA = (*ait)->get("chrgA").as< double >();}
                      catch ( std::runtime_error ) { anion = false; }

                      try { chrgN = (*ait)->get("chrgN").as< double >(); }
                      catch ( std::runtime_error ) { neutral = false; }
                      chrgN = (*ait)->get("chrgN").as< double >();

                      try { chrgC = (*ait)->get("chrgC").as< double >();}
                      catch ( std::runtime_error ) { cation = false; }

                      if (anion) { chrg[-1] = chrgA; }
                      if (neutral) { chrg[0] = chrgN; }
                      if (cation) { chrg[1] = chrgC; }

                      // Retrieve dipole polarizability tensor
                      vector<double> polV =
                         (*ait)->get("ptensor").as< vector<double> >();

                      if ( polV.size() != 9 ) {
                         cout << "... ... Invalid input for polarizability of "
                              << "atom " << atmname
                              << "in fragment " << fragname
                              << "in segment " << segname
                              << "." << endl;
                         cout << "... ... Format is: xx yx zx xy yy zy xz yz zz"
                              << "." << endl;
                         throw std::runtime_error("Error in xml file");
                      }

                      matrix polT = matrix( vec(polV[0], polV[3], polV[6]),
                                            vec(polV[1], polV[4], polV[7]),
                                            vec(polV[2], polV[5], polV[8]) );
                      
                      // Add data to tensor and charge maps
                      ptensors[ namekey ] = polT;
                      chrgs[ namekey ] = chrg;

                } /* exit loop over atoms */
          } /* exit loop over residues */

          // Log charge states processed for this segment
          hasChrg[segname][-1] = anion;
          hasChrg[segname][0]  = neutral;
          hasChrg[segname][+1] = cation;

          cout << endl << "... ... Segment " << segname << ": "
               << "Found data for ";
          if (anion)  { cout << "negative "; }
          if (neutral){ cout << "neutral "; }
          if (cation) { cout << "positive "; }
          cout << "state. ";

    } /* exit loop over molecules */

    if ( ptensors.size() != chrgs.size() ) {
         cout << "Inconsistent data set size for charges and polarizabilities."
              << endl;
         throw std::runtime_error("Error in xml file.");
    }

    cout << endl << "... ... Found multipole data for "
         << ptensors.size() << " atom types.";

    // TODO Move out to atom class
    _ptensors  = ptensors;
    _ptit  = _ptensors.begin();
    _chrgs  = chrgs;
    _chrgit = _chrgs.begin();
    _hasChrg = hasChrg;
}

/**
 * \brief Equip topology with atomistic charges
 * @param top
 * @return
 */
void EMultipole::EquipTop(Topology *top) {
    cout << endl << "... ... Equip topology.";

    map< string, bool > skip;

    // Explore charge states for each segment
    vector < Segment * > ::iterator segit;
    for (segit = top->Segments().begin();
            segit < top->Segments().end();
            segit++) {
        try {
            bool anion = _hasChrg.at((*segit)->getName()).at(-1);
            bool neutral = _hasChrg.at((*segit)->getName()).at(0);
            bool cation = _hasChrg.at((*segit)->getName()).at(+1);
            (*segit)->AddChrgState(-1, anion);
            (*segit)->AddChrgState(0, neutral);
            (*segit)->AddChrgState(+1, cation);
            /*
            cout << endl << "... ... " << (*segit)->getId() << " "
                    << (*segit)->getName() << " " << anion << " "
                    << neutral << " " << cation << ".";
            */
        }
        catch (out_of_range) {
            if ( ! skip.count((*segit)->getName())) {
                 cout << endl << "... ... WARNING No charge input for segment "
                     << (*segit)->getName() << ". Skipping all ... ";
                 skip[(*segit)->getName()] = true;
            }
            (*segit)->AddChrgState(-1, false);
            (*segit)->AddChrgState(0, false);
            (*segit)->AddChrgState(+1, false); 
        }
    }

    // Equip atoms with charges
    for (segit = top->Segments().begin();
            segit < top->Segments().end();
            segit++) {

        bool anion = (*segit)->hasChrgState(-1);
        bool neutral = (*segit)->hasChrgState(0);
        bool cation = (*segit)->hasChrgState(+1);

        vector < Atom * >::iterator atmit;
        for ( atmit = (*segit)->Atoms().begin();
              atmit != (*segit)->Atoms().end();
              atmit++ ) {             

              if ( anion || neutral || cation) {
                   string keyname =  (*atmit)->getName()
                                   + "_" + (*atmit)->getFragment()->getName()
                                   + "_" + (*atmit)->getSegment()->getName();
                   (*atmit)->setQ(_chrgs.at(keyname));
                   (*atmit)->setPTensor(_ptensors.at(keyname));
              }
              else {
                  map< int, double > zero;
                  zero[-1] = zero[0] = zero[1] = 0;
                  matrix zeroes;
                  zeroes.ZeroMatrix();
                  (*atmit)->setQ(zero);
                  (*atmit)->setPTensor(zeroes);
              }

              /*
              cout << endl << "... ... " << (*atmit)->getId() << " "
                      << (*atmit)->getQ(-1) << " "
                      << (*atmit)->getQ(0) << " "
                      << (*atmit)->getQ(1);
              cout << endl << (*atmit)->getPTensor().getCol(0);
              cout << endl << (*atmit)->getPTensor().getCol(1);
              cout << endl << (*atmit)->getPTensor().getCol(2);
              */
        }
    }
}

/**
 * \brief Compute electrostatic contribution to site energies
 * @param top
 * @return
 */
bool EMultipole::EvaluateFrame(Topology *top) {
    cout << endl << "... ... EvaluateFrame ";

    // Equip topology with charges; have to do this here, because
    // on the '::Initialize' stage the topology was not yet populated
    this->EquipTop(top);

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
    // Hence multiply all moments calculated via multipole interaction tensors
    // by this factor => yields ind. dpl. moments in C*nm
    double _A3toNM3 = 1.000e-3;
    
    // Initialize vector to store induced dipole moments,
    // primary field (from permanent multipoles) and secondary field
    if ( !_muInd.size() ) {
         vector < Atom *>::iterator ait;
         for ( ait = top->Atoms().begin();
               ait != top->Atoms().end();
               ait++ ) {
               _muInd.push_back( vec(0,0,0) );
               _pField.push_back( vec(0,0,0) );
               _iField.push_back( vec(0,0,0) );
               (*ait)->chrg(0);
        }
    }

    // TODO Starting here, substitute <Molecule> with QMCrgunit;
    // TODO top->CrgUnits() instead of top->Molecules()
    // TODO Also change from Atom to atom class.

    // TODO Parallelise starting here

    vector < Segment * >::iterator sit;
    for ( sit = top->Segments().begin();
          sit != top->Segments().end();
          sit++ ) {

          Segment *seg = *sit;

          cout << endl << "... ... Segment " << seg->getName() << seg->getId();

          // Zero out energy contributions for all states
          _pInter[-1] = _pInter[0] = _pInter[1] = 0.0;
          _iInter[-1] = _iInter[0] = _iInter[1] = 0.0;
          _pIntra[-1] = _pIntra[0] = _pIntra[1] = 0.0;
          _iIntra[-1] = _iIntra[0] = _iIntra[1] = 0.0;

          int state;
          // state -1 <=> anion
          // ...

          if (seg->hasChrgState(-1)) {
              state = -1;
              seg->chrg(state);
              Induce(top);
              CalcIntEnergy(top, state);
              Depolarize(top);
          }          
          if (seg->hasChrgState(1)) {
              state = 1;
              seg->chrg(state);
              Induce(top);
              CalcIntEnergy(top, state);
              Depolarize(top);
          }          
          if (seg->hasChrgState(0)) {
              state = 0;
              seg->chrg(state);
              Induce(top);
              CalcIntEnergy(top, state);
              Depolarize(top);
          }
          
          // Output results

          if (seg->hasChrgState(+1)) {
              double eNC = (_pInter[+1] + _iInter[+1])
                         - (_pInter[0] + _iInter[0]);
              
              cout << endl 
                   << "... ... ... E(\"Cation\") - E(\"Neutral\"): " << eNC;
                   
              // TODO Call members of charge unit to set energies
          }
          
          if (seg->hasChrgState(-1)) {
              double eNA = (_pInter[-1] + _iInter[-1])
                         - (_pInter[0] + _iInter[0]); 

              cout << endl
                   << "... ... ... E(\"Anion\")  - E(\"Neutral\"): " << eNA;
                   
              // TODO Call members of charge unit to set energies
          }

          cout << endl
               << "... ... ... e_P_Inter " << _pInter[-1] << " | " << _pInter[0] << " | " << _pInter[1];
          cout << endl
               << "... ... ... e_I_Inter " << _iInter[-1] << " | " << _iInter[0] << " | " << _iInter[1];
          cout << endl
               << "... ... ... e_P_Intra " << _pIntra[-1] << " | " << _pIntra[0] << " | " << _pIntra[1];
          cout << endl
               << "... ... ... e_I_Intra " << _iIntra[-1] << " | " << _iIntra[0] << " | " << _iIntra[1];



          return 1; // OVERRIDE

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

    for ( int i = 0; i < mol->Atoms().size(); i++ ) {
        Atom *atm = mol->getAtom(i);
        atm->chrg(state);

    }
}


/**
 * \brief Calculate induced dipole moments on the basis of Thole's dipole
 *        interaction model
 * @param top
 */
void EMultipole::Induce(Topology *top) {

    cout << endl
         << "... ... ... Induce *** RETURN ***";
    return;
    
    // Calculate field generated by charges, perm. multipoles

    vector < Atom *>::iterator ait;
    vector < Atom *>::iterator bit;
    for ( ait = top->Atoms().begin();
          ait != top->Atoms().end();
          ait++ ) {
          Atom *atm = *ait;
          
          for ( bit = ait+1;
                bit != top->Atoms().end();
                bit++ ) {
                Atom *btm = *bit;

                if (btm->getSegment() == atm->getSegment()) { continue; }

                // Add field seen by a due to b (and vice versa) to field vector
                // Pay attention to directionality: R_ab != R_ba
                vec R = top->PbShortestConnect(
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
    for ( ait = top->Atoms().begin();
          ait != top->Atoms().end();
          ait++ ) {

          Atom *atm = *ait;

          for ( bit = ait+1;
                bit != top->Atoms().end();
                bit++ ) {

                Atom *btm = *bit;

                if (btm == atm) { continue; }

                int a = atm->getId();
                int b = btm->getId();

                vec R = top->PbShortestConnect(
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

    for ( ait = top->Atoms().begin();
          ait != top->Atoms().end();
          ait++ ) {

          Atom *atm = *ait;
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
void EMultipole::Depolarize(Topology *top) {
    cout << endl
         << "... ... ... Depolarize";

    vector < vec >::iterator vit;
    for ( vit = _muInd.begin();
          vit < _muInd.end();
          vit++ ) {
          (*vit) = vec(0,0,0);
    }
}


/**
 * \brief Calculate energy for specific charge state of system
 * @param top
 * @param state
 */
void EMultipole::CalcIntEnergy(Topology *top, int &state) {

    cout << endl
         << "... ... ... Interaction energy ";


    // Progress bar
    int full = top->Atoms().size();
    int cent = int(full / 100);
    cout << endl;

    // pick Atom from molecule, pick Atom from topology
    // check: intra- or inter- ?
    // if intra: apply 1-n scaling
    vector < Atom *>::iterator ait;
    for ( ait = top->Atoms().begin();
          ait != top->Atoms().end();
          ait++ ) {

          Atom *atm = *ait;

          if (atm->getId() % cent == 0) {
              this->PrintProgBar( int(atm->getId() / cent) ) ;
          }

          vector < Atom *>::iterator bit;
          for ( bit = ait+1;
                bit != top->Atoms().end();
                bit++ ) {


                Atom *btm = *bit;

               //cout << endl
                    //<< atm->getId() << " | " << btm->getId() << ": ";
               
                // In principal, intramolecular interactions are already
                // included in QM energies. Hence, focus on intermolecular
                // energy and keep track of intramolecular contributions just
                // to document. Accordingly:

                // It matters whether R = R_ba or = R_ab;
                // to work with equations further below, R = R_b - R_a.
                vec R = top->PbShortestConnect(atm->getPos(), btm->getPos());
                double absR = abs(R);
                if (_useCutoff && absR > _cutoff) { continue; }

                // "intra-action" or interaction ?
                bool intra = false;
                if ( atm->getSegment() == btm->getSegment() ) { intra = true; }
                
                //cout << "Intra " << intra << "; ";

                double u = absR / pow(
                    atm->getPTensor().get(0,0) *
                    btm->getPTensor().get(0,0) , 1/6);
                
                double ePerm = btm->getQ() * this->T0(R) * atm->getQ();

                //cout << "ePerm " << ePerm << "; ";
                
                // Note the factor 1/2, which accounts for the work necessary
                // to induce multipoles.
                double eInd = 0.5*(
                    _muInd[btm->getId()-1] *   this->T1(R,u) * atm->getQ() +
                    _muInd[atm->getId()-1] * - this->T1(R,u) * btm->getQ() );

                //cout << "eInd " << eInd << "; ";

                // TODO Apply 1-n scaling to intramolecular interactions.
                //      This scaling, however, is in this context only
                //      important for the induction process in ::Induce.

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

    cout << endl
	 << "++++++++++++++++++++++++++++++++++++++++++++++++++++++ " << endl;
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
    for ( _ptit = _ptensors.begin(); _ptit != _ptensors.end(); _ptit ++ ) {
          cout << _ptit->first << endl;
          cout << _ptit->second;
    }
    cout << endl;
  }
}

/**
 * \brief Progress bar with percent indicator
 * @param percent
 */
void EMultipole::PrintProgBar(int percent) {
      string bar;

      for (int i = 0; i < 25; i++) {
           if( i < (percent/4)) {
               bar.replace(i,1,"=");
           }
           else if ( i == (percent/4)) {
               bar.replace(i,1,">");
           }
           else {
               bar.replace(i,1," ");
           }
      }

      cout<< "\r" "... ... ... |" << bar << "| ";
      cout.width(3);
      cout<< percent << "%     " << flush;
}




}} /* exit namespaces votca, ctp */









