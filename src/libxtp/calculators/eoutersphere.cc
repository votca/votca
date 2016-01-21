/*
 *            Copyright 2009-2016 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */


#include "eoutersphere.h"



namespace votca { namespace xtp {
    
// +++++++++++++++++++++++++++++ //
// EOUTERSPHERE MEMBER FUNCTIONS //
// +++++++++++++++++++++++++++++ //

void EOutersphere::Initialize(Property *opt) {

    // update options with the VOTCASHARE defaults   
    UpdateWithDefaults( opt );

    string key = "options." + Identify();

    /* ---- OPTIONS.XML Structure -----
     *
     * <eoutersphere>
     *
     *      <multipoles>MULTIPOLES.XML</multipoles>
     *
     *      <method></method>
     *      <pekar></pekar>
     *      <cutoff></cutoff>
     *
     * </eoutersphere>
     *
     */

    _method = opt->get(key+".method").as<string> ();

    if (_method == "constant") {
        _lambdaConstant = opt->get(key+".lambdaconst").as<double> ();
    }
    else if (_method == "spheres") {
        _pekarFactor = opt->get(key+".pekar").as<double> ();
        list< Property* > typeinfo = opt->Select(key+".segment");
        list< Property* > ::iterator rit;
        for (rit = typeinfo.begin(); rit != typeinfo.end(); ++rit) {
            string type = (*rit)->get("type").as<string>();
            double radius = (*rit)->get("radius").as<double>();
            _radius[type] = radius;
        }
    }
    else if (_method == "dielectric") {
        _pekarFactor = opt->get(key+".pekar").as<double> ();
        _lambdaCutoff = opt->get(key+".cutoff").as<double> ();
    }
    else {
        cout << "ERROR: Typo? Method for reorg. energy calc.: " << _method;
        throw std::runtime_error("Unrecognized <method> in options file.");
    }
    

    this->EStatify(NULL, opt);
}


void EOutersphere::EStatify(Topology *top, Property *options) {

    cout << endl << "... ... Estatify system: ";

    string key = "options.eoutersphere";
    string allocFile = options->get(key+".multipoles").as<string> ();

    // ++++++++++++++++++++++++++++++++ //
    // Load polar-site indices from XML //
    // ++++++++++++++++++++++++++++++++ //

    // => Output to maps:
    map<string, vector<int> > alloc_frag_mpoleIdx;
    map<string, vector<string> > alloc_frag_mpoleName;
    map<string, string > alloc_seg_mpoleFiles_n;
    map<string, string > alloc_seg_mpoleFiles_e;
    map<string, string > alloc_seg_mpoleFiles_h;

    Property allocation; // <- Which polar sites are part of which fragment?
    load_property_from_xml(allocation, allocFile.c_str());


    /* --- MULTIPOLES.XML Structure ---
     *
     * <topology>
     *
     *     <molecules>
     *          <molecule>
     *          <name></name>
     *
     *          <segments>
     *
     *              <segment>
     *              <name>DCV</name>
     *
     *              <multipoles_n></multipoles_n>
     *              <multipoles_e></multipoles_e>
     *              <multipoles_h></multipoles_h>
     *
     *              <map2md></map2md>
     *
     *              <fragments>
     *                  <fragment>
     *                  <name></name>
     *                  <mpoles></mpoles>
     *                  </fragment>
     *              </fragments>
     *              ...
     *              ...
     */


    key = "topology.molecules.molecule";
    list<Property *> mols = allocation.Select(key);
    list<Property *>::iterator molit;
    for (molit = mols.begin(); molit != mols.end(); molit++) {

        string molName = (*molit)->get("name").as<string> ();

        key = "segments.segment";
        list<Property *> segs = (*molit)->Select(key);
        list<Property *>::iterator segit;

        for (segit = segs.begin(); segit != segs.end(); segit++) {

            string segName = (*segit)->get("name").as<string> ();

            // GDMA filenames for cation (h), neutral (n), anion (e) state
            string mpoleFile_n = (*segit)->get("multipoles_n").as<string> ();
            alloc_seg_mpoleFiles_n[segName] = mpoleFile_n;
            if ( (*segit)->exists("multipoles_e")) {
                string mpoleFile_e = (*segit)->get("multipoles_e").as<string>();
                alloc_seg_mpoleFiles_e[segName] = mpoleFile_e;
            }
            if ( (*segit)->exists("multipoles_h")) {
                string mpoleFile_h = (*segit)->get("multipoles_h").as<string>();
                alloc_seg_mpoleFiles_h[segName] = mpoleFile_h;
            }

            // Default: Project multipoles onto rigidified coordinates
            if ( (*segit)->exists("map2md")) {
                int map2md = (*segit)->get("map2md").as<int>();
                _map2md[segName] = (map2md == 0) ? false : true;
            }
            else {
                _map2md[segName] = false; // i.e. map to rigidified coordinates
            }

            key = "fragments.fragment";
            list<Property *> frags = (*segit)->Select(key);
            list<Property *>::iterator fragit;

            for (fragit = frags.begin(); fragit != frags.end(); fragit++) {

                string fragName = (*fragit)->get("name").as<string> ();
                string mapKeyName = fragName + segName + molName;

                string mpoles = (*fragit)->get("mpoles").as<string> ();

                Tokenizer tokPoles(mpoles, " \t\n");
                vector<string> mpoleInfo;
                tokPoles.ToVector(mpoleInfo);

                vector<int> mpoleIdcs;
                vector<string> mpoleNames;

                vector<string> ::iterator strit;
                for (strit=mpoleInfo.begin(); strit<mpoleInfo.end(); strit++) {

                    Tokenizer tokPoleInfo( (*strit), " :");
                    vector<string> poleInfo;
                    tokPoleInfo.ToVector(poleInfo);

                    int mpoleIdx = boost::lexical_cast<int>(poleInfo[0]);
                    string mpoleName = poleInfo[1];

                    mpoleIdcs.push_back(mpoleIdx);
                    mpoleNames.push_back(mpoleName);

                }

                alloc_frag_mpoleIdx[mapKeyName] = mpoleIdcs;
                alloc_frag_mpoleName[mapKeyName] = mpoleNames;
            }
        }
    }

    // ++++++++++++++++++++++ //
    // Parse GDMA punch files //
    // ++++++++++++++++++++++ //

    // => Output to PolarSite Container (template container)
    map<string, vector<PolarSite*> > map_seg_polarSites;
    map<string, vector<bool> >       map_seg_chrgStates;


    // Multipoles for neutral state
    int state = 0;
    map<string, string > ::iterator strit;
    for (strit = alloc_seg_mpoleFiles_n.begin();
         strit != alloc_seg_mpoleFiles_n.end();
         strit++) {

        string segName = strit->first;
        string filename = strit->second;

        vector< PolarSite* > poles = ParseGdmaFile(filename, state);

        map_seg_polarSites[segName] = poles;
        map_seg_chrgStates[segName].resize(3);
        map_seg_chrgStates[segName][0] = false; // <- negative
        map_seg_chrgStates[segName][1] = true;  // <- neutral
        map_seg_chrgStates[segName][2] = false; // <- positive
        poles.clear();

    }

    // Multipoles for negative state
    state = -1;
    if ( alloc_seg_mpoleFiles_e.size() ) {
        for (strit = alloc_seg_mpoleFiles_e.begin();
             strit != alloc_seg_mpoleFiles_e.end();
             strit++) {
            string segName = strit->first;
            string filename = strit->second;

            vector< PolarSite* > polesAnion = ParseGdmaFile(filename, state);
            map_seg_chrgStates[segName][state+1] = true;

            // Merge with polar sites for neutral state
            vector< PolarSite* > polesNeutral = map_seg_polarSites[segName];

            assert(polesAnion.size() == polesNeutral.size());
            for (unsigned int i = 0; i < polesNeutral.size(); i++) {

                polesNeutral[i]->setQs( polesAnion[i]->getQs(state), state );
                polesNeutral[i]->setPs( polesAnion[i]->getPs(state), state );
                delete polesAnion[i];
            }
            polesAnion.clear();
        }
    }

    // Multipoles for positive state
    state = +1;
    if ( alloc_seg_mpoleFiles_h.size() ) {
        for (strit = alloc_seg_mpoleFiles_h.begin();
             strit != alloc_seg_mpoleFiles_h.end();
             strit++) {
            string segName = strit->first;
            string filename = strit->second;

            vector< PolarSite* > polesCation = ParseGdmaFile(filename, state);
            map_seg_chrgStates[segName][state+1] = true;

            // Merge with polar sites for neutral state
            vector< PolarSite* > polesNeutral = map_seg_polarSites[segName];

            assert(polesCation.size() == polesNeutral.size());
            for (unsigned int i = 0; i < polesNeutral.size(); i++) {

                polesNeutral[i]->setQs( polesCation[i]->getQs(state), state );
                polesNeutral[i]->setPs( polesCation[i]->getPs(state), state );
                delete polesCation[i];
            }
            polesCation.clear();
        }
    }

    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
    // Forward information on polar-site templates to topology //
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++ //

    // Containers to use for mapping
    // map<string, vector<PolarSite*> >     map_seg_polarSites;
    // map<string, vector<int> >            alloc_frag_mpoleIdx;
    // map<string, vector<string> >         alloc_frag_mpoleName;

    // TO CONTINUE FROM HERE ...
    // Loop over segments, loop over fragments in segments.
    // => For each segment, look up associated polar sites using map 1;
    //      => For each fragment, look up associated
    //         polar sites idcs. using map 2;
    //          => For each index in retrieved container, create new polar site,
    //             importing information from polar site container of segment
    //          => Translate positions to MD frame, copy local frame from
    //             fragment

    _map_seg_polarSites = map_seg_polarSites;
    _map_seg_chrgStates = map_seg_chrgStates;
    _alloc_frag_mpoleIdx =  alloc_frag_mpoleIdx;
    _alloc_frag_mpoleName = alloc_frag_mpoleName;
}


vector<PolarSite*> EOutersphere::ParseGdmaFile(string filename, int state) {

    int poleCount = 1;
    double Q0_total = 0.0;
    string units = "";
    bool useDefaultPs = true;
    bool warn_anisotropy = false;

    vector<PolarSite*> poles;
    PolarSite *thisPole = NULL;

    vector<double> Qs; // <- multipole moments
    double         P1; // <- dipole polarizability

    std::string line;
    std::ifstream intt;
    intt.open(filename.c_str());

    if (intt.is_open() ) {
        while ( intt.good() ) {

            std::getline(intt, line);
            vector<string> split;
            Tokenizer toker(line, " \t");
            toker.ToVector(split);

            if ( !split.size()      ||
                  split[0] == "!"   ||
                  split[0].substr(0,1) == "!" ) { continue; }

    // ! Interesting information here, e.g.
    // ! DCV2T opt
    // ! SP        RB3LYP          6-311+G(d,p)
    // Units bohr
    //
    // C          -4.2414603400   -3.8124751600    0.0017575736    Rank  2
    //  -0.3853409355
    //  -0.0002321905   0.2401559510   0.6602334308
    //  -0.7220625314   0.0004894995  -0.0003833545   0.4526409813  -0.50937399
    //  P 1.75


            // Units used
            if ( split[0] == "Units") {
                units = split[1];
                if (units != "bohr" && units != "angstrom") {
                    throw std::runtime_error( "Unit " + units + " in file "
                                            + filename + " not supported.");
                }
            }

            // element,  position,  rank limit
            else if ( split.size() == 6 && split[4] == "Rank" ) {

                Qs.clear();
                P1 = -1.;

                int id = poleCount++;  // <- starts from 1
                string name = split[0];

                double BOHR2NM = 0.0529189379;
                double ANGSTROM2NM = 0.1;
                double x, y, z;

                if (units == "bohr") {
                    x = BOHR2NM * boost::lexical_cast<double>(split[1]);
                    y = BOHR2NM * boost::lexical_cast<double>(split[2]);
                    z = BOHR2NM * boost::lexical_cast<double>(split[3]);
                }
                else if (units == "angstrom") {
                    x = ANGSTROM2NM * boost::lexical_cast<double>(split[1]);
                    y = ANGSTROM2NM * boost::lexical_cast<double>(split[2]);
                    z = ANGSTROM2NM * boost::lexical_cast<double>(split[3]);
                }
                else {
                    throw std::runtime_error( "Unit " + units + " in file "
                                            + filename + " not supported.");
                }
                
                vec pos = vec(x,y,z);

                int rank = boost::lexical_cast<int>(split[5]);

                PolarSite *newPole = new PolarSite(id, name);
                newPole->setRank(rank);
                newPole->setPos(pos);
                poles.push_back(newPole);
                thisPole = newPole;

            }

            // 'P', dipole polarizability
            else if ( split[0] == "P" && (split.size() == 2 || split.size() == 7) ) {
                if (split.size() == 7) {
                    warn_anisotropy = true;
                }
                P1 = 1e-3 * boost::lexical_cast<double>(split[1]);
                thisPole->setPs(P1, state);
                useDefaultPs = false;
            }

            // multipole line
            else {

                int lineRank = int( sqrt(thisPole->getQs(state).size()) + 0.5 );

                if (lineRank == 0) {
                    Q0_total += boost::lexical_cast<double>(split[0]);
                }

                for (unsigned int i = 0; i < split.size(); i++) {

                    double qXYZ = boost::lexical_cast<double>(split[i]);

                    // Convert e*(a_0)^k to e*(nm)^k where k = rank
                    double BOHR2NM = 0.0529189379;
                    qXYZ *= pow(BOHR2NM, lineRank); // OVERRIDE

                    Qs.push_back(qXYZ);

                }
                thisPole->setQs(Qs, state);
            }

        } /* Exit loop over lines */
    }
    else { cout << endl << "ERROR: No such file " << filename << endl; }

    cout << endl << "... ... ... Reading " << filename <<
                    ": Q0(Total) = " << Q0_total << flush;

    if (useDefaultPs) {

        cout << endl << "... ... ... NOTE Using default Thole polarizabilities "
             << "for charge state " << state << ". ";

        vector< PolarSite* > ::iterator pol;
        for (pol = poles.begin(); pol < poles.end(); ++pol) {
            string elem = (*pol)->getName();
            double alpha = 0.0;
            // Original set of Thole polarizabilites
            if      (elem == "C") { alpha = 1.75e-3;  } // <- conversion from
            else if (elem == "H") { alpha = 0.696e-3; } //    A³ to nm³ = 10⁻³
            else if (elem == "N") { alpha = 1.073e-3; }
            else if (elem == "O") { alpha = 0.837e-3; }
            else if (elem == "S") { alpha = 2.926e-3; }
            // Different set of Thole polarizabilities
            //if      (elem == "C") { alpha = 1.334e-3; } // <- conversion from
            //else if (elem == "H") { alpha = 0.496e-3; } //    A³ to nm³ = 10⁻³
            //else if (elem == "N") { alpha = 1.073e-3; }
            //else if (elem == "O") { alpha = 0.837e-3; }
            //else if (elem == "S") { alpha = 3.300e-3; }
            else { throw runtime_error("No polarizability given "
                                       "for polar site type " + elem + ". "); }
            (*pol)->setPs(alpha, state);
        }
    }
    
    if (warn_anisotropy)
    cout << endl << endl
         << "WARNING '" << filename << "': EOutersphere does not support "
         << "tensorial polarizabilities." 
         << endl;

    return poles;
}


void EOutersphere::DistributeMpoles(Topology *top) {

    // +++++++++++++++++++++++++++++++++++++ //
    // Equip TOP with distributed multipoles //
    // +++++++++++++++++++++++++++++++++++++ //

    vector<Segment*> ::iterator sit;
    for (sit = top->Segments().begin();
         sit < top->Segments().end();
         ++sit) {

        Segment *seg = *sit;
        vector<PolarSite*> poleSites = _map_seg_polarSites.at(seg->getName());
        seg->setChrgStates(_map_seg_chrgStates[seg->getName()]);

        bool map2md = _map2md[seg->getName()];

        vector<Fragment*> ::iterator fit;
        for (fit = seg->Fragments().begin();
             fit < seg->Fragments().end();
             ++fit) {

            Fragment *frag = *fit;

            string idkey = frag->getName() + seg->getName()
                         + seg->getMolecule()->getName();
            vector<int> polesInFrag = _alloc_frag_mpoleIdx.at(idkey);
            vector<string> namesInFrag = _alloc_frag_mpoleName.at(idkey);

            if (map2md) {
                if (polesInFrag.size() != frag->Atoms().size()) {
                    cout << endl << "ERROR: Segment " << seg->getName() <<
                            " Fragment " << frag->getName() <<
                            ": MAP2MD = TRUE requires same number of polar sites "
                            " as there are atoms to perform mapping. " << endl;
                    throw runtime_error("User not paying attention. ");
                }
            }


            for (unsigned int i = 0; i < polesInFrag.size(); i++) {

                string name = namesInFrag[i];
                int poleId = polesInFrag[i];

                PolarSite *templ = poleSites[poleId-1];
                PolarSite *newSite = top->AddPolarSite(name);
                newSite->ImportFrom(templ);
                seg->AddPolarSite(newSite);
                frag->AddPolarSite(newSite);

                // Shift + rotate
                if (!map2md) {
                    newSite->Translate(frag->getTransQM2MD());
                    newSite->Rotate(frag->getRotQM2MD(), frag->getCoMD());
                }
                else {
                    vec mdpos = frag->Atoms()[i]->getPos();
                    newSite->setPos(mdpos);
                    if (newSite->getRank() > 0) {
                        cout << endl << "ERROR: MAP2MD = TRUE "
                        " prevents use of higher-rank multipoles. " << endl;
                        throw runtime_error("User not paying attention. ");
                    }
                }
            }
        }
    }
}


QMPair *EOutersphere::RequestNextPair(int opId, Topology *top) {
    
    _nextPairMutex.Lock();

    QMPair *workOnThis;

    if (_nextPair == top->NBList().end()) {
        workOnThis = NULL;
    }
    else {
        workOnThis = *_nextPair;
        ++_nextPair;
        cout << "\r... ... OP " << opId
             << " evaluating energy for pair "
             << workOnThis->getId() << flush;
    }

    _nextPairMutex.Unlock();

    return workOnThis;
}


bool EOutersphere::EvaluateFrame(Topology *top) {

    // ++++++++++++++++++++ //
    // Ridigidfy + Estatify //
    // ++++++++++++++++++++ //

    // Rigidify if (a) not rigid yet (b) rigidification at all possible
    if (!top->isRigid()) {
        bool isRigid = top->Rigidify();
        if (!isRigid) { return 0; }
    }
    else { cout << endl << "... ... System is already rigidified."; }

    // Forward multipoles to topology
    if (_method != "constant" && !top->isEStatified()) {
        this->DistributeMpoles(top);
        cout << endl << "... ... Created " << top->PolarSites().size()
                 << " multipole sites.";
    }
    else if (_method == "constant") {
        cout << endl << "... ... Method = 'constant' => Do not estatify." << flush;
        cout << endl << "... ... Using constant lambda = " << _lambdaConstant << flush;
    }
    else { cout << endl << "... ... System is already estatified."; }


    // +++++++++++++++++++++++++++++++++ //
    // Create + start threads (Site Ops) //
    // +++++++++++++++++++++++++++++++++ //

    vector<PairOpOutersphere*> pairOps;
    _nextPair = top->NBList().begin();
    cout << endl;

    // Forward to first segment specified in options
    //for ( ; (*_nextPair)->getId() != this->_firstSeg &&
    //          _nextSite < top->Segments().end(); ++_nextSite) { ; }

    for (unsigned int id = 0; id < _nThreads; id++) {
        PairOpOutersphere *newOp = new PairOpOutersphere(id, top, this);
        pairOps.push_back(newOp);
    }

    for (unsigned int id = 0; id < _nThreads; id++) {
        pairOps[id]->InitSlotData(top);
    }

    for (unsigned int id = 0; id < _nThreads; id++) {
        pairOps[id]->Start();
    }

    for (unsigned int id = 0; id < _nThreads; id++) {
        pairOps[id]->WaitDone();
    }

    for (unsigned int id = 0; id < _nThreads; id++) {
        delete pairOps[id];
    }

    pairOps.clear();

    return true;

    //if      (_method == "constant")   { this->ConstLambda(top); }
    //else if (_method == "spheres")    { this->SpheresLambda(top);}
    //else if (_method == "dielectric") { this->DielectricLambda(top); }
}


// +++++++++++++++++++++++++++ //
// PAIR OPERATOR MEMBER FUNCTS //
// +++++++++++++++++++++++++++ //

void EOutersphere::PairOpOutersphere::EvalSite(Topology *top, QMPair *qmpair) {

    if (_master->_method == "dielectric") {
        this->_segsSphere.clear();
        this->_polsSphere.clear();
        this->_polsPair.clear();

        vector<Segment*> ::iterator sit;
        for (sit = top->Segments().begin(); sit < top->Segments().end(); ++sit) {

            // TODO Introduce volume-correction when using cutoff
            //if ( abs(_top->PbShortestConnect((*sit)->getPos(), qmpair->getPos()))
            //        > _master->_lambdaCutoff) { continue; }

            if ( (*sit)->getId() == qmpair->Seg1()->getId() ||
                 (*sit)->getId() == qmpair->Seg2()->getId() ) {
                 _polsPair.push_back( _polarSites[(*sit)->getId() - 1] );
            }

            else {
                _segsSphere.push_back(*sit);
                _polsSphere.push_back( _polarSites[(*sit)->getId() - 1] );
            }
        }

        bool hasCation = qmpair->first->hasChrgState(+1)
                      && qmpair->second->hasChrgState(+1);
        bool hasAnion  = qmpair->first->hasChrgState(-1)
                      && qmpair->second->hasChrgState(-1);

        if (hasCation) {
            int state = +1;
            this->ResetFields();
            this->ChargeDeltaPair(state);
            double lOut = this->CalcOuter(top, qmpair, state);
            qmpair->setLambdaO(lOut, state);
            this->ResetFields();
        }

        if (hasAnion) {
            int state = -1;
            this->ResetFields();
            this->ChargeDeltaPair(state);
            double lOut = this->CalcOuter(top, qmpair, state);
            qmpair->setLambdaO(lOut, state);
            this->ResetFields();
        }        
    }
    else if (_master->_method == "constant") {
        
        bool hasCation = qmpair->first->hasChrgState(+1)
                      && qmpair->second->hasChrgState(+1);
        bool hasAnion  = qmpair->first->hasChrgState(-1)
                      && qmpair->second->hasChrgState(-1);
        
        if (hasCation) {
            int state = +1;
            double lOut = _master->_lambdaConstant;
            qmpair->setLambdaO(lOut, state);
        }

        if (hasAnion) {
            int state = -1;
            double lOut = _master->_lambdaConstant;
            qmpair->setLambdaO(lOut, state);
        }  
    }
    else if (_master->_method == "spheres") {
        double e = 1.602176487e-19;
        double EPS0 = 8.85418782e-12;
        double NM = 1e-09;
        // TODO extract radii from input
        _master->_radius.at("c60");
        double R1 = _master->_radius.at(qmpair->Seg1()->getName());
        double R2 = _master->_radius.at(qmpair->Seg2()->getName());
        double lambda = _master->_pekarFactor * e / (4.*M_PI*EPS0) *
                                     (  1. / ( 2.*R1*NM )
                                      + 1. / ( 2.*R2*NM )
                                      - 1. / ( qmpair->Dist()*NM ));
        qmpair->setLambdaO(lambda, -1);
        qmpair->setLambdaO(lambda, +1);
    }
    else {
        assert(false); // No such method
    }
}


double EOutersphere::PairOpOutersphere::CalcOuter(Topology *top, QMPair *qmpair, int state) {

    vector< vector<PolarSite*> > ::iterator sit1;
    vector< PolarSite* > ::iterator pit1;
    vector< vector<PolarSite*> > ::iterator sit2;
    vector< PolarSite* > ::iterator pit2;    


    double dD2dV = 0.0;
    double int2displ = 1/(4*M_PI) * 1.602176487e-19 / 1.000e-18;
    //double int2eV   = 8.854187817e-12 *int2N_C*int2N_C /1.602176487e-19;

    int polCount = 0;

    for (sit1 = _polsSphere.begin(); sit1 < _polsSphere.end(); ++sit1) {
    for (pit1 = (*sit1).begin(); pit1 < (*sit1).end(); ++pit1) {

        ++polCount;

    for (sit2 = _polsPair.begin(); sit2 < _polsPair.end(); ++sit2) {
    for (pit2 = (*sit2).begin(); pit2 < (*sit2).end(); ++pit2) {

        this->_actor.Field_AT_DUE(*(*pit1), *(*pit2));

    }}
        dD2dV += (*pit1)->FPx * (*pit1)->FPx +
                 (*pit1)->FPy * (*pit1)->FPy +
                 (*pit1)->FPz * (*pit1)->FPz;
    }}

    dD2dV *= int2displ * int2displ;

    return 0.5 * _master->_pekarFactor/8.854187817e-12
               * top->BoxVolume()*1.000e-27/top->PolarSites().size()
               * dD2dV /1.602176487e-19;
}


void EOutersphere::PairOpOutersphere::ChargeDeltaPair(int state) {

    vector< vector<PolarSite*> > ::iterator sit;
    vector< PolarSite* > ::iterator pit;

    sit = _polsPair.begin();

    for (pit = (*sit).begin(); pit < (*sit).end(); ++pit) {
        (*pit)->ChargeDelta(state, 0);
    }

    ++sit;

    for (pit = (*sit).begin(); pit < (*sit).end(); ++pit) {
        (*pit)->ChargeDelta(0, state);
    }
}


void EOutersphere::PairOpOutersphere::ResetFields() {

    vector< vector<PolarSite*> > ::iterator sit;
    vector< PolarSite* > ::iterator pit;

    for (sit = _polsSphere.begin(); sit < _polsSphere.end(); ++sit) {
        for (pit = (*sit).begin(); pit < (*sit).end(); ++pit) {
            (*pit)->ResetFieldP();
        }
    }
}


void EOutersphere::PairOpOutersphere::Run(void) {

    while (true) {

        _pair = _master->RequestNextPair(_id, _top);

        if (_pair == NULL) { break; }
        else { this->EvalSite(_top, _pair); }
    }
}


void EOutersphere::PairOpOutersphere::InitSlotData(Topology *top) {

    vector< Segment* > ::iterator sitRef;
    vector< vector<PolarSite*> > ::iterator sitNew;
    vector< PolarSite* > ::iterator pitRef;
    vector< PolarSite* > ::iterator pitNew;

    _polarSites.resize(top->Segments().size());
    assert(top->Segments().size() == _polarSites.size());

    for (sitRef = top->Segments().begin(), sitNew = _polarSites.begin();
         sitRef < top->Segments().end();
         ++sitRef, ++sitNew) {

        (*sitNew).resize((*sitRef)->PolarSites().size());

        for (pitRef = (*sitRef)->PolarSites().begin(),
             pitNew = (*sitNew).begin();
             pitRef < (*sitRef)->PolarSites().end();
             ++pitRef, ++ pitNew) {

            *pitNew = new PolarSite();
            (*pitNew)->ImportFrom(*pitRef, "full");
            (*pitNew)->Charge(0);
        }
    }
}


EOutersphere::PairOpOutersphere::~PairOpOutersphere() {

    vector< vector<PolarSite*> > ::iterator sit;
    vector<PolarSite*> ::iterator pol;

    for (sit = _polarSites.begin(); sit < _polarSites.end(); ++sit) {
        for (pol = (*sit).begin(); pol < (*sit).end(); ++pol) {
            delete *pol;
        }
        (*sit).clear();
    }
    _polarSites.clear();
}


// +++++++++++++++++++++++++++ //
// CHEAPER METHODS (DISCARDED) //
// +++++++++++++++++++++++++++ //

void EOutersphere::ConstLambda(Topology *top) {

    QMNBList ::iterator pit;
    for (pit = top->NBList().begin(); pit != top->NBList().end(); pit++) {
        QMPair *qpair = *pit;
        qpair->setLambdaO(_lambdaConstant, -1);
        qpair->setLambdaO(_lambdaConstant, +1);
    }
}


void EOutersphere::SpheresLambda(Topology *top) {
    
    //double e = 1.602176487e-19;
    //double EPS0 = 8.85418782e-12;
    //double NM = 1e-09;    
    
    QMNBList ::iterator pit;
    for (pit = top->NBList().begin(); pit != top->NBList().end(); pit++) {
        //QMPair *pair = *pit;

        // TODO extract radii from input
        throw std::runtime_error(" EOutersphere -> Need to read in radii");
        //double R1 = 1;
        //double R2 = 2;
        /*double lambda = _pekarFactor * e / (4.*M_PI*EPS0) *
                                     (  1. / ( 2.*R1*NM )
                                      + 1. / ( 2.*R2*NM )
                                      - 1. / ( pair->Dist()*NM ));*/
        assert (false); // TODO e, h calc. // pair->setLambdaO(lambda);
    }
}


void EOutersphere::DielectricLambda(Topology *top) {

    double e = 1.602176487e-19;
    double EPS0 = 8.85418782e-12;
    double NM = 1.e-09;
    double NM3 = 1.e-27;

    QMNBList ::iterator pit;
    for (pit = top->NBList().begin(); pit != top->NBList().end(); pit++) {

        QMPair *pair = *pit;
        double lambda = 0.;
        vec seg1pos = pair->Seg1()->getPos();
        vec seg2pos = pair->Seg2()->getPos();

        vector< Segment* > ::iterator sit;
        for (sit = top->Segments().begin();
             sit < top->Segments().end();
             sit++) {
             Segment *seg = *sit;

             // Segment external to pair? Cut-off obeyed?
             if ( seg == pair->Seg1() || seg == pair->Seg2() ) { continue; }

             vec dR1bc = top->PbShortestConnect( seg->getPos(), seg1pos );
             vec dR2bc = top->PbShortestConnect( seg->getPos(), seg2pos );
             if ( abs( (dR1bc+dR2bc)*0.5 ) > _lambdaCutoff ) { continue; }

             vec dR1   = seg1pos - seg->getPos();
             vec dR2   = seg2pos - seg->getPos();
             vec shift1 = dR1bc - dR1;
             vec shift2 = dR2bc - dR2;

             vector< Atom* > ::iterator ait;
             vector< Atom* > ::iterator bit;
             for (ait = seg->Atoms().begin(); ait < seg->Atoms().end(); ait++) {

                 Atom *Ext = *ait;
                 // Electric induction field
                 vec D = vec(0.,0.,0.);

                 for (bit = pair->Seg1()->Atoms().begin();
                         bit < pair->Seg1()->Atoms().end();
                         bit++) {

                     Atom *Int = *bit;
                     
                     double dQ = Int->getQ(-1) - Int->getQ(0);
                     
                     vec R = Ext->getPos() - Int->getPos() - shift1;
                     double dr = abs(R);
                     double dr3 = dr*dr*dr;

                     D += e * dQ / (4.*M_PI*dr3*NM3) * R * NM;

                 }

                 for (bit = pair->Seg2()->Atoms().begin();
                         bit < pair->Seg2()->Atoms().end();
                         bit++) {
                     Atom *Int = *bit;

                     double dQ = Int->getQ(0) - Int->getQ(-1);
                     
                     vec R = Ext->getPos() - Int->getPos() - shift2;
                     double dr = abs(R);
                     double dr3 = dr*dr*dr;

                     D += e * dQ / (4.*M_PI*dr3*NM3) * R * NM;

                 }
                 
                 lambda += 1/e * _pekarFactor/(2*EPS0)
                         * top->BoxVolume()*NM3 / top->Atoms().size();   
             } /* exit loop over external atoms in segment */
        } /* exit loop over external segments */

        assert (false); // TODO e, h calculation // pair->setLambdaO(lambda);

    } /* exit loop over pairs */

}


// +++++++++++++++++++++++++++ //
// INTERACTOR MEMBER FUNCTIONS //
// +++++++++++++++++++++++++++ //

void EOutersphere::InteractorMod::Field_AT_DUE(PolarSite &pol1, PolarSite &pol) {

    // NOTE >>> e12 points from polar site 1 to polar site 2 <<< NOTE //
    e12  = _top->PbShortestConnect(pol1.getPos(), pol.getPos());
    R    = 1/abs(e12);
    R2   = R*R;
    R3   = R2*R;
    R4   = R3*R;
    R5   = R4*R;
    e12 *= R;

        rax = e12.getX();
        ray = e12.getY();
        raz = e12.getZ();
        rbx = - rax;
        rby = - ray;
        rbz = - raz;

        cxx = 1;
        cxy = 0;
        cxz = 0;
        cyx = 0;
        cyy = 1;
        cyz = 0;
        czx = 0;
        czy = 0;
        czz = 1;

    double Fx = 0.0;
    double Fy = 0.0;
    double Fz = 0.0;

    // Field generated by rank-0 m'pole
        Fx += T1x_00() * pol.Q00;
        Fy += T1y_00() * pol.Q00;
        Fz += T1z_00() * pol.Q00;

    // Field generated by rank-1 m'pole
    if (pol._rank > 0) {
        Fx += T1x_1x() * pol.Q1x;
        Fx += T1x_1y() * pol.Q1y;
        Fx += T1x_1z() * pol.Q1z;

        Fy += T1y_1x() * pol.Q1x;
        Fy += T1y_1y() * pol.Q1y;
        Fy += T1y_1z() * pol.Q1z;

        Fz += T1z_1x() * pol.Q1x;
        Fz += T1z_1y() * pol.Q1y;
        Fz += T1z_1z() * pol.Q1z;
    }

    // Field generated by rank-2 m'pole
    if (pol._rank > 1) {
        Fx += T1x_20()  * pol.Q20;
        Fx += T1x_21c() * pol.Q21c;
        Fx += T1x_21s() * pol.Q21s;
        Fx += T1x_22c() * pol.Q22c;
        Fx += T1x_22s() * pol.Q22s;

        Fy += T1y_20()  * pol.Q20;
        Fy += T1y_21c() * pol.Q21c;
        Fy += T1y_21s() * pol.Q21s;
        Fy += T1y_22c() * pol.Q22c;
        Fy += T1y_22s() * pol.Q22s;

        Fz += T1z_20()  * pol.Q20;
        Fz += T1z_21c() * pol.Q21c;
        Fz += T1z_21s() * pol.Q21s;
        Fz += T1z_22c() * pol.Q22c;
        Fz += T1z_22s() * pol.Q22s;
    }

    pol1.FPx += Fx;
    pol1.FPy += Fy;
    pol1.FPz += Fz;
}



}}



