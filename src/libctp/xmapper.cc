#include <votca/ctp/xmapper.h>


namespace votca { namespace ctp {

    
void XMpsMap::GenerateMap(string xml_file, 
                          string alloc_table, 
                          Topology *top, 
                          vector<XJob*> &xjobs) {
    
    this->CollectMapFromXML(xml_file);
    this->CollectSegMpsAlloc(alloc_table, top);
    this->CollectSitesFromMps(xjobs);
    
}
    
    
void XMpsMap::CollectMapFromXML(string xml_file) {

    cout << endl 
         << "... ... ... Allocate polar sites to fragments. "
         << flush;

    string allocFile = xml_file;

    // ++++++++++++++++++++++++++++++++ //
    // Load polar-site indices from XML //
    // ++++++++++++++++++++++++++++++++ //

    // => Output to maps:
    map<string, vector<int> > alloc_frag_mpoleIdx;
    map<string, vector<string> > alloc_frag_mpoleName;
    map<string, vector<int> > alloc_frag_trihedron;
    map<string, vector<double> > alloc_frag_weights;
    map<string, vector<int> > alloc_frag_isVirtualMp;

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


    string key = "topology.molecules.molecule";
    list<Property *> mols = allocation.Select(key);
    list<Property *>::iterator molit;
    for (molit = mols.begin(); molit != mols.end(); molit++) {

        string molName = (*molit)->get("name").as<string> ();

        key = "segments.segment";
        list<Property *> segs = (*molit)->Select(key);
        list<Property *>::iterator segit;

        for (segit = segs.begin(); segit != segs.end(); segit++) {

            string segName = (*segit)->get("name").as<string> ();

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

                // Local frame for polar sites
                vector<int> trihedron_mps;
                if ((*fragit)->exists("localframe_mps")) {
                   cout << endl
                         << "... ... ... ... " << segName << ": "
                         << "Defining distinct local frame for polar sites."
                         << flush;
                   trihedron_mps = (*fragit)->get("localframe_mps")
                                         .as< vector<int> >();
                }
                else {
                   trihedron_mps = (*fragit)->get("localframe")
                                         .as< vector<int> >();
                }
                
                // Mapping weights for polar sites
                vector<double> weights_mps;
                if ((*fragit)->exists("weights_mps")) {
                    cout << endl
                         << "... ... ... ... " << segName << ": "
                         << "Using distinct weights for polar sites."
                         << flush;
                   weights_mps = (*fragit)->get("weights_mps")
                                       .as< vector<double> >();
                }
                else {
                   weights_mps = (*fragit)->get("weights")
                                       .as< vector<double> >();
                }
                
                // Virtual vs real (= atom-centered) polar sites
                vector<int> isVirtualMp;
                if ((*fragit)->exists("virtual_mps")) {
                    cout << endl
                         << "... ... ... ... " << segName << ": "
                         << "Checking for virtual expansion sites."
                         << flush;
                    vector<int> isVirtual_ints = (*fragit)->get("virtual_mps")
                                        .as< vector<int> >();
                    for (int i = 0; i < isVirtual_ints.size(); ++i) {
                        bool isVirtual = (isVirtual_ints[i] == 0) ? false:true;
                        isVirtualMp.push_back(isVirtual);
                    }
                }
                else {
                    for (int i = 0; i < weights_mps.size(); ++i) {
                        isVirtualMp.push_back(false);
                    }
                }

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

                if ( (mpoleIdcs.size() != mpoleNames.size())
                   ||(mpoleIdcs.size() != weights_mps.size()) 
                   ||(mpoleIdcs.size() != isVirtualMp.size()) ) {
                    
                    cout << endl 
                         << "ERROR: Bad multipole definition in fragment '"
                         << mapKeyName
                         << "'. Misscounted?"
                         << endl;
                    throw std::runtime_error("Revise input.");                    
                }
                
                alloc_frag_mpoleIdx[mapKeyName]         = mpoleIdcs;
                alloc_frag_mpoleName[mapKeyName]        = mpoleNames;
                alloc_frag_trihedron[mapKeyName]        = trihedron_mps;
                alloc_frag_weights[mapKeyName]          = weights_mps;
                alloc_frag_isVirtualMp[mapKeyName]      = isVirtualMp;
            }
        }
    }

    _alloc_frag_mpoleIdx    = alloc_frag_mpoleIdx;
    _alloc_frag_mpoleName   = alloc_frag_mpoleName;
    _alloc_frag_trihedron   = alloc_frag_trihedron;
    _alloc_frag_weights     = alloc_frag_weights;
    _alloc_frag_isVirtualMp = alloc_frag_isVirtualMp;
}
    
    
void XMpsMap::CollectSegMpsAlloc(string alloc_table, Topology *top) {

    _alloc_table = alloc_table;
    
    std::string line;
    std::ifstream intt;
    intt.open(alloc_table.c_str());

    if (intt.is_open() ) {
        while ( intt.good() ) {

            std::getline(intt, line);
            vector<string> split;
            Tokenizer toker(line, " ");
            toker.ToVector(split);

            if ( !split.size()      ||
                  split[0] == "#"   ||
                  split[0].substr(0,1) == "#" ) { continue; }

            // Line sample:
            // 1022 C60 C60_n.mps C60_e.mps C60_h.mps

            int segId       = boost::lexical_cast<int>(split[0]);
            string segName  = boost::lexical_cast<string>(split[1]);

            // Some compliance checks
            assert( split.size() == 5 );                                        // <- File format correct?
            assert( top->getSegment(segId)->getName() == segName );             // <- Input matches topology?

            _segId_mpsFile_n[segId] = split[2];
            _segId_mpsFile_e[segId] = split[3];
            _segId_mpsFile_h[segId] = split[4];

        } // Exit loop over lines
    }
    else { cout << endl << "ERROR: No such file " << alloc_table << endl;
           throw std::runtime_error("Input file missing.");               }

    assert( _segId_mpsFile_n.size() == top->Segments().size() );                // <- Input for all segments?

}


void XMpsMap::CollectSitesFromMps(vector<XJob*> &xjobs) {

    // +++++++++++++ //
    // Parse + Store //
    // +++++++++++++ //

    map<int,string> ::iterator misit;

    // Neutral species
    for (misit = _segId_mpsFile_n.begin();
         misit != _segId_mpsFile_n.end();
         ++misit) {
        if (_mpsFile_pSites.count(misit->second) > 0 ) { continue; }
        else { _mpsFile_pSites[misit->second] = Parse_GDMA(misit->second, 0); }
    }

    // Anion
    for (misit = _segId_mpsFile_e.begin();
         misit != _segId_mpsFile_e.end();
         ++misit) {
        if (_mpsFile_pSites.count(misit->second) > 0 ) { continue; }
        else { _mpsFile_pSites[misit->second] = Parse_GDMA(misit->second, -1); }
    }

    // Cation
    for (misit = _segId_mpsFile_h.begin();
         misit != _segId_mpsFile_h.end();
         ++misit) {
        if (_mpsFile_pSites.count(misit->second) > 0 ) { continue; }
        else { _mpsFile_pSites[misit->second] = Parse_GDMA(misit->second, +1); }
    }
    
    // Job Seg1 Seg2
    vector<XJob*> :: iterator jit;
    for (jit = xjobs.begin();
         jit < xjobs.end();
         ++jit) {

        vector<string> mpsfiles = (*jit)->getSegMps();
        
        for (int i = 0; i < mpsfiles.size(); ++i) {
            string mpsfile = mpsfiles[i];
            if (_mpsFile_pSites_job.count(mpsfile) > 0 ) { ; }
            else { _mpsFile_pSites_job[mpsfile] = Parse_GDMA(mpsfile, 0); }
        }
    }

    cout << endl
         << "... ... ..."
         << " Parsed " << _mpsFile_pSites.size()
         << " mps-files from " << _alloc_table
         << ", " << _mpsFile_pSites_job.size()
         << " mps-files from XJobs."
         << flush;
}


vector<APolarSite*> XMpsMap::Parse_GDMA(string filename, int state) {

    return APS_FROM_MPS(filename,state);
    
}


void XMpsMap::EquipWithPolSites(Topology *top) {
    
    // Warning: Direct mapping of k > 0 multipoles to MD coordinates
    bool print_huge_map2md_warning = false;
    
    // Log warning: Symmetry = 1 and k > 0 multipoles.
    map<string,bool> warned_symm_idkey;

    // +++++++++++++++++++++++++++++++++++++ //
    // Equip TOP with distributed multipoles //
    // +++++++++++++++++++++++++++++++++++++ //

    vector<Segment*> ::iterator sit;
    for (sit = top->Segments().begin();
         sit < top->Segments().end();
         ++sit) {

        Segment *seg                = *sit;
        int segId                   = seg->getId();

        bool map2md                 = _map2md[seg->getName()];

        string mps_n                = _segId_mpsFile_n[segId];
        string mps_h                = _segId_mpsFile_h[segId];
        string mps_e                = _segId_mpsFile_e[segId];

        vector<APolarSite*> pols_n   = _mpsFile_pSites[mps_n];
        vector<APolarSite*> pols_h   = _mpsFile_pSites[mps_h];
        vector<APolarSite*> pols_e   = _mpsFile_pSites[mps_e];

        // Merge polar sites
        assert(pols_n.size() == pols_h.size());
        assert(pols_n.size() == pols_e.size());

        for (int i = 0; i < pols_n.size(); i++) {
            pols_n[i]->setQs( pols_h[i]->getQs(+1), +1 );
            pols_n[i]->setPs( pols_h[i]->getPs(+1), +1 );
        }
        for (int i = 0; i < pols_n.size(); ++i) {
            pols_n[i]->setQs(pols_e[i]->getQs(-1), -1 );
            pols_n[i]->setPs(pols_e[i]->getPs(-1), -1 );
        }

        vector<Fragment*> ::iterator fit;
        for (fit = seg->Fragments().begin();
             fit < seg->Fragments().end();
             ++fit) {

            Fragment *frag = *fit;

            // Retrieve polar-site data for this fragment
            string idkey                 = frag->getName() + seg->getName()
                                         + seg->getMolecule()->getName();
            vector<int> polesInFrag      = _alloc_frag_mpoleIdx.at(idkey);
            vector<string> namesInFrag   = _alloc_frag_mpoleName.at(idkey);
            vector<double> weightsInFrag = _alloc_frag_weights.at(idkey);
            vector<int> isVirtualMp      = _alloc_frag_isVirtualMp.at(idkey);

            if (map2md && polesInFrag.size() != frag->Atoms().size()) {
                cout << endl
                     << "ERROR: Segment " << seg->getName()
                     << " Fragment " << frag->getName()
                     << ": MAP2MD = TRUE requires same number of polar "
                     << "sites as there are atoms to perform mapping. "
                     << endl;
                throw runtime_error("Check mapping or switch to map2md = 0");
            }

            matrix rotateMP2MD;
            vec translateMP2MD;

            // Determine transformation matrices to execute mapping
            if (!map2md) {

                vector<APolarSite*> trihedron_pol;
                vector<Atom*>      trihedron_atm;

                vector<int> trihedron_ints  = _alloc_frag_trihedron.at(idkey);
                vector<int> ::iterator iit;
                for (iit = trihedron_ints.begin();
                     iit < trihedron_ints.end();
                     ++iit) {
                    trihedron_pol.push_back(pols_n[(*iit)-1]);
                }
                
                trihedron_ints  = frag->getTrihedron();
                for (iit = trihedron_ints.begin();
                     iit < trihedron_ints.end();
                     ++iit) {
                    vector< Atom* > ::iterator ait;
                    for (ait = frag->Atoms().begin();
                         ait < frag->Atoms().end();
                         ++ait) {
                        if ((*ait)->getQMId() == (*iit)) {
                            trihedron_atm.push_back(*ait);
                        }
                    }
                }


                int symmetry = trihedron_pol.size();
                assert (trihedron_pol.size() <= trihedron_atm.size() );               

                vec xMD, yMD, zMD;
                vec xQM, yQM, zQM;

                if (symmetry == 3) {
                    vec r1MD = trihedron_atm[0]->getPos();
                    vec r2MD = trihedron_atm[1]->getPos();
                    vec r3MD = trihedron_atm[2]->getPos();
                    vec r1QM = trihedron_pol[0]->getPos();
                    vec r2QM = trihedron_pol[1]->getPos();
                    vec r3QM = trihedron_pol[2]->getPos();

                    xMD = r2MD - r1MD;
                    yMD = r3MD - r1MD;
                    xQM = r2QM - r1QM;
                    yQM = r3QM - r1QM;

                    zMD = xMD ^ yMD;
                    zQM = xQM ^ yQM;

                    yMD = zMD ^ xMD;
                    yQM = zQM ^ xQM;

                    xMD = xMD.normalize();
                    yMD = yMD.normalize();
                    zMD = zMD.normalize();
                    xQM = xQM.normalize();
                    yQM = yQM.normalize();
                    zQM = zQM.normalize();
                }

                else if (symmetry == 2) {

                    vec r1MD = trihedron_atm[0]->getPos();
                    vec r2MD = trihedron_atm[1]->getPos();
                    vec r1QM = trihedron_pol[0]->getPos();
                    vec r2QM = trihedron_pol[1]->getPos();

                    xMD = r2MD - r1MD;
                    xQM = r2QM - r1QM;

                    // Normalising not necessary, but recommendable, when doing...
                    xMD = xMD.normalize();
                    xQM = xQM.normalize();

                    vec yMDtmp = vec(0,0,0);
                    vec yQMtmp = vec(0,0,0);

    // ... this: Check whether one of the components is equal or close to
    // zero. If so, this easily gives a second leg for the trihedron.
    if      ( xMD.getX()*xMD.getX() < 1e-6 ) { yMDtmp = vec(1,0,0); }
    else if ( xMD.getY()*xMD.getY() < 1e-6 ) { yMDtmp = vec(0,1,0); }
    else if ( xMD.getZ()*xMD.getZ() < 1e-6 ) { yMDtmp = vec(0,0,1); }
    if      ( xQM.getX()*xQM.getX() < 1e-6 ) { yQMtmp = vec(1,0,0); }
    else if ( xQM.getY()*xQM.getY() < 1e-6 ) { yQMtmp = vec(0,1,0); }
    else if ( xQM.getZ()*xQM.getZ() < 1e-6 ) { yQMtmp = vec(0,0,1); }

    if ( abs(yMDtmp) < 0.5 ) {
       // All components of xMD are unequal to zero => division is safe.
       // Choose vector from plane with xMD * inPlaneVec = 0:
       double tmp_x = 1.;
       double tmp_y = 1.;
       double tmp_z = 1/xMD.getZ() * (-xMD.getX()*tmp_x - xMD.getY()*tmp_y);
       yMDtmp = vec(tmp_x, tmp_y, tmp_z);
       yMDtmp.normalize();
    }
    if ( abs(yQMtmp) < 0.5 ) {
       double tmp_x = 1.;
       double tmp_y = 1.;
       double tmp_z = 1/xQM.getZ() * (-xQM.getX()*tmp_x - xQM.getY()*tmp_y);
       yQMtmp = vec(tmp_x, tmp_y, tmp_z);
       yQMtmp.normalize();
    }

                    // Now proceed as for symmetry 3
                    zMD = xMD ^ yMDtmp;
                    yMD = zMD ^ xMD;
                    zQM = xQM ^ yQMtmp;
                    yQM = zQM ^ xQM;

                    xMD.normalize();
                    yMD.normalize();
                    zMD.normalize();
                    xQM.normalize();
                    yQM.normalize();
                    zQM.normalize();
                }

                else if (symmetry == 1) {

                    if (!warned_symm_idkey[idkey]) {
                        cout << endl << "... ... ... "
                         << "WARNING: Symmetry = 1 for fragment "
                         << frag->getName() << ": This will generate artifacts "
                         << "when mapping higher-rank multipoles (dipoles, ..)."
                         << flush;
                        warned_symm_idkey[idkey] = true;
                    }

                    xMD = vec(1,0,0);
                    yMD = vec(0,1,0);
                    zMD = vec(0,0,1);
                    xQM = vec(1,0,0);
                    yQM = vec(0,1,0);
                    zQM = vec(0,0,1);
                }

                else {
                    cout << endl
                         << "NOTE: Invalid definition of local frame in fragment "
                         << frag->getName();
                    cout << ". Assuming point particle for mapping. "
                         << endl;
                    cout << endl
                         << "WARNING: Symmetry = 1 for fragment "
                         << frag->getName() << ": This will generate artifacts "
                         << "when mapping higher-rank multipoles (dipoles, ..)."
                         << endl;

                    xMD = vec(1,0,0);
                    yMD = vec(0,1,0);
                    zMD = vec(0,0,1);
                    xQM = vec(1,0,0);
                    yQM = vec(0,1,0);
                    zQM = vec(0,0,1);
                }

                matrix rotMD = matrix(xMD, yMD, zMD);
                matrix rotMP = matrix(xQM, yQM, zQM);
                
                rotateMP2MD = rotMD * rotMP.Transpose();


                // ++++++++++++++++++ //
                // Transform fragment //
                // ++++++++++++++++++ //

                
                vec CoMP = vec(0.,0.,0.);
                double W = 0.0;
                for (int i = 0; i < polesInFrag.size(); ++i) {

                    double weight = weightsInFrag[i];
                    
                    vec pos = pols_n[polesInFrag[i]-1]->getPos();
                    
                    CoMP += weight*pos;
                    W += weight;

                }
                CoMP /= W;

                translateMP2MD = frag->getCoMD() - CoMP;

            }            

            // Create polar sites 
            for (int i = 0; i < polesInFrag.size(); i++) {

                string name             = namesInFrag[i];
                int poleId              = polesInFrag[i];

                APolarSite *templ        = pols_n[poleId-1];
                APolarSite *newSite      = top->AddAPolarSite(name);
                newSite->ImportFrom(templ);
                newSite->setIsVirtual(isVirtualMp[i]);
                seg->AddAPolarSite(newSite);
                frag->AddAPolarSite(newSite);

                // Shift + rotate
                if (!map2md) {
                    newSite->Translate(translateMP2MD);
                    newSite->Rotate(rotateMP2MD, frag->getCoMD());
                }
                else {
                    vec mdpos = frag->Atoms()[i]->getPos();
                    newSite->setPos(mdpos);
                    if (newSite->getRank() > 0) {
                        print_huge_map2md_warning = true;
                    }
                }
            }
        }
    }

    if (print_huge_map2md_warning) {
        cout << endl << endl
             << "**************************************************************"
             << "WARNING: MAP2MD = TRUE while using higher-rank multipoles can "
             << "mess up the orientation of those multipoles if the coordinate "
             << "frame used in the .mps file does not agree with the global MD "
             << "frame. If you know what you are doing - proceed ... "
             << "**************************************************************"
             << endl;
    }

    top->setIsEStatified(true);

}


vector<APolarSite*> XMpsMap::MapPolSitesToSeg(const vector<APolarSite*> &pols_n, Segment *seg) {

    bool print_huge_map2md_warning = false;

    vector<APolarSite*> return_pols;
    return_pols.reserve(pols_n.size());

    int segId                   = seg->getId();
    bool map2md                 = _map2md[seg->getName()];

    vector<Fragment*> ::iterator fit;
    for (fit = seg->Fragments().begin();
         fit < seg->Fragments().end();
         ++fit) {

        Fragment *frag = *fit;

        // Retrieve polar-site data for this fragment
        string idkey                = frag->getName() + seg->getName()
                                    + seg->getMolecule()->getName();
        vector<int> polesInFrag     = _alloc_frag_mpoleIdx.at(idkey);
        vector<string> namesInFrag  = _alloc_frag_mpoleName.at(idkey);
        vector<double> weightsInFrag= _alloc_frag_weights.at(idkey);

        if (map2md && polesInFrag.size() != frag->Atoms().size()) {
            cout << endl
                 << "ERROR: Segment " << seg->getName()
                 << " Fragment " << frag->getName()
                 << ": MAP2MD = TRUE requires same number of polar "
                 << "sites as there are atoms to perform mapping. "
                 << endl;
            throw runtime_error("Check mapping or switch to map2md = 0");
        }

        matrix rotateMP2MD;
        vec translateMP2MD;

        // Determine transformation matrices to execute mapping
        if (!map2md) {

            vector<APolarSite*> trihedron_pol;
            vector<Atom*>      trihedron_atm;

            vector<int> trihedron_ints  = _alloc_frag_trihedron.at(idkey);
            vector<int> ::iterator iit;
            for (iit = trihedron_ints.begin();
                 iit < trihedron_ints.end();
                 ++iit) {
                trihedron_pol.push_back(pols_n[(*iit)-1]);
            }

            trihedron_ints  = frag->getTrihedron();
            for (iit = trihedron_ints.begin();
                 iit < trihedron_ints.end();
                 ++iit) {
                vector< Atom* > ::iterator ait;
                for (ait = frag->Atoms().begin();
                     ait < frag->Atoms().end();
                     ++ait) {
                    if ((*ait)->getQMId() == (*iit)) {
                        trihedron_atm.push_back(*ait);
                    }
                }
            }


            int symmetry = trihedron_pol.size();
            assert (trihedron_pol.size() == trihedron_atm.size() );

            vec xMD, yMD, zMD;
            vec xQM, yQM, zQM;

            if (symmetry == 3) {
                vec r1MD = trihedron_atm[0]->getPos();
                vec r2MD = trihedron_atm[1]->getPos();
                vec r3MD = trihedron_atm[2]->getPos();
                vec r1QM = trihedron_pol[0]->getPos();
                vec r2QM = trihedron_pol[1]->getPos();
                vec r3QM = trihedron_pol[2]->getPos();

                xMD = r2MD - r1MD;
                yMD = r3MD - r1MD;
                xQM = r2QM - r1QM;
                yQM = r3QM - r1QM;

                zMD = xMD ^ yMD;
                zQM = xQM ^ yQM;

                yMD = zMD ^ xMD;
                yQM = zQM ^ xQM;

                xMD = xMD.normalize();
                yMD = yMD.normalize();
                zMD = zMD.normalize();
                xQM = xQM.normalize();
                yQM = yQM.normalize();
                zQM = zQM.normalize();
            }

            else if (symmetry == 2) {

                vec r1MD = trihedron_atm[0]->getPos();
                vec r2MD = trihedron_atm[1]->getPos();
                vec r1QM = trihedron_pol[0]->getPos();
                vec r2QM = trihedron_pol[1]->getPos();

                xMD = r2MD - r1MD;
                xQM = r2QM - r1QM;

                // Normalising not necessary, but recommendable, when doing...
                xMD = xMD.normalize();
                xQM = xQM.normalize();

                vec yMDtmp = vec(0,0,0);
                vec yQMtmp = vec(0,0,0);

    // ... this: Check whether one of the components is equal or close to
    // zero. If so, this easily gives a second leg for the trihedron.
    if      ( xMD.getX()*xMD.getX() < 1e-6 ) { yMDtmp = vec(1,0,0); }
    else if ( xMD.getY()*xMD.getY() < 1e-6 ) { yMDtmp = vec(0,1,0); }
    else if ( xMD.getZ()*xMD.getZ() < 1e-6 ) { yMDtmp = vec(0,0,1); }
    if      ( xQM.getX()*xQM.getX() < 1e-6 ) { yQMtmp = vec(1,0,0); }
    else if ( xQM.getY()*xQM.getY() < 1e-6 ) { yQMtmp = vec(0,1,0); }
    else if ( xQM.getZ()*xQM.getZ() < 1e-6 ) { yQMtmp = vec(0,0,1); }

    if ( abs(yMDtmp) < 0.5 ) {
    // All components of xMD are unequal to zero => division is safe.
    // Choose vector from plane with xMD * inPlaneVec = 0:
    double tmp_x = 1.;
    double tmp_y = 1.;
    double tmp_z = 1/xMD.getZ() * (-xMD.getX()*tmp_x - xMD.getY()*tmp_y);
    yMDtmp = vec(tmp_x, tmp_y, tmp_z);
    yMDtmp.normalize();
    }
    if ( abs(yQMtmp) < 0.5 ) {
    double tmp_x = 1.;
    double tmp_y = 1.;
    double tmp_z = 1/xQM.getZ() * (-xQM.getX()*tmp_x - xQM.getY()*tmp_y);
    yQMtmp = vec(tmp_x, tmp_y, tmp_z);
    yQMtmp.normalize();
    }

                // Now proceed as for symmetry 3
                zMD = xMD ^ yMDtmp;
                yMD = zMD ^ xMD;
                zQM = xQM ^ yQMtmp;
                yQM = zQM ^ xQM;

                xMD.normalize();
                yMD.normalize();
                zMD.normalize();
                xQM.normalize();
                yQM.normalize();
                zQM.normalize();
            }

            else if (symmetry == 1) {

                //cout << endl
                //     << "WARNING: Symmetry = 1 for fragment "
                //     << frag->getName() << ": This will generate artifacts "
                //     << "when mapping higher-rank multipoles (dipoles, ..)."
                //     << endl;

                xMD = vec(1,0,0);
                yMD = vec(0,1,0);
                zMD = vec(0,0,1);
                xQM = vec(1,0,0);
                yQM = vec(0,1,0);
                zQM = vec(0,0,1);
            }

            else {
                cout << endl
                     << "NOTE: Invalid definition of local frame in fragment "
                     << frag->getName();
                cout << ". Assuming point particle for mapping. "
                     << endl;
                cout << endl
                     << "WARNING: Symmetry = 1 for fragment "
                     << frag->getName() << ": This will generate artifacts "
                     << "when mapping higher-rank multipoles (dipoles, ..)."
                     << endl;

                xMD = vec(1,0,0);
                yMD = vec(0,1,0);
                zMD = vec(0,0,1);
                xQM = vec(1,0,0);
                yQM = vec(0,1,0);
                zQM = vec(0,0,1);
            }

            matrix rotMD = matrix(xMD, yMD, zMD);
            matrix rotMP = matrix(xQM, yQM, zQM);

            rotateMP2MD = rotMD * rotMP.Transpose();


            // ++++++++++++++++++ //
            // Transform fragment //
            // ++++++++++++++++++ //


            vec CoMP = vec(0.,0.,0.);
            double W = 0.0;
            for (int i = 0; i < polesInFrag.size(); ++i) {

                double weight = weightsInFrag[i];

                vec pos = pols_n[polesInFrag[i]-1]->getPos();

                CoMP += weight*pos;
                W += weight;

            }
            CoMP /= W;

            translateMP2MD = frag->getCoMD() - CoMP;

        }

        // Create polar sites
        for (int i = 0; i < polesInFrag.size(); i++) {

            string name             = namesInFrag[i];
            int poleId              = polesInFrag[i];

            APolarSite *templ        = pols_n[poleId-1];
            APolarSite *newSite      = new APolarSite(-1, name);
            newSite->ImportFrom(templ);

            // Shift + rotate
            if (!map2md) {
                newSite->Translate(translateMP2MD);
                newSite->Rotate(rotateMP2MD, frag->getCoMD());
            }
            else {
                vec mdpos = frag->Atoms()[i]->getPos();
                newSite->setPos(mdpos);
                if (newSite->getRank() > 0) {
                    print_huge_map2md_warning = true;
                }
            }

            newSite->Charge(0);
            
            return_pols.push_back(newSite);

        }
    } // End loop over fragments

    if (print_huge_map2md_warning) {
        cout << endl << endl
             << "**************************************************************"
             << "WARNING: MAP2MD = TRUE while using higher-rank multipoles can "
             << "mess up the orientation of those multipoles if the coordinate "
             << "frame used in the .mps file does not agree with the global MD "
             << "frame. If you know what you are doing - proceed ... "
             << "**************************************************************"
             << endl;
    }

    return return_pols;
}


}}