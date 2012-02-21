#ifndef EMULTIPOLE2_H
#define EMULTIPOLE2_H


#include <votca/ctp/qmcalculator2.h>

namespace votca { namespace ctp {


class EMultipole2 : public QMCalculator2
{

public:

    EMultipole2() {};
   ~EMultipole2() {};

    string   Identify() { return "EMultipole (Parallel)"; }

    void     Initialize(Topology *top, Property *options);
    void     EStatify(Topology *top, Property *options);
    vector<PolarSite*> ParseGdmaFile(string filename, int state);

    bool     EvaluateFrame(Topology *top);
    void     EvalSite(Topology *top, Segment *seg, int slot);
    void     Induce() { return; }
    void     CalcIntEnergy() { return; }
    
    void     InitSlotData(Topology *top) { ; }
    void     PostProcess(Topology *top) { ; }
    Segment *RequestNextSite(int opId, Topology *top);
    void     LockCout() { _coutMutex.Lock(); }
    void     UnlockCout() { _coutMutex.Unlock(); }


    // ++++++++++++++++++++++++++++++++++++++ //
    // Site workers (i.e. individual threads) //
    // ++++++++++++++++++++++++++++++++++++++ //

    class SiteOperator : public Thread
    {
    public:

        SiteOperator(int id, Topology *top,
                     EMultipole2 *master)
                   : _id(id), _top(top), _seg(NULL),
                     _master(master)      {};

       ~SiteOperator() {};

        int  getId() { return _id; }
        void setId(int id) { _id = id; }

        void Run(void);


    protected:

        int                      _id;
        Topology                *_top;
        Segment                 *_seg;
        EMultipole2             *_master;
    };


private:


    // Allocation of polar sites to fragments and segments
    map<string, vector<PolarSite*> >     _map_seg_polarSites;
    map<string, vector<int> >            _alloc_frag_mpoleIdx;
    map<string, vector<string> >         _alloc_frag_mpoleName;

    // Thread management
    vector<Segment*> ::iterator _nextSite;
    Mutex                       _nextSiteMutex;
    Mutex                       _coutMutex;

};



void EMultipole2::Initialize(Topology *top, Property *options) {

    cout << endl << "... ... Initialize with " << _nThreads << " threads.";
    cout << endl << "... ... Parameterize model";

    if (!top->isEStatified()) { this->EStatify(top, options); }



}


void EMultipole2::EStatify(Topology *top, Property *options) {

    cout << endl << "... ... Estatify system";

    string key = "options.emultipole";
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
            if ( (*segit)->exists("multipole_e")) {
                string mpoleFile_e = (*segit)->get("multipoles_e").as<string>();
                alloc_seg_mpoleFiles_e[segName] = mpoleFile_e;
            }
            if ( (*segit)->exists("multipole_H")) {
                string mpoleFile_h = (*segit)->get("multipoles_h").as<string>();
                alloc_seg_mpoleFiles_h[segName] = mpoleFile_h;
            }                       

            key = "fragments.fragment";
            list<Property *> frags = (*segit)->Select(key);
            list<Property *>::iterator fragit;

            for (fragit = frags.begin(); fragit != frags.end(); fragit++) {

                string fragName = (*fragit)->get("name").as<string> ();
                string mapKeyName = fragName + segName + molName;

                string mpoles = (*fragit)->get("mpoles").as<string> ();

                Tokenizer tokPoles(mpoles, " \t");
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


    int state = 0; // <- neutral
    map<string, string > ::iterator strit;
    for (strit = alloc_seg_mpoleFiles_n.begin();
         strit != alloc_seg_mpoleFiles_n.end();
         strit++) {

        string segName = strit->first;
        string filename = strit->second;

        vector< PolarSite* > poles = ParseGdmaFile(filename, state);

        map_seg_polarSites[segName] = poles;
        poles.clear();

    } /* Exit loop over files */

    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
    // Forward information on polar-site templates to topology //
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++ //

    // Containers to use for mapping
    // map<string, vector<PolarSite*> >     map_seg_polarSites;
    // map<string, vector<int> >            alloc_frag_mpoleIdx;
    // map<string, vector<string> >         alloc_frag_mpoleName;

    
    map<string, vector<PolarSite*> > ::iterator it1;
    for (it1 = map_seg_polarSites.begin();
            it1 != map_seg_polarSites.end();
            it1++) {

        vector<PolarSite*> ::iterator it2;
        vector<PolarSite*> poles = it1->second;

        for( it2 = poles.begin();
                it2 < poles.end();
                it2++) {

            PolarSite *pol = *it2;
            cout << endl << "ID " << pol->_id;
            cout << " RANK " << pol->_rank;

            vector<double> ::iterator dit;
            for (dit = pol->_Qs[1].begin();
                    dit < pol->_Qs[1].end();
                    dit++) {

                cout << " " << *dit;
            }
        }
    }

    // TO CONTINUE...
    // Loop over segments, loop over fragments in segments.
    // => For each segment, look up associated polar sites using map 1;
    //      => For each fragment, look up associated
    //         polar sites idcs. using map 2;
    //          => For each index in retrieved container, create new polar site,
    //             importing information from polar site container of segment
    //          => Translate positions to MD frame, copy local frame from
    //             fragment

    _map_seg_polarSites = map_seg_polarSites;
    _alloc_frag_mpoleIdx =  alloc_frag_mpoleIdx;
    _alloc_frag_mpoleName = alloc_frag_mpoleName;
    
}



vector<PolarSite*> EMultipole2::ParseGdmaFile(string filename, int state) {

    int poleCount = 1;
    string units = "";

    vector<PolarSite*> poles;
    PolarSite *thisPole = NULL;

    vector<double> Qs;

    std::string line;
    std::ifstream intt;
    intt.open(filename.c_str());

    if (intt.is_open() ) {
        while ( intt.good() ) {

            std::getline(intt, line);
            vector<string> split;
            Tokenizer toker(line, " ");
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


            // Units used
            if ( split[0] == "Units") {
                units = split[1];
                if (units != "bohr") {
                    throw std::runtime_error( "Unit " + units + " in file "
                                            + filename + " not supported.");
                }
            }

            // element,  position,  rank limit
            else if ( split.size() == 6 ) {

                Qs.clear();

                int id = poleCount++;  // <- starts from 1
                string name = split[0];

                double x = boost::lexical_cast<double>(split[1]);
                double y = boost::lexical_cast<double>(split[2]);
                double z = boost::lexical_cast<double>(split[3]);
                vec pos = vec(x,y,z);

                int rank = boost::lexical_cast<int>(split[5]);

                PolarSite *newPole = new PolarSite(id, name);
                newPole->setRank(rank);
                poles.push_back(newPole);
                thisPole = newPole;

            }

            // multipole line
            else {

                for (int i = 0; i < split.size(); i++) {

                    double qXYZ = boost::lexical_cast<double>(split[i]);
                    Qs.push_back(qXYZ);

                }
                thisPole->setQs(Qs, state);
            }

        } /* Exit loop over lines */
    }
    else { cout << endl << "ERROR: No file " << filename << endl; }

    return poles;
    
}













bool EMultipole2::EvaluateFrame(Topology *top) {

    // Rigidify if (a) not rigid yet (b) rigidification at all possible
    if (!top->isRigid()) {
        bool isRigid = top->Rigidify();
        if (!isRigid) { return 0; }
    }
    else { cout << endl << "... ... System is already rigidified."; }
    cout << endl;

    // +++++++++++++++++++++++++++++++++++++ //
    // Equip TOP with distributed multipoles //
    // +++++++++++++++++++++++++++++++++++++ //


    vector<Segment*> ::iterator sit;
    for (sit = top->Segments().begin();
         sit < top->Segments().end();
         ++sit) {

        Segment *seg = *sit;
        vector<PolarSite*> poleSites = _map_seg_polarSites.at(seg->getName());

        vector<Fragment*> ::iterator fit;
        for (fit = seg->Fragments().begin();
             fit < seg->Fragments().end();
             ++fit) {

            Fragment *frag = *fit;
            string idkey = frag->getName() + seg->getName()
                         + seg->getMolecule()->getName();
            vector<int> polesInFrag = _alloc_frag_mpoleIdx.at(idkey);
            vector<string> namesInFrag = _alloc_frag_mpoleName.at(idkey);

            for (int i = 0; i < polesInFrag.size(); i++) {

                string name = namesInFrag[i];
                int poleId = polesInFrag[i];
                PolarSite *templ = poleSites[poleId];

                PolarSite *newSite = top->AddPolarSite(name);
                newSite->ImportFrom(templ);



            }






        }



    }





    vector<SiteOperator*> siteOps;
    this->InitSlotData(top);

    _nextSite = top->Segments().begin();

    for (int id = 0; id < _nThreads; id++) {
        SiteOperator *newOp = new SiteOperator(id, top, this);
        siteOps.push_back(newOp);
    }

    for (int id = 0; id < _nThreads; id++) {
        siteOps[id]->Start();
    }

    for (int id = 0; id < _nThreads; id++) {
        siteOps[id]->WaitDone();
    }

    for (int id = 0; id < _nThreads; id++) {
        delete siteOps[id];
    }

    siteOps.clear();

    this->PostProcess(top);
    return 1;
}


void EMultipole2::EvalSite(Topology *top, Segment *seg, int slot) {

    this->LockCout();
    cout << "\r... ... Evaluate site " << seg->getId();
    this->UnlockCout();

}



// +++++++++++++++++ //
// Thread Management //
// +++++++++++++++++ //

Segment *EMultipole2::RequestNextSite(int opId, Topology *top) {

    _nextSiteMutex.Lock();

    Segment *workOnThis;

    if (_nextSite == top->Segments().end()) {
        workOnThis = NULL;
    }
    else {
        workOnThis = *_nextSite;
        _nextSite++;
    }

    _nextSiteMutex.Unlock();

    return workOnThis;
}


// +++++++++++++++++++++++++++++ //
// SiteOperator Member Functions //
// +++++++++++++++++++++++++++++ //

void EMultipole2::SiteOperator::Run(void) {

    while (true) {

        Segment *seg = _master->RequestNextSite(_id, _top);

        if (seg == NULL) { break; }
        else { this->_master->EvalSite(_top, seg, _id); }
    }
}

}}









#endif