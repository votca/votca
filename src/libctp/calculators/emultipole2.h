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
    void     DistributeMpoles(Topology *top);

    bool     EvaluateFrame(Topology *top);
    void     EvalSite(Topology *top, Segment *seg, int slot);
    void     Induce() { return; }
    void     CalcIntEnergy() { return; }
    
    void     InitSlotData(Topology *top) { ; }
    void     PostProcess(Topology *top) { ; }
    Segment *RequestNextSite(int opId, Topology *top);
    void     LockCout() { _coutMutex.Lock(); }
    void     UnlockCout() { _coutMutex.Unlock(); }


    // +++++++++++++++++++++++++++ //
    // Multipole Interaction Class //
    // +++++++++++++++++++++++++++ //

    class Interactor
    {
    public:

        Interactor(Topology *top) : _top(top) {};
        Interactor() {};
       ~Interactor() {};

        inline double Energy(PolarSite &pol1, PolarSite &pol2);

    private:

        vec    e12;
        double R;       //  |
        double R2;      //  |
        double R3;      //  |-> NOTE: reciprocal, i.e. e.g. R3 = 1/(R*R*R)
        double R4;      //  |
        double R5;      //  |

        double rax, ray, raz;
        double rbx, rby, rbz;
        double cxx, cxy, cxz;
        double cyx, cyy, cyz;
        double czx, czy, czz;

        inline double T00_00() { return R; }

        inline double T1x_00() { return R2 * rax; }
        inline double T1y_00() { return R2 * ray; }
        inline double T1z_00() { return R2 * raz; }
        inline double T00_1x() { return R2 * rbx; }
        inline double T00_1y() { return R2 * rby; }
        inline double T00_1z() { return R2 * rbz; }

        inline double T20_00()  { return R3 * 0.5 * (3 * raz*raz - 1); }
        inline double T21c_00() { return R3 * sqrt(3) * rax * raz; }
        inline double T21s_00() { return R3 * sqrt(3) * ray * raz; }
        inline double T22c_00() { return R3 * 0.5 * sqrt(3) * (rax*rax - ray*ray); }
        inline double T22s_00() { return R3 * sqrt(3) * rax*ray; }
        inline double T00_20()  { return R3 * 0.5 * (3 * rbz*rbz - 1); }
        inline double T00_21c() { return R3 * sqrt(3) * rbx * rbz; }
        inline double T00_21s() { return R3 * sqrt(3) * rby * rbz; }
        inline double T00_22c() { return R3 * 0.5 * sqrt(3) * (rbx*rbx - rby*rby); }
        inline double T00_22s() { return R3 * sqrt(3) * rbx*rby; }

        inline double T1x_1x() { return R3 * (3 * rax*rbx + cxx); }
        inline double T1x_1y() { return R3 * (3 * rax*rby + cxy); }
        inline double T1x_1z() { return R3 * (3 * rax*rbz + cxz); }
        inline double T1y_1x() { return R3 * (3 * ray*rbx + cyx); }
        inline double T1y_1y() { return R3 * (3 * ray*rby + cyy); }
        inline double T1y_1z() { return R3 * (3 * ray*rbz * cyz); }
        inline double T1z_1x() { return R3 * (3 * raz*rbx * czx); }
        inline double T1z_1y() { return R3 * (3 * raz*rby * czy); }
        inline double T1z_1z() { return R3 * (3 * raz*rbz * czz); }

        inline double T20_1x()  { return R4 * 0.5 * (15*raz*raz*rbx + 6*raz*czx - 3*rbx); }
        inline double T20_1y()  { return R4 * 0.5 * (15*raz*raz*rby + 6*raz*czy - 3*rby); }
        inline double T20_1z()  { return R4 * 0.5 * (15*raz*raz*rbz + 6*raz*czz - 3*rbz); }
        inline double T21c_1x() { return R4 * sqrt(3) * (rax*czx + cxx*raz + 5*rax*raz*rbx); }
        inline double T21c_1y() { return R4 * sqrt(3) * (rax*czy + cxy*raz + 5*rax*raz*rby); }
        inline double T21c_1z() { return R4 * sqrt(3) * (rax*czz + cxz*raz + 5*rax*raz*rbz); }
        inline double T21s_1x() { return R4 * sqrt(3) * (ray*czx + cyx*raz + 5*ray*raz*rbx); }
        inline double T21s_1y() { return R4 * sqrt(3) * (ray*czy + cyy*raz + 5*ray*raz*rby); }
        inline double T21s_1z() { return R4 * sqrt(3) * (ray*czz + cyz*raz + 5*ray*raz*rbz); }
        inline double T22c_1x() { return R4 * 0.5 * sqrt(3) * ( 5*(rax*rax-ray*ray)*rbx + 2*rax*cxx - 2*ray*cyx); }
        inline double T22c_1y() { return R4 * 0.5 * sqrt(3) * ( 5*(rax*rax-ray*ray)*rby + 2*rax*cxy - 2*ray*cyy); }
        inline double T22c_1z() { return R4 * 0.5 * sqrt(3) * ( 5*(rax*rax-ray*ray)*rbz + 2*rax*cxz - 2*ray*cyz); }
        inline double T22s_1x() { return R4 * sqrt(3) * ( 5*rax*ray*rbx + rax*cyx + ray*cxx ); }
        inline double T22s_1y() { return R4 * sqrt(3) * ( 5*rax*ray*rby + rax*cyy * ray*cxy ); }
        inline double T22s_1z() { return R4 * sqrt(3) * ( 5*rax*ray*rbz + rax*cyz * ray*cxz ); }
        
        inline double T1x_20()  { return R4 * 0.5 * (15*rbz*rbz*rax + 6*rbz*cxz - 3*rax); }
        inline double T1y_20()  { return R4 * 0.5 * (15*rbz*rbz*ray + 6*rbz*cyz - 3*ray); }
        inline double T1z_20()  { return R4 * 0.5 * (15*rbz*rbz*raz + 6*rbz*czz - 3*raz); }
        inline double T1x_21c() { return R4 * sqrt(3) * (rbx*cxz + cxx*rbz + 5*rbx*rbz*rax); }
        inline double T1y_21c() { return R4 * sqrt(3) * (rbx*cyz + cyx*rbz + 5*rbx*rbz*ray); }
        inline double T1z_21c() { return R4 * sqrt(3) * (rbx*czz + czx*rbz + 5*rbx*rbz*raz); }
        inline double T1x_21s() { return R4 * sqrt(3) * (rby*cxz + cxy*rbz + 5*rby*rbz*rax); }
        inline double T1y_21s() { return R4 * sqrt(3) * (rby*cyz + cyy*rbz + 5*rby*rbz*ray); }
        inline double T1z_21s() { return R4 * sqrt(3) * (rby*czz + czy*rbz + 5*rby*rbz*raz); }
        inline double T1x_22c() { return R4 * 0.5 * sqrt(3) * ( 5*(rbx*rbx-rby*rby)*rax + 2*rbx*cxx - 2*rby*cxy); }
        inline double T1y_22c() { return R4 * 0.5 * sqrt(3) * ( 5*(rbx*rbx-rby*rby)*ray + 2*rbx*cyx - 2*rby*cyy); }
        inline double T1z_22c() { return R4 * 0.5 * sqrt(3) * ( 5*(rbx*rbx-rby*rby)*raz + 2*rbx*czx - 2*rby*czy); }
        inline double T1x_22s() { return R4 * sqrt(3) * ( 5*rbx*rby*rax + rbx*cxy + rby*cxx ); }
        inline double T1y_22s() { return R4 * sqrt(3) * ( 5*rbx*rby*ray + rbx*cyy * rby*cyx ); }
        inline double T1z_22s() { return R4 * sqrt(3) * ( 5*rbx*rby*raz + rbx*czy * rby*czx ); }

        inline double T20_20()  { return R5 * 0.75 * (35*raz*raz*rbz*rbz - 5*raz*raz - 5*rbz*rbz + 20*raz*rbz*czz + 2*czz*czz + 1); }
        inline double T20_21c() { return R5 * 0.5 * sqrt(3) * (35*raz*raz*rbx*rbz - 5*rbx*rbz + 10*raz*rbx*czz + 10*raz*rbz*czx + 2*czx*czz); }
        inline double T20_21s() { return R5 * 0.5 * sqrt(3) * (35*raz*raz*rby*rbz - 5*rby*rbz + 10*raz*rby*czz + 10*raz*rbz*czy + 2*czy*czz); }
        inline double T20_22c() { return R5 * 0.25 * sqrt(3) * (35*raz*raz*rbx*rbx - 35*raz*raz*rby*rby - 5*rbx*rbx - 5*rby*rby + 20*raz*rbx*czx - 20*raz*rby*czy + 2*czx*czx - 2*czy*czy); }
        inline double T20_22s() { return R5 * 0.5 * sqrt(3) * (35*raz*raz*rbx*rby - 5*rbx*rby + 10*raz*rbx*czy + 10*raz*rby*czx + 2*czx*czy); }
        inline double T21c_21c() { return R5 * (35*rax*raz*rbx*rbz + 5*rax*rbx*czz + 5*rax*rbz*czx + 5*raz*rbx*cxz + 5*raz*rbz*cxx + cxx*czz + cxz*czx); }
        inline double T21c_21s() { return R5 * (35*rax*raz*rby*rbx + 5*rax*rby*czz + 5*rax*rbz*czy + 5*raz*rby*cxz + 5*raz*rbz*cxy + cxy*czz + cxz*czy); }
        inline double T21c_22c() { return R5 * 0.5 * (35*rax*raz*rbx*rbx - 35*rax*raz*rby*rby + 10*rax*rbx*czx - 10*rax*rby*czy + 10*raz*rbx*cxx - 10*raz*rby*cxy + 2*cxx*czx - 2*cxy*czy); }
        inline double T21c_22s() { return R5 * (35*rax*raz*rbx*rby * 5*rax*rbx*czy + 5*rax*rby*czx + 5*raz*rbx*cxy + 5*raz*rby*cxx + cxx*czy + cxy*czx); }
        inline double T21s_21s() { return R5 * (35*ray*raz*rby*rbz + 5*ray*rby*czz + 5*ray*rbz*czy + 5*raz*rby*cyz + 5*raz*rbz*cyy + cyy*czz + cyz*czy); }
        inline double T21s_22c() { return R5 * 0.5 * (35*ray*raz*rbx*rbx - 35*ray*raz*rby*rby + 10*ray*rbx*czx - 10*ray*rby*czy + 10*raz*rbx*cyx - 10*raz*rby*cyy + 2*cyx*czx - 2*cyy*czy); }
        inline double T21s_22s() { return R5 * (35*ray*raz*rbx*rby + 5*ray*rbx*czy + 5*ray*rby*czx + 5*raz*rbx*cyy + 5*raz*rby*cyx + cyx*czy + cyy*czx); }
        inline double T22c_22c() { return R5 * 0.25 * (35*rax*rax*rbx*rbx - 35*rax*rax*rby*rby - 35*ray*ray*rbx*rbx + 35*ray*ray*rby*rby + 20*rax*rbx*cxx - 20*rax*rby*cxy - 20*ray*rbx*cyx + 20*ray*rby*cyy + 2*cxx*cxx - 2*cxy*cxy - 2*cyx*cyx + 2*cyy*cyy); }
        inline double T22c_22s() { return R5 * 0.5 * (35*rax*rax*rbx*rby - 35*ray*ray*rbx*rby + 10*rax*rbx*cxy + 10*rax*rby*cxx - 10*ray*rbx*cyy - 10*ray*rby*cyx + 2*cxx*cxy - 2*cyx*cyy); }
        inline double T22s_22s() { return R5 * (35*rax*ray*rbx*rby + 5*rax*rbx*cyy + 5*rax*rby*cyx + 5*ray*rbx*cxy + 5*ray*rby*cxx + cxx*cyy + cxy*cyx); }

        inline double T21c_20() { return R5 * 0.5 * sqrt(3) * (35*rbz*rbz*rax*raz - 5*rax*raz + 10*rbz*rax*czz + 10*rbz*raz*czx + 2*czx*czz); }
        inline double T21s_20() { return R5 * 0.5 * sqrt(3) * (35*rbz*rbz*ray*raz - 5*ray*raz + 10*rbz*ray*czz + 10*rbz*raz*cyz + 2*cyz*czz); }
        inline double T22c_20() { return R5 * 0.25 * sqrt(3) * (35*rbz*rbz*rax*rax - 35*rbz*rbz*ray*ray - 5*rax*rax - 5*ray*ray + 20*rbz*rax*cxz - 20*rbz*ray*cyz + 2*cxz*cxz - 2*cyz*cyz); }
        inline double T22s_20() { return R5 * 0.5 * sqrt(3) * (35*rbz*rbz*rax*ray - 5*rax*ray + 10*rbz*rax*cyz + 10*rbz*ray*cxz + 2*cxz*cyz); }
        inline double T21s_21c() { return R5 * (35*rbx*rbz*ray*rax + 5*rbx*ray*czz + 5*rbx*raz*cyz + 5*rbz*ray*czx + 5*rbz*raz*cyx + cyx*czz + czx*cyz); }
        inline double T22c_21c() { return R5 * 0.5 * (35*rbx*rbz*rax*rax - 35*rbx*rbz*ray*ray + 10*rbx*rax*cxz - 10*rbx*ray*cyz + 10*rbz*rax*cxx - 10*rbz*ray*cyx + 2*cxx*cxz - 2*cyx*cyz); }
        inline double T22s_21c() { return R5 * (35*rbx*rbz*rax*ray * 5*rbx*rax*cyz + 5*rbx*ray*cxz + 5*rbz*rax*cyx + 5*rbz*ray*cxx + cxx*cyz + cyx*cxz); }
        inline double T22c_21s() { return R5 * 0.5 * (35*rby*rbz*rax*rax - 35*rby*rbz*ray*ray + 10*rby*rax*cxz - 10*rby*ray*cyz + 10*rbz*rax*cxy - 10*rbz*ray*cyy + 2*cxy*cxz - 2*cyy*cyz); }
        inline double T22s_21s() { return R5 * (35*rby*rbz*rax*ray + 5*rby*rax*cyz + 5*rby*ray*cxz + 5*rbz*rax*cyy + 5*rbz*ray*cxy + cxy*cyz + cyy*cxz); }
        inline double T22s_22c() { return R5 * 0.5 * (35*rbx*rbx*rax*ray - 35*rby*rby*rax*ray + 10*rbx*rax*cyx + 10*rbx*ray*cxx - 10*rby*rax*cyy - 10*rby*ray*cxy + 2*cxx*cyx - 2*cxy*cyy); }
        
        Topology *_top;

    };


    // ++++++++++++++++++++++++++++++++++++++ //
    // Site workers (i.e. individual threads) //
    // ++++++++++++++++++++++++++++++++++++++ //

    class SiteOpMultipole : public Thread
    {
    public:

        SiteOpMultipole(int id, Topology *top,
                     EMultipole2 *master)
                   : _id(id), _top(top), _seg(NULL),
                     _master(master)      { _actor = Interactor(top); };

       ~SiteOpMultipole();

        int  getId() { return _id; }
        void setId(int id) { _id = id; }

        void InitSlotData(Topology *top);
        void Run(void);

        vector<PolarSite*> &PolarSites() { return _polarSites; }


    protected:

        int                      _id;
        Topology                *_top;
        Segment                 *_seg;
        EMultipole2             *_master;

        vector<PolarSite*>       _polarSites;
        Interactor               _actor;
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

    // interaction parameters
    bool            _useCutoff;
    double          _cutoff;
    bool            _useExp;
    double          _aDamp;
    bool            _useScaling;
    vector<double>  _scale1;

    // convergence parameters
    float           _omegSOR;
    double          _epsTol;
    int             _maxIter;




};



void EMultipole2::Initialize(Topology *top, Property *opt) {

    cout << endl << "... ... Initialize with " << _nThreads << " threads.";

    if (!top->isEStatified()) { this->EStatify(top, opt); }


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

            // Merge with polar sites for neutral state
            vector< PolarSite* > polesNeutral = map_seg_polarSites[segName];

            assert(polesAnion.size() == polesNeutral.size());
            for (int i = 0; i < polesNeutral.size(); i++) {

                polesNeutral[i]->setQs( polesAnion[i]->getQs(state), state );
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

            // Merge with polar sites for neutral state
            vector< PolarSite* > polesNeutral = map_seg_polarSites[segName];

            assert(polesCation.size() == polesNeutral.size());
            for (int i = 0; i < polesNeutral.size(); i++) {

                polesNeutral[i]->setQs( polesCation[i]->getQs(state), state );
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

                double BOHR2NM = 0.0529189379;
                double x = BOHR2NM * boost::lexical_cast<double>(split[1]);
                double y = BOHR2NM * boost::lexical_cast<double>(split[2]);
                double z = BOHR2NM * boost::lexical_cast<double>(split[3]);
                vec pos = vec(x,y,z);

                int rank = boost::lexical_cast<int>(split[5]);

                PolarSite *newPole = new PolarSite(id, name);
                newPole->setRank(rank);
                newPole->setPos(pos);
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


void EMultipole2::DistributeMpoles(Topology *top) {

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

                PolarSite *templ = poleSites[poleId-1];
                PolarSite *newSite = top->AddPolarSite(name);
                newSite->ImportFrom(templ);
                seg->AddPolarSite(newSite);
                frag->AddPolarSite(newSite);

                // Shift + rotate
                newSite->Translate(frag->getTransQM2MD());
                newSite->Rotate(frag->getRotQM2MD(), frag->getCoMD());
            }
        }
    }
}










bool EMultipole2::EvaluateFrame(Topology *top) {

    // Rigidify if (a) not rigid yet (b) rigidification at all possible
    if (!top->isRigid()) {
        bool isRigid = top->Rigidify();
        if (!isRigid) { return 0; }
    }
    else { cout << endl << "... ... System is already rigidified."; }

    // Forward multipoles to topology
    this->DistributeMpoles(top);

    cout << endl << "... ... Created " << top->PolarSites().size()
                 << " multipole sites." << endl;

    /*
    string mpNAME = "mp_system.pdb";
    FILE *mpPDB = NULL;
    mpPDB = fopen(mpNAME.c_str(), "w");
    vector<Segment*>::iterator sit;
    for (sit = top->Segments().begin(); sit < top->Segments().end(); ++sit) {
        (*sit)->WritePDB(mpPDB, "Multipoles", "MD");
    }
    fclose(mpPDB);
    */

    /*
    vector<PolarSite*>::iterator pol;
    for (pol = top->PolarSites().begin();
            pol < top->PolarSites().end();
            ++pol) {

        (*pol)->PrintInfo(cout);
    }
    cout << "SIZEOF " << sizeof(PolarSite) << endl;
    */


    vector<Segment*> ::iterator sit1;
    vector<Segment*> ::iterator sit2;
    vector<PolarSite*> ::iterator pol1;
    vector<PolarSite*> ::iterator pol2;

    int count = 0;
    Interactor actor = Interactor(top);

    for (sit1 = top->Segments().begin();
         sit1 < top->Segments().end();
         ++sit1) {
    for (sit2 = sit1 + 1;
         sit2 < top->Segments().end();
         ++sit2) {

        Segment *seg1 = *sit1;
        Segment *seg2 = *sit2;


        if ( abs( top->PbShortestConnect(seg1->getPos(),
                                         seg2->getPos()) ) > _cutoff) {
            continue;
        }

        ++count;
        /*
        cout << count << " " << seg1->getId() << " "<< seg2->getId() << " "
             << abs( top->PbShortestConnect(seg1->getPos(),
                                         seg2->getPos()) )
             << endl; */



        
        for (pol1 = seg1->PolarSites().begin();
             pol1 < seg1->PolarSites().end();
             ++pol1) {
        for (pol2 = seg2->PolarSites().begin();
             pol2 < seg2->PolarSites().end();
             ++pol2) {

             double E = actor.Energy(*(*pol1), *(*pol2));




        }} 
    }}





    // Create + start threads
    vector<SiteOpMultipole*> siteOps;
    

    _nextSite = top->Segments().begin();

    for (int id = 0; id < _nThreads; id++) {
        SiteOpMultipole *newOp = new SiteOpMultipole(id, top, this);
        siteOps.push_back(newOp);
    }

    for (int id = 0; id < _nThreads; id++) {
        siteOps[id]->InitSlotData(top);
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

void EMultipole2::SiteOpMultipole::Run(void) {

    while (true) {

        Segment *seg = _master->RequestNextSite(_id, _top);

        if (seg == NULL) { break; }
        else { this->_master->EvalSite(_top, seg, _id); }
    }
}


void EMultipole2::SiteOpMultipole::InitSlotData(Topology *top) {

    _polarSites.reserve(top->PolarSites().size());

    vector<PolarSite*> ::iterator pol;
    for (pol = top->PolarSites().begin();
         pol < top->PolarSites().end();
         ++pol) {
        
        PolarSite *threadPole = new PolarSite();
        threadPole->ImportFrom(*pol, "full");
        _polarSites.push_back(threadPole);
    }
}


EMultipole2::SiteOpMultipole::~SiteOpMultipole() {

    vector<PolarSite*> ::iterator pol;
    for (pol = _polarSites.begin(); pol < _polarSites.end(); ++pol) {
        delete *pol;
    }

    _polarSites.clear();

}

// +++++++++++++++++++++++++++ //
// Interactor Member Functions //
// +++++++++++++++++++++++++++ /&/

inline double EMultipole2::Interactor::Energy(PolarSite &pol1,PolarSite &pol2) {

    double E = 0.0;

    e12  = _top->PbShortestConnect(pol1.getPos(), pol2.getPos());
    R    = 1/abs(e12);
    R2   = R*R;
    R3   = R2*R;
    R4   = R2*R2;
    R5   = R2*R3;
    e12 /= R;
    // NOTE >>> e12 points from pole 1 to pole 2 <<< NOTE

    

    if (pol1._rank > 0 || pol2._rank > 0) {
        rax =   pol1._locX * e12;
        ray =   pol1._locY * e12;
        raz =   pol1._locZ * e12;
        rbx = - pol2._locX * e12;
        rby = - pol2._locY * e12;
        rbz = - pol2._locZ * e12;

        cxx = pol1._locX * pol2._locX;
        cxy = pol1._locX * pol2._locY;
        cxz = pol1._locX * pol2._locZ;
        cyx = pol1._locY * pol2._locX;
        cyy = pol1._locY * pol2._locY;
        cyz = pol1._locY * pol2._locZ;
        czx = pol1._locZ * pol2._locX;
        czy = pol1._locZ * pol2._locY;
        czz = pol1._locZ * pol2._locZ;
    }



    

    
    

    return E;
}


}}









#endif