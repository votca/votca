#ifndef DENSITY_H
#define DENSITY_H


namespace votca { namespace ctp {


#include <votca/ctp/qmcalculator.h>

class Density : public QMCalculator
{
public:

    string      Identify() { return "Density"; }
    void        Initialize(Topology *top, Property *options);
    bool        EvaluateFrame(Topology *top);

private:

    vec         _axis;
    double      _resolution;
    string      _outfile;
    string      _outfile_EA_IP;
    string      _density_type;
    bool        _auto_bin;
    double      _min;
    double      _max;
};


void Density::Initialize(Topology *top, Property *options) {

    string key      = "options.density";
    _axis           = options->get(key+".axis").as< vec >();
    _resolution     = options->get(key+".resolution").as< double >();
    _outfile        = options->get(key+".output").as< string >();
    _outfile_EA_IP  = options->get(key+".output_e").as< string >();
    _density_type   = options->get(key+".density_type").as< string >();

    int autobin     = options->get(key+".auto_bin").as< int >();
    _auto_bin       = (autobin == 1) ? true : false;

    if (!_auto_bin) {
        _min            = options->get(key+".min").as< double >();
        _max            = options->get(key+".max").as< double >();
    }
    
    // Normalize axis
    _axis           = _axis / sqrt(_axis*_axis);
}



bool Density::EvaluateFrame(Topology *top) {

    map< string, vector< double > > map_seg_zs; // for atomistic number density
    map< string, vector< double > > map_com_zs; // for segment number density
    map< string, vector< Segment* > > map_com_seg;
    map< string, bool > set_seg;

    vector< Segment* > ::iterator sit;
    vector< Atom* > ::iterator ait;

    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
    // Collect profile positions from atoms in system, order by segment name //
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //

    double MAX = -1e100;
    double MIN =  1e100;
    double RES = _resolution;

    for (sit = top->Segments().begin();
         sit < top->Segments().end();
         ++sit) {

        Segment *seg = *sit;
        if (!set_seg.count(seg->getName())) { set_seg[seg->getName()] = true; }

        double z_com = seg->getPos() * _axis;

        map_com_zs[seg->getName()].push_back(z_com);
        map_com_seg[seg->getName()].push_back(seg);

        // Include this segment in density profile?
        // if (!wildcmp(seg->getName().c_str(),_seg_pattern.c_str())) {
        //     continue;
        // }

        for (ait = seg->Atoms().begin();
             ait < seg->Atoms().end();
             ++ait) {

            double z = (*ait)->getPos() * _axis;

            MAX = (z > MAX) ? z : MAX;
            MIN = (z < MIN) ? z : MIN;

            map_seg_zs[seg->getName()].push_back(z);
        }
    }

    if (!_auto_bin) {
        MAX = _max;
        MIN = _min;
    }

    int BIN = int( (MAX-MIN)/_resolution + 0.5 ) + 1;
    int SEG = set_seg.size();

    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
    // Calculate density profile from z-list, segment-name based //
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //

    vector< vector< int > > seg_hist_Ns;
    vector< vector< int > > com_hist_Ns;
    vector< vector< double > > com_profile_ea;
    vector< vector< double > > com_profile_ec;
    map< string, bool > ::iterator setit;

    for (setit = set_seg.begin(); setit != set_seg.end(); ++setit) {

        // Retrieve appropriate z list
        string segName = setit->first;
        vector< double > seg_zs = map_seg_zs[segName];
        vector< double > com_zs = map_com_zs[segName];
        vector< Segment* > com_seg = map_com_seg[segName];
        vector< vector< double > > hist_zs;
        vector< vector< double > > hist_zs_com;
        vector< vector< double > > binned_ea; // anion site energies  (<> EA)
        vector< vector< double > > binned_ec; // cation site energies (<> IP)
        hist_zs.resize(BIN);
        hist_zs_com.resize(BIN);
        binned_ea.resize(BIN);
        binned_ec.resize(BIN);
        
        // Perform binning
        vector< double > ::iterator zit;
        vector< Segment* > ::iterator sit;
        for (zit = seg_zs.begin(); zit < seg_zs.end(); ++zit) {

            int bin = int( ((*zit)-MIN)/_resolution + 0.5 );
            hist_zs[bin].push_back((*zit));
        }

        for (zit = com_zs.begin(), sit = com_seg.begin();
             zit < com_zs.end();
             ++zit, ++sit) {

            int bin = int( ((*zit)-MIN)/_resolution + 0.5 );
            hist_zs_com[bin].push_back((*zit));
            binned_ea[bin].push_back((*sit)->getSiteEnergy(-1));
            binned_ec[bin].push_back((*sit)->getSiteEnergy(+1));
        }

        // Reduce bins
        vector< int > hist_Ns;
        vector< int > hist_Ns_com;
        vector< double > binned_avg_ea;
        vector< double > binned_avg_ec;
        hist_Ns.resize(BIN);
        hist_Ns_com.resize(BIN);
        binned_avg_ea.resize(BIN);
        binned_avg_ec.resize(BIN);
        // Atomistic number density
        for (int bin = 0; bin < BIN; ++bin) {
            hist_Ns[bin] = hist_zs[bin].size();
        }
        // Segment-based number density
        for (int bin = 0; bin < BIN; ++bin) {
            hist_Ns_com[bin] = hist_zs_com[bin].size();
        }
        // IP and EA
        for (int bin = 0; bin < BIN; ++bin) {
            double avg_ea = 0.0;
            double avg_ec = 0.0;

            for (int entry = 0; entry < binned_ea[bin].size(); ++entry) {
                avg_ea += binned_ea[bin][entry];
                avg_ec += binned_ec[bin][entry];
            }
            avg_ea /= binned_ea[bin].size();
            avg_ec /= binned_ec[bin].size();

            binned_avg_ea[bin] = avg_ea;
            binned_avg_ec[bin] = avg_ec;
        }
        
        seg_hist_Ns.push_back(hist_Ns);
        com_hist_Ns.push_back(hist_Ns_com);
        com_profile_ea.push_back(binned_avg_ea);
        com_profile_ec.push_back(binned_avg_ec);
    }


    // +++++++++++++++++++++++ //
    // Output number densities //
    // +++++++++++++++++++++++ //

    FILE *out;
    out = fopen(_outfile.c_str(), "w");

    fprintf(out, "# DENSITY PROFILE ALONG AXIS z = %4.7f %4.7f %4.7f \n",
            _axis.getX(), _axis.getY(), _axis.getZ());
    fprintf(out, "# MIN z %4.7f MAX z %4.7f \n", MIN, MAX);


    fprintf(out, "# z");
    for (setit = set_seg.begin(); setit != set_seg.end(); ++setit) {
        fprintf(out, " N(%1s,z) ", setit->first.c_str());
    }
    fprintf(out, " N(z) \n");

    for (int bin = 0; bin < BIN; ++bin) {

        double z = MIN + bin*_resolution;
        fprintf(out, "%4.7f", z);
        int N_z = 0;
        if (_density_type == "atoms") {
            for (int s = 0; s < SEG; ++s) {
                fprintf(out, " %7d ", seg_hist_Ns[s][bin]);
                N_z += seg_hist_Ns[s][bin];
            }
        }
        else if (_density_type == "segments") {
            for (int s = 0; s < SEG; ++s) {
                fprintf(out, " %7d ", com_hist_Ns[s][bin]);
                N_z += com_hist_Ns[s][bin];
            }
        }
        else {
            cout << endl << "ERROR: Invalid density type " << _density_type
                 << ". Should be either atoms or segments" << flush;
            break;
        }

        fprintf(out, " %7d \n", N_z);
    }

    fclose(out);


    // ++++++++++++++++++++++++++++++++ //
    // Output spatially-resolved EA, IP //
    // ++++++++++++++++++++++++++++++++ //

    out = fopen(_outfile_EA_IP.c_str(), "w");

    fprintf(out, "# ENERGY LANDSCAPE ALONG AXIS z = %4.7f %4.7f %4.7f \n",
            _axis.getX(), _axis.getY(), _axis.getZ());
    fprintf(out, "# MIN z %4.7f MAX z %4.7f \n", MIN, MAX);

    fprintf(out, "# z");
    for (setit = set_seg.begin(); setit != set_seg.end(); ++setit) {
        fprintf(out, " EA(%1s,z) ", setit->first.c_str());
        fprintf(out, " IP(%1s,z) ", setit->first.c_str());
    }
    fprintf(out, " \n");

    for (int bin = 0; bin < BIN; ++bin) {

        double z = MIN + bin*_resolution;
        fprintf(out, "%4.7f", z);

        for (int s = 0; s < SEG; ++s) {
            fprintf(out, " %4.7f ", com_profile_ea[s][bin]);
            fprintf(out, " %4.7f ", com_profile_ec[s][bin]);
        }

        fprintf(out, " \n");
    }
    fclose(out);

    return true;
}





}}

#endif

