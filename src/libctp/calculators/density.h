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
};


void Density::Initialize(Topology *top, Property *options) {

    string key      = "options.density";
    _axis           = options->get(key+".axis").as< vec >();
    _resolution     = options->get(key+".resolution").as< double >();
    _outfile        = options->get(key+".output").as< string >();

    // Normalize axis
    _axis           = _axis / sqrt(_axis*_axis);
}



bool Density::EvaluateFrame(Topology *top) {

    map< string, vector< double > > map_seg_zs;
    map< string, bool > set_seg;

    vector< Segment* > ::iterator sit;
    vector< Atom* > ::iterator ait;

    // Collect profile positions from atoms in system, order by segment name
    double MAX = -1e100;
    double MIN =  1e100;
    double RES = _resolution;

    for (sit = top->Segments().begin();
         sit < top->Segments().end();
         ++sit) {

        Segment *seg = *sit;
        if (!set_seg.count(seg->getName())) { set_seg[seg->getName()] = true; }

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
    
    int BIN = int( (MAX-MIN)/_resolution + 0.5 ) + 1;
    int SEG = set_seg.size();

    // Calculate density profile from z-list
    vector< vector< int > > seg_hist_Ns;
    map< string, bool > ::iterator setit;

    for (setit = set_seg.begin(); setit != set_seg.end(); ++setit) {

        // Retrieve appropriate z list
        string segName = setit->first;
        vector< double > seg_zs = map_seg_zs[segName];
        vector< vector< double > > hist_zs;
        hist_zs.resize(BIN);
        
        // Perform binning
        vector< double > ::iterator zit;
        for (zit = seg_zs.begin(); zit < seg_zs.end(); ++zit) {

            int bin = int( ((*zit)-MIN)/_resolution + 0.5 );
            hist_zs[bin].push_back((*zit));
        }

        // Reduce bins
        vector< int > hist_Ns;
        hist_Ns.resize(BIN);
        for (int bin = 0; bin < BIN; ++bin) {
            hist_Ns[bin] = hist_zs[bin].size();
        }

        seg_hist_Ns.push_back(hist_Ns);
    }   


    FILE *out;
    out = fopen(_outfile.c_str(), "w");

    fprintf(out, "# DENSITY PROFILE ALONG AXIS z = %4.7f %4.7f %4.7f \n",
            _axis.getX(), _axis.getY(), _axis.getZ());
    fprintf(out, "# MIN z %4.7f MAX z %4.7f \n", MIN, MAX);


    fprintf(out, "# z");
    for (setit = set_seg.begin(); setit != set_seg.end(); ++setit) {
        fprintf(out, " N(%5s,z) ", setit->first.c_str());
    }
    fprintf(out, " N(z) \n");

    for (int bin = 0; bin < BIN; ++bin) {

        double z = MIN + bin*_resolution;
        fprintf(out, "%4.7f", z);
        int N_z = 0;
        for (int s = 0; s < SEG; ++s) {
            fprintf(out, " %7d ", seg_hist_Ns[s][bin]);
            N_z += seg_hist_Ns[s][bin];
        }
        fprintf(out, " %7d \n", N_z);
    }
    fclose(out);

    return true;
}





}}

#endif

