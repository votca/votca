#include <votca/xtp/poissongrid.h>


namespace votca {
namespace xtp {
namespace POI {

PoissonGrid::PoissonGrid(Topology *top, vector<PolarSeg*> &fg, vector<PolarSeg*> &bg, 
    Logger *log) : _log(log) {
    
    
    
    
    _center = votca::tools::vec(0,0,0);
    for (vector<PolarSeg*>::iterator sit = fg.begin(); sit < fg.end();
        ++sit) _center += (*sit)->getPos();
    if (fg.size()> 0) _center /= fg.size();
    // Shift center to front lower left corner of central cell
    //_center = _center - 0.5*vec(dx0, dy0, dz0);
    
    LOG(logDEBUG,*_log) << flush;
    LOG(logDEBUG,*_log) << "Setup <PoissonGrid> @ " << _center << flush;
    
    
    // GRID PARAMETERS
    double dx0 = 0.5;
    double dx1 = 1.0;
    int threshold_dx0_dx1 = 4;
    double x0 = _center.getX();
    double x01 = threshold_dx0_dx1*dx0;
    double x1 = 6.;
    
    //double dy0 = 0.5;
    //double dy1 = 1.0;
    //int threshold_dy0_dy1 = 4;
    
    //double dz0 = 0.5;
    //double dz1 = 0.5;
    //int threshold_dz0_dz1 = 4;
    
    
    
    
    // x dimension
    vector< double > array_dx;
    double xi = x0;
    double dx = dx0;
    array_dx.push_back(dx);
    for (int i = 1; true; ++i) {
        
        if (xi <= x01) {
            dx = dx0;            
        }
        else if (xi <= x1) {
            dx =  (xi - x01) / (x1 - x01) * dx1
                + (x1 - xi)  / (x1 - x01) * dx0;
        }
        else {
            dx = dx1;
        }
        
        xi += dx;
        array_dx.push_back(dx);

        cout << endl << dx << " => " << xi << flush;
        if (xi > x1) break;
    }
    
    
    
    
    
    


}

}}}
