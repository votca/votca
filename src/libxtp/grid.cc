#include <votca/xtp/grid.h>
#include <math.h>       /* ceil */

using namespace std;
using namespace votca::tools;


namespace votca { namespace xtp {
    namespace ub = boost::numeric::ublas;
    
Grid::Grid(const Grid &obj)
    :_cutoff(obj._cutoff),_gridspacing(obj._gridspacing),_cutoff_inside(obj._cutoff_inside),_shift_cutoff(obj._shift_cutoff),
    _shift_cutoff_inside(obj._shift_cutoff_inside),_useVdWcutoff(obj._useVdWcutoff),_useVdWcutoff_inside(obj._useVdWcutoff_inside),
    _cubegrid(obj._cubegrid),_padding(obj._padding),_createpolarsites(obj._createpolarsites),
    _xsteps(obj._xsteps),_ysteps(obj._ysteps),_zsteps(obj._zsteps) {
    _lowerbound=obj._lowerbound;
    _gridpoints=obj._gridpoints;
    std::vector<APolarSite*>::const_iterator pit;
    for(pit=obj._all_gridsites.begin();pit<obj._all_gridsites.end();++pit){
       APolarSite *apolarsite= new APolarSite(*pit,false);
       if(!apolarsite->getIsVirtual()) _gridsites.push_back(apolarsite);
       _all_gridsites.push_back(apolarsite);   
    }     
     _sites_seg = new PolarSeg(0, _gridsites);
     _atomlist=obj._atomlist;
    };
        
        
Grid::~Grid() {
        std::vector<APolarSite*>::iterator pit;
        for(pit=_all_gridsites.begin();pit!=_all_gridsites.end();++pit){
             delete *pit;
        }
        _all_gridsites.clear();
        if (_sites_seg != NULL) delete _sites_seg;
    }

Grid &Grid::operator=(const Grid & obj){
    _cutoff=obj._cutoff;
    _gridspacing=obj._gridspacing;
    _cutoff_inside=obj._cutoff_inside;
    _shift_cutoff=obj._shift_cutoff;
    _shift_cutoff_inside=obj._shift_cutoff_inside;
    _useVdWcutoff=obj._useVdWcutoff;
    _useVdWcutoff_inside=obj._useVdWcutoff_inside;
    _cubegrid=obj._cubegrid;
    _padding=obj._padding;
    _createpolarsites=obj._createpolarsites;
    _gridpoints=obj._gridpoints;
    _lowerbound=obj._lowerbound;
    _xsteps=obj._xsteps;
    _ysteps=obj._ysteps;
    _zsteps=obj._zsteps;
    std::vector<APolarSite*>::const_iterator pit;
    for(pit=obj._all_gridsites.begin();pit<obj._all_gridsites.end();++pit){
       APolarSite *apolarsite= new APolarSite(*pit,false);
       if(!apolarsite->getIsVirtual()) _gridsites.push_back(apolarsite);
       _all_gridsites.push_back(apolarsite);   
    }     
     _sites_seg = new PolarSeg(0, _gridsites);
     _atomlist=obj._atomlist;
     return *this;
}
    
void Grid::printGridtoxyzfile(const char* _filename){
        //unit is Angstrom in xyz file 
        double NmtoA=10.0; 
        ofstream points;
        points.open(_filename, ofstream::out);
        points << _gridpoints.size() << endl;
        points << endl;
        for ( unsigned i = 0 ; i < _gridpoints.size(); i++){
            points << "X " << _gridpoints[i](0)*NmtoA << " " << _gridpoints[i](1)*NmtoA << " " << _gridpoints[i](2)*NmtoA << endl;

        }
        points.close();
    }    


void Grid::readgridfromCubeFile(string filename, bool ignore_zeros){
        Elements _elements;
        _cubegrid=true;   
        double Bohr2Nm=1.0/18.897259886; 
        double Bohr2A=1.0/1.8897259886; 
        if(_gridpoints.size()>0) throw std::runtime_error("Grid object already has points.");
        ifstream in1;
        string s;
        in1.open(filename.c_str(), ios::in);
        getline(in1, s);
        getline(in1, s);
        int natoms;
        double xstart,ystart,zstart;
        double xincr,yincr,zincr;

        double tempdouble;
        string name="H";
        in1 >> natoms;
        in1 >> xstart;
        in1 >> ystart;
        in1 >> zstart;
        in1 >> _xsteps;
        in1 >> xincr;
        in1 >> tempdouble;
        in1 >> tempdouble;
        in1 >> _ysteps;
        in1 >> tempdouble;
        in1 >> yincr;
        in1 >> tempdouble;
        in1 >> _zsteps;
        in1 >> tempdouble;
        in1 >> tempdouble;
        in1 >> zincr;          
        
        if(xincr==yincr && yincr==zincr && zincr==xincr) _gridspacing=xincr*Bohr2A;
        else throw std::runtime_error("Gridspacing in x,y,z is different, currently not implemented, loading aborted");
        _lowerbound=vec(xstart*Bohr2A,ystart*Bohr2A,zstart*Bohr2A);
        
        
        _atomlist= new vector< QMAtom* >;
        for (int iatom =0; iatom < std::abs(natoms); iatom++) {
                 // get center coordinates in Bohr
                 double x ;
                 double y ;
                 double z ;
                 int atnum ;
                 double crg ;

                 // get from first cube
                 in1 >> atnum;
                 in1 >> crg;
                 in1 >> x;
                 in1 >> y;
                 in1 >> z;
                 
                 QMAtom *qmatom=new QMAtom(_elements.getEleName(atnum),x*Bohr2A,y*Bohr2A,z*Bohr2A,crg,false);
                 _atomlist->push_back(qmatom);
        }
        double potential=0.0;
        //cout << "File has " << xsteps << " : "<< ysteps << " : "<<zsteps << " gridpoints"<<endl;   
        for (int _ix = 0; _ix < _xsteps; _ix++) {
            double posx=(xstart+_ix*xincr)*Bohr2Nm;

           for (int _iy = 0; _iy < _ysteps; _iy++) {
              double posy=(ystart+_iy*yincr)*Bohr2Nm;


              for (int _iz = 0; _iz < _zsteps; _iz++) {
                double posz=(zstart+_iz*zincr)*Bohr2Nm;
                in1 >> potential;
                vec temp=vec(posx,posy,posz);
                ub::vector<double> temppos=temp.converttoub();
                APolarSite *apolarsite= new APolarSite(0,name);
                apolarsite->setRank(0);        
                apolarsite->setQ00(0,0); // <- charge state 0 <> 'neutral'
                apolarsite->setIsoP(0.0);
                apolarsite->setPos(temp);
                apolarsite->setPhi(potential,0);
                if(potential!=0.0 || !ignore_zeros){
                _gridsites.push_back(apolarsite);
                _gridpoints.push_back(temppos);
                }
                else {apolarsite->setIsVirtual(true);}
                 _all_gridsites.push_back(apolarsite);


              }}}
        if (_sites_seg != NULL) delete _sites_seg;
        _sites_seg = new PolarSeg(0, _gridsites);

        }         

void Grid::printgridtoCubefile(string filename){
            double A2Bohr=1.8897259886;
            //Creates Cube file of Grid in Angstrom and 
            if(!_cubegrid){
                throw std::runtime_error("Grid cannot be written to cube file as grid is not regular");
            }
            if(this->getTotalSize()<1){
                throw std::runtime_error("Grid object is empty. Setup grid first!");
            }
            
            //cout << _upperbound*A2Bohr<<endl;
            //cout << _lowerbound*A2Bohr<<endl;
            
            Elements _elements;
            FILE *out;
            out = fopen(filename.c_str(), "w");
            
            fprintf(out, "Electrostatic potential around molecule \n" );
            fprintf(out, "Created by VOTCA-XTP \n");
            fprintf(out, "%lu %f %f %f \n", _atomlist->size(), _lowerbound.getX()*A2Bohr, _lowerbound.getY()*A2Bohr,_lowerbound.getZ()*A2Bohr);
            fprintf(out, "%d %f 0.0 0.0 \n", _xsteps, _gridspacing*A2Bohr); 
            fprintf(out, "%d 0.0 %f 0.0 \n",  _ysteps, _gridspacing*A2Bohr);
            fprintf(out, "%d 0.0 0.0 %f \n", _zsteps, _gridspacing*A2Bohr);
            
            vector<QMAtom* >::const_iterator ait;
            for (ait=_atomlist->begin(); ait != _atomlist->end(); ++ait) {
                    
                    double x = (*ait)->x*A2Bohr;
                    double y = (*ait)->y*A2Bohr;
                    double z = (*ait)->z*A2Bohr;

                    string element = (*ait)->type;
                    int atnum = _elements.getEleNum (element);
                    double crg =_elements.getNucCrgECP(element);

                    fprintf(out, "%d %f %f %f %f\n", atnum, crg, x, y, z);
                }
            vector< APolarSite* >::iterator pit;
            int Nrecord=0.0;
            for(pit=_all_gridsites.begin();pit!=_all_gridsites.end();++pit){
                Nrecord++;

                double _potential=(*pit)->getPhi();
                if (Nrecord == 6) {
                    fprintf(out, "%E \n", _potential);
                    Nrecord = 0;
                } else {
                    fprintf(out, "%E ", _potential);
                }
                                
                
            }             
        
        fclose(out);        
        }    

//Create a 12^depth geodesic grid for a sphere of a given radius/cutoff

void Grid::subdivide(const vec &v1, const vec &v2, const vec &v3, std::vector<vec> &spherepoints, const int depth) {
    if(depth == 0) {
        spherepoints.push_back(v1);
        spherepoints.push_back(v2);
        spherepoints.push_back(v3);
        return;
    }
    const vec v12 = (v1 + v2).normalize();
    const vec v23 = (v2 + v3).normalize();
    const vec v31 = (v3 + v1).normalize();
    subdivide(v1, v12, v31, spherepoints, depth - 1);
    subdivide(v2, v23, v12, spherepoints, depth - 1);
    subdivide(v3, v31, v23, spherepoints, depth - 1);
    subdivide(v12, v23, v31,spherepoints, depth - 1);
}

void Grid::initialize_sphere(std::vector<vec> &spherepoints, const int depth) {
    const double X = 0.525731112119133606;
    const double Z = 0.850650808352039932;
    std::vector<vec> vdata;
    vdata.push_back(vec(-X, 0.0, Z));
    vdata.push_back(vec(X, 0.0, Z));
    vdata.push_back(vec(-X, 0.0, -Z));
    vdata.push_back(vec(X, 0.0, -Z));
    vdata.push_back(vec(0.0, Z, X));
    vdata.push_back(vec(0.0, Z, -X));
    vdata.push_back(vec(0.0, -Z, X));
    vdata.push_back(vec(0.0, -Z, -X));
    vdata.push_back(vec(Z, X, 0.0));
    vdata.push_back(vec(-Z, X, 0.0));
    vdata.push_back(vec(Z, -X, 0.0));
    vdata.push_back(vec(-Z, -X, 0.0));
    int tindices[20][3] = {
        {0, 4, 1}, { 0, 9, 4 }, { 9, 5, 4 }, { 4, 5, 8 }, { 4, 8, 1 },
        { 8, 10, 1 }, { 8, 3, 10 }, { 5, 3, 8 }, { 5, 2, 3 }, { 2, 7, 3 },
        { 7, 10, 3 }, { 7, 6, 10 }, { 7, 11, 6 }, { 11, 0, 6 }, { 0, 1, 6 },
        { 6, 1, 10 }, { 9, 0, 11 }, { 9, 11, 2 }, { 9, 2, 5 }, { 7, 2, 11 }
    };
    for(int i = 0; i < 20; i++)
        subdivide(vdata[tindices[i][0]], vdata[tindices[i][1]], vdata[tindices[i][2]], spherepoints, depth);
}



void Grid::setupradialgrid(const int depth) {
    _cubegrid = false;
    std::vector<vec> spherepoints;
    initialize_sphere(spherepoints, depth);
    double A2Nm = 0.1;
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    vector<QMAtom* >::const_iterator ait;
    for (ait = _atomlist->begin(); ait != _atomlist->end(); ++ait) {
        x += (*ait)->x;
        y += (*ait)->y;
        z += (*ait)->z;
    }

    double numofatoms = double(_atomlist->size());
    vec centerofmolecule = vec(x / numofatoms, y / numofatoms, z / numofatoms);
    //cout <<endl;
    //cout << "Center of molecule is "<< centerofmolecule<<"."<<endl;
    
    bool cutoff_smaller_molecule=false;
    double temp_cutoff=_cutoff;
    for (ait = _atomlist->begin(); ait != _atomlist->end(); ++ait) {
        x = (*ait)->x;
        y = (*ait)->y;
        z = (*ait)->z;
        
        vec temppos=vec(x,y,z);
        double dist=abs(centerofmolecule-temppos)+2.0; // 2.0 because that is approximately the VdW radius
        if (dist>temp_cutoff){
            cutoff_smaller_molecule=true;
            temp_cutoff=dist;
        }
    }
    
    if (cutoff_smaller_molecule){       
        cout <<endl;
        cout << "Specified cutoff of "<< _cutoff<<" to small. Cutoff is set to "<< temp_cutoff <<"Angstroem."<<endl;
        _cutoff=temp_cutoff;
    }
        
    //till here everything was Angstroem
        
    vector<vec>::const_iterator git;
    for (git = spherepoints.begin(); git != spherepoints.end(); ++git) {
        vec position=((*git)*_cutoff+centerofmolecule)*A2Nm;
        _gridpoints.push_back(position.converttoub());
        if(_createpolarsites){                   
            string name="H";
            APolarSite *apolarsite= new APolarSite(0,name);
            apolarsite->setRank(0);        
            apolarsite->setQ00(0,0); // <- charge state 0 <> 'neutral'
            apolarsite->setIsoP(0.0);
            apolarsite->setPos(position);
            _gridsites.push_back(apolarsite);
            _all_gridsites.push_back(apolarsite);       
        }
    }
    if (_sites_seg != NULL) delete _sites_seg;
           _sites_seg = new PolarSeg(0, _gridsites);

}

void Grid::setupgrid(){
           
            
    double A2Nm=0.1;
    Elements _elements;
    double xmin=std::numeric_limits<double>::max();
    double ymin=xmin;
    double zmin=xmin;

    double xmax=std::numeric_limits<double>::min();
    double ymax=xmax;
    double zmax=xmax;
    double xtemp,ytemp,ztemp;

    if(_useVdWcutoff){
        _padding=0.0;
        for (vector<QMAtom* >::const_iterator atom = _atomlist->begin(); atom != _atomlist->end(); ++atom ){
            if(_elements.getVdWChelpG((*atom)->type)+_shift_cutoff>_padding) _padding=_elements.getVdWChelpG((*atom)->type)+_shift_cutoff; 
        }
    } 


    for (vector<QMAtom* >::const_iterator atom = _atomlist->begin(); atom != _atomlist->end(); ++atom ) {
        xtemp=(*atom)->x;
        ytemp=(*atom)->y;
        ztemp=(*atom)->z;
        if (xtemp<xmin) xmin=xtemp;
        if (xtemp>xmax) xmax=xtemp;
        if (ytemp<ymin) ymin=ytemp;
        if (ytemp>ymax)  ymax=ytemp;
        if (ztemp<zmin) zmin=ztemp;
        if (ztemp>zmax)  zmax=ztemp;
    }    

    _lowerbound=vec(xmin-_padding,ymin-_padding,zmin-_padding);
    vec _upperbound=vec(xmax+_padding,ymax+_padding,zmax+_padding);
    vec steps=(_upperbound-_lowerbound)/_gridspacing;
    _xsteps=int(ceil(steps.getX()));
    _ysteps=int(ceil(steps.getY()));
    _zsteps=int(ceil(steps.getZ()));
    
    // needed to symmetrize grid around molecule
    double padding_x=(steps.getX()-_xsteps)*_gridspacing*0.5+_padding;
    double padding_y=(steps.getY()-_ysteps)*_gridspacing*0.5+_padding;
    double padding_z=(steps.getZ()-_zsteps)*_gridspacing*0.5+_padding;

    ub::vector<double> temppos= ub::zero_vector<double>(3);
    for(int i=0;i<=_xsteps;i++){
        double x=xmin-padding_x+i*_gridspacing; 
        for(int j=0;j<=_ysteps;j++){
            double y=ymin-padding_y+j*_gridspacing; 
            for(int k=0;k<=_zsteps;k++){
                double z=zmin-padding_z+k*_gridspacing; 
                bool _is_valid = false;
                    for (vector<QMAtom* >::const_iterator atom = _atomlist->begin(); atom != _atomlist->end(); ++atom ) {
                        //cout << "Punkt " << x <<":"<< y << ":"<<z << endl;
                        xtemp=(*atom)->x;
                        ytemp=(*atom)->y;
                        ztemp=(*atom)->z;
                        double distance2=pow((x-xtemp),2)+pow((y-ytemp),2)+pow((z-ztemp),2);
                        if(_useVdWcutoff) _cutoff=_elements.getVdWChelpG((*atom)->type)+_shift_cutoff;
                        if(_useVdWcutoff_inside)_cutoff_inside=_elements.getVdWChelpG((*atom)->type)+_shift_cutoff_inside;
                        //cout << "Punkt " << x <<":"<< y << ":"<<z << ":"<< distance2 << ":"<< (*atom)->type <<":"<<pow(VdW,2)<< endl;
                        if ( distance2<pow(_cutoff_inside,2)){
                            _is_valid = false;
                            break;
                            }
                        else if ( distance2<pow(_cutoff,2))  _is_valid = true;
                    }
                    if (_is_valid || _cubegrid){
                        temppos(0)=A2Nm*x;
                        temppos(1)=A2Nm*y;        
                        temppos(2)=A2Nm*z;   
                        if(_createpolarsites){
                            // APolarSite are in nm so convert
                            vec temp=vec(temppos);
                            string name="H";
                            APolarSite *apolarsite= new APolarSite(0,name);
                            apolarsite->setRank(0);        
                            apolarsite->setQ00(0,0); // <- charge state 0 <> 'neutral'
                            apolarsite->setIsoP(0.0);
                            apolarsite->setPos(temp);
                            if(_is_valid){
                                _gridsites.push_back(apolarsite);
                                _gridpoints.push_back(temppos);
                                }
                            else {apolarsite->setIsVirtual(true);}
                            _all_gridsites.push_back(apolarsite);
                            }
                        else if(!_createpolarsites){_gridpoints.push_back(temppos);}
                    }                    
                }                          
            }                  
        }
    if (_sites_seg != NULL) delete _sites_seg;
    _sites_seg = new PolarSeg(0, _gridsites);
}
    
    
}}