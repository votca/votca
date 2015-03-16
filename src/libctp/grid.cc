
#include <votca/ctp/grid.h>

using namespace std;
using namespace votca::tools;


namespace votca { namespace ctp {
    namespace ub = boost::numeric::ublas;
    
void Grid::printGridtofile(const char* _filename){
            //unit is nm
            double A2nm=10;
            ofstream points;
            points.open(_filename, ofstream::out);
            points << _gridpoints.size() << endl;
            points << endl;
            for ( int i = 0 ; i < _gridpoints.size(); i++){
                points << "X " << A2nm*_gridpoints[i](0) << " " << A2nm*_gridpoints[i](1) << " " << A2nm*_gridpoints[i](2) << endl;

            }
            points.close();
        }    


void Grid::readgridfromCubeFile(string filename, bool ignore_zeros){
        _cubegrid=true;   
        double Bohr2Nm=1.0/18.897259886; 
        if(_gridpoints.size()>0) throw std::runtime_error("Grid object already has points.");
        ifstream in1;
        string s;
        in1.open(filename.c_str(), ios::in);
        getline(in1, s);
        getline(in1, s);
        int natoms;
        double xstart,ystart,zstart;
        double xincr,yincr,zincr;
        double xsteps,ysteps,zsteps;
        double tempdouble;
        string name="H";
        in1 >> natoms;
        in1 >> xstart;
        in1 >> ystart;
        in1 >> zstart;
        in1 >> xsteps;
        in1 >> xincr;
        in1 >> tempdouble;
        in1 >> tempdouble;
        in1 >> ysteps;
        in1 >> tempdouble;
        in1 >> yincr;
        in1 >> tempdouble;
        in1 >> zsteps;
        in1 >> tempdouble;
        in1 >> tempdouble;
        in1 >> zincr;          
        
        if(xincr==yincr && yincr==zincr && zincr==xincr) _gridspacing==xincr*Bohr2Nm;
        else throw std::runtime_error("Gridspacing in x,y,z is different, currently not implemented, loading aborted");
        
        
        double potential;
            
        for (int _ix = 0; _ix < xsteps; _ix++) {
            double posx=(xstart+_ix*xincr)*Bohr2Nm;

           for (int _iy = 0; _iy < ysteps; _iy++) {
              double posy=(ystart+_iy*yincr)*Bohr2Nm;


              for (int _iz = 0; _iz < zsteps; _iz++) {
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
            if(_gridpoints.size()<1){
                throw std::runtime_error("Grid object is empty. Setup grid first!");
            }
            vec steps=(_upperbound-_lowerbound)/_gridspacing;
            //cout << _upperbound*A2Bohr<<endl;
            //cout << _lowerbound*A2Bohr<<endl;
            //cout << steps<<endl;
            
            Elements _elements;
            FILE *out;
            out = fopen(filename.c_str(), "w");
            
            fprintf(out, "Electrostatic potential around molecule \n" );
            fprintf(out, "Created by VOTCA-CTP \n");
            fprintf(out, "%d %f %f %f \n", _atomlist->size(), _lowerbound.getX()*A2Bohr, _lowerbound.getY()*A2Bohr,_lowerbound.getZ()*A2Bohr);
            fprintf(out, "%d %f 0.0 0.0 \n", int(steps.getX())+1, _gridspacing*A2Bohr);
            fprintf(out, "%d 0.0 %f 0.0 \n",  int(steps.getY())+1, _gridspacing*A2Bohr);
            fprintf(out, "%d 0.0 0.0 %f \n", int(steps.getZ())+1, _gridspacing*A2Bohr);
            
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
void Grid::setupradialgrid(int depth){}

void Grid::setupgrid(){
           
            
            double AtoNm=0.1;
            Elements _elements;
            double xmin=std::numeric_limits<double>::max();
            double ymin=xmin;
            double zmin=xmin;

            double xmax=std::numeric_limits<double>::min();
            double ymax=xmax;
            double zmax=xmax;
            double xtemp,ytemp,ztemp;
            //setup one polarsite and use copy constructor later
            
            
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

                double boxdimx=xmax-xmin+2*_padding;
               

                double x=xmin-_padding;
                _lowerbound=vec(xmin-_padding,ymin-_padding,zmin-_padding);
                _upperbound=vec(xmax+_padding,ymax+_padding,zmax+_padding);
                

                ub::vector<double> temppos= ub::zero_vector<double>(3);
                while(x< xmax+_padding){
                   double y=ymin-_padding;
                   while(y< ymax+_padding){
                        double z=zmin-_padding;
                        while(z< zmax+_padding){
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
                                temppos(0)=AtoNm*x;
                                temppos(1)=AtoNm*y;        
                                temppos(2)=AtoNm*z;   
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
                            }
                            z+=_gridspacing; 
                        }
                        y+=_gridspacing;
                     //cout << "Punkt " << x  <<":"<<  xmax+padding <<":"<< y << ":"<<z << endl;
                    }
                  x+=_gridspacing;
                  //cout << (x<xmax+padding) << endl;     
                }
              
                if (_sites_seg != NULL) delete _sites_seg;
                _sites_seg = new PolarSeg(0, _gridsites);
        }
    
    
}}