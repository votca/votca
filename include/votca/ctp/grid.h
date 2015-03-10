/* 
 *            Copyright 2009-2012 The VOTCA Development Team
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

#ifndef __CTP_GRID__H
#define	__CTP_GRID__H


#include <votca/ctp/elements.h>
#include <string>
#include <map>
#include <vector>
#include <votca/tools/property.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/multi_array.hpp>
#include <votca/ctp/qmatom.h>
#include <votca/ctp/logger.h>
#include <votca/ctp/apolarsite.h>
/**
* \brief Takes a list of atoms, and creates different grids for it. Right now only CHELPG grid.
*
* 
* 
*/
using namespace std;
using namespace votca::tools;


namespace votca { namespace ctp {
    namespace ub = boost::numeric::ublas;
    
  class Grid{
    public:
        
        
        Grid( bool createpolarsites, bool useVdWcutoff, bool useVdWcutoff_inside)
            :_cutoff(3),_gridspacing(0.3),_cutoff_inside(1.5),_shift_cutoff(0.0),_shift_cutoff_inside(0.0),
             _useVdWcutoff(useVdWcutoff),_useVdWcutoff_inside(useVdWcutoff_inside),_cubegrid(false),_padding(3.0),
             _createpolarsites(createpolarsites), _sites_seg(NULL) {};
           
        
        Grid()
            :_cutoff(3),_gridspacing(0.3),_cutoff_inside(1.5),_shift_cutoff(0.0),_shift_cutoff_inside(0.0),
             _useVdWcutoff(false),_useVdWcutoff_inside(false),_cubegrid(false),_padding(3.0),
             _createpolarsites(false), _sites_seg(NULL) {};
           
        
        ~Grid() {};
        
        std::vector< ub::vector<double> > &getGrid() {return _gridpoints;}
        std::vector< APolarSite* > &Sites() {return _gridsites;}
        std::vector< APolarSite*>* getSites() {return &_gridsites;} 
        PolarSeg* getSeg(){return _sites_seg;}

        
        void setCutoffs(double cutoff, double cutoff_inside){_cutoff=cutoff;_cutoff_inside=cutoff_inside;}
        void setCutoffshifts(double shift_cutoff, double shift_cutoff_inside){_shift_cutoff=shift_cutoff;_shift_cutoff_inside=shift_cutoff_inside;}
        void setSpacing(double spacing){_gridspacing=spacing;}
        void setPadding(double padding){_padding=padding;}
        void setCubegrid(bool cubegrid){_cubegrid=cubegrid;_createpolarsites=true;}
        void setAtomlist(vector< QMAtom* >* Atomlist){_atomlist=Atomlist;}
        int  getTotalSize(){return _gridpoints.size();}
        
        int getsize(){
            int size=0.0;
            if(_cubegrid){size=_gridsites.size();}
            else{size=_gridpoints.size();}

            return size; 
        }

        void printGridtofile(const char* _filename){
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
        
        
       
        void readgridfromCubeFile(string filename, bool ignore_zeros){
           
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
        
        double val1;
            
        for (int _ix = 0; _ix < xsteps; _ix++) {
            double posx=(xstart+_ix*xincr)*Bohr2Nm;

           for (int _iy = 0; _iy < ysteps; _iy++) {
              double posy=(ystart+_iy*yincr)*Bohr2Nm;


              for (int _iz = 0; _iz < zsteps; _iz++) {
                double posz=(zstart+_iz*zincr)*Bohr2Nm;
                in1 >> val1;
                vec temp=vec(posx,posy,posz);
                ub::vector<double> temppos=temp.converttoub();
                APolarSite *apolarsite= new APolarSite(0,name);
                apolarsite->setRank(0);        
                apolarsite->setQ00(0,0); // <- charge state 0 <> 'neutral'
                apolarsite->setIsoP(0.0);
                apolarsite->setPos(temp);
                if(val1!=0.0 || !ignore_zeros){
                _gridsites.push_back(apolarsite);
                _gridpoints.push_back(temppos);
                }
                else {apolarsite->setIsVirtual(true);}
                 _all_gridsites.push_back(apolarsite);


              }}}
        if (_sites_seg != NULL) delete _sites_seg;
        _sites_seg = new PolarSeg(0, _gridsites);

        }         
        
        
        void printgridtoCubefile(string filename){
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
        
        
        
        
        
        
        
  
        
        //setup will return a grid in nm not in A, although setupgrid internally uses A.
        void setupgrid(){
           
            
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
       
        void setupCHELPgrid(){
            

            //_padding=2.8; // Additional distance from molecule to set up grid according to CHELPG paper [Journal of Computational Chemistry 11, 361, 1990]
            _gridspacing=0.3; // Grid spacing according to same paper 
            _cutoff=2.8;
            _useVdWcutoff_inside=true;
            _shift_cutoff_inside=0.0;
            _useVdWcutoff=false;
            setupgrid();
        }
        
        
            

               
                        
                        
        
      
  private:
      std::vector< ub::vector<double> > _gridpoints;
      std::vector< APolarSite* > _gridsites;
      std::vector< APolarSite* > _all_gridsites;
      PolarSeg *_sites_seg;
      const vector< QMAtom* >* _atomlist;
      double _gridspacing;
      double _cutoff;
      double _cutoff_inside;
      double _shift_cutoff_inside;
      double _shift_cutoff;
      double _padding;
      bool   _createpolarsites;
      bool   _useVdWcutoff;
      bool   _useVdWcutoff_inside;
      bool   _cubegrid;
      vec _upperbound;
      vec _lowerbound;
      
        
    };   
    
 
    
}}

#endif	/* GRID_H */