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
             _useVdWcutoff(useVdWcutoff),_useVdWcutoff_inside(useVdWcutoff_inside),
             _createpolarsites(createpolarsites), _sites_seg(NULL) {};
           
        
        Grid()
            :_cutoff(3),_gridspacing(0.3),_cutoff_inside(1.5),_shift_cutoff(0.0),_shift_cutoff_inside(0.0),
             _useVdWcutoff(false),_useVdWcutoff_inside(false),
             _createpolarsites(false), _sites_seg(NULL) {};
           
        
        ~Grid() {};
        
        std::vector< ub::vector<double> > &getGrid() {return _gridpoints;}
        std::vector< APolarSite* > &Sites() {return _gridsites;}
        std::vector< APolarSite*>* getSites() {return &_gridsites;} 
        PolarSeg* getSeg(){return _sites_seg;}
        
        
        void setCutoffs(double cutoff, double cutoff_inside){_cutoff=cutoff;_cutoff_inside=cutoff_inside;}
        void setCutoffshifts(double shift_cutoff, double shift_cutoff_inside){_shift_cutoff=shift_cutoff;_shift_cutoff_inside=shift_cutoff_inside;}
        void setSpacing(double spacing){_gridspacing=spacing;}
    
      
        
        
        int getsize(){ return _gridpoints.size(); }
        
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
        
       
        void generateCubefile(const char* _filename){
        
        
        
        
        
        }
        
        
        
        
  
        
        //setup will return a grid in nm not in A, although setupgrid internally uses A.
        void setupgrid(const vector< QMAtom* >& Atomlist){
           
            double _padding=0.0;
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
            for (vector<QMAtom* >::const_iterator atom = Atomlist.begin(); atom != Atomlist.end(); ++atom ){
                if(_elements.getVdWChelpG((*atom)->type)+_shift_cutoff>_padding) _padding=_elements.getVdWChelpG((*atom)->type)+_shift_cutoff; 
            }
        } 
            else _padding=_cutoff;
                
                
         for (vector<QMAtom* >::const_iterator atom = Atomlist.begin(); atom != Atomlist.end(); ++atom ) {
                xtemp=(*atom)->x;
                ytemp=(*atom)->y;
                ztemp=(*atom)->z;
                if (xtemp<xmin)
                    xmin=xtemp;
                if (xtemp>xmax)
                    xmax=xtemp;
                 if (ytemp<ymin)
                    ymin=ytemp;
                if (ytemp>ymax)
                    ymax=ytemp;
                 if (ztemp<zmin)
                    zmin=ztemp;
                if (ztemp>zmax)
                    zmax=ztemp;

            }    

                double boxdimx=xmax-xmin+2*_padding;
               

                double x=xmin-_padding;


                ub::vector<double> temppos= ub::zero_vector<double>(3);
                while(x< xmax+_padding){
                   double y=ymin-_padding;
                   while(y< ymax+_padding){
                        double z=zmin-_padding;
                        while(z< zmax+_padding){
                            bool _is_valid = false;
                                for (vector<QMAtom* >::const_iterator atom = Atomlist.begin(); atom != Atomlist.end(); ++atom ) {
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
                            if (_is_valid){
                                temppos(0)=AtoNm*x;
                                temppos(1)=AtoNm*y;        
                                temppos(2)=AtoNm*z;
                                _gridpoints.push_back(temppos);
                                if(_createpolarsites){
                                    // APolarSite are in nm so convert
                                    vec temp=vec(temppos);
                                    string name="H";
                                    APolarSite *apolarsite= new APolarSite(0,name);
                                    apolarsite->setRank(0);        
                                    apolarsite->setQ00(0,0); // <- charge state 0 <> 'neutral'
                                    apolarsite->setIsoP(0.0);
                                    apolarsite->setPos(temp);
                                    _gridsites.push_back(apolarsite);
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
       
        void setupCHELPgrid(const vector< QMAtom* >& Atomlist){
            

            //_padding=2.8; // Additional distance from molecule to set up grid according to CHELPG paper [Journal of Computational Chemistry 11, 361, 1990]
            _gridspacing=0.3; // Grid spacing according to same paper 
            _cutoff=2.8;
            _useVdWcutoff_inside=true;
            _shift_cutoff_inside=0.0;
            _useVdWcutoff=false;
            setupgrid(Atomlist);
        }
        
        
            

               
                        
                        
        
      
  private:
      std::vector< ub::vector<double> > _gridpoints;
      std::vector< APolarSite* > _gridsites;
      PolarSeg *_sites_seg;
      
      double _gridspacing;
      double _cutoff;
      double _cutoff_inside;
      double _shift_cutoff_inside;
      double _shift_cutoff;
      bool   _createpolarsites;
      bool   _useVdWcutoff;
      bool   _useVdWcutoff_inside;

      
        
    };   
    
 
    
}}

#endif	/* GRID_H */