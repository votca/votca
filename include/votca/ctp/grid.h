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
        
        
        Grid(bool outsidemolecule, bool createpolarsites)
            :_cutoff(1),_padding(1),_gridspacing(1),_outsidemolecule(true),_createpolarsites(false) {};
           
        
        Grid()
            :_cutoff(1),_padding(1),_gridspacing(1),_outsidemolecule(true),_createpolarsites(false) {};
           
        
        ~Grid() {};
        
        std::vector< ub::vector<double> > &getGrid() {return _gridpoints;}
        std::vector< APolarSite > &getSites() {return _gridsites;}
        void setCutoff(double cutoff){_cutoff=cutoff;}
        void setSpacing(double spacing){_gridspacing=spacing;}
        void setPadding(double padding){_padding=padding;}
    
      
        
        
        int getsize(){ return _gridpoints.size(); }
        
        void printGridtofile(const char* _filename){
            ofstream points;
            points.open(_filename, ofstream::out);
            points << _gridpoints.size() << endl;
            points << endl;
            for ( int i = 0 ; i < _gridpoints.size(); i++){
                points << "X " << _gridpoints[i](0) << " " << _gridpoints[i](1) << " " << _gridpoints[i](2) << endl;

            }
            points.close();
        }
        
       
        
        
        void setupgrid(const vector< QMAtom* >& Atomlist){
            Elements _elements;
            double xmin=1000;
            double ymin=1000;
            double zmin=1000;

            double xmax=-1000;
            double ymax=-1000;
            double zmax=-1000;
            double xtemp,ytemp,ztemp;
            //setup one polarsite and use copy constructor later
         
        
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
                                    double VdW=_elements.getVdWChelpG((*atom)->type);
                                    //cout << "Punkt " << x <<":"<< y << ":"<<z << ":"<< distance2 << ":"<< (*atom)->type <<":"<<pow(VdW,2)<< endl;
                                    if ( _outsidemolecule && distance2<pow(VdW,2)){
                                        //cout << "Punkt" << x <<":"<< y << ":"<<z << "rejected" << endl;

                                        _is_valid = false;
                                        break;
                                        }
                                    else if (distance2<pow(_cutoff,2)){
                                        //cout << "hier" << endl;
                                        _is_valid = true;
                                    }



                                }
                            if (_is_valid){
                                temppos(0)=x;
                                temppos(1)=y;        
                                temppos(2)=z;
                                _gridpoints.push_back(temppos);
                                if(_createpolarsites){
                                    vec temp=vec(x,y,z);
                                    APolarSite apolarsite;
                                    apolarsite.setPos(temp);
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
           
        
        
        
        
        }
       
        void setupCHELPgrid(const vector< QMAtom* >& Atomlist){
            

            _padding=2.8; // Additional distance from molecule to set up grid according to CHELPG paper [Journal of Computational Chemistry 11, 361, 1990]
            _gridspacing=0.3; // Grid spacing according to same paper 
            _cutoff=2.8;
            _outsidemolecule=true;
            _createpolarsites=false;
            
            setupgrid(Atomlist);
        }
        
        
            

               
                        
                        
        
      
  private:
      std::vector< ub::vector<double> > _gridpoints;
      std::vector< APolarSite > _gridsites;
      
      double _gridspacing;
      double _cutoff;
      double _padding;
      bool   _outsidemolecule;
      bool   _createpolarsites;
      
        
    };   
    
 
    
}}

#endif	/* GRID_H */