/* 
 *            Copyright 2009-2018 The VOTCA Development Team
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

#ifndef __XTP_GRID__H
#define	__XTP_GRID__H


#include <votca/tools/elements.h>
#include <string>
#include <vector>
#include <votca/xtp/qmatom.h>
#include <votca/xtp/logger.h>
#include <votca/xtp/apolarsite.h>
#include <votca/xtp/polarseg.h>
/**
* \brief Takes a list of atoms, and creates different grids for it. Right now only CHELPG grid.
*
* 
* 
*/


namespace votca { namespace xtp {
   
  class Grid{
    public:
        
        
        Grid( bool createpolarsites, bool useVdWcutoff, bool useVdWcutoff_inside)
            :_cutoff(3),_gridspacing(0.3),_cutoff_inside(1.5),_shift_cutoff(0.0),_shift_cutoff_inside(0.0),
             _useVdWcutoff(useVdWcutoff),_useVdWcutoff_inside(useVdWcutoff_inside),_cubegrid(false),_padding(3.0),
             _createpolarsites(createpolarsites), _atomlist(NULL), 
            _lowerbound(tools::vec(0,0,0)), _xsteps(0),_ysteps(0),_zsteps(0) {};
           
        
        Grid()
            :_cutoff(3),_gridspacing(0.3),_cutoff_inside(1.5),_shift_cutoff(0.0),_shift_cutoff_inside(0.0),
             _useVdWcutoff(false),_useVdWcutoff_inside(false),_cubegrid(false),_padding(3.0),
             _createpolarsites(false), _atomlist(NULL),
             _lowerbound(tools::vec(0,0,0)),_xsteps(0),_ysteps(0),_zsteps(0) {};
             
        Grid(std::vector< tools::vec > points)
             :_gridpoints(points){};
           
        
        ~Grid();
        
        Grid(const Grid &obj);
        
        Grid& operator=(const Grid &obj);
        
        const std::vector< tools::vec > &getGridPositions() const {return _gridpoints;}
        Eigen::VectorXd &getGridValues(){return _gridvalues;}
        const Eigen::VectorXd &getGridValues() const{return _gridvalues;}
        std::vector< xtp::APolarSite* > &Sites() {return _gridsites;}
        
        void setCutoffs(double cutoff, double cutoff_inside){_cutoff=cutoff;_cutoff_inside=cutoff_inside;}
        
        void setAtomlist(std::vector< QMAtom* >* Atomlist){_atomlist=Atomlist;}
        unsigned getsize(){return _gridpoints.size();}
        
        int getTotalSize(){
            int size=0.0;
            if(_cubegrid){size=_all_gridsites.size();}
            else{size=_gridpoints.size();}
            return size; 
        }

        void printGridtoxyzfile(std::string filename);
        
        void readgridfromCubeFile(std::string filename, bool ignore_zeros=true);
       
        void printgridtoCubefile(std::string filename);
        
        
        void setupgrid();
       
        void setupCHELPgrid(){
            _padding=3*tools::conv::ang2bohr; // Additional distance from molecule to set up grid according to CHELPG paper [Journal of Computational Chemistry 11, 361, 1990]
            _gridspacing=0.3*tools::conv::ang2bohr; // Grid spacing according to same paper 
            _cutoff=2.8*tools::conv::ang2bohr;
            _useVdWcutoff_inside=true;
            _shift_cutoff_inside=0.0;
            _useVdWcutoff=false;
            setupgrid();
        }
       
      
  private:
     
      std::vector< tools::vec > _gridpoints;
      Eigen::VectorXd _gridvalues;
      std::vector< xtp::APolarSite* > _gridsites;
      std::vector< xtp::APolarSite* > _all_gridsites;
      
      double _cutoff;
      double _gridspacing;
      double _cutoff_inside;
      double _shift_cutoff;
      double _shift_cutoff_inside;
      bool   _useVdWcutoff;
      bool   _useVdWcutoff_inside;
      bool   _cubegrid;
      double _padding;
      bool   _createpolarsites; 
      std::vector< QMAtom* >* _atomlist;
      tools::vec _lowerbound;
      int _xsteps, _ysteps, _zsteps;
      
 
    };   
    
 
    
}}

#endif	/* GRID_H */
