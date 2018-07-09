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

#include <votca/xtp/grid.h>
#include <math.h>       /* ceil */
#include <votca/tools/constants.h>




namespace votca { namespace xtp {
  using namespace tools;
    
Grid::Grid(const Grid &obj)
    :_cutoff(obj._cutoff),_gridspacing(obj._gridspacing),_cutoff_inside(obj._cutoff_inside),_shift_cutoff(obj._shift_cutoff),
    _shift_cutoff_inside(obj._shift_cutoff_inside),_useVdWcutoff(obj._useVdWcutoff),_useVdWcutoff_inside(obj._useVdWcutoff_inside),
    _cubegrid(obj._cubegrid),_padding(obj._padding),_createpolarsites(obj._createpolarsites),
    _xsteps(obj._xsteps),_ysteps(obj._ysteps),_zsteps(obj._zsteps) {
    _lowerbound=obj._lowerbound;
    _gridpoints=obj._gridpoints;
    std::vector<xtp::APolarSite*>::const_iterator pit;
    for(pit=obj._all_gridsites.begin();pit<obj._all_gridsites.end();++pit){
       xtp::APolarSite *apolarsite= new xtp::APolarSite(*pit,false);
       if(!apolarsite->getIsVirtual()) _gridsites.push_back(apolarsite);
       _all_gridsites.push_back(apolarsite);   
    }     
     _atomlist=obj._atomlist;
    };
        
        
Grid::~Grid() {
        std::vector<xtp::APolarSite*>::iterator pit;
        for(pit=_all_gridsites.begin();pit!=_all_gridsites.end();++pit){
             delete *pit;
        }
        _all_gridsites.clear();
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
    std::vector<xtp::APolarSite*>::const_iterator pit;
    for(pit=obj._all_gridsites.begin();pit<obj._all_gridsites.end();++pit){
       xtp::APolarSite *apolarsite= new xtp::APolarSite(*pit,false);
       if(!apolarsite->getIsVirtual()) _gridsites.push_back(apolarsite);
       _all_gridsites.push_back(apolarsite);   
    }     
     _atomlist=obj._atomlist;
     return *this;
}
    
void Grid::printGridtoxyzfile(std::string filename){
        //unit is Angstrom in xyz file 
        std::ofstream points;
        points.open(filename.c_str(), std::ofstream::out);
        points << _gridpoints.size() << endl;
        points << endl;
        for ( const auto& point:_gridpoints){
            points << "X " << point.getX()*conv::bohr2ang << " " 
                    << point.getY()*conv::bohr2ang << " " 
                    << point.getZ()*conv::bohr2ang << endl;

        }
        points.close();
        return;
    }    


void Grid::readgridfromCubeFile(std::string filename, bool ignore_zeros){
        Elements _elements;
        _cubegrid=true;   

        if(_gridpoints.size()>0) throw std::runtime_error("Grid object already has points.");
        std::ifstream in1;
        std::string s;
        in1.open(filename.c_str(), std::ios::in);
        getline(in1, s);
        getline(in1, s);
        int natoms;
        double xstart,ystart,zstart;
        double xincr,yincr,zincr;

        double tempdouble;
        std::string name="H";
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
        
        
        if(xincr==yincr && yincr==zincr && zincr==xincr){
          _gridspacing=xincr*conv::bohr2ang;
        }
        else{
          throw std::runtime_error("Gridspacing in x,y,z is different, currently not implemented, loading aborted");
        }
        _lowerbound=vec(xstart*conv::bohr2ang,ystart*conv::bohr2ang,zstart*conv::bohr2ang);
        
        _atomlist= new std::vector< QMAtom* >;
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
                 
                 QMAtom *qmatom=new QMAtom(iatom,_elements.getEleName(atnum),x,y,z);
                 _atomlist->push_back(qmatom);
        }
        double potential=0.0;
        //cout << "File has " << xsteps << " : "<< ysteps << " : "<<zsteps << " gridpoints"<<endl;   
        for (int _ix = 0; _ix < _xsteps; _ix++) {
            double posx=(xstart+_ix*xincr);

           for (int _iy = 0; _iy < _ysteps; _iy++) {
              double posy=(ystart+_iy*yincr);


              for (int _iz = 0; _iz < _zsteps; _iz++) {
                double posz=(zstart+_iz*zincr);
                in1 >> potential;
                vec pos=vec(posx,posy,posz);
                xtp::APolarSite *apolarsite= new xtp::APolarSite(0,name);
                apolarsite->setRank(0);        
                apolarsite->setQ00(0,0); // <- charge state 0 <> 'neutral'
                apolarsite->setIsoP(0.0);
                vec temp=pos*tools::conv::bohr2nm;
                apolarsite->setPos(temp);
                apolarsite->setPhi(potential,0);
                if(potential!=0.0 || !ignore_zeros){
                _gridsites.push_back(apolarsite);
                _gridpoints.push_back(pos);
                }
                else {apolarsite->setIsVirtual(true);}
                 _all_gridsites.push_back(apolarsite);


              }}}
        
        return;
        }         

void Grid::printgridtoCubefile(std::string filename){

            //Creates Cube file of Grid in Angstrom and 
            if(!_cubegrid){
                throw std::runtime_error("Grid cannot be written to cube file as grid is not regular");
            }
            if(this->getTotalSize()<1){
                throw std::runtime_error("Grid object is empty. Setup grid first!");
            }
          
            
            Elements _elements;
            FILE *out;
            out = fopen(filename.c_str(), "w");
            
            fprintf(out, "Electrostatic potential around molecule \n" );
            fprintf(out, "Created by VOTCA-XTP \n");
            fprintf(out, "%lu %f %f %f \n", _atomlist->size(), _lowerbound.getX(),
                    _lowerbound.getY(),_lowerbound.getZ());
            fprintf(out, "%d %f 0.0 0.0 \n", _xsteps, _gridspacing); 
            fprintf(out, "%d 0.0 %f 0.0 \n",  _ysteps, _gridspacing);
            fprintf(out, "%d 0.0 0.0 %f \n", _zsteps, _gridspacing);
            
            std::vector<QMAtom* >::const_iterator ait;
            for (ait=_atomlist->begin(); ait != _atomlist->end(); ++ait) {
              const tools::vec& pos=(*ait)->getPos();
                   
                    std::string element = (*ait)->getType();
                    int atnum = _elements.getEleNum (element);
                    int crg=(*ait)->getNuccharge();

                    fprintf(out, "%d %d %f %f %f\n", atnum, crg, pos.getX(), pos.getY(), pos.getZ());
                }
            std::vector< xtp::APolarSite* >::iterator pit;
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
        return;
        }    

void Grid::setupgrid(){
           
            
    Elements _elements;
    std::vector<QMAtom* >::const_iterator atom;
    if (_useVdWcutoff) {
        _padding = 0.0;
        for (atom = _atomlist->begin(); atom != _atomlist->end(); ++atom) {
          double VdWcutoff=_elements.getVdWChelpG((*atom)->getType())*tools::conv::ang2bohr;
          if ( VdWcutoff+ _shift_cutoff > _padding) {
            _padding = VdWcutoff + _shift_cutoff;
          }
        }
      }
    
    double xmin=std::numeric_limits<double>::max();
    double ymin=xmin;
    double zmin=xmin;

    double xmax=std::numeric_limits<double>::min();
    double ymax=xmax;
    double zmax=xmax;
    double xtemp,ytemp,ztemp;

    for (atom = _atomlist->begin(); atom != _atomlist->end(); ++atom ) {
      const tools::vec & pos=(*atom)->getPos();
        xtemp=pos.getX();
        ytemp=pos.getY();
        ztemp=pos.getZ();
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
    
    
    for(int i=0;i<=_xsteps;i++){
        double x=xmin-padding_x+i*_gridspacing; 
        for(int j=0;j<=_ysteps;j++){
            double y=ymin-padding_y+j*_gridspacing; 
            for(int k=0;k<=_zsteps;k++){
                double z=zmin-padding_z+k*_gridspacing; 
                bool _is_valid = false;
                vec gridpos=vec(x,y,z);
                    for (atom = _atomlist->begin(); atom != _atomlist->end(); ++atom ) {
                        vec atompos=(*atom)->getPos();
                        double distance2=(gridpos-atompos)*(gridpos-atompos);
                        if(_useVdWcutoff){
                          _cutoff=_elements.getVdWChelpG((*atom)->getType())*tools::conv::ang2bohr+_shift_cutoff;
                        }
                        if(_useVdWcutoff_inside){
                          _cutoff_inside=_elements.getVdWChelpG((*atom)->getType())*tools::conv::ang2bohr+_shift_cutoff_inside;
                        }
                        if ( distance2<(_cutoff_inside*_cutoff_inside)){
                            _is_valid = false;
                            break;
                            }
                        else if ( distance2<(_cutoff*_cutoff))  _is_valid = true;
                    }
                    if (_is_valid || _cubegrid){
                        
                        if(_createpolarsites){
                          
                            std::string name="X";
                            xtp::APolarSite *apolarsite= new xtp::APolarSite(0,name);
                            apolarsite->setRank(0);        
                            apolarsite->setQ00(0,0); // <- charge state 0 <> 'neutral'
                            apolarsite->setIsoP(0.0);
                            vec temp=gridpos*conv::bohr2nm;
                            apolarsite->setPos(temp);
                            if(_is_valid){
                                _gridsites.push_back(apolarsite);
                                _gridpoints.push_back(gridpos);
                                }
                            else {apolarsite->setIsVirtual(true);}
                            _all_gridsites.push_back(apolarsite);
                            }
                        else if(!_createpolarsites){_gridpoints.push_back(gridpos);}
                    }                    
                }                          
            }                  
        }
    
    _gridvalues=Eigen::VectorXd::Zero(_gridpoints.size());
    return;
}
    
    
}}
