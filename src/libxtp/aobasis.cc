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
#include "votca/xtp/aobasis.h"
#include "votca/xtp/aoshell.h"
#include "votca/xtp/qmatom.h"
#include "votca/tools/elements.h"
#include "votca/xtp/aomatrix.h"
#include <votca/tools/constants.h>




namespace votca { namespace xtp {

 AOBasis::~AOBasis() {
        for (AOShell* shell:_aoshells) delete shell;
        _aoshells.clear();
         }

AOShell* AOBasis::addShell( const Shell& shell, const QMAtom& atom, int startIndex ){
        AOShell* aoshell = new AOShell( shell, atom, startIndex );
        _aoshells.push_back(aoshell);
        return aoshell;
        }

AOShell* AOBasis::addECPShell( const Shell& shell, const QMAtom& atom, int startIndex, bool nonlocal ){
        AOShell* aoshell = new AOShell( shell, atom, startIndex,nonlocal );
        _aoshells.push_back(aoshell);
        return aoshell;
        }

void AOBasis::ReorderMOs(Eigen::MatrixXd &v, const std::string& start, const std::string& target) {

    if (start == target) {
        return;
    }
    
    if(target=="orca" || target=="nwchem"){
        std::vector<int> multiplier = getMultiplierVector(target,start);
        // and reorder rows of _orbitals->_mo_coefficients() accordingly
        MultiplyMOs(v, multiplier);
    }

    // get reordering vector _start -> target
    std::vector<int> order = getReorderVector(start, target);
    
    // Sanity check
    if (v.rows() != int(order.size())) {
        std::cerr << "Size mismatch in ReorderMOs" << v.rows() << ":" << order.size() << std::endl;
        throw std::runtime_error("Abort!");
    }

    // actual swapping of coefficients
    for (unsigned _i_orbital = 0; _i_orbital < v.cols(); _i_orbital++) {
        for (unsigned s = 1, d; s < order.size(); ++s) {
            for (d = order[s]; d < s; d = order[d]) {
                ;
            }
            if (d == s) while (d = order[d], d != s) std::swap(v(s,_i_orbital), v(d,_i_orbital));
        }
    }

    // NWChem has some strange minus in d-functions
    if (start == "nwchem" || start == "orca") {
        std::vector<int> multiplier = getMultiplierVector(start, target);
        // and reorder rows of _orbitals->_mo_coefficients() accordingly
        MultiplyMOs(v, multiplier);
    }
    return;
}

void AOBasis::ReorderMatrix(Eigen::MatrixXd &v,const std::string& start,const std::string& target ){
    if (start==target){
        return;
    }
    std::vector<int> order = getReorderVector(start, target);
    std::vector<int> multiplier=getMultiplierVector(start,target);
    
     if (v.cols() != int(order.size())) {
        std::cerr << "Size mismatch in ReorderMatrix" << v.cols() << ":" << order.size() << std::endl;
        throw std::runtime_error("Abort!");
    }

    if (start != "xtp") {
      std::vector<int>newmultiplier=std::vector<int>(multiplier.size());
      for(unsigned i=0;i<newmultiplier.size();i++){
        newmultiplier[i]=multiplier[order[i]];
      }
    multiplier=newmultiplier;
    }
    
    Eigen::MatrixXd temp=v;
    for(unsigned i=0;i<temp.cols();i++){
        int i_index=order[i];
        for(unsigned j=0;j<temp.rows();j++){
            int j_index=order[j];
            v(i_index,j_index)=multiplier[i]*multiplier[j]*temp(i,j);
        }
    }
    
    return;
}


void AOBasis::MultiplyMOs(Eigen::MatrixXd &v, std::vector<int> const &multiplier )  {
          // Sanity check
          if ( v.cols() != int(multiplier.size()) ) {
              std::cerr << "Size mismatch in MultiplyMOs" << v.cols() << ":" << multiplier.size() << std::endl;
              throw std::runtime_error( "Abort!");
          }
          for ( unsigned _i_basis = 0; _i_basis < v.cols(); _i_basis++ ){
            for ( unsigned _i_orbital = 0; _i_orbital < v.rows(); _i_orbital++ ){               
                   v(_i_basis ,_i_orbital) = multiplier[_i_basis] * v(_i_basis,_i_orbital  );
               }
           }
          return;
    }


    //this is for gaussian only to transform from gaussian ordering cartesian to gaussian spherical
    Eigen::MatrixXd AOBasis::getTransformationCartToSpherical(const std::string& package) {
      Eigen::MatrixXd _trafomatrix;
      if (package != "gaussian") {
        std::cout << " I should not have been called, will do nothing! " << std::endl;
      } else {
        // go through basisset, determine function sizes
        int _dim_sph = 0;
        int _dim_cart = 0;
        for (const AOShell* shell:(*this)) {
          const std::string& _type = shell->getType();

          _dim_sph += NumFuncShell(_type);
          _dim_cart += NumFuncShell_cartesian(_type);

        }
        _trafomatrix = Eigen::MatrixXd::Zero(_dim_sph, _dim_cart);

        int _row_start = 0;
        int _col_start = 0;
         for (const AOShell* shell:(*this)) {
          const std::string& _type = shell->getType();
          int _row_end = _row_start + NumFuncShell(_type);
          int _col_end = _col_start + NumFuncShell_cartesian(_type);
          Eigen::Block<Eigen::MatrixXd> block = _trafomatrix.block(_row_start, _col_start, NumFuncShell(_type), NumFuncShell_cartesian(_type));
          addTrafoCartShell(shell, block);
          _row_start = _row_end;
          _col_start = _col_end;

        }
      }
      return _trafomatrix;
    }


void AOBasis::addTrafoCartShell( const AOShell* shell , Eigen::Block<Eigen::MatrixXd>& submatrix ){

    // fill _local according to _lmax;
    int lmax = shell->getLmax();
    std::string type = shell->getType();

    int sph_size =NumFuncShell( type ) + OffsetFuncShell( type );
    int cart_size = NumFuncShell_cartesian( type ) + OffsetFuncShell_cartesian( type )  ;
    Eigen::MatrixXd local = Eigen::MatrixXd::Zero(sph_size,cart_size);

    // s-functions
    local(0,0) = 1.0; // s
    // p-functions
    if ( lmax > 0 ){
        local(1,1) = 1.0;
        local(2,2) = 1.0;
        local(3,3) = 1.0;
    }
    // d-functions
    if ( lmax > 1 ){
        local(4,4) = -0.5;             // d3z2-r2 (dxx)
        local(4,5) = -0.5;             // d3z2-r2 (dyy)
        local(4,6) =  1.0;             // d3z2-r2 (dzz)
        local(5,8) =  1.0;             // dxz
        local(6,9) =  1.0;             // dyz
        local(7,4) = 0.5*sqrt(3.0);    // dx2-y2 (dxx)
        local(7,5) = -local(7,4);      // dx2-y2 (dyy)
        local(8,7) = 1.0;              // dxy
     }
    if ( lmax > 2 ){
        std::cerr << " Gaussian input with f- functions or higher not yet supported!" << std::endl;
        exit(1);
    }
    // now copy to _trafo
    for ( int i_sph = 0 ; i_sph < NumFuncShell( type ) ; i_sph++ ){
        for  ( int i_cart = 0 ; i_cart < NumFuncShell_cartesian( type ) ; i_cart++ ){
            submatrix( i_sph , i_cart ) = local( i_sph + OffsetFuncShell( type ) , i_cart +  OffsetFuncShell_cartesian( type ) );
        }
    }
    return;
}

std::vector<int> AOBasis::getMultiplierVector( const std::string& start, const std::string& target){
    std::vector<int> multiplier;
    multiplier.reserve(_AOBasisSize);
    std::string s;
    std::string t;
    if(start=="xtp"){
      s=target;
      t=start;
    }else{
      s=start;
      t=target;
    }
    // go through basisset
    for (const AOShell* shell:(*this)) { 
        addMultiplierShell(  s, t, shell->getType(), multiplier );
    }
    return multiplier;
    }

void AOBasis::addMultiplierShell(const std::string& start, const std::string& target, const std::string& shell_type, std::vector<int>& multiplier) {
//multipliers were all found using code, hard to establish

    if (target == "xtp") {
        // current length of vector
        //int _cur_pos = multiplier.size() - 1;

        // single type shells defined here
        if (shell_type.length() == 1) {
            if (shell_type == "S") {
                multiplier.push_back(1);
            }
            else if (shell_type == "P") {
                multiplier.push_back(1);
                multiplier.push_back(1);
                multiplier.push_back(1);
            }
            else if (shell_type == "D") {
                if (start == "nwchem") {
                    multiplier.push_back(1);
                    multiplier.push_back(1);
                    multiplier.push_back(-1);
                    multiplier.push_back(1);
                    multiplier.push_back(1);
                } else if (start == "orca"){
                    multiplier.push_back(1);
                    multiplier.push_back(1);
                    multiplier.push_back(1);
                    multiplier.push_back(1);
                    multiplier.push_back(1);
                } else {
                    std::cerr << "Tried to get multipliers d-functions from package " << start << ".";
                    throw std::runtime_error("Multiplication not implemented yet!");
                }
            }
            else if (shell_type == "F") {
                if ( start == "orca" ){
                    multiplier.push_back(1);
                    multiplier.push_back(1);
                    multiplier.push_back(1);
                    multiplier.push_back(1);
                    multiplier.push_back(1);
                    multiplier.push_back(-1);
                    multiplier.push_back(-1);

                }else if (start == "nwchem"){
                    multiplier.push_back(1);
                    multiplier.push_back(1);
                    multiplier.push_back(-1);
                    multiplier.push_back(+1);
                    multiplier.push_back(+1);
                    multiplier.push_back(+1);
                    multiplier.push_back(-1);
                } else {
                std::cerr << "Tried to get multipliers f-functions from package " << start << ".";
                throw std::runtime_error("Multiplication not implemented yet!");
                }
            }
            else if (shell_type == "G") {
                if ( start == "orca" ){
                    //Not checked yet
                    multiplier.push_back(1);
                    multiplier.push_back(1);
                    multiplier.push_back(1);
                    multiplier.push_back(1);
                    multiplier.push_back(1);
                    multiplier.push_back(-1);
                    multiplier.push_back(-1);
                    multiplier.push_back(-1);
                    multiplier.push_back(-1);

                } else {
                std::cerr << "Tried to get multipliers g-functions from package " << start << ".";
                throw std::runtime_error("Multiplication not implemented yet!");
                }
            }else{
                std::cerr << "Tried to get multipliers h-functions . ";
                throw std::runtime_error("Multiplication not implemented yet!");
            }
        } else {
            // for combined shells, iterate over all contributions
            //_nbf = 0;
            for (unsigned i = 0; i < shell_type.length(); ++i) {
                std::string local_shell = std::string(shell_type, i, 1);
                addMultiplierShell(start, target, local_shell, multiplier);
            }
        }
    } else {

        std::cerr << "Tried to reorder functions (multiplier)  from " << start << " to " << target << std::endl;
        throw std::runtime_error("Reordering not implemented yet!");


    }
    return;
}


std::vector<int>  AOBasis::getReorderVector(const std::string& start,const std::string& target){
    std::vector<int> neworder;
    neworder.reserve(_AOBasisSize);
    std::string s;
    std::string t;
    if(start=="xtp"){
      s=target;
      t=start;
    }else{
      s=start;
      t=target;
    }
    // go through basisset
     for (const AOShell* shell:(*this)) {
        addReorderShell( s, t, shell->getType(), neworder );
    }
     if(start=="xtp"){
       neworder=invertOrder(neworder);
     } 
    return neworder;
}

std::vector<int> AOBasis::invertOrder(const std::vector<int>& order ){
 
    std::vector<int>neworder=std::vector<int>(order.size());
    for(unsigned i=0;i<order.size();i++){
        neworder[order[i]]=int(i);
    }
    return neworder;
    }

    void AOBasis::addReorderShell(const std::string& start, const std::string& target, const std::string& shell_type, std::vector<int>& order) {
        //Reordering is given by email from gaussian, orca output MOs, and http://www.nwchem-sw.org/index.php/Release66:Basis for nwchem
        
      // current length of vector

      int cur_pos = order.size() - 1;

      if (target == "xtp") {
        // single type shells defined here
        if (shell_type.length() == 1) {
          if (shell_type == "S") {
            order.push_back(cur_pos + 1);
          }//for S

            //votca order is z,y,x e.g. Y1,0 Y1,-1 Y1,1
          else if (shell_type == "P") {
            if (start == "orca") {
                //orca order is z,x,y Y1,0,Y1,1,Y1,-1
              order.push_back(cur_pos + 1);
              order.push_back(cur_pos + 3);
              order.push_back(cur_pos + 2);
            } else if (start == "gaussian" || start == "nwchem") {
                //nwchem gaussian x,y,z Y1,1 Y1,-1 Y1,0
              order.push_back(cur_pos + 3);
              order.push_back(cur_pos + 2);
              order.push_back(cur_pos + 1);
            } else if (start == "votca") {//for usage with old orb files
                 //old votca x,y,z Y1,1 Y1,-1 Y1,0
              order.push_back(cur_pos + 3);
              order.push_back(cur_pos + 2);
              order.push_back(cur_pos + 1);
            } else {
              std::cerr << "Tried to reorder p-functions from package " << start << ".";
              throw std::runtime_error("Reordering not implemented yet!");
            }
          }//for P
            //votca order is d3z2-r2 dyz dxz dxy dx2-y2 e.g. Y2,0 Y2,-1 Y2,1 Y2,-2 Y2,2
          else if (shell_type == "D") {
            //orca order is d3z2-r2 dxz dyz dx2-y2 dxy e.g. Y2,0 Y2,1 Y2,-1 Y2,2 Y2,-2
            if (start == "gaussian" || start == "orca") {
              order.push_back(cur_pos + 1);
              order.push_back(cur_pos + 3);
              order.push_back(cur_pos + 2);
              order.push_back(cur_pos + 5);
              order.push_back(cur_pos + 4);
            } else if (start == "nwchem") {
              // nwchem order is dxy dyz d3z2-r2 -dxz dx2-y2, e.g. Y2,-2 Y2,-1 Y2,0 Y2,1 Y2,2 
              order.push_back(cur_pos + 4);
              order.push_back(cur_pos + 2);
              order.push_back(cur_pos + 1);
              order.push_back(cur_pos + 3);
              order.push_back(cur_pos + 5);
            } else if (start == "votca") { //for usage with old orb files
              order.push_back(cur_pos + 3);
              order.push_back(cur_pos + 2);
              order.push_back(cur_pos + 4);
              order.push_back(cur_pos + 1);
              order.push_back(cur_pos + 5);
            } else {
              std::cerr << "Tried to reorder d-functions from package " << start << ".";
              throw std::runtime_error("Reordering not implemented yet!");
            }
          } else if (shell_type == "F") {
               //ordering for votca is Yl,0 Yl,-1 Yl,1 ......Yl,-m Yl,m
            if (start == "gaussian" || start == "orca") {
                //ordering for gaussian and orca is Yl,0 Yl,1 Yl,-1 ......Yl,m Yl,-m
              order.push_back(cur_pos + 1);
              order.push_back(cur_pos + 3);
              order.push_back(cur_pos + 2);
              order.push_back(cur_pos + 5);
              order.push_back(cur_pos + 4);
              order.push_back(cur_pos + 7);
              order.push_back(cur_pos + 6);
            } else if (start == "nwchem") {
                //ordering for nwchem is fxxy-yyy, fxyz,fyzz-xxy-yyy,fzzz-xxz-yyz,f-xzz+xxx+xyy,fxxz-yyz,fxyy-xxx
                // e.g. Y3,-3 Y3,-2 Y3,-1 Y3,0 Y3,1 Y3,2 Y3,3
              order.push_back(cur_pos + 6);
              order.push_back(cur_pos + 4);
              order.push_back(cur_pos + 2);
              order.push_back(cur_pos + 1);
              order.push_back(cur_pos + 3);
              order.push_back(cur_pos + 5);
              order.push_back(cur_pos + 7);
            } else {
              std::cerr << "Tried to reorder f-functions from package " << start << ".";
              throw std::runtime_error("Reordering not implemented yet!");
            }
          }else if (shell_type == "G") {
               //ordering for votca is Yl,0 Yl,-1 Yl,1 ......Yl,-m Yl,m
            if (start == "gaussian" || start == "orca") {
                 //ordering for gaussian and orca is Yl,0 Yl,1 Yl,-1 ......Yl,m Yl,-m
              order.push_back(cur_pos + 1);
              order.push_back(cur_pos + 3);
              order.push_back(cur_pos + 2);
              order.push_back(cur_pos + 5);
              order.push_back(cur_pos + 4);
              order.push_back(cur_pos + 7);
              order.push_back(cur_pos + 6);
              order.push_back(cur_pos + 9);
              order.push_back(cur_pos + 8);
            }else  {
              std::cerr << "Tried to reorder g-functions from package " << start << ".";
              throw std::runtime_error("Reordering not implemented");
            }
          }else {
            std::cerr << "Tried to reorder functions  of shell type " << shell_type << std::endl;
            throw std::runtime_error("Reordering not implemented");
          }
        } else {
          // for combined shells, iterate over all contributions
          //_nbf = 0;
          for (unsigned i = 0; i < shell_type.length(); ++i) {
            std::string local_shell = std::string(shell_type, i, 1);
            this->addReorderShell(start, target, local_shell, order);
          }
        }

      } else {
        std::cerr << "Tried to reorder functions (neworder) from " << start << " to " << target << std::endl;
        throw std::runtime_error("Reordering not implemented yet!");
      }
      return;
    }

    const std::vector<const AOShell*> AOBasis::getShellsperAtom(int AtomId)const {
      std::vector<const AOShell*> result;
      for (const auto& aoshell : _aoshells) {
        if (aoshell->getAtomIndex() == AtomId) {
          result.push_back(aoshell);
        }
      }
      return result;
    }

    int AOBasis::getFuncperAtom(int AtomId) const {
      int number = 0;
      for (const auto& aoshell : _aoshells) {
        if (aoshell->getAtomIndex() == AtomId) {
          number += aoshell->_numFunc;
        }
      }
      return number;
    }

    void AOBasis::AOBasisFill(const BasisSet& bs, std::vector<QMAtom* >& atoms, int fragbreak) {
      tools::Elements elementinfo;
      std::vector<QMAtom* > ::iterator ait;
      _AOBasisSize = 0;
      _AOBasisFragA = 0;
      _AOBasisFragB = 0;
      // loop over atoms
      for (QMAtom* atom : atoms) {
        const std::string& name = atom->getType();
        atom->_nuccharge = elementinfo.getNucCrg(name);
        const Element& element = bs.getElement(name);
        for (const Shell& shell:element) {
          int numfuncshell = NumFuncShell(shell.getType());
          AOShell* aoshell = addShell(shell, *atom, _AOBasisSize);
          _AOBasisSize += numfuncshell;
          for (const GaussianPrimitive& gaussian:shell) {
            aoshell->addGaussian(gaussian);
          }
          aoshell->CalcMinDecay();
          aoshell->normalizeContraction();
        }
        if (atom->getAtomID() < fragbreak) _AOBasisFragA = _AOBasisSize;
      }

      if (fragbreak < 0) {
        _AOBasisFragA = _AOBasisSize;
        _AOBasisFragB = 0;
      } else {
        _AOBasisFragB = _AOBasisSize - _AOBasisFragA;
      }
      return;
    }

    void AOBasis::ECPFill(const BasisSet& bs, std::vector<QMAtom* >& atoms) {

      _AOBasisSize = 0;
      for (QMAtom* atom : atoms) {
        std::string name = atom->getType();
        if (name == "H" || name == "He") {
          continue;
        }
        const Element& element = bs.getElement(name);
        atom->_ecpcharge = element.getNcore();
        int lmax = element.getLmax();
        for (const Shell& shell:element) {
          if (shell.getType().size() > 1) {
            throw std::runtime_error("In ecps no combined shells e.g. SP are allowed");
          }
          //Local part is with L=Lmax
          bool nonlocal = false;
          if (shell.getLmax() < lmax) {
            nonlocal = true;
          }

          AOShell* aoshell = addECPShell(shell, *atom, _AOBasisSize, nonlocal);
          _AOBasisSize += NumFuncShell(shell.getType());
         for (const GaussianPrimitive& gaussian:shell) {
            aoshell->addGaussian(gaussian);
          }
          aoshell->CalcMinDecay();
        }
      }
      return;
    }











}}
