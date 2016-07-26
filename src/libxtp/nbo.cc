/*
 *            Copyright 2009-2016 The VOTCA Development Team
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


#include <votca/xtp/nbo.h>
#include <votca/xtp/aomatrix.h>
namespace votca { namespace xtp {

void NBO::EvaluateNBO(std::vector< QMAtom* >& _atomlist,  ub::matrix<double> &_dmat,AOBasis &basis,BasisSet &bs){
    AOOverlap _overlap;
    // initialize overlap matrix
    _overlap.Initialize(basis._AOBasisSize);
    // Fill overlap
    _overlap.Fill(&basis);
    
    ub::matrix<double> _prodmat = ub::prod( _dmat, _overlap._aomatrix );
    

    
    ub::matrix<double> P=ub::prod(_overlap._aomatrix,_prodmat);
   
    IntercenterOrthogonalisation(P,_overlap._aomatrix, _atomlist ,bs);
    
    
   
         //cout << id << " "<< id+nooffunc << endl;
      
    }
    //cout << id << " " << _dmat.size1() << endl;





void NBO::IntercenterOrthogonalisation(ub::matrix<double> &P,ub::matrix<double> &Overlap,vector< QMAtom* >& _atomlist ,BasisSet &bs){
    
    
    vector < QMAtom* > :: iterator atom;

// make an array which stores for each atom the the starting point for all s functions, p functions, d functions etc...   
    
    
    //atom/shell/starting 
    std::vector< std::map< string,std::vector< std::pair<int,int> > > >sorting;
    
    int functionindex=0;
    for (atom = _atomlist.begin(); atom < _atomlist.end(); ++atom){
        std::map< string,std::vector< std::pair<int,int> > > shellsort;
        Element* element = bs.getElement((*atom)->type);
        for (Element::ShellIterator its = element->firstShell(); its != element->lastShell(); its++) {
            
           // for loop because shells can also consist of SP shells or alike
            for(unsigned i = 0; i <(*its)->getType().length(); ++i) {
            string local_shell = string( (*its)->getType(), i, 1 );
            std::vector< std::pair<int,int> >& index =shellsort[local_shell];
            int end=functionindex+(*its)->FindnumofFunc(local_shell);
            index.push_back(std::pair<int,int>(functionindex,end));
            
            functionindex=end;
            }
        }
        sorting.push_back(shellsort);
         }
    
    
    /*
    for (unsigned i=0;i<sorting.size();i++){
        cout << "Atom "<<i<<endl;
        typedef std::map< string,std::vector< std::pair<int,int> > >::iterator it_type;
        for(it_type iterator = sorting[i].begin(); iterator != sorting[i].end(); iterator++) {
    cout<< "Shell "<<iterator->first<<":";
    std::vector< std::pair<int,int> >& index=iterator->second;
    for (unsigned j=0;j<index.size();j++){
        cout <<" ["<<index[j].first<<":"<<index[j].second<<"]";
    }
    cout<<endl;
}
    }
    */
    
    
   
    
    
    
    return;
}


void NBO::TransformMatrixtoNewBasis(ub::matrix<double>& Matrix,const ub::matrix<double>& transformation){
    ub::matrix<double> temp=ub::prod(Matrix,ub::trans(transformation));
    Matrix=ub::prod(transformation,temp);
    return;
}



void NBO::LoadMatrices(std::string fn_projectionMatrix, std::string fn_overlapMatrix){

    //TODO: Yuriy, fill this in
    
    return;
    }

}}