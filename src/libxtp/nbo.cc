/*
 *            Copyright 2009-2017 The VOTCA Development Team
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
#include <votca/tools/linalg.h>
namespace votca { namespace xtp {

void NBO::EvaluateNBO(std::vector< ctp::QMAtom* >& _atomlist,const  ub::matrix<double> &_dmat,const AOBasis &basis, BasisSet &bs){
    AOOverlap _overlap;
    // initialize overlap matrix
    _overlap.Initialize(basis.AOBasisSize());
    // Fill overlap
    _overlap.Fill(basis);
    
    ub::matrix<double> _prodmat = ub::prod( _dmat, _overlap.Matrix() );
    
    
    
    ub::matrix<double> P=ub::prod(_overlap.Matrix(),_prodmat);
   
    ub::matrix<double> PNAOs_trans=IntercenterOrthogonalisation(P,_overlap.Matrix(), _atomlist ,bs);
    
    cout<<P<<endl;
    cout<<_overlap.Matrix()<<endl;
    
    vector < ctp::QMAtom* > :: iterator atom;
    for (atom = _atomlist.begin(); atom < _atomlist.end(); ++atom){
    
    std::vector<int> func=_elements.getMinimalBasis((*atom)->type,_ECP);
  
    
    
   
         //cout << id << " "<< id+nooffunc << endl;
      
    }
    //cout << id << " " << _dmat.size1() << endl;

    return;
}


ub::matrix<double>NBO::IntercenterOrthogonalisation(ub::matrix<double> &P, ub::matrix<double> &Overlap,vector< ctp::QMAtom* >& _atomlist , BasisSet &bs){
    
    
    vector< ctp::QMAtom* >::iterator atom;

// make an array which stores for each atom the the starting point for all s functions, p functions, d functions etc...   
    
    
    //atom/shell/starting 
    std::vector< std::map< int,std::vector< int > > >sorting;
    
    int functionindex=0;
    for (atom = _atomlist.begin(); atom < _atomlist.end(); ++atom){
        std::map< int,std::vector< int > > shellsort;
        Element* element = bs.getElement((*atom)->type);
        for (Element::ShellIterator its = element->firstShell(); its != element->lastShell(); its++) {
            
           // for loop because shells can also consist of SP shells or alike
            for(unsigned i = 0; i <(*its)->getType().length(); ++i) {
            string local_shell = string( (*its)->getType(), i, 1 );
            int l=FindLmax(local_shell);
            std::vector< int >& index =shellsort[l];
           
            index.push_back(functionindex);
            
            functionindex+=NumFuncShell(local_shell);
            }
        }
        sorting.push_back(shellsort);
         }
    typedef std::map< int,std::vector< int > >::iterator it_type;
    
    
    for (unsigned i=0;i<sorting.size();i++){
        cout << "Atom "<<i<<endl;
        
        for(it_type iterator = sorting[i].begin(); iterator != sorting[i].end(); iterator++) {
    cout<< "Shell "<<iterator->first<<":";
    std::vector< int >& index=iterator->second;
    for (unsigned j=0;j<index.size();j++){
        cout <<" "<<index[j];
    }
    cout<<endl;
}
    }
    
    
    ub::matrix<double> transformation=ub::zero_matrix<double>(P.size1(),P.size2());
    ub::vector<double> occupancies=ub::zero_vector<double>(P.size1());
    // Stting up atomic transformations
    
    //atoms
    for (unsigned i=0;i<sorting.size();i++){
        //shell
        for(it_type iterator = sorting[i].begin(); iterator != sorting[i].end(); iterator++) {
            int ll=2*(iterator->first)+1;
            std::vector< int >& index=iterator->second;
            unsigned size=index.size();
            ub::matrix<double> P_L=ub::zero_matrix<double>(size,size);
            ub::matrix<double> S_L=ub::zero_matrix<double>(size,size);
            
            //filling P_L and S_L
            for (unsigned i=0;i<size;i++){
                for (unsigned j=0;j<size;j++){
                    //summing over m
                    for (int m=0;m<ll;m++){
                        P_L(i,j)+=P(index[i]+m,index[j]+m);
                        S_L(i,j)+=Overlap(index[i]+m,index[j]+m);
                    } 
                }
            }
            cout << P_L<< endl;
            cout << S_L<<endl;
            
            //Setup generalized eigenproblem for each P_L/S_L
            ub::vector<double> eigenvalues;
            ub::matrix<double> eigenvectors;
            bool check=linalg_eigenvalues_general( P_L,S_L, eigenvalues, eigenvectors);
            if (!check){
                throw runtime_error("Diagonalisation failed");
            }
            cout << eigenvalues<<endl;
            cout <<eigenvectors<<endl;
            
            
            for (unsigned i=0;i<size;i++){
                for (int m=0;m<ll;m++){
                occupancies(index[i]+m)=eigenvalues(i);
                }
                for (unsigned j=0;j<size;j++){
                    for (int m=0;m<ll;m++){
                    transformation(index[i]+m,index[j]+m)=eigenvectors(i,j);
                }
            }
        }
        }
    }
        cout<<transformation<<endl;
        cout<<occupancies<<endl;
    
        TransformMatrixtoNewBasis(P,transformation);
        TransformMatrixtoNewBasis(Overlap,transformation);
    
    return transformation;
}


void NBO::TransformMatrixtoNewBasis(ub::matrix<double>& Matrix, const ub::matrix<double>& transformation){
    ub::matrix<double> temp=ub::prod(Matrix,ub::trans(transformation));
    Matrix=ub::prod(transformation,temp);
    return;
}



void NBO::LoadMatrices(std::string fn_projectionMatrix, std::string fn_overlapMatrix){

    //TODO: Yuriy, fill this in
    
    return;
    }

}}