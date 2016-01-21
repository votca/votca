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
 *Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */


#include <votca/xtp/ERIs.h>


using namespace std;
using namespace votca::tools;

namespace votca {
    namespace xtp {
        namespace ub = boost::numeric::ublas;
        
        
        
        void ERIs::Initialize (AOBasis &_dftbasis, AOBasis &_auxbasis,  AOOverlap &_auxAOverlap, AOCoulomb &_auxAOcoulomb){

           
           
            _threecenter.Fill( _auxbasis, _dftbasis );
            /* cout << "\n" << endl;
             for (int k=0; k<_auxbasis._AOBasisSize; k++){
                cout << " s-s " << _threecenter.getDatamatrix(k)(0,0) << "  s-px  " << _threecenter.getDatamatrix(k)(0,4) << " sym " << _threecenter.getDatamatrix(k)(4,0) << endl;
                       
             }
            exit(0); */
            
            //ub::matrix<double> _inverse=ub::zero_matrix<double>( _auxAOverlap.Dimension(), _auxAOverlap.Dimension());
            
            /* AOOverlap _auxoverlap_inverse;               
            AOOverlap _auxoverlap_cholesky_inverse;      
            _auxoverlap_inverse.Initialize( _auxbasis._AOBasisSize);
            _auxAOcoulomb.Symmetrize(_auxAOverlap , _auxbasis, _auxoverlap_inverse , _auxoverlap_cholesky_inverse);
            */
            //linalg_invert( _auxAOverlap.Matrix() , _inverse);
            //ub::matrix<double> _temp=ub::prod(_auxAOcoulomb.Matrix(),_inverse);
            _Vcoulomb=_auxAOcoulomb.Matrix();
            
            
            //cout << "Vcoulomb"<< _Vcoulomb<< endl;
 /*
            int size4c=_dftbasis.AOBasisSize();

            // reconstructing 4center-integrals
	  
            for (int k=0; k<_auxbasis._AOBasisSize; k++){
                
                ub::matrix<double> shit = _threecenter.getDatamatrix(k);
                
               for (int alpha=0; alpha<size4c; alpha++){
                    for (int beta=0; beta<size4c ;beta++){
                        
                        
                        
                        cout << "SHIT " << k << " " <<  alpha << " " << beta << " " << shit(alpha,beta) << " vs " << shit(beta,alpha) << endl;
                        
                        
                    }
               }
                
            }
            
             for (int k=0; k<_auxbasis._AOBasisSize; k++){
				for (int l=0; l<_auxbasis._AOBasisSize; l++){
                                 
                                    
                                    cout << " VCOU " << k << " " << l << " " << scientific << _Vcoulomb(k,l) << " vs " << _Vcoulomb(l,k) << endl; 
                                    
                                }
                                }
            
            
            
            
            //exit(0);
            
            
        
            

           
            
            typedef boost::multi_array<double, 4> fourdim;
            fourdim  fourcenter(boost::extents[size4c][size4c][size4c][size4c]);
            cout <<endl;
       
                          for (int alpha=0; alpha<size4c; alpha++){
                    for (int beta=0; beta<size4c ;beta++){
           
                       for (int mu=0; mu<size4c ;mu++){
                            for (int nu=0; nu<size4c; nu++){
			      fourcenter[alpha][beta][mu][nu]=0.0;
			      for (int k=0; k<_auxbasis._AOBasisSize; k++){
				for (int l=0; l<_auxbasis._AOBasisSize; l++){
                                                 //cout<<_threecenter.getDatamatrix(k)(alpha,beta)<<"  "<<_threecenter.getDatamatrix(l)(mu,nu)<<"  "<<_Vcoulomb(k,l)<<endl;
            fourcenter[alpha][beta][mu][nu]+=_Vcoulomb(k,l)*_threecenter.getDatamatrix(k)(alpha,beta)*_threecenter.getDatamatrix(l)(mu,nu);
             //cout<<   fourcenter[alpha][beta][mu][nu]<< endl;                            
                                    }}
                                    

            cout << "4c("<<alpha+1<<":"<<beta+1<<":"<<mu+1<<":"<<nu+1<<") = "<< fourcenter[alpha][beta][mu][nu]<< endl;
            }}
            }}

exit(0);
	    
            
            
            
          
            AOOverlap _auxoverlap_inverse;               // will also be needed in PPM itself
            AOOverlap _auxoverlap_cholesky_inverse;      // will also be needed in PPM itself
            _auxoverlap_inverse.Initialize( _auxbasis._AOBasisSize);
            _auxAOcoulomb.Symmetrize(_auxAOverlap , _auxbasis, _auxoverlap_inverse , _auxoverlap_cholesky_inverse);

            std::vector< ub::matrix<double> > I_times_sqrtV;
            for (int i=0; i< _auxbasis._AOBasisSize; i++){
                 I_times_sqrtV.push_back(ub::zero_matrix<double>(_dftbasis._AOBasisSize, _dftbasis._AOBasisSize));        
            }
            
            cout <<"Coulomb symm"<< _auxAOcoulomb.Matrix() <<endl;
            
            for (int m=0;m<_auxbasis._AOBasisSize; m++){
                for (int l=0;l<_auxbasis._AOBasisSize; l++){
                    I_times_sqrtV[m]+=_auxAOcoulomb.Matrix()(m,l)*_threecenter.getDatamatrix(l);
                }
            }
            
            
            int size4c=_dftbasis.AOBasisSize();
            typedef boost::multi_array<double, 4> fourdim;
            fourdim  fourcenter(boost::extents[size4c][size4c][size4c][size4c]);
            
            for (int alpha=0;alpha<size4c;alpha++){
                for (int beta=0;beta<size4c;beta++){
                   for (int mu=0;mu<size4c;mu++){
                        for (int nu=0;nu<size4c;nu++){
                    fourcenter[alpha][beta][mu][nu]=0.0;
                    for (int m=0;m<_auxbasis._AOBasisSize;m++){
                        fourcenter[alpha][beta][mu][nu]+= I_times_sqrtV[m](alpha,beta)*I_times_sqrtV[m](mu,nu);
                    }
                    cout << "4c("<<alpha+1<<":"<<beta+1<<":"<<mu+1<<":"<<nu+1<<")="<< fourcenter[alpha][beta][mu][nu]<< endl;
                    
                        }}}exit(0);}
                    
          
            */
            
            
            
         
        }
        
        

        
        
        
        
        
        void ERIs::Initialize_Symmetric (AOBasis &_dftbasis, AOBasis &_auxbasis,  AOOverlap &_auxAOverlap, AOCoulomb &_auxAOcoulomb){

           
           
            _threecenter.Fill( _auxbasis, _dftbasis );
            
            
            ub::matrix<double> _inverse=ub::zero_matrix<double>( _auxAOverlap.Dimension(), _auxAOverlap.Dimension());
            
            AOOverlap _auxoverlap_inverse;               
            AOOverlap _auxoverlap_cholesky_inverse;      
            _auxoverlap_inverse.Initialize( _auxbasis._AOBasisSize);
            _auxAOcoulomb.Symmetrize_DFT(_auxAOverlap , _auxbasis, _auxoverlap_inverse , _auxoverlap_cholesky_inverse);
            
           // linalg_invert( _auxAOverlap.Matrix() , _inverse);
           // ub::matrix<double> _temp=ub::prod(_auxAOcoulomb.Matrix(),_inverse);
           // _Vcoulomb=ub::prod(_inverse,_temp);
            
            //cout << "Vcoulomb"<< _Vcoulomb<< endl;


          // multiply each element of _threecenters with symmetrized Coulomb
          // K(l';a,b) = sum_l{ sqrt{V}(l',l) I(l;a,b) }
            
            std::vector< ub::matrix<double> > K ;
            for (int i=0; i< _auxbasis._AOBasisSize; i++){
                 K.push_back(ub::zero_matrix<double>(_dftbasis._AOBasisSize, _dftbasis._AOBasisSize));        
            }

            for (int m=0;m<_auxbasis._AOBasisSize; m++){
                for (int l=0;l<_auxbasis._AOBasisSize; l++){
                    K[m]+=_auxAOcoulomb.Matrix()(m,l)*_threecenter.getDatamatrix(l);
                }
            }
            
            // update _threecenters
            
            for (int m=0;m<_auxbasis._AOBasisSize; m++){
                    _threecenter.getDatamatrix(m)=K[m];
            }
            
            
          
            int size4c=_dftbasis.AOBasisSize();
            typedef boost::multi_array<double, 4> fourdim;
            fourdim  fourcenter(boost::extents[size4c][size4c][size4c][size4c]);
            
            for (int alpha=0; alpha<size4c; alpha++){
                for (int beta=0; beta<size4c; beta++){
                   for (int mu=0; mu<size4c; mu++){
                        for (int nu=0; nu<size4c; nu++){
                    fourcenter[alpha][beta][mu][nu]=0.0;
                    for (int m=0;m<_auxbasis._AOBasisSize;m++){
                        fourcenter[alpha][beta][mu][nu]+= K[m](alpha,beta)*K[m](mu,nu);
                    }
                    cout << "4c("<<alpha+1<<":"<<beta+1<<":"<<mu+1<<":"<<nu+1<<")="<< fourcenter[alpha][beta][mu][nu]<< endl;
                    
                        }}}}
            
            exit(0);
                    
          
            
            
            
            
         

        
        
        }
        
        
        void ERIs::CalculateERIs (ub::matrix<double> &DMAT, AOOverlap &_auxAOverlap, ub::matrix<double>& _AOIntegrals ){
            
            
            
            /***** TEST VIA RECONSTRUCTED 4center ERIs ***************/
            
/*
            int size4c=DMAT.size1(); // _dftbasis.AOBasisSize();
            int auxsize = _Vcoulomb.size1();
            
            typedef boost::multi_array<double, 4> fourdim;
            fourdim  fourcenter(boost::extents[size4c][size4c][size4c][size4c]);
            cout <<endl;
       
           
            for (int alpha=0;alpha<size4c;alpha++){
                    for (int beta=0;beta<size4c;beta++){
                       for (int mu=0;mu<size4c;mu++){
                            for (int nu=0;nu<size4c;nu++){
                        fourcenter[alpha][beta][mu][nu]=0.0;

                        
                        for (int k=0;k<auxsize;k++){
                                    for (int l=0;l<auxsize;l++){
                                                 //cout<<_threecenter.getDatamatrix(k)(alpha,beta)<<"  "<<_threecenter.getDatamatrix(l)(mu,nu)<<"  "<<_Vcoulomb(k,l)<<endl;
            fourcenter[alpha][beta][mu][nu]+=_Vcoulomb(k,l)*_threecenter.getDatamatrix(k)(alpha,beta)*_threecenter.getDatamatrix(l)(mu,nu);
             //cout<<   fourcenter[alpha][beta][mu][nu]<< endl;                            
                                    }}
                                    

            //cout << "4c("<<alpha+1<<":"<<beta+1<<":"<<mu+1<<":"<<nu+1<<")="<< fourcenter[alpha][beta][mu][nu]<< endl;
            }}
            }}
            
            

            */
            
            
      /*      
            // TEST Number of fitted electrons
            ub::matrix<double> _FITMAT = ub::zero_matrix<double>(DMAT.size1(),DMAT.size2());
            ub::matrix<double> _inverse=ub::zero_matrix<double>( _auxAOverlap.Dimension(), _auxAOverlap.Dimension());
            
           //_auxAOverlap.Print("AUX-OL");
            linalg_invert( _auxAOverlap._aomatrix , _inverse);
            
           ub:: matrix<double> _testS = ub::prod(_auxAOverlap._aomatrix , _inverse);
                      for (int _i = 0; _i < _threecenter.getSize(); _i++) {
                for (int _j = 0; _j < _threecenter.getSize(); _j++) {
                    cout << " S*S-1(" << _i << "," << _j << ") " << setprecision (9) << _testS(_i,_j) << endl; 
                      }
        }
           
          */  
            //cout <<  "3c size " << _threecenter.getSize() << " _auxAODim " << _auxAOverlap.Dimension() << " DMAT1 " << DMAT.size1() << " DMAT2 " << DMAT.size2() << endl;
            
            
            //cout << _inverse << endl;
       /*     
            cout << " 3C size " << _threecenter.getSize() << endl;
            
            for (int _i = 0; _i < _threecenter.getSize(); _i++) {
                for (int _j = 0; _j < _threecenter.getSize(); _j++) {

                //cout << " S-1(" << _i << "," << _j << ") " << setprecision (9) << _inverse(_i,_j) << endl; 
                
                    ub::matrix<double>& _3c = _threecenter.getDatamatrix(_j);

                    for ( int alpha = 0 ; alpha < DMAT.size1(); alpha++ ){
                        for ( int beta = 0 ; beta < DMAT.size2(); beta++ ){
                        
                            _FITMAT(alpha,beta) += _inverse(_i,_j)*_3c(alpha,beta)*_AOIntegrals(0,_i);
                            cout << " inV " << _inverse(_i,_j) << " 3C " << _3c(alpha,beta) << endl;
                            
                        }
                    }
                    
                    
                }
            }
            
           
  
            
            
            
            
//            ub::vector<double> Fitarray=_FITMAT.data();

            double _Nfit = 0.0;
            
            for ( int alpha = 0 ; alpha < DMAT.size1(); alpha++ ){
                        for ( int beta = 0 ; beta < DMAT.size1(); beta++ ){
                    
                        cout << alpha << "  " << beta << "  " << _FITMAT(alpha,beta) << endl;  
                        
                        _Nfit += DMAT(alpha,beta) * _FITMAT(alpha,beta);
                }
                }
            

            
            
            cout << " \nNumber of electrons in Density Fit  " << _Nfit << endl;
            
            
            //exit(0);
            exit(0);
            
            
            _ERIs=ub::zero_matrix<double>(DMAT.size1(),DMAT.size2());
            ub::vector<double> dmatasarray=DMAT.data();
            
            ub::matrix<double> Itilde=ub::zero_matrix<double>(_threecenter.getSize(),1);
            //cout << _threecenter.getSize() << " Size-Threecenter"<<endl;
            //check Efficiency !!!! someday 
            for ( int _i=0; _i<_threecenter.getSize();_i++){
                ub::vector<double>threecenterasarray=(_threecenter.getDatamatrix(_i)).data();
                //cout << _threecenter.getDatamatrix(_i).size1() << "x"<< _threecenter.getDatamatrix(_i).size2() <<" Size-Threecenter,matrix"<<endl;
                // Trace over prod::DMAT,I(l)=componentwise product over 
                for ( int _j=0; _j<threecenterasarray.size();_j++){
                    Itilde(_i,0)+=dmatasarray[_j]*threecenterasarray[_j];
                }
            }
            //cout << "Itilde " <<Itilde << endl;
            ub::matrix<double>K=ub::prod(_Vcoulomb,Itilde);
            //cout << "K " << K << endl;
            for ( int _i=0; _i<K.size1(); _i++){
                
            _ERIs+=_threecenter.getDatamatrix(_i)*K(_i,0);    
            //cout << "I " << _threecenter.getDatamatrix(_i) << endl;
            //cout<< "ERIs " <<_ERIs<< endl;
            }
            
            
            
                          
            for (int alpha=0;alpha<size4c;alpha++){
                    for (int beta=0;beta<size4c;beta++){

                        double localERI = 0.0;
                        
                        for (int mu=0;mu<size4c;mu++){
                            for (int nu=0;nu<size4c;nu++){
                                
                                localERI += DMAT(mu,nu) * fourcenter[alpha][beta][mu][nu];
                                
                         }
                    }
                        
                        
                        cout << alpha+1<<":"<<beta+1<< " : " << localERI << " vs " << _ERIs(alpha,beta) << endl;
                        
                 }
            }
            
            exit(0); 
            
            
            
            
           
           // CalculateEnergy(dmatasarray);
        }
        
        
        
        
        
        // This is simply Trace(prod(ERIs*Dmat))=Totenergy
        void ERIs::CalculateEnergy(ub::vector<double> &dmatasarray){
            _ERIsenergy=0;
            ub::vector<double> ERIsasarray=_ERIs.data();
            for (unsigned _i=0;_i<ERIsasarray.size();_i++){
                _ERIsenergy+=dmatasarray[_i]*ERIsasarray[_i];
                
            }
            
            
          */  
            
        }
        
        
        void ERIs::printERIs(){
          for (unsigned i=0; i< _ERIs.size1(); i++){
                for (unsigned j=0; j< _ERIs.size2();j++){
                    cout << "ERIs [" << i<<":"<<j<<"]="<<_ERIs(i,j)<<endl;
                }
            }
        }
        
        
        
        
        
    }}
