/*
 *            Copyright 2009-2019 The VOTCA Development Team
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
 */

#pragma once
#ifndef STERNHEIMER_H
#define STERNHEIMER_H

#include <votca/tools/property.h>
#include <fstream>
#include <votca/xtp/eigen.h>
#include <votca/xtp/logger.h>
#include <votca/xtp/aobasis.h>
#include <votca/xtp/orbitals.h>

namespace votca{
namespace xtp{
class Orbitals;
class AOBasis;
class Sternheimer{
    public:
        Sternheimer(Orbitals& orbitals, Logger& log) : _orbitals(orbitals), _log(log){};
        
        //Calculates the dielectric matrix via the non selfconsistent Sternheimer equation for frequency w 
        Eigen::MatrixXcd CalculateDielectricMatrix(std::complex<double> w);
        
        //Calculates the screened coulomb interaction matrix from the dielectric matrix at frequency w
        Eigen::MatrixXcd CalculateScreenedCoulombOS(std::complex<double> w);
        Eigen::MatrixXcd CalculateScreenedCoulombSC(std::complex<double> w);
        //Calculates the Greens function of system
        Eigen::MatrixXcd CalculateGreensfunction(Orbitals& orbitals, double w, double eta);
        //Calculates the selfenergy of the system
        Eigen::MatrixXcd CalculateSelfEnergy(Orbitals& orbitals, double w, double eta, Eigen::VectorXd grid);
        //Calculates the spectral function
        double CalculateSpectralfunction(Eigen::MatrixXcd Selfenergy, Eigen::MatrixXcd deltaSE, double w);
        
        Eigen::VectorXcd CalculateqpCorrection(Eigen::MatrixXcd Pertubation, Orbitals& orb);
        
        Eigen::VectorXcd CalculateqpCorrection(Eigen::MatrixXcd SelfEnergy, Eigen::MatrixXcd VXC, Orbitals& orb);
        
        Eigen::MatrixXd CalculateCoulombMatrix();
        
        Eigen::MatrixXd DipoleMatrix();
        
        bool evaluate();
        
    private:
        Logger& _log;
        
        Orbitals& _orbitals;
        
        //returns the overlap matrix for all occupied states
        Eigen::MatrixXd CalculateOverlapMatrix();
        //returns the density matrix for all occupied states
        Eigen::MatrixXd CalculateDensityMatrix();
        //returns the hamiltonian matrix for all occupied states
        Eigen::MatrixXd CalculateHamiltonian();
        //Calculates coulomb matrix 
        
        //sets up the left hand side of the sternheimer equation 
        Eigen::MatrixXcd SternheimerLHS(Eigen::MatrixXcd hamiltonian, Eigen::MatrixXcd overlap, double eps, std::complex<double> omega, bool pm);
        //sets up the right hand side of the sternheimer equation
        Eigen::VectorXcd SternheimerRHS(Eigen::MatrixXcd overlap, Eigen::MatrixXcd density, Eigen::MatrixXcd pertubation, Eigen::VectorXcd coeff);
        //solves the sternheimer equation via the Biconjugate gradient method
        Eigen::VectorXcd SternheimerSolve(Eigen::MatrixXcd& LHS, Eigen::VectorXcd& RHS);
        //Calculates the matrix elements for the pertubation used in the Sternheimer equation
        Eigen::MatrixXcd DeltaVMatrix(Eigen::MatrixXcd deltaV);
        //Calculates the pertubation of the electron density using the non self-consistent one shot approach
        Eigen::MatrixXcd CalculateDeltaN(std::complex<double> w);
        //Pade Appoximation on complex grid
        Eigen::VectorXcd PadeAppox(Eigen::VectorXcd grid);
        //LHS of the linear System for analytic Greens function
        Eigen::MatrixXcd AnalyticGreensLHS(Orbitals& orbitals, double w, double eta);
        //RHS of the linear system for the analytic Greens function
        Eigen::MatrixXcd AnalyticGreensRHS(Orbitals& orbitals, double w);
        //Calculates the analytic part of the Greens function
        Eigen::MatrixXcd AnalyticGreensfunction(Orbitals& orbitals, double w, double eta);
        //Calculates the nonanalytic part of the Greens function
        Eigen::MatrixXcd NonanalyticGreensfunction(Orbitals& orbitals, double w);
        //Solves the linear system for the analytic Greens function
        Eigen::MatrixXcd AnalyticGreensolver(Eigen::MatrixXcd A, Eigen::MatrixXcd b);
        
        Eigen::MatrixXcd CalculateAnalyticSelfEnergy(Orbitals& orbitals, double w, double eta, double delta, Eigen::VectorXd grid);
        
        Eigen::MatrixXcd CalculateNonAnalyticSelfEnergy(Orbitals& orbitals, double w, double delta);

        Eigen::MatrixXcd SelfEnergyC(Orbitals orbitals, Eigen::VectorXd grid, double w, double eta);
        
        Eigen::MatrixXcd SelfEnergyEx(Orbitals orbitals);
        
        Eigen::MatrixXcd CalculateKSselfenergie(Eigen::MatrixXcd Selfenergie, Eigen::MatrixXd XC);
        
        Eigen::MatrixXd Polarisability(Eigen::MatrixXd delta_n, int x, int y, double omega);
};
}//namespace xtp
}//namespace votca
#endif /* STERNHEIMER_H */


//Eigen::MatrixXcd SternheimerW::Polarisability(std::complex<double> w, std::string gridtype){
//        
//        clock_t begin = clock();
//
//        
//        initializeMultishift(_orbitals.getBasisSetSize());
//        
//        AOBasis dftbasis=_orbitals.SetupDftBasis();
//        NumericalIntegration numint;
//        numint.GridSetup(gridtype, _orbitals.QMAtoms(), dftbasis); //For now use medium grid
//        const int& basis_size = _orbitals.getBasisSetSize();
//        std::vector<std::pair<double,const Eigen::Vector3d*> > grid = numint.getGridpoints();
//        
//        Eigen::MatrixXcd Polar = Eigen::MatrixXcd::Zero(3,3);
//        
//        double r_x = 0;
//        double rp_x= 0;
//        
//        std::cout<<"Started calculating deltaN"<<std::endl;
//        std::vector<Eigen::MatrixXcd> DeltaN=DeltaNOS(w, gridtype);
//        
//        Eigen::Vector3d gridcont=Eigen::Vector3d::Zero();
//        
//        std::complex<double> rpint;
//        std::complex<double> rint;
//        
//        Calculate second integral
//        
//        AOBasis basis=_orbitals.SetupDftBasis();
//        AODipole dipole;
//        dipole.Fill(basis);
//        
//        std::cout<<"Dipole integral complete"<<std::endl;
//        
//        Calculate first integral
//        for(int i=0;i<3;i++){
//            for(int j=0;j<3;j++){
//                for(int n=0; n<basis_size; n++){
//                    for(int m=0; m<basis_size; m++){
//                        for(int r=0; r<grid.size(); r++){  
//                            gridcont = *grid.at(r).second;
//                            r_x=gridcont(i);
//                            Polar(i,j)=Polar(i,j)+DeltaN.at(r)(n,m)*r_x*grid.at(r).first;
//                        }
//                        Polar(i,j)=Polar(i,j)+dipole.Matrix()[j](n,m);
//                    }
//                }
//                std::cout<<"Done with second index direction: "<<i<<std::endl;
//            }
//            std::cout<<"Done with first index direction: "<<i<<std::endl;
//        }
//        std::cout<<"Polarisability: "<<std::endl<<Polar<<std::endl;
//        
//        std::cout<<"Time elapsed: "<<double(clock()-begin)/CLOCKS_PER_SEC<<std::endl;
//        
//        return Polar;
//    }
//
//    std::vector<Eigen::MatrixXcd> SternheimerW::DeltaNOS(std::complex<double> w, std::string gridtype){
//        
//        
//        Setting up grid for evaluation
//        AOBasis dftbasis=_orbitals.SetupDftBasis();
//        NumericalIntegration numint;
//        numint.GridSetup(gridtype, _orbitals.QMAtoms(), dftbasis); //For now use medium grid
//        
//        std::vector<std::pair<double,const Eigen::Vector3d*> > gridpoints = numint.getGridpoints();
//        
//        std::cout<<"gridsize= "<<numint.getGridSize()<<std::endl;
//        Setting up constants
//        const int& num_occ_lvls = _orbitals.getNumberOfAlphaElectrons();
//        const int& basis_size = _orbitals.getBasisSetSize();
//        
//        Setting up constant matrices and vectors needed for the Sternheimer equation     
//        const Eigen::MatrixXd H=Hamiltonian();
//        const Eigen::MatrixXd S=OverlapMatrix();
//        const Eigen::MatrixXd p=DensityMatrix();      
//        Eigen::MatrixXcd V=Eigen::MatrixXcd::Zero(basis_size, basis_size);
//        
//        const Eigen::MatrixXd& mo_coefficients = _orbitals.MOCoefficients();
//        const Eigen::VectorXd& mo_energies = _orbitals.MOEnergies();
//        
//        Setting up container for solutions in each gridpoint 
//        
//        Eigen::MatrixXcd delta_c_p=Eigen::MatrixXcd::Zero(basis_size, basis_size);
//        Eigen::MatrixXcd delta_c_m=Eigen::MatrixXcd::Zero(basis_size, basis_size);
//        
//        Setting up container for LHS and RHS of the system
//        std::vector<Eigen::MatrixXcd> LHS_P;
//        std::vector<Eigen::MatrixXcd> LHS_M;
//        Eigen::VectorXcd RHS=Eigen::VectorXcd::Zero(basis_size);
//        
//        H_new is used if omega=0 to shift the eigenvalues away from 0
//        Eigen::MatrixXcd H_new=Eigen::MatrixXcd::Zero(basis_size,basis_size);
//        double alpha=2*(mo_energies(num_occ_lvls)-mo_energies(2));
//        
//        double pi=tools::conv::Pi;
//        Initilazing solution vector and container for solution
//        
//        std::vector<Eigen::MatrixXcd> result;
//        Eigen::MatrixXcd delta_n=Eigen::MatrixXcd::Zero(basis_size, basis_size);
//        
//        Iterating over all occupied states and saving the LHS
//        for(int v=0; v<num_occ_lvls; v++){
//            
//            Setting up LHS since it is not dependent on r
//            if(w==0.0){                      
//                H_new=H+alpha*S*p.transpose();        
//                LHS_P.push_back(SternheimerLHS(H_new,S,mo_energies(v),w,true));
//                LHS_M.push_back(SternheimerLHS(H_new,S,mo_energies(v),w,false));
//            }
//            else{
//                LHS_P.push_back(SternheimerLHS(H,S,mo_energies(v),w,true));
//                LHS_M.push_back(SternheimerLHS(H,S,mo_energies(v),w,false));
//            }
//        }    
//            Iterate over all spacial gridpoints
//        for(int r=0; r<gridpoints.size(); r++){
//            
//            
//            std::cout<<(*gridpoints[r])<<std::endl;
//            V=CoulombMatrix(*gridpoints.at(r).second);
//            
//            for(int v=0; v<num_occ_lvls; v++){
//    
//                RHS=SternheimerRHS(S,p,V,mo_coefficients.col(v));
//                
//                delta_c_p.col(v)=SternheimerSolve(LHS_P.at(v), RHS);
//                delta_c_m.col(v)=SternheimerSolve(LHS_M.at(v), RHS);
//                
//            }
//            
//            for(int v=0; v<basis_size; v++){
//                for(int i=0;i<basis_size; i++){
//                    for(int j=0;j<basis_size; j++){
//                        delta_n(i,j)=delta_n(i,j)+2*mo_coefficients(i,v)*delta_c_p(j,v)+2*mo_coefficients(i,v)*delta_c_m(j,v);                        
//                    }
//                }
//            }
//        
//            result.push_back(delta_n);
//            delta_n=Eigen::MatrixXcd::Zero(basis_size,basis_size);
//            
//            if(r%100==0){
//                
//                std::cout<<"I'm in iteration "<<r<<std::endl;
//                
//            }
//        }
//        
//        clock_t end = clock();
//        
//        
//        return result;
//        
//    }