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

#ifndef __XTP_GRIDBOX__H
#define	__XTP_GRIDBOX__H

#include <votca/tools/vec.h>
#include <votca/xtp/aobasis.h>

#include <votca/tools/linalg.h>
#include <votca/xtp/grid_containers.h>


namespace votca { namespace xtp {

    namespace ub = boost::numeric::ublas;

        class GridBox {
            
        public: 
            
            const std::vector<tools::vec>& getGridPoints() const{return grid_pos;}
            
            const std::vector<double>& getGridWeights() const{return weights;}
            
            const std::vector<const AOShell* >& getShells() const{return significant_shells;}
            
            const std::vector<ub::range>& getAOranges() const{return aoranges;}
            
            unsigned size() const{return grid_pos.size();}
            
            unsigned Shellsize() const{return significant_shells.size();}
            
            void addGridPoint(const GridContainers::integration_grid& point){
                grid_pos.push_back(point.grid_pos);
                weights.push_back(point.grid_weight);
            };
            
            void addShell(const AOShell* shell){
              significant_shells.push_back(shell);  
            };
            
            
            void FillAOandAOGrad(ub::matrix<double>& aomatrix,ub::matrix<double>& aogradmatrix){
                
                
                return;
            }
            
            void PrepareForIntegration(){
                matrix_size=0;
                if
                std::vector<unsigned> start;
                std::vector<unsigned> end;
                
                for (const AOShell* shell: significant_shells){
                    aoranges.push_back(ub::range(matrix_size,shell->getNumFunc());
                    matrix_size+=shell->getNumFunc();
                    start.push_back(shell->getStartIndex());
                    end.push_back(shell->getStartIndex()+shell->getNumFunc());
                    
                }
                std::vector<unsigned> startindex;
                std::vector<unsigned> endindex;
                bool before=false;
                for(unsigned i=0;i<start.size()-1;++i){
                    if(before){continue;}
                    startindex.push_back(start[i]);
                    if(end[i]==start[i+1]){
                        before=true;
                        endindex.push_back(end[i+1]);
                    }
                    else{
                    endindex.push_back(end[i]);
                    before=false;
                    }
                }
                unsigned shellstart=0;
                for(unsigned i=0;i<startindex.size();++i){
                       ranges.push_back(ub::range(startindex[i],endindex[i]));
                       unsigned size=startindex[i]-endindex[i];
                       inv_ranges.push_back(ub::range(shellstart,size));
                       shellstart+=size;
                    }
                }

                return;
            }
            
            ub::matrix<double> ReadFromBigMatrix(const ub::matrix<double>& bigmatrix){
                ub::matrix<double> matrix=ub::matrix<double>(matrix_size);
                for(unsigned i=0;i<ranges.size();i++){
                    ub::project(matrix,inv_ranges[i],inv_ranges[i])=ub::project(bigmatrix,ranges[i],ranges[i]);
                }    
                
                return matrix;
            }
            
            
            void AddtoBigMatrix(ub::matrix<double>& bigmatrix,const ub::matrix<double>& smallmatrix){
                ub::matrix<double> matrix=ub::matrix<double>(matrix_size);
                for(unsigned i=0;i<ranges.size();i++){
                    ub::project(bigmatrix,ranges[i],ranges[i])+=ub::project(smallmatrix,inv_ranges[i],inv_ranges[i]);
                }    
                return;
            }
            
            
            
            private
            
                bool is_small;   
                unsigned matrix_size;
                std::vector<ub::range> aoranges;
                std::vector<ub::range> ranges;
                std::vector<ub::range> inv_ranges;
                std::vector< tools::vec > grid_pos;
                std::vector<const AOShell* > significant_shells;
                std::vector< double > weights;
                std::vector< double > densities;
                std::vector< ub::matrix<double> > dens_grad;
                
            };

        };

    }}
#endif	/* NUMERICAL_INTEGRATION_H */