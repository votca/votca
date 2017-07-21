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

#include <votca/xtp/gridbox.h>





namespace votca { namespace xtp {
    namespace ub = boost::numeric::ublas;
    
            void GridBox::AddtoBigMatrix(ub::matrix<double>& bigmatrix, const ub::matrix<double>& smallmatrix) {


            for (unsigned i = 0; i < ranges.size(); i++) {
                for (unsigned j = 0; j < ranges.size(); j++) {
                    ub::project(bigmatrix, ranges[i], ranges[j]) += ub::project(smallmatrix, inv_ranges[i], inv_ranges[j]);
                }
            }
            return;
        }

            

        ub::matrix<double> GridBox::ReadFromBigMatrix(const ub::matrix<double>& bigmatrix) {
            
            ub::matrix<double> _matrix = ub::zero_matrix<double>(matrix_size);
            for (unsigned i = 0; i < ranges.size(); i++) {
                for (unsigned j = 0; j < ranges.size(); j++) {
                    
                    ub::project(_matrix, inv_ranges[i], inv_ranges[j]) = ub::project(bigmatrix, ranges[i], ranges[j]);
                }
            }
            return _matrix;
        }

        void GridBox::PrepareForIntegration() {
            matrix_size = 0;

            std::vector<unsigned> start;
            std::vector<unsigned> end;

            for (unsigned i=0;i< significant_shells.size();++i) {
                const AOShell* shell=significant_shells[i];
                aoranges.push_back(ub::range(matrix_size, matrix_size+shell->getNumFunc()));
                matrix_size += shell->getNumFunc();
                start.push_back(shell->getStartIndex());
                end.push_back(shell->getStartIndex() + shell->getNumFunc());
            }
            std::vector<unsigned> startindex;
            std::vector<unsigned> endindex;
            
            if(start.size()>1){
            startindex.push_back(start[0]);
            
            for (unsigned i = 0; i < start.size() - 1; ++i) {
                
                if(end[i]!=start[i+1]){
                    startindex.push_back(start[i+1]);
                    endindex.push_back(end[i]);
                }
            }
            endindex.push_back(end[end.size()-1]);
            }
            else{
                startindex=start;
                endindex=end;
             }
            unsigned shellstart = 0;
            for (unsigned i = 0; i < startindex.size(); ++i) {
                ranges.push_back(ub::range(startindex[i], endindex[i]));
                
                unsigned size = endindex[i]-startindex[i];
                
                inv_ranges.push_back(ub::range(shellstart, shellstart + size));
                shellstart += size;
            }

            
            return;
        }
    
    
}}