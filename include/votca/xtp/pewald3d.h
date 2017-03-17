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

#ifndef VOTCA_XTP_PEWALD3D_H
#define VOTCA_XTP_PEWALD3D_H

#include <votca/csg/boundarycondition.h>
#include <votca/xtp/polartop.h>
#include <votca/xtp/xjob.h>
#include <votca/xtp/qmthread.h>
#include <votca/xtp/ewaldnd.h>

namespace CSG = votca::csg;

namespace votca { namespace xtp {
    
    // NOTE: This is not a conventional 3D Ewald summation, so use carefully
    //       (tuned for the purpose of site-energy calculations)
    // NOTE: PolarTop should be set-up with three containers: FGC, FGN, BGN
    //       MGN is set-up in constructor using the real-space c/o (from input)
    //       The topology is used to retrieve information on the sim. box (PB)
    //       All polar segments should be positioned as nearest images of the
    //       foreground charge density (FGC, FGN).
    //       All polar segments should be appropriately charged (Q00).
    // NOTE: The k-shell grouping algorithm can fail for strongly skewed boxes.
    //        
    
    class PEwald3D3D : public Ewald3DnD
    {
        
    public:
        
        PEwald3D3D(Topology *top, PolarTop *ptop, Property *opt, Logger *log);
       ~PEwald3D3D();
        
        std::string IdentifyMethod() { return "Polar 3D x 3D"; }
        void GenerateKVectors(vector<PolarSeg*> &ps1, vector<PolarSeg*> &ps2);
        
        EWD::triple<> ConvergeRealSpaceSum(vector<PolarSeg*> &target);
        EWD::triple<> ConvergeReciprocalSpaceSum(vector<PolarSeg*> &target);
        EWD::triple<> CalculateForegroundCorrection(vector<PolarSeg*> &target);
        EWD::triple<> CalculateShapeCorrection(vector<PolarSeg*> &target);
        EWD::triple<> CalculateHigherRankCorrection(vector<PolarSeg*> &target) { return EWD::triple<>(0,0,0); }        
        EWD::triple<> CalculateK0Correction(vector<PolarSeg*> &target) { return EWD::triple<>(0,0,0); }
        
        void Field_ConvergeRealSpaceSum();
        void Field_ConvergeReciprocalSpaceSum();
        void Field_CalculateForegroundCorrection();
        void Field_CalculateShapeCorrection();
        
        void Potential_ConvergeRealSpaceSum(vector<PolarSeg*> &target);
        void Potential_ConvergeReciprocalSpaceSum(vector<PolarSeg*> &target);
        void Potential_CalculateForegroundCorrection(vector<PolarSeg*> &target);
        void Potential_CalculateShapeCorrection(vector<PolarSeg*> &target);
        
        void ScanCutoff();
        
    private:
        
    };
    
    struct TinyNeighbour
    {
        TinyNeighbour(PolarSeg *seg, vec L) : _nb(seg), _L(L) { ; }
        PolarSeg *_nb;
        vec _L;
    };

}}


#endif
