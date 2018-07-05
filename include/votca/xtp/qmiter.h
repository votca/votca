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

#ifndef __QMITER__H
#define	__QMITER__H


#include <votca/ctp/xjob.h>
#include <votca/ctp/xinductor.h>
#include <votca/xtp/qminterface.h>
// add gwbse header for excited state support
#include <votca/xtp/gwbse.h>
#include <votca/xtp/qmpackagefactory.h>
#include <votca/xtp/orbitals.h>


namespace votca { namespace xtp {


    
class QMMIter
{

public:

    QMMIter(int id) : _id(id), _hasdRdQ(false), _hasQM(false), _hasMM(false)  { ; }
   ~QMMIter() { ; }

   void ConvertPSitesToQMAtoms(std::vector< ctp::PolarSeg* > &, std::vector< QMAtom* > &);
   void ConvertQMAtomsToPSites(std::vector< QMAtom* > &, std::vector< ctp::PolarSeg* > &);
   void UpdatePosChrgFromQMAtoms(std::vector< QMAtom* > &, std::vector< ctp::PolarSeg* > &);  
   void UpdateMPSFromGDMA( std::vector<std::vector<double> > &multipoles,  std::vector< ctp::PolarSeg* > &psegs);
   void GenerateQMAtomsFromPolarSegs(ctp::PolarTop *ptop, Orbitals &orb, bool split_dpl, double dpl_spacing);   

   void setdRdQ(double dR_RMS, double dQ_RMS, double dQ_SUM);
   void setQMSF(double energy_QM, double energy_SF, double energy_GWBSE);
   void setE_FM(double ef00, double ef01, double ef02, 
                  double ef11, double ef12, double em0,
                  double em1,  double em2, double efm);

   double getRMSdR() { return _dR_RMS; }
   double getRMSdQ() { return _dQ_RMS; }
   double getSUMdQ() { return _dQ_SUM; }
   int getId() { return _id;}

   double getSFEnergy() { assert(_hasQM); return _e_SF; }
   double getFMEnergy() { assert(_hasMM); return _e_fm_; }
   double getQMEnergy() { assert(_hasQM); return _e_QM; }
   double getGWBSEEnergy() { assert(_hasGWBSE); return _e_GWBSE; }
   double getMMEnergy();
   double getQMMMEnergy();


private:

    int    _id;

    bool   _hasdRdQ;
    bool   _hasQM;
    bool   _hasMM;
    bool   _hasGWBSE;

    double _dR_RMS;
    double _dQ_RMS;
    double _dQ_SUM;       

    double _e_QM;
    double _e_SF;
    double _e_GWBSE;

    double _ef_00;
    double _ef_01;
    double _ef_02;
    double _ef_11;
    double _ef_12;
    double _em_0_;
    double _em_1_;
    double _em_2_;
    double _e_fm_;


};
    
 

}}

#endif
