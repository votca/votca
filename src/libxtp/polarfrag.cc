#include <votca/xtp/polarfrag.h>
#include <votca/xtp/dmaspace.h>

namespace votca {
namespace xtp {


PolarFrag::~PolarFrag()  { 
    clear(); /* Polar sites cleaned by PolarSeg */
    if (_indu_cg_site != NULL) delete _indu_cg_site;
    if (_perm_cg_site != NULL) delete _perm_cg_site;
}
    
    
const votca::tools::vec &PolarFrag::CalcPosCenterOfGeom() {    
    _pos = votca::tools::vec(0,0,0);    
    for (unsigned int i = 0; i < this->size(); ++i) {
        _pos += (*this)[i]->getPos();
    }
    if (this->size() > 0) _pos /= double(this->size());
    return _pos;
}
    
    
const votca::tools::vec &PolarFrag::CalcPosPolarWeights() {
    // Establish total weights
    double sum_iso_p = 0.0;
    double sum_abs_q = 0.0;
    for (unsigned int i = 0; i < this->size(); ++i) {
        sum_iso_p += (*this)[i]->getIsoP();
        sum_abs_q += std::abs((*this)[i]->getQ00());
    }

    _pos = votca::tools::vec(0,0,0);
    // Any polar sites in here? If not, return (0,0,0)
    if (this->size() < 1) ;
    // Noteworthy charge and polarizability? If not, return CoG
    else if (sum_abs_q <= 1e-2 && sum_iso_p <= 1e-2) {
        _pos = this->CalcPosCenterOfGeom();
    }
    // Else: return weighted 0.5*(center-of-polarity + center-of-charge)
    else {
        for (unsigned int i = 0; i < this->size(); ++i) {
            double weight = 0.0;
            if (sum_abs_q < 1e-2) {
                assert(sum_iso_p >= 1e-2 && "<CalcPosPolarWeights> P-ERROR");
                weight = (*this)[i]->getIsoP()/sum_iso_p;
            }
            else if (sum_iso_p < 1e-2) {
                assert(sum_abs_q >= 1e-2 && "<CalcPosPolarWeights> Q-ERROR");
                weight = std::abs((*this)[i]->getQ00())/sum_abs_q;
            }
            else {
                weight = 0.5 * (
                    (*this)[i]->getIsoP()/sum_iso_p
                  + std::abs((*this)[i]->getQ00())/sum_abs_q);
            }
            _pos += (*this)[i]->getPos() * weight;
        }
    }
    return _pos;
}
    
    
void PolarFrag::GeneratePermInduCgSite(bool do_cg_polarizabilities) {
    // ATTENTION The same method appears in <PolarSeg>
    assert(!do_cg_polarizabilities && "NOT IMPLEMENTED, NOT NEEDED?");
    // Collapse multipole moments : position, rank L
    vec target_pos = _pos;
    int state = 0;
    int L = 2;
    vector<double> QCG(L*L+2*L+1, 0.0);  // permanent
    vector<double> uQCG(L*L+2*L+1, 0.0); // induced

    for (PolarFrag::iterator pit = begin();
        pit < end(); ++pit) {
        // PERMANENT MOMENTS
        // Convert real to complex moments            
        vector<double> Qlm = (*pit)->getQs(0);
        DMA::ComplexSphericalMoments Xlm(Qlm);
        // Shift moments
        DMA::MomentShift mshift;
        vec shift = target_pos - (*pit)->getPos();
        DMA::RegularSphericalHarmonics Clm(-shift);
        vector<DMA::cmplx> Xlm_shifted = mshift.Shift(Xlm, Clm);            
        // Convert complex to real moments & add to base
        DMA::RealSphericalMoments Qlm_shifted(Xlm_shifted);
        Qlm_shifted.AddToVector(QCG);

        // INDUCED MOMENTS
        // Convert real to complex moments
        vec u1 = (*pit)->getU1();
        vector<double> uQlm(L*L+2*L+1, 0.0);
        uQlm[1] = u1.getZ(); // NOTE order is z-x-y == 10-11c-11s
        uQlm[2] = u1.getX();
        uQlm[3] = u1.getY();
        DMA::ComplexSphericalMoments uXlm(uQlm);
        // Shift moments
        DMA::RegularSphericalHarmonics uClm(-shift);
        vector<DMA::cmplx> uXlm_shifted = mshift.Shift(uXlm, uClm);
        // Convert complex to real moments & add to base
        DMA::RealSphericalMoments uQlm_shifted(uXlm_shifted);
        uQlm_shifted.AddToVector(uQCG);
    }
        
    // Collapse polarizabilities
    votca::tools::matrix PCG;
    PCG.ZeroMatrix();
    
    // Zero induced dipole
    vec u1_cg_red = vec(0,0,0);
        
    // Generate new coarse-grained site from the above
    APolarSite *indu_cg_site = new APolarSite(this->getId(), "FGU");
    APolarSite *perm_cg_site = new APolarSite(this->getId(), "FGP");
    
    indu_cg_site->setResolution(APolarSite::coarsegrained);
    indu_cg_site->setPos(target_pos);
    indu_cg_site->setRank(L);
    
    perm_cg_site->setResolution(APolarSite::coarsegrained);
    perm_cg_site->setPos(target_pos);
    perm_cg_site->setRank(L);
    
    // ATTENTION Save INDUCED   moments as PERMANENT moments (<indu_cg_site>)
    // ATTENTION Save PERMANENT moments as PERMANENT moments (<perm_cg_site>)
    indu_cg_site->setQs(uQCG, state);
    indu_cg_site->setPs(PCG, state);
    indu_cg_site->setU1(u1_cg_red);
    indu_cg_site->Charge(state);
    
    perm_cg_site->setQs(QCG, state);
    perm_cg_site->setPs(PCG, state);
    perm_cg_site->setU1(u1_cg_red);
    perm_cg_site->Charge(state);
    
    // Deallocate previously allocated sites
    if (_indu_cg_site != NULL) delete _indu_cg_site;
    if (_perm_cg_site != NULL) delete _perm_cg_site;
    _indu_cg_site = indu_cg_site;
    _perm_cg_site = perm_cg_site;
    return;
}
    
    
    
    
    
    
    
    
}}
