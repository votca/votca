#include <votca/xtp/ewald/xmapper.h>


namespace votca { namespace xtp {

/*
void XMpsMap::Gen_BGN(Topology *top, PolarTop *new_ptop, QMThread *thread) {
    // Generates (periodic) background charge distribution for Ewald summations
    // 'NEW' instances of polar sites are not registered in the topology.
    // Modifies the polar topology (*PolarTop*) passed during function call
    
    // DECLARE TARGET CONTAINERS
    std::vector<PolarSeg*> bgN;
    std::vector<Segment*> segs_bgN;
    std::vector<Segment*>::iterator sit;
    
    // PARTITION SEGMENTS ONTO BACKGROUND + FOREGROUND
    segs_bgN.reserve(top->Segments().size());
    for (sit = top->Segments().begin();
         sit < top->Segments().end();
         ++sit) {        
        Segment *seg = *sit;        
        segs_bgN.push_back(seg);      
    }
    
    // CREATE POLAR SITES FOR FOREGROUND + BACKGROUND
    // Background
    //bool only_active_sites = true;
    bgN.reserve(segs_bgN.size());
    for (sit = segs_bgN.begin(); sit < segs_bgN.end(); ++sit) {
        Segment *seg = *sit;
        // Look up appropriate set of polar sites
        std::string mps = _segId_mpsFile_n[seg->getId()];
        std::vector<APolarSite*> psites_raw  = _mpsFile_pSites[mps];
        PolarSeg *new_pseg = this->std::mapPolSitesToSeg(psites_raw, seg);
        bgN.push_back(new_pseg);        
    }
    
    // PROPAGATE SHELLS TO POLAR TOPOLOGY
    new_ptop->setBGN(bgN);    
    new_ptop->setSegsBGN(segs_bgN);
    return;
}
*/

void XMapper::Gen_FGC_FGN_BGN(std::string mapfile, const Topology &top, XJob *xjob) {
    // Generates foreground/background charge distribution for Ewald summations
    // 'NEW' instances of polar sites are not registered in the topology.
    // Stores the resulting 'polar topology' with the XJob class.
    
XTP_LOG(Log::info, *log_) << "Mps-Mapper: Generate FGC FGN BGN" << std::flush;
    BackgroundRegion BGN(0, *log_);
    BackgroundRegion FGN(1, *log_);
    BackgroundRegion FGC(2, *log_);

    PolarMapper polmap(*log_);
    XTP_LOG(Log::info, *log_) << "Reading from" << mapfile << std::flush;


  std::vector<std::string> xjobmps =
      xjob->getSegMps();  // this is a nasty hack to get the state of the
                         // foreground
  std::cout << "MPS FILE OF THIS JOB: " << xjobmps[0] << std::endl;
  std::cout << "!!!!! STILL NEED TO DETERMINE QMSTATE FROM HERE !!!!! "
            << std::endl;
  tools::Tokenizer toker(xjobmps[0], "_.");
  std::vector<std::string> split = toker.ToVector();

  std::string state_string = *(split.rbegin() + 1);



    polmap.LoadMappingFile(mapfile);
    Index seg_index = 0;
    for (auto segment : top.Segments()) {
      PolarSegment mol = polmap.map(segment, SegId(seg_index, "n"));
      // segments in the foreground
      if (xjob->isInCenter(segment.getId())) {
        FGN.push_back(mol);
        PolarSegment mol_x =
            polmap.map(segment, SegId(seg_index, state_string));
        FGC.push_back(mol_x);
      } else {
        BGN.push_back(mol);
      }
      seg_index++;
    }

    XTP_LOG(Log::info, *log_) << "Generated BGN, FGN and FGC" << std::flush;

    // DECLARE TARGET CONTAINERS
    // Convert this to old PolarTop
    Topology new_top = top;
    PolarTop *new_ptop = new PolarTop(&new_top);   
    std::vector<PolarSeg *> fgC;
    std::vector<PolarSeg *> fgN;
    std::vector<PolarSeg *> bgN;
    std::vector<Segment *> segs_fgC;
    std::vector<Segment *> segs_fgN;
    std::vector<Segment *> segs_bgN;

    segs_fgC.reserve(xjob->getSegments().size());
    segs_fgN.reserve(xjob->getSegments().size());
    segs_bgN.reserve(top.Segments().size() - xjob->getSegments().size());
    for (auto segment : top.Segments()) {
      // Foreground
      if (xjob->isInCenter(segment.getId())) {
        segs_fgN.push_back(&segment);
        segs_fgC.push_back(&segment);
      }
      // Background
      else {
        segs_bgN.push_back(&segment);
      }
    }
    XTP_LOG(Log::info, *log_)
        << "Converted BGN, FGN and FGC segments" << std::flush;

    bool only_active_sites = false;
    fgN.reserve(segs_fgN.size());
    // get all NEW PolarSites of this segment and convert them to OLD PolarSites
    int state = 0;  // neutral ground state
    for (auto &segment : FGN) {
      std::vector<APolarSite *> psites;
      psites.reserve(segment.size());
      for (auto site : segment) {
        APolarSite *psite = new APolarSite();
        psite->ConvertFromPolarSite(site, state);
        psite->Charge(state);
        psites.push_back(psite);
      }
      // now make an OLD PolarSeg from the new PolarSegment
      PolarSeg *new_pseg = new PolarSeg(int(segment.getId()), psites);
      fgN.push_back(new_pseg);
    }

    XTP_LOG(Log::info, *log_) << "Created FGN Polar segments" << std::flush;

    // now the excited one
    fgC.reserve(segs_fgC.size());
    if (state_string == "n") state = 0;
    if (state_string == "e") state = -1;
    if (state_string == "h") state = 1;
    for (auto &segment : FGC) {
      std::vector<APolarSite *> psites;
      psites.reserve(segment.size());
      for (auto site : segment) {
        APolarSite *psite = new APolarSite();
        psite->ConvertFromPolarSite(site, state);
        psite->Charge(state);
        psites.push_back(psite);
      }

      if (state != 0) {
        //// FGC segments need APolarSites for neutral variant as well
        int segid = segment.getId();
        PolarSegment mol_n =
            polmap.map(top.Segments()[segid], SegId(segid, "n"));
        int pidx = 0;
        for (auto site_n : mol_n) {
          psites[pidx]->addNeutral(site_n);
          pidx++;
        }
      }
      // now make an OLD PolarSeg from the new PolarSegment
      PolarSeg *new_pseg = new PolarSeg(int(segment.getId()), psites);
      fgC.push_back(new_pseg);
    }
    XTP_LOG(Log::info, *log_) << "Created FGC Polar segments" << std::flush;

    // segments in BGN are NEW POLARSEGMENTS
    state = 0;  // all neutral again
    for (auto segment : BGN) {
      // get all NEW PolarSites of this segment and convert them to OLD
      // PolarSites
      std::vector<APolarSite *> psites;
      psites.reserve(segment.size());
      for (auto site : segment) {
        APolarSite *psite = new APolarSite();
        psite->ConvertFromPolarSite(site, state);
        psite->Charge(state);  // set state of this site: ground state
        psites.push_back(psite);
      }

      // now make an OLD PolarSeg from the new PolarSegment
      PolarSeg *new_pseg = new PolarSeg(int(segment.getId()), psites);
      bgN.push_back(new_pseg);
    }
    XTP_LOG(Log::info, *log_) << "Created BGN Polar segments" << std::flush;

    // PROPAGATE SHELLS TO POLAR TOPOLOGY
    new_ptop->setFGC(fgC);
    new_ptop->setFGN(fgN);
    new_ptop->setBGN(bgN);
    new_ptop->setSegsFGC(segs_fgC);
    new_ptop->setSegsFGN(segs_fgN);
    new_ptop->setSegsBGN(segs_bgN);

    XTP_LOG(Log::info, *log_) << "Propagated to polar topology" << std::flush;

    // Center polar topology
    vec center = xjob->Center();
    new_ptop->CenterAround(votca::tools::conv::bohr2nm* center);
    xjob->setPolarTop(new_ptop);
    XTP_LOG(Log::info, *log_) << "Centered polar topology" << std::flush;
 

    /// OLD STUFF


/*

    // DECLARE TARGET CONTAINERS
    PolarTop *new_ptop = new PolarTop(top);    
    std::vector<PolarSeg*> fgC;
    std::vector<PolarSeg*> fgN;
    std::vector<PolarSeg*> bgN;
    std::vector<Segment*> segs_fgC;
    std::vector<Segment*> segs_fgN;
    std::vector<Segment*> segs_bgN;
    std::vector<Segment*>::iterator sit;
    
    // PARTITION SEGMENTS ONTO BACKGROUND + FOREGROUND
    segs_fgC.reserve(job->getSegments().size());
    segs_fgN.reserve(job->getSegments().size());
    segs_bgN.reserve(top->Segments().size()-job->getSegments().size());
    for (sit = top->Segments().begin();
         sit < top->Segments().end();
         ++sit) {        
        Segment *seg = *sit;        
        // Foreground
        if (job->isInCenter(seg->getId())) {
            segs_fgN.push_back(seg);
            segs_fgC.push_back(seg);
        }
        // Background
        else {
            segs_bgN.push_back(seg);
        }        
    }
    
    // CREATE POLAR SITES FOR FOREGROUND + BACKGROUND
    // Foreground
    bool only_active_sites = false;
    fgC.reserve(segs_fgC.size());
    fgN.reserve(segs_fgN.size());
    for (unsigned int i = 0; i < job->getSegments().size(); ++i) {        
        Segment *seg = job->getSegments()[i];
        // Charged => mps-file from job
        std::string mps_C = job->getSegMps()[i];
        std::vector<APolarSite*> psites_raw_C 
            = this->GetOrCreateRawSites(mps_C,thread);
        PolarSeg *psegC 
            = this->std::mapPolSitesToSeg(psites_raw_C, seg, only_active_sites);        
        fgC.push_back(psegC);
        // Neutral => look up mps file
        std::string mps_N = _segId_mpsFile_n[seg->getId()];
        std::vector<APolarSite*> psites_raw_N  = _mpsFile_pSites[mps_N];
        PolarSeg *psegN
            = this->std::mapPolSitesToSeg(psites_raw_N, seg, only_active_sites);
        fgN.push_back(psegN);
    }
    // Background
    only_active_sites = true;
    bgN.reserve(segs_bgN.size());
    for (sit = segs_bgN.begin(); sit < segs_bgN.end(); ++sit) {
        Segment *seg = *sit;
        // Look up appropriate set of polar sites
        std::string mps = _segId_mpsFile_n[seg->getId()];
        std::vector<APolarSite*> psites_raw  = _mpsFile_pSites[mps];
        PolarSeg *psegN = this->std::mapPolSitesToSeg(psites_raw, seg);
        bgN.push_back(psegN);        
    }
    
    // PROPAGATE SHELLS TO POLAR TOPOLOGY
    new_ptop->setFGC(fgC);
    new_ptop->setFGN(fgN);
    new_ptop->setBGN(bgN);    
    new_ptop->setSegsFGC(segs_fgC);
    new_ptop->setSegsFGN(segs_fgN);
    new_ptop->setSegsBGN(segs_bgN);    
    // Center polar topology
    vec center = job->Center();
    new_ptop->CenterAround(center);
    job->setPolarTop(new_ptop);
    */
    
}

/*
void XMpsMap::Gen_FGC_Load_FGN_BGN(Topology *top, XJob *job, std::string archfile, 
    QMThread *thread) {
    
    // LOAD BACKGROUND POLARIZATION STATE
    PolarTop bgp_ptop = PolarTop(top);
    bgp_ptop.LoadFromDrive(archfile);
    
    // SANITY CHECKS I
    if (bgp_ptop.QM0().size() || bgp_ptop.MM1().size() || bgp_ptop.MM2().size()
        || bgp_ptop.FGC().size() || bgp_ptop.FGN().size()) {
        std::cout << std::endl;
        std::cout << "ERROR The polar topology from '" << archfile 
            << "' contains more than just background. ";
        std::cout << std::endl;
        throw std::runtime_error
            ("Sanity checks I in XMpsMap::Gen_FGC_Load_FGN_BGN failed.");
    }
    
    // DECLARE TARGET CONTAINERS
    PolarTop *new_ptop = new PolarTop(top);    
    std::vector<PolarSeg*> fgC;
    std::vector<PolarSeg*> fgN;
    std::vector<PolarSeg*> bgN;
    std::vector<PolarSeg*>::iterator psit;    
    std::vector<Segment*> segs_fgC;
    std::vector<Segment*> segs_fgN;
    std::vector<Segment*> segs_bgN;
    std::vector<Segment*>::iterator sit;
    
    // PARTITION SEGMENTS ONTO BACKGROUND + FOREGROUND
    segs_fgC.reserve(job->getSegments().size());
    segs_fgN.reserve(job->getSegments().size());
    segs_bgN.reserve(top->Segments().size()-job->getSegments().size());
    for (sit = top->Segments().begin();
         sit < top->Segments().end();
         ++sit) {        
        Segment *seg = *sit;        
        // Foreground
        if (job->isInCenter(seg->getId())) {
            segs_fgN.push_back(seg);
            segs_fgC.push_back(seg);
        }
        // Background
        else {
            segs_bgN.push_back(seg);
        }        
    }
    
    // CREATE POLAR SITES FOR FGC
    bool only_active_sites = false;
    fgC.reserve(segs_fgC.size());
    for (unsigned int i = 0; i < job->getSegments().size(); ++i) {        
        Segment *seg = job->getSegments()[i];
        // Charged => mps-file from job
        std::string mps_C = job->getSegMps()[i];
        std::vector<APolarSite*> psites_raw_C 
            = this->GetOrCreateRawSites(mps_C,thread);
        PolarSeg *psegC
            = this->std::mapPolSitesToSeg(psites_raw_C, seg, only_active_sites);        
        fgC.push_back(psegC);
    }
    
    // DIVIDE POLAR SEGMENTS FROM RESURRECTED BACKGROUND ONTO FGN, BGN
    for (psit = bgp_ptop.BGN().begin(); psit < bgp_ptop.BGN().end(); ++psit) {
        // Move to (neutral) foreground?
        if (job->isInCenter((*psit)->getId())) {
            fgN.push_back(*psit);
        }
        // Move to (neutral) background?
        else {
            bgN.push_back(*psit);
        }
    }
    
    // SANITY CHECKS II
    if ((fgN.size() != fgC.size())
        || (fgN.size() + bgN.size() != top->Segments().size())
        || (bgp_ptop.BGN().size() != top->Segments().size())) {
        std::cout << std::endl;
        std::cout << "ERROR Is the background binary compatible with this system? ";
        std::cout << "(archive = '" << archfile << "')";
        std::cout << std::endl;
        throw std::runtime_error
            ("Sanity checks II in XMpsMap::Gen_FGC_Load_FGN_BGN failed.");
    }
    
    // PROPAGATE SHELLS TO POLAR TOPOLOGY
    new_ptop->setFGC(fgC);
    new_ptop->setFGN(fgN);
    new_ptop->setBGN(bgN);    
    new_ptop->setSegsFGC(segs_fgC);
    new_ptop->setSegsFGN(segs_fgN);
    new_ptop->setSegsBGN(segs_bgN);
    // Remember to remove ownership from temporary bgp_ptop
    bgp_ptop.RemoveAllOwnership();
    // Center polar topology
    vec center = job->Center();
    new_ptop->CenterAround(center);
    job->setPolarTop(new_ptop);
    return;
}


void XMpsMap::Gen_QM_MM1_MM2(Topology *top, XJob *job, double co1, double co2, QMThread *thread) {
    // Generates QM MM1 MM2, centered around job->Center().
    // 'NEW' instances of polar sites are not registered in the topology.
    // Stores the resulting 'polar topology' with the XJob class.    
    
    // TARGET CONTAINERS
    PolarTop *new_ptop = new PolarTop(top);
    
    std::vector<PolarSeg*> qm0;
    std::vector<PolarSeg*> mm1;
    std::vector<PolarSeg*> mm2;
    
    std::vector<Segment*> segs_qm0;
    std::vector<Segment*> segs_mm1;
    std::vector<Segment*> segs_mm2;    
    std::vector<Segment*> ::iterator sit;
    
    // PARTITION SEGMENTS ONTO SHELLS    
    for (sit = top->Segments().begin();
         sit < top->Segments().end();
         ++sit) {
        
        Segment *seg = *sit;
        
        // QM0 SHELL
        if (job->isInCenter(seg->getId())) {
            segs_qm0.push_back(seg);
        } 
        // MM1 SHELL
        else if (job->isWithinDist(seg->getPos(),co1,top)) {
            segs_mm1.push_back(seg);
        }
        // MM2 SHELL
        else if (job->isWithinDist(seg->getPos(),co2,top)) {
            segs_mm2.push_back(seg);
        }
        // BEYOND
        else {
            ;
        }
    }
    
    // CREATE POLAR SEGMENTS FROM SHELLS    
    // ... QM0 SHELL
    bool only_active_sites = false;
    qm0.reserve(segs_qm0.size());
    for (unsigned int i = 0; i < job->getSegments().size(); ++i) {        
        Segment *seg = job->getSegments()[i];
        std::vector<APolarSite*> psites_raw 
                = this->GetOrCreateRawSites(job->getSegMps()[i],thread);
        PolarSeg *psites_std::mapped
                = this->std::mapPolSitesToSeg(psites_raw, seg, only_active_sites);        
        qm0.push_back(psites_std::mapped);        
    }
    // ... MM1 SHELL
    only_active_sites = true;
    mm1.reserve(segs_mm1.size());
    for (sit = segs_mm1.begin(); sit < segs_mm1.end(); ++sit) {
        Segment *seg = *sit;
        // Look up appropriate set of polar sites
        std::string mps = _segId_mpsFile_n[seg->getId()];
        std::vector<APolarSite*> psites_raw  = _mpsFile_pSites[mps];
        PolarSeg *psites_std::mapped
                = this->std::mapPolSitesToSeg(psites_raw, seg, only_active_sites);
        mm1.push_back(psites_std::mapped);        
    }
    // ... MM2 SHELL
    only_active_sites = true;
    mm2.reserve(segs_mm2.size());
    for (sit = segs_mm2.begin(); sit < segs_mm2.end(); ++sit) {
        Segment *seg = *sit;
        // Look up appropriate set of polar sites
        std::string mps = _segId_mpsFile_n[seg->getId()];
        std::vector<APolarSite*> psites_raw  = _mpsFile_pSites[mps];
        PolarSeg *psites_std::mapped
                = this->std::mapPolSitesToSeg(psites_raw, seg, only_active_sites);
        mm2.push_back(psites_std::mapped);
    }
    
    // PROPAGATE SHELLS TO POLAR TOPOLOGY
    new_ptop->setQM0(qm0);
    new_ptop->setMM1(mm1);
    new_ptop->setMM2(mm2);
    
    new_ptop->setSegsQM0(segs_qm0);
    new_ptop->setSegsMM1(segs_mm1);
    new_ptop->setSegsMM2(segs_mm2);
    
    // CENTER AROUND QM REGION
    new_ptop->CenterAround(job->Center());
    job->setPolarTop(new_ptop);
    
}*/


}}
