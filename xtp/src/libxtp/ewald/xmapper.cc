#include <votca/xtp/ewald/xmapper.h>

namespace votca {
namespace xtp {

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

void XMapper::Gen_FGC_FGN_BGN(std::string mapfile, const Topology &top,
                              XJob *xjob) {
  // Generates foreground/background charge distribution for Ewald summations
  // 'NEW' instances of polar sites are not registered in the topology.
  // Stores the resulting 'polar topology' with the XJob class.

  XTP_LOG(Log::info, *log_) << "Mps-Mapper: Generate FGC FGN BGN" << std::flush;
  BackgroundRegion BGN(0, *log_);
  BackgroundRegion FGN(1, *log_);
  BackgroundRegion FGC(2, *log_);

  PolarMapper polmap(*log_);

  std::vector<std::string> xjobmps =
      xjob->getSegMps();  // this is a nasty hack to get the state of the
                          // foreground
  std::cout << "MPS FILE OF THIS JOB: " << xjobmps[0] << std::endl;
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
      PolarSegment mol_x = polmap.map(segment, SegId(seg_index, state_string));
      FGC.push_back(mol_x);
    } else {
      BGN.push_back(mol);
    }
    seg_index++;
  }

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
      PolarSegment mol_n = polmap.map(top.Segments()[segid], SegId(segid, "n"));
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

  // PROPAGATE SHELLS TO POLAR TOPOLOGY
  new_ptop->setFGC(fgC);
  new_ptop->setFGN(fgN);
  new_ptop->setBGN(bgN);
  new_ptop->setSegsFGC(segs_fgC);
  new_ptop->setSegsFGN(segs_fgN);
  new_ptop->setSegsBGN(segs_bgN);

  // Center polar topology
  vec center = xjob->Center();
  new_ptop->CenterAround(votca::tools::conv::bohr2nm * center);
  xjob->setPolarTop(new_ptop);
}

/*
void XMpsMap::Gen_QM_MM1_MM2(Topology *top, XJob *job, double co1, double co2,
QMThread *thread) {
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
                = this->std::mapPolSitesToSeg(psites_raw, seg,
only_active_sites); qm0.push_back(psites_std::mapped);
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
                = this->std::mapPolSitesToSeg(psites_raw, seg,
only_active_sites); mm1.push_back(psites_std::mapped);
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
                = this->std::mapPolSitesToSeg(psites_raw, seg,
only_active_sites); mm2.push_back(psites_std::mapped);
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

}  // namespace xtp
}  // namespace votca
