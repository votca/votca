#ifndef _dlp_io_layer_H
#define _dlp_io_layer_H

struct FullCellT {
  double uc[9];
  double rc[9];
  double ac[9];
  double pc[10];
};

struct FrameSiteT {
  char name[9];
  char type[9];

  double r[3];
  double v[3];
  double f[3];

  double m,q;
  int    id,im,ig;
};

struct FrameSpecsT {
  char  title[81];
  double cell[9];
  double tstep,energy;
  int  nstep,nsites,imcon,keytrj;
};

struct FieldSiteT {
  char name[9];
  char type[9];
  double m,q;
  int    idmol,idgrp,ifrzn,nrept;
};

struct MolecSpecsT {
  char name[81];
  int  id,nrept,nsites,ngroups;
};

struct FieldSpecsT {
  char title[81];
  char units[9];
  int  ineut;
  int  nmols,natms;
  int  nbonds,nangles,ndhdrs,nexcls;
};

struct char81 {
  char str[81];
};

struct char9 {
  char str[9];
};

extern "C"
{
  void field_scan_(int*, int*, int*, int*);

  void field_read_(int*, FieldSpecsT*, MolecSpecsT*, FieldSiteT*);

  void field_write_(int*, FieldSpecsT*, MolecSpecsT*, FieldSiteT*);

  void traj_read_(int*, FrameSpecsT*, FrameSiteT*);

  void traj_write_(int*, FrameSpecsT*, FrameSiteT*);

  void conf_read_(int*, FrameSpecsT*, FrameSiteT*);

  void conf_write_(int*, FrameSpecsT*, FrameSiteT*);
}

#endif

