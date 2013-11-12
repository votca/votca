
/***********************************************************************
 *     
 *     DL_POLY wrapper-routines needed for minimalistic I/O for 
 *     CONFIG, FIELD & HISTORY files (wrappers for VOTCA calls) 
 *     NOTE: maximum length of strings is defined by lenrec=255 (parse_module.f)
 *     
 * copyright (c) 2013, daresbury laboratory
 *     author    - andrey brukhno  july 2013
 * all rights reserved.
 * 
 * redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 
 * 1. redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 2. redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. neither the name of the <organization> nor the names of its
 *    contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 * 
 * this software is provided by the copyright holders and contributors "as is"
 * and any express or implied warranties, including, but not limited to, the
 * implied warranties of merchantability and fitness for a particular purpose
 * are disclaimed. in no event shall the copyright owner or contributors be
 * liable for any direct, indirect, incidental, special, exemplary, or
 * consequential damages (including, but not limited to, procurement of
 * substitute goods or services; loss of use, data, or profits; or business
 * interruption) however caused and on any theory of liability, whether in
 * contract, strict liability, or tort (including negligence or otherwise)
 * arising in any way out of the use of this software, even if advised of the
 * possibility of such damage.
 *     
 */

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

