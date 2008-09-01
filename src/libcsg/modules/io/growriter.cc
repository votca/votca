// 
// File:   growriter.cc
// Author: ruehle
//
// Created on January 8, 2008, 18:09 PM
//

#include <stdio.h>
#include <string>
#include "growriter.h"

using namespace std;

void GROWriter::Open(string file, bool bAppend)
{
    _out = fopen(file.c_str(), bAppend ? "at" : "wt");
}

void GROWriter::Close()
{
    fclose(_out);
}

void GROWriter::Write(Configuration *conf)
{
    char nm[6],format[100];
  int  ai,i,resnr,l,vpr;
  Topology *top = conf->getTopology();

  fprintf (_out,"%s\n","what a nice title");
  fprintf (_out,"%5d\n",conf->getTopology()->BeadCount());
  
  bool v = false; // we don't write velocities!  
  int pr = 3; // precision of writeout
  
  /* build format string for printing, 
     something like "%8.3f" for x and "%8.4f" for v */
  /*if (pr<0)
    pr=0;
  if (pr>30)
    pr=30;*/
  l=pr+5;
  vpr=pr+1;
  if (v)
    sprintf(format,"%%%d.%df%%%d.%df%%%d.%df%%%d.%df%%%d.%df%%%d.%df\n",
	    l,pr,l,pr,l,pr,l,vpr,l,vpr,l,vpr);
  else
    sprintf(format,"%%%d.%df%%%d.%df%%%d.%df\n",l,pr,l,pr,l,pr);
  
  for (i=0; i< conf->getTopology()->BeadCount(); i++) {
    resnr=top->getBead(i)->getResnr();
    string resname = top->getResidue(resnr)->getName();
    string atomname = top->getBead(i)->getName();
    
    fprintf(_out,"%5d%-5.5s%5.5s%5d",
            (resnr+1)%100000,resname.c_str(),atomname.c_str(),(i+1)%100000);
    /* next fprintf uses built format string */
    vec r = conf->Pos(i);
    vec vv = vec(0,0,0);
    
    if (v)
      fprintf(_out,format,
	      r.getX(), r.getY(), r.getZ(), vv.getX(), vv.getY(), vv.getZ());
    else
      fprintf(_out,format,
	      r.getX(), r.getY(), r.getZ());
  }

  // write the boy
  matrix box = conf->getBox();
  
  if (pr<5) 
    pr=5;
  l=pr+5;
  
  if (box[0][1] || box[0][2] || box[1][0] || box[1][2] ||
      box[2][0] || box[2][1]) {
    sprintf(format,"%%%d.%df%%%d.%df%%%d.%df"
	    "%%%d.%df%%%d.%df%%%d.%df%%%d.%df%%%d.%df%%%d.%df\n",
	    l,pr,l,pr,l,pr,l,pr,l,pr,l,pr,l,pr,l,pr,l,pr);
    fprintf(_out,format,
	    box[0][0],box[1][1],box[2][2],
	    box[0][1],box[0][2],box[1][0],
	    box[1][2],box[2][0],box[2][1]);
  } else {
    sprintf(format,"%%%d.%df%%%d.%df%%%d.%df\n",l,pr,l,pr,l,pr);
    fprintf(_out,format,
	    box[0][0],box[1][1],box[2][2]);
  }
  fflush(_out);
}
