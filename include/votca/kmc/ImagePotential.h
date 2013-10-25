/*
 * Copyright 2009-2013 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

void Graph::LoadGraph {
    
    vector<NodeMultiple*> node;
    
    // Load nodes
    votca::tools::Database db;
    db.Open( _filename );
    if(verbose >= 1) {cout << "LOADING GRAPH" << endl << "database file: " << _filename << endl; }
    votca::tools::Statement *stmt = db.Prepare("SELECT _id-1, name, posX, posY, posZ, UnCnN"+_carriertype+", UcNcC"+_carriertype+",eAnion,eNeutral,eCation,ucCnN"+_carriertype+" FROM segments;");

    int i=0;
    while (stmt->Step() != SQLITE_DONE)
    {
        NodeMultiple *newNode = new NodeMultiple();
        node.push_back(newNode);

        int newid = stmt->Column<int>(0);
        string name = stmt->Column<string>(1);
        node[i]->id = newid;
        myvec nodeposition = myvec(stmt->Column<double>(2)*1E-9, stmt->Column<double>(3)*1E-9, stmt->Column<double>(4)*1E-9); // converted from nm to m
        node[i]->position = nodeposition;
        node[i]->reorg_intorig = stmt->Column<double>(5); // UnCnN
        node[i]->reorg_intdest = stmt->Column<double>(6); // UcNcC
        double eAnion = stmt->Column<double>(7);
        double eNeutral = stmt->Column<double>(8);
        double eCation = stmt->Column<double>(9);
        double internalenergy = stmt->Column<double>(10); // UcCnN
        double siteenergy = 0;
        if(_carriertype == "e")
        {
            siteenergy = eAnion + internalenergy;
        }
        else if(_carriertype == "h")
        {
            siteenergy = eCation + internalenergy;
        }
        node[i]->siteenergy = siteenergy;
        if (votca::tools::wildcmp(_injection_name.c_str(), name.c_str()))
        {
            node[i]->injectable = 1;
        }
        else
        {
            node[i]->injectable = 0;
        }
        i++;
    }
    delete stmt;
    if(verbose >= 1) { cout << "segments: " << node.size() << endl; }
    
    // Load pairs and rates
    int numberofpairs = 0;
    stmt = db.Prepare("SELECT seg1-1 AS 'segment1', seg2-1 AS 'segment2', rate12"+_carriertype+" AS 'rate', drX, drY, drZ, Jeff2"+_carriertype+", lO"+_carriertype+" FROM pairs UNION SELECT seg2-1 AS 'segment1', seg1-1 AS 'segment2', rate21"+_carriertype+" AS 'rate', -drX AS 'drX', -drY AS 'drY', -drZ AS 'drZ', Jeff2"+_carriertype+", lO"+_carriertype+" FROM pairs ORDER BY segment1;");
    while (stmt->Step() != SQLITE_DONE)
    {
        int seg1 = stmt->Column<int>(0);
        int seg2 = stmt->Column<int>(1);

        double rate12 = stmt->Column<double>(2);

        myvec dr = myvec(stmt->Column<double>(3)*1E-9, stmt->Column<double>(4)*1E-9, stmt->Column<double>(5)*1E-9); // converted from nm to m
        double Jeff2 = stmt->Column<double>(6);
        double reorg_out = stmt->Column<double>(7); 
        node[seg1]->AddEvent(seg2,rate12,dr,Jeff2,reorg_out);
        numberofpairs ++;
    }    
    delete stmt;

    if(verbose >= 1) { cout << "pairs: " << numberofpairs/2 << endl; }
    
    // Calculate initial escape rates !!!THIS SHOULD BE MOVED SO THAT IT'S NOT DONE TWICE IN CASE OF COULOMB INTERACTION!!!
    for(unsigned int i=0; i<node.size(); i++)
    {
        node[i]->InitEscapeRate();
    }
    return node;
}    
    
    
    
}

#ifndef __VOTCA_KMC_IMAGEPOTENTIAL_H_
#define __VOTCA_KMC_IMAGEPOTENTIAL_H_

namespace votca { namespace kmc {
  
using namespace std;

struct s_imagepot {
  s_imagepot(int mx_, int dx_, int dy_, int dz_);
  int mx; // x-coordinate of point "M" for which the potential is computed
  int dx, dy, dz; // coordinates of point "P" relative to point "M"
};

s_imagepot::s_imagepot(int mx_, int dx_, int dy_, int dz_) {
  mx = mx_; dx = dx_; dy = dy_; dz = dz_;
}

class ImagePotential {
public:
  ImagePotential(long N); // N = amount of periodic images we want to use (N=0 is no image charges, N=1 is first two image charges etc.)
  double get_pot(s_imagepot &s);
private:
  double image_potential(s_imagepot &s, long N);
  double precalc[L-1][2*RC+1][2*RC+1][2*RC+1]; // Multi-dimensional array with precalculated image potentials
};

double ImagePotential::get_pot(s_imagepot &s) {
  int smx = s.mx;
  double result = 0;
  if (smx>0 && smx<L) { // Don't bother the caller with this check
    result = precalc[smx-1][s.dx+RC][s.dy+RC][s.dz+RC];
  }
  return result;
}

ImagePotential::ImagePotential(long N) {
  // Precalculate the image potential for all relevant configurations
  for (int il=1; il<L; il++) {
    for (int ix=-RC; ix<=RC; ix++) {
      for (int iy=-RC; iy<=RC; iy++) {
        for (int iz=-RC; iz<=RC; iz++) {
          s_imagepot s(il,ix,iy,iz);
          precalc[il-1][ix+RC][iy+RC][iz+RC] = image_potential(s,N);
        }
      }
    }
  }
}

double ImagePotential::image_potential(s_imagepot &s,long N) {
  long distsyz = s.dy*s.dy + s.dz*s.dz;
  long dists = s.dx*s.dx + distsyz;
  double result = 0;
  if (dists!=0) { // dists==0 is special case: charge interacts with its own image charges
    result = 1.0/sqrt(double(dists));
  }
	for (int i=0; i<100; i++) {
    long dist1;
		long dist2;
		long sign;
    if (div(i,2).rem==0) { // even generation (x-position of image charges is -p.x + 2*j*L, j=...,-1,0,1,...)
			sign = -1;
      dist1 = i*L + 2*s.mx + s.dx;
      dist2 = i*L + 2*L - 2*s.mx - s.dx;
		}
		else { // odd generation (x-position of image charges is -p.x + 2*j*L, j=...,-1,0,1,...)
			sign = 1; 
      dist1 = (i+1)*L + s.dx;
			dist2 = (i+1)*L - s.dx;
		}
		long dists1 = dist1*dist1 + distsyz;
		long dists2 = dist2*dist2 + distsyz;
    result += sign*(1.0/sqrt(double(dists1)) + 1.0/sqrt(double(dists2)));
  }
  return result;
}
