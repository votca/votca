/*
 * Copyright 2009-2011 The VOTCA Development Team (http://www.votca.org)
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


#include "eoutersphere2.h"


namespace votca { namespace ctp {

void Eoutersphere2::Initialize(Topology *top, Property *opt) {

    string key = "options.eoutersphere.method";

    if (opt->exists(key+".method")) {
        _method = opt->get(key+".method").as<string> ();
    }
    else {
        cout << "ERROR: No method for calculation of outersph. reorg. energy "
                "specified.";
        throw std::runtime_error("Missing input in options file.");
    }
    if (_method == "constant") {
        _lambdaConstant = opt->get(key+".lambdaconst").as<double> ();
    }
    else if (_method == "spheres") {
        _pekarFactor = opt->get(key+".pekar").as<double> ();
    }
    else if (_method == "dielectric") {
        _pekarFactor = opt->get(key+".pekar").as<double> ();
        _lambdaCutoff = opt->get(key+".cutoff").as<double> ();
    }
    else {
        cout << "ERROR: Unidentified method for outer-reorg. energy calc.: "
             << _method;
        throw std::runtime_error("Unrecognized option in options file.");
    }
}

bool Eoutersphere2::EvaluateFrame(Topology *top) {

    if (_method == "constant") { this->ConstLambda(top); }
    else if (_method == "spheres") { this->SpheresLambda(top);}
    else if (_method == "dielectric") { this->DielectricLambda(top); }
}

void Eoutersphere2::ConstLambda(Topology *top) {

    QMNBList2 ::iterator pit;
    for (pit = top->NBList().begin(); pit != top->NBList().end(); pit++) {
        QMPair2 *pair = *pit;
        pair->setLambdaO(_lambdaConstant);
    }
}

void Eoutersphere2::SpheresLambda(Topology *top) {
    
    double e = 1.602176487e-19;
    double EPS0 = 8.85418782e-12;
    double NM = 1e-09;    
    
    QMNBList2 ::iterator pit;
    for (pit = top->NBList().begin(); pit != top->NBList().end(); pit++) {
        QMPair2 *pair = *pit;

        // TODO extract radii from input
        throw std::runtime_error(" EOutersphere2 -> Need to read in radii");
        double R1 = 1;
        double R2 = 2;
        double lambda = _pekarFactor * e / (4.*M_PI*EPS0) *
                                     (  1. / ( 2.*R1*NM )
                                      + 1. / ( 2.*R2*NM )
                                      - 1. / ( pair->Dist()*NM ));
        pair->setLambdaO(lambda);
    }
}

void Eoutersphere2::DielectricLambda(Topology *top) {

    double e = 1.602176487e-19;
    double EPS0 = 8.85418782e-12;
    double NM = 1.e-09;
    double NM3 = 1.e-27;

    QMNBList2 ::iterator pit;
    for (pit = top->NBList().begin(); pit != top->NBList().end(); pit++) {

        QMPair2 *pair = *pit;
        double lambda = 0.;
        vec seg1pos = pair->Seg1()->getPos();
        vec seg2pos = pair->Seg2()->getPos();

        vector< Segment* > ::iterator sit;
        for (sit = top->Segments().begin();
             sit < top->Segments().end();
             sit++) {
             Segment *seg = *sit;

             // Segment external to pair? Cut-off obeyed?
             if ( seg == pair->Seg1() || seg == pair->Seg2() ) { continue; }

             vec dR1bc = top->PbShortestConnect( seg->getPos(), seg1pos );
             vec dR2bc = top->PbShortestConnect( seg->getPos(), seg2pos );
             if ( abs( (dR1bc+dR2bc)*0.5 ) > _lambdaCutoff ) { continue; }

             vec dR1   = seg1pos - seg->getPos();
             vec dR2   = seg2pos - seg->getPos();
             vec shift1 = dR1bc - dR1;
             vec shift2 = dR2bc - dR2;

             vector< Atom* > ::iterator ait;
             vector< Atom* > ::iterator bit;
             for (ait = seg->Atoms().begin(); ait < seg->Atoms().end(); ait++) {

                 Atom *Ext = *ait;
                 // Electric induction field
                 vec D = vec(0.,0.,0.);

                 for (bit = pair->Seg1()->Atoms().begin();
                         bit < pair->Seg1()->Atoms().end();
                         bit++) {

                     Atom *Int = *bit;
                     
                     double dQ = Int->getQ(-1) - Int->getQ(0);
                     
                     vec R = Ext->getPos() - Int->getPos() - shift1;
                     double dr = abs(R);
                     double dr3 = dr*dr*dr;

                     D += e * dQ / (4.*M_PI*dr3*NM3) * R * NM;

                 }

                 for (bit = pair->Seg2()->Atoms().begin();
                         bit < pair->Seg2()->Atoms().end();
                         bit++) {
                     Atom *Int = *bit;

                     double dQ = Int->getQ(0) - Int->getQ(-1);
                     
                     vec R = Ext->getPos() - Int->getPos() - shift2;
                     double dr = abs(R);
                     double dr3 = dr*dr*dr;

                     D += e * dQ / (4.*M_PI*dr3*NM3) * R * NM;

                 }
                 
                 lambda += 1/e * _pekarFactor/(2*EPS0)
                         * top->BoxVolume()*NM3 / top->Atoms().size();   
             } /* exit loop over external atoms in segment */
        } /* exit loop over external segments */

        pair->setLambdaO(lambda);

    } /* exit loop over pairs */

}

}}



