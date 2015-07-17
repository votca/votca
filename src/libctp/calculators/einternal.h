/*
 *            Copyright 2009-2012 The VOTCA Development Team
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


#ifndef _VOTCA_CTP_EINTERNAL_H
#define _VOTCA_CTP_EINTERNAL_H

#include <votca/ctp/qmcalculator.h>

namespace votca { namespace ctp {


class EInternal : public QMCalculator
{
public:

    EInternal() { };
   ~EInternal() { };

    string Identify() { return "einternal"; }
    void Initialize(Property *options);
    void ParseEnergiesXML(Property *options);
    bool EvaluateFrame(Topology *top);

private:

    map<string, double> _seg_U_cC_nN_e;
    map<string, double> _seg_U_nC_nN_e;
    map<string, double> _seg_U_cN_cC_e;

    map<string, double> _seg_U_cC_nN_h;
    map<string, double> _seg_U_nC_nN_h;
    map<string, double> _seg_U_cN_cC_h;

    map<string, bool>   _seg_has_e;
    map<string, bool>   _seg_has_h;

    map<string, bool>   _has_seg;

};

void EInternal::Initialize(Property *options) {

    /* ---- OPTIONS.XML Structure -----
     *
     * <einternal>
     *
     *      <energiesXML>ENERGIES.XML</energiesXML>
     *
     * </einternal>
     *
     */

    this->ParseEnergiesXML(options);
}

void EInternal::ParseEnergiesXML(Property *opt) {

    // update options with the VOTCASHARE defaults   
    UpdateWithDefaults( opt );
    string key = "options." + Identify();

    string energiesXML = opt->get(key+".energiesXML").as<string> ();

    cout << endl
         << "... ... Site, reorg. energies from " << energiesXML << ". "
         << flush;

    Property alloc;
    load_property_from_xml(alloc, energiesXML.c_str());

    /* --- ENERGIES.XML Structure ---
     *
     * <topology>
     *
     *     <molecules>
     *          <molecule>
     *          <name></name>
     *
     *          <segments>
     *
     *              <segment>
     *              <name></name>
     *
     *              <!-- U_sG_sG, s->state, G->geometry !-->
     *
     *              <U_cC_nN_e></U_cC_nN_e>
     *              <U_cC_nN_h></U_cC_nN_h>
     *
     *              <U_nC_nN_e></U_nC_nN_e>
     *              <U_nC_nN_h></U_nC_nN_h>
     *
     *              <U_cN_cC_e></U_cN_cC_e>
     *              <U_cN_cC_h></U_cN_cC_h>
     *
     *              </segment>
     *
     *              <segment>
     *                  ...
     *
     */

    key = "topology.molecules.molecule";
    list<Property*> mols = alloc.Select(key);
    list<Property*> ::iterator molit;
    for (molit = mols.begin(); molit != mols.end(); ++molit) {

        key = "segments.segment";
        list<Property*> segs = (*molit)->Select(key);
        list<Property*> ::iterator segit;

        for (segit = segs.begin(); segit != segs.end(); ++segit) {

            string segName = (*segit)->get("name").as<string> ();

            bool has_seg = true;
            bool has_e = false;
            bool has_h = false;

            double U_cC_nN_e = 0.0;
            double U_cC_nN_h = 0.0;
            double U_nC_nN_e = 0.0;
            double U_nC_nN_h = 0.0;
            double U_cN_cC_e = 0.0;
            double U_cN_cC_h = 0.0;

            if ( (*segit)->exists("U_cC_nN_e") &&
                 (*segit)->exists("U_nC_nN_e") &&
                 (*segit)->exists("U_cN_cC_e")    ) {

                U_cC_nN_e = (*segit)->get("U_cC_nN_e").as< double > ();
                U_nC_nN_e = (*segit)->get("U_nC_nN_e").as< double > ();
                U_cN_cC_e = (*segit)->get("U_cN_cC_e").as< double > ();

                has_e = true;
            }

            if ( (*segit)->exists("U_cC_nN_h") &&
                 (*segit)->exists("U_nC_nN_h") &&
                 (*segit)->exists("U_cN_cC_h")    ) {

                U_cC_nN_h = (*segit)->get("U_cC_nN_h").as< double > ();
                U_nC_nN_h = (*segit)->get("U_nC_nN_h").as< double > ();
                U_cN_cC_h = (*segit)->get("U_cN_cC_h").as< double > ();

                has_h = true;
            }

            _seg_U_cC_nN_e[segName] = U_cC_nN_e;
            _seg_U_nC_nN_e[segName] = U_nC_nN_e;
            _seg_U_cN_cC_e[segName] = U_cN_cC_e;
            _seg_has_e[segName] = has_e;

            _seg_U_cC_nN_h[segName] = U_cC_nN_h;
            _seg_U_nC_nN_h[segName] = U_nC_nN_h;
            _seg_U_cN_cC_h[segName] = U_cN_cC_h;
            _seg_has_h[segName] = has_h;
            _has_seg[segName] = has_seg;
        }
    }
}

bool EInternal::EvaluateFrame(Topology *top) {

    vector< Segment* > ::iterator sit;
    int count = 0;
    for (sit = top->Segments().begin(); sit < top->Segments().end(); ++sit) {

        string segName = (*sit)->getName();

        try {
            //bool has_seg = _has_seg.at(segName);
        }
        catch (out_of_range) {
            cout << endl << "... ... WARNING: No energy information for seg ["
                         << segName << "]. Skipping... ";
            continue;
        }

        ++count;

        if (_seg_has_e[segName]) {

            double u  = _seg_U_cC_nN_e[segName];
            double l1 = _seg_U_nC_nN_e[segName];
            double l2 = _seg_U_cN_cC_e[segName];
            bool has_e = true;

            (*sit)->setU_cC_nN(u, -1);
            (*sit)->setU_nC_nN(l1, -1);
            (*sit)->setU_cN_cC(l2, -1);
            (*sit)->setHasState(has_e, -1);
        }

        if (_seg_has_h[segName]) {

            double u  = _seg_U_cC_nN_h[segName];
            double l1 = _seg_U_nC_nN_h[segName];
            double l2 = _seg_U_cN_cC_h[segName];
            bool has_h = true;

            (*sit)->setU_cC_nN(u, +1);
            (*sit)->setU_nC_nN(l1, +1);
            (*sit)->setU_cN_cC(l2, +1);
            (*sit)->setHasState(has_h, +1);
        }
    }

    cout << endl
         << "... ... Read in site, reorg. energies for " 
         << count << " segments. " << flush;


    return 1;
}


}}

#endif //_VOTCA_CTP_EINTERNAL_H
