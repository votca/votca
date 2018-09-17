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


#ifndef _VOTCA_XTP_EINTERNAL_H
#define _VOTCA_XTP_EINTERNAL_H

#include <votca/ctp/qmcalculator.h>

namespace votca { namespace xtp {


class EInternal : public ctp::QMCalculator
{
public:

    EInternal() { };
   ~EInternal() { };

    std::string Identify() { return "einternal"; }
    void Initialize(tools::Property *options);
    void ParseEnergiesXML(tools::Property *options);
    bool EvaluateFrame(ctp::Topology *top);

private:

    std::map<std::string, double> _seg_U_cC_nN_e;
    std::map<std::string, double> _seg_U_nC_nN_e;
    std::map<std::string, double> _seg_U_cN_cC_e;

    std::map<std::string, double> _seg_U_cC_nN_h;
    std::map<std::string, double> _seg_U_nC_nN_h;
    std::map<std::string, double> _seg_U_cN_cC_h;
    
    std::map<std::string, double> _seg_U_xX_nN_s;
    std::map<std::string, double> _seg_U_nX_nN_s;
    std::map<std::string, double> _seg_U_xN_xX_s;
    
    std::map<std::string, double> _seg_U_xX_nN_t;
    std::map<std::string, double> _seg_U_nX_nN_t;
    std::map<std::string, double> _seg_U_xN_xX_t;
    

    std::map<std::string, bool>   _seg_has_e;
    std::map<std::string, bool>   _seg_has_h;
    std::map<std::string, bool>   _seg_has_s;
    std::map<std::string, bool>   _seg_has_t;

    std::map<std::string, bool>   _has_seg;

};

void EInternal::Initialize(tools::Property *options) {

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

void EInternal::ParseEnergiesXML(tools::Property *opt) {

    // update options with the VOTCASHARE defaults   
    UpdateWithDefaults( opt, "xtp" );
    std::string key = "options." + Identify();

    std::string energiesXML = opt->get(key+".energiesXML").as<std::string> ();

    std::cout << std::endl
         << "... ... Site, reorg. energies from " << energiesXML << ". "
         << std::flush;

    tools::Property alloc;
    tools::load_property_from_xml(alloc, energiesXML.c_str());

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
    std::list<tools::Property*> mols = alloc.Select(key);
    for (tools::Property* molprop:mols) {

        key = "segments.segment";
        std::list<tools::Property*> segs = molprop->Select(key);
        for (tools::Property* segprop:segs) {

            std::string segName =segprop->get("name").as<std::string> ();

            bool has_seg = true;
            bool has_e = false;
            bool has_h = false;
            bool has_s = false;
            bool has_t = false;

            double U_cC_nN_e = 0.0;
            double U_cC_nN_h = 0.0;
            double U_nC_nN_e = 0.0;
            double U_nC_nN_h = 0.0;
            double U_cN_cC_e = 0.0;
            double U_cN_cC_h = 0.0;
            
            double U_xX_nN_s = 0.0;
            double U_xX_nN_t = 0.0;
            double U_nX_nN_s = 0.0;
            double U_nX_nN_t = 0.0;
            double U_xN_xX_s = 0.0;
            double U_xN_xX_t = 0.0;
            

            if ( segprop->exists("U_cC_nN_e") &&
                segprop->exists("U_nC_nN_e") &&
                 segprop->exists("U_cN_cC_e")    ) {

                U_cC_nN_e = segprop->get("U_cC_nN_e").as< double > ();
                U_nC_nN_e = segprop->get("U_nC_nN_e").as< double > ();
                U_cN_cC_e = segprop->get("U_cN_cC_e").as< double > ();

                has_e = true;
            }
            
            if ( segprop->exists("U_cC_nN_h") &&
                 segprop->exists("U_nC_nN_h") &&
                 segprop->exists("U_cN_cC_h")    ) {

                U_cC_nN_h = segprop->get("U_cC_nN_h").as< double > ();
                U_nC_nN_h = segprop->get("U_nC_nN_h").as< double > ();
                U_cN_cC_h = segprop->get("U_cN_cC_h").as< double > ();

                has_h = true;
            }
            
            if ( segprop->exists("U_xX_nN_s") &&
                 segprop->exists("U_nX_nN_s") &&
                 segprop->exists("U_xN_xX_s")    ) {

                U_xX_nN_s = segprop->get("U_xX_nN_s").as< double > ();
                U_nX_nN_s = segprop->get("U_nX_nN_s").as< double > ();
                U_xN_xX_s = segprop->get("U_xN_xX_s").as< double > ();

                has_s = true;
            }
            if ( segprop->exists("U_xX_nN_t") &&
                 segprop->exists("U_nX_nN_t") &&
                 segprop->exists("U_xN_xX_t")    ) {

                U_xX_nN_t = segprop->get("U_xX_nN_t").as< double > ();
                U_nX_nN_t = segprop->get("U_nX_nN_t").as< double > ();
                U_xN_xX_t = segprop->get("U_xN_xX_t").as< double > ();

                has_t = true;
            }
            _seg_U_cC_nN_e[segName] = U_cC_nN_e;
            _seg_U_nC_nN_e[segName] = U_nC_nN_e;
            _seg_U_cN_cC_e[segName] = U_cN_cC_e;
            _seg_has_e[segName] = has_e;

            _seg_U_cC_nN_h[segName] = U_cC_nN_h;
            _seg_U_nC_nN_h[segName] = U_nC_nN_h;
            _seg_U_cN_cC_h[segName] = U_cN_cC_h;
            _seg_has_h[segName] = has_h;
                        
            _seg_U_xX_nN_s[segName] = U_xX_nN_s;
            _seg_U_nX_nN_s[segName] = U_nX_nN_s;
            _seg_U_xN_xX_s[segName] = U_xN_xX_s;
            _seg_has_s[segName] = has_s;
                        
            _seg_U_xX_nN_t[segName] = U_xX_nN_t;
            _seg_U_nX_nN_t[segName] = U_nX_nN_t;
            _seg_U_xN_xX_t[segName] = U_xN_xX_t;
            _seg_has_t[segName] = has_t;
            
            _has_seg[segName] = has_seg;
           
        }
    }
}

bool EInternal::EvaluateFrame(ctp::Topology *top) {

    int count = 0;
    for (ctp::Segment* seg:top->Segments()) {

        std::string segName = seg->getName();
        
        try {
            _has_seg.at(segName);
        }
        catch (const std::exception& out_of_range) {
            std::cout << std::endl << "... ... WARNING: No energy information for seg ["
                         << segName << "]. Skipping... ";
            continue;
        }

        ++count;

        if (_seg_has_e[segName]) {

            double u  = _seg_U_cC_nN_e[segName];
            double l1 = _seg_U_nC_nN_e[segName];
            double l2 = _seg_U_cN_cC_e[segName];
            bool has_e = true;

            seg->setU_cC_nN(u, -1);
            seg->setU_nC_nN(l1, -1);
            seg->setU_cN_cC(l2, -1);
            seg->setHasState(has_e, -1);
        }

        if (_seg_has_h[segName]) {

            double u  = _seg_U_cC_nN_h[segName];
            double l1 = _seg_U_nC_nN_h[segName];
            double l2 = _seg_U_cN_cC_h[segName];
            bool has_h = true;

            seg->setU_cC_nN(u, +1);
            seg->setU_nC_nN(l1, +1);
            seg->setU_cN_cC(l2, +1);
            seg->setHasState(has_h, +1);
        }
         if (_seg_has_s[segName]) {

            double u  = _seg_U_xX_nN_s[segName];
            double l1 = _seg_U_nX_nN_s[segName];
            double l2 = _seg_U_xN_xX_s[segName];
            bool has_s = true;
            
            seg->setU_xX_nN(u, +2);
            seg->setU_nX_nN(l1, +2);
            seg->setU_xN_xX(l2, +2);
            seg->setHasState(has_s, +2);
        }
        if (_seg_has_t[segName]) {

            double u  = _seg_U_xX_nN_t[segName];
            double l1 = _seg_U_nX_nN_t[segName];
            double l2 = _seg_U_xN_xX_t[segName];
            bool has_t = true;

            seg->setU_xX_nN(u, +3);
            seg->setU_nX_nN(l1, +3);
            seg->setU_xN_xX(l2, +3);
            seg->setHasState(has_t, +3);
        }
    }

    std::cout << std::endl
         << "... ... Read in site, reorg. energies for " 
         << count << " segments. " << std::flush;


    return 1;
}


}}

#endif //_VOTCA_XTP_EINTERNAL_H
