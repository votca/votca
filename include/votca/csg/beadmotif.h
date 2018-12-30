/*
 * Copyright 2009-2018 The VOTCA Development Team (http://www.votca.org)
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

#ifndef _VOTCA_CSG_BEADMOTIF_H
#define _VOTCA_CSG_BEADMOTIF_H

#include <list>
#include <map>
#include <memory>
#include <set>
#include <string>

#include <votca/csg/basebead.h>

#include <votca/tools/graph.h>

namespace votca {
  namespace csg {

    /**
     * \brief Designed determine what kind of structure a beadstructure has
     *
     * Wants a beadstructure is created there are several different forms it can
     * have, this class helps to classify and break structures up into the 
     * appropriate sub class. The possible classes include:
     *
     * 1. Single 
     * 2. Line
     * 3. Loop
     * 4. Fused Ring
     * 5. Other
     * 6. Undefined
     *
     * The Single, line, loop and Fused Ring types are all elementary types that 
     * represent a fundamental structural unit. 
     *
     * Other represents a type that has not been broken up into its fundamental
     * component.
     * Undefined means the structure has not yet been categorized. 
     *
     * Though the Beadmotif inherits from the Beadstructure its methods are kept
     * private the reason is that the type might be incorrect if an  AddBead
     * from the original class is called. 
     **/

    class BeadMotif : private BeadStructure{
      public:
        BeadMotif() : type_(motif_type::undefined) {};
        ~BeadMotif() {}

        enum motif_type {
          single;
          line;
          loop;
          fused_ring;
          other;
          undefined;
        };

        motif_type getType() const;

        /**
         * \brief Calculates the type of the motif
         **/
        void CalculateType();
      private:

        motif_type type_;
        bool isSingle_();
        bool isLoop_();
        bool isFusedRing_();
        bool isLine_();

    };
  }
}

#endif // _VOTCA_CSG_BEADMOTIF_H
