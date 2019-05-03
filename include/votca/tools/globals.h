/*
 * Copyright 2009-2019 The VOTCA Development Team (http://www.votca.org)
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

#pragma once
#ifndef __VOTCA_TOOLS_GLOBALS_H
#define __VOTCA_TOOLS_GLOBALS_H

#include <string>
#include <votca/tools/votca_config.h>

namespace votca {
namespace tools {

/**
    \brief class to store global variables

    This class is used to access global variables
*/

struct globals {
  /// be loud and noisy
  static bool verbose;
  /// web of the package
  static std::string url;
  /// email address of the developers
  static std::string email;

  /// If Eigen is overloaded with MKL
  static bool VOTCA_MKL;

  /// man pages format strings
  struct man {
    static std::string option;
    static std::string header;
    static std::string name;
    static std::string authors;
    static std::string copyright;
    static std::string synopsis;
    static std::string description;
    static std::string options;
  };

  /// TEX pages format strings
  struct tex {
    static std::string section;
    static std::string label;
    static std::string description;
    static std::string options;
    static std::string option;
  };

  // constants
  struct constants {
    static const double pi;
    static const double kB;    // eV/K
    static const double hbar;  // eV*s
  };

  // conversion factors
  struct conversion {
    static const double Bohr2nm;
    static const double nm2Bohr;
    static const double Ang2Bohr;
    static const double Ryd2eV;
    static const double Hrt2eV;
    static const double int2eV;  // ewald internal to eV conversion
  };
};

}  // namespace tools
}  // namespace votca

#endif /* __VOTCA_TOOLS_GLOBALS_H */
