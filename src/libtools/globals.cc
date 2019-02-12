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

#include <votca/tools/globals.h>

namespace votca {
namespace tools {

bool globals::verbose = false;
std::string globals::url = "http://www.votca.org";
std::string globals::email = "devs@votca.org";

#if defined(MKL)
bool globals::VOTCA_MKL = true;
#else
bool globals::VOTCA_MKL = false;
#endif

std::string globals::man::option(".TP\n\\fB%1%\\fR\n%2%\n");

std::string globals::man::header(".TH \"%1%\" 1 \"\" \"Version: %2%\"\n\n");

std::string globals::man::name(
    ".SH NAME\n"
    "\n.P\n"
    "%1% \\- Part of the VOTCA package\n"
    "\n.P\n"
    "For more info please visit %2%\n\n");

std::string globals::man::authors(
    "\n.SH AUTHORS\n"
    "\n.P\n"
    "Written and maintained by the VOTCA Development Team <%1%>\n");

std::string globals::man::copyright(
    "\n.SH COPYRIGHT\n"
    "\n.P\n\n"
    "Copyright 2009\\-2019 The VOTCA Development Team (%1%).\n"
    "\n.P\n"
    "Licensed under the Apache License, Version 2.0 (the \"License\") "
    "you may not use this file except in compliance with the License. "
    "You may obtain a copy of the License at"
    "\n.P\n"
    "http://www.apache.org/licenses/LICENSE\\-2.0\n"
    "\n.P\n"
    "Unless required by applicable law or agreed to in writing, software "
    "distributed under the License is distributed on an \"AS IS\" BASIS, "
    "WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. "
    "See the License for the specific language governing permissions and "
    "limitations under the License.");

std::string globals::man::synopsis(
    "\n.SH SYNOPSIS\n"
    "\n.P\n\\fB%1%\\fR [\\fIOPTION\\fR] [\\fIARGUMENT\\fR]\n");

std::string globals::man::description(
    "\n.SH DESCRIPTION\n"
    "\n.P\n%1%\n");

std::string globals::man::options("\n.SH OPTIONS\n");

std::string globals::tex::section("\n\\subsection{%1%}");
std::string globals::tex::label("\n\\label{prog:%1%}");
std::string globals::tex::description("\n%1%");
std::string globals::tex::option("\n\\item[ ] \\texttt{%1%} %2%");
std::string globals::tex::options(
    "\n\\begin{compactitem}%1%\n\\end{compactitem}\n");

// constants
const double globals::constants::pi = 3.1415926535897;
const double globals::constants::kB = 8.617332478E-5;      // eV/K
const double globals::constants::hbar = 6.5821192815E-16;  // eV*s

// conversion factors
const double globals::conversion::Bohr2nm = 0.052917721092;
const double globals::conversion::nm2Bohr = 18.897259886;
const double globals::conversion::Ang2Bohr = 1.8897259886;
const double globals::conversion::Ryd2eV = 13.6058;
const double globals::conversion::Hrt2eV = 27.21138386;
const double globals::conversion::int2eV =
    1 / (4 * constants::pi * 8.854187817e-12) * 1.602176487e-19 / 1.000e-9;

}  // namespace tools
}  // namespace votca
