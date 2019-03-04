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

/**

  \mainpage VOTCA C++ reference

  \section intro_sec Introduction

  This page is the C++ code documentation of the VOTCA package
  (http://www.votca.org). The main target of VOTCA is the implementation of
  systematic coarse-graining techniques. However, it offers a powerful,
  object-oriented framework to develop analysis tools for particle based
  molecular simulations.

  \section started_sec Getting started
  To start developing custom analysis tools, a good place to start is the
  csgapps repository:

  https://github.com/votca/csgapps

  It contains several small analysis tools which were implemented based upon the
  VOTCA framework. We highly recomment to use an IDE such as Netbeans for
  development since it offers lots of guides to get started with new code (code
  completion, code documentation popups, navigation thourh code, ...).

  The main container for the whole structure is the Topology, so it is a good
  advise to get comfortable with this class. Also the standard applications in
  csg/src/tools might help.

  \section beginner_sec For beginners: how to avoid frustration

  For those not familiar with object oriented code: don't try to dig into every
  single function in order to understand what exactly is going on. This strategy
  only works for very small projects and is not intended for oject oriented
  programs. Think about the code in layers of abstraction! Your main focus
  should be on the global structure and understand how objects relate to each
  other. The code was designed that you don't have to redo and understand all
  the nasty details!

 */

#ifndef _VOTCA_CSG_VERSION_H
#define _VOTCA_CSG_VERSION_H

#include <string>

namespace votca {
namespace csg {
const std::string &CsgVersionStr();
void               HelpTextHeader(const std::string &tool_name);
}  // namespace csg
}  // namespace votca

#endif /* _VOTCA_CSG_VERSION_H */
