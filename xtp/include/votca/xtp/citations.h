/*
 *            Copyright 2009-2026 The VOTCA Development Team
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

#pragma once
#ifndef VOTCA_XTP_CITATIONS_H
#define VOTCA_XTP_CITATIONS_H

// Standard includes
#include <string>
#include <vector>

namespace votca {
namespace xtp {

/**
 * \brief References for the methods used in a run.
 *
 * The code that runs a method calls Cite with the key of its reference
 * (the keys of the bibliography of the manual). The application prints the
 * references used, in the order of first use, at the end of the run. The
 * table of calculators gives the references every run of a calculator uses
 * and those that depend on its options (printed by --cite).
 */
namespace citations {

struct Reference {
  std::string key;
  std::string text;
  std::string doi;  // may be empty
};

/// a reference that depends on the options of a calculator
struct Optional {
  std::string key;
  std::string when;
};

struct CalculatorEntry {
  std::vector<std::string> keys;  // every run
  std::vector<Optional> optional;
};

/// marks a reference as used (thread safe; unknown keys are printed as such)
void Cite(const std::string& key);
/// the keys used so far, in the order of first use
std::vector<std::string> Used();
void Reset();

/// all references
const std::vector<Reference>& References();
bool Exists(const std::string& key);
/// one line: text and DOI
std::string Format(const std::string& key);

/// the block printed at the end of a run (empty if nothing was cited)
std::string Report();

/// the references of a calculator or tool (empty for unknown names)
CalculatorEntry ForCalculator(const std::string& name);
/// all calculators and tools of the table
std::vector<std::string> Calculators();
/// the text printed by --cite
std::string CalculatorReport(const std::string& name);

}  // namespace citations
}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_CITATIONS_H
