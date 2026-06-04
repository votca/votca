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

#include <votca/xtp/extended_hueckel.h>

namespace votca {
namespace xtp {

namespace {
constexpr double ev_to_ha = 0.03674932217565499;
}

ExtendedHuckelParameters::ExtendedHuckelParameters() {
  auto add = [&](const std::string& el, int l, double ev) {
    eps_[{el, l}] = ev * ev_to_ha;
  };

  // H
  add("H", 0, -13.60);

  // Second row
  add("B", 0, -15.20);
  add("B", 1, -8.50);

  add("C", 0, -21.40);
  add("C", 1, -11.40);

  add("N", 0, -26.00);
  add("N", 1, -13.40);

  add("O", 0, -32.30);
  add("O", 1, -14.80);

  add("F", 0, -40.00);
  add("F", 1, -18.10);

  // Third row
  add("Al", 0, -12.30);
  add("Al", 1, -6.50);

  add("Si", 0, -17.30);
  add("Si", 1, -9.20);

  add("P", 0, -18.60);
  add("P", 1, -10.70);
  add("P", 2, -14.00);

  add("S", 0, -20.00);
  add("S", 1, -11.00);
  add("S", 2, -14.20);

  add("Cl", 0, -30.00);
  add("Cl", 1, -15.00);
  add("Cl", 2, -18.00);

  // Halogens often relevant in XTP
  add("Br", 0, -28.00);
  add("Br", 1, -14.00);
  add("Br", 2, -16.00);

  add("I", 0, -25.00);
  add("I", 1, -12.00);
  add("I", 2, -14.00);

  // Optional: common heavier main-group
  add("Ga", 0, -12.10);
  add("Ga", 1, -6.70);

  add("Ge", 0, -16.00);
  add("Ge", 1, -8.90);

  add("As", 0, -18.90);
  add("As", 1, -11.00);
  add("As", 2, -14.50);

  add("Se", 0, -20.80);
  add("Se", 1, -12.00);
  add("Se", 2, -14.80);
}

bool ExtendedHuckelParameters::Has(const std::string& element, int l) const {
  return eps_.find({element, l}) != eps_.end();
}

bool ExtendedHuckelParameters::HasElement(const std::string& element) const {
  for (const auto& kv : eps_) {
    if (kv.first.first == element) {
      return true;
    }
  }
  return false;
}

double ExtendedHuckelParameters::Get(const std::string& element, int l) const {
  auto it = eps_.find({element, l});
  if (it == eps_.end()) {
    throw std::runtime_error("No exact EHT parameter for element " + element +
                             " and l=" + std::to_string(l));
  }
  return it->second;
}

double ExtendedHuckelParameters::GetWithFallback(const std::string& element,
                                                 int l, int* used_l) const {
  // 1. Exact match
  auto exact = eps_.find({element, l});
  if (exact != eps_.end()) {
    if (used_l != nullptr) {
      *used_l = l;
    }
    return exact->second;
  }

  // 2. Element missing entirely
  if (!HasElement(element)) {
    throw std::runtime_error("No EHT parameters available for element " +
                             element);
  }

  // 3. Fallback logic
  // For higher-l shells: prefer d, then p, then s
  if (l >= 2) {
    if (Has(element, 2)) {
      if (used_l != nullptr) {
        *used_l = 2;
      }
      return Get(element, 2);
    }
    if (Has(element, 1)) {
      if (used_l != nullptr) {
        *used_l = 1;
      }
      return Get(element, 1);
    }
    if (Has(element, 0)) {
      if (used_l != nullptr) {
        *used_l = 0;
      }
      return Get(element, 0);
    }
  }

  // For p shells: fall back to s
  if (l == 1) {
    if (Has(element, 0)) {
      if (used_l != nullptr) {
        *used_l = 0;
      }
      return Get(element, 0);
    }
  }

  // For completeness: if some weird lower-l case comes in
  // choose the highest available shell on that element
  for (int ll = 6; ll >= 0; --ll) {
    if (Has(element, ll)) {
      if (used_l != nullptr) {
        *used_l = ll;
      }
      return Get(element, ll);
    }
  }

  throw std::runtime_error("No usable EHT fallback parameter for element " +
                           element + " and l=" + std::to_string(l));
}

}  // namespace xtp
}  // namespace votca