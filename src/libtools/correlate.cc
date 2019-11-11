/*
 * Copyright 2009-2019 The VOTCA Development Team (http://www.votca.org)
 *
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

#include <votca/tools/correlate.h>
#include <votca/tools/eigen.h>

namespace votca {
namespace tools {

void Correlate::CalcCorrelations(DataCollection<double>::selection &data) {
  Index N = Index(data[0].size());
  double Nd = (double)N;
  Eigen::Map<Eigen::ArrayXd> m0(data[0].data(), N);
  double xm = m0.sum();
  xm /= Nd;
  double xsq = m0.abs2().sum();

  for (Index v = 1; v < data.size(); v++) {
    Eigen::Map<Eigen::ArrayXd> m_v(data[v].data(), N);
    double ym = m_v.sum();
    double ysq = m_v.abs2().sum();
    double p = (m_v * m0).sum();
    ym /= Nd;
    double norm = std::sqrt((xsq - Nd * xm * xm) * (ysq - Nd * ym * ym));
    p = (p - Nd * xm * ym) / norm;
    _corr.push_back(p);
  }
}

}  // namespace tools
}  // namespace votca
