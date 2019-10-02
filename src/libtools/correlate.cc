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

#include <math.h>
#include <votca/tools/correlate.h>

namespace votca {
namespace tools {

/**
    \todo clean implementation!!!
*/
void Correlate::CalcCorrelations(DataCollection<double>::selection *data) {
  size_t N;
  double xm(0), xsq(0);

  N = (*data)[0].size();
  for (size_t i = 0; i < N; i++) {
    xm += (*data)[0][i];
    xsq += (*data)[0][i] * (*data)[0][i];
  }
  xm /= (double)N;

  for (size_t v = 1; v < data->size(); v++) {
    std::pair<std::string, double> p("do_names", 0);
    double ym(0), ysq(0);

    for (size_t i = 0; i < N; i++) {
      ym += (*data)[v][i];
      ysq += (*data)[v][i] * (*data)[v][i];
      p.second += (*data)[v][i] * (*data)[0][i];
    }
    ym /= (double)N;
    double norm = (xsq - ((double)N) * xm * xm) * (ysq - ((double)N) * ym * ym);
    p.second = (p.second - ((double)N) * xm * ym) / sqrt(norm);
    _corr.push_back(p);
  }
}

}  // namespace tools
}  // namespace votca
