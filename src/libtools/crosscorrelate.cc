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

#include <votca/tools/crosscorrelate.h>
#include <votca/tools/votca_config.h>

#ifndef NOFFTW
#include <fftw3.h>
#endif

namespace votca {
namespace tools {

/**
    \todo clean implementation!!!
*/

#ifdef NOFFTW
void CrossCorrelate::AutoCorrelate(DataCollection<double>::selection&) {

  throw std::runtime_error(
      "CrossCorrelate::AutoCorrelate is not compiled-in due to disabling of "
      "FFTW -recompile Votca Tools with FFTW3 support ");
}
#else
void CrossCorrelate::AutoCorrelate(DataCollection<double>::selection& data) {
  size_t N = data[0].size();
  if (N > (size_t)std::numeric_limits<int>::max()) {
    throw std::runtime_error("CrossCorrelate::AutoCorrelate: size is too big");
  }

  _corrfunc.resize(N);

  fftw_complex* tmp;
  fftw_plan fft, ifft;

  tmp = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (N / 2 + 1));

  fft = fftw_plan_dft_r2c_1d((int)N, &data[0][0], tmp, FFTW_ESTIMATE);
  ifft = fftw_plan_dft_c2r_1d((int)N, tmp, &_corrfunc[0], FFTW_ESTIMATE);
  fftw_execute(fft);

  tmp[0][0] = tmp[0][1] = 0;
  for (size_t i = 1; i < N / 2 + 1; i++) {
    tmp[i][0] = tmp[i][0] * tmp[i][0] + tmp[i][1] * tmp[i][1];
    tmp[i][1] = 0;
  }
  fftw_execute(ifft);

  double d = _corrfunc[0];
  for (size_t i = 0; i < N; i++) {
    _corrfunc[i] = _corrfunc[i] / d;
  }
  fftw_destroy_plan(fft);
  fftw_destroy_plan(ifft);
  fftw_free(tmp);
}
#endif

}  // namespace tools
}  // namespace votca
