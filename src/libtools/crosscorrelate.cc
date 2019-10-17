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
void CrossCorrelate::AutoCorrelate(DataCollection<double>::selection* data,
                                   bool average) {
#ifdef NOFFTW
  throw std::runtime_error(
      "CrossCorrelate::AutoCorrelate is not compiled-in due to disabling of "
      "FFTW -recompile Votca Tools with FFTW3 support ");
#else
  size_t N = (*data)[0].size();
  _corrfunc.resize(N);

  fftw_complex* tmp;
  fftw_plan fft, ifft;

  tmp = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (N / 2 + 1));

  fft = fftw_plan_dft_r2c_1d(N, &(*data)[0][0], tmp, FFTW_ESTIMATE);
  ifft = fftw_plan_dft_c2r_1d(N, tmp, &_corrfunc[0], FFTW_ESTIMATE);
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
#endif
}

void CrossCorrelate::AutoFourier(std::vector<double>& ivec) {
#ifdef NOFFTW
  throw std::runtime_error(
      "CrossCorrelate::AutoFourier is not compiled-in due to disabling of FFTW "
      "-recompile Votca Tools with FFTW3 support ");
#else
  size_t N = ivec.size();
  _corrfunc.resize(N);

  fftw_complex* tmp;
  fftw_plan fft;

  tmp = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (N / 2 + 1));

  fft = fftw_plan_dft_r2c_1d(N, &ivec[0], tmp, FFTW_ESTIMATE);
  fftw_execute(fft);

  tmp[0][0] = tmp[0][1] = 0;
  for (size_t i = 1; i < N / 2 + 1; i++) {
    tmp[i][0] = tmp[i][0] * tmp[i][0] + tmp[i][1] * tmp[i][1];
    tmp[i][1] = 0;
  }

  // copy the real component of temp to the _corrfunc vector
  for (size_t i = 0; i < N; i++) {
    _corrfunc[i] = tmp[i][0];
  }

  fftw_destroy_plan(fft);
  fftw_free(tmp);
#endif
}

void CrossCorrelate::FFTOnly(std::vector<double>& ivec) {
#ifdef NOFFTW
  throw std::runtime_error(
      "CrossCorrelate::FFTOnly is not compiled-in due to disabling of FFTW "
      "-recompile Votca Tools with FFTW3 support ");
#else
  size_t N = ivec.size();
  _corrfunc.resize(N);

  fftw_complex* tmp;
  fftw_plan fft;

  tmp = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (N / 2 + 1));

  fft = fftw_plan_dft_r2c_1d(N, &ivec[0], tmp, FFTW_ESTIMATE);
  fftw_execute(fft);

  // copy the real component of temp to the _corrfunc vector
  for (size_t i = 0; i < N / 2 + 1; i++) {
    _corrfunc[i] = tmp[i][0];
  }

  fftw_destroy_plan(fft);
  fftw_free(tmp);
#endif
}

void CrossCorrelate::DCTOnly(std::vector<double>& ivec) {
#ifdef NOFFTW
  throw std::runtime_error(
      "CrossCorrelate::DCTOnly is not compiled-in due to disabling of FFTW "
      "-recompile Votca Tools with FFTW3 support ");
#else
  size_t N = ivec.size();
  _corrfunc.resize(N);

  std::vector<double> tmp(N);

  fftw_plan fft;

  // do real to real discrete cosine trafo
  fft = fftw_plan_r2r_1d(N, &ivec[0], &tmp[0], FFTW_REDFT10, FFTW_ESTIMATE);
  fftw_execute(fft);
  fftw_destroy_plan(fft);
  _corrfunc = tmp;

#endif
}

void CrossCorrelate::AutoCosine(std::vector<double>& ivec) {
#ifdef NOFFTW
  throw std::runtime_error(
      "CrossCorrelate::AutoCosine is not compiled-in due to disabling of FFTW "
      "-recompile Votca Tools with FFTW3 support ");
#else
  size_t N = ivec.size();
  _corrfunc.resize(N);

  std::vector<double> tmp(N);
  fftw_plan fft;

  // do real to real discrete cosine trafo
  fft = fftw_plan_r2r_1d(N, &ivec[0], &tmp[0], FFTW_REDFT10, FFTW_ESTIMATE);
  fftw_execute(fft);
  fftw_destroy_plan(fft);
  // compute autocorrelation
  tmp[0] = 0;
  for (size_t i = 1; i < N; i++) {
    tmp[i] = tmp[i] * tmp[i];
  }

  _corrfunc = tmp;

#endif
}

void CrossCorrelate::AutoCorr(std::vector<double>& ivec) {
#ifdef NOFFTW
  throw std::runtime_error(
      "CrossCorrelate::AutoCorr is not compiled-in due to disabling of FFTW "
      "-recompile Votca Tools with FFTW3 support ");
#else
  size_t N = ivec.size();
  _corrfunc.resize(N);

  fftw_complex* tmp;
  fftw_plan fft, ifft;

  tmp = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (N / 2 + 1));

  fft = fftw_plan_dft_r2c_1d(N, &ivec[0], tmp, FFTW_ESTIMATE);
  ifft = fftw_plan_dft_c2r_1d(N, tmp, &_corrfunc[0], FFTW_ESTIMATE);
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
#endif
}

}  // namespace tools
}  // namespace votca
