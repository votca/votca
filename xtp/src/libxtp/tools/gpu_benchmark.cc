/*
 * Copyright 2009-2021 The VOTCA Development Team (http://www.votca.org)
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

// Standard includes
#include <chrono>
#include <iomanip>
#include <numeric>
#include <string>

#include "votca/tools/property.h"
#include "votca/xtp/bse_operator.h"
#include "votca/xtp/eigen.h"
#include "votca/xtp/logger.h"
#include "votca/xtp/openmp_cuda.h"
#include "votca/xtp/orbitals.h"
#include "votca/xtp/rpa.h"
#include "votca/xtp/threecenter.h"

// Local private VOTCA includes
#include "gpu_benchmark.h"

namespace votca {
namespace xtp {

void GPUBenchmark::ParseOptions(const tools::Property& options) {

  repetitions_ = options.get("repetitions").as<Index>();
}

std::pair<double, double> CalcStatistics(const std::vector<double>& v) {
  double mean = std::reduce(v.begin(), v.end()) / ((double)v.size());
  double sq_sum = std::inner_product(v.begin(), v.end(), v.begin(), 0.0);
  double stdev = std::sqrt(sq_sum / ((double)v.size()) - mean * mean);
  return {mean, stdev};
}

template <class T>
tools::Property RunPart(T&& payload, const std::string& name,
                        Index repetitions) {
  std::vector<double> individual_timings;
  individual_timings.reserve(repetitions);
  Index count = 0;
  std::cout << name << std::endl;
  for (Index i = 0; i < repetitions; i++) {
    std::chrono::time_point<std::chrono::steady_clock> start =
        std::chrono::steady_clock::now();
    count += payload();
    std::chrono::time_point<std::chrono::steady_clock> end =
        std::chrono::steady_clock::now();
    individual_timings.push_back(
        std::chrono::duration_cast<std::chrono::duration<double> >(end - start)
            .count());
  }
  // this is more to prevent optimisation of the whole loop if the payload is
  // small
  if (count != repetitions) {
    std::cout << "sth went wrong with the count" << std::endl;
  }
  auto [mean1, std1] = CalcStatistics(individual_timings);
  std::cout << "avg:" << mean1 << " std:" << std1 << std::endl;
  tools::Property output(name, "", "");
  output.add("avg", std::to_string(mean1));
  output.add("std", std::to_string(std1));
  tools::Property& runs = output.add("runs", "");
  for (double time : individual_timings) {
    runs.add("timing", std::to_string(time));
  }
  return output;
}

bool GPUBenchmark::Run() {

  Orbitals orb;
  orb.ReadFromCpt(job_name_ + ".orb");
  std::cout << "\n\nCreating benchmark for " << job_name_ << ".orb"
            << std::endl;
  AOBasis basis = orb.SetupDftBasis();
  AOBasis auxbasis = orb.SetupAuxBasis();

  std::cout << "Repetitions:" << repetitions_ << std::endl;
  std::cout << "Number of CPUs:" << OPENMP::getMaxThreads()
            << " \nNumber gpus:" << OpenMP_CUDA::UsingGPUs() << std::endl;
  std::cout << "Using MKL for Eigen overload:" << XTP_HAS_MKL_OVERLOAD()
            << std::endl;

  TCMatrix_gwbse Mmn;
  Index max_3c = std::max(orb.getBSEcmax(), orb.getGWAmax());
  Mmn.Initialize(auxbasis.AOBasisSize(), orb.getRPAmin(), max_3c,
                 orb.getRPAmin(), orb.getRPAmax());
  std::cout << "BasisSet:" << basis.Name() << " size:" << basis.AOBasisSize()
            << std::endl;
  std::cout << "AuxBasisSet:" << auxbasis.Name()
            << " size:" << auxbasis.AOBasisSize() << std::endl;
  std::cout << "rpamin:" << orb.getRPAmin() << " rpamax:" << orb.getRPAmax()
            << std::endl;

  tools::Property output("GPU_Benchmark", "", "");
  output.add("Repetitions", std::to_string(repetitions_));
  output.add("CPUs", std::to_string(OPENMP::getMaxThreads()));

  output.add("GPUs", std::to_string(OpenMP_CUDA::UsingGPUs()));
  output.add("MKL_overload", std::to_string(XTP_HAS_MKL_OVERLOAD()));
  output.add("Basisset", basis.Name());
  output.add("Basissetsize", std::to_string(basis.AOBasisSize()));
  output.add("AuxBasisset", auxbasis.Name());
  output.add("AuxBasissetsize", std::to_string(auxbasis.AOBasisSize()));
  output.add("rpamin", std::to_string(orb.getRPAmin()));
  output.add("rpamax", std::to_string(orb.getRPAmax()));
  output.add(RunPart(
      [&]() {
        Mmn.Fill(auxbasis, basis, orb.MOs().eigenvectors());
        return 1;
      },
      "Filling Three Center", repetitions_));

  Eigen::MatrixXd op =
      Eigen::MatrixXd::Identity(auxbasis.AOBasisSize(), auxbasis.AOBasisSize());
  RunPart(
      [&]() {
        Mmn.MultiplyRightWithAuxMatrix(op);
        return 1;
      },
      "Multiplication of tensor with matrix", repetitions_);

  Logger log;
  RPA rpa(log, Mmn);
  rpa.configure(orb.getHomo(), orb.getRPAmin(), orb.getRPAmax());
  rpa.setRPAInputEnergies(orb.MOs().eigenvalues().segment(
      orb.getRPAmin(), orb.getRPAmax() - orb.getRPAmax() + 1));
  double frequency = 0.5;
  Eigen::MatrixXd result =
      Eigen::MatrixXd::Zero(auxbasis.AOBasisSize(), auxbasis.AOBasisSize());
  output.add(RunPart(
      [&]() {
        result += rpa.calculate_epsilon_i(frequency);
        result += rpa.calculate_epsilon_r(-frequency);
        return 1;
      },
      "RPA evaluation", repetitions_));
  frequency += result.sum();
  Index hqp_size = orb.getBSEcmax() - orb.getBSEvmin() + 1;
  Eigen::MatrixXd Hqp_fake = Eigen::MatrixXd::Random(hqp_size, hqp_size);
  BSEOperator_Options opt;
  opt.cmax = orb.getBSEcmax();
  opt.homo = orb.getHomo();
  opt.qpmin = orb.getGWAmin();
  opt.rpamin = orb.getRPAmin();
  opt.vmin = orb.getBSEvmin();
  SingletOperator_TDA s_op(result, Mmn, Hqp_fake);
  s_op.configure(opt);
  TripletOperator_TDA t_op(result, Mmn, Hqp_fake);
  t_op.configure(opt);
  SingletOperator_BTDA_B sbtda_op(result, Mmn, Hqp_fake);
  sbtda_op.configure(opt);
  HxOperator hx_op(result, Mmn, Hqp_fake);
  hx_op.configure(opt);
  Index spacesize = 100;  // A searchspace of 100 is pretty decent
  Eigen::MatrixXd state = Eigen::MatrixXd::Random(s_op.size(), spacesize);
  Eigen::MatrixXd result_op = Eigen::MatrixXd::Zero(s_op.size(), spacesize);

  output.add(RunPart(
      [&]() {
        result_op += s_op * state;
        return 1;
      },
      "SingletOperator_TDA", repetitions_));

  output.add(RunPart(
      [&]() {
        result_op += t_op * state;
        return 1;
      },
      "TripletOperator_TDA", repetitions_));

  output.add(RunPart(
      [&]() {
        result_op += sbtda_op * state;
        return 1;
      },
      "SingletOperator_BTDA_B", repetitions_));

  output.add(RunPart(
      [&]() {
        result_op += hx_op * state;
        return 1;
      },
      "HxOperator", repetitions_));

  std::cout << output << std::endl;
  frequency += result_op.sum();
  return true;
}

}  // namespace xtp
}  // namespace votca
