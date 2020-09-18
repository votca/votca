/*
 * Copyright 2009-2020 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
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
#define BOOST_TEST_MAIN

#define BOOST_TEST_MODULE aobasis_test

// Standard includes
#include <fstream>

// Third party includes
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/xtp/aobasis.h"
#include "votca/xtp/aomatrix.h"
#include "votca/xtp/convergenceacc.h"
#include "votca/xtp/orbitals.h"

using namespace votca::xtp;
using namespace votca;

BOOST_AUTO_TEST_SUITE(aobasis_test)

BOOST_AUTO_TEST_CASE(FillNormBasis_test) {

  Orbitals orbitals;
  orbitals.QMAtoms().LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                                  "/aobasis/Al.xyz");
  BasisSet basis;
  basis.Load(std::string(XTP_TEST_DATA_FOLDER) + "/aobasis/notnormalized.xml");
  AOBasis aobasis;
  aobasis.Fill(basis, orbitals.QMAtoms());

  const AOShell& shell = aobasis.getShell(0);
  std::vector<double> ref_results = {0.1831079647, 0.9155398233};
  Index i = 0;
  bool check_norm = true;
  for (const AOGaussianPrimitive& gaussian : shell) {
    if (std::abs(ref_results[i] - gaussian.getContraction()) > 1e-7) {
      check_norm = false;
      break;
    }
    i++;
  }

  i = 0;
  if (!check_norm) {
    for (const AOGaussianPrimitive& gaussian : shell) {
      std::cout << "Ref:" << ref_results[i]
                << " result:" << gaussian.getContraction() << std::endl;
      i++;
    }
  }
  BOOST_CHECK_EQUAL(check_norm, 1);
}

BOOST_AUTO_TEST_CASE(Serializing) {

  Orbitals orbitals;
  orbitals.QMAtoms().LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                                  "/aobasis/molecule.xyz");
  BasisSet basis;
  basis.Load(std::string(XTP_TEST_DATA_FOLDER) + "/aobasis/3-21G.xml");
  AOBasis aobasis;
  aobasis.Fill(basis, orbitals.QMAtoms());

  CheckpointFile ff("aobasis.hdf5");
  CheckpointWriter ww = ff.getWriter();
  aobasis.WriteToCpt(ww);

  CheckpointReader rr = ff.getReader();
  AOBasis aobasis2;
  aobasis2.ReadFromCpt(rr);

  // no real way to test if two aobasis are equal so we check if the matrices
  // are the same
  AOOverlap overlap1;
  overlap1.Fill(aobasis);
  AOOverlap overlap2;
  overlap2.Fill(aobasis2);
  bool check = overlap1.Matrix().isApprox(overlap2.Matrix(), 1e-8);
  BOOST_CHECK_EQUAL(check, 1);

  auto func_per_atom1 = aobasis.getFuncPerAtom();
  auto func_per_atom2 = aobasis2.getFuncPerAtom();

  for (size_t i = 0; i < func_per_atom1.size(); i++) {
    BOOST_CHECK_EQUAL(func_per_atom1[i], func_per_atom2[i]);
  }
}

BOOST_AUTO_TEST_CASE(Adding_a_basis) {

  QMMolecule mol("a", 0);
  mol.LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) + "/aobasis/molecule.xyz");

  QMMolecule mol2 = mol;
  Eigen::Vector3d shift = {1, 2, 3};
  mol2.Translate(shift);
  BasisSet basis;
  basis.Load(std::string(XTP_TEST_DATA_FOLDER) + "/aobasis/3-21G.xml");
  AOBasis aobasis;
  aobasis.Fill(basis, mol);

  AOBasis aobasis2;
  aobasis2.Fill(basis, mol2);

  AOBasis combined = aobasis;
  combined.add(aobasis2);

  // no real way to test if two aobasis are equal so we check if the matrices
  // are the same
  AOOverlap overlap1;
  overlap1.Fill(aobasis);
  AOOverlap overlap2;
  overlap2.Fill(aobasis2);
  AOOverlap overlapcombined;
  overlapcombined.Fill(combined);
  bool check = overlap1.Matrix().isApprox(
      overlapcombined.Matrix().topLeftCorner(aobasis.AOBasisSize(),
                                             aobasis.AOBasisSize()),
      1e-8);
  BOOST_CHECK_EQUAL(check, 1);

  bool check2 = overlap2.Matrix().isApprox(
      overlapcombined.Matrix().bottomRightCorner(aobasis2.AOBasisSize(),
                                                 aobasis2.AOBasisSize()),
      1e-8);
  BOOST_CHECK_EQUAL(check2, 1);

  auto func_per_atom1 = aobasis.getFuncPerAtom();
  auto func_per_atom2 = aobasis2.getFuncPerAtom();
  func_per_atom1.insert(func_per_atom1.end(), func_per_atom2.begin(),
                        func_per_atom2.end());
  auto func_per_atom3 = combined.getFuncPerAtom();

  for (size_t i = 0; i < func_per_atom3.size(); i++) {
    BOOST_CHECK_EQUAL(func_per_atom1[i], func_per_atom3[i]);
  }
}

BOOST_AUTO_TEST_SUITE_END()
