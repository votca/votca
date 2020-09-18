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

#include "votca/xtp/qmmolecule.h"
#define BOOST_TEST_MAIN

#define BOOST_TEST_MODULE ecpaobasis_test

// Standard includes
#include <fstream>

// Third party includes
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/xtp/aopotential.h"
#include "votca/xtp/ecpaobasis.h"
#include "votca/xtp/orbitals.h"

using namespace votca::xtp;
using namespace votca;

BOOST_AUTO_TEST_SUITE(ecpaobasis_test)

BOOST_AUTO_TEST_CASE(Serializing) {

  Orbitals orbitals;

  orbitals.QMAtoms().LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                                  "/ecpaobasis/molecule.xyz");

  QMMolecule mol = orbitals.QMAtoms();

  ECPBasisSet ecp;
  ecp.Load(std::string(XTP_TEST_DATA_FOLDER) + "/ecpaobasis/ecp.xml");
  ECPAOBasis aoecp;
  aoecp.Fill(ecp, orbitals.QMAtoms());

  CheckpointFile ff("aoecp.hdf5");
  CheckpointWriter ww = ff.getWriter();
  aoecp.WriteToCpt(ww);

  CheckpointReader rr = ff.getReader();
  ECPAOBasis aoecp2;
  aoecp2.ReadFromCpt(rr);

  BasisSet s;
  s.Load(std::string(XTP_TEST_DATA_FOLDER) + "/ecpaobasis/3-21G.xml");
  AOBasis b;
  b.Fill(s, orbitals.QMAtoms());

  // no real way to test if two ecps are equal so we check if the matrices
  // are the same

  AOECP ecpmat1;
  ecpmat1.FillPotential(b, aoecp);
  AOECP ecpmat2;
  ecpmat2.FillPotential(b, aoecp2);
  bool check = ecpmat1.Matrix().isApprox(ecpmat2.Matrix(), 1e-8);
  BOOST_CHECK_EQUAL(check, 1);

  aoecp2.AddECPChargeToMolecule(mol);

  for (Index i = 0; i < mol.size(); i++) {
    BOOST_CHECK_EQUAL(mol[i].getNuccharge(),
                      orbitals.QMAtoms()[i].getNuccharge());
  }
}

BOOST_AUTO_TEST_SUITE_END()
