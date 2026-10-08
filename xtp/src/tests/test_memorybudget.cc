/*
 * Copyright 2009-2026 The VOTCA Development Team (http://www.votca.org)
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

#define BOOST_TEST_MODULE memorybudget_test

// Third party includes
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/xtp/memorybudget.h"

using namespace votca::xtp;

BOOST_AUTO_TEST_SUITE(memorybudget_test)

BOOST_AUTO_TEST_CASE(budget_shares_and_defaults) {
  // no budget: the fixed allowance, whatever the process holds
  MemoryBudget::Set(0.0);
  BOOST_CHECK(!MemoryBudget::Given());
  BOOST_CHECK_EQUAL(MemoryBudget::AvailableBytes(),
                    MemoryBudget::DefaultAllowanceBytes());

  // a budget: its share minus the resident memory and a margin
  MemoryBudget::Set(100.0);
  BOOST_CHECK(MemoryBudget::Given());
  BOOST_CHECK_CLOSE(MemoryBudget::JobBudgetBytes(), 100e9, 1e-12);
  const double resident = MemoryBudget::ResidentBytes();
  if (resident >= 0) {
    BOOST_CHECK_CLOSE(MemoryBudget::AvailableBytes(), 100e9 - resident - 5e9,
                      1e-3);
  }
  // shared by concurrent jobs
  MemoryBudget::SetConcurrentJobs(4);
  BOOST_CHECK_CLOSE(MemoryBudget::JobBudgetBytes(), 25e9, 1e-12);
  BOOST_CHECK(MemoryBudget::AvailableBytes() < 25e9);
  // a budget smaller than what is in use leaves nothing, not a negative
  MemoryBudget::Set(1e-3);
  BOOST_CHECK_EQUAL(MemoryBudget::AvailableBytes(), 0.0);
  MemoryBudget::Set(0.0);
}

BOOST_AUTO_TEST_SUITE_END()
