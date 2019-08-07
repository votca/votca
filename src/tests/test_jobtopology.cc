/*
 * Copyright 2009-2019 The VOTCA Development Team (http://www.votca.org)
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

#define BOOST_TEST_MODULE jobtopology_test
#include <boost/test/unit_test.hpp>
#include <votca/xtp/jobtopology.h>

using namespace votca::xtp;
using namespace std;

BOOST_AUTO_TEST_SUITE(jobtopology_test)

BOOST_AUTO_TEST_CASE(constructor) {

  ofstream jobstream("job.xml");
  jobstream << "	<job>" << std::endl;
  jobstream << "		<id>0</id>" << std::endl;
  jobstream << "		<tag>seg0:n</tag>" << std::endl;
  jobstream << "		<input>" << std::endl;
  jobstream << "			<site_energy>0:n</site_energy>"
            << std::endl;
  jobstream << "			<regions>" << std::endl;
  jobstream << "				<region>" << std::endl;
  jobstream << "					<id>0</id>" << std::endl;
  jobstream
      << "					<segments>0:n</segments>"
      << std::endl;
  jobstream << "				</region>" << std::endl;
  jobstream << "			</regions>" << std::endl;
  jobstream << "		</input>" << std::endl;
  jobstream << "		<status>AVAILABLE</status>" << std::endl;
  jobstream << "	</job>" << std::endl;

  votca::tools::Property prop;
  votca::tools::load_property_from_xml(prop, "job.xml");
  std::string workdir = ".";
  Logger log;
  Job job(prop.get("job"));
  JobTopology top(job, log, workdir);
}

BOOST_AUTO_TEST_SUITE_END()
