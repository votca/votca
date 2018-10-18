/*
 * Copyright 2009-2018 The VOTCA Development Team (http://www.votca.org)
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

#define BOOST_TEST_MODULE basisset_test
#include <boost/test/unit_test.hpp>
#include <votca/xtp/basisset.h>
#include <fstream>
#include <iostream>
#include <votca/xtp/orbitals.h>
#include <votca/xtp/aobasis.h>
#include <votca/xtp/aoshell.h>
using namespace votca::xtp;
using namespace std;
BOOST_AUTO_TEST_SUITE(basisset_test)

BOOST_AUTO_TEST_CASE(Contraction_test) {   
std::ofstream basisfile("contracted.xml");
basisfile<<"<basis name=\"cc-pVTZ\">"<<std::endl;
basisfile<<"<element name=\"C\">"<<std::endl;
basisfile<<"		<shell scale=\"1.0\" type=\"S\">"<<std::endl;
basisfile<<"			<constant decay=\"8236.0\">"<<std::endl;
basisfile<<"				<contractions factor=\"0.0005424302\" type=\"S\"/>"<<std::endl;
basisfile<<"			</constant>"<<std::endl;
basisfile<<"			<constant decay=\"1235.0\">"<<std::endl;
basisfile<<"				<contractions factor=\"0.0041964279\" type=\"S\"/>"<<std::endl;
basisfile<<"			</constant>"<<std::endl;
basisfile<<"			<constant decay=\"280.8\">"<<std::endl;
basisfile<<"				<contractions factor=\"0.0215409141\" type=\"S\"/>"<<std::endl;
basisfile<<"			</constant>"<<std::endl;
basisfile<<"			<constant decay=\"79.27\">"<<std::endl;
basisfile<<"				<contractions factor=\"0.0836149496\" type=\"S\"/>"<<std::endl;
basisfile<<"			</constant>"<<std::endl;
basisfile<<"			<constant decay=\"25.59\">"<<std::endl;
basisfile<<"				<contractions factor=\"0.2398716189\" type=\"S\"/>"<<std::endl;
basisfile<<"			</constant>"<<std::endl;
basisfile<<"			<constant decay=\"8.997\">"<<std::endl;
basisfile<<"				<contractions factor=\"0.4437518201\" type=\"S\"/>"<<std::endl;
basisfile<<"			</constant>"<<std::endl;
basisfile<<"			<constant decay=\"3.319\">"<<std::endl;
basisfile<<"				<contractions factor=\"0.3535796965\" type=\"S\"/>"<<std::endl;
basisfile<<"			</constant>"<<std::endl;
basisfile<<"			<constant decay=\"0.3643\">"<<std::endl;
basisfile<<"				<contractions factor=\"-0.0091763661\" type=\"S\"/>"<<std::endl;
basisfile<<"			</constant>"<<std::endl;
basisfile<<"		</shell>"<<std::endl;
basisfile<<"		<shell scale=\"1.0\" type=\"S\">"<<std::endl;
basisfile<<"			<constant decay=\"8236.0\">"<<std::endl;
basisfile<<"				<contractions factor=\"-0.0001963922\" type=\"S\"/>"<<std::endl;
basisfile<<"			</constant>"<<std::endl;
basisfile<<"			<constant decay=\"1235.0\">"<<std::endl;
basisfile<<"				<contractions factor=\"-0.0015259503\" type=\"S\"/>"<<std::endl;
basisfile<<"			</constant>"<<std::endl;
basisfile<<"			<constant decay=\"280.8\">"<<std::endl;
basisfile<<"				<contractions factor=\"-0.007890449\" type=\"S\"/>"<<std::endl;
basisfile<<"			</constant>"<<std::endl;
basisfile<<"			<constant decay=\"79.27\">"<<std::endl;
basisfile<<"				<contractions factor=\"-0.0315148705\" type=\"S\"/>"<<std::endl;
basisfile<<"			</constant>"<<std::endl;
basisfile<<"			<constant decay=\"25.59\">"<<std::endl;
basisfile<<"				<contractions factor=\"-0.0969100083\" type=\"S\"/>"<<std::endl;
basisfile<<"			</constant>"<<std::endl;
basisfile<<"			<constant decay=\"8.997\">"<<std::endl;
basisfile<<"				<contractions factor=\"-0.2205415263\" type=\"S\"/>"<<std::endl;
basisfile<<"			</constant>"<<std::endl;
basisfile<<"			<constant decay=\"3.319\">"<<std::endl;
basisfile<<"				<contractions factor=\"-0.2960691129\" type=\"S\"/>"<<std::endl;
basisfile<<"			</constant>"<<std::endl;
basisfile<<"			<constant decay=\"0.3643\">"<<std::endl;
basisfile<<"				<contractions factor=\"1.0405034329\" type=\"S\"/>"<<std::endl;
basisfile<<"			</constant>"<<std::endl;
basisfile<<"		</shell>"<<std::endl;
basisfile<<"		<shell scale=\"1.0\" type=\"S\">"<<std::endl;
basisfile<<"			<constant decay=\"0.9059\">"<<std::endl;
basisfile<<"				<contractions factor=\"1.0\" type=\"S\"/>"<<std::endl;
basisfile<<"			</constant>"<<std::endl;
basisfile<<"		</shell>"<<std::endl;
basisfile<<"		<shell scale=\"1.0\" type=\"S\">"<<std::endl;
basisfile<<"			<constant decay=\"0.1285\">"<<std::endl;
basisfile<<"				<contractions factor=\"1.0\" type=\"S\"/>"<<std::endl;
basisfile<<"			</constant>"<<std::endl;
basisfile<<"		</shell>"<<std::endl;
basisfile<<"		<shell scale=\"1.0\" type=\"P\">"<<std::endl;
basisfile<<"			<constant decay=\"18.71\">"<<std::endl;
basisfile<<"				<contractions factor=\"0.0394263872\" type=\"P\"/>"<<std::endl;
basisfile<<"			</constant>"<<std::endl;
basisfile<<"			<constant decay=\"4.133\">"<<std::endl;
basisfile<<"				<contractions factor=\"0.2440889849\" type=\"P\"/>"<<std::endl;
basisfile<<"			</constant>"<<std::endl;
basisfile<<"			<constant decay=\"1.2\">"<<std::endl;
basisfile<<"				<contractions factor=\"0.8154920089\" type=\"P\"/>"<<std::endl;
basisfile<<"			</constant>"<<std::endl;
basisfile<<"		</shell>"<<std::endl;
basisfile<<"		<shell scale=\"1.0\" type=\"P\">"<<std::endl;
basisfile<<"			<constant decay=\"0.3827\">"<<std::endl;
basisfile<<"				<contractions factor=\"1.0\" type=\"P\"/>"<<std::endl;
basisfile<<"			</constant>"<<std::endl;
basisfile<<"		</shell>"<<std::endl;
basisfile<<"		<shell scale=\"1.0\" type=\"P\">"<<std::endl;
basisfile<<"			<constant decay=\"0.1209\">"<<std::endl;
basisfile<<"				<contractions factor=\"1.0\" type=\"P\"/>"<<std::endl;
basisfile<<"			</constant>"<<std::endl;
basisfile<<"		</shell>"<<std::endl;
basisfile<<"		<shell scale=\"1.0\" type=\"D\">"<<std::endl;
basisfile<<"			<constant decay=\"1.097\">"<<std::endl;
basisfile<<"				<contractions factor=\"1.0\" type=\"D\"/>"<<std::endl;
basisfile<<"			</constant>"<<std::endl;
basisfile<<"		</shell>"<<std::endl;
basisfile<<"		<shell scale=\"1.0\" type=\"D\">"<<std::endl;
basisfile<<"			<constant decay=\"0.318\">"<<std::endl;
basisfile<<"				<contractions factor=\"1.0\" type=\"D\"/>"<<std::endl;
basisfile<<"			</constant>"<<std::endl;
basisfile<<"		</shell>"<<std::endl;
basisfile<<"		<shell scale=\"1.0\" type=\"F\">"<<std::endl;
basisfile<<"			<constant decay=\"0.761\">"<<std::endl;
basisfile<<"				<contractions factor=\"1.0\" type=\"F\"/>"<<std::endl;
basisfile<<"			</constant>"<<std::endl;
basisfile<<"		</shell>"<<std::endl;
basisfile<<"	</element>"<<std::endl;
basisfile<<"<element name=\"O\">"<<std::endl;
basisfile<<"		<shell scale=\"1.0\" type=\"S\">"<<std::endl;
basisfile<<"			<constant decay=\"15330.0\">"<<std::endl;
basisfile<<"				<contractions factor=\"0.0005201983\" type=\"S\"/>"<<std::endl;
basisfile<<"			</constant>"<<std::endl;
basisfile<<"			<constant decay=\"2299.0\">"<<std::endl;
basisfile<<"				<contractions factor=\"0.0040233448\" type=\"S\"/>"<<std::endl;
basisfile<<"			</constant>"<<std::endl;
basisfile<<"			<constant decay=\"522.4\">"<<std::endl;
basisfile<<"				<contractions factor=\"0.0207290833\" type=\"S\"/>"<<std::endl;
basisfile<<"			</constant>"<<std::endl;
basisfile<<"			<constant decay=\"147.3\">"<<std::endl;
basisfile<<"				<contractions factor=\"0.0810823271\" type=\"S\"/>"<<std::endl;
basisfile<<"			</constant>"<<std::endl;
basisfile<<"			<constant decay=\"47.55\">"<<std::endl;
basisfile<<"				<contractions factor=\"0.2362263521\" type=\"S\"/>"<<std::endl;
basisfile<<"			</constant>"<<std::endl;
basisfile<<"			<constant decay=\"16.76\">"<<std::endl;
basisfile<<"				<contractions factor=\"0.4435182094\" type=\"S\"/>"<<std::endl;
basisfile<<"			</constant>"<<std::endl;
basisfile<<"			<constant decay=\"6.207\">"<<std::endl;
basisfile<<"				<contractions factor=\"0.3586705887\" type=\"S\"/>"<<std::endl;
basisfile<<"			</constant>"<<std::endl;
basisfile<<"			<constant decay=\"0.6882\">"<<std::endl;
basisfile<<"				<contractions factor=\"-0.0083497972\" type=\"S\"/>"<<std::endl;
basisfile<<"			</constant>"<<std::endl;
basisfile<<"		</shell>"<<std::endl;
basisfile<<"		<shell scale=\"1.0\" type=\"S\">"<<std::endl;
basisfile<<"			<constant decay=\"15330.0\">"<<std::endl;
basisfile<<"				<contractions factor=\"-0.000197236\" type=\"S\"/>"<<std::endl;
basisfile<<"			</constant>"<<std::endl;
basisfile<<"			<constant decay=\"2299.0\">"<<std::endl;
basisfile<<"				<contractions factor=\"-0.0015350107\" type=\"S\"/>"<<std::endl;
basisfile<<"			</constant>"<<std::endl;
basisfile<<"			<constant decay=\"522.4\">"<<std::endl;
basisfile<<"				<contractions factor=\"-0.0079511839\" type=\"S\"/>"<<std::endl;
basisfile<<"			</constant>"<<std::endl;
basisfile<<"			<constant decay=\"147.3\">"<<std::endl;
basisfile<<"				<contractions factor=\"-0.0321134529\" type=\"S\"/>"<<std::endl;
basisfile<<"			</constant>"<<std::endl;
basisfile<<"			<constant decay=\"47.55\">"<<std::endl;
basisfile<<"				<contractions factor=\"-0.100269643\" type=\"S\"/>"<<std::endl;
basisfile<<"			</constant>"<<std::endl;
basisfile<<"			<constant decay=\"16.76\">"<<std::endl;
basisfile<<"				<contractions factor=\"-0.2340471118\" type=\"S\"/>"<<std::endl;
basisfile<<"			</constant>"<<std::endl;
basisfile<<"			<constant decay=\"6.207\">"<<std::endl;
basisfile<<"				<contractions factor=\"-0.3014109278\" type=\"S\"/>"<<std::endl;
basisfile<<"			</constant>"<<std::endl;
basisfile<<"			<constant decay=\"0.6882\">"<<std::endl;
basisfile<<"				<contractions factor=\"1.0349196495\" type=\"S\"/>"<<std::endl;
basisfile<<"			</constant>"<<std::endl;
basisfile<<"		</shell>"<<std::endl;
basisfile<<"		<shell scale=\"1.0\" type=\"S\">"<<std::endl;
basisfile<<"			<constant decay=\"1.752\">"<<std::endl;
basisfile<<"				<contractions factor=\"1.0\" type=\"S\"/>"<<std::endl;
basisfile<<"			</constant>"<<std::endl;
basisfile<<"		</shell>"<<std::endl;
basisfile<<"		<shell scale=\"1.0\" type=\"S\">"<<std::endl;
basisfile<<"			<constant decay=\"0.2384\">"<<std::endl;
basisfile<<"				<contractions factor=\"1.0\" type=\"S\"/>"<<std::endl;
basisfile<<"			</constant>"<<std::endl;
basisfile<<"		</shell>"<<std::endl;
basisfile<<"		<shell scale=\"1.0\" type=\"P\">"<<std::endl;
basisfile<<"			<constant decay=\"34.46\">"<<std::endl;
basisfile<<"				<contractions factor=\"0.0411634896\" type=\"P\"/>"<<std::endl;
basisfile<<"			</constant>"<<std::endl;
basisfile<<"			<constant decay=\"7.749\">"<<std::endl;
basisfile<<"				<contractions factor=\"0.2577628359\" type=\"P\"/>"<<std::endl;
basisfile<<"			</constant>"<<std::endl;
basisfile<<"			<constant decay=\"2.28\">"<<std::endl;
basisfile<<"				<contractions factor=\"0.8024192744\" type=\"P\"/>"<<std::endl;
basisfile<<"			</constant>"<<std::endl;
basisfile<<"		</shell>"<<std::endl;
basisfile<<"		<shell scale=\"1.0\" type=\"P\">"<<std::endl;
basisfile<<"			<constant decay=\"0.7156\">"<<std::endl;
basisfile<<"				<contractions factor=\"1.0\" type=\"P\"/>"<<std::endl;
basisfile<<"			</constant>"<<std::endl;
basisfile<<"		</shell>"<<std::endl;
basisfile<<"		<shell scale=\"1.0\" type=\"P\">"<<std::endl;
basisfile<<"			<constant decay=\"0.214\">"<<std::endl;
basisfile<<"				<contractions factor=\"1.0\" type=\"P\"/>"<<std::endl;
basisfile<<"			</constant>"<<std::endl;
basisfile<<"		</shell>"<<std::endl;
basisfile<<"		<shell scale=\"1.0\" type=\"D\">"<<std::endl;
basisfile<<"			<constant decay=\"2.314\">"<<std::endl;
basisfile<<"				<contractions factor=\"1.0\" type=\"D\"/>"<<std::endl;
basisfile<<"			</constant>"<<std::endl;
basisfile<<"		</shell>"<<std::endl;
basisfile<<"		<shell scale=\"1.0\" type=\"D\">"<<std::endl;
basisfile<<"			<constant decay=\"0.645\">"<<std::endl;
basisfile<<"				<contractions factor=\"1.0\" type=\"D\"/>"<<std::endl;
basisfile<<"			</constant>"<<std::endl;
basisfile<<"		</shell>"<<std::endl;
basisfile<<"		<shell scale=\"1.0\" type=\"F\">"<<std::endl;
basisfile<<"			<constant decay=\"1.428\">"<<std::endl;
basisfile<<"				<contractions factor=\"1.0\" type=\"F\"/>"<<std::endl;
basisfile<<"			</constant>"<<std::endl;
basisfile<<"		</shell>"<<std::endl;
basisfile<<"	</element>"<<std::endl;
basisfile<<"</basis>"<<std::endl;
basisfile.close();
   
std::ofstream xyzfile("CO.xyz");
xyzfile << " 2" << std::endl;
xyzfile << " CO" << std::endl;
xyzfile << " C            .000000     .000000     .000000" << std::endl;
xyzfile << " O            1.000000     .000000     .000000" << std::endl;
xyzfile.close();

Orbitals orbitals;
orbitals.QMAtoms().LoadFromXYZ("CO.xyz");
BasisSet basis;
basis.LoadBasisSet("contracted.xml");

}

BOOST_AUTO_TEST_SUITE_END()
