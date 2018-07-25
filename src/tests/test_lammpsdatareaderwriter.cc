/*
 * Copyright 2009-2018 The VOTCA Development Team (http://www.votca.org)
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

#define BOOST_TEST_MAIN

#define BOOST_TEST_MODULE lammpdatatopologyreaderwriter_test
#include <boost/test/unit_test.hpp>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <votca/tools/matrix.h>
#include <votca/tools/types.h>
#include <votca/tools/elements.h>
#include <votca/csg/orthorhombicbox.h>
#include <votca/csg/bead.h>
#include <votca/csg/trajectorywriter.h>
#include <votca/csg/trajectoryreader.h>

using namespace std;
using namespace votca::csg;
using namespace votca::tools;

// used for rounding doubles so we can compare them
double round_(double v, int p);

// Check if file exists
bool fexists_(const string filename);

void printTestFile_(const string filename);

BOOST_AUTO_TEST_SUITE(lammpsdatareaderwriter_test)

	/**
	 * \brief Test the trajectory reader
	 *
	 * This test is designed to test the trajectory reader this is done by
	 * creating a small lammps data file. A topology object is created with
	 * some default values. The file is then read in with the 
	 * trajectory reader and the values in the top object are then examined
	 * to ensure they no longer represent the default state but the values 
	 * from the file. 
	 */
	BOOST_AUTO_TEST_CASE(test_trajectoryreader){

		// Create a .data file with (polymer strand)
		// and read from it. Create a topology object with the same 
		// molecule to enable the ability to read in the trajectory 
		// file
		string lammpsdatafilename = "test_polymer.data";
		if(fexists_(lammpsdatafilename)){
			remove(lammpsdatafilename.c_str());
		}
		printTestFile_(lammpsdatafilename);

		Topology top;


	}

/** 
 * \brief Testing trajectory writer
 *
 * This test first creates a topology object and assigns default values to
 * it. It then writes the topology info to a lammps data file. The data 
 * file is then read into the topology file and the values are compared. 
 */
BOOST_AUTO_TEST_CASE(test_trajectorywriter) {

	// Create a topology object with a simple system (polymer strand)
	// and write it to a lammps data file
	Topology top;

}

BOOST_AUTO_TEST_SUITE_END()

/*****************************************************************************
 * Internal test functions                                                   *
 *****************************************************************************/

// used for rounding doubles so we can compare them
double round_(double v, int p) {
	v *= pow(10, p);
	v = round(v);
	v /= pow(10, p);
	return v;
}

// Check if file exists
bool fexists_(const string filename) {
	std::ifstream ifile(filename);
	return (bool)ifile;
}

void printTestFile_(const string filename){

	ofstream outfile(filename);
	outfile << "LAMMPS data file via write_data, version 17 Nov 2015, timestep = 961" << endl;
	outfile << "" << endl;
	outfile << "100 atoms" << endl;
	outfile << "1 atom types" << endl;
	outfile << "99 bonds" << endl;
	outfile << "1 bond types" << endl;
	outfile << "98 angles" << endl;
	outfile << "1 angle types" << endl;
	outfile << "97 dihedrals" << endl;
	outfile << "1 dihedral types" << endl;
	outfile << "" << endl;
	outfile << "0.0000000000000000e+00 1.5850000000000000e+02 xlo xhi" << endl;
	outfile << "0.0000000000000000e+00 1.5850000000000000e+02 ylo yhi" << endl;
	outfile << "0.0000000000000000e+00 1.0000000000000000e+02 zlo zhi" << endl;
	outfile << "" << endl;
	outfile << "Masses" << endl;
	outfile << "" << endl;
	outfile << "1 14.02" << endl;
	outfile << "" << endl;
	outfile << "Pair Coeffs # lj/cut" << endl;
	outfile << "" << endl;
	outfile << "1 0.112 4.01" << endl;
	outfile << "" << endl;
	outfile << "Bond Coeffs # harmonic" << endl;
	outfile << "" << endl;
	outfile << "1 350 1.53" << endl;
	outfile << "" << endl;
	outfile << "Angle Coeffs # harmonic" << endl;
	outfile << "" << endl;
	outfile << "1 60 109.5" << endl;
	outfile << "" << endl;
	outfile << "Atoms # molecular" << endl;
	outfile << "" << endl;
	outfile << "5 1 1 5.8137734438047275e+01 5.2820532540049605e+01 5.1610735195303207e+01 0 0 0" << endl;
	outfile << "6 1 1 5.7593847093164598e+01 5.3264325139761141e+01 5.0249644669646678e+01 0 0 0" << endl;
	outfile << "7 1 1 5.6511930501583258e+01 5.4329619289202022e+01 5.0444704115951488e+01 0 0 0" << endl;
	outfile << "9 1 1 5.8047091053480429e+01 5.6222029048615425e+01 4.9874845165148727e+01 0 0 0" << endl;
	outfile << "8 1 1 5.7157290110423276e+01 5.5616832485184531e+01 5.0964522757346906e+01 0 0 0" << endl;
	outfile << "10 1 1 5.8796606242207154e+01 5.7431388649541120e+01 5.0436740373291755e+01 0 0 0" << endl;
	outfile << "13 1 1 6.1470904287728864e+01 5.8661356795663160e+01 5.0956738288440924e+01 0 0 0" << endl;
	outfile << "14 1 1 6.2318874922711110e+01 5.9820763044333191e+01 5.1482733521136190e+01 0 0 0" << endl;
	outfile << "11 1 1 5.9740399827120577e+01 5.6964684624717705e+01 5.1547015047030342e+01 0 0 0" << endl;
	outfile << "12 1 1 6.0550755372027595e+01 5.8152824805563078e+01 5.2068213576080311e+01 0 0 0" << endl;
	outfile << "15 1 1 6.3220534569625151e+01 5.9309964725714572e+01 5.2607828752412175e+01 0 0 0" << endl;
	outfile << "16 1 1 6.4122151012703725e+01 6.0440684531421851e+01 5.3106224006392388e+01 0 0 0" << endl;
	outfile << "17 1 1 6.4989383763573784e+01 5.9916464190367087e+01 5.4251845650219749e+01 0 0 0" << endl;
	outfile << "18 1 1 6.4098633759524319e+01 5.9598860291991549e+01 5.5454314244516617e+01 0 0 0" << endl;
	outfile << "19 1 1 6.4944652589155282e+01 5.8985111104358268e+01 5.6571073918330995e+01 0 0 0" << endl;
	outfile << "20 1 1 6.5994561083409252e+01 5.9993009252428024e+01 5.7042894764699355e+01 0 0 0" << endl;
	outfile << "22 1 1 6.7926851198445149e+01 6.0321059085483142e+01 5.8608478010470662e+01 0 0 0" << endl;
	outfile << "23 1 1 6.8901563116380558e+01 6.0573575982606258e+01 5.7456282631605447e+01 0 0 0" << endl;
	outfile << "21 1 1 6.6831300248918623e+01 5.9354803574005928e+01 5.8152660598019651e+01 0 0 0" << endl;
	outfile << "1 1 1 6.2806039595490141e+01 5.2512717926489422e+01 4.9887303784980617e+01 0 0 0" << endl;
	outfile << "3 1 1 6.0466882959009020e+01 5.2652416089059159e+01 5.0746419497822217e+01 0 0 0" << endl;
	outfile << "4 1 1 5.9321955572837354e+01 5.1875910671417962e+01 5.1400092405450614e+01 0 0 0" << endl;
	outfile << "2 1 1 6.1665374508442042e+01 5.1725903855301063e+01 5.0535634777284422e+01 0 0 0" << endl;
	outfile << "24 1 1 7.0088507511783732e+01 6.1393197682645742e+01 5.7966976153016141e+01 0 0 0" << endl;
	outfile << "25 1 1 7.0848284090923443e+01 6.0554827815202422e+01 5.8996738916525899e+01 0 0 0" << endl;
	outfile << "26 1 1 7.2074450966338077e+01 6.1312416880982234e+01 5.9509767238047807e+01 0 0 0" << endl;
	outfile << "27 1 1 7.1640887270430980e+01 6.2567844751612334e+01 6.0268733565300963e+01 0 0 0" << endl;
	outfile << "28 1 1 7.2872102996450593e+01 6.3190209979261759e+01 6.0930409638305782e+01 0 0 0" << endl;
	outfile << "30 1 1 7.5169466126282174e+01 6.4075591986631721e+01 6.0538420613075047e+01 0 0 0" << endl;
	outfile << "31 1 1 7.4966311016695201e+01 6.5312955853856934e+01 6.1415663004489325e+01 0 0 0" << endl;
	outfile << "29 1 1 7.3857038025118229e+01 6.3671680434633977e+01 5.9861961539904001e+01 0 0 0" << endl;
	outfile << "32 1 1 7.4876121418370914e+01 6.6569610690469432e+01 6.0546437243887951e+01 0 0 0" << endl;
	outfile << "33 1 1 7.4715912082279033e+01 6.7788173718226062e+01 6.1457064952540193e+01 0 0 0" << endl;
	outfile << "34 1 1 7.4830560097661390e+01 6.9072909689161762e+01 6.0634855859033337e+01 0 0 0" << endl;
	outfile << "35 1 1 7.4709645502788135e+01 7.0274258093193097e+01 6.1576173418651145e+01 0 0 0" << endl;
	outfile << "36 1 1 7.5010835174024820e+01 7.1575619503013996e+01 6.0825198644724438e+01 0 0 0" << endl;
	outfile << "37 1 1 7.6493197887939786e+01 7.1625085032176173e+01 6.0437841872973479e+01 0 0 0" << endl;
	outfile << "77 1 1 8.9325614908665898e+01 8.0680262816921143e+01 6.0577003363101213e+01 0 0 0" << endl;
	outfile << "82 1 1 9.4253138814091812e+01 7.9091853937814449e+01 6.1540762438737737e+01 0 0 0" << endl;
	outfile << "83 1 1 9.5504280116567031e+01 7.8339602252550620e+01 6.1085185014966697e+01 0 0 0" << endl;
	outfile << "81 1 1 9.3547457379129199e+01 7.9682975603340012e+01 6.0319447707117909e+01 0 0 0" << endl;
	outfile << "84 1 1 9.6206032664762077e+01 7.7734909971149520e+01 6.2302004664526756e+01 0 0 0" << endl;
	outfile << "80 1 1 9.2279019213154214e+01 8.0411859731792688e+01 6.0765384931122867e+01 0 0 0" << endl;
	outfile << "85 1 1 9.7445478326252115e+01 7.6964192967750179e+01 6.1840617631852759e+01 0 0 0" << endl;
	outfile << "87 1 1 9.8786905710954187e+01 7.7464497107296893e+01 6.3900829432248614e+01 0 0 0" << endl;
	outfile << "86 1 1 9.8158347846061886e+01 7.6356873786335953e+01 6.3051475555508198e+01 0 0 0" << endl;
	outfile << "79 1 1 9.1572544651542529e+01 8.0991873739757821e+01 5.9538798993591456e+01 0 0 0" << endl;
	outfile << "38 1 1 7.7366931968969823e+01 7.1712215563035130e+01 6.1695176668606017e+01 0 0 0" << endl;
	outfile << "46 1 1 8.2867326178284316e+01 7.1164859884934685e+01 6.4670965608452619e+01 0 0 0" << endl;
	outfile << "43 1 1 8.1431831883962957e+01 7.3576180437134539e+01 6.3315300572605416e+01 0 0 0" << endl;
	outfile << "45 1 1 8.2938976862641297e+01 7.1563572608818035e+01 6.3194752323319385e+01 0 0 0" << endl;
	outfile << "44 1 1 8.1590912386546293e+01 7.2152844394267660e+01 6.2773327974293338e+01 0 0 0" << endl;
	outfile << "68 1 1 8.4531862597954628e+01 8.1281414571148204e+01 6.7827593132784131e+01 0 0 0" << endl;
	outfile << "98 1 1 1.0415802336592691e+02 7.6358596582309673e+01 6.1198416436081907e+01 0 0 0" << endl;
	outfile << "99 1 1 1.0370780457168856e+02 7.6827422710965806e+01 5.9814071418582913e+01 0 0 0" << endl;
	outfile << "100 1 1 1.0278495052273817e+02 7.8038808732039655e+01 5.9962906483300024e+01 0 0 0" << endl;
	outfile << "39 1 1 7.7252209503495536e+01 7.3093591786513144e+01 6.2352746390331660e+01 0 0 0" << endl;
	outfile << "40 1 1 7.8044827779464242e+01 7.4142034261801641e+01 6.1563624170662820e+01 0 0 0" << endl;
	outfile << "42 1 1 7.9970226879064057e+01 7.3999725933328548e+01 6.3158991319547482e+01 0 0 0" << endl;
	outfile << "41 1 1 7.9545330382415571e+01 7.3860433702007583e+01 6.1694333446766457e+01 0 0 0" << endl;
	outfile << "70 1 1 8.4409527738959412e+01 8.2578008698560396e+01 6.5691216835626065e+01 0 0 0" << endl;
	outfile << "69 1 1 8.3861556555824635e+01 8.1369588056649022e+01 6.6454424119676872e+01 0 0 0" << endl;
	outfile << "71 1 1 8.5885220100259431e+01 8.2353311829452366e+01 6.5351706991771749e+01 0 0 0" << endl;
	outfile << "73 1 1 8.7492743787320066e+01 8.0988615573905847e+01 6.4007439355582278e+01 0 0 0" << endl;
	outfile << "74 1 1 8.7620370400056515e+01 7.9876788057036194e+01 6.2961936151310034e+01 0 0 0" << endl;
	outfile << "72 1 1 8.6013679805442990e+01 8.1232506066144936e+01 6.4316860660638241e+01 0 0 0" << endl;
	outfile << "75 1 1 8.7073062888496963e+01 8.0357863304214618e+01 6.1614635633869128e+01 0 0 0" << endl;
	outfile << "76 1 1 8.8029170099206496e+01 8.1379831276688407e+01 6.0992553341540230e+01 0 0 0" << endl;
	outfile << "89 1 1 1.0057695670950783e+02 7.9209901405693046e+01 6.3981720510380036e+01 0 0 0" << endl;
	outfile << "88 1 1 9.9943044449021414e+01 7.8107969370881960e+01 6.3131608378479854e+01 0 0 0" << endl;
	outfile << "90 1 1 1.0119235454425196e+02 7.8589240972391963e+01 6.5239048021103315e+01 0 0 0" << endl;
	outfile << "91 1 1 1.0193393760015982e+02 7.9662098183193820e+01 6.6041355942009332e+01 0 0 0" << endl;
	outfile << "92 1 1 1.0316170728832537e+02 8.0121411076907719e+01 6.5252168429990675e+01 0 0 0" << endl;
	outfile << "97 1 1 1.0500806213866807e+02 7.7451009503315490e+01 6.1848913257186148e+01 0 0 0" << endl;
	outfile << "96 1 1 1.0540558800723898e+02 7.7014614332020813e+01 6.3261977383557557e+01 0 0 0" << endl;
	outfile << "94 1 1 1.0527976752008665e+02 7.9340855792933525e+01 6.4177356744590838e+01 0 0 0" << endl;
	outfile << "93 1 1 1.0415771620812434e+02 7.8964132433286139e+01 6.5148020705178652e+01 0 0 0" << endl;
	outfile << "95 1 1 1.0619145951604089e+02 7.8131710417899569e+01 6.3954571381514654e+01 0 0 0" << endl;
	outfile << "78 1 1 9.0290559831794539e+01 8.1702985417906305e+01 5.9975016466002131e+01 0 0 0" << endl;
	outfile << "52 1 1 8.1440133017846961e+01 6.9168533917762332e+01 7.0554093637393166e+01 0 0 0" << endl;
	outfile << "51 1 1 8.1951420104865321e+01 6.9708300491886376e+01 6.9216901470055106e+01 0 0 0" << endl;
	outfile << "53 1 1 8.1139395734809028e+01 7.0335404559593755e+01 7.1495882334663577e+01 0 0 0" << endl;
	outfile << "59 1 1 8.2915237922444760e+01 7.3812293512714973e+01 7.3941050702264917e+01 0 0 0" << endl;
	outfile << "60 1 1 8.3087946679492347e+01 7.5066710190083526e+01 7.3083449112750031e+01 0 0 0" << endl;
	outfile << "47 1 1 8.1816168669082302e+01 7.0062603425354112e+01 6.4825384398723614e+01 0 0 0" << endl;
	outfile << "49 1 1 8.1229273779412352e+01 7.0880874629970194e+01 6.7114733225531353e+01 0 0 0" << endl;
	outfile << "48 1 1 8.1663518653716537e+01 6.9664910326201991e+01 6.6294323906729460e+01 0 0 0" << endl;
	outfile << "50 1 1 8.0796780356164817e+01 7.0418582496791984e+01 6.8507340711124343e+01 0 0 0" << endl;
	outfile << "62 1 1 8.3185343092380009e+01 7.7560123403753167e+01 7.3138264265003045e+01 0 0 0" << endl;
	outfile << "63 1 1 8.2343260445774973e+01 7.7648127351366156e+01 7.1864318159643261e+01 0 0 0" << endl;
	outfile << "64 1 1 8.2786695955482301e+01 7.8875771790109852e+01 7.1067026862460082e+01 0 0 0" << endl;
	outfile << "65 1 1 8.4253521712145044e+01 7.8702543438875878e+01 7.0667440332348178e+01 0 0 0" << endl;
	outfile << "66 1 1 8.4739038846207052e+01 7.9943435119395744e+01 6.9916690024880992e+01 0 0 0" << endl;
	outfile << "67 1 1 8.4013699524302055e+01 8.0051869847146762e+01 6.8575119432467091e+01 0 0 0" << endl;
	outfile << "61 1 1 8.2808477726655866e+01 7.6308383286601369e+01 7.3932351774246925e+01 0 0 0" << endl;
	outfile << "54 1 1 8.0364287775279905e+01 6.9824924188188220e+01 7.2712240255197784e+01 0 0 0" << endl;
	outfile << "55 1 1 8.0172100064076687e+01 7.0972440593323071e+01 7.3704282647795253e+01 0 0 0" << endl;
	outfile << "56 1 1 8.1486234702804140e+01 7.1219005843720581e+01 7.4448774864110874e+01 0 0 0" << endl;
	outfile << "57 1 1 8.1356040603794312e+01 7.2464100885425907e+01 7.5328025439832388e+01 0 0 0" << endl;
	outfile << "58 1 1 8.1477273459174981e+01 7.3715469993842632e+01 7.4455793910286005e+01 0 0 0" << endl;
	outfile << "" << endl;
	outfile << "Velocities" << endl;
	outfile << "" << endl;
	outfile << "5 2.4762633191726287e-02 3.7712101178784003e-04 1.5769358250450456e-02" << endl;
	outfile << "6 4.9428431380202226e-03 1.0804180320459038e-02 2.4204258774847094e-03" << endl;
	outfile << "7 1.9452054479231866e-02 -1.9406319411225625e-02 -2.1385143854993158e-03" << endl;
	outfile << "9 -7.9893801910391332e-03 -6.0622972972037645e-03 -2.5306997868292510e-03" << endl;
	outfile << "8 -3.7990815273581538e-03 -3.1915707136970888e-03 -6.8016787572024130e-03" << endl;
	outfile << "10 -4.3317455563978251e-03 -3.5291411347403760e-02 -2.4748603833764364e-02" << endl;
	outfile << "13 1.1852189664289304e-03 -1.9173483195200052e-02 -1.9480571746689759e-03" << endl;
	outfile << "14 -9.3903507227086937e-03 -3.1679159764559831e-02 -7.0400810544051827e-03" << endl;
	outfile << "11 -1.4152172709569304e-02 -1.7020259039630063e-03 1.1237102729935737e-02" << endl;
	outfile << "12 -1.8301425172796027e-02 -8.8266569617304943e-03 -2.2500338235775699e-02" << endl;
	outfile << "15 2.0039395269746380e-02 3.1027580396968457e-02 -1.3614213866368465e-03" << endl;
	outfile << "16 -7.3085072941391155e-03 2.6535433528437791e-02 -6.1521562509386559e-03" << endl;
	outfile << "17 1.6709374465406234e-03 -6.2697199111260028e-03 -8.3895047544977458e-04" << endl;
	outfile << "18 -9.2999862554835507e-03 -1.9871852313778610e-02 9.3787240244730700e-04" << endl;
	outfile << "19 1.3384583039511622e-02 3.3431456140446576e-02 5.5249038095012892e-03" << endl;
	outfile << "20 9.2699698605564068e-03 3.4691500714413528e-02 -7.6639466981432624e-03" << endl;
	outfile << "22 -1.1241072808184854e-02 1.6664784344516594e-02 -6.9689951103898296e-03" << endl;
	outfile << "23 2.2367913250104277e-02 -2.3213698670117343e-03 1.0259233005471137e-02" << endl;
	outfile << "21 1.4549935701929731e-02 1.2995173042938948e-02 1.2507109711396788e-02" << endl;
	outfile << "1 -2.6641145036433292e-02 2.6190072117130156e-03 -1.9437903334184874e-02" << endl;
	outfile << "3 -9.7874601716667481e-03 2.4608048576392576e-02 -6.3318129601802765e-03" << endl;
	outfile << "4 -2.1360438661634912e-02 -2.6661667544608404e-02 -7.7821647704330319e-03" << endl;
	outfile << "2 -2.7130308256301864e-02 -1.3878973024088717e-03 -2.8057024579053831e-03" << endl;
	outfile << "24 -8.6078753300249874e-03 1.1593854471737771e-02 8.3089011212495890e-03" << endl;
	outfile << "25 1.9649313476471506e-02 -7.0021088596946123e-03 -8.0204710702006833e-03" << endl;
	outfile << "26 -1.4234890484802198e-02 -1.8286751528428895e-02 1.3647265223365740e-02" << endl;
	outfile << "27 -2.0280590607577743e-02 -1.6379697812307598e-02 3.4137198505099801e-03" << endl;
	outfile << "28 2.0074216429186720e-02 -4.9256614939011522e-03 -5.4266005495092725e-03" << endl;
	outfile << "30 -7.8725194102736082e-03 -9.4031638775822889e-03 -1.5798531573843918e-02" << endl;
	outfile << "31 8.9089839623491610e-03 2.2022802264774731e-03 3.0642400579465767e-02" << endl;
	outfile << "29 3.9212356654378826e-02 5.2091264431473189e-03 -4.0623234656550940e-03" << endl;
	outfile << "32 8.3925235498200981e-03 1.9917559147257211e-02 -1.7900208622364161e-02" << endl;
	outfile << "33 -8.0264441590715066e-03 3.3064233234334134e-02 1.1441442546196201e-02" << endl;
	outfile << "34 2.1128000569418524e-02 -4.6312987016410845e-03 -9.4292550851877224e-03" << endl;
	outfile << "35 -1.0761442895669767e-02 1.5640066866490363e-02 1.0922409385097188e-02" << endl;
	outfile << "36 -5.4260102161066857e-03 -2.8891165013975702e-02 3.0758902974153700e-02" << endl;
	outfile << "37 1.4879727341226747e-03 2.0199055148914548e-02 1.4738987481173291e-02" << endl;
	outfile << "77 -1.6573276895111278e-02 -1.1421948599648337e-03 -6.7826546496229621e-03" << endl;
	outfile << "82 2.8010730689331853e-02 -4.1000406681648414e-02 5.8622195616645752e-03" << endl;
	outfile << "83 -1.5479045896950589e-02 7.8872304747784269e-04 6.6513453037488901e-03" << endl;
	outfile << "81 -1.6049086642226859e-02 1.0221670734401533e-03 6.3003654927853113e-03" << endl;
	outfile << "84 -2.1928531948100308e-02 -5.5067164697141825e-03 5.4234407281079505e-03" << endl;
	outfile << "80 -1.6897643170552402e-02 2.2609651403303747e-03 -9.4200179691879028e-03" << endl;
	outfile << "85 -2.3612843763106413e-02 -1.7825877752803955e-02 4.5127649813374478e-04" << endl;
	outfile << "87 -2.9095159853415792e-02 -1.6343295508658748e-02 2.1341302865212630e-02" << endl;
	outfile << "86 3.3544927202013384e-03 9.2564707554206333e-03 2.6534784988163890e-02" << endl;
	outfile << "79 2.2516623677300861e-02 7.7632069780299044e-03 2.6844678746668932e-02" << endl;
	outfile << "38 -1.2203074573956540e-02 7.9997465117948221e-03 -9.8448119878856086e-05" << endl;
	outfile << "46 5.2863504055621173e-03 -1.2772642351490133e-02 1.9904269205190091e-03" << endl;
	outfile << "43 1.1130986475597845e-02 -1.0534977479942793e-02 -2.1215695600093031e-02" << endl;
	outfile << "45 2.0556493709178971e-02 2.1966213688943640e-02 -1.8384905728320491e-02" << endl;
	outfile << "44 -3.3146790019430942e-03 1.3453890798779244e-02 -5.2493172584168260e-03" << endl;
	outfile << "68 1.5211952137003628e-02 -5.6400759859683378e-03 1.0147559772200568e-02" << endl;
	outfile << "98 -1.3958291473051412e-02 2.7429880017782425e-03 1.4534048370076172e-02" << endl;
	outfile << "99 1.5371611827224103e-03 -1.7683122440084961e-02 2.2417001899310805e-02" << endl;
	outfile << "100 1.2956789664940821e-02 7.2651264839478837e-03 3.1712584943495317e-04" << endl;
	outfile << "39 -3.0858990017561268e-02 -2.1869068911718442e-02 2.4334894467195708e-02" << endl;
	outfile << "40 1.4725489000154482e-03 -2.7847485303198131e-02 3.6934368794878908e-02" << endl;
	outfile << "42 -7.8211015687012931e-03 -1.9906163415893100e-03 1.7771316111664690e-02" << endl;
	outfile << "41 -1.5563090096177526e-02 2.5391725443022921e-02 -2.2838116127881648e-02" << endl;
	outfile << "70 2.2574907751353439e-03 5.2720137114613125e-03 -5.6502101066269175e-03" << endl;
	outfile << "69 -2.0782298544526144e-02 -1.9123089612201313e-02 -5.5495421451414086e-03" << endl;
	outfile << "71 1.5199190772362998e-02 -2.0194965464532446e-02 -8.2940543155542550e-03" << endl;
	outfile << "73 -1.6733194200098458e-03 -2.5426522352323434e-02 7.1947801547689550e-03" << endl;
	outfile << "74 2.9385979581141593e-02 1.9753246802476903e-03 -1.5547553205971563e-02" << endl;
	outfile << "72 1.7419455036659887e-02 -7.4221784851914484e-03 -6.3724935552474087e-03" << endl;
	outfile << "75 1.3218074846576994e-02 -8.4692500140810013e-03 -1.3957714645733560e-02" << endl;
	outfile << "76 -9.7094440363709560e-03 -4.7481291676793585e-03 -1.3457755620754540e-02" << endl;
	outfile << "89 1.2803236068771474e-02 2.1431319014253332e-02 -5.7602042720647020e-03" << endl;
	outfile << "88 -5.5086145008018100e-03 -3.4188552760361666e-04 1.2851883465185576e-02" << endl;
	outfile << "90 1.4651065729978933e-02 -2.1659397768623492e-02 -1.3429783326948224e-02" << endl;
	outfile << "91 1.5538972728427542e-02 1.9221554970218064e-02 -6.2764525864740161e-03" << endl;
	outfile << "92 7.5070794796031825e-03 8.0673236809167010e-04 7.7067548495011072e-03" << endl;
	outfile << "97 -2.4115252194968975e-03 -2.3279059035975283e-03 1.4346519252967949e-02" << endl;
	outfile << "96 3.2646203255364041e-02 -4.4536795353391770e-03 2.2538156091833172e-03" << endl;
	outfile << "94 7.0694385399962222e-03 -1.8288484994808916e-02 -2.6458708662829903e-02" << endl;
	outfile << "93 -1.3940127631482876e-02 -1.5723190649898598e-02 4.9622714912770406e-03" << endl;
	outfile << "95 2.4833272288747690e-03 1.8041308633154535e-02 7.3921118765592972e-03" << endl;
	outfile << "78 -2.4054043783725516e-02 5.9694545475322582e-03 5.1009346878962761e-03" << endl;
	outfile << "52 -2.3333680454902868e-02 7.2426822107082801e-03 -1.1748386992241778e-04" << endl;
	outfile << "51 -2.2543560369082227e-02 -2.5852976015415714e-02 7.3427833507981521e-03" << endl;
	outfile << "53 -1.1217888805950368e-02 8.1302434447706705e-03 -4.6563963543899776e-03" << endl;
	outfile << "59 5.4303110363990350e-03 2.1461242335297798e-04 -3.2611206453245761e-03" << endl;
	outfile << "60 6.2607296660893022e-03 -1.8646978578999347e-02 -4.0975371835725062e-03" << endl;
	outfile << "47 2.3742877322270890e-02 1.5060613317396630e-02 -1.7888690536014551e-02" << endl;
	outfile << "49 -2.1119124274878238e-02 4.8418744490236650e-03 -6.9244336697025554e-03" << endl;
	outfile << "48 9.7188868453367081e-03 -9.6147820781278884e-03 1.9755353283227351e-02" << endl;
	outfile << "50 -2.3032122827548228e-03 -1.0584973567871944e-03 -4.0808953207440799e-03" << endl;
	outfile << "62 3.1433822359340002e-02 9.8461551003507560e-03 -1.3452413274615443e-02" << endl;
	outfile << "63 1.1733858605872626e-02 -3.1906921305811513e-02 -4.7401808899959477e-03" << endl;
	outfile << "64 -1.2110188149774643e-02 -4.4640083174155342e-03 -1.6107405102854154e-02" << endl;
	outfile << "65 -1.4450073083055176e-02 3.0230706268667210e-02 -4.6823830531080284e-04" << endl;
	outfile << "66 4.2081467182798057e-03 2.9028903065311237e-02 -2.0106880786399578e-02" << endl;
	outfile << "67 -2.5857871253387107e-04 -2.8253509334472220e-02 -3.5700373330979050e-03" << endl;
	outfile << "61 1.8360372621294543e-02 1.8030597281484120e-02 -1.9022806091684543e-03" << endl;
	outfile << "54 -5.7431940873534587e-03 1.0872788768227556e-02 1.3025488593651966e-02" << endl;
	outfile << "55 -2.9248565279124894e-02 2.0601220669770508e-03 -4.7066384051876022e-03" << endl;
	outfile << "56 -2.6726844290510664e-02 9.3528816199405160e-03 -2.7501515672069238e-02" << endl;
	outfile << "57 -1.7980703138665951e-02 -2.0041585094997282e-02 1.4030560841356713e-02" << endl;
	outfile << "58 -1.6968971591429328e-02 5.3087226038875824e-03 1.4873355151351859e-02" << endl;
	outfile << "" << endl;
	outfile << "Bonds" << endl;
	outfile << "" << endl;
	outfile << "1 1 5 6" << endl;
	outfile << "2 1 6 7" << endl;
	outfile << "3 1 7 8" << endl;
	outfile << "4 1 9 10" << endl;
	outfile << "5 1 8 9" << endl;
	outfile << "6 1 10 11" << endl;
	outfile << "7 1 13 14" << endl;
	outfile << "8 1 14 15" << endl;
	outfile << "9 1 11 12" << endl;
	outfile << "10 1 12 13" << endl;
	outfile << "11 1 15 16" << endl;
	outfile << "12 1 16 17" << endl;
	outfile << "13 1 17 18" << endl;
	outfile << "14 1 18 19" << endl;
	outfile << "15 1 19 20" << endl;
	outfile << "16 1 20 21" << endl;
	outfile << "17 1 22 23" << endl;
	outfile << "18 1 23 24" << endl;
	outfile << "19 1 21 22" << endl;
	outfile << "20 1 1 2" << endl;
	outfile << "21 1 3 4" << endl;
	outfile << "22 1 4 5" << endl;
	outfile << "23 1 2 3" << endl;
	outfile << "24 1 24 25" << endl;
	outfile << "25 1 25 26" << endl;
	outfile << "26 1 26 27" << endl;
	outfile << "27 1 27 28" << endl;
	outfile << "28 1 28 29" << endl;
	outfile << "29 1 30 31" << endl;
	outfile << "30 1 31 32" << endl;
	outfile << "31 1 29 30" << endl;
	outfile << "32 1 32 33" << endl;
	outfile << "33 1 33 34" << endl;
	outfile << "34 1 34 35" << endl;
	outfile << "35 1 35 36" << endl;
	outfile << "36 1 36 37" << endl;
	outfile << "37 1 37 38" << endl;
	outfile << "38 1 77 78" << endl;
	outfile << "39 1 82 83" << endl;
	outfile << "40 1 83 84" << endl;
	outfile << "41 1 81 82" << endl;
	outfile << "42 1 84 85" << endl;
	outfile << "43 1 80 81" << endl;
	outfile << "44 1 85 86" << endl;
	outfile << "45 1 87 88" << endl;
	outfile << "46 1 86 87" << endl;
	outfile << "47 1 79 80" << endl;
	outfile << "48 1 38 39" << endl;
	outfile << "49 1 46 47" << endl;
	outfile << "50 1 43 44" << endl;
	outfile << "51 1 45 46" << endl;
	outfile << "52 1 44 45" << endl;
	outfile << "53 1 68 69" << endl;
	outfile << "54 1 98 99" << endl;
	outfile << "55 1 99 100" << endl;
	outfile << "56 1 39 40" << endl;
	outfile << "57 1 40 41" << endl;
	outfile << "58 1 42 43" << endl;
	outfile << "59 1 41 42" << endl;
	outfile << "60 1 70 71" << endl;
	outfile << "61 1 69 70" << endl;
	outfile << "62 1 71 72" << endl;
	outfile << "63 1 73 74" << endl;
	outfile << "64 1 74 75" << endl;
	outfile << "65 1 72 73" << endl;
	outfile << "66 1 75 76" << endl;
	outfile << "67 1 76 77" << endl;
	outfile << "68 1 89 90" << endl;
	outfile << "69 1 88 89" << endl;
	outfile << "70 1 90 91" << endl;
	outfile << "71 1 91 92" << endl;
	outfile << "72 1 92 93" << endl;
	outfile << "73 1 97 98" << endl;
	outfile << "74 1 96 97" << endl;
	outfile << "75 1 94 95" << endl;
	outfile << "76 1 93 94" << endl;
	outfile << "77 1 95 96" << endl;
	outfile << "78 1 78 79" << endl;
	outfile << "79 1 52 53" << endl;
	outfile << "80 1 51 52" << endl;
	outfile << "81 1 53 54" << endl;
	outfile << "82 1 59 60" << endl;
	outfile << "83 1 60 61" << endl;
	outfile << "84 1 47 48" << endl;
	outfile << "85 1 49 50" << endl;
	outfile << "86 1 48 49" << endl;
	outfile << "87 1 50 51" << endl;
	outfile << "88 1 62 63" << endl;
	outfile << "89 1 63 64" << endl;
	outfile << "90 1 64 65" << endl;
	outfile << "91 1 65 66" << endl;
	outfile << "92 1 66 67" << endl;
	outfile << "93 1 67 68" << endl;
	outfile << "94 1 61 62" << endl;
	outfile << "95 1 54 55" << endl;
	outfile << "96 1 55 56" << endl;
	outfile << "97 1 56 57" << endl;
	outfile << "98 1 57 58" << endl;
	outfile << "99 1 58 59" << endl;
	outfile << "" << endl;
	outfile << "Angles" << endl;
	outfile << "" << endl;
	outfile << "1 1 4 5 6" << endl;
	outfile << "2 1 5 6 7" << endl;
	outfile << "3 1 6 7 8" << endl;
	outfile << "4 1 8 9 10" << endl;
	outfile << "5 1 7 8 9" << endl;
	outfile << "6 1 9 10 11" << endl;
	outfile << "7 1 12 13 14" << endl;
	outfile << "8 1 13 14 15" << endl;
	outfile << "9 1 10 11 12" << endl;
	outfile << "10 1 11 12 13" << endl;
	outfile << "11 1 14 15 16" << endl;
	outfile << "12 1 15 16 17" << endl;
	outfile << "13 1 16 17 18" << endl;
	outfile << "14 1 17 18 19" << endl;
	outfile << "15 1 18 19 20" << endl;
	outfile << "16 1 19 20 21" << endl;
	outfile << "17 1 21 22 23" << endl;
	outfile << "18 1 22 23 24" << endl;
	outfile << "19 1 20 21 22" << endl;
	outfile << "20 1 2 3 4" << endl;
	outfile << "21 1 3 4 5" << endl;
	outfile << "22 1 1 2 3" << endl;
	outfile << "23 1 23 24 25" << endl;
	outfile << "24 1 24 25 26" << endl;
	outfile << "25 1 25 26 27" << endl;
	outfile << "26 1 26 27 28" << endl;
	outfile << "27 1 27 28 29" << endl;
	outfile << "28 1 29 30 31" << endl;
	outfile << "29 1 30 31 32" << endl;
	outfile << "30 1 28 29 30" << endl;
	outfile << "31 1 31 32 33" << endl;
	outfile << "32 1 32 33 34" << endl;
	outfile << "33 1 33 34 35" << endl;
	outfile << "34 1 34 35 36" << endl;
	outfile << "35 1 35 36 37" << endl;
	outfile << "36 1 36 37 38" << endl;
	outfile << "37 1 76 77 78" << endl;
	outfile << "38 1 81 82 83" << endl;
	outfile << "39 1 82 83 84" << endl;
	outfile << "40 1 80 81 82" << endl;
	outfile << "41 1 83 84 85" << endl;
	outfile << "42 1 79 80 81" << endl;
	outfile << "43 1 84 85 86" << endl;
	outfile << "44 1 86 87 88" << endl;
	outfile << "45 1 85 86 87" << endl;
	outfile << "46 1 78 79 80" << endl;
	outfile << "47 1 37 38 39" << endl;
	outfile << "48 1 45 46 47" << endl;
	outfile << "49 1 42 43 44" << endl;
	outfile << "50 1 44 45 46" << endl;
	outfile << "51 1 43 44 45" << endl;
	outfile << "52 1 67 68 69" << endl;
	outfile << "53 1 97 98 99" << endl;
	outfile << "54 1 98 99 100" << endl;
	outfile << "55 1 38 39 40" << endl;
	outfile << "56 1 39 40 41" << endl;
	outfile << "57 1 41 42 43" << endl;
	outfile << "58 1 40 41 42" << endl;
	outfile << "59 1 69 70 71" << endl;
	outfile << "60 1 68 69 70" << endl;
	outfile << "61 1 70 71 72" << endl;
	outfile << "62 1 72 73 74" << endl;
	outfile << "63 1 73 74 75" << endl;
	outfile << "64 1 71 72 73" << endl;
	outfile << "65 1 74 75 76" << endl;
	outfile << "66 1 75 76 77" << endl;
	outfile << "67 1 88 89 90" << endl;
	outfile << "68 1 87 88 89" << endl;
	outfile << "69 1 89 90 91" << endl;
	outfile << "70 1 90 91 92" << endl;
	outfile << "71 1 91 92 93" << endl;
	outfile << "72 1 96 97 98" << endl;
	outfile << "73 1 95 96 97" << endl;
	outfile << "74 1 93 94 95" << endl;
	outfile << "75 1 92 93 94" << endl;
	outfile << "76 1 94 95 96" << endl;
	outfile << "77 1 77 78 79" << endl;
	outfile << "78 1 51 52 53" << endl;
	outfile << "79 1 50 51 52" << endl;
	outfile << "80 1 52 53 54" << endl;
	outfile << "81 1 58 59 60" << endl;
	outfile << "82 1 59 60 61" << endl;
	outfile << "83 1 46 47 48" << endl;
	outfile << "84 1 48 49 50" << endl;
	outfile << "85 1 47 48 49" << endl;
	outfile << "86 1 49 50 51" << endl;
	outfile << "87 1 61 62 63" << endl;
	outfile << "88 1 62 63 64" << endl;
	outfile << "89 1 63 64 65" << endl;
	outfile << "90 1 64 65 66" << endl;
	outfile << "91 1 65 66 67" << endl;
	outfile << "92 1 66 67 68" << endl;
	outfile << "93 1 60 61 62" << endl;
	outfile << "94 1 53 54 55" << endl;
	outfile << "95 1 54 55 56" << endl;
	outfile << "96 1 55 56 57" << endl;
	outfile << "97 1 56 57 58" << endl;
	outfile << "98 1 57 58 59" << endl;
	outfile << "" << endl;
	outfile << "Dihedrals" << endl;
	outfile << "" << endl;
	outfile << "1 1 4 5 6 7" << endl;
	outfile << "2 1 5 6 7 8" << endl;
	outfile << "3 1 6 7 8 9" << endl;
	outfile << "4 1 8 9 10 11" << endl;
	outfile << "5 1 7 8 9 10" << endl;
	outfile << "6 1 9 10 11 12" << endl;
	outfile << "7 1 12 13 14 15" << endl;
	outfile << "8 1 13 14 15 16" << endl;
	outfile << "9 1 10 11 12 13" << endl;
	outfile << "10 1 11 12 13 14" << endl;
	outfile << "11 1 14 15 16 17" << endl;
	outfile << "12 1 15 16 17 18" << endl;
	outfile << "13 1 16 17 18 19" << endl;
	outfile << "14 1 17 18 19 20" << endl;
	outfile << "15 1 18 19 20 21" << endl;
	outfile << "16 1 19 20 21 22" << endl;
	outfile << "17 1 21 22 23 24" << endl;
	outfile << "18 1 22 23 24 25" << endl;
	outfile << "19 1 20 21 22 23" << endl;
	outfile << "20 1 2 3 4 5" << endl;
	outfile << "21 1 3 4 5 6" << endl;
	outfile << "22 1 1 2 3 4" << endl;
	outfile << "23 1 23 24 25 26" << endl;
	outfile << "24 1 24 25 26 27" << endl;
	outfile << "25 1 25 26 27 28" << endl;
	outfile << "26 1 26 27 28 29" << endl;
	outfile << "27 1 27 28 29 30" << endl;
	outfile << "28 1 29 30 31 32" << endl;
	outfile << "29 1 30 31 32 33" << endl;
	outfile << "30 1 28 29 30 31" << endl;
	outfile << "31 1 31 32 33 34" << endl;
	outfile << "32 1 32 33 34 35" << endl;
	outfile << "33 1 33 34 35 36" << endl;
	outfile << "34 1 34 35 36 37" << endl;
	outfile << "35 1 35 36 37 38" << endl;
	outfile << "36 1 36 37 38 39" << endl;
	outfile << "37 1 76 77 78 79" << endl;
	outfile << "38 1 81 82 83 84" << endl;
	outfile << "39 1 82 83 84 85" << endl;
	outfile << "40 1 80 81 82 83" << endl;
	outfile << "41 1 83 84 85 86" << endl;
	outfile << "42 1 79 80 81 82" << endl;
	outfile << "43 1 84 85 86 87" << endl;
	outfile << "44 1 86 87 88 89" << endl;
	outfile << "45 1 85 86 87 88" << endl;
	outfile << "46 1 78 79 80 81" << endl;
	outfile << "47 1 37 38 39 40" << endl;
	outfile << "48 1 45 46 47 48" << endl;
	outfile << "49 1 42 43 44 45" << endl;
	outfile << "50 1 44 45 46 47" << endl;
	outfile << "51 1 43 44 45 46" << endl;
	outfile << "52 1 67 68 69 70" << endl;
	outfile << "53 1 97 98 99 100" << endl;
	outfile << "54 1 38 39 40 41" << endl;
	outfile << "55 1 39 40 41 42" << endl;
	outfile << "56 1 41 42 43 44" << endl;
	outfile << "57 1 40 41 42 43" << endl;
	outfile << "58 1 69 70 71 72" << endl;
	outfile << "59 1 68 69 70 71" << endl;
	outfile << "60 1 70 71 72 73" << endl;
	outfile << "61 1 72 73 74 75" << endl;
	outfile << "62 1 73 74 75 76" << endl;
	outfile << "63 1 71 72 73 74" << endl;
	outfile << "64 1 74 75 76 77" << endl;
	outfile << "65 1 75 76 77 78" << endl;
	outfile << "66 1 88 89 90 91" << endl;
	outfile << "67 1 87 88 89 90" << endl;
	outfile << "68 1 89 90 91 92" << endl;
	outfile << "69 1 90 91 92 93" << endl;
	outfile << "70 1 91 92 93 94" << endl;
	outfile << "71 1 96 97 98 99" << endl;
	outfile << "72 1 95 96 97 98" << endl;
	outfile << "73 1 93 94 95 96" << endl;
	outfile << "74 1 92 93 94 95" << endl;
	outfile << "75 1 94 95 96 97" << endl;
	outfile << "76 1 77 78 79 80" << endl;
	outfile << "77 1 51 52 53 54" << endl;
	outfile << "78 1 50 51 52 53" << endl;
	outfile << "79 1 52 53 54 55" << endl;
	outfile << "80 1 58 59 60 61" << endl;
	outfile << "81 1 59 60 61 62" << endl;
	outfile << "82 1 46 47 48 49" << endl;
	outfile << "83 1 48 49 50 51" << endl;
	outfile << "84 1 47 48 49 50" << endl;
	outfile << "85 1 49 50 51 52" << endl;
	outfile << "86 1 61 62 63 64" << endl;
	outfile << "87 1 62 63 64 65" << endl;
	outfile << "88 1 63 64 65 66" << endl;
	outfile << "89 1 64 65 66 67" << endl;
	outfile << "90 1 65 66 67 68" << endl;
	outfile << "91 1 66 67 68 69" << endl;
	outfile << "92 1 60 61 62 63" << endl;
	outfile << "93 1 53 54 55 56" << endl;
	outfile << "94 1 54 55 56 57" << endl;
	outfile << "95 1 55 56 57 58" << endl;
	outfile << "96 1 56 57 58 59" << endl;
	outfile << "97 1 57 58 59 60" << endl;
	outfile.close();
}


