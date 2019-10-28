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

#define BOOST_TEST_MODULE gaussian_test
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <votca/xtp/orbitals.h>
#include <votca/xtp/qmpackagefactory.h>
using namespace votca::xtp;
using namespace std;

BOOST_AUTO_TEST_SUITE(gaussian_test)

BOOST_AUTO_TEST_CASE(polar_test) {
  std::ofstream polar("polar_gaussian.log");
  polar << " Set GAUSS_MDEF to \"2gb\"." << std::endl;
  polar << "          MDV=     134217728 using IRadAn=       1." << std::endl;
  polar << "          Solving linear equations simultaneously, MaxMat=       0."
        << std::endl;
  polar << "          There are     3 degrees of freedom in the 1st order "
           "CPHF.  IDoFFX=0 NUNeed=     3."
        << std::endl;
  polar << "      3 vectors produced by pass  0 Test12= 2.00D-15 3.33D-08 "
           "XBig12= 5.14D+00 1.22D+00."
        << std::endl;
  polar << " AX will form     3 AO Fock derivatives at one time." << std::endl;
  polar << "      3 vectors produced by pass  1 Test12= 2.00D-15 3.33D-08 "
           "XBig12= 7.21D-02 1.18D-01."
        << std::endl;
  polar << "      3 vectors produced by pass  2 Test12= 2.00D-15 3.33D-08 "
           "XBig12= 7.97D-05 6.02D-03."
        << std::endl;
  polar << "      3 vectors produced by pass  3 Test12= 2.00D-15 3.33D-08 "
           "XBig12= 9.00D-08 2.03D-04."
        << std::endl;
  polar << "      3 vectors produced by pass  4 Test12= 2.00D-15 3.33D-08 "
           "XBig12= 2.65D-09 2.03D-05."
        << std::endl;
  polar << "      3 vectors produced by pass  5 Test12= 2.00D-15 3.33D-08 "
           "XBig12= 1.35D-11 1.46D-06."
        << std::endl;
  polar << "      1 vectors produced by pass  6 Test12= 2.00D-15 3.33D-08 "
           "XBig12= 1.71D-14 4.77D-08."
        << std::endl;
  polar << " InvSVY:  IOpt=1 It=  1 EMax= 4.44D-16" << std::endl;
  polar << " Solved reduced A of dimension    19 with     3 vectors."
        << std::endl;
  polar << " FullF1:  Do perturbations      1 to       3." << std::endl;
  polar << " SCF Polarizability for W=    0.000000:" << std::endl;
  polar << "                1             2             3 " << std::endl;
  polar << "      1  0.114169D+02" << std::endl;
  polar << "      2 -0.270043D-03  0.114172D+02" << std::endl;
  polar << "      3  0.830290D-04 -0.221681D-03  0.114175D+02" << std::endl;
  polar << " Isotropic polarizability for W=    0.000000       11.42 Bohr**3."
        << std::endl;
  polar << " SCF Static Hyperpolarizability:" << std::endl;
  polar << " K=  1 block:" << std::endl;
  polar << "                1 " << std::endl;
  polar << "      1 -0.285540D+02" << std::endl;
  polar << " K=  2 block:" << std::endl;
  polar << "                1             2 " << std::endl;
  polar << "      1 -0.238877D+01" << std::endl;
  polar << "      2  0.122480D+02 -0.680382D+01" << std::endl;
  polar << " K=  3 block:" << std::endl;
  polar << "                1             2             3 " << std::endl;
  polar << "      1  0.111937D+02" << std::endl;
  polar << "      2  0.323093D+01 -0.229864D+02" << std::endl;
  polar << "      3  0.163031D+02  0.919020D+01  0.117877D+02" << std::endl;
  polar << " End of Minotr F.D. properties file   721 does not exist."
        << std::endl;
  polar << " End of Minotr F.D. properties file   722 does not exist."
        << std::endl;
  polar << " End of Minotr F.D. properties file   788 does not exist."
        << std::endl;
  polar << " Leave Link 1002 at Fri Aug  9 23:54:14 2019, MaxMem=   134217728 "
           "cpu:               0.7 elap:               0.8"
        << std::endl;
  polar << " (Enter /curc/sw/gaussian/16/g16/l601.exe)" << std::endl;
  polar
      << " Copying SCF densities to generalized density rwf, IOpCl= 0 IROHF=0."
      << std::endl;
  polar << "" << std::endl;
  polar << " ******************************************************************"
           "****"
        << std::endl;
  polar << "" << std::endl;
  polar << "            Population analysis using the SCF density."
        << std::endl;
  polar << "" << std::endl;
  polar << " ******************************************************************"
           "****"
        << std::endl;
  polar << "" << std::endl;
  polar << " Alpha  occ. eigenvalues --  -10.14484  -0.71198  -0.40074  "
           "-0.40071  -0.40071"
        << std::endl;
  polar << " Alpha virt. eigenvalues --    0.16282   0.20865   0.20866   "
           "0.20867   0.72527"
        << std::endl;
  polar << " Alpha virt. eigenvalues --    0.72528   0.72529   1.06018   "
           "1.09140   1.09142"
        << std::endl;
  polar << " Alpha virt. eigenvalues --    1.09144   1.71498" << std::endl;
  polar << "          Condensed to atoms (all electrons):" << std::endl;
  polar << "               1          2          3          4          5"
        << std::endl;
  polar << "     1  H    0.491489   0.365636  -0.024611  -0.024619  -0.024618"
        << std::endl;
  polar << "     2  C    0.365636   5.404352   0.365627   0.365623   0.365625"
        << std::endl;
  polar << "     3  H   -0.024611   0.365627   0.491522  -0.024623  -0.024625"
        << std::endl;
  polar << "     4  H   -0.024619   0.365623  -0.024623   0.491532  -0.024623"
        << std::endl;
  polar << "     5  H   -0.024618   0.365625  -0.024625  -0.024623   0.491522"
        << std::endl;
  polar << " Mulliken charges with hydrogens summed into heavy atoms:"
        << std::endl;
  polar << "               1" << std::endl;
  polar << "     2  C   -0.000000" << std::endl;
  polar << " Electronic spatial extent (au):  <R**2>=             35.2117"
        << std::endl;
  polar << " Charge=             -0.0000 electrons" << std::endl;
  polar << " Dipole moment (field-independent basis, Debye):" << std::endl;
  polar << "    X=             -0.0001    Y=             -0.0000    Z=         "
           "    -0.0001  Tot=              0.0001"
        << std::endl;
  polar << " Quadrupole moment (field-independent basis, Debye-Ang):"
        << std::endl;
  polar << "   XX=             -8.2204   YY=             -8.2204   ZZ=         "
           "    -8.2201"
        << std::endl;
  polar << "   XY=             -0.0000   XZ=              0.0002   YZ=         "
           "    -0.0001"
        << std::endl;
  polar << " Traceless Quadrupole moment (field-independent basis, Debye-Ang):"
        << std::endl;
  polar << "   XX=             -0.0001   YY=             -0.0001   ZZ=         "
           "     0.0002"
        << std::endl;
  polar << "   XY=             -0.0000   XZ=              0.0002   YZ=         "
           "    -0.0001"
        << std::endl;
  polar << " Octapole moment (field-independent basis, Debye-Ang**2):"
        << std::endl;
  polar << "  XXX=             -0.8568  YYY=             -0.2040  ZZZ=         "
           "     0.3544  XYY=              0.3676"
        << std::endl;
  polar << "  XXY=             -0.0716  XXZ=              0.3361  XZZ=         "
           "     0.4895  YZZ=              0.2759"
        << std::endl;
  polar << "  YYZ=             -0.6898  XYZ=              0.0970" << std::endl;
  polar << " Hexadecapole moment (field-independent basis, Debye-Ang**3):"
        << std::endl;
  polar << " XXXX=            -14.1430 YYYY=            -14.3584 ZZZZ=         "
           "   -14.5865 XXXY=              0.0729"
        << std::endl;
  polar << " XXXZ=             -0.1404 YYYX=             -0.1812 YYYZ=         "
           "    -0.0147 ZZZX=              0.3923"
        << std::endl;
  polar << " ZZZY=             -0.0487 XXYY=             -5.1040 XXZZ=         "
           "    -4.8757 YYZZ=             -4.6606"
        << std::endl;
  polar << " XXYZ=              0.0636 YYXZ=             -0.2518 ZZXY=         "
           "     0.1084"
        << std::endl;
  polar << " N-N= 1.347284297650D+01 E-N=-1.194830087342D+02  KE= "
           "3.987780184382D+01"
        << std::endl;
  polar << "  Exact polarizability:  11.417  -0.000  11.417   0.000  -0.000  "
           "11.417"
        << std::endl;
  polar << " Approx polarizability:  12.844  -0.000  12.844   0.000  -0.000  "
           "12.845"
        << std::endl;
  polar << " No NMR shielding tensors so no spin-rotation constants."
        << std::endl;
  polar << " Leave Link  601 at Fri Aug  9 23:54:14 2019, MaxMem=   134217728 "
           "cpu:               0.1 elap:               0.1"
        << std::endl;
  polar << " (Enter /curc/sw/gaussian/16/g16/l9999.exe)" << std::endl;
  polar << "" << std::endl;
  polar << " ------------------------------------------------------------------"
           "----"
        << std::endl;
  polar << "" << std::endl;
  polar << " Electric dipole moment (input orientation):" << std::endl;
  polar << " (Debye = 10**-18 statcoulomb cm , SI units = C m)" << std::endl;
  polar << "                  (au)            (Debye)         (10**-30 SI)"
        << std::endl;
  polar << "   Tot        0.390098D-04      0.991531D-04      0.330739D-03"
        << std::endl;
  polar << "   x         -0.241725D-04     -0.614403D-04     -0.204943D-03"
        << std::endl;
  polar << "   y         -0.187198D-04     -0.475810D-04     -0.158713D-03"
        << std::endl;
  polar << "   z         -0.242286D-04     -0.615830D-04     -0.205419D-03"
        << std::endl;
  polar << "" << std::endl;
  polar << " Dipole polarizability, Alpha (input orientation)." << std::endl;
  polar << " (esu units = cm**3 , SI units = C**2 m**2 J**-1)" << std::endl;
  polar << " Alpha(0;0):" << std::endl;
  polar << "               (au)            (10**-24 esu)      (10**-40 SI)"
        << std::endl;
  polar << "   iso        0.114172D+02      0.169185D+01      0.188244D+01"
        << std::endl;
  polar << "   aniso      0.770584D-03      0.114189D-03      0.127052D-03"
        << std::endl;
  polar << "   xx         0.114169D+02      0.169182D+01      0.188240D+01"
        << std::endl;
  polar << "   yx        -0.270043D-03     -0.400162D-04     -0.445241D-04"
        << std::endl;
  polar << "   yy         0.114172D+02      0.169185D+01      0.188244D+01"
        << std::endl;
  polar << "   zx         0.830290D-04      0.123036D-04      0.136896D-04"
        << std::endl;
  polar << "   zy        -0.221681D-03     -0.328498D-04     -0.365503D-04"
        << std::endl;
  polar << "   zz         0.114175D+02      0.169189D+01      0.188249D+01"
        << std::endl;
  polar << "" << std::endl;
  polar << " First dipole hyperpolarizability, Beta (input orientation)."
        << std::endl;
  polar << " ||, _|_  parallel and perpendicular components, (z) with respect "
           "to z axis,"
        << std::endl;
  polar << " vector components x,y,z.  Values do not include the 1/n! factor "
           "of 1/2."
        << std::endl;
  polar << " (esu units = statvolt**-1 cm**4 , SI units = C**3 m**3 J**-2)"
        << std::endl;
  polar << " Beta(0;0,0):" << std::endl;
  polar << "               (au)            (10**-30 esu)      (10**-50 SI)"
        << std::endl;
  polar << "   || (z)    -0.300405D-02     -0.259527D-04     -0.963209D-05"
        << std::endl;
  polar << "   _|_(z)    -0.100135D-02     -0.865090D-05     -0.321070D-05"
        << std::endl;
  polar << "   x         -0.861972D-02     -0.744677D-04     -0.276379D-04"
        << std::endl;
  polar << "   y         -0.716324D-02     -0.618848D-04     -0.229679D-04"
        << std::endl;
  polar << "   z         -0.150203D-01     -0.129763D-03     -0.481604D-04"
        << std::endl;
  polar << "   ||         0.374817D-02      0.323813D-04      0.120180D-04"
        << std::endl;
  polar << "   xxx       -0.285540D+02     -0.246684D+00     -0.915545D-01"
        << std::endl;
  polar << "   xxy       -0.238877D+01     -0.206371D-01     -0.765925D-02"
        << std::endl;
  polar << "   yxy        0.122480D+02      0.105813D+00      0.392716D-01"
        << std::endl;
  polar << "   yyy       -0.680382D+01     -0.587797D-01     -0.218155D-01"
        << std::endl;
  polar << "   xxz        0.111937D+02      0.967050D-01      0.358911D-01"
        << std::endl;
  polar << "   yxz        0.323093D+01      0.279128D-01      0.103595D-01"
        << std::endl;
  polar << "   yyz       -0.229864D+02     -0.198585D+00     -0.737028D-01"
        << std::endl;
  polar << "   zxz        0.163031D+02      0.140846D+00      0.522736D-01"
        << std::endl;
  polar << "   zyz        0.919020D+01      0.793961D-01      0.294671D-01"
        << std::endl;
  polar << "   zzz        0.117877D+02      0.101837D+00      0.377956D-01"
        << std::endl;
  polar << "" << std::endl;
  polar << " ------------------------------------------------------------------"
           "----"
        << std::endl;
  polar << "" << std::endl;
  polar << " Dipole orientation:" << std::endl;
  polar << "     1         -0.38912045          0.77097464         -1.86366599"
        << std::endl;
  polar << "     6          0.00000000          0.00000000         -0.00000000"
        << std::endl;
  polar << "     1          0.80477824         -1.88022178         -0.19092896"
        << std::endl;
  polar << "     1          1.33171252          1.21830236          0.98062434"
        << std::endl;
  polar << "     1         -1.74724566         -0.10910912          1.07441313"
        << std::endl;
  polar << "" << std::endl;
  polar << " Electric dipole moment (dipole orientation):" << std::endl;
  polar << " (Debye = 10**-18 statcoulomb cm , SI units = C m)" << std::endl;
  polar << "                  (au)            (Debye)         (10**-30 SI)"
        << std::endl;
  polar << "   Tot        0.390098D-04      0.991531D-04      0.330739D-03"
        << std::endl;
  polar << "   x          0.000000D+00      0.000000D+00      0.000000D+00"
        << std::endl;
  polar << "   y          0.000000D+00      0.000000D+00      0.000000D+00"
        << std::endl;
  polar << "   z          0.390098D-04      0.991531D-04      0.330739D-03"
        << std::endl;
  polar << "" << std::endl;
  polar << " Dipole polarizability, Alpha (dipole orientation)." << std::endl;
  polar << " (esu units = cm**3 , SI units = C**2 m**2 J**-1)" << std::endl;
  polar << " Alpha(0;0):" << std::endl;
  polar << "               (au)            (10**-24 esu)      (10**-40 SI)"
        << std::endl;
  polar << "   iso        0.114172D+02      0.169185D+01      0.188244D+01"
        << std::endl;
  polar << "   aniso      0.770584D-03      0.114189D-03      0.127052D-03"
        << std::endl;
  polar << "   xx         0.114171D+02      0.169184D+01      0.188243D+01"
        << std::endl;
  polar << "   yx        -0.108358D-04     -0.160570D-05     -0.178658D-05"
        << std::endl;
  polar << "   yy         0.114175D+02      0.169190D+01      0.188249D+01"
        << std::endl;
  polar << "   zx         0.191662D-03      0.284014D-04      0.316009D-04"
        << std::endl;
  polar << "   zy        -0.281776D-03     -0.417548D-04     -0.464585D-04"
        << std::endl;
  polar << "   zz         0.114170D+02      0.169182D+01      0.188240D+01"
        << std::endl;
  polar << "" << std::endl;
  polar << " First dipole hyperpolarizability, Beta (dipole orientation)."
        << std::endl;
  polar << " ||, _|_  parallel and perpendicular components, (z) with respect "
           "to z axis,"
        << std::endl;
  polar << " vector components x,y,z.  Values do not include the 1/n! factor "
           "of 1/2."
        << std::endl;
  polar << " (esu units = statvolt**-1 cm**4 , SI units = C**3 m**3 J**-2)"
        << std::endl;
  polar << " Beta(0;0,0):" << std::endl;
  polar << "               (au)            (10**-30 esu)      (10**-50 SI)"
        << std::endl;
  polar << "   || (z)     0.362152D-02      0.312872D-04      0.116119D-04"
        << std::endl;
  polar << "   _|_(z)     0.120717D-02      0.104291D-04      0.387064D-05"
        << std::endl;
  polar << "   x          0.404320D-02      0.349301D-04      0.129639D-04"
        << std::endl;
  polar << "   y         -0.264327D-02     -0.228358D-04     -0.847526D-05"
        << std::endl;
  polar << "   z          0.181076D-01      0.156436D-03      0.580596D-04"
        << std::endl;
  polar << "   ||         0.374817D-02      0.323813D-04      0.120180D-04"
        << std::endl;
  polar << "   xxx       -0.103540D+02     -0.894502D-01     -0.331986D-01"
        << std::endl;
  polar << "   xxy        0.299662D+01      0.258885D-01      0.960826D-02"
        << std::endl;
  polar << "   yxy        0.188452D+02      0.162808D+00      0.604244D-01"
        << std::endl;
  polar << "   yyy       -0.180699D+02     -0.156110D+00     -0.579386D-01"
        << std::endl;
  polar << "   xxz        0.190252D+02      0.164363D+00      0.610018D-01"
        << std::endl;
  polar << "   yxz        0.109034D+02      0.941973D-01      0.349604D-01"
        << std::endl;
  polar << "   yyz       -0.129760D+01     -0.112102D-01     -0.416056D-02"
        << std::endl;
  polar << "   zxz       -0.848985D+01     -0.733457D-01     -0.272215D-01"
        << std::endl;
  polar << "   zyz        0.150724D+02      0.130214D+00      0.483276D-01"
        << std::endl;
  polar << "   zzz       -0.177216D+02     -0.153101D+00     -0.568218D-01"
        << std::endl;
  polar << "" << std::endl;
  polar << " ------------------------------------------------------------------"
           "----"
        << std::endl;
  polar << "" << std::endl;
  polar << " Test job not archived." << std::endl;
  polar << " 1\\1\\GINC-SHAS0345\\SP\\RPBE1PBE\\3-21G\\C1H4\\JOBR9774\\09-Aug-"
           "2019\\0\\\\#p po"
        << std::endl;
  polar << " lar PBE1PBE/3-21g punch=mo nosymm test\\\\CH4_chelpg "
           "PBE0\\\\0,1\\H,0,0.528"
        << std::endl;
  polar << " 8,0.161,0.9359\\C,0,0.,0.,0.\\H,0,0.2051,0.824,-0.6786\\H,0,0."
           "3345,-0.931"
        << std::endl;
  polar << " 4,-0.4496\\H,0,-1.0685,-0.0537,0.1921\\\\Version=ES64L-G16RevA."
           "03\\HF=-40."
        << std::endl;
  polar << " 2407948\\RMSD=2.129e-09\\Dipole=-0.0000242,-0.0000187,-0."
           "0000242\\Polar=1"
        << std::endl;
  polar << " 1.416942,-0.00027,11.4171738,0.000083,-0.0002217,11."
           "4174661\\HyperPolar"
        << std::endl;
  polar << " =-28.5540132,-2.3887674,12.2480366,-6.8038165,11.1937143,3."
           "2309343,-22"
        << std::endl;
  polar << " .9864197,16.3031034,9.1901961,11.7876986\\Quadrupole=-0.0000506,-"
           "0.0001"
        << std::endl;
  polar << " 1,0.0001606,-0.0000309,0.000121,-0.0000456\\PG=C01 [X(C1H4)]\\\\@"
        << std::endl;
  polar << "" << std::endl;
  polar << " Normal termination of Gaussian 16 at Fri Aug  9 23:54:15 2019."
        << std::endl;
  polar.close();

  QMPackageFactory::RegisterAll();
  std::unique_ptr<QMPackage> gaussian =
      std::unique_ptr<QMPackage>(QMPackages().Create("gaussian"));
  Logger log;
  gaussian->setLog(&log);
  gaussian->setRunDir(".");
  gaussian->setLogFileName("polar_gaussian.log");
  Eigen::Matrix3d polar_mat = gaussian->GetPolarizability();

  Eigen::Matrix3d polar_ref = Eigen::Matrix3d::Zero();

  polar << "   xx         0.114169D+02      0.169182D+01      0.188240D+01"
        << std::endl;
  polar << "   yx        -0.270043D-03     -0.400162D-04     -0.445241D-04"
        << std::endl;
  polar << "   yy         0.114172D+02      0.169185D+01      0.188244D+01"
        << std::endl;
  polar << "   zx         0.830290D-04      0.123036D-04      0.136896D-04"
        << std::endl;
  polar << "   zy        -0.221681D-03     -0.328498D-04     -0.365503D-04"
        << std::endl;
  polar << "   zz         0.114175D+02      0.169189D+01      0.188249D+01"
        << std::endl;

  polar_ref << 0.114169e+02, -0.270043e-03, 0.830290e-04, -0.270043e-03,
      0.114172e+02, -0.221681e-03, 0.830290e-04, -0.221681e-03, 0.114175e+02;

  bool polar_check = polar_ref.isApprox(polar_mat, 1e-5);
  if (!polar_check) {
    std::cout << "res" << std::endl;
    std::cout << polar_mat << std::endl;
    std::cout << "ref " << std::endl;
    std::cout << polar_ref << std::endl;
  }
}

BOOST_AUTO_TEST_CASE(read_withexternal) {

  std::ofstream ext("gaussian_ext.log");
  ext << " Entering Gaussian System, Link 0=/home/apps/gaussian_6_9/g09/g09"
      << std::endl;
  ext << " Initial command:" << std::endl;
  ext << " /home/apps/gaussian_6_9/g09/l1.exe /scratch/Gau-28211.inp "
         "-scrdir=/scratch/"
      << std::endl;
  ext << " ******************************************" << std::endl;
  ext << " Gaussian 09:  EM64L-G09RevA.01  8-May-2009" << std::endl;
  ext << "                30-Jul-2019 " << std::endl;
  ext << " ******************************************" << std::endl;
  ext << " %mem=1Gb" << std::endl;
  ext << " %nprocshared=1" << std::endl;
  ext << " Will use up to    1 processors via shared memory." << std::endl;
  ext << " --------------------------------------------------------"
      << std::endl;
  ext << " #p pop=minimal PBE1PBE/3-21g punch=mo nosymm test charge"
      << std::endl;
  ext << " --------------------------------------------------------"
      << std::endl;
  ext << " 1/18=10,28=-2,38=1,125=101/1;" << std::endl;
  ext << " 2/12=2,15=1,17=6,18=5,40=1/2;" << std::endl;
  ext << " 3/5=5,11=2,16=1,25=1,30=1,74=-13/1,2,3;" << std::endl;
  ext << " 4//1;" << std::endl;
  ext << " 5/5=2,38=5/2;" << std::endl;
  ext << " 6/7=2,8=2,9=2,10=2,28=1/1;" << std::endl;
  ext << " 99/5=1,9=1,10=32/99;" << std::endl;
  ext << " Leave Link    1 at Tue Jul 30 12:36:51 2019, MaxMem=  134217728 "
         "cpu:       0.0"
      << std::endl;
  ext << " (Enter /home/apps/gaussian_6_9/g09/l101.exe)" << std::endl;
  ext << " ------------" << std::endl;
  ext << " CH4_ext PBE0" << std::endl;
  ext << " ------------" << std::endl;
  ext << " Symbolic Z-matrix:" << std::endl;
  ext << " Charge =  0 Multiplicity = 1" << std::endl;
  ext << " H                     0.5288    0.161     0.9359 " << std::endl;
  ext << " C                     0.        0.        0. " << std::endl;
  ext << " H                     0.2051    0.824    -0.6786 " << std::endl;
  ext << " H                     0.3345   -0.9314   -0.4496 " << std::endl;
  ext << " H                    -1.0685   -0.0537    0.1921 " << std::endl;
  ext << " " << std::endl;
  ext << " NAtoms=      5 NQM=      5 NQMF=      0 NMic=      0 NMicF=      0 "
         "NTot=      5."
      << std::endl;
  ext << "                    Isotopes and Nuclear Properties:" << std::endl;
  ext << " (Nuclear quadrupole moments (NQMom) in fm**2, nuclear magnetic "
         "moments (NMagM)"
      << std::endl;
  ext << "  in nuclear magnetons)" << std::endl;
  ext << "" << std::endl;
  ext << "  Atom         1           2           3           4           5"
      << std::endl;
  ext << " IAtWgt=           1          12           1           1           1"
      << std::endl;
  ext << " AtmWgt=   1.0078250  12.0000000   1.0078250   1.0078250   1.0078250"
      << std::endl;
  ext << " NucSpn=           1           0           1           1           1"
      << std::endl;
  ext << " AtZEff=   0.0000000   0.0000000   0.0000000   0.0000000   0.0000000"
      << std::endl;
  ext << " NQMom=    0.0000000   0.0000000   0.0000000   0.0000000   0.0000000"
      << std::endl;
  ext << " NMagM=    2.7928460   0.0000000   2.7928460   2.7928460   2.7928460"
      << std::endl;
  ext << " Background charge distribution read from input stream:" << std::endl;
  ext << " Point Charges:" << std::endl;
  ext << " XYZ=    5.0000    0.0000    0.0000 Q=    1.0000 A=    0.0000 R=    "
         "0.0000 C=    0.0000"
      << std::endl;
  ext << " XYZ=   -5.0000    0.0000    0.0000 Q=   -1.0000 A=    0.0000 R=    "
         "0.0000 C=    0.0000"
      << std::endl;
  ext << " Sum of input charges=            0.000000" << std::endl;
  ext << "              Moments=           10.000000            0.000000       "
         "     0.000000"
      << std::endl;
  ext << " Leave Link  101 at Tue Jul 30 12:36:51 2019, MaxMem=  134217728 "
         "cpu:       0.1"
      << std::endl;
  ext << " (Enter /home/apps/gaussian_6_9/g09/l202.exe)" << std::endl;
  ext << "                          Input orientation:                         "
         " "
      << std::endl;
  ext << " --------------------------------------------------------------------"
         "-"
      << std::endl;
  ext << " Center     Atomic      Atomic             Coordinates (Angstroms)"
      << std::endl;
  ext << " Number     Number       Type             X           Y           Z"
      << std::endl;
  ext << " --------------------------------------------------------------------"
         "-"
      << std::endl;
  ext << "      1          1           0        0.528800    0.161000    "
         "0.935900"
      << std::endl;
  ext << "      2          6           0        0.000000    0.000000    "
         "0.000000"
      << std::endl;
  ext << "      3          1           0        0.205100    0.824000   "
         "-0.678600"
      << std::endl;
  ext << "      4          1           0        0.334500   -0.931400   "
         "-0.449600"
      << std::endl;
  ext << "      5          1           0       -1.068500   -0.053700    "
         "0.192100"
      << std::endl;
  ext << " --------------------------------------------------------------------"
         "-"
      << std::endl;
  ext << "                    Distance matrix (angstroms):" << std::endl;
  ext << "                    1          2          3          4          5"
      << std::endl;
  ext << "     1  H    0.000000" << std::endl;
  ext << "     2  C    1.086950   0.000000" << std::endl;
  ext << "     3  H    1.775095   1.086987   0.000000" << std::endl;
  ext << "     4  H    1.775021   1.086985   1.774997   0.000000" << std::endl;
  ext << "     5  H    1.775022   1.086958   1.774974   1.774978   0.000000"
      << std::endl;
  ext << " Symmetry turned off by external request." << std::endl;
  ext << " Stoichiometry    CH4" << std::endl;
  ext << " Framework group  C1[X(CH4)]" << std::endl;
  ext << " Deg. of freedom     9" << std::endl;
  ext << " Full point group                 C1      NOp   1" << std::endl;
  ext << " Rotational constants (GHZ):    159.1656622    159.1551086    "
         "159.1530548"
      << std::endl;
  ext << " Leave Link  202 at Tue Jul 30 12:36:52 2019, MaxMem=  134217728 "
         "cpu:       0.0"
      << std::endl;
  ext << " (Enter /home/apps/gaussian_6_9/g09/l301.exe)" << std::endl;
  ext << " Standard basis: 3-21G (6D, 7F)" << std::endl;
  ext << " Ernie: Thresh=  0.10000D-02 Tol=  0.10000D-05 Strict=F."
      << std::endl;
  ext << " Integral buffers will be    131072 words long." << std::endl;
  ext << " Raffenetti 2 integral format." << std::endl;
  ext << " Two-electron integral symmetry is turned off." << std::endl;
  ext << "    17 basis functions,    27 primitive gaussians,    17 cartesian "
         "basis functions"
      << std::endl;
  ext << "     5 alpha electrons        5 beta electrons" << std::endl;
  ext << "       nuclear repulsion energy        13.4728429172 Hartrees."
      << std::endl;
  ext << " IExCor= 1009 DFT=T Ex+Corr=PBE1PBE ExCW=0 ScaHFX=  0.250000"
      << std::endl;
  ext << " ScaDFX=  0.750000  0.750000  1.000000  1.000000 ScalE2=  1.000000  "
         "1.000000"
      << std::endl;
  ext << " IRadAn=      0 IRanWt=     -1 IRanGd=            0 ICorTp=0"
      << std::endl;
  ext << " NAtoms=    5 NActive=    5 NUniq=    5 SFac= 7.50D-01 NAtFMM=   80 "
         "NAOKFM=F Big=F"
      << std::endl;
  ext << " Background charge distribution read from rwf." << std::endl;
  ext << " FoFCou: FMM=F IPFlag=           0 FMFlag=      100000 FMFlg1=       "
         "    0"
      << std::endl;
  ext << "         NFxFlg=           0 DoJE=F BraDBF=T KetDBF=T FulRan=T"
      << std::endl;
  ext << "         Omega=  0.000000  0.000000  1.000000  0.000000  0.000000 "
         "ICntrl=    1100 IOpCl=  0"
      << std::endl;
  ext << "         NMat0=    1 NMatS0=    1 NMatT0=    0 NMatD0=    1 NMtDS0=  "
         "  0 NMtDT0=    0"
      << std::endl;
  ext << "         I1Cent=    10001006 NGrid=           7." << std::endl;
  ext << " Symmetry not used in FoFCou." << std::endl;
  ext << " Self energy of the charges =        -0.0529177209 a.u." << std::endl;
  ext << " MM   energy of the charges =         0.0000000000 a.u." << std::endl;
  ext << " Nuclei-charges interaction =        -0.0043517531 a.u." << std::endl;
  ext << " Nuclear repulsion after external point charges =       "
         "13.4155734432 Hartrees."
      << std::endl;
  ext << " Leave Link  301 at Tue Jul 30 12:36:52 2019, MaxMem=  134217728 "
         "cpu:       0.1"
      << std::endl;
  ext << " (Enter /home/apps/gaussian_6_9/g09/l302.exe)" << std::endl;
  ext << " NPDir=0 NMtPBC=     1 NCelOv=     1 NCel=       1 NClECP=     1 "
         "NCelD=      1"
      << std::endl;
  ext << "         NCelK=      1 NCelE2=     1 NClLst=     1 CellRange=     "
         "0.0."
      << std::endl;
  ext << " One-electron integrals computed using PRISM." << std::endl;
  ext << " NBasis=    17 RedAO= T  NBF=    17" << std::endl;
  ext << " NBsUse=    17 1.00D-06 NBFU=    17" << std::endl;
  ext << " Precomputing XC quadrature grid using" << std::endl;
  ext << " IXCGrd= 2 IRadAn=           0 IRanWt=          -1 IRanGd=           "
         "0 AccXCQ= 1.00D-10."
      << std::endl;
  ext << " NRdTot=     303 NPtTot=       38382 NUsed=       40581 NTot=       "
         "40613"
      << std::endl;
  ext << " NSgBfM=    17    17    17    17    17 NAtAll=     5     5."
      << std::endl;
  ext << " Leave Link  302 at Tue Jul 30 12:36:52 2019, MaxMem=  134217728 "
         "cpu:       0.3"
      << std::endl;
  ext << " (Enter /home/apps/gaussian_6_9/g09/l303.exe)" << std::endl;
  ext << " DipDrv:  MaxL=1." << std::endl;
  ext << " Leave Link  303 at Tue Jul 30 12:36:52 2019, MaxMem=  134217728 "
         "cpu:       0.0"
      << std::endl;
  ext << " (Enter /home/apps/gaussian_6_9/g09/l401.exe)" << std::endl;
  ext << " Harris functional with IExCor= 1009 diagonalized for initial guess."
      << std::endl;
  ext << " Harris En= -40.3536082458206    " << std::endl;
  ext << " Leave Link  401 at Tue Jul 30 12:36:52 2019, MaxMem=  134217728 "
         "cpu:       0.3"
      << std::endl;
  ext << " (Enter /home/apps/gaussian_6_9/g09/l502.exe)" << std::endl;
  ext << " Closed shell SCF:" << std::endl;
  ext << " Requested convergence on RMS density matrix=1.00D-08 within 128 "
         "cycles."
      << std::endl;
  ext << " Requested convergence on MAX density matrix=1.00D-06." << std::endl;
  ext << " Requested convergence on             energy=1.00D-06." << std::endl;
  ext << " No special actions if energy rises." << std::endl;
  ext << " Using DIIS extrapolation, IDIIS=  1040." << std::endl;
  ext << " Two-electron integral symmetry not used." << std::endl;
  ext << "        40523 words used for storage of precomputed grid."
      << std::endl;
  ext << " Keep R1 ints in memory in canonical form, NReq=929941." << std::endl;
  ext << " IEnd=       60254 IEndB=       60254 NGot=   134217728 MDV=   "
         "134162768"
      << std::endl;
  ext << " LenX=   134162768 LenY=   134161886" << std::endl;
  ext << " Symmetry not used in FoFDir." << std::endl;
  ext << " MinBra= 0 MaxBra= 1 Meth= 1." << std::endl;
  ext << " IRaf=       0 NMat=   1 IRICut=       1 DoRegI=T DoRafI=F ISym2E= 0 "
         "JSym2E=0."
      << std::endl;
  ext << " Integral accuracy reduced to 1.0D-05 until final iterations."
      << std::endl;
  ext << " SCF Done:  E(RPBE1PBE) =  -40.2974255693     A.U. after   10 cycles"
      << std::endl;
  ext << "             Convg  =    0.1387D-08             -V/T =  2.0105"
      << std::endl;
  ext << " KE= 3.987794385638D+01 PE=-1.194848521547D+02 EE= 2.589390928582D+01"
      << std::endl;
  ext << " Leave Link  502 at Tue Jul 30 12:36:54 2019, MaxMem=  134217728 "
         "cpu:       1.2"
      << std::endl;
  ext << " (Enter /home/apps/gaussian_6_9/g09/l601.exe)" << std::endl;
  ext << " Copying SCF densities to generalized density rwf, IOpCl= 0 IROHF=0."
      << std::endl;
  ext << "" << std::endl;
  ext << " ********************************************************************"
         "**"
      << std::endl;
  ext << "" << std::endl;
  ext << "            Population analysis using the SCF density." << std::endl;
  ext << "" << std::endl;
  ext << " ********************************************************************"
         "**"
      << std::endl;
  ext << "" << std::endl;
  ext << " Alpha  occ. eigenvalues --  -10.14431  -0.71244  -0.40545  -0.40329 "
         " -0.39207"
      << std::endl;
  ext << " Alpha virt. eigenvalues --    0.15786   0.19739   0.20238   0.23990 "
         "  0.71815"
      << std::endl;
  ext << " Alpha virt. eigenvalues --    0.72761   0.72921   1.05788   1.08164 "
         "  1.08590"
      << std::endl;
  ext << " Alpha virt. eigenvalues --    1.11584   1.71628" << std::endl;
  ext << " N-N= 1.341557344318D+01 E-N=-1.194848521052D+02  KE= "
         "3.987794385638D+01"
      << std::endl;
  ext << " No NMR shielding tensors so no spin-rotation constants."
      << std::endl;
  ext << " Leave Link  601 at Tue Jul 30 12:36:54 2019, MaxMem=  134217728 "
         "cpu:       0.1"
      << std::endl;
  ext << " (Enter /home/apps/gaussian_6_9/g09/l9999.exe)" << std::endl;
  ext << "" << std::endl;
  ext << " Test job not archived." << std::endl;
  ext << " 1\\1\\GINC-N41\\SP\\RPBE1PBE\\3-21G\\C1H4\\SIM_JOBR9774\\30-Jul-"
         "2019\\0\\\\#p pop"
      << std::endl;
  ext << " =minimal PBE1PBE/3-21g punch=mo nosymm test charge\\\\CH4_ext "
         "PBE0\\\\0,1\\"
      << std::endl;
  ext << " H,0,0.5288,0.161,0.9359\\C,0,0.,0.,0.\\H,0,0.2051,0.824,-0.6786\\H,"
         "0,0.33"
      << std::endl;
  ext << " 45,-0.9314,-0.4496\\H,0,-1.0685,-0.0537,0.1921\\\\Version=EM64L-"
         "G09RevA.0"
      << std::endl;
  ext << " 1\\HF=-40.2974256\\RMSD=1.387e-09\\Dipole=-0.2538894,-0.0007213,-0."
         "000152"
      << std::endl;
  ext << " 5\\Quadrupole=0.2196626,-0.0932575,-0.1264051,0.0180435,-0.0874911,-"
         "0.0"
      << std::endl;
  ext << " 249764\\PG=C01 [X(C1H4)]\\\\@" << std::endl;
  ext << "" << std::endl;
  ext << " Job cpu time:  0 days  0 hours  0 minutes  2.4 seconds."
      << std::endl;
  ext << " Normal termination of Gaussian 09 at Tue Jul 30 12:36:54 2019."
      << std::endl;
  ext.close();

  QMPackageFactory::RegisterAll();
  std::unique_ptr<QMPackage> gaussian =
      std::unique_ptr<QMPackage>(QMPackages().Create("gaussian"));
  Logger log;
  Orbitals orb;
  gaussian->setLog(&log);
  gaussian->setRunDir(".");
  gaussian->setLogFileName("gaussian_ext.log");
  gaussian->ParseLogFile(orb);
  const QMMolecule& seg = orb.QMAtoms();
  double ang2bohr = votca::tools::conv::ang2bohr;
  QMMolecule ref("ref", 0);
  Eigen::Vector3d pos1 = {0.528800, 0.161000, 0.935900};
  QMAtom s1(0, "H", pos1 * ang2bohr);
  pos1 = {0.000000, 0.000000, 0.000000};
  QMAtom s2(1, "C", pos1 * ang2bohr);
  pos1 = {0.205100, 0.824000, -0.678600};
  QMAtom s3(2, "H", pos1 * ang2bohr);
  pos1 = {0.334500, -0.931400, -0.449600};
  QMAtom s4(3, "H", pos1 * ang2bohr);
  pos1 = {-1.068500, -0.053700, 0.192100};
  QMAtom s5(4, "H", pos1 * ang2bohr);
  ref.push_back(s1);
  ref.push_back(s2);
  ref.push_back(s3);
  ref.push_back(s4);
  ref.push_back(s5);
  BOOST_CHECK_CLOSE(orb.getScaHFX(), 0.25, 1e-5);

  BOOST_CHECK_EQUAL(seg.size(), ref.size());
  for (Index i = 0; i < seg.size(); i++) {
    bool check_pos = ref[i].getPos().isApprox(seg[i].getPos(), 1e-5);
    BOOST_CHECK_EQUAL(check_pos, true);
    if (!check_pos) {
      std::cout << "result " << i << std::endl;
      std::cout << seg[i].getPos().transpose() << std::endl;
      std::cout << "ref " << i << std::endl;
      std::cout << ref[i].getPos().transpose() << std::endl;
    }
    BOOST_CHECK_EQUAL(ref[i].getElement(), seg[i].getElement());
  }
  double ref_tot =
      -40.2974256 - (-0.0529177209);  // HF - Self energy ext charges
  BOOST_CHECK_CLOSE(orb.getDFTTotalEnergy(), ref_tot, 1e-5);
  std::ofstream basisfile("3-21G_small.xml");
  basisfile << "<basis name=\"3-21G\">" << endl;
  basisfile << "  <element name=\"H\">" << endl;
  basisfile << "    <shell scale=\"1.0\" type=\"S\">" << endl;
  basisfile << "      <constant decay=\"5.447178e+00\">" << endl;
  basisfile << "        <contractions factor=\"1.562850e-01\" type=\"S\"/>"
            << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "      <constant decay=\"8.245470e-01\">" << endl;
  basisfile << "        <contractions factor=\"9.046910e-01\" type=\"S\"/>"
            << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "    </shell>" << endl;
  basisfile << "    <shell scale=\"1.0\" type=\"S\">" << endl;
  basisfile << "      <constant decay=\"1.831920e-01\">" << endl;
  basisfile << "        <contractions factor=\"1.000000e+00\" type=\"S\"/>"
            << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "    </shell>" << endl;
  basisfile << "  </element>" << endl;
  basisfile << "  <element name=\"C\">" << endl;
  basisfile << "    <shell scale=\"1.0\" type=\"S\">" << endl;
  basisfile << "      <constant decay=\"1.722560e+02\">" << endl;
  basisfile << "        <contractions factor=\"6.176690e-02\" type=\"S\"/>"
            << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "      <constant decay=\"2.591090e+01\">" << endl;
  basisfile << "        <contractions factor=\"3.587940e-01\" type=\"S\"/>"
            << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "      <constant decay=\"5.533350e+00\">" << endl;
  basisfile << "        <contractions factor=\"7.007130e-01\" type=\"S\"/>"
            << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "    </shell>" << endl;
  basisfile << "    <shell scale=\"1.0\" type=\"SP\">" << endl;
  basisfile << "      <constant decay=\"3.664980e+00\">" << endl;
  basisfile << "        <contractions factor=\"-3.958970e-01\" type=\"S\"/>"
            << endl;
  basisfile << "        <contractions factor=\"2.364600e-01\" type=\"P\"/>"
            << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "      <constant decay=\"7.705450e-01\">" << endl;
  basisfile << "        <contractions factor=\"1.215840e+00\" type=\"S\"/>"
            << endl;
  basisfile << "        <contractions factor=\"8.606190e-01\" type=\"P\"/>"
            << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "    </shell>" << endl;
  basisfile << "    <shell scale=\"1.0\" type=\"SP\">" << endl;
  basisfile << "      <constant decay=\"1.958570e-01\">" << endl;
  basisfile << "        <contractions factor=\"1.000000e+00\" type=\"S\"/>"
            << endl;
  basisfile << "        <contractions factor=\"1.000000e+00\" type=\"P\"/>"
            << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "    </shell>" << endl;
  basisfile << "  </element>" << endl;
  basisfile << "</basis>" << endl;
  basisfile.close();

  std::ofstream fort("gaussian_ext_mos.fort");
  fort << "(5D15.8)" << std::endl;
  fort << "    1 Alpha MO OE=-0.10144311D+02" << std::endl;
  fort << "-0.94219105D-03 0.15611639D-01 0.98598140D+00 0.10850707D+00 "
          "0.31392292D-03"
       << std::endl;
  fort << "-0.15169346D-06-0.12143039D-05-0.80313324D-01-0.20760220D-03 "
          "0.77748626D-06"
       << std::endl;
  fort << " 0.33748015D-05-0.92812281D-03 0.15567021D-01-0.93315546D-03 "
          "0.15584960D-01"
       << std::endl;
  fort << "-0.87355191D-03 0.15379376D-01" << std::endl;
  fort << "    2 Alpha MO OE=-0.71244389D+00" << std::endl;
  fort << " 0.12613249D+00 0.24272106D-01-0.21183035D+00 0.20143373D+00 "
          "0.13332601D-01"
       << std::endl;
  fort << "-0.42946479D-04 0.11443711D-03 0.62514672D+00 "
          "0.22841992D-01-0.10645715D-04"
       << std::endl;
  fort << "-0.31122632D-03 0.12291206D+00 0.20288994D-01 0.12417714D+00 "
          "0.21753987D-01"
       << std::endl;
  fort << " 0.11119334D+00 0.86877638D-02" << std::endl;
  fort << "    3 Alpha MO OE=-0.40545269D+00" << std::endl;
  fort << " 0.25928463D+00 0.26968603D+00 0.28167933D-02-0.24307550D-02 "
          "0.71720725D-01"
       << std::endl;
  fort << " 0.12381960D+00 0.35036144D+00-0.17190646D-01 0.53920816D-01 "
          "0.92758986D-01"
       << std::endl;
  fort << " 0.26222208D+00-0.82776446D-01-0.79058005D-01-0.16908962D+00-0."
          "16689111D+00"
       << std::endl;
  fort << "-0.11989174D-01-0.67834275D-02" << std::endl;
  fort << "    4 Alpha MO OE=-0.40329176D+00" << std::endl;
  fort << "-0.51220638D-01-0.53815473D-01-0.19692920D-02 "
          "0.17389758D-02-0.50613524D-01"
       << std::endl;
  fort << " 0.35842910D+00-0.11625696D+00 0.11866250D-01-0.39387570D-01 "
          "0.27838461D+00"
       << std::endl;
  fort << "-0.90190602D-01 0.24389876D+00 "
          "0.23442493D+00-0.19866182D+00-0.19756859D+00"
       << std::endl;
  fort << " 0.92300353D-02 0.53129869D-02" << std::endl;
  fort << "    5 Alpha MO OE=-0.39207040D+00" << std::endl;
  fort << " 0.71865134D-01 0.80273518D-01 0.14054719D-01-0.12466442D-01 "
          "0.37746750D+00"
       << std::endl;
  fort << " 0.25470180D-01-0.86715554D-01-0.85943143D-01 0.34931272D+00 "
          "0.23944257D-01"
       << std::endl;
  fort << "-0.81790302D-01 0.93479632D-01 0.93470861D-01 0.83683230D-01 "
          "0.87223023D-01"
       << std::endl;
  fort << "-0.27151150D+00-0.17516474D+00" << std::endl;
  fort << "    6 Alpha MO OE= 0.15785693D+00" << std::endl;
  fort << "-0.26081477D-01-0.12019517D+01-0.17647168D+00 0.75568862D-01 "
          "0.11526700D+00"
       << std::endl;
  fort << "-0.19881171D-02 0.16711564D-01 0.23409194D+01 "
          "0.38240117D+00-0.63814349D-02"
       << std::endl;
  fort << " 0.60071236D-01-0.19271229D-01-0.92644675D+00-0.21389762D-01-0."
          "10180868D+01"
       << std::endl;
  fort << "-0.46172422D-02-0.39204127D+00" << std::endl;
  fort << "    7 Alpha MO OE= 0.19739292D+00" << std::endl;
  fort << " 0.97951471D-01 0.14950488D+01-0.15547299D-01 "
          "0.93622041D-02-0.46010731D-01"
       << std::endl;
  fort << "-0.12137815D+00-0.30794625D+00 "
          "0.20171164D+00-0.18486834D+00-0.48603330D+00"
       << std::endl;
  fort << "-0.12346877D+01-0.31084876D-01-0.54906587D+00-0.71875519D-01-0."
          "11949437D+01"
       << std::endl;
  fort << " 0.32357460D-03-0.60272726D-01" << std::endl;
  fort << "    8 Alpha MO OE= 0.20237932D+00" << std::endl;
  fort << "-0.23956631D-01-0.31459567D+00 0.99432488D-02-0.62155803D-02 "
          "0.34937275D-01"
       << std::endl;
  fort << "-0.30843744D+00 0.11451188D+00-0.12862346D+00 "
          "0.14212665D+00-0.12523101D+01"
       << std::endl;
  fort << " 0.46553785D+00 0.10642221D+00 "
          "0.16035912D+01-0.79702138D-01-0.11345203D+01"
       << std::endl;
  fort << " 0.13082230D-04 0.43066766D-01" << std::endl;
  fort << "    9 Alpha MO OE= 0.23990368D+00" << std::endl;
  fort << "-0.53090158D-01-0.26263860D+00 0.51369664D-01-0.41646206D-01 "
          "0.28015341D+00"
       << std::endl;
  fort << " 0.16229376D-01-0.57168698D-01-0.65204622D+00 0.12779616D+01 "
          "0.71605295D-01"
       << std::endl;
  fort << "-0.25180915D+00-0.56901387D-01-0.31357213D+00-0.55153722D-01-0."
          "29106502D+00"
       << std::endl;
  fort << " 0.17285938D+00 0.18866764D+01" << std::endl;
  fort << "   10 Alpha MO OE= 0.71815361D+00" << std::endl;
  fort << "-0.19384154D+00-0.46903444D-01-0.89715286D-02-0.81232373D-02-0."
          "79759482D+00"
       << std::endl;
  fort << "-0.54348878D-01 0.18529582D+00 0.10540291D+00 0.12813025D+01 "
          "0.87035466D-01"
       << std::endl;
  fort << "-0.29627897D+00-0.24049804D+00-0.54987099D-01-0.21924396D+00-0."
          "51324066D-01"
       << std::endl;
  fort << " 0.54376136D+00 0.64184874D-01" << std::endl;
  fort << "   11 Alpha MO OE= 0.72760999D+00" << std::endl;
  fort << " 0.11450024D+00 0.18190414D-01 0.12736017D-02 0.11020663D-02 "
          "0.10368982D+00"
       << std::endl;
  fort << "-0.73687827D+00 0.23200870D+00-0.14638325D-01-0.17568831D+00 "
          "0.12563721D+01"
       << std::endl;
  fort << "-0.39581679D+00-0.53308095D+00-0.71884014D-01 0.44919032D+00 "
          "0.64559993D-01"
       << std::endl;
  fort << "-0.14586852D-01 0.10892251D-02" << std::endl;
  fort << "   12 Alpha MO OE= 0.72920758D+00" << std::endl;
  fort << "-0.58741428D+00-0.81827538D-01-0.18418108D-02-0.15965091D-02-0."
          "14886884D+00"
       << std::endl;
  fort << "-0.24630094D+00-0.71580126D+00 0.20945683D-01 0.25469979D+00 "
          "0.42423663D+00"
       << std::endl;
  fort << " 0.12335484D+01 0.18060215D+00 0.21243247D-01 0.36463708D+00 "
          "0.45784892D-01"
       << std::endl;
  fort << " 0.18643071D-01-0.18799114D-02" << std::endl;
  fort << "   13 Alpha MO OE= 0.10578832D+01" << std::endl;
  fort << " 0.79541904D+00-0.54613802D+00 "
          "0.11975888D+00-0.16164347D+00-0.22335450D+00"
       << std::endl;
  fort << " 0.59259401D-02-0.47439053D-01-0.14349828D+00 "
          "0.29719264D+00-0.71098835D-02"
       << std::endl;
  fort << " 0.57437031D-01 0.64973974D+00-0.30410734D+00 "
          "0.69704415D+00-0.38349810D+00"
       << std::endl;
  fort << " 0.43680119D+00 0.78901278D-01" << std::endl;
  fort << "   14 Alpha MO OE= 0.10816366D+01" << std::endl;
  fort << "-0.81572475D+00 0.13443208D+01 0.10770186D-01-0.89581534D-02 "
          "0.12353705D+00"
       << std::endl;
  fort << " 0.28872343D+00 "
          "0.74805805D+00-0.20608035D-01-0.14934314D+00-0.35438692D+00"
       << std::endl;
  fort << "-0.91799654D+00 0.32329465D+00-0.43875395D+00 "
          "0.66741740D+00-0.97837580D+00"
       << std::endl;
  fort << " 0.65430269D-01-0.31238599D-01" << std::endl;
  fort << "   15 Alpha MO OE= 0.10859019D+01" << std::endl;
  fort << " 0.16097408D+00-0.29337911D+00-0.66685359D-02 "
          "0.49433532D-02-0.92739437D-01"
       << std::endl;
  fort << " 0.75429672D+00-0.27121740D+00 0.13630485D-01 "
          "0.11278072D+00-0.92861199D+00"
       << std::endl;
  fort << " 0.33380473D+00-0.89531494D+00 0.13592800D+01 "
          "0.62996054D+00-0.10280768D+01"
       << std::endl;
  fort << "-0.45282790D-01 0.26497861D-01" << std::endl;
  fort << "   16 Alpha MO OE= 0.11158366D+01" << std::endl;
  fort << " 0.13376575D+00-0.38647497D+00-0.26632675D-01-0.16503724D-02-0."
          "72590849D+00"
       << std::endl;
  fort << "-0.44796048D-01 0.15612316D+00 0.81926043D-01 0.91131389D+00 "
          "0.56677662D-01"
       << std::endl;
  fort << "-0.19775176D+00 0.18462554D+00-0.47039976D+00 "
          "0.16212333D+00-0.43357818D+00"
       << std::endl;
  fort << "-0.11126856D+01 0.15526004D+01" << std::endl;
  fort << "   17 Alpha MO OE= 0.17162770D+01" << std::endl;
  fort << "-0.15420845D+00-0.76284712D+00-0.16433265D-01-0.19974513D+01 "
          "0.18392824D-01"
       << std::endl;
  fort << " 0.12658087D-03-0.30262457D-03 "
          "0.39472871D+01-0.42942443D-01-0.39733891D-03"
       << std::endl;
  fort << " 0.72253262D-03-0.15146383D+00-0.77424534D+00-0.15256700D+00-0."
          "76986091D+00"
       << std::endl;
  fort << "-0.13562257D+00-0.83084497D+00" << std::endl;
  fort.close();
  gaussian->setMOsFileName("gaussian_ext_mos.fort");
  orb.setDFTbasisName("3-21G_small.xml");
  gaussian->ParseMOsFile(orb);
  Eigen::VectorXd MOs_energy_ref = Eigen::VectorXd::Zero(17);
  MOs_energy_ref << -10.1443, -0.712444, -0.405453, -0.403292, -0.39207,
      0.157857, 0.197393, 0.202379, 0.239904, 0.718154, 0.72761, 0.729208,
      1.05788, 1.08164, 1.0859, 1.11584, 1.71628;
  bool check_eng = MOs_energy_ref.isApprox(orb.MOs().eigenvalues(), 1e-5);
  BOOST_CHECK_EQUAL(check_eng, true);
  if (!check_eng) {
    std::cout << "result eng" << std::endl;
    std::cout << orb.MOs().eigenvalues() << std::endl;
    std::cout << "ref eng" << std::endl;
    std::cout << MOs_energy_ref << std::endl;
  }
  Eigen::MatrixXd MOs_coeff_ref = Eigen::MatrixXd::Zero(17, 17);
  MOs_coeff_ref << -0.000942191, 0.126132, 0.259285, -0.0512206, 0.0718651,
      -0.0260815, 0.0979515, -0.0239566, -0.0530902, -0.193842, 0.1145,
      -0.587414, 0.795419, -0.815725, 0.160974, 0.133766, -0.154208, 0.0156116,
      0.0242721, 0.269686, -0.0538155, 0.0802735, -1.20195, 1.49505, -0.314596,
      -0.262639, -0.0469034, 0.0181904, -0.0818275, -0.546138, 1.34432,
      -0.293379, -0.386475, -0.762847, 0.985981, -0.21183, 0.00281679,
      -0.00196929, 0.0140547, -0.176472, -0.0155473, 0.00994325, 0.0513697,
      -0.00897153, 0.0012736, -0.00184181, 0.119759, 0.0107702, -0.00666854,
      -0.0266327, -0.0164333, 0.108507, 0.201434, -0.00243076, 0.00173898,
      -0.0124664, 0.0755689, 0.0093622, -0.00621558, -0.0416462, -0.00812324,
      0.00110207, -0.00159651, -0.161643, -0.00895815, 0.00494335, -0.00165037,
      -1.99745, -1.2143e-06, 0.000114437, 0.350361, -0.116257, -0.0867156,
      0.0167116, -0.307946, 0.114512, -0.0571687, 0.185296, 0.232009, -0.715801,
      -0.0474391, 0.748058, -0.271217, 0.156123, -0.000302625, -1.51693e-07,
      -4.29465e-05, 0.12382, 0.358429, 0.0254702, -0.00198812, -0.121378,
      -0.308437, 0.0162294, -0.0543489, -0.736878, -0.246301, 0.00592594,
      0.288723, 0.754297, -0.044796, 0.000126581, 0.000313923, 0.0133326,
      0.0717207, -0.0506135, 0.377468, 0.115267, -0.0460107, 0.0349373,
      0.280153, -0.797595, 0.10369, -0.148869, -0.223355, 0.123537, -0.0927394,
      -0.725908, 0.0183928, -0.0803133, 0.625147, -0.0171906, 0.0118663,
      -0.0859431, 2.34092, 0.201712, -0.128623, -0.652046, 0.105403, -0.0146383,
      0.0209457, -0.143498, -0.020608, 0.0136305, 0.081926, 3.94729, 3.3748e-06,
      -0.000311226, 0.262222, -0.0901906, -0.0817903, 0.0600712, -1.23469,
      0.465538, -0.251809, -0.296279, -0.395817, 1.23355, 0.057437, -0.917997,
      0.333805, -0.197752, 0.000722533, 7.77486e-07, -1.06457e-05, 0.092759,
      0.278385, 0.0239443, -0.00638143, -0.486033, -1.25231, 0.0716053,
      0.0870355, 1.25637, 0.424237, -0.00710988, -0.354387, -0.928612,
      0.0566777, -0.000397339, -0.000207602, 0.022842, 0.0539208, -0.0393876,
      0.349313, 0.382401, -0.184868, 0.142127, 1.27796, 1.2813, -0.175688,
      0.2547, 0.297193, -0.149343, 0.112781, 0.911314, -0.0429424, -0.000928123,
      0.122912, -0.0827764, 0.243899, 0.0934796, -0.0192712, -0.0310849,
      0.106422, -0.0569014, -0.240498, -0.533081, 0.180602, 0.64974, 0.323295,
      -0.895315, 0.184626, -0.151464, 0.015567, 0.020289, -0.079058, 0.234425,
      0.0934709, -0.926447, -0.549066, 1.60359, -0.313572, -0.0549871,
      -0.071884, 0.0212432, -0.304107, -0.438754, 1.35928, -0.4704, -0.774245,
      -0.000933155, 0.124177, -0.16909, -0.198662, 0.0836832, -0.0213898,
      -0.0718755, -0.0797021, -0.0551537, -0.219244, 0.44919, 0.364637,
      0.697044, 0.667417, 0.629961, 0.162123, -0.152567, 0.015585, 0.021754,
      -0.166891, -0.197569, 0.087223, -1.01809, -1.19494, -1.13452, -0.291065,
      -0.0513241, 0.06456, 0.0457849, -0.383498, -0.978376, -1.02808, -0.433578,
      -0.769861, -0.000873552, 0.111193, -0.0119892, 0.00923004, -0.271512,
      -0.00461724, 0.000323575, 1.30822e-05, 0.172859, 0.543761, -0.0145869,
      0.0186431, 0.436801, 0.0654303, -0.0452828, -1.11269, -0.135623,
      0.0153794, 0.00868776, -0.00678343, 0.00531299, -0.175165, -0.392041,
      -0.0602727, 0.0430668, 1.88668, 0.0641849, 0.00108923, -0.00187991,
      0.0789013, -0.0312386, 0.0264979, 1.5526, -0.830845;
  bool check_coeff = MOs_coeff_ref.isApprox(orb.MOs().eigenvectors(), 1e-5);
  BOOST_CHECK_EQUAL(check_coeff, true);
  if (!check_coeff) {
    std::cout << "result coeff" << std::endl;
    std::cout << orb.MOs().eigenvectors() << std::endl;
    std::cout << "ref coeff" << std::endl;
    std::cout << MOs_coeff_ref << std::endl;
  }
}

BOOST_AUTO_TEST_CASE(getcharges) {

  std::ofstream chelpg("gaussian_charges.log");
  chelpg << " Entering Gaussian System, Link 0=/home/apps/gaussian_6_9/g09/g09"
         << std::endl;
  chelpg << " This is part of the Gaussian(R) 09 program.  It is based on"
         << std::endl;
  chelpg << "  " << std::endl;
  chelpg << "" << std::endl;
  chelpg << " Cite this work as:" << std::endl;
  chelpg << " Gaussian 09, Revision A.01," << std::endl;
  chelpg << " ******************************************" << std::endl;
  chelpg << " Gaussian 09:  EM64L-G09RevA.01  8-May-2009" << std::endl;
  chelpg << "                30-Jul-2019 " << std::endl;
  chelpg << " ******************************************" << std::endl;
  chelpg << " %mem=1Gb" << std::endl;
  chelpg << " %nprocshared=1" << std::endl;
  chelpg << " Will use up to    1 processors via shared memory." << std::endl;
  chelpg << " ------------------------------------------------" << std::endl;
  chelpg << " #p pop=CHELPG PBE1PBE/3-21g punch=mo nosymm test" << std::endl;
  chelpg << " ------------------------------------------------" << std::endl;
  chelpg << " 1/38=1/1;" << std::endl;
  chelpg << " 2/12=2,15=1,17=6,18=5,40=1/2;" << std::endl;
  chelpg << " 3/5=5,11=2,16=1,25=1,30=1,74=-13/1,2,3;" << std::endl;
  chelpg << " 4//1;" << std::endl;
  chelpg << " 5/5=2,38=5/2;" << std::endl;
  chelpg << " 6/7=2,8=2,9=2,10=2,15=8,20=3,28=1/1,2;" << std::endl;
  chelpg << " 99/5=1,9=1,10=32/99;" << std::endl;
  chelpg << " Leave Link    1 at Tue Jul 30 12:35:31 2019, MaxMem=  134217728 "
            "cpu:       0.0"
         << std::endl;
  chelpg << " (Enter /home/apps/gaussian_6_9/g09/l101.exe)" << std::endl;
  chelpg << " ---------------" << std::endl;
  chelpg << " CH4_chelpg PBE0" << std::endl;
  chelpg << " ---------------" << std::endl;
  chelpg << " Symbolic Z-matrix:" << std::endl;
  chelpg << " Charge =  0 Multiplicity = 1" << std::endl;
  chelpg << " H                     0.5288    0.161     0.9359 " << std::endl;
  chelpg << " C                     0.        0.        0. " << std::endl;
  chelpg << " H                     0.2051    0.824    -0.6786 " << std::endl;
  chelpg << " H                     0.3345   -0.9314   -0.4496 " << std::endl;
  chelpg << " H                    -1.0685   -0.0537    0.1921 " << std::endl;
  chelpg << " Symmetry turned off by external request." << std::endl;
  chelpg << " Stoichiometry    CH4" << std::endl;
  chelpg << " Framework group  C1[X(CH4)]" << std::endl;
  chelpg << " Deg. of freedom     9" << std::endl;
  chelpg << " Full point group                 C1      NOp   1" << std::endl;
  chelpg << " Rotational constants (GHZ):    159.1656622    159.1551086    "
            "159.1530548"
         << std::endl;
  chelpg << " Leave Link  202 at Tue Jul 30 12:35:31 2019, MaxMem=  134217728 "
            "cpu:       0.0"
         << std::endl;
  chelpg << " (Enter /home/apps/gaussian_6_9/g09/l301.exe)" << std::endl;
  chelpg << " Standard basis: 3-21G (6D, 7F)" << std::endl;
  chelpg << " Ernie: Thresh=  0.10000D-02 Tol=  0.10000D-05 Strict=F."
         << std::endl;
  chelpg << " Integral buffers will be    131072 words long." << std::endl;
  chelpg << " Raffenetti 2 integral format." << std::endl;
  chelpg << " Two-electron integral symmetry is turned off." << std::endl;
  chelpg << "    17 basis functions,    27 primitive gaussians,    17 "
            "cartesian basis functions"
         << std::endl;
  chelpg << "     5 alpha electrons        5 beta electrons" << std::endl;
  chelpg << "       nuclear repulsion energy        13.4728429172 Hartrees."
         << std::endl;
  chelpg << " Symmetry not used in FoFCou." << std::endl;
  chelpg << " Harris En= -40.2961308964506    " << std::endl;
  chelpg << " Leave Link  401 at Tue Jul 30 12:35:32 2019, MaxMem=  134217728 "
            "cpu:       0.3"
         << std::endl;
  chelpg << " (Enter /home/apps/gaussian_6_9/g09/l502.exe)" << std::endl;
  chelpg << " Closed shell SCF:" << std::endl;
  chelpg << " Requested convergence on RMS density matrix=1.00D-08 within 128 "
            "cycles."
         << std::endl;
  chelpg << " Requested convergence on MAX density matrix=1.00D-06."
         << std::endl;
  chelpg << " Requested convergence on             energy=1.00D-06."
         << std::endl;
  chelpg << " No special actions if energy rises." << std::endl;
  chelpg << " Using DIIS extrapolation, IDIIS=  1040." << std::endl;
  chelpg << " Two-electron integral symmetry not used." << std::endl;
  chelpg << "        40523 words used for storage of precomputed grid."
         << std::endl;
  chelpg << " Keep R1 ints in memory in canonical form, NReq=929941."
         << std::endl;
  chelpg << " IEnd=       60254 IEndB=       60254 NGot=   134217728 MDV=   "
            "134162768"
         << std::endl;
  chelpg << " LenX=   134162768 LenY=   134161886" << std::endl;
  chelpg << " Symmetry not used in FoFDir." << std::endl;
  chelpg << " MinBra= 0 MaxBra= 1 Meth= 1." << std::endl;
  chelpg << " IRaf=       0 NMat=   1 IRICut=       1 DoRegI=T DoRafI=F "
            "ISym2E= 0 JSym2E=0."
         << std::endl;
  chelpg << " Integral accuracy reduced to 1.0D-05 until final iterations."
         << std::endl;
  chelpg << "" << std::endl;
  chelpg
      << " SCF Done:  E(RPBE1PBE) =  -40.2407953770     A.U. after    9 cycles"
      << std::endl;
  chelpg << "             Convg  =    0.1115D-08             -V/T =  2.0091"
         << std::endl;
  chelpg
      << " KE= 3.987780150847D+01 PE=-1.194830079285D+02 EE= 2.589156812589D+01"
      << std::endl;
  chelpg << " Leave Link  502 at Tue Jul 30 12:35:33 2019, MaxMem=  134217728 "
            "cpu:       1.2"
         << std::endl;
  chelpg << " (Enter /home/apps/gaussian_6_9/g09/l601.exe)" << std::endl;
  chelpg
      << " Copying SCF densities to generalized density rwf, IOpCl= 0 IROHF=0."
      << std::endl;
  chelpg << "" << std::endl;
  chelpg << " *****************************************************************"
            "*****"
         << std::endl;
  chelpg << "" << std::endl;
  chelpg << "            Population analysis using the SCF density."
         << std::endl;
  chelpg << "" << std::endl;
  chelpg << " *****************************************************************"
            "*****"
         << std::endl;
  chelpg << "" << std::endl;
  chelpg << " Alpha  occ. eigenvalues --  -10.14484  -0.71198  -0.40074  "
            "-0.40071  -0.40071"
         << std::endl;
  chelpg << " Alpha virt. eigenvalues --    0.16282   0.20865   0.20866   "
            "0.20867   0.72527"
         << std::endl;
  chelpg << " Alpha virt. eigenvalues --    0.72528   0.72529   1.06018   "
            "1.09140   1.09142"
         << std::endl;
  chelpg << " Alpha virt. eigenvalues --    1.09144   1.71498" << std::endl;
  chelpg << "          Condensed to atoms (all electrons):" << std::endl;
  chelpg << "              1          2          3          4          5"
         << std::endl;
  chelpg << "     1  H    0.491488   0.365637  -0.024611  -0.024619  -0.024618"
         << std::endl;
  chelpg << "     2  C    0.365637   5.404354   0.365627   0.365623   0.365625"
         << std::endl;
  chelpg << "     3  H   -0.024611   0.365627   0.491521  -0.024623  -0.024625"
         << std::endl;
  chelpg << "     4  H   -0.024619   0.365623  -0.024623   0.491532  -0.024624"
         << std::endl;
  chelpg << "     5  H   -0.024618   0.365625  -0.024625  -0.024624   0.491521"
         << std::endl;
  chelpg << " Mulliken atomic charges:" << std::endl;
  chelpg << "              1" << std::endl;
  chelpg << "     1  H    0.216724" << std::endl;
  chelpg << "     2  C   -0.866866" << std::endl;
  chelpg << "     3  H    0.216711" << std::endl;
  chelpg << "     4  H    0.216711" << std::endl;
  chelpg << "     5  H    0.216720" << std::endl;
  chelpg << " Sum of Mulliken atomic charges =   0.00000" << std::endl;
  chelpg << " Mulliken charges with hydrogens summed into heavy atoms:"
         << std::endl;
  chelpg << "              1" << std::endl;
  chelpg << "     2  C    0.000000" << std::endl;
  chelpg << " Sum of Mulliken charges with hydrogens summed into heavy atoms = "
            "  0.00000"
         << std::endl;
  chelpg << " Electronic spatial extent (au):  <R**2>=             35.2117"
         << std::endl;
  chelpg << " Charge=              0.0000 electrons" << std::endl;
  chelpg << " Leave Link  601 at Tue Jul 30 12:35:33 2019, MaxMem=  134217728 "
            "cpu:       0.1"
         << std::endl;
  chelpg << " (Enter /home/apps/gaussian_6_9/g09/l602.exe)" << std::endl;
  chelpg << " Breneman (CHELPG) radii used." << std::endl;
  chelpg << " Generate Potential Derived Charges using the Breneman model, "
            "NDens= 1."
         << std::endl;
  chelpg << " Grid spacing= 0.300 Box extension= 2.800" << std::endl;
  chelpg
      << " NStep X,Y,Z=   25     26     26   Total possible points=       16900"
      << std::endl;
  chelpg << " Number of Points to Fit=    5365" << std::endl;
  chelpg << "" << std::endl;
  chelpg << " *****************************************************************"
            "*****"
         << std::endl;
  chelpg << "" << std::endl;
  chelpg << "            Electrostatic Properties Using The SCF Density"
         << std::endl;
  chelpg << "" << std::endl;
  chelpg << " *****************************************************************"
            "*****"
         << std::endl;
  chelpg << "" << std::endl;
  chelpg << "       Atomic Center    1 is at   0.528800  0.161000  0.935900"
         << std::endl;
  chelpg << "       Atomic Center    2 is at   0.000000  0.000000  0.000000"
         << std::endl;
  chelpg << "       Atomic Center    3 is at   0.205100  0.824000 -0.678600"
         << std::endl;
  chelpg << "       Atomic Center    4 is at   0.334500 -0.931400 -0.449600"
         << std::endl;
  chelpg << "       Atomic Center    5 is at  -1.068500 -0.053700  0.192100"
         << std::endl;
  chelpg << "    5365 points will be used for fitting atomic charges"
         << std::endl;
  chelpg << " Fitting point charges to electrostatic potential" << std::endl;
  chelpg << " Charges from ESP fit, RMS=   0.00087 RRMS=   0.35760:"
         << std::endl;
  chelpg << " Charge=   0.00000 Dipole=    -0.0002     0.0016     0.0002 Tot=  "
            "   0.0016"
         << std::endl;
  chelpg << "              1" << std::endl;
  chelpg << "     1  H    0.132498" << std::endl;
  chelpg << "     2  C   -0.529733" << std::endl;
  chelpg << "     3  H    0.132585" << std::endl;
  chelpg << "     4  H    0.132207" << std::endl;
  chelpg << "     5  H    0.132443" << std::endl;
  chelpg << " -----------------------------------------------------------------"
         << std::endl;
  chelpg << "" << std::endl;
  chelpg << "              Electrostatic Properties (Atomic Units)"
         << std::endl;
  chelpg << "" << std::endl;
  chelpg << " -----------------------------------------------------------------"
         << std::endl;
  chelpg << "    Center     Electric         -------- Electric Field --------"
         << std::endl;
  chelpg << "               Potential          X             Y             Z"
         << std::endl;
  chelpg << " -----------------------------------------------------------------"
         << std::endl;
  chelpg << "    1 Atom     -1.111171" << std::endl;
  chelpg << "    2 Atom    -14.682105" << std::endl;
  chelpg << "    3 Atom     -1.111172" << std::endl;
  chelpg << "    4 Atom     -1.111174" << std::endl;
  chelpg << "    5 Atom     -1.111176" << std::endl;
  chelpg << " -----------------------------------------------------------------"
         << std::endl;
  chelpg << " Leave Link  602 at Tue Jul 30 12:35:33 2019, MaxMem=  134217728 "
            "cpu:       0.1"
         << std::endl;
  chelpg << " (Enter /home/apps/gaussian_6_9/g09/l9999.exe)" << std::endl;
  chelpg << "" << std::endl;
  chelpg << " Test job not archived." << std::endl;
  chelpg << " 1\\1\\GINC-N41\\SP\\RPBE1PBE\\3-21G\\C1H4\\SIM_JOBR9774\\30-Jul-"
            "2019\\0\\\\#p pop"
         << std::endl;
  chelpg << " =CHELPG PBE1PBE/3-21g punch=mo nosymm test\\\\CH4_chelpg "
            "PBE0\\\\0,1\\H,0,0"
         << std::endl;
  chelpg << " .5288,0.161,0.9359\\C,0,0.,0.,0.\\H,0,0.2051,0.824,-0.6786\\H,0,"
            "0.3345,-0"
         << std::endl;
  chelpg << " .9314,-0.4496\\H,0,-1.0685,-0.0537,0.1921\\\\Version=EM64L-"
            "G09RevA.01\\HF="
         << std::endl;
  chelpg << " -40.2407954\\RMSD=1.115e-09\\Dipole=-0.0000227,-0.0000173,-0."
            "0000229\\Qua"
         << std::endl;
  chelpg << " drupole=-0.00005,-0.0001103,0.0001603,-0.0000307,0.0001211,-0."
            "0000455\\"
         << std::endl;
  chelpg << " PG=C01 [X(C1H4)]\\\\@" << std::endl;
  chelpg << " Normal termination of Gaussian 09 at Tue Jul 30 12:35:34 2019."
         << std::endl;
  chelpg.close();

  QMPackageFactory::RegisterAll();
  std::unique_ptr<QMPackage> gaussian =
      std::unique_ptr<QMPackage>(QMPackages().Create("gaussian"));
  Logger log;
  gaussian->setLog(&log);
  gaussian->setRunDir(".");
  gaussian->setLogFileName("gaussian_charges.log");
  StaticSegment seg = gaussian->GetCharges();
  double ang2bohr = votca::tools::conv::ang2bohr;
  StaticSegment ref("ref", 0);
  Eigen::Vector3d pos1 = {0.528800, 0.161000, 0.935900};
  StaticSite s1(0, "H", pos1 * ang2bohr);
  s1.setCharge(0.132498);
  pos1 = {0.000000, 0.000000, 0.000000};
  StaticSite s2(1, "C", pos1 * ang2bohr);
  s2.setCharge(-0.529733);
  pos1 = {0.205100, 0.824000, -0.678600};
  StaticSite s3(2, "H", pos1 * ang2bohr);
  s3.setCharge(0.132585);
  pos1 = {0.334500, -0.931400, -0.449600};
  StaticSite s4(3, "H", pos1 * ang2bohr);
  s4.setCharge(0.132207);
  pos1 = {-1.068500, -0.053700, 0.192100};
  StaticSite s5(4, "H", pos1 * ang2bohr);
  s5.setCharge(0.132443);
  ref.push_back(s1);
  ref.push_back(s2);
  ref.push_back(s3);
  ref.push_back(s4);
  ref.push_back(s5);

  BOOST_CHECK_EQUAL(seg.size(), ref.size());
  for (Index i = 0; i < seg.size(); i++) {
    BOOST_CHECK_EQUAL(ref[i].Q().isApprox(seg[i].Q(), 1e-5), true);
    bool check_pos = ref[i].getPos().isApprox(seg[i].getPos(), 1e-5);
    BOOST_CHECK_EQUAL(check_pos, true);
    if (!check_pos) {
      std::cout << "result " << i << std::endl;
      std::cout << seg[i].getPos().transpose() << std::endl;
      std::cout << "ref " << i << std::endl;
      std::cout << ref[i].getPos().transpose() << std::endl;
    }
    BOOST_CHECK_EQUAL(ref[i].getElement(), seg[i].getElement());
  }
}

BOOST_AUTO_TEST_CASE(read_optgeometry) {

  std::ofstream opt("gaussian_opt.log");

  opt << " Entering Gaussian System, Link 0=/home/apps/gaussian_6_9/g09/g09"
      << std::endl;
  opt << " Initial command:" << std::endl;
  opt << " /home/apps/gaussian_6_9/g09/l1.exe /scratch/Gau-28293.inp "
         "-scrdir=/scratch/"
      << std::endl;
  opt << " Entering Link 1 = /home/apps/gaussian_6_9/g09/l1.exe PID=     28294."
      << std::endl;
  opt << "  " << std::endl;
  opt << " ******************************************" << std::endl;
  opt << " Gaussian 09:  EM64L-G09RevA.01  8-May-2009" << std::endl;
  opt << "                30-Jul-2019 " << std::endl;
  opt << " ******************************************" << std::endl;
  opt << " %mem=1Gb" << std::endl;
  opt << " %nprocshared=1" << std::endl;
  opt << " Will use up to    1 processors via shared memory." << std::endl;
  opt << " -----------------------------------------------------" << std::endl;
  opt << " #p opt pop=minimal PBE1PBE/3-21g punch=mo nosymm test" << std::endl;
  opt << " -----------------------------------------------------" << std::endl;
  opt << " Ernie: Thresh=  0.10000D-02 Tol=  0.10000D-05 Strict=F."
      << std::endl;
  opt << " Integral buffers will be    131072 words long." << std::endl;
  opt << " N-N= 1.339011150672D+01 E-N=-1.192924236182D+02  KE= "
         "3.984200803650D+01"
      << std::endl;
  opt << " No NMR shielding tensors so no spin-rotation constants."
      << std::endl;
  opt << " Leave Link  601 at Tue Jul 30 12:37:58 2019, MaxMem=  134217728 "
         "cpu:       0.1"
      << std::endl;
  opt << " (Enter /home/apps/gaussian_6_9/g09/l9999.exe)" << std::endl;
  opt << "" << std::endl;
  opt << " Test job not archived." << std::endl;
  opt << " 1\\1\\GINC-N41\\FOpt\\RPBE1PBE\\3-21G\\C1H4\\SIM_JOBR9774\\30-Jul-"
         "2019\\0\\\\#p o"
      << std::endl;
  opt << " pt pop=minimal PBE1PBE/3-21g punch=mo nosymm test\\\\CH4_opt "
         "PBE0\\\\0,1\\H"
      << std::endl;
  opt << " ,0.5320667034,0.1620160329,0.9416501581\\C,-0.0000199957,-0."
         "0000124141,"
      << std::endl;
  opt << " -0.0000290191\\H,0.2063990471,0.8290973723,-0.6827572201\\H,0."
         "3365791865"
      << std::endl;
  opt << " ,-0.9371634348,-0.4523731748\\H,-1.0751249412,-0.0540375562,0."
         "193309255"
      << std::endl;
  opt << " 9\\\\Version=EM64L-G09RevA.01\\HF=-40.240904\\RMSD=9.625e-10\\RMSF="
         "8.181e-0"
      << std::endl;
  opt << " 5\\Dipole=0.0000023,-0.0000026,-0.0000049\\Quadrupole=0.000016,-0."
         "000036"
      << std::endl;
  opt << " 8,0.0000208,0.0000019,0.0000228,-0.0000136\\PG=C01 [X(C1H4)]\\\\@"
      << std::endl;
  opt << "" << std::endl;
  opt << "" << std::endl;
  opt << " EVERYTHING'S GOT A MORAL, IF ONLY YOU CAN FIND IT." << std::endl;
  opt << "" << std::endl;
  opt << "                                    -- LEWIS CARROL, ALICE IN "
         "WONDERLAND"
      << std::endl;
  opt << " Job cpu time:  0 days  0 hours  0 minutes  5.4 seconds."
      << std::endl;
  opt << " File lengths (MBytes):  RWF=      5 Int=      0 D2E=      0 Chk=    "
         "  1 Scr=      1"
      << std::endl;
  opt << " Normal termination of Gaussian 09 at Tue Jul 30 12:37:58 2019."
      << std::endl;

  opt.close();
  QMPackageFactory::RegisterAll();
  std::unique_ptr<QMPackage> gaussian =
      std::unique_ptr<QMPackage>(QMPackages().Create("gaussian"));
  Logger log;
  gaussian->setLog(&log);
  gaussian->setRunDir(".");
  gaussian->setLogFileName("gaussian_opt.log");
  Orbitals orb;
  gaussian->ParseLogFile(orb);

  double ref_tot = -40.240904;  // HF - Self energy ext charges
  BOOST_CHECK_CLOSE(ref_tot, orb.getDFTTotalEnergy(), 1e-5);
  const QMMolecule& seg = orb.QMAtoms();
  double ang2bohr = votca::tools::conv::ang2bohr;
  QMMolecule ref("ref", 0);
  Eigen::Vector3d pos1 = {0.5320667034, 0.1620160329, 0.9416501581};
  QMAtom s1(0, "H", pos1 * ang2bohr);
  pos1 = {-0.0000199957, -0.0000124141, -0.0000290191};
  QMAtom s2(1, "C", pos1 * ang2bohr);
  pos1 = {0.2063990471, 0.8290973723, -0.6827572201};
  QMAtom s3(2, "H", pos1 * ang2bohr);
  pos1 = {0.3365791865, -0.9371634348, -0.4523731748};
  QMAtom s4(3, "H", pos1 * ang2bohr);
  pos1 = {-1.0751249412, -0.0540375562, 0.1933092559};
  QMAtom s5(4, "H", pos1 * ang2bohr);
  ref.push_back(s1);
  ref.push_back(s2);
  ref.push_back(s3);
  ref.push_back(s4);
  ref.push_back(s5);
  BOOST_CHECK_EQUAL(seg.size(), ref.size());
  for (Index i = 0; i < seg.size(); i++) {
    BOOST_CHECK_EQUAL(ref[i].getPos().isApprox(seg[i].getPos(), 1e-5), true);
    BOOST_CHECK_EQUAL(ref[i].getElement(), seg[i].getElement());
  }
}

BOOST_AUTO_TEST_SUITE_END()
