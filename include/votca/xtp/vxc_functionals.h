/*
 *            Copyright 2009-2019 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#pragma once
#ifndef VOTCA_XTP_VXCFUNCTIONALS_H
#define VOTCA_XTP_VXCFUNCTIONALS_H

#include <map>
#include <string>

namespace votca {
namespace xtp {

/**
    \brief conversion of functional string into integer



*/
class Vxc_Functionals {
 public:
  Vxc_Functionals() { FillMaps(); };
  ~Vxc_Functionals(){};

  const int &getID(std::string name) const {
    if (_stringtoID.count(name) == 0) {
      throw std::runtime_error("Functional " + name + " not supported");
    }
    return _stringtoID.at(name);
  }

 private:
  std::map<std::string, int> _stringtoID;

  inline void FillMaps() { FillstringtoID(); }

  inline void FillstringtoID() {

    _stringtoID["XC_LDA_X"] = 1;
    _stringtoID["XC_LDA_C_WIGNER"] = 2;
    _stringtoID["XC_LDA_C_RPA"] = 3;
    _stringtoID["XC_LDA_C_HL"] = 4;
    _stringtoID["XC_LDA_C_GL"] = 5;
    _stringtoID["XC_LDA_C_XALPHA"] = 6;
    _stringtoID["XC_LDA_C_VWN"] = 7;
    _stringtoID["XC_LDA_C_VWN_RPA"] = 8;
    _stringtoID["XC_LDA_C_PZ"] = 9;
    _stringtoID["XC_LDA_C_PZ_MOD"] = 10;
    _stringtoID["XC_LDA_C_OB_PZ"] = 11;
    _stringtoID["XC_LDA_C_PW"] = 12;
    _stringtoID["XC_LDA_C_PW_MOD"] = 13;
    _stringtoID["XC_LDA_C_OB_PW"] = 14;
    _stringtoID["XC_LDA_C_2D_AMGB"] = 15;
    _stringtoID["XC_LDA_C_2D_PRM"] = 16;
    _stringtoID["XC_LDA_C_vBH"] = 17;
    _stringtoID["XC_LDA_C_1D_CSC"] = 18;
    _stringtoID["XC_LDA_X_2D"] = 19;
    _stringtoID["XC_LDA_XC_TETER93"] = 20;
    _stringtoID["XC_LDA_X_1D"] = 21;
    _stringtoID["XC_LDA_C_ML1"] = 22;
    _stringtoID["XC_LDA_C_ML2"] = 23;
    _stringtoID["XC_LDA_C_GOMBAS"] = 24;
    _stringtoID["XC_LDA_C_PW_RPA"] = 25;
    _stringtoID["XC_LDA_C_1D_LOOS"] = 26;
    _stringtoID["XC_LDA_C_RC04"] = 27;
    _stringtoID["XC_LDA_C_VWN_1"] = 28;
    _stringtoID["XC_LDA_C_VWN_2"] = 29;
    _stringtoID["XC_LDA_C_VWN_3"] = 30;
    _stringtoID["XC_LDA_C_VWN_4"] = 31;
    _stringtoID["XC_LDA_K_TF"] = 50;
    _stringtoID["XC_LDA_K_LP"] = 51;
    _stringtoID["XC_GGA_C_Q2D"] = 47;
    _stringtoID["XC_GGA_X_Q2D"] = 48;
    _stringtoID["XC_GGA_X_PBE_MOL"] = 49;
    _stringtoID["XC_GGA_K_TFVW"] = 52;
    _stringtoID["XC_GGA_K_REVAPBEINT"] = 53;
    _stringtoID["XC_GGA_K_APBEINT"] = 54;
    _stringtoID["XC_GGA_K_REVAPBE"] = 55;
    _stringtoID["XC_GGA_X_AK13"] = 56;
    _stringtoID["XC_GGA_K_MEYER"] = 57;
    _stringtoID["XC_GGA_X_LV_RPW86"] = 58;
    _stringtoID["XC_GGA_X_PBE_TCA"] = 59;
    _stringtoID["XC_GGA_X_PBEINT"] = 60;
    _stringtoID["XC_GGA_C_ZPBEINT"] = 61;
    _stringtoID["XC_GGA_C_PBEINT"] = 62;
    _stringtoID["XC_GGA_C_ZPBESOL"] = 63;
    _stringtoID["XC_GGA_XC_OPBE_D"] = 65;
    _stringtoID["XC_GGA_XC_OPWLYP_D"] = 66;
    _stringtoID["XC_GGA_XC_OBLYP_D"] = 67;
    _stringtoID["XC_GGA_X_VMT84_GE"] = 68;
    _stringtoID["XC_GGA_X_VMT84_PBE"] = 69;
    _stringtoID["XC_GGA_X_VMT_GE"] = 70;
    _stringtoID["XC_GGA_X_VMT_PBE"] = 71;
    _stringtoID["XC_GGA_C_N12_SX"] = 79;
    _stringtoID["XC_GGA_C_N12"] = 80;
    _stringtoID["XC_GGA_X_N12"] = 82;
    _stringtoID["XC_GGA_C_VPBE"] = 83;
    _stringtoID["XC_GGA_C_OP_XALPHA"] = 84;
    _stringtoID["XC_GGA_C_OP_G96"] = 85;
    _stringtoID["XC_GGA_C_OP_PBE"] = 86;
    _stringtoID["XC_GGA_C_OP_B88"] = 87;
    _stringtoID["XC_GGA_C_FT97"] = 88;
    _stringtoID["XC_GGA_C_SPBE"] = 89;
    _stringtoID["XC_GGA_X_SSB_SW"] = 90;
    _stringtoID["XC_GGA_X_SSB"] = 91;
    _stringtoID["XC_GGA_X_SSB_D"] = 92;
    _stringtoID["XC_GGA_XC_HCTH_407P"] = 93;
    _stringtoID["XC_GGA_XC_HCTH_P76"] = 94;
    _stringtoID["XC_GGA_XC_HCTH_P14"] = 95;
    _stringtoID["XC_GGA_XC_B97_GGA1"] = 96;
    _stringtoID["XC_GGA_XC_HCTH_A"] = 97;
    _stringtoID["XC_GGA_X_BPCCAC"] = 98;
    _stringtoID["XC_GGA_C_REVTCA"] = 99;
    _stringtoID["XC_GGA_C_TCA"] = 100;
    _stringtoID["XC_GGA_X_PBE"] = 101;
    _stringtoID["XC_GGA_X_PBE_R"] = 102;
    _stringtoID["XC_GGA_X_B86"] = 103;
    _stringtoID["XC_GGA_X_HERMAN"] = 104;
    _stringtoID["XC_GGA_X_B86_MGC"] = 105;
    _stringtoID["XC_GGA_X_B88"] = 106;
    _stringtoID["XC_GGA_X_G96"] = 107;
    _stringtoID["XC_GGA_X_PW86"] = 108;
    _stringtoID["XC_GGA_X_PW91"] = 109;
    _stringtoID["XC_GGA_X_OPTX"] = 110;
    _stringtoID["XC_GGA_X_DK87_R1"] = 111;
    _stringtoID["XC_GGA_X_DK87_R2"] = 112;
    _stringtoID["XC_GGA_X_LG93"] = 113;
    _stringtoID["XC_GGA_X_FT97_A"] = 114;
    _stringtoID["XC_GGA_X_FT97_B"] = 115;
    _stringtoID["XC_GGA_X_PBE_SOL"] = 116;
    _stringtoID["XC_GGA_X_RPBE"] = 117;
    _stringtoID["XC_GGA_X_WC"] = 118;
    _stringtoID["XC_GGA_X_MPW91"] = 119;
    _stringtoID["XC_GGA_X_AM05"] = 120;
    _stringtoID["XC_GGA_X_PBEA"] = 121;
    _stringtoID["XC_GGA_X_MPBE"] = 122;
    _stringtoID["XC_GGA_X_XPBE"] = 123;
    _stringtoID["XC_GGA_X_2D_B86_MGC"] = 124;
    _stringtoID["XC_GGA_X_BAYESIAN"] = 125;
    _stringtoID["XC_GGA_X_PBE_JSJR"] = 126;
    _stringtoID["XC_GGA_X_2D_B88"] = 127;
    _stringtoID["XC_GGA_X_2D_B86"] = 128;
    _stringtoID["XC_GGA_X_2D_PBE"] = 129;
    _stringtoID["XC_GGA_C_PBE"] = 130;
    _stringtoID["XC_GGA_C_LYP"] = 131;
    _stringtoID["XC_GGA_C_P86"] = 132;
    _stringtoID["XC_GGA_C_PBE_SOL"] = 133;
    _stringtoID["XC_GGA_C_PW91"] = 134;
    _stringtoID["XC_GGA_C_AM05"] = 135;
    _stringtoID["XC_GGA_C_XPBE"] = 136;
    _stringtoID["XC_GGA_C_LM"] = 137;
    _stringtoID["XC_GGA_C_PBE_JRGX"] = 138;
    _stringtoID["XC_GGA_X_OPTB88_VDW"] = 139;
    _stringtoID["XC_GGA_X_PBEK1_VDW"] = 140;
    _stringtoID["XC_GGA_X_OPTPBE_VDW"] = 141;
    _stringtoID["XC_GGA_X_RGE2"] = 142;
    _stringtoID["XC_GGA_C_RGE2"] = 143;
    _stringtoID["XC_GGA_X_RPW86"] = 144;
    _stringtoID["XC_GGA_X_KT1"] = 145;
    _stringtoID["XC_GGA_XC_KT2"] = 146;
    _stringtoID["XC_GGA_C_WL"] = 147;
    _stringtoID["XC_GGA_C_WI"] = 148;
    _stringtoID["XC_GGA_X_MB88"] = 149;
    _stringtoID["XC_GGA_X_SOGGA"] = 150;
    _stringtoID["XC_GGA_X_SOGGA11"] = 151;
    _stringtoID["XC_GGA_C_SOGGA11"] = 152;
    _stringtoID["XC_GGA_C_WI0"] = 153;
    _stringtoID["XC_GGA_XC_TH1"] = 154;
    _stringtoID["XC_GGA_XC_TH2"] = 155;
    _stringtoID["XC_GGA_XC_TH3"] = 156;
    _stringtoID["XC_GGA_XC_TH4"] = 157;
    _stringtoID["XC_GGA_X_C09X"] = 158;
    _stringtoID["XC_GGA_C_SOGGA11_X"] = 159;
    _stringtoID["XC_GGA_X_LB"] = 160;
    _stringtoID["XC_GGA_XC_HCTH_93"] = 161;
    _stringtoID["XC_GGA_XC_HCTH_120"] = 162;
    _stringtoID["XC_GGA_XC_HCTH_147"] = 163;
    _stringtoID["XC_GGA_XC_HCTH_407"] = 164;
    _stringtoID["XC_GGA_XC_EDF1"] = 165;
    _stringtoID["XC_GGA_XC_XLYP"] = 166;
    _stringtoID["XC_GGA_XC_B97"] = 167;
    _stringtoID["XC_GGA_XC_B97_1"] = 168;
    _stringtoID["XC_GGA_XC_B97_2"] = 169;
    _stringtoID["XC_GGA_XC_B97_D"] = 170;
    _stringtoID["XC_GGA_XC_B97_K"] = 171;
    _stringtoID["XC_GGA_XC_B97_3"] = 172;
    _stringtoID["XC_GGA_XC_PBE1W"] = 173;
    _stringtoID["XC_GGA_XC_MPWLYP1W"] = 174;
    _stringtoID["XC_GGA_XC_PBELYP1W"] = 175;
    _stringtoID["XC_GGA_XC_SB98_1a"] = 176;
    _stringtoID["XC_GGA_XC_SB98_1b"] = 177;
    _stringtoID["XC_GGA_XC_SB98_1c"] = 178;
    _stringtoID["XC_GGA_XC_SB98_2a"] = 179;
    _stringtoID["XC_GGA_XC_SB98_2b"] = 180;
    _stringtoID["XC_GGA_XC_SB98_2c"] = 181;
    _stringtoID["XC_GGA_X_LBM"] = 182;
    _stringtoID["XC_GGA_X_OL2"] = 183;
    _stringtoID["XC_GGA_X_APBE"] = 184;
    _stringtoID["XC_GGA_K_APBE"] = 185;
    _stringtoID["XC_GGA_C_APBE"] = 186;
    _stringtoID["XC_GGA_K_TW1"] = 187;
    _stringtoID["XC_GGA_K_TW2"] = 188;
    _stringtoID["XC_GGA_K_TW3"] = 189;
    _stringtoID["XC_GGA_K_TW4"] = 190;
    _stringtoID["XC_GGA_X_HTBS"] = 191;
    _stringtoID["XC_GGA_X_AIRY"] = 192;
    _stringtoID["XC_GGA_X_LAG"] = 193;
    _stringtoID["XC_GGA_XC_MOHLYP"] = 194;
    _stringtoID["XC_GGA_XC_MOHLYP2"] = 195;
    _stringtoID["XC_GGA_XC_TH_FL"] = 196;
    _stringtoID["XC_GGA_XC_TH_FC"] = 197;
    _stringtoID["XC_GGA_XC_TH_FCFO"] = 198;
    _stringtoID["XC_GGA_XC_TH_FCO"] = 199;
    _stringtoID["XC_GGA_C_OPTC"] = 200;
    _stringtoID["XC_GGA_K_VW"] = 500;
    _stringtoID["XC_GGA_K_GE2"] = 501;
    _stringtoID["XC_GGA_K_GOLDEN"] = 502;
    _stringtoID["XC_GGA_K_YT65"] = 503;
    _stringtoID["XC_GGA_K_BALTIN"] = 504;
    _stringtoID["XC_GGA_K_LIEB"] = 505;
    _stringtoID["XC_GGA_K_ABSP1"] = 506;
    _stringtoID["XC_GGA_K_ABSP2"] = 507;
    _stringtoID["XC_GGA_K_GR"] = 508;
    _stringtoID["XC_GGA_K_LUDENA"] = 509;
    _stringtoID["XC_GGA_K_GP85"] = 510;
    _stringtoID["XC_GGA_K_PEARSON"] = 511;
    _stringtoID["XC_GGA_K_OL1"] = 512;
    _stringtoID["XC_GGA_K_OL2"] = 513;
    _stringtoID["XC_GGA_K_FR_B88"] = 514;
    _stringtoID["XC_GGA_K_FR_PW86"] = 515;
    _stringtoID["XC_GGA_K_DK"] = 516;
    _stringtoID["XC_GGA_K_PERDEW"] = 517;
    _stringtoID["XC_GGA_K_VSK"] = 518;
    _stringtoID["XC_GGA_K_VJKS"] = 519;
    _stringtoID["XC_GGA_K_ERNZERHOF"] = 520;
    _stringtoID["XC_GGA_K_LC94"] = 521;
    _stringtoID["XC_GGA_K_LLP"] = 522;
    _stringtoID["XC_GGA_K_THAKKAR"] = 523;
    _stringtoID["XC_GGA_X_WPBEH"] = 524;
    _stringtoID["XC_GGA_X_HJS_PBE"] = 525;
    _stringtoID["XC_GGA_X_HJS_PBE_SOL"] = 526;
    _stringtoID["XC_GGA_X_HJS_B88"] = 527;
    _stringtoID["XC_GGA_X_HJS_B97X"] = 528;
    _stringtoID["XC_GGA_X_ITYH"] = 529;
    _stringtoID["XC_GGA_X_SFAT"] = 530;
    _stringtoID["XC_HYB_GGA_X_N12_SX"] = 81;
    _stringtoID["XC_HYB_GGA_XC_B3PW91"] = 401;
    _stringtoID["XC_HYB_GGA_XC_B3LYP"] = 402;
    _stringtoID["XC_HYB_GGA_XC_B3P86"] = 403;
    _stringtoID["XC_HYB_GGA_XC_O3LYP"] = 404;
    _stringtoID["XC_HYB_GGA_XC_mPW1K"] = 405;
    _stringtoID["XC_HYB_GGA_XC_PBEH"] = 406;
    _stringtoID["XC_HYB_GGA_XC_B97"] = 407;
    _stringtoID["XC_HYB_GGA_XC_B97_1"] = 408;
    _stringtoID["XC_HYB_GGA_XC_B97_2"] = 410;
    _stringtoID["XC_HYB_GGA_XC_X3LYP"] = 411;
    _stringtoID["XC_HYB_GGA_XC_B1WC"] = 412;
    _stringtoID["XC_HYB_GGA_XC_B97_K"] = 413;
    _stringtoID["XC_HYB_GGA_XC_B97_3"] = 414;
    _stringtoID["XC_HYB_GGA_XC_MPW3PW"] = 415;
    _stringtoID["XC_HYB_GGA_XC_B1LYP"] = 416;
    _stringtoID["XC_HYB_GGA_XC_B1PW91"] = 417;
    _stringtoID["XC_HYB_GGA_XC_mPW1PW"] = 418;
    _stringtoID["XC_HYB_GGA_XC_MPW3LYP"] = 419;
    _stringtoID["XC_HYB_GGA_XC_SB98_1a"] = 420;
    _stringtoID["XC_HYB_GGA_XC_SB98_1b"] = 421;
    _stringtoID["XC_HYB_GGA_XC_SB98_1c"] = 422;
    _stringtoID["XC_HYB_GGA_XC_SB98_2a"] = 423;
    _stringtoID["XC_HYB_GGA_XC_SB98_2b"] = 424;
    _stringtoID["XC_HYB_GGA_XC_SB98_2c"] = 425;
    _stringtoID["XC_HYB_GGA_X_SOGGA11_X"] = 426;
    _stringtoID["XC_HYB_GGA_XC_HSE03"] = 427;
    _stringtoID["XC_HYB_GGA_XC_HSE06"] = 428;
    _stringtoID["XC_HYB_GGA_XC_HJS_PBE"] = 429;
    _stringtoID["XC_HYB_GGA_XC_HJS_PBE_SOL"] = 430;
    _stringtoID["XC_HYB_GGA_XC_HJS_B88"] = 431;
    _stringtoID["XC_HYB_GGA_XC_HJS_B97X"] = 432;
    _stringtoID["XC_HYB_GGA_XC_CAM_B3LYP"] = 433;
    _stringtoID["XC_HYB_GGA_XC_TUNED_CAM_B3LYP"] = 434;
    _stringtoID["XC_HYB_GGA_XC_BHANDH"] = 435;
    _stringtoID["XC_HYB_GGA_XC_BHANDHLYP"] = 436;
    _stringtoID["XC_HYB_GGA_XC_MB3LYP_RC04"] = 437;
    _stringtoID["XC_HYB_GGA_XC_MPWLYP1M"] = 453;
    _stringtoID["XC_HYB_GGA_XC_REVB3LYP"] = 454;
    _stringtoID["XC_HYB_GGA_XC_CAMY_BLYP"] = 455;
    _stringtoID["XC_HYB_GGA_XC_PBE0_13"] = 456;
    _stringtoID["XC_MGGA_XC_OTPSS_D"] = 64;
    _stringtoID["XC_MGGA_C_CS"] = 72;
    _stringtoID["XC_MGGA_C_MN12_SX"] = 73;
    _stringtoID["XC_MGGA_C_MN12_L"] = 74;
    _stringtoID["XC_MGGA_C_M11_L"] = 75;
    _stringtoID["XC_MGGA_C_M11"] = 76;
    _stringtoID["XC_MGGA_C_M08_SO"] = 77;
    _stringtoID["XC_MGGA_C_M08_HX"] = 78;
    _stringtoID["XC_MGGA_X_LTA"] = 201;
    _stringtoID["XC_MGGA_X_TPSS"] = 202;
    _stringtoID["XC_MGGA_X_M06_L"] = 203;
    _stringtoID["XC_MGGA_X_GVT4"] = 204;
    _stringtoID["XC_MGGA_X_TAU_HCTH"] = 205;
    _stringtoID["XC_MGGA_X_BR89"] = 206;
    _stringtoID["XC_MGGA_X_BJ06"] = 207;
    _stringtoID["XC_MGGA_X_TB09"] = 208;
    _stringtoID["XC_MGGA_X_RPP09"] = 209;
    _stringtoID["XC_MGGA_X_2D_PRHG07"] = 210;
    _stringtoID["XC_MGGA_X_2D_PRHG07_PRP10"] = 211;
    _stringtoID["XC_MGGA_X_REVTPSS"] = 212;
    _stringtoID["XC_MGGA_X_PKZB"] = 213;
    _stringtoID["XC_MGGA_X_M05"] = 214;
    _stringtoID["XC_MGGA_X_M05_2X"] = 215;
    _stringtoID["XC_MGGA_X_M06_HF"] = 216;
    _stringtoID["XC_MGGA_X_M06"] = 217;
    _stringtoID["XC_MGGA_X_M06_2X"] = 218;
    _stringtoID["XC_MGGA_X_M08_HX"] = 219;
    _stringtoID["XC_MGGA_X_M08_SO"] = 220;
    _stringtoID["XC_MGGA_X_MS0"] = 221;
    _stringtoID["XC_MGGA_X_MS1"] = 222;
    _stringtoID["XC_MGGA_X_MS2"] = 223;
    _stringtoID["XC_MGGA_X_MS2H"] = 224;
    _stringtoID["XC_MGGA_X_M11_L"] = 226;
    _stringtoID["XC_MGGA_X_MN12_L"] = 227;
    _stringtoID["XC_MGGA_X_MN12_SX"] = 228;
    _stringtoID["XC_MGGA_C_CC06"] = 229;
    _stringtoID["XC_MGGA_X_MK00"] = 230;
    _stringtoID["XC_MGGA_C_TPSS"] = 231;
    _stringtoID["XC_MGGA_C_VSXC"] = 232;
    _stringtoID["XC_MGGA_C_M06_L"] = 233;
    _stringtoID["XC_MGGA_C_M06_HF"] = 234;
    _stringtoID["XC_MGGA_C_M06"] = 235;
    _stringtoID["XC_MGGA_C_M06_2X"] = 236;
    _stringtoID["XC_MGGA_C_M05"] = 237;
    _stringtoID["XC_MGGA_C_M05_2X"] = 238;
    _stringtoID["XC_MGGA_C_PKZB"] = 239;
    _stringtoID["XC_MGGA_C_BC95"] = 240;
    _stringtoID["XC_MGGA_C_REVTPSS"] = 241;
    _stringtoID["XC_MGGA_XC_TPSSLYP1W"] = 242;
    _stringtoID["XC_MGGA_X_MK00B"] = 243;
    _stringtoID["XC_MGGA_X_BLOC"] = 244;
    _stringtoID["XC_MGGA_X_MODTPSS"] = 245;
    _stringtoID["XC_HYB_MGGA_X_M11"] = 225;
    _stringtoID["XC_HYB_MGGA_XC_M05"] = 438;
    _stringtoID["XC_HYB_MGGA_XC_M05_2X"] = 439;
    _stringtoID["XC_HYB_MGGA_XC_B88B95"] = 440;
    _stringtoID["XC_HYB_MGGA_XC_B86B95"] = 441;
    _stringtoID["XC_HYB_MGGA_XC_PW86B95"] = 442;
    _stringtoID["XC_HYB_MGGA_XC_BB1K"] = 443;
    _stringtoID["XC_HYB_MGGA_XC_M06_HF"] = 444;
    _stringtoID["XC_HYB_MGGA_XC_MPW1B95"] = 445;
    _stringtoID["XC_HYB_MGGA_XC_MPWB1K"] = 446;
    _stringtoID["XC_HYB_MGGA_XC_X1B95"] = 447;
    _stringtoID["XC_HYB_MGGA_XC_XB1K"] = 448;
    _stringtoID["XC_HYB_MGGA_XC_M06"] = 449;
    _stringtoID["XC_HYB_MGGA_XC_M06_2X"] = 450;
    _stringtoID["XC_HYB_MGGA_XC_PW6B95"] = 451;
    _stringtoID["XC_HYB_MGGA_XC_PWB6K"] = 452;
    _stringtoID["XC_HYB_MGGA_XC_TPSSH"] = 457;
    _stringtoID["XC_HYB_MGGA_XC_REVTPSSH"] = 458;
  }
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_VXCFUNCTIONALS_H
