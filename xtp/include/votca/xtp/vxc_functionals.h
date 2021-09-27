/*
 *            Copyright 2009-2020 The VOTCA Development Team
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
#ifndef VOTCA_XTP_VXC_FUNCTIONALS_H
#define VOTCA_XTP_VXC_FUNCTIONALS_H

// Standard includes
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

  int getID(std::string name) const {
    if (!stringtoID_.count(name)) {
      throw std::runtime_error("Functional " + name + " not supported");
    }
    return stringtoID_.at(name);
  }

 private:
  std::map<std::string, int> stringtoID_;

  inline void FillMaps() { FillstringtoID(); }

  inline void FillstringtoID() {

    stringtoID_["XC_LDA_X"] = 1;
    stringtoID_["XC_LDA_C_WIGNER"] = 2;
    stringtoID_["XC_LDA_C_RPA"] = 3;
    stringtoID_["XC_LDA_C_HL"] = 4;
    stringtoID_["XC_LDA_C_GL"] = 5;
    stringtoID_["XC_LDA_C_XALPHA"] = 6;
    stringtoID_["XC_LDA_C_VWN"] = 7;
    stringtoID_["XC_LDA_C_VWN_RPA"] = 8;
    stringtoID_["XC_LDA_C_PZ"] = 9;
    stringtoID_["XC_LDA_C_PZ_MOD"] = 10;
    stringtoID_["XC_LDA_C_OB_PZ"] = 11;
    stringtoID_["XC_LDA_C_PW"] = 12;
    stringtoID_["XC_LDA_C_PW_MOD"] = 13;
    stringtoID_["XC_LDA_C_OB_PW"] = 14;
    stringtoID_["XC_LDA_C_2D_AMGB"] = 15;
    stringtoID_["XC_LDA_C_2D_PRM"] = 16;
    stringtoID_["XC_LDA_C_vBH"] = 17;
    stringtoID_["XC_LDA_C_1D_CSC"] = 18;
    stringtoID_["XC_LDA_X_2D"] = 19;
    stringtoID_["XC_LDA_XC_TETER93"] = 20;
    stringtoID_["XC_LDA_X_1D"] = 21;
    stringtoID_["XC_LDA_C_ML1"] = 22;
    stringtoID_["XC_LDA_C_ML2"] = 23;
    stringtoID_["XC_LDA_C_GOMBAS"] = 24;
    stringtoID_["XC_LDA_C_PW_RPA"] = 25;
    stringtoID_["XC_LDA_C_1D_LOOS"] = 26;
    stringtoID_["XC_LDA_C_RC04"] = 27;
    stringtoID_["XC_LDA_C_VWN_1"] = 28;
    stringtoID_["XC_LDA_C_VWN_2"] = 29;
    stringtoID_["XC_LDA_C_VWN_3"] = 30;
    stringtoID_["XC_LDA_C_VWN_4"] = 31;
    stringtoID_["XC_LDA_K_TF"] = 50;
    stringtoID_["XC_LDA_K_LP"] = 51;
    stringtoID_["XC_GGA_C_Q2D"] = 47;
    stringtoID_["XC_GGA_X_Q2D"] = 48;
    stringtoID_["XC_GGA_X_PBE_MOL"] = 49;
    stringtoID_["XC_GGA_K_TFVW"] = 52;
    stringtoID_["XC_GGA_K_REVAPBEINT"] = 53;
    stringtoID_["XC_GGA_K_APBEINT"] = 54;
    stringtoID_["XC_GGA_K_REVAPBE"] = 55;
    stringtoID_["XC_GGA_X_AK13"] = 56;
    stringtoID_["XC_GGA_K_MEYER"] = 57;
    stringtoID_["XC_GGA_X_LV_RPW86"] = 58;
    stringtoID_["XC_GGA_X_PBE_TCA"] = 59;
    stringtoID_["XC_GGA_X_PBEINT"] = 60;
    stringtoID_["XC_GGA_C_ZPBEINT"] = 61;
    stringtoID_["XC_GGA_C_PBEINT"] = 62;
    stringtoID_["XC_GGA_C_ZPBESOL"] = 63;
    stringtoID_["XC_GGA_XC_OPBE_D"] = 65;
    stringtoID_["XC_GGA_XC_OPWLYP_D"] = 66;
    stringtoID_["XC_GGA_XC_OBLYP_D"] = 67;
    stringtoID_["XC_GGA_X_VMT84_GE"] = 68;
    stringtoID_["XC_GGA_X_VMT84_PBE"] = 69;
    stringtoID_["XC_GGA_X_VMT_GE"] = 70;
    stringtoID_["XC_GGA_X_VMT_PBE"] = 71;
    stringtoID_["XC_GGA_C_N12_SX"] = 79;
    stringtoID_["XC_GGA_C_N12"] = 80;
    stringtoID_["XC_GGA_X_N12"] = 82;
    stringtoID_["XC_GGA_C_VPBE"] = 83;
    stringtoID_["XC_GGA_C_OP_XALPHA"] = 84;
    stringtoID_["XC_GGA_C_OP_G96"] = 85;
    stringtoID_["XC_GGA_C_OP_PBE"] = 86;
    stringtoID_["XC_GGA_C_OP_B88"] = 87;
    stringtoID_["XC_GGA_C_FT97"] = 88;
    stringtoID_["XC_GGA_C_SPBE"] = 89;
    stringtoID_["XC_GGA_X_SSB_SW"] = 90;
    stringtoID_["XC_GGA_X_SSB"] = 91;
    stringtoID_["XC_GGA_X_SSB_D"] = 92;
    stringtoID_["XC_GGA_XC_HCTH_407P"] = 93;
    stringtoID_["XC_GGA_XC_HCTH_P76"] = 94;
    stringtoID_["XC_GGA_XC_HCTH_P14"] = 95;
    stringtoID_["XC_GGA_XC_B97_GGA1"] = 96;
    stringtoID_["XC_GGA_XC_HCTH_A"] = 97;
    stringtoID_["XC_GGA_X_BPCCAC"] = 98;
    stringtoID_["XC_GGA_C_REVTCA"] = 99;
    stringtoID_["XC_GGA_C_TCA"] = 100;
    stringtoID_["XC_GGA_X_PBE"] = 101;
    stringtoID_["XC_GGA_X_PBE_R"] = 102;
    stringtoID_["XC_GGA_X_B86"] = 103;
    stringtoID_["XC_GGA_X_HERMAN"] = 104;
    stringtoID_["XC_GGA_X_B86_MGC"] = 105;
    stringtoID_["XC_GGA_X_B88"] = 106;
    stringtoID_["XC_GGA_X_G96"] = 107;
    stringtoID_["XC_GGA_X_PW86"] = 108;
    stringtoID_["XC_GGA_X_PW91"] = 109;
    stringtoID_["XC_GGA_X_OPTX"] = 110;
    stringtoID_["XC_GGA_X_DK87_R1"] = 111;
    stringtoID_["XC_GGA_X_DK87_R2"] = 112;
    stringtoID_["XC_GGA_X_LG93"] = 113;
    stringtoID_["XC_GGA_X_FT97_A"] = 114;
    stringtoID_["XC_GGA_X_FT97_B"] = 115;
    stringtoID_["XC_GGA_X_PBE_SOL"] = 116;
    stringtoID_["XC_GGA_X_RPBE"] = 117;
    stringtoID_["XC_GGA_X_WC"] = 118;
    stringtoID_["XC_GGA_X_MPW91"] = 119;
    stringtoID_["XC_GGA_X_AM05"] = 120;
    stringtoID_["XC_GGA_X_PBEA"] = 121;
    stringtoID_["XC_GGA_X_MPBE"] = 122;
    stringtoID_["XC_GGA_X_XPBE"] = 123;
    stringtoID_["XC_GGA_X_2D_B86_MGC"] = 124;
    stringtoID_["XC_GGA_X_BAYESIAN"] = 125;
    stringtoID_["XC_GGA_X_PBE_JSJR"] = 126;
    stringtoID_["XC_GGA_X_2D_B88"] = 127;
    stringtoID_["XC_GGA_X_2D_B86"] = 128;
    stringtoID_["XC_GGA_X_2D_PBE"] = 129;
    stringtoID_["XC_GGA_C_PBE"] = 130;
    stringtoID_["XC_GGA_C_LYP"] = 131;
    stringtoID_["XC_GGA_C_P86"] = 132;
    stringtoID_["XC_GGA_C_PBE_SOL"] = 133;
    stringtoID_["XC_GGA_C_PW91"] = 134;
    stringtoID_["XC_GGA_C_AM05"] = 135;
    stringtoID_["XC_GGA_C_XPBE"] = 136;
    stringtoID_["XC_GGA_C_LM"] = 137;
    stringtoID_["XC_GGA_C_PBE_JRGX"] = 138;
    stringtoID_["XC_GGA_X_OPTB88_VDW"] = 139;
    stringtoID_["XC_GGA_X_PBEK1_VDW"] = 140;
    stringtoID_["XC_GGA_X_OPTPBE_VDW"] = 141;
    stringtoID_["XC_GGA_X_RGE2"] = 142;
    stringtoID_["XC_GGA_C_RGE2"] = 143;
    stringtoID_["XC_GGA_X_RPW86"] = 144;
    stringtoID_["XC_GGA_X_KT1"] = 145;
    stringtoID_["XC_GGA_XC_KT2"] = 146;
    stringtoID_["XC_GGA_C_WL"] = 147;
    stringtoID_["XC_GGA_C_WI"] = 148;
    stringtoID_["XC_GGA_X_MB88"] = 149;
    stringtoID_["XC_GGA_X_SOGGA"] = 150;
    stringtoID_["XC_GGA_X_SOGGA11"] = 151;
    stringtoID_["XC_GGA_C_SOGGA11"] = 152;
    stringtoID_["XC_GGA_C_WI0"] = 153;
    stringtoID_["XC_GGA_XC_TH1"] = 154;
    stringtoID_["XC_GGA_XC_TH2"] = 155;
    stringtoID_["XC_GGA_XC_TH3"] = 156;
    stringtoID_["XC_GGA_XC_TH4"] = 157;
    stringtoID_["XC_GGA_X_C09X"] = 158;
    stringtoID_["XC_GGA_C_SOGGA11_X"] = 159;
    stringtoID_["XC_GGA_X_LB"] = 160;
    stringtoID_["XC_GGA_XC_HCTH_93"] = 161;
    stringtoID_["XC_GGA_XC_HCTH_120"] = 162;
    stringtoID_["XC_GGA_XC_HCTH_147"] = 163;
    stringtoID_["XC_GGA_XC_HCTH_407"] = 164;
    stringtoID_["XC_GGA_XC_EDF1"] = 165;
    stringtoID_["XC_GGA_XC_XLYP"] = 166;
    stringtoID_["XC_GGA_XC_B97"] = 167;
    stringtoID_["XC_GGA_XC_B97_1"] = 168;
    stringtoID_["XC_GGA_XC_B97_2"] = 169;
    stringtoID_["XC_GGA_XC_B97_D"] = 170;
    stringtoID_["XC_GGA_XC_B97_K"] = 171;
    stringtoID_["XC_GGA_XC_B97_3"] = 172;
    stringtoID_["XC_GGA_XC_PBE1W"] = 173;
    stringtoID_["XC_GGA_XC_MPWLYP1W"] = 174;
    stringtoID_["XC_GGA_XC_PBELYP1W"] = 175;
    stringtoID_["XC_GGA_XC_SB98_1a"] = 176;
    stringtoID_["XC_GGA_XC_SB98_1b"] = 177;
    stringtoID_["XC_GGA_XC_SB98_1c"] = 178;
    stringtoID_["XC_GGA_XC_SB98_2a"] = 179;
    stringtoID_["XC_GGA_XC_SB98_2b"] = 180;
    stringtoID_["XC_GGA_XC_SB98_2c"] = 181;
    stringtoID_["XC_GGA_X_LBM"] = 182;
    stringtoID_["XC_GGA_X_OL2"] = 183;
    stringtoID_["XC_GGA_X_APBE"] = 184;
    stringtoID_["XC_GGA_K_APBE"] = 185;
    stringtoID_["XC_GGA_C_APBE"] = 186;
    stringtoID_["XC_GGA_K_TW1"] = 187;
    stringtoID_["XC_GGA_K_TW2"] = 188;
    stringtoID_["XC_GGA_K_TW3"] = 189;
    stringtoID_["XC_GGA_K_TW4"] = 190;
    stringtoID_["XC_GGA_X_HTBS"] = 191;
    stringtoID_["XC_GGA_X_AIRY"] = 192;
    stringtoID_["XC_GGA_X_LAG"] = 193;
    stringtoID_["XC_GGA_XC_MOHLYP"] = 194;
    stringtoID_["XC_GGA_XC_MOHLYP2"] = 195;
    stringtoID_["XC_GGA_XC_TH_FL"] = 196;
    stringtoID_["XC_GGA_XC_TH_FC"] = 197;
    stringtoID_["XC_GGA_XC_TH_FCFO"] = 198;
    stringtoID_["XC_GGA_XC_TH_FCO"] = 199;
    stringtoID_["XC_GGA_C_OPTC"] = 200;
    stringtoID_["XC_GGA_K_VW"] = 500;
    stringtoID_["XC_GGA_K_GE2"] = 501;
    stringtoID_["XC_GGA_K_GOLDEN"] = 502;
    stringtoID_["XC_GGA_K_YT65"] = 503;
    stringtoID_["XC_GGA_K_BALTIN"] = 504;
    stringtoID_["XC_GGA_K_LIEB"] = 505;
    stringtoID_["XC_GGA_K_ABSP1"] = 506;
    stringtoID_["XC_GGA_K_ABSP2"] = 507;
    stringtoID_["XC_GGA_K_GR"] = 508;
    stringtoID_["XC_GGA_K_LUDENA"] = 509;
    stringtoID_["XC_GGA_K_GP85"] = 510;
    stringtoID_["XC_GGA_K_PEARSON"] = 511;
    stringtoID_["XC_GGA_K_OL1"] = 512;
    stringtoID_["XC_GGA_K_OL2"] = 513;
    stringtoID_["XC_GGA_K_FR_B88"] = 514;
    stringtoID_["XC_GGA_K_FR_PW86"] = 515;
    stringtoID_["XC_GGA_K_DK"] = 516;
    stringtoID_["XC_GGA_K_PERDEW"] = 517;
    stringtoID_["XC_GGA_K_VSK"] = 518;
    stringtoID_["XC_GGA_K_VJKS"] = 519;
    stringtoID_["XC_GGA_K_ERNZERHOF"] = 520;
    stringtoID_["XC_GGA_K_LC94"] = 521;
    stringtoID_["XC_GGA_K_LLP"] = 522;
    stringtoID_["XC_GGA_K_THAKKAR"] = 523;
    stringtoID_["XC_GGA_X_WPBEH"] = 524;
    stringtoID_["XC_GGA_X_HJS_PBE"] = 525;
    stringtoID_["XC_GGA_X_HJS_PBE_SOL"] = 526;
    stringtoID_["XC_GGA_X_HJS_B88"] = 527;
    stringtoID_["XC_GGA_X_HJS_B97X"] = 528;
    stringtoID_["XC_GGA_X_ITYH"] = 529;
    stringtoID_["XC_GGA_X_SFAT"] = 530;
    stringtoID_["XC_HYB_GGA_X_N12_SX"] = 81;
    stringtoID_["XC_HYB_GGA_XC_B3PW91"] = 401;
    stringtoID_["XC_HYB_GGA_XC_B3LYP"] = 402;
    stringtoID_["XC_HYB_GGA_XC_B3P86"] = 403;
    stringtoID_["XC_HYB_GGA_XC_O3LYP"] = 404;
    stringtoID_["XC_HYB_GGA_XC_mPW1K"] = 405;
    stringtoID_["XC_HYB_GGA_XC_PBEH"] = 406;
    stringtoID_["XC_HYB_GGA_XC_B97"] = 407;
    stringtoID_["XC_HYB_GGA_XC_B97_1"] = 408;
    stringtoID_["XC_HYB_GGA_XC_B97_2"] = 410;
    stringtoID_["XC_HYB_GGA_XC_X3LYP"] = 411;
    stringtoID_["XC_HYB_GGA_XC_B1WC"] = 412;
    stringtoID_["XC_HYB_GGA_XC_B97_K"] = 413;
    stringtoID_["XC_HYB_GGA_XC_B97_3"] = 414;
    stringtoID_["XC_HYB_GGA_XC_MPW3PW"] = 415;
    stringtoID_["XC_HYB_GGA_XC_B1LYP"] = 416;
    stringtoID_["XC_HYB_GGA_XC_B1PW91"] = 417;
    stringtoID_["XC_HYB_GGA_XC_mPW1PW"] = 418;
    stringtoID_["XC_HYB_GGA_XC_MPW3LYP"] = 419;
    stringtoID_["XC_HYB_GGA_XC_SB98_1a"] = 420;
    stringtoID_["XC_HYB_GGA_XC_SB98_1b"] = 421;
    stringtoID_["XC_HYB_GGA_XC_SB98_1c"] = 422;
    stringtoID_["XC_HYB_GGA_XC_SB98_2a"] = 423;
    stringtoID_["XC_HYB_GGA_XC_SB98_2b"] = 424;
    stringtoID_["XC_HYB_GGA_XC_SB98_2c"] = 425;
    stringtoID_["XC_HYB_GGA_X_SOGGA11_X"] = 426;
    stringtoID_["XC_HYB_GGA_XC_HSE03"] = 427;
    stringtoID_["XC_HYB_GGA_XC_HSE06"] = 428;
    stringtoID_["XC_HYB_GGA_XC_HJS_PBE"] = 429;
    stringtoID_["XC_HYB_GGA_XC_HJS_PBE_SOL"] = 430;
    stringtoID_["XC_HYB_GGA_XC_HJS_B88"] = 431;
    stringtoID_["XC_HYB_GGA_XC_HJS_B97X"] = 432;
    stringtoID_["XC_HYB_GGA_XC_CAM_B3LYP"] = 433;
    stringtoID_["XC_HYB_GGA_XC_TUNED_CAM_B3LYP"] = 434;
    stringtoID_["XC_HYB_GGA_XC_BHANDH"] = 435;
    stringtoID_["XC_HYB_GGA_XC_BHANDHLYP"] = 436;
    stringtoID_["XC_HYB_GGA_XC_MB3LYP_RC04"] = 437;
    stringtoID_["XC_HYB_GGA_XC_MPWLYP1M"] = 453;
    stringtoID_["XC_HYB_GGA_XC_REVB3LYP"] = 454;
    stringtoID_["XC_HYB_GGA_XC_CAMY_BLYP"] = 455;
    stringtoID_["XC_HYB_GGA_XC_PBE0_13"] = 456;
  }
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_VXC_FUNCTIONALS_H
