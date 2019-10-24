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
#ifndef VOTCA_XTP_AOMATRIX_H
#define VOTCA_XTP_AOMATRIX_H

#include <votca/xtp/aobasis.h>
#include <votca/xtp/aoshell.h>
#include <votca/xtp/multiarray.h>
#include <votca/xtp/staticsite.h>

namespace votca {
namespace xtp {

namespace Cart {
enum Index {
  s,  // s
  x,
  y,
  z,  // p
  xx,
  xy,
  xz,
  yy,
  yz,
  zz,  // d
  xxx,
  xxy,
  xxz,
  xyy,
  xyz,
  xzz,
  yyy,
  yyz,
  yzz,
  zzz,  // f
  xxxx,
  xxxy,
  xxxz,
  xxyy,
  xxyz,
  xxzz,
  xyyy,
  xyyz,
  xyzz,
  xzzz,
  yyyy,
  yyyz,
  yyzz,
  yzzz,
  zzzz,  // g
  xxxxx,
  xxxxy,
  xxxxz,
  xxxyy,
  xxxyz,
  xxxzz,
  xxyyy,
  xxyyz,
  xxyzz,
  xxzzz,
  xyyyy,
  xyyyz,
  xyyzz,
  xyzzz,
  xzzzz,
  yyyyy,
  yyyyz,
  yyyzz,
  yyzzz,
  yzzzz,
  zzzzz,  // h

  xxxxxx,
  xxxxxy,
  xxxxxz,
  xxxxyy,
  xxxxyz,
  xxxxzz,
  xxxyyy,
  xxxyyz,
  xxxyzz,
  xxxzzz,
  xxyyyy,
  xxyyyz,
  xxyyzz,
  xxyzzz,  // i
  xxzzzz,
  xyyyyy,
  xyyyyz,
  xyyyzz,
  xyyzzz,
  xyzzzz,
  xzzzzz,
  yyyyyy,
  yyyyyz,
  yyyyzz,
  yyyzzz,
  yyzzzz,
  yzzzzz,
  zzzzzz,

  xxxxxxx,
  xxxxxxy,
  xxxxxxz,
  xxxxxyy,
  xxxxxyz,
  xxxxxzz,
  xxxxyyy,
  xxxxyyz,
  xxxxyzz,
  xxxxzzz,
  xxxyyyy,
  xxxyyyz,  // j
  xxxyyzz,
  xxxyzzz,
  xxxzzzz,
  xxyyyyy,
  xxyyyyz,
  xxyyyzz,
  xxyyzzz,
  xxyzzzz,
  xxzzzzz,
  xyyyyyy,
  xyyyyyz,
  xyyyyzz,
  xyyyzzz,
  xyyzzzz,
  xyzzzzz,
  xzzzzzz,
  yyyyyyy,
  yyyyyyz,
  yyyyyzz,
  yyyyzzz,
  yyyzzzz,
  yyzzzzz,
  yzzzzzz,
  zzzzzzz,

  xxxxxxxx,
  xxxxxxxy,
  xxxxxxxz,
  xxxxxxyy,
  xxxxxxyz,
  xxxxxxzz,
  xxxxxyyy,
  xxxxxyyz,
  xxxxxyzz,
  xxxxxzzz,
  xxxxyyyy,
  xxxxyyyz,
  xxxxyyzz,
  xxxxyzzz,
  xxxxzzzz,  // k
  xxxyyyyy,
  xxxyyyyz,
  xxxyyyzz,
  xxxyyzzz,
  xxxyzzzz,
  xxxzzzzz,
  xxyyyyyy,
  xxyyyyyz,
  xxyyyyzz,
  xxyyyzzz,
  xxyyzzzz,
  xxyzzzzz,
  xxzzzzzz,
  xyyyyyyy,
  xyyyyyyz,
  xyyyyyzz,
  xyyyyzzz,
  xyyyzzzz,
  xyyzzzzz,
  xyzzzzzz,
  xzzzzzzz,
  yyyyyyyy,
  yyyyyyyz,
  yyyyyyzz,
  yyyyyzzz,
  yyyyzzzz,
  yyyzzzzz,
  yyzzzzzz,
  yzzzzzzz,
  zzzzzzzz,
};
}

/* "superclass" AOSuperMatrix contains all common functionality for
 * atomic orbital matrix types
 */
class AOSuperMatrix {
 public:
  static int getBlockSize(int _lmax);
  static Eigen::MatrixXd getTrafo(const AOGaussianPrimitive& gaussian);
};

// base class for 1D atomic orbital matrix types (overlap, Coulomb, ESP)
template <class T>
class AOMatrix : public AOSuperMatrix {
 public:
  // Access functions
  int Dimension() { return _aomatrix.rows(); }
  const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& Matrix() const {
    return _aomatrix;
  }
  void Fill(const AOBasis& aobasis);
  // integrate F
  static std::vector<double> XIntegrate(int size, double U);

 protected:
  virtual void FillBlock(
      Eigen::Block<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>& matrix,
      const AOShell& shell_row, const AOShell& shell_col) const = 0;
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> _aomatrix;
};

/* base class class for 3D atomic orbital matrix types
 * (in principle, we could make nD and 1D and 3D are just special types)
 */
class AOMatrix3D : public AOSuperMatrix {
 public:
  const std::array<Eigen::MatrixXd, 3>& Matrix() const { return _aomatrix; }
  void Fill(const AOBasis& aobasis);

 protected:
  std::array<Eigen::MatrixXd, 3> _aomatrix;
  virtual void FillBlock(std::vector<Eigen::Block<Eigen::MatrixXd>>& matrix,
                         const AOShell& shell_row,
                         const AOShell& shell_col) const = 0;
};

/* derived class for atomic orbital gradient matrices, required for
 * momentum transition dipoles
 */
class AOMomentum : public AOMatrix3D {
 protected:
  void FillBlock(std::vector<Eigen::Block<Eigen::MatrixXd>>& matrix,
                 const AOShell& shell_row, const AOShell& shell_col) const;
};

/* derived class for atomic orbital electrical dipole matrices, required for
 * electical transition dipoles
 */
class AODipole : public AOMatrix3D {
 public:
  void setCenter(const Eigen::Vector3d& r) {
    _r = r;
  }  // definition of a center around which the moment should be calculated
 protected:
  void FillBlock(std::vector<Eigen::Block<Eigen::MatrixXd>>& matrix,
                 const AOShell& shell_row, const AOShell& shell_col) const;

 private:
  Eigen::Vector3d _r = Eigen::Vector3d::Zero();
};

// derived class for atomic orbital nuclear potential
class AOESP : public AOMatrix<double> {
 public:
  void Fillnucpotential(const AOBasis& aobasis, const QMMolecule& atoms);
  void Fillextpotential(
      const AOBasis& aobasis,
      const std::vector<std::unique_ptr<StaticSite>>& externalsites);
  const Eigen::MatrixXd& getNuclearpotential() const {
    return _nuclearpotential;
  }
  const Eigen::MatrixXd& getExternalpotential() const {
    return _externalpotential;
  }
  void setPosition(const Eigen::Vector3d& r) { _r = r; };

 protected:
  void FillBlock(Eigen::Block<Eigen::MatrixXd>& matrix,
                 const AOShell& shell_row, const AOShell& shell_col) const;

 private:
  Eigen::Vector3d _r = Eigen::Vector3d::Zero();
  Eigen::MatrixXd _nuclearpotential;
  Eigen::MatrixXd _externalpotential;
};

// derived class for Effective Core Potentials
class AOECP : public AOMatrix<double> {
 public:
  void setECP(const AOBasis* ecp) { _ecp = ecp; }

 protected:
  void FillBlock(Eigen::Block<Eigen::MatrixXd>& matrix,
                 const AOShell& shell_row, const AOShell& shell_col) const;

 private:
  const AOBasis* _ecp;
  Eigen::MatrixXd calcVNLmatrix(
      int lmax_ecp, const Eigen::Vector3d& posC,
      const AOGaussianPrimitive& g_row, const AOGaussianPrimitive& g_col,
      const Eigen::Matrix<int, 4, 5>& power_ecp,
      const Eigen::Matrix<double, 4, 5>& gamma_ecp,
      const Eigen::Matrix<double, 4, 5>& pref_ecp) const;

  void getBLMCOF(int lmax_ecp, int lmax_dft, const Eigen::Vector3d& pos,
                 tensor3d& BLC, tensor3d& C) const;
  Eigen::VectorXd CalcNorms(double decay, int size) const;
  Eigen::VectorXd CalcInt_r_exp(int nmax, double decay) const;
};

// derived class for kinetic energy
class AOKinetic : public AOMatrix<double> {
 protected:
  void FillBlock(Eigen::Block<Eigen::MatrixXd>& matrix,
                 const AOShell& shell_row, const AOShell& shell_col) const;
};

// derived class for atomic orbital overlap
class AOOverlap : public AOMatrix<double> {
 public:
  Eigen::MatrixXd FillShell(const AOShell& shell) const;
  int Removedfunctions() const { return removedfunctions; }
  double SmallestEigenValue() const { return smallestEigenvalue; }

  Eigen::MatrixXd Pseudo_InvSqrt(double etol);
  Eigen::MatrixXd Sqrt();

 protected:
  void FillBlock(Eigen::Block<Eigen::MatrixXd>& matrix,
                 const AOShell& shell_row, const AOShell& shell_col) const;

 private:
  int removedfunctions;
  double smallestEigenvalue;
};

class AODipole_Potential : public AOMatrix<double> {
 public:
  void Fillextpotential(
      const AOBasis& aobasis,
      const std::vector<std::unique_ptr<StaticSite>>& externalsites);
  const Eigen::MatrixXd& getExternalpotential() const {
    return _externalpotential;
  }

 protected:
  void FillBlock(Eigen::Block<Eigen::MatrixXd>& matrix,
                 const AOShell& shell_row, const AOShell& shell_col) const;

 private:
  void setSite(const StaticSite* site) { _site = site; };
  const StaticSite* _site;
  Eigen::MatrixXd _externalpotential;
};

class AOQuadrupole_Potential : public AOMatrix<double> {
 public:
  void Fillextpotential(
      const AOBasis& aobasis,
      const std::vector<std::unique_ptr<StaticSite>>& externalsites);
  const Eigen::MatrixXd& getExternalpotential() const {
    return _externalpotential;
  }

 protected:
  void FillBlock(Eigen::Block<Eigen::MatrixXd>& matrix,
                 const AOShell& shell_row, const AOShell& shell_col) const;

 private:
  void setSite(const StaticSite* site) { _site = site; };

  const StaticSite* _site;
  Eigen::MatrixXd _externalpotential;
};

// derived class for atomic orbital Coulomb interaction
class AOCoulomb : public AOMatrix<double> {
 public:
  Eigen::MatrixXd Pseudo_InvSqrt_GWBSE(const AOOverlap& auxoverlap,
                                       double etol);
  Eigen::MatrixXd Pseudo_InvSqrt(double etol);
  int Removedfunctions() const { return removedfunctions; }

 protected:
  void FillBlock(Eigen::Block<Eigen::MatrixXd>& matrix,
                 const AOShell& shell_row, const AOShell& shell_col) const;

 private:
  int removedfunctions;
};

class AOPlanewave : public AOMatrix<std::complex<double>> {
 public:
  void Fillextpotential(const AOBasis& aobasis,
                        const std::vector<Eigen::Vector3d>& kpoints);
  Eigen::MatrixXd getExternalpotential() const {
    return _externalpotential.real();
  }

 protected:
  void FillBlock(Eigen::Block<Eigen::MatrixXcd>& matrix,
                 const AOShell& shell_row, const AOShell& shell_col) const;

 private:
  void setkVector(const Eigen::Vector3d& k) { _k = k; };
  Eigen::Vector3d _k;
  Eigen::MatrixXcd _externalpotential;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_AOMATRIX_H
