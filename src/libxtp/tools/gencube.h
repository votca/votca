/* 
 *            Copyright 2009-2018 The VOTCA Development Team
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

#ifndef _VOTCA_XTP_GENCUBE_H
#define _VOTCA_XTP_GENCUBE_H
#include <boost/progress.hpp>
#include <stdio.h>
#include <votca/xtp/aobasis.h>
#include <boost/format.hpp>
#include <votca/tools/elements.h>
#include <votca/xtp/logger.h>
#include <votca/tools/constants.h>

namespace votca {
  namespace xtp {

    using namespace std;

    class GenCube : public QMTool {
      public:

        GenCube() {
        };

        ~GenCube() {
        };

        string Identify() {
          return "gencube";
        }

        void Initialize(Property *options);
        bool Evaluate();


      private:

        Eigen::VectorXd EvaluateBasisAtPosition(const AOBasis& dftbasis,const tools::vec& pos);

        void calculateCube();
        void subtractCubes();

        string _orbfile;
        string _output_file;
        string _infile1;
        string _infile2;

        bool _dostateonly;

        double _padding;
        int _xsteps;
        int _ysteps;
        int _zsteps;
        QMState _state;
        string _mode;
        Logger _log;

    };

    void GenCube::Initialize(Property* options) {

      // update options with the VOTCASHARE defaults   
      UpdateWithDefaults(options, "xtp");


      string key = "options." + Identify();
      _orbfile = options->get(key + ".input").as<string> ();
      _output_file = options->get(key + ".output").as<string> ();

      // padding
      _padding = options->get(key + ".padding").as<double> ();

      // steps
      _xsteps = options->get(key + ".xsteps").as<int> ();
      _ysteps = options->get(key + ".ysteps").as<int> ();
      _zsteps = options->get(key + ".zsteps").as<int> ();

      std::string statestring=options->get(key + ".state").as<string> ();
      _state.FromString(statestring);
      _dostateonly = options->ifExistsReturnElseReturnDefault<bool>(key + ".diff2gs",false);

      _mode = options->get(key + ".mode").as<string> ();
      if (_mode == "subtract") {
        _infile1 = options->get(key + ".infile1").as<string> ();
        _infile2 = options->get(key + ".infile2").as<string> ();
      }

      // get the path to the shared folders with xml files
      char *votca_share = getenv("VOTCASHARE");
      if (votca_share == NULL) throw std::runtime_error("VOTCASHARE not set, cannot open help files.");

    }

    void GenCube::calculateCube() {

      XTP_LOG(logDEBUG, _log) << "Reading serialized QM data from " << _orbfile << flush;

      Orbitals orbitals;
      XTP_LOG(logDEBUG, _log) << " Loading QM data from " << _orbfile << flush;
      orbitals.ReadFromCpt(_orbfile);

      std::vector<QMAtom*> atoms = orbitals.QMAtoms();

      // determine min and max in each cartesian direction
      double xmin = std::numeric_limits<double>::max();
      double xmax = std::numeric_limits<double>::min();
      double ymin = xmin;
      double ymax = xmax;
      double zmin = xmin;
      double zmax = xmax;

      for (const QMAtom* atom:atoms) {
        const tools::vec& pos = atom->getPos();
        // get center coordinates in Bohr
        double x = pos.getX();
        double y = pos.getY();
        double z = pos.getZ();
        if (x > xmax) xmax = x;
        if (x < xmin) xmin = x;
        if (y > ymax) ymax = y;
        if (y < ymin) ymin = y;
        if (z > zmax) zmax = z;
        if (z < zmin) zmin = z;
      }
      // generate cube grid
      double xstart = xmin - _padding;
      double xstop = xmax + _padding;
      double ystart = ymin - _padding;
      double ystop = ymax + _padding;
      double zstart = zmin - _padding;
      double zstop = zmax + _padding;

      double xincr = (xstop - xstart) / double(_xsteps);
      double yincr = (ystop - ystart) / double(_ysteps);
      double zincr = (zstop - zstart) / double(_zsteps);

      std::ofstream out(_output_file);
      if (!out.is_open()) {
        throw std::runtime_error("Bad file handle: " + _output_file);
      }

      // write cube header
      if(_state.isTransition()){
        out<<boost::format("Transition state: %1$s \n")% _state.ToString();
      }

      bool do_amplitude=(_state.Type().isSingleParticleState());

      if(do_amplitude){
        out<<boost::format("%1$s with energy %2$f eV \n") %_state.ToLongString() % (orbitals.getExcitedStateEnergy(_state)*tools::conv::hrt2ev);
      }else{
        if(_dostateonly){
          out<<boost::format( "Difference electron density of excited state %1$s \n")% _state.ToString();
        }else{
          out<<boost::format("Total electron density of %1$s state\n")% _state.ToLongString();
        }
      }


      out<< "Created by VOTCA-XTP \n";
      if (do_amplitude) {
        out<<boost::format("-%1$lu %2$f %3$f %4$f \n") %atoms.size() %xstart %ystart %zstart;
      } else {
        out<<boost::format("%1$lu %2$f %3$f %4$f \n") %atoms.size() %xstart %ystart %zstart;
      }

      out<<boost::format( "%1$d %2$f 0.0 0.0 \n") %(_xsteps + 1) % xincr;
      out<<boost::format( "%1$d 0.0 %2$f 0.0 \n") %(_ysteps + 1) % yincr;
      out<<boost::format( "%1$d 0.0 0.0 %2$f \n") %(_zsteps + 1) % zincr;
      Elements _elements;
      for (const QMAtom* atom:atoms) {
        const tools::vec& pos = atom->getPos();
        // get center coordinates in Bohr
        double x = pos.getX();
        double y = pos.getY();
        double z = pos.getZ();
        string element = atom->getElement();
        int atnum = _elements.getEleNum(element);
        double crg = atom->getNuccharge();
        out<<boost::format( "%1$d %2$f %3$f %4$f %5$f\n") %atnum %crg %x %y %z;
      }

      if (do_amplitude) {
        out<<boost::format("  1 %1$d \n") %(_state.Index()+1);
      }

      // load DFT basis set (element-wise information) from xml file
      BasisSet dftbs;
      dftbs.LoadBasisSet(orbitals.getDFTbasis());
      XTP_LOG(logDEBUG, _log) << " Loaded DFT Basis Set " << orbitals.getDFTbasis() << flush;

      // fill DFT AO basis by going through all atoms 
      AOBasis dftbasis;
      dftbasis.AOBasisFill(dftbs, orbitals.QMAtoms());

      Eigen::MatrixXd mat=Eigen::MatrixXd::Zero(dftbasis.AOBasisSize(), dftbasis.AOBasisSize());
      if(_dostateonly){
        if(_state.Type().isExciton()){
          std::vector< Eigen::MatrixXd > DMAT=orbitals.DensityMatrixExcitedState(_state);
          mat=DMAT[1] - DMAT[0];
        }
      }else{
        mat=orbitals.DensityMatrixFull(_state);
      }
      int amplitudeindex=0;
      if(do_amplitude){
        if(_state.Type()==QMStateType::DQPstate){
          mat=orbitals.CalculateQParticleAORepresentation();
          amplitudeindex=_state.Index()-orbitals.getGWAmin();
        }else{
          mat=orbitals.MOCoefficients();
          amplitudeindex=_state.Index();
        }
      }


      XTP_LOG(logDEBUG, _log) << " Calculating cube data ... \n" << flush;
      _log.setPreface(logDEBUG, (format(" ... ...")).str());

      boost::progress_display progress(_xsteps);
      // eval density at cube grid points
      for (int ix = 0; ix <= _xsteps; ix++) {
        double x = xstart + double(ix) * xincr;
        for (int iy = 0; iy <= _ysteps; iy++) {
          double y = ystart + double(iy) * yincr;
          int Nrecord = 0;
          for (int iz = 0; iz <= _zsteps; iz++) {
            double z = zstart + double(iz) * zincr;
            Nrecord++;
            vec pos = vec(x, y, z);
            Eigen::VectorXd tmat=EvaluateBasisAtPosition(dftbasis,pos);
            double value=0.0;
            if(do_amplitude){
              value = (mat.col(amplitudeindex).transpose() * tmat).value();
            }else{  
              value = (tmat.transpose() * mat * tmat).value();
            }
            if (Nrecord == 6 || iz == _zsteps) {
              out<<boost::format("%1$E \n")% value;
              Nrecord = 0;
            } else {
              out<<boost::format("%1$E ")% value;
            }
          }// z-component
        }// y-component
        ++progress;
      } // x-component


      out.close();
      XTP_LOG(logDEBUG, _log) << "Wrote cube data to " << _output_file << flush;
      return;
    }

    Eigen::VectorXd GenCube::EvaluateBasisAtPosition(const AOBasis& dftbasis,const tools::vec& pos){

      // get value of orbitals at each gridpoint
      Eigen::VectorXd tmat = Eigen::VectorXd::Zero(dftbasis.AOBasisSize());
      for (const AOShell* shell:dftbasis) {
        const double decay =shell->getMinDecay();
        const tools::vec& shellpos = shell->getPos();
        tools::vec dist = shellpos - pos;
        double distsq = dist*dist;
        // if contribution is smaller than -ln(1e-10), calc density
        if ((decay * distsq) < 20.7) {
          Eigen::VectorBlock<Eigen::VectorXd> tmat_block = tmat.segment(shell->getStartIndex(), shell->getNumFunc());
          shell->EvalAOspace(tmat_block, pos);
        }
      }
      return tmat;
    }

    void GenCube::subtractCubes() {

      // open infiles for reading
      ifstream in1;
      XTP_LOG(logDEBUG, _log) << " Reading first cube from " << _infile1 << flush;
      in1.open(_infile1.c_str(), ios::in);
      ifstream in2;
      XTP_LOG(logDEBUG, _log) << " Reading second cube from " << _infile2 << flush;
      in2.open(_infile2.c_str(), ios::in);
      string s;

      FILE *out;
      out = fopen(_output_file.c_str(), "w");

      // first two lines of header are garbage
      getline(in1, s);
      fprintf(out, "%s\n", s.c_str());
      getline(in1, s);
      fprintf(out, "%s subtraction \n", s.c_str());
      getline(in2, s);
      getline(in2, s);

      // read rest from header
      int natoms;
      double xstart;
      double ystart;
      double zstart;
      // first line
      in1 >> natoms;
      bool do_amplitude=false;
      if (natoms < 0) do_amplitude = true;
      in1 >> xstart;
      in1 >> ystart;
      in1 >> zstart;
      // check from second file
      int tempint;
      double tempdouble;
      in2 >> tempint;
      if (tempint != natoms) {
        cerr << "Atom numbers do not match";
        exit(1);
      }
      in2 >> tempdouble;
      if (tempdouble != xstart) {
        cerr << "Xstart does not match";
        exit(1);
      }
      in2 >> tempdouble;
      if (tempdouble != ystart) {
        cerr << "Ystart does not match";
        exit(1);
      }
      in2 >> tempdouble;
      if (tempdouble != zstart) {
        cerr << "Zstart does not match";
        exit(1);
      }

      fprintf(out, "%d %f %f %f \n", natoms, xstart, ystart, zstart);

      // grid information from first cube
      double xincr;
      double yincr;
      double zincr;
      in1 >> _xsteps;
      in1 >> xincr;
      in1 >> tempdouble;
      in1 >> tempdouble;
      in1 >> _ysteps;
      in1 >> tempdouble;
      in1 >> yincr;
      in1 >> tempdouble;
      in1 >> _zsteps;
      in1 >> tempdouble;
      in1 >> tempdouble;
      in1 >> zincr;

      // check second cube
      in2 >> tempint;
      if (tempint != _xsteps) {
        cerr << "xsteps does not match";
        exit(1);
      }
      in2 >> tempdouble;
      if (tempdouble != xincr) {
        cerr << "xincr does not match";
        exit(1);
      }
      in2 >> tempdouble;
      in2 >> tempdouble;
      in2 >> tempint;
      if (tempint != _ysteps) {
        cerr << "ysteps does not match";
        exit(1);
      }
      in2 >> tempdouble;
      in2 >> tempdouble;
      if (tempdouble != yincr) {
        cerr << "yincr does not match";
        exit(1);
      }
      in2 >> tempdouble;
      in2 >> tempint;
      if (tempint != _zsteps) {
        cerr << "zsteps does not match";
        exit(1);
      }
      in2 >> tempdouble;
      in2 >> tempdouble;
      in2 >> tempdouble;
      if (tempdouble != zincr) {
        cerr << "zincr does not match";
        exit(1);
      }

      fprintf(out, "%d %f 0.0 0.0 \n", _xsteps, xincr);
      fprintf(out, "%d 0.0 %f 0.0 \n", _ysteps, yincr);
      fprintf(out, "%d 0.0 0.0 %f \n", _zsteps, zincr);
      // atom information

      for (int iatom = 0; iatom < std::abs(natoms); iatom++) {
        // get center coordinates in Bohr
        double x;
        double y;
        double z;
        int atnum;
        double crg;

        // get from first cube
        in1 >> atnum;
        in1 >> crg;
        in1 >> x;
        in1 >> y;
        in1 >> z;

        // check second cube
        in2 >> tempint;
        if (tempint != atnum) {
          cerr << "atnum does not match";
          exit(1);
        }
        in2 >> tempdouble;
        if (tempdouble != crg) {
          cerr << "crg does not match";
          exit(1);
        }
        in2 >> tempdouble;
        if (tempdouble != x) {
          cerr << "x does not match";
          exit(1);
        }
        in2 >> tempdouble;
        if (tempdouble != y) {
          cerr << "y does not match";
          exit(1);
        }
        in2 >> tempdouble;
        if (tempdouble != z) {
          cerr << "z does not match";
          exit(1);
        }

        fprintf(out, "%d %f %f %f %f\n", atnum, crg, x, y, z);


      }

      if (do_amplitude) {
        int ntotal;
        int nis;
        in1 >> ntotal;
        in1 >> nis;

        in2 >> tempint;
        if (tempint != ntotal) {
          cerr << "ntotal does not match";
          exit(1);
        }
        in2 >> tempint;
        if (tempint != nis) {
          cerr << "nis does not match";
          exit(1);
        }
        fprintf(out, "  1 %d \n", nis);
      }
      // now read data
      double val1;
      double val2;
      for (int _ix = 0; _ix < _xsteps; _ix++) {
        for (int _iy = 0; _iy < _ysteps; _iy++) {
          int Nrecord = 0;
          for (int _iz = 0; _iz < _zsteps; _iz++) {
            Nrecord++;
            in1 >> val1;
            in2 >> val2;
            if (Nrecord == 6 || _iz == _zsteps - 1) {
              fprintf(out, "%E \n", val1 - val2);
              Nrecord = 0;
            } else {
              fprintf(out, "%E ", val1 - val2);
            }
          }
        }
      }

      fclose(out);
      XTP_LOG(logDEBUG, _log) << "Wrote subtracted cube data to " << _output_file << flush;
    }

    bool GenCube::Evaluate() {

      _log.setReportLevel(logDEBUG);
      _log.setMultithreading(true);

      _log.setPreface(logINFO, "\n... ...");
      _log.setPreface(logERROR, "\n... ...");
      _log.setPreface(logWARNING, "\n... ...");
      _log.setPreface(logDEBUG, "\n... ...");

      // calculate new cube
      if (_mode == "new") {
        calculateCube();
      } else if (_mode == "subtract") {
        subtractCubes();
      }

      return true;
    }
  }
}


#endif
