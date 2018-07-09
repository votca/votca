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

        class GenCube : public xtp::QMTool {
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



            void calculateCube();
            void subtractCubes();

            string _orbfile;
            string _output_file;
            string _infile1;
            string _infile2;
            bool _do_groundstate;
            bool _do_bse;
            bool _do_qp;
            bool _do_transition;
            bool _do_singlet;
            bool _do_triplet;
            bool _do_ks;

            double _padding;
            int _xsteps;
            int _ysteps;
            int _zsteps;
            int _state;
            string _spin;
            string _type;
            string _mode;
            xtp::Logger _log;

        };

        void GenCube::Initialize(Property* options) {

            _do_groundstate = false;
            _do_bse = false;
            _do_qp = false;
            _do_transition = false;
            _do_singlet = false;
            _do_triplet = false;
            _do_ks = false;

            // update options with the VOTCASHARE defaults   
            UpdateWithDefaults(options, "xtp");


            string key = "options." + Identify();
            // _jobfile = options->get(key + ".file").as<string>();

            // key = "options." + Identify();


            // orbitals file or pure DFT output
            _orbfile = options->get(key + ".input").as<string> ();
            _output_file = options->get(key + ".output").as<string> ();

            // padding
            _padding = options->get(key + ".padding").as<double> ();

            // steps
            _xsteps = options->get(key + ".xsteps").as<int> ();
            _ysteps = options->get(key + ".ysteps").as<int> ();
            _zsteps = options->get(key + ".zsteps").as<int> ();

            _state = options->get(key + ".state").as<int> ();
            _spin = options->get(key + ".spin").as<string> ();
            if (_spin == "singlet") {
                _do_singlet = true;
            } else if (_spin == "triplet") {
                _do_triplet = true;
            } else {
                throw std::runtime_error("Spin not known, only singlet or triplet possible");
            }

            _type = options->get(key + ".type").as<string> ();

            _mode = options->get(key + ".mode").as<string> ();

            if (_mode == "subtract") {
                _infile1 = options->get(key + ".infile1").as<string> ();
                _infile2 = options->get(key + ".infile2").as<string> ();

            }


            if (_type == "ground" || _state == 0) _do_groundstate = true;
            else if (_type == "transition") _do_transition = true;
            else if (_type == "excited") {
                _do_bse = true;
                _do_groundstate = true;
            } else if (_type == "excited-gs") _do_bse = true;
            else if (_type == "qp") _do_qp = true;
            else if (_type == "ks") _do_ks = true;
            else {
                throw std::runtime_error("Option for type not known");
            }

            // get the path to the shared folders with xml files
            char *votca_share = getenv("VOTCASHARE");
            if (votca_share == NULL) throw std::runtime_error("VOTCASHARE not set, cannot open help files.");

        }

        void GenCube::calculateCube() {

            XTP_LOG(xtp::logDEBUG, _log) << "Reading serialized QM data from " << _orbfile << flush;

            Orbitals _orbitals;
            XTP_LOG(xtp::logDEBUG, _log) << " Loading QM data from " << _orbfile << flush;
            _orbitals.ReadFromCpt(_orbfile);

            if (_do_qp && !_orbitals.hasQPdiag()) {
                throw std::runtime_error("Orbitals file does not contain QP coefficients");
            }

            if (_do_bse && _do_singlet && !_orbitals.hasBSESinglets()) {
                throw std::runtime_error("Orbitals file does not contain Singlet BSE coefficients");
            }
            if (_do_bse && _do_triplet && !_orbitals.hasBSETriplets()) {
                throw std::runtime_error("Orbitals file does not contain Triplet BSE coefficients");
            }

            // load the QM data from serialized orbitals object
            // get atoms
            std::vector<QMAtom*> _atoms = _orbitals.QMAtoms();

            // determine min and max in each cartesian direction
            double xmin = std::numeric_limits<double>::max();
            double xmax = std::numeric_limits<double>::min();
            double ymin = xmin;
            double ymax = xmax;
            double zmin = xmin;
            double zmax = xmax;

            vector< QMAtom* > ::iterator ait;
            for (ait = _atoms.begin(); ait != _atoms.end(); ++ait) {
                const tools::vec& pos = (*ait)->getPos();
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

            // open output stream
            FILE *out;
            out = fopen(_output_file.c_str(), "w");

            // write cube header
            if (_do_groundstate && !_do_bse) {
                fprintf(out, "Electron density of neutral state \n");
            } else if (_do_bse) {
                if (_do_groundstate) {
                    fprintf(out, "Total electron density of excited state  %i spin %s \n", _state, _spin.c_str());
                } else {
                    fprintf(out, "Difference electron density of excited state  %i spin %s \n", _state, _spin.c_str());
                }
            } else if (_do_qp) {
                fprintf(out, "Quasiparticle state %i with energy %f eV \n", _state, _orbitals.QPdiagEnergies()[_state - 1 - _orbitals.getGWAmin()] * tools::conv::hrt2ev);
            } else if (_do_ks) {
                fprintf(out, "Kohn-Sham state %i with energy %f eV \n", _state, _orbitals.MOEnergies()[_state - 1] * tools::conv::hrt2ev);
            } else if (_do_transition) {
                fprintf(out, "Transition state  between Groundstate and state %i \n", _state);
            }
            fprintf(out, "Created by VOTCA-XTP \n");
            if (_do_qp) {
                fprintf(out, "-%lu %f %f %f \n", _atoms.size(), xstart, ystart, zstart);
            } else {
                fprintf(out, "%lu %f %f %f \n", _atoms.size(), xstart, ystart, zstart);
            }

            fprintf(out, "%d %f 0.0 0.0 \n", _xsteps + 1, xincr);
            fprintf(out, "%d 0.0 %f 0.0 \n", _ysteps + 1, yincr);
            fprintf(out, "%d 0.0 0.0 %f \n", _zsteps + 1, zincr);
            Elements _elements;
            for (ait = _atoms.begin(); ait != _atoms.end(); ++ait) {
                const tools::vec& pos = (*ait)->getPos();
                // get center coordinates in Bohr
                double x = pos.getX();
                double y = pos.getY();
                double z = pos.getZ();

                string element = (*ait)->getType();
                int atnum = _elements.getEleNum(element);
                double crg = (*ait)->getNuccharge();
                fprintf(out, "%d %f %f %f %f\n", atnum, crg, x, y, z);
            }

            if (_do_qp || _do_ks) {
                fprintf(out, "  1 %d \n", _state);
            }

            // load DFT basis set (element-wise information) from xml file
            BasisSet dftbs;
            dftbs.LoadBasisSet(_orbitals.getDFTbasis());
            XTP_LOG(xtp::logDEBUG, _log) << " Loaded DFT Basis Set " << _orbitals.getDFTbasis() << flush;

            // fill DFT AO basis by going through all atoms 
            AOBasis dftbasis;
            dftbasis.AOBasisFill(&dftbs, _orbitals.QMAtoms());

            // now depending on the type of cube
            if (_do_groundstate || _do_bse || _do_transition) {
                Eigen::MatrixXd DMAT_tot = Eigen::MatrixXd::Zero(dftbasis.AOBasisSize(), dftbasis.AOBasisSize());
                // ground state only if requested
                if (_do_groundstate) {
                    Eigen::MatrixXd DMATGS = _orbitals.DensityMatrixGroundState();
                    DMAT_tot = DMATGS; // Ground state + hole_contribution + electron contribution
                    XTP_LOG(xtp::logDEBUG, _log) << " Calculated ground state density matrix " << flush;
                }

                if (_state > 0) {
                    if (_do_transition) {
                        DMAT_tot = _orbitals.TransitionDensityMatrix(_spin, _state - 1);
                        XTP_LOG(xtp::logDEBUG, _log) << " Calculated transition state density matrix " << flush;
                    }                        // excited state if requested
                    else if (_do_bse) {
                        std::vector< Eigen::MatrixXd > DMAT = _orbitals.DensityMatrixExcitedState(_spin, _state - 1);

                        DMAT_tot = DMAT_tot + DMAT[1] - DMAT[0]; // Ground state + hole_contribution + electron contribution
                        XTP_LOG(xtp::logDEBUG, _log) << " Calculated excited state density matrix " << flush;
                    }
                }

                XTP_LOG(xtp::logDEBUG, _log) << " Calculating cube data ... \n" << flush;
                _log.setPreface(xtp::logDEBUG, (format(" ... ...")).str());

                boost::progress_display progress(_xsteps);
                // eval density at cube grid points
                for (int _ix = 0; _ix <= _xsteps; _ix++) {
                    double _x = xstart + double(_ix) * xincr;
                    for (int _iy = 0; _iy <= _ysteps; _iy++) {
                        double _y = ystart + double(_iy) * yincr;
                        int Nrecord = 0;
                        for (int _iz = 0; _iz <= _zsteps; _iz++) {
                            double _z = zstart + double(_iz) * zincr;
                            Nrecord++;
                            vec pos = vec(_x, _y, _z);
                            // get value of orbitals at each gridpoint
                            Eigen::VectorXd tmat = Eigen::VectorXd::Zero(dftbasis.AOBasisSize());
                            for (AOBasis::AOShellIterator _row = dftbasis.firstShell(); _row != dftbasis.lastShell(); _row++) {
                                const double decay = (*_row)->getMinDecay();
                                const tools::vec& shellpos = (*_row)->getPos();
                                tools::vec dist = shellpos - pos;
                                double distsq = dist*dist;
                                // if contribution is smaller than -ln(1e-10), calc density
                                if ((decay * distsq) < 20.7) {
                                    Eigen::VectorBlock<Eigen::VectorXd> tmat_block = tmat.segment((*_row)->getStartIndex(), (*_row)->getNumFunc());
                                    (*_row)->EvalAOspace(tmat_block, pos);
                                }
                            }
                            double density_at_grid = (tmat.transpose() * DMAT_tot * tmat).value();
                            if (Nrecord == 6 || _iz == _zsteps) {
                                fprintf(out, "%E \n", density_at_grid);
                                Nrecord = 0;
                            } else {
                                fprintf(out, "%E ", density_at_grid);
                            }
                        }// z-component
                    }// y-component
                    ++progress;
                } // x-component
            } // ground or excited state
            _log.setPreface(xtp::logDEBUG, (format("\n ... ...")).str());

            // diagonalized QP, if requested
            if ((_do_ks || _do_qp) && _state > 0) {
                Eigen::VectorXd Ftemp;
                if (_do_qp) {
                    Eigen::VectorXd QPcoefs = _orbitals.QPdiagCoefficients().col(_state - 1);
                    // get QPdiag coefficients for the requested state
                    Eigen::MatrixXd MOs = _orbitals.MOCoefficients().block(_orbitals.getGWAmin(), 0, _orbitals.getGWAtot(), dftbasis.AOBasisSize());
                    // get DFT MO coefficients
                    Ftemp = MOs.transpose() * QPcoefs;
                }
                if (_do_ks) {
                    Ftemp = _orbitals.MOCoefficients().col(_state - 1);
                }
                for (int _ix = 0; _ix <= _xsteps; _ix++) {
                    double _x = xstart + double(_ix) * xincr;
                    for (int _iy = 0; _iy <= _ysteps; _iy++) {
                        double _y = ystart + double(_iy) * yincr;
                        int Nrecord = 0;
                        for (int _iz = 0; _iz <= _zsteps; _iz++) {
                            double _z = zstart + double(_iz) * zincr;
                            Nrecord++;
                            // get value of orbitals at each gridpoint
                            Eigen::VectorXd tmat = Eigen::VectorXd::Zero(dftbasis.AOBasisSize());
                            vec pos = vec(_x, _y, _z);
                            for (AOBasis::AOShellIterator _row = dftbasis.firstShell(); _row != dftbasis.lastShell(); _row++) {
                                const double decay = (*_row)->getMinDecay();
                                const tools::vec& shellpos = (*_row)->getPos();
                                tools::vec dist = shellpos - pos;
                                double distsq = dist*dist;
                                // if contribution is smaller than -ln(1e-10), calc density
                                if ((decay * distsq) < 20.7) {
                                    Eigen::VectorBlock<Eigen::VectorXd> tmat_block = tmat.segment((*_row)->getStartIndex(), (*_row)->getNumFunc());
                                    (*_row)->EvalAOspace(tmat_block, pos);
                                }
                            }
                            double QP_at_grid = (Ftemp.transpose() * tmat).value();
                            if (Nrecord == 6 || _iz == _zsteps) {
                                fprintf(out, "%E \n", QP_at_grid);
                                Nrecord = 0;
                            } else {
                                fprintf(out, "%E ", QP_at_grid);
                            }
                        }// z-component
                    }// y-component
                } // x-component
            }

            fclose(out);
            XTP_LOG(xtp::logDEBUG, _log) << "Wrote cube data to " << _output_file << flush;

            return;
        }

        void GenCube::subtractCubes() {

            // open infiles for reading
            ifstream in1;
            XTP_LOG(xtp::logDEBUG, _log) << " Reading first cube from " << _infile1 << flush;
            in1.open(_infile1.c_str(), ios::in);
            ifstream in2;
            XTP_LOG(xtp::logDEBUG, _log) << " Reading second cube from " << _infile2 << flush;
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
            if (natoms < 0) _do_qp = true;
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

            if (_do_qp) {
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
            XTP_LOG(xtp::logDEBUG, _log) << "Wrote subtracted cube data to " << _output_file << flush;

        }

        bool GenCube::Evaluate() {

            _log.setReportLevel(xtp::logDEBUG);
            _log.setMultithreading(true);

            _log.setPreface(xtp::logINFO, "\n... ...");
            _log.setPreface(xtp::logERROR, "\n... ...");
            _log.setPreface(xtp::logWARNING, "\n... ...");
            _log.setPreface(xtp::logDEBUG, "\n... ...");

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
