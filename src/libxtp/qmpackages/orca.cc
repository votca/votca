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

#include "orca.h"
#include <votca/xtp/qminterface.h>
#include <votca/ctp/segment.h>
#include <votca/xtp/basisset.h>
#include <votca/tools/elements.h>
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <stdio.h>
#include <iomanip>





namespace votca {
    namespace xtp {
      using namespace std;

        void Orca::Initialize(tools::Property &options) {

            //good luck

            // Orca file names
            std::string fileName = "system";

            _input_file_name = fileName + ".inp";
            _log_file_name = fileName + ".log";
            _shell_file_name = fileName + ".sh";
            _orb_file_name = fileName + ".gbw";

            std::string key = "package";
            std::string _name = options.get(key + ".name").as<std::string> ();

            if (_name != "orca") {
                cerr << "Tried to use " << _name << " package. ";
                throw std::runtime_error("Wrong options file");
            }

            _executable = options.get(key + ".executable").as<std::string> ();
            _charge = options.get(key + ".charge").as<int> ();
            _spin = options.get(key + ".spin").as<int> ();
            _options = options.get(key + ".options").as<std::string> ();
            _memory = options.get(key + ".memory").as<std::string> ();
            _threads = options.get(key + ".threads").as<int> ();
            _scratch_dir = options.get(key + ".scratch").as<std::string> ();
            _cleanup = options.get(key + ".cleanup").as<std::string> ();
            _auxbasisset_name = options.get(key + ".auxbasisset").as<std::string> ();



            if (options.exists(key + ".outputVxc")) {
                _output_Vxc = options.get(key + ".outputVxc").as<bool> ();
            } else _output_Vxc = false;
            if (_output_Vxc) {
                throw std::runtime_error("Sorry " + _name + " does not support Vxc output");
            }

            if ( _write_pseudopotentials )  _ecp_name = options.get(key + ".ecp").as<std::string> ();
            _basisset_name = options.get(key + ".basisset").as<std::string> ();
            _write_basis_set = options.get(key + ".writebasisset").as<bool> ();
            _write_pseudopotentials = options.get(key + ".writepseudopotentials").as<bool> ();
            if ( _write_pseudopotentials )  _ecp_name = options.get(key + ".ecp").as<std::string> ();

            // check if the optimize keyword is present, if yes, read updated coords
            std::string::size_type iop_pos = _options.find(" Opt"); /*optimization word in orca*/
            if (iop_pos != std::string::npos) {
                _is_optimization = true;
            } else {
                _is_optimization = false;
            }

            // check if the esp keyword is present, if yes, get the charges and save them
            iop_pos = _options.find(" chelpg"); /*electrostatic potential related to partial atomic charges I guess is chelpg in orca but check */
            if (iop_pos != std::string::npos ||  _options.find(" CHELPG")!= std::string::npos) {
                _get_charges = true;
            } else {
                _get_charges = false;
            }

            // check if the guess should be prepared, if yes, append the guess later
            _write_guess = false;
            iop_pos = _options.find("Guess MORead");
            if (iop_pos != std::string::npos) _write_guess = true;
            iop_pos = _options.find("Guess MORead\n");
            if (iop_pos != std::string::npos) _write_guess = true;
        }


/* Custom basis sets are written on a per-element basis to
     * the system.bas/aux file(s), which are then included in the
     * Orca input file using GTOName = "system.bas/aux"
     */
    void Orca::WriteBasisset(std::vector<QMAtom*>& qmatoms, std::string& bs_name, std::string& el_file_name) {

      
      std::vector<std::string> UniqueElements= FindUniqueElements(qmatoms);
      
      tools::Elements elementInfo;
      BasisSet bs;
      bs.LoadBasisSet(bs_name);
      CTP_LOG(ctp::logDEBUG, *_pLog) << "Loaded Basis Set " << bs_name << flush;
      ofstream el_file;

      el_file.open(el_file_name.c_str());
      el_file << "$DATA" << endl;

        for (const std::string& element_name:UniqueElements) {
          const Element& element = bs.getElement(element_name);
          el_file << elementInfo.getEleFull(element_name) << endl;
          for (const Shell& shell:element) {
            string type = shell.getType();
            // check combined shells
            for (unsigned i = 0; i < type.size(); ++i) {
              string subtype = string(type, i, 1);
              el_file << subtype << " " << shell.getSize() << endl;
              int sh_idx = 0;
              for (const GaussianPrimitive& gaussian:shell) {
                sh_idx++;
                el_file << " " << sh_idx << " " << indent(gaussian._decay);
                el_file << " " << indent(gaussian._contraction[FindLmax(subtype)]);
                el_file << endl;
              }

            }
          }

        }
      el_file << "STOP\n";
      el_file.close();

      return;
    }

/* Coordinates are written in standard Element,x,y,z format to the
     * input file.
     */
    void Orca::WriteCoordinates(std::ofstream& inp_file, std::vector<QMAtom*>& qmatoms) {

      for (QMAtom* atom : qmatoms) {
        tools::vec pos = atom->getPos() * tools::conv::bohr2ang;
        inp_file << setw(3) << atom->getType().c_str()
                << setw(12) << setiosflags(ios::fixed) << setprecision(5) << pos.getX()
                << setw(12) << setiosflags(ios::fixed) << setprecision(5) << pos.getY()
                << setw(12) << setiosflags(ios::fixed) << setprecision(5) << pos.getZ()
                << endl;
      }
      inp_file << "* \n" << endl;
      return;
    }


/* If custom ECPs are used, they need to be specified in the input file
     * in a section following the basis set includes.
     */
    void Orca::WriteECP(std::ofstream& inp_file, std::vector<QMAtom*>& qmatoms) {

      inp_file << endl;
      std::vector<std::string> UniqueElements= FindUniqueElements(qmatoms);
           
      BasisSet ecp;
      ecp.LoadPseudopotentialSet(_ecp_name);

      CTP_LOG(ctp::logDEBUG, *_pLog) << "Loaded Pseudopotentials " << _ecp_name << flush;

      for (const std::string& element_name:UniqueElements) {
          try{    
           ecp.getElement(element_name);
          }catch(std::runtime_error& error){
            CTP_LOG(ctp::logDEBUG, *_pLog) << "No pseudopotential for " << element_name<<" available" << flush;
            continue;
          }
          const Element& element = ecp.getElement(element_name);
 
          inp_file << "\n" << "NewECP" << " " << element_name << endl;
          inp_file << "N_core" << " " << element.getNcore() << endl;
          inp_file << "lmax" << " " << getLName(element.getLmax()) << endl;
          //For Orca the order doesn't matter but let's write it in ascending order
          // write remaining shells in ascending order s,p,d...
          for (int i = 0; i <= element.getLmax(); i++) {
              for (const Shell& shell:element) {
                if (shell.getLmax() == i) {
                // shell type, number primitives, scale factor
                inp_file << shell.getType() << " " << shell.getSize() << endl;
                int sh_idx = 0;
                 for (const GaussianPrimitive& gaussian:shell) {
                  sh_idx++;
                  inp_file << sh_idx << " " << gaussian._decay << " " << gaussian._contraction[0] << " " << gaussian._power << endl;
                }
              }
            }
          }
          inp_file << "end\n " << "\n" << endl;
        }
      return;
    }
        
        void Orca::WriteChargeOption(){
           std::string::size_type iop_pos = _options.find("pointcharges");
              if (iop_pos == std::string::npos) {
                _options = _options + "\n %pointcharges \"background.crg\"";
              }
        }

        /* For QM/MM the molecules in the MM environment are represented by
         * their atomic partial charge distributions. ORCA expects them in
         * q,x,y,z format in a separate file "background.crg"
         */
        void Orca::WriteBackgroundCharges() {
          
            std::ofstream crg_file;
            std::string _crg_file_name_full = _run_dir + "/background.crg";
            crg_file.open(_crg_file_name_full.c_str());
            int total_background = 0;

            for (std::shared_ptr<ctp::PolarSeg> seg:_PolarSegments) {
                for (ctp::APolarSite* site:*seg) {
                    if (site->getQ00() != 0.0) total_background++;
                    if (site->getRank() > 0 || _with_polarization ) {
                        std::vector<std::vector<double>> _split_multipoles = SplitMultipoles(site);
                        total_background+= _split_multipoles.size();
                    }
                }
            } //counting only
            
            crg_file << total_background << endl;
            boost::format fmt("%1$+1.7f %2$+1.7f %3$+1.7f %4$+1.7f");
            //now write
            for (std::shared_ptr<ctp::PolarSeg> seg:_PolarSegments) {
                for (ctp::APolarSite* site:*seg) {
                    string sitestring=boost::str(fmt % site->getQ00() % ((site->getPos().getX())*votca::tools::conv::nm2ang) 
                            % (site->getPos().getY()*votca::tools::conv::nm2ang) 
                            % (site->getPos().getZ()*votca::tools::conv::nm2ang) 
                            );
                    if (site->getQ00() != 0.0) crg_file << sitestring << endl;
                    if (site->getRank() > 0 || _with_polarization ) {
                        std::vector< std::vector<double> > _split_multipoles = SplitMultipoles(site);
                        for (const auto& mpoles:_split_multipoles){
                           string multipole=boost::str( fmt % mpoles[3] % mpoles[0] % mpoles[1] % mpoles[2] );
                            crg_file << multipole << endl;
                        }
                    }
                }
            }
            
            return;
        }

        /**
         * Prepares the *.inp file from a vector of segments
         * Appends a guess constructed from monomer orbitals if supplied, Not implemented yet
         */
        bool Orca::WriteInputFile(Orbitals& orbitals) {

            std::vector<std::string> results;
            std::string temp_suffix = "/id";
            std::string scratch_dir_backup = _scratch_dir;
            std::ofstream inp_file;
            std::string inp_file_name_full = _run_dir + "/" + _input_file_name;
            inp_file.open(inp_file_name_full.c_str());
            // header
            inp_file << "* xyz  " << _charge << " " << _spin << endl;

            std::vector< QMAtom* > qmatoms = orbitals.QMAtoms();
            // put coordinates
            WriteCoordinates(inp_file, qmatoms);
            // add parallelization info
            inp_file << "%pal\n " << "nprocs " << _threads << "\nend" << "\n" << endl;
            // basis set info
            if (_write_basis_set) {
                std::string el_file_name = _run_dir + "/" + "system.bas";
                WriteBasisset(qmatoms, _basisset_name, el_file_name);
                inp_file << "%basis\n " << endl;
                inp_file << "GTOName" << " " << "=" << "\"system.bas\";" << endl;
                if (_auxbasisset_name != "") {
                    std::string aux_file_name = _run_dir + "/" + "system.aux";
                    WriteBasisset(qmatoms, _auxbasisset_name, aux_file_name);
                    inp_file << "GTOAuxName" << " " << "=" << "\"system.aux\";" << endl;
                }
            } // write_basis set

            // ECPs
            /* WRITING ECP INTO system.inp FILE for ORCA**/
            if (_write_pseudopotentials) {
                WriteECP(inp_file, qmatoms);
            } // write pseudopotentials
            /* END   OF WRITING BASISSET/ECP INTO system.inp FILE for ORCA*************/
            inp_file << "end\n " << "\n" << endl; //This end is for the basis set block
            if (_write_charges) {
                WriteBackgroundCharges();
            }

            inp_file << _options << "\n";
            inp_file << endl;
            inp_file.close();
            // and now generate a shell script to run both jobs, if neccessary
            CTP_LOG(ctp::logDEBUG, *_pLog) << "Setting the scratch dir to " << _scratch_dir + temp_suffix << flush;
            _scratch_dir = scratch_dir_backup + temp_suffix;
            WriteShellScript();
            _scratch_dir = scratch_dir_backup;
            return true;
        }

        bool Orca::WriteShellScript() {
            ofstream shell_file;
            std::string shell_file_name_full = _run_dir + "/" + _shell_file_name;
            shell_file.open(shell_file_name_full.c_str());
            shell_file << "#!/bin/bash" << endl;
            shell_file << "mkdir -p " << _scratch_dir << endl;
            
            if(_write_guess){
              if(!(boost::filesystem::exists( _run_dir + "/molA.gbw" ) && boost::filesystem::exists( _run_dir + "/molB.gbw" )  )){
              throw runtime_error("Using guess relies on a molA.gbw and a molB.gbw file being in the directory.");
            }
              shell_file<<_executable<<"_mergefrag molA.gbw molB.gbw dimer.gbw > merge.log"<<endl;
            }
            shell_file << _executable << " " << _input_file_name << " > " << _log_file_name << endl; //" 2> run.error" << endl;
            shell_file.close();
            return true;
        }

        /**
         * Runs the Orca job.
         */
        bool Orca::Run() {

            CTP_LOG(ctp::logDEBUG, *_pLog) << "Running Orca job" << flush;

            if (std::system(NULL)) {
                std::string _command = "cd " + _run_dir + "; sh " + _shell_file_name;
                int check = std::system(_command.c_str());
                if (check == -1) {
                    CTP_LOG(ctp::logERROR, *_pLog) << _input_file_name << " failed to start" << flush;
                    return false;
                }
                if (CheckLogFile()) {
                    CTP_LOG(ctp::logDEBUG, *_pLog) << "Finished Orca job" << flush;
                    return true;
                } else {
                    CTP_LOG(ctp::logDEBUG, *_pLog) << "Orca job failed" << flush;
                }
            } else {
                CTP_LOG(ctp::logERROR, *_pLog) << _input_file_name << " failed to start" << flush;
                return false;
            }

            return true;
        }

        /**
         * Cleans up after the Orca job
         */
        void Orca::CleanUp() {

          if(_write_guess){
            remove((_run_dir + "/" +"molA.gbw").c_str());
            remove((_run_dir + "/" +"molB.gbw").c_str());
            remove((_run_dir + "/" +"dimer.gbw").c_str());           
          }
            // cleaning up the generated files
            if (_cleanup.size() != 0) {
                tools::Tokenizer tok_cleanup(_cleanup, ",");
                std::vector <std::string> cleanup_info;
                tok_cleanup.ToVector(cleanup_info);
                for (const std::string& substring:cleanup_info) {
                    if (substring == "inp") {
                        std::string file_name = _run_dir + "/" + _input_file_name;
                        remove(file_name.c_str());
                    }

                    if (substring == "bas") {
                        std::string file_name = _run_dir + "/system.bas";
                        remove(file_name.c_str());
                    }

                    if (substring == "log") {
                        std::string file_name = _run_dir + "/" + _log_file_name;
                        remove(file_name.c_str());
                    }

                    if (substring == "gbw") {
                        std::string file_name = _run_dir + "/" + _orb_file_name;
                        remove(file_name.c_str());
                    }

                    if (substring == "ges") {
                        std::string file_name = _run_dir + "/system.ges";
                        remove(file_name.c_str());
                    }
                    if (substring == "prop") {
                        std::string file_name = _run_dir + "/system.prop";
                        remove(file_name.c_str());
                    }
                }
            }
            return;
        }

        bool Orca::ParseLogFile(Orbitals& orbitals) {
            bool found_success=false;
            orbitals.setQMpackage("orca");
            orbitals.setDFTbasisName(_basisset_name);
            if (_write_pseudopotentials) {
                orbitals.setECPName(_ecp_name);
            } 
            CTP_LOG(ctp::logDEBUG, *_pLog) << "Parsing " << _log_file_name << flush;
            std::string log_file_name_full = _run_dir + "/" + _log_file_name;
            // check if LOG file is complete
            if (!CheckLogFile()) return false;
            std::map <int, double> energies;
            std::map <int, double> occupancy;

            std::string line;
            unsigned levels = 0;
            int number_of_electrons = 0;
            std::vector<std::string> results;

            std::ifstream input_file(log_file_name_full.c_str());

            if (input_file.fail()) {
                CTP_LOG(ctp::logERROR, *_pLog) << "File " << log_file_name_full << " not found " << flush;
                return false;
            } else {
                CTP_LOG(ctp::logDEBUG, *_pLog) << "Reading Coordinates and occupationnumbers and energies from " << log_file_name_full << flush;
            }
            //Coordinates of the final configuration depending on whether it is an optimization or not
            while (input_file) {
                getline(input_file, line);
                boost::trim(line);

                if (_is_optimization) {
                    throw runtime_error("Not implemented yet!");
                }
                bool found_optimization = true;
                std::string::size_type coordinates_pos = line.find("CARTESIAN COORDINATES (ANGSTROEM)");

                if (found_optimization && coordinates_pos != std::string::npos) {
                    CTP_LOG(ctp::logDEBUG, *_pLog) << "Getting the coordinates" << flush;
                    bool has_QMAtoms = orbitals.hasQMAtoms();
                    // three garbage lines
                    getline(input_file, line);
                    // now starts the data in format
                    // _id type Qnuc x y z
                    std::vector<std::string> row=GetLineAndSplit(input_file, "\t ");
                    int nfields = row.size();
                    int atom_id = 0;
                    while (nfields == 4) {
                        std::string atom_type = row.at(0);
                        double x = boost::lexical_cast<double>(row.at(1));
                        double y = boost::lexical_cast<double>(row.at(2));
                        double z = boost::lexical_cast<double>(row.at(3));
                        row=GetLineAndSplit(input_file, "\t ");
                        nfields = row.size();
                        tools::vec pos=tools::vec(x,y,z);
                        pos*=tools::conv::ang2bohr;
                        if (has_QMAtoms == false) {
                            orbitals.AddAtom(atom_id,atom_type, pos);
                        } else {
                            QMAtom* pAtom = orbitals.QMAtoms().at(atom_id);
                            pAtom->setPos(pos);
                        }
                        atom_id++;
                    }
                }

                std::string::size_type energy_pos = line.find("FINAL SINGLE");
                if (energy_pos != std::string::npos) {
                    boost::algorithm::split(results, line, boost::is_any_of(" "), boost::algorithm::token_compress_on);
                    std::string energy = results[4];
                    boost::trim(energy);
                    orbitals.setQMEnergy(boost::lexical_cast<double>(energy));
                    CTP_LOG(ctp::logDEBUG, *_pLog) << (boost::format("QM energy[Hrt]: %4.8f ") % orbitals.getQMEnergy()).str() << flush;
                }

                std::string::size_type HFX_pos = line.find("Fraction HF Exchange ScalHFX");
                if (HFX_pos != std::string::npos) {
                    boost::algorithm::split(results, line, boost::is_any_of(" "), boost::algorithm::token_compress_on);
                    double ScaHFX = boost::lexical_cast<double>(results.back());
                    orbitals.setScaHFX(ScaHFX);
                    CTP_LOG(ctp::logDEBUG, *_pLog) << "DFT with " << ScaHFX << " of HF exchange!" << flush;
                }

                std::string::size_type dim_pos = line.find("Basis Dimension");
                if (dim_pos != std::string::npos) {
                    boost::algorithm::split(results, line, boost::is_any_of(" "), boost::algorithm::token_compress_on);
                    std::string dim = results[4]; //The 4th element of results vector is the Basis Dim
                    boost::trim(dim);
                    levels = boost::lexical_cast<int>(dim);
                    CTP_LOG(ctp::logDEBUG, *_pLog) << "Basis Dimension: " << levels << flush;
                    CTP_LOG(ctp::logDEBUG, *_pLog) << "Energy levels: " << levels << flush;
                }

                std::string::size_type OE_pos = line.find("ORBITAL ENERGIES");
                if (OE_pos != std::string::npos) {
                    number_of_electrons=0;
                    getline(input_file, line);
                    getline(input_file, line);
                    getline(input_file, line);
                    if (line.find("E(Eh)") == std::string::npos) {
                        CTP_LOG(ctp::logDEBUG, *_pLog) << "Warning: Orbital Energies not found in log file" << flush;
                    }
                    for (unsigned i = 0; i < levels; i++) {
                        results=GetLineAndSplit(input_file, " ");
                        std::string no = results[0];
                        boost::trim(no);
                        unsigned levelnumber = boost::lexical_cast<unsigned>(no);
                        if (levelnumber != i) {
                            CTP_LOG(ctp::logDEBUG, *_pLog) << "Have a look at the orbital energies something weird is going on" << flush;
                        }
                        std::string oc = results[1];
                        boost::trim(oc);
                        double occ = boost::lexical_cast<double>(oc);
                        // We only count alpha electrons, each orbital must be empty or doubly occupied
                        if (occ == 2 || occ == 1) {
                            number_of_electrons++;
                            occupancy[i] = occ;
                        } else if (occ == 0) {
                            occupancy[i] = occ;
                        } else {
                            if (occ == 1){
                                CTP_LOG(ctp::logDEBUG, *_pLog) << "Watch out! No distinction between alpha and beta electrons. Check if occ = 1 is suitable for your calculation " << flush;
                                number_of_electrons++;
                                occupancy[i] = occ;
                            } else {
                            throw runtime_error("Only empty or doubly occupied orbitals are allowed not running the right kind of DFT calculation");
                            }
                        }
                        std::string e = results[2];
                        boost::trim(e);
                        energies [i] = boost::lexical_cast<double>(e);
                    }
                }
                /*
                 *  Partial charges from the input file
                 */
                std::string::size_type charge_pos = line.find("CHELPG Charges");

                if (charge_pos != std::string::npos && _get_charges) {
                    CTP_LOG(ctp::logDEBUG, *_pLog) << "Getting charges" << flush;
                    getline(input_file, line);
                    std::vector<std::string> row=GetLineAndSplit(input_file, "\t ");
                    int nfields = row.size();
                    while (nfields == 4) {
                        int atom_id = boost::lexical_cast< int >(row.at(0));
                        atom_id++;
                        std::string atom_type = row.at(1);
                        double atom_charge = boost::lexical_cast< double >(row.at(3));
                        row=GetLineAndSplit(input_file, "\t ");
                        nfields = row.size();
                        QMAtom* pAtom;
                        if (!orbitals.hasQMAtoms()) {
                            pAtom =orbitals.AddAtom(atom_id - 1,atom_type, tools::vec(0.0));
                        } else {
                            pAtom = orbitals.QMAtoms().at(atom_id - 1);
                        }
                        pAtom->setPartialcharge(atom_charge);
                    }
                }


                std::string::size_type success = line.find("*                     SUCCESS                       *");
                if (success != std::string::npos ) {
                    found_success=true;
                }
            }

            CTP_LOG(ctp::logDEBUG, *_pLog) << "Alpha electrons: " << number_of_electrons << flush;
            int occupied_levels = number_of_electrons;
            int unoccupied_levels = levels - occupied_levels;
            CTP_LOG(ctp::logDEBUG, *_pLog) << "Occupied levels: " << occupied_levels << flush;
            CTP_LOG(ctp::logDEBUG, *_pLog) << "Unoccupied levels: " << unoccupied_levels << flush;

            /************************************************************/

            // copying information to the orbitals object
       
            orbitals.setBasisSetSize(levels);
            orbitals.setNumberOfAlphaElectrons(number_of_electrons);
            orbitals.setNumberOfOccupiedLevels(occupied_levels);
            orbitals.setSelfEnergy(0.0);

            // copying energies to a vector
            orbitals.MOEnergies().resize(levels);
            //_level = 1;
            for (int i = 0; i < orbitals.MOEnergies().size(); i++) {
                orbitals.MOEnergies()[i] = energies[ i ];
            }



            CTP_LOG(ctp::logDEBUG, *_pLog) << "Done reading Log file" << flush;
            return found_success;
        }

        bool Orca::CheckLogFile() {
            // check if the log file exists
            ifstream input_file((_run_dir + "/" + _log_file_name).c_str());

            if (input_file.fail()) {
                CTP_LOG(ctp::logERROR, *_pLog) << "Orca LOG is not found" << flush;
                return false;
            };

            std::string line;
            while (input_file) {
                getline(input_file, line);
                boost::trim(line);
                std::string::size_type error = line.find("FATAL ERROR ENCOUNTERED");
                if (error != std::string::npos) {
                    CTP_LOG(ctp::logERROR, *_pLog) << "ORCA encountered a fatal error, maybe a look in the log file may help." << flush;
                    return false;
                }
                error = line.find("mpirun detected that one or more processes exited with non-zero status");
                if (error != std::string::npos) {
                    CTP_LOG(ctp::logERROR, *_pLog) << "ORCA had an mpi problem, maybe your openmpi version is not good." << flush;
                    return false;
                }
            }
            return true;
        }

        // Parses the Orca gbw file and stores data in the Orbitals object

        bool Orca::ParseOrbitalsFile(Orbitals& orbitals) {
            if (!CheckLogFile()) return false;
            std::vector<double> coefficients;
            int basis_size = orbitals.getBasisSetSize();
            int levels = orbitals.getBasisSetSize();
            if (basis_size == 0 || levels == 0) {
                throw runtime_error("Basis size not set, calculator does not parse log file first");
            }
            CTP_LOG(ctp::logDEBUG, *_pLog) << "Reading the gbw file, this may or may not work so be careful: " << flush;
            ifstream infile;
            infile.open((_run_dir + "/" + _orb_file_name).c_str(), ios::binary | ios::in);
            if (!infile) {
                throw runtime_error("Could not open " + _orb_file_name + " file");
            }
            infile.seekg(24, ios::beg);
            char* buffer = new char [8];
            infile.read(buffer, 8);
            long int offset = *((long int*) buffer);

            infile.seekg(offset, ios::beg);
            infile.read(buffer, 4);
            int op_read = *((int*) buffer);
            infile.seekg(offset + 4, ios::beg);
            infile.read(buffer, 4);
            int dim_read = *((int*) buffer);
            infile.seekg(offset + 8, ios::beg);
            CTP_LOG(ctp::logDEBUG, *_pLog) << "Number of operators: " << op_read << " Basis dimension: " << dim_read << flush;
            int n = op_read * dim_read*dim_read;
            delete[] buffer;
            buffer = new char [8];
            for (int i = 0; i < n; i++) {
                infile.read(buffer, 8);
                double mocoeff = *((double*) buffer);
                //CTP_LOG(logDEBUG,*_pLog) << mocoeff<< flush ;
                coefficients.push_back(mocoeff);
            }
            delete[] buffer;

            infile.close();
            // i -> MO, j -> AO
            (orbitals.MOCoefficients()).resize(levels, basis_size);
            for (int i = 0; i < orbitals.MOCoefficients().rows(); i++) {
                for (int j = 0; j < orbitals.MOCoefficients().cols(); j++) {
                    orbitals.MOCoefficients()(j, i) = coefficients[j * basis_size + i];
                   
                }
            }  
           ReorderOutput(orbitals);
            CTP_LOG(ctp::logDEBUG, *_pLog) << "Done parsing" << flush;
            return true;
        }

        std::string Orca::getLName(int lnum) {
            if (lnum == 0) {
                return "S";
            } else if (lnum == 1) {
                return "P";
            } else if (lnum == 2) {
                return "D";
            } else if (lnum == 3) {
                return "F";
            } else {
                throw runtime_error("Orca::getLName functions higher than F not implemented");
            }
            return "0";
        }

        std::string Orca::indent(const double &number) {
            std::stringstream ssnumber;
            if (number >= 0) {
                ssnumber << "    ";
            } else {
                ssnumber << "   ";
            }
            ssnumber << setiosflags(ios::fixed) << setprecision(15) << std::scientific << number;
            std::string snumber = ssnumber.str();
            return snumber;
        }




    }
}
