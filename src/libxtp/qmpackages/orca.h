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

#ifndef __VOTCA_XTP_ORCA_H
#define	__VOTCA_XTP_ORCA_H


#include <votca/ctp/apolarsite.h>
#include <votca/xtp/qmpackage.h>

#include <string>



namespace votca { namespace xtp {
/**
    \brief Wrapper for the ORCA program

    The ORCA class executes the ORCA package
    and extracts information from its log and io files

*/
class Orca : public QMPackage
{
public:

   std::string getPackageName() { return "orca"; }

   void Initialize( tools::Property &options );

   bool WriteInputFile( Orbitals& orbitals);

   bool WriteShellScript();

   bool Run();

   void CleanUp();

   bool CheckLogFile();

   bool ParseLogFile( Orbitals& orbitals );

   bool ParseOrbitalsFile( Orbitals& orbitals );


   std::string getScratchDir( ) { return _scratch_dir; }

private:

    std::string                              _shell_file_name;
    std::string                              _scratch_dir;
    bool                                _is_optimization;

    std::string                              _cleanup;

    std::string indent( const double &number );
    std::string getLName(int lnum);

    void WriteBasisset(std::vector<QMAtom*>& qmatoms, std::string& _bs_name, std::string& _el_file_name);
    void WriteCoordinates(std::ofstream& _com_file, std::vector<QMAtom*>& qmatoms);
    void WriteECP(std::ofstream& _com_file, std::vector<QMAtom*>& qmatoms);
    void WriteBackgroundCharges();
    
    void WriteChargeOption();
};


}}

#endif	/* __VOTCA_XTP_ORCA_H */
