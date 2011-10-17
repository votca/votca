/*
 * Copyright 2009-2011 The VOTCA Development Team (http://www.votca.org)
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

#include <votca/csg/trajectoryreader.h>
#include <votca/csg/trajectorywriter.h>
#include <votca/csg/topologyreader.h>
#include <votca/ctp/calculatorfactory.h>
#include <votca/ctp/qmapplication.h>
#include <votca/ctp/version.h>

namespace votca { namespace ctp {

    QMApplication::QMApplication()
{
    CalculatorFactory::RegisterAll();
}

QMApplication::~QMApplication()
{}

void QMApplication::Initialize(void)
{
    // register all io plugins
    TrajectoryWriter::RegisterPlugins();
    TrajectoryReader::RegisterPlugins();
    TopologyReader::RegisterPlugins();
    CalculatorFactory::RegisterAll();

    AddProgramOptions()
        ("segments,s", boost::program_options::value<string>(), "  description of conjugated segments")
        ("options,o", boost::program_options::value<string>(), "  program and calculator options")
        ("file,f", boost::program_options::value<string>(), "  sqlite state file")
        ("first-frame,i", boost::program_options::value<int>()->default_value(1), "  start frame (first frame is 1)")
        ("nframes,n", boost::program_options::value<int>()->default_value(-1), "  process so many frames")
        //  this is shit, move it out!
        //("nnnames", boost::program_options::value<string>()->default_value("*"), "  List of strings that the concatenation of the two molnames must match to be analyzed")
        ;
}

bool QMApplication::EvaluateOptions(void)
{
    CheckRequired("segments", "please provide an xml file with the description of segments.");
    CheckRequired("options", "please provide an xml file with the program options");
    CheckRequired("file", "no database file specified");
    return true;
}

void QMApplication::Run()
{
    _qmtop.LoadListCharges(_op_vm["segments"].as<string>());
    // read in program options from options.xml
    load_property_from_xml(_options, _op_vm["options"].as<string>());

    int first_frame = OptionsMap()["first-frame"].as<int>(); /// starting frame
    if(first_frame == 0) throw std::runtime_error("error, first frame is 0 but we start counting with 1");
    first_frame--;

    int nframes = OptionsMap()["nframes"].as<int>(); /// number of frames to be processed

    BeginEvaluate();

    /// load qmtop from state saver
    cout << "Loading qmtopology via state saver." << endl;
    string statefile = OptionsMap()["file"].as<string>();
    StateSaverSQLite loader;
    loader.Open(_qmtop, statefile);
    if(loader.FramesInDatabase() != 1)
        throw std::runtime_error("database contains none or more than one frame which is not supported yet.");

    while(loader.NextFrame()) {
        EvaluateFrame();
        loader.WriteFrame();
    }
    loader.Close();
    EndEvaluate();
}

void QMApplication::ShowHelpText(std::ostream &out)
{
    string name =  ProgramName();
    if(VersionString() != "")
         name = name + ", version " + VersionString();

    votca::ctp::HelpTextHeader(name);
    HelpText(out);
    out << "\n\n" << OptionsDesc() << endl;
}

/*void QMApplication::PrintNbs(string filename){
    ofstream out_nbl;
    out_nbl.open(filename.c_str());
    QMNBList &nblist = _qmtop.nblist();
    if(out_nbl!=0){
        out_nbl << "Neighbours, J(0), J_eff, rate, r_ij, abs(r_ij) [nm]" << endl;
        QMNBList::iterator iter;
        for ( iter  = nblist.begin(); iter != nblist.end() ; ++iter){
            ///Hack!
            if ((*iter)->Js().size() > 0 ){
                out_nbl << "(" << (*iter)->first->getId() << "," << (*iter)->second->getId() << "): ";
                out_nbl << (*iter)->Js()[0] << " " << sqrt((*iter)->calcJeff2()) << " " << (*iter)->rate12() << " ";
                out_nbl << (*iter)->r().getX() << " " << (*iter)->r().getY() << " " << (*iter)->r().getZ() << " ";
                out_nbl << " " << (*iter)->dist() << endl;
            }
        }
    }
    out_nbl.close();
}*/

void QMApplication::AddCalculator(QMCalculator* calculator){
    _calculators.push_back(calculator);
}

void QMApplication::BeginEvaluate(){
    list<QMCalculator *>::iterator iter;
    for (iter = _calculators.begin(); iter != _calculators.end(); ++iter){
        (*iter)->Initialize(&_qmtop, &_options);
    }
}

bool QMApplication::EvaluateFrame(){
    list<QMCalculator *>::iterator iter;
    for (iter = _calculators.begin(); iter != _calculators.end(); ++iter){
        (*iter)->EvaluateFrame(&_qmtop);
    }
}

void QMApplication::EndEvaluate(){
    list<QMCalculator *>::iterator iter;
    for (iter = _calculators.begin(); iter != _calculators.end(); ++iter){
        (*iter)->EndEvaluate(&_qmtop);
    }
}
}}
