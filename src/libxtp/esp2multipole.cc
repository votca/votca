/*
 *            Copyright 2009-2016 The VOTCA Development Team
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


#include <votca/xtp/esp2multipole.h>
#include <boost/format.hpp>

namespace votca { namespace xtp {

void Esp2multipole::Initialize(Property* options) {
    string key = Identify();
    _use_ecp=false;
    _do_svd=false;
    
    _use_mulliken=false;
    _use_CHELPG=false;
    _use_GDMA=false;
    _use_CHELPG_SVD=false;
         
    _state     = options->get(key + ".state").as<string> (); 
    _state_no     = options->get(key + ".statenumber").as<int> ();
    _spin     = options->get(key + ".spin").as<string> ();
    if ( options->exists(key+".ecp")) {
       _use_ecp=options->get(key + ".ecp").as<bool> ();
    }
    if ( options->exists(key+".method")) {
         _method = options->get(key + ".method").as<string> ();
         if (_method=="Mulliken")_use_mulliken=true; 
         else if(_method=="mulliken")_use_mulliken=true; 
         else if(_method=="CHELPG")_use_CHELPG=true; 
         else if(_method=="GDMA") throw std::runtime_error("GDMA not implemented yet");
         else if(_method=="CHELPG_SVD") throw std::runtime_error("CHELPG_SVD not implemented yet"); 
         else  throw std::runtime_error("Method not recognized. Only Mulliken and CHELPG implemented");
         }
    else _use_CHELPG=true;
    if (!_use_mulliken){
         _integrationmethod     = options->get(key + ".integrationmethod").as<string> ();
    }
    

    if (!(_integrationmethod=="numeric" || _integrationmethod=="analytic")){
        std::runtime_error("Method not recognized. Only numeric and analytic available");
    }
    if ( options->exists(key+".gridsize")) {
         _gridsize = options->get(key+".gridsize").as<string>();
         }
    else _gridsize="medium";
    if ( options->exists(key+".openmp")) {
         _openmp_threads = options->get(key+".openmp").as<int>();
         }
    else _openmp_threads=0;
    if ( options->exists(key+".svd")) {
         _do_svd = options->get(key+".svd.do_svd").as<bool>();
         _conditionnumber = options->get(key+".svd.conditionnumber").as<double>();
         }
    
              
    
  
    
    // get the path to the shared folders with xml files
    char *votca_share = getenv("VOTCASHARE");    
    if(votca_share == NULL) throw std::runtime_error("VOTCASHARE not set, cannot open help files.");

  
   
}
string  Esp2multipole::GetIdentifier(){
    string identifier;
    if(_state=="transition" && _spin=="singlet" ){
        identifier=(boost::format("n2s%i") % _state_no).str();
    } 
    else if(_state=="transition" && _spin=="triplet" ){
        identifier=(boost::format("n2t%i") % _state_no).str();
    } 
    else if(_state=="excited" && _spin=="triplet" ){
        identifier=(boost::format("t%i") % _state_no).str();
    } 
    else if(_state=="excited" && _spin=="singlet" ){
        identifier=(boost::format("s%i") % _state_no).str();
    } 
    else if(_state=="excited" && _spin=="singlet" ){
        identifier=(boost::format("s%i") % _state_no).str();
    } 
    else if(_state=="ground"  ){
        identifier="n";
    } 
    else{
        throw std::runtime_error("Esp2multipole::GetIdentifier did not recognize config");
    }
    return identifier;
}

void Esp2multipole::WritetoFile(string _output_file, string identifier){
    _use_mps=false;
    _use_pdb=false;
    string data_format  = boost::filesystem::extension(_output_file);  
    //cout << data_format << endl;
    if (data_format==".mps")_use_mps=true; 
    else if(data_format==".pdb")_use_pdb=true;    
    else  throw std::runtime_error("Outputfile format not recognized. Export only to .pdb and .mps");
     string tag="TOOL:"+Identify()+"_"+_state+"_"+_spin+boost::lexical_cast<string>(_state_no);
    if(_use_mps){
        QMMInterface Converter;
        PolarSeg* result=Converter.Convert(_Atomlist);
        
        result->WriteMPS(_output_file,tag);
        }
    else if(_use_pdb){
        FILE *out;
        Orbitals _orbitals;
        std::vector< QMAtom* >::iterator at;
         for (at=_Atomlist.begin();at<_Atomlist.end();++at){
            
            _orbitals.AddAtom(*(*at));
        }
        out = fopen(_output_file.c_str(), "w");
       _orbitals.WritePDB(out, tag); 
    }
}




void Esp2multipole::Extractingcharges( Orbitals& _orbitals ){
    int threads=1;
#ifdef _OPENMP
            if ( _openmp_threads > 0 ) omp_set_num_threads(_openmp_threads); 
            threads=omp_get_max_threads();
#endif
   LOG(logDEBUG, *_log) << "===== Running on "<< threads << " threads ===== " << flush;

        vector< QMAtom* > Atomlist =_orbitals.QMAtoms();
        std::vector< QMAtom* >::iterator at;
        for (at=Atomlist.begin();at<Atomlist.end();++at){
            QMAtom * atom=new QMAtom(*(*at));
            _Atomlist.push_back(atom);
        }
        ub::matrix<double> DMAT_tot;
        BasisSet bs;
        bs.LoadBasisSet(_orbitals.getDFTbasis());
        AOBasis basis;
        basis.AOBasisFill(&bs, _Atomlist );
        
        
        ub::matrix<double> _MO_Coefficients = *(_orbitals.getOrbitals()); // this is a copy?
        
        //basis.ReorderMOs(_orbitals.MOCoefficients(), _orbitals.getQMpackage(), "votca" );  
        basis.ReorderMOs(_MO_Coefficients, _orbitals.getQMpackage(), "votca" );  
        bool _do_transition=false;
        if(_state=="transition"){
            _do_transition=true;
            if (_spin=="singlet"){
                //DMAT_tot=_orbitals.TransitionDensityMatrix(_orbitals.MOCoefficients() , _orbitals.BSESingletCoefficients(), _state_no-1);
                DMAT_tot=_orbitals.TransitionDensityMatrix(_MO_Coefficients, _orbitals.BSESingletCoefficients(), _state_no-1);
            }
            else if (_spin=="triplet"){
                //DMAT_tot=_orbitals.TransitionDensityMatrix(_orbitals.MOCoefficients() , _orbitals.BSETripletCoefficients(), _state_no-1); 
                DMAT_tot=_orbitals.TransitionDensityMatrix(_MO_Coefficients, _orbitals.BSETripletCoefficients(), _state_no-1); 
            }
            else throw std::runtime_error("Spin entry not recognized");
        }
        else if (_state=="ground" || _state=="excited"){
            
        
            //ub::matrix<double> &DMATGS=_orbitals.DensityMatrixGroundState(_orbitals.MOCoefficients());
            ub::matrix<double> &DMATGS=_orbitals.DensityMatrixGroundState(_MO_Coefficients);
            DMAT_tot=DMATGS;
            if ( _state_no > 0 && _state=="excited"){
                std::vector<ub::matrix<double> > DMAT;
                if (_spin=="singlet"){
                    //DMAT = _orbitals.DensityMatrixExcitedState( _orbitals.MOCoefficients() , _orbitals.BSESingletCoefficients(), _state_no-1);
                    DMAT = _orbitals.DensityMatrixExcitedState( _MO_Coefficients , _orbitals.BSESingletCoefficients(), _state_no-1);
                }
                else if (_spin=="triplet"){
                    //DMAT = _orbitals.DensityMatrixExcitedState( _orbitals.MOCoefficients() , _orbitals.BSETripletCoefficients(), _state_no-1);
                    DMAT = _orbitals.DensityMatrixExcitedState( _MO_Coefficients , _orbitals.BSETripletCoefficients(), _state_no-1);
                }
                else throw std::runtime_error("Spin entry not recognized");
                DMAT_tot=DMAT_tot-DMAT[0]+DMAT[1];
            }            
	   // Ground state + hole_contribution + electron contribution
	}
        else throw std::runtime_error("State entry not recognized");
        
        
        if (_use_mulliken) {
            Mulliken mulliken;
            mulliken.setUseECPs(_use_ecp);
            mulliken.EvaluateMulliken(_Atomlist, DMAT_tot, basis, bs, _do_transition);
                
        }
        else if (_use_CHELPG){         
            Espfit esp=Espfit(_log);
            esp.setUseECPs(_use_ecp);
            if(_do_svd){
                esp.setUseSVD(_do_svd,_conditionnumber);
            }
            if (_integrationmethod=="numeric")  {
                esp.Fit2Density(_Atomlist, DMAT_tot, basis,bs,_gridsize); 
            }
            else if (_integrationmethod=="analytic")  esp.Fit2Density_analytic(_Atomlist,DMAT_tot,basis);
        }
}       

}}
