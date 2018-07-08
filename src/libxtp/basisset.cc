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
#include "votca/xtp/basisset.h"
#include <votca/tools/property.h>


namespace votca { namespace xtp {
    

       
int FindLmax(const std::string& _type ){

        int _lmax;
        // single type shells
        if ( _type.length() == 1 ){
            if ( _type == "S" ){ _lmax = 0;}
            else if ( _type == "P" ){ _lmax = 1;}
            else if ( _type == "D" ){ _lmax = 2;}
            else if ( _type == "F" ){ _lmax = 3;}
            else if ( _type == "G" ){ _lmax = 4;}
            else if ( _type == "H" ){ _lmax = 5;}
            else if ( _type == "I" ){ _lmax = 6;}
            else{
                throw std::runtime_error("FindLmax: Shelltype not known");
            }
        } else {
            _lmax = -1;
            for(unsigned i = 0; i < _type.length(); ++i) {
                std::string local_shell = std::string( _type, i, 1 );
                int _test = FindLmax( local_shell  );
                if ( _test > _lmax ) { _lmax = _test;} 
            }   
        }

        return _lmax;
        
    }
 
 int FindLmin(const std::string& _type ){

        int _lmin;
        // single type shells
        if ( _type.length() == 1 ){
            if ( _type == "S" ){ _lmin = 0;}
            else if ( _type == "P" ){ _lmin = 1;}
            else if ( _type == "D" ){ _lmin = 2;}
            else if ( _type == "F" ){ _lmin = 3;}
            else if ( _type == "G" ){ _lmin = 4;}
            else if ( _type == "H" ){ _lmin = 5;}
            else if ( _type == "I" ){ _lmin = 6;}
            else{
                throw std::runtime_error("FindLmax: Shelltype not known");
            }
        } else {
            _lmin = 10;
            for(unsigned i = 0; i < _type.length(); ++i) {
                std::string local_shell = std::string( _type, i, 1 );
                int _test = FindLmin( local_shell  );
                if(_test==0){return 0;}
                if ( _test < _lmin ) { _lmin = _test;} 
            }   
        }

        return _lmin;
        
    }
 
  int OffsetFuncShell(const std::string& shell_type ) {
    int _nbf;
    // single type shells
    if ( shell_type.length() == 1 ){
       if ( shell_type == "S" ){ _nbf = 0;}
       else if ( shell_type == "P" ){ _nbf = 1;}
       else if ( shell_type == "D" ){ _nbf = 4;}
       else if ( shell_type == "F" ){ _nbf = 9;}
       else if ( shell_type == "G" ){ _nbf = 16;}
       else if ( shell_type == "H" ){ _nbf = 25;}
       else if ( shell_type == "I" ){ _nbf = 36;}
       else{
                throw std::runtime_error("OffsetFuncShell: Shelltype not known");
            }
    } else {
        // for combined shells, go over all contributions and find minimal offset
        _nbf = 1000;
        for(unsigned i = 0; i < shell_type.length(); ++i) {
            std::string local_shell = std::string( shell_type, i, 1 );
            int _test = OffsetFuncShell( local_shell  );
            if ( _test < _nbf ) { _nbf = _test;} 
        }   
    }
    return _nbf;
        } 
  
  
    int NumFuncShell( const std::string& shell_type) {
        int _nbf = 0;
        // single type shells
        if ( shell_type.length() == 1 ){
            if ( shell_type == "S" ){ _nbf = 1;}
            else if ( shell_type == "P" ){ _nbf = 3;}
            else if ( shell_type == "D" ){ _nbf = 5;}
            else if ( shell_type == "F" ){ _nbf = 7;}
            else if ( shell_type == "G" ){ _nbf = 9;}
            else if ( shell_type == "H" ){ _nbf = 11;}
            else if ( shell_type == "I" ){ _nbf = 13;}
            else{
                throw std::runtime_error("FindnumofFunc: Shelltype not known");
            }
    } else {
        // for combined shells, go over all contributions and add functions
        _nbf = 0;
        
           for (unsigned i = 0; i < shell_type.length(); ++i) {
                    std::string local_shell = std::string(shell_type, i, 1);
                    _nbf += NumFuncShell(local_shell);
            }   
    }
    return _nbf;
        
    }
    
    
      std::vector<int> NumFuncSubShell(const std::string& shell_type) {
          std::vector <int> subshells;
        // single type shells
        if ( shell_type.length() == 1 ){
           subshells.push_back( NumFuncShell(shell_type));
        // for combined shells, go over all contributions and add functions
        }else{
           for (unsigned i = 0; i < shell_type.length(); ++i) {
                    std::string local_shell = std::string(shell_type, i, 1);
                   subshells.push_back( NumFuncShell(local_shell));
            }   
    }
    return subshells;        
    }
    
    int NumFuncShell_cartesian(const std::string& shell_type ) {
    int _nbf;
    // single type shells defined here
    if ( shell_type.length() == 1 ){
       if ( shell_type == "S" ){ _nbf = 1;}
       else if ( shell_type == "P" ){ _nbf = 3;}
       else if ( shell_type == "D" ){ _nbf = 6;}
       else if ( shell_type == "F" ){ _nbf = 10;}
       else if ( shell_type == "G" ){ _nbf = 15;}
       else if ( shell_type == "H" ){ _nbf = 21;}
       else if ( shell_type == "I" ){ _nbf = 28;}
       else {
                    throw std::runtime_error("NumFuncShell_cartesian shell_type not known");
                }
    } else {
        // for combined shells, sum over all contributions
        _nbf = 0;
        for( unsigned i = 0; i < shell_type.length(); ++i) {
            std::string local_shell =    std::string( shell_type, i, 1 );
            _nbf += NumFuncShell_cartesian( local_shell  );
        }
    }

    return _nbf;
}
    
    
    
       

int OffsetFuncShell_cartesian(const std::string& shell_type ) {
    int _nbf;
    // single type shells
    if ( shell_type.length() == 1 ){
       if ( shell_type == "S" ){ _nbf = 0;}
       else if ( shell_type == "P" ){ _nbf = 1;}
       else if ( shell_type == "D" ){ _nbf = 4;}
       else if ( shell_type == "F" ){ _nbf = 10;}
       else if ( shell_type == "G" ){ _nbf = 20;}
       else if ( shell_type == "H" ){ _nbf = 35;}
       else if ( shell_type == "I" ){ _nbf = 56;}
       else {
                    throw std::runtime_error("OffsetFuncShell_cartesian shell_type not known");
                }
    } else {
        // for combined shells, go over all contributions and find minimal offset
        _nbf = 1000;
        for(unsigned i = 0; i < shell_type.length(); ++i) {
            std::string local_shell = std::string( shell_type, i, 1 );
            int _test = OffsetFuncShell_cartesian( local_shell  );
            if ( _test < _nbf ) { _nbf = _test;} 
        }   
    }
    return _nbf;
}    
    

void BasisSet::LoadBasisSet ( std::string name ) 
{    
    tools::Property basis_property;
    _name=name;
    // if name contains .xml, assume a basisset .xml file is located in the working directory
    std::size_t found_xml = name.find(".xml");
    std::string xmlFile;
    if (found_xml!=std::string::npos) {
      xmlFile = name;
    } else {
      // get the path to the shared folders with xml files
      char *votca_share = getenv("VOTCASHARE");
      if(votca_share == NULL) throw std::runtime_error("VOTCASHARE not set, cannot open help files.");
      xmlFile = std::string(getenv("VOTCASHARE")) + std::string("/xtp/basis_sets/") + name + std::string(".xml");
    }
    bool success = load_property_from_xml(basis_property, xmlFile);
    
    if ( !success ) {; }
    
    std::list<tools::Property*> elementProps = basis_property.Select("basis.element");
        
    for (std::list<tools::Property*> ::iterator  ite = elementProps.begin(); ite != elementProps.end(); ++ite) 
    {       
        std::string elementName = (*ite)->getAttribute<std::string>("name");
        Element *element = addElement( elementName );
        //cout << "\nElement " << elementName;
        
        std::list<tools::Property*> shellProps = (*ite)->Select("shell");
        for (std::list<tools::Property*> ::iterator  its = shellProps.begin(); its != shellProps.end(); ++its) 
        {            
            std::string shellType = (*its)->getAttribute<std::string>("type");
            double shellScale = (*its)->getAttribute<double>("scale");
            
            Shell* shell = element->addShell( shellType, shellScale );
            //cout << "\n\tShell " << shellType;
            
            std::list<tools::Property*> constProps = (*its)->Select("constant");
            for (std::list<tools::Property*> ::iterator  itc = constProps.begin(); itc != constProps.end(); ++itc) {
                double decay = (*itc)->getAttribute<double>("decay");
                // << " decay "<<decay<<endl;
                std::vector<double> contraction=std::vector<double>(shell->getLmax()+1,0.0);
                std::list<tools::Property*> contrProps = (*itc)->Select("contractions");
                for (std::list<tools::Property*> ::iterator itcont = contrProps.begin(); itcont != contrProps.end(); ++itcont){
                    std::string contrType = (*itcont)->getAttribute<std::string>("type");
                    double contrFactor = (*itcont)->getAttribute<double>("factor");
                    //cout << " factor " << contrFactor << endl;
                    if ( contrType == "S" ) contraction[0] = contrFactor;
                    else if ( contrType == "P" ) contraction[1] = contrFactor;
                    else if ( contrType == "D" ) contraction[2] = contrFactor;
                    else if ( contrType == "F" ) contraction[3] = contrFactor;
                    else if ( contrType == "G" ) contraction[4] = contrFactor;
                    else if ( contrType == "H" ) contraction[5] = contrFactor;
                    else if ( contrType == "I" ) contraction[6] = contrFactor;
                    else{
                        throw std::runtime_error("LoadBasiset:Contractiontype not known");
                    }
                }    
                 
                shell->addGaussian(decay, contraction);
            }
            
        }
       
    }
 return;  
}


void BasisSet::LoadPseudopotentialSet ( std::string name ) 
{    
    tools::Property basis_property;
    _name=name;
  // if name contains .xml, assume a ecp .xml file is located in the working directory
    std::size_t found_xml = name.find(".xml");
    std::string xmlFile;
    if (found_xml!=std::string::npos) {
      xmlFile = name;
    } else {
      // get the path to the shared folders with xml files
      char *votca_share = getenv("VOTCASHARE");
      if(votca_share == NULL) throw std::runtime_error("VOTCASHARE not set, cannot open help files.");
      xmlFile = std::string(getenv("VOTCASHARE")) + std::string("/xtp/ecps/") + name + std::string(".xml");
    }
    bool success = load_property_from_xml(basis_property, xmlFile);
    
    if ( !success ) {; }
    
    std::list<tools::Property*> elementProps = basis_property.Select("pseudopotential.element");
        
    for (std::list<tools::Property*> ::iterator  ite = elementProps.begin(); ite != elementProps.end(); ++ite) 
    {       
        std::string elementName = (*ite)->getAttribute<std::string>("name");
        int lmax = (*ite)->getAttribute<int>("lmax");
        int ncore = (*ite)->getAttribute<int>("ncore");
        
        Element *element = addElement( elementName, lmax, ncore );
        //cout << "\nElement " << elementName;
        
        std::list<tools::Property*> shellProps = (*ite)->Select("shell");
        for (std::list<tools::Property*> ::iterator  its = shellProps.begin(); its != shellProps.end(); ++its) 
        {            
            std::string shellType = (*its)->getAttribute<std::string>("type");
            double shellScale = 1.0;
            
            Shell* shell = element->addShell( shellType, shellScale );
            //cout << "\n\tShell " << shellType;
            
            std::list<tools::Property*> constProps = (*its)->Select("constant");
            for (std::list<tools::Property*> ::iterator  itc = constProps.begin(); itc != constProps.end(); ++itc) 
            {
                int power = (*itc)->getAttribute<int>("power");
                double decay = (*itc)->getAttribute<double>("decay");
                //double contraction = (*itc)->getAttribute<double>("contraction");
                std::vector<double> contraction;
                // just testing here with single value 
                contraction.push_back((*itc)->getAttribute<double>("contraction"));
                //shell->addGaussian(decay, contraction);
                shell->addGaussian(power, decay, contraction);
                //cout << "\n\t\t" << decay << " " << contraction << endl;
            }
            
        }
       
    }
    return;
}

// adding an Element to a Basis Set
Element* BasisSet::addElement( std::string elementType ) {
    Element *element = new Element( elementType );
    _elements[elementType] = element;
    return element;
};

// adding an Element to a Pseudopotential Library
Element* BasisSet::addElement( std::string elementType, int lmax, int ncore ) {
    Element *element = new Element( elementType, lmax, ncore );
    _elements[elementType] = element;
    return element;
};

// cleanup the basis set
BasisSet::~BasisSet() {
    
    for ( std::map< std::string,Element* >::iterator it = _elements.begin(); it !=  _elements.end(); it++ ) {
         delete (*it).second;
     }
    
    _elements.clear();
};
        



}}
