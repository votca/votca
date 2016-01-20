#ifndef VOTCA_XTP_PTOPREADER_H
#define VOTCA_XTP_PTOPREADER_H

#include <votca/xtp/qmtool.h>
#include <votca/xtp/polartop.h>
#include <votca/xtp/ewdspace.h>

namespace votca {
namespace xtp {


class PtopReader : public QMTool
{
public:
    
    PtopReader() { };
   ~PtopReader() { };

    string Identify() { return "ptopreader"; }
    void   Initialize(Property *options);
    bool   Evaluate();
    
private:

    string _ptop_file;
    string _ptop_key;
    string _ptop_xml_file;
    string _ptop_tab_file;

};


void PtopReader::Initialize(Property *opt) {
    
    string key = "options.ptopreader";
    _ptop_file = opt->get(key+".ptop_file").as<string>();
    if (opt->exists(key+".ptop_key")) {
        _ptop_key = opt->get(key+".ptop_key").as<string>();
    }
    else _ptop_key = "bgn";
    _ptop_xml_file = _ptop_file+".xml";
    _ptop_tab_file = _ptop_file+".tab";
}


bool PtopReader::Evaluate() {
    
    string pfx = "... ... ";
    using boost::format;
    
    // LOAD BACKGROUND POLARIZATION STATE
    PolarTop ptop = PolarTop(NULL);
    ptop.LoadFromDrive(_ptop_file);
    ptop.RemoveAllOwnership();
    
    vector<PolarSeg*> bgn;
    if (_ptop_key == "bgn")
        bgn = ptop.BGN();
    else if (_ptop_key == "fgc")
        bgn = ptop.FGC();
    else if (_ptop_key == "fgn")
        bgn = ptop.FGN();
    else if (_ptop_key == "qm0")
        bgn = ptop.QM0();
    
    
    cout << endl << pfx << "Background size = " << bgn.size() << flush;
    
    
    std::ofstream ofs;
    ofs.open(_ptop_xml_file.c_str(), ofstream::out);
    ofs << "<ptop>\n";
    
    
    for (vector<PolarSeg*>::iterator sit = bgn.begin();
        sit != bgn.end(); ++sit) {
        vec pos = (*sit)->getPos();
        vec u1_tot = vec(0,0,0);        
        for (PolarSeg::iterator pit = (*sit)->begin();
            pit != (*sit)->end(); ++pit) {            
            u1_tot += (*pit)->getU1();     
        }        
        ofs << "\t<pseg>\n";
        ofs << (format("\t\t<id>%1$d</id>\n") % (*sit)->getId());
        ofs << (format("\t\t<size>%1$d</size>\n") % (*sit)->size());
        ofs << (format("\t\t<pos>%1$1.4f %2$1.4f %3$1.4f</pos>\n") 
            % pos.getX() % pos.getY() % pos.getZ());
        ofs << (format("\t\t<dpl>%1$1.7e %2$1.7e %3$1.7e</dpl>\n") 
            % u1_tot.getX() % u1_tot.getY() % u1_tot.getZ());
        
        for (PolarSeg::iterator pit = (*sit)->begin();
            pit != (*sit)->end(); ++pit) {
            vec pos = (*pit)->getPos();
            vec u1 = (*pit)->getU1();
            ofs << "\t\t<psit>\n";
            ofs << (format("\t\t\t<pos>%1$1.4f %2$1.4f %3$1.4f</pos>\n") 
                % pos.getX() % pos.getY() % pos.getZ());
            ofs << (format("\t\t\t<dpl>%1$1.7e %2$1.7e %3$1.7e</dpl>\n") 
                % u1.getX() % u1.getY() % u1.getZ());            
            ofs << "\t\t</psit>\n";
        }
        
        ofs << "\t</pseg>\n";
    }
    
    ofs << "</ptop>\n";
    ofs.close();
    
    
    ofs.open(_ptop_tab_file.c_str(), ofstream::out);    
    for (vector<PolarSeg*>::iterator sit = bgn.begin();
        sit != bgn.end(); ++sit) {
        vec pos = (*sit)->getPos();
        vec d1_tot = (*sit)->CalcTotD();
        vec u1_tot = vec(0,0,0);        
        for (PolarSeg::iterator pit = (*sit)->begin();
            pit != (*sit)->end(); ++pit) {            
            u1_tot += (*pit)->getU1();     
        }

        ofs << (format("sxyz %1$1.4f %2$1.4f %3$1.4f ") 
            % pos.getX() % pos.getY() % pos.getZ());
        ofs << (format("su1  %1$+1.7e %2$+1.7e %3$+1.7e ") 
            % u1_tot.getX() % u1_tot.getY() % u1_tot.getZ());
        ofs << (format(" sd1  %1$+1.7e %2$+1.7e %3$+1.7e\n") 
            % d1_tot.getX() % d1_tot.getY() % d1_tot.getZ());

        for (PolarSeg::iterator pit = (*sit)->begin();
            pit != (*sit)->end(); ++pit) {
            vec pos = (*pit)->getPos();
            vec u1 = (*pit)->getU1();
            ofs << (format("xyz %1$1.4f %2$1.4f %3$1.4f ") 
                % pos.getX() % pos.getY() % pos.getZ());
            ofs << (format("dpl %1$1.7e %2$1.7e %3$1.7e\n") 
                % u1.getX() % u1.getY() % u1.getZ());            
        }        
    }
    ofs.close();
    return true;
    
    
}





}}

#endif
