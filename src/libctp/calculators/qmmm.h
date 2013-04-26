#ifndef __QMMMCALC__H
#define	__QMMMCALC__H



namespace votca { namespace ctp {
    
    
class QMMM : public QMCalculator
{

public:
    
    QMMM() {};
   ~QMMM() {};
   
    string  Identify() { return "QMMM"; }
    void    Initialize(Topology *top, Property *options);
    bool    EvaluateFrame(Topology *top);
   
   

private:



};




void QMMM::Initialize(Topology *top, Property *options) {
    
    /*
     
     <qmmm>
     *  <mapper>
     *    <mpsmapxml>
     *    <segmpstab>
     *  </mapper>
     * 
     
     
     
     */
    ;
}


bool QMMM::EvaluateFrame(Topology *top) {
    
    return true;
}

    
    
    
    
}}

#endif /* __QMMM__H */