#ifndef _VOTCA_XTP_EIMPORT_H
#define _VOTCA_XTP_EIMPORT_H

#include <votca/xtp/qmcalculator.h>


namespace votca { namespace xtp {

class EImport : public QMCalculator
{
public:

    EImport() {};
   ~EImport() {};

   string Identify() { return "eimport"; }

   void   Initialize(Property *options);
   bool   EvaluateFrame(Topology *top);
   void   StochasticEnergies(Topology *top, string &_probabilityfile, int state);
   double SphereIntersection(double radius, double distance);
   double Correlation(double distance, vector<double> distances, vector<double>bcoeff, double sigma);


private:

    string _energiesFile;
    bool   _reset;
    string      _probabilityfile_h;
    string      _probabilityfile_e;
    double      _sigma_h;
    double      _sigma_e;
    double      _avgestatic_h;
    double      _avgestatic_e;
    double      _cutoff;
    bool        _stochastic;
    int         _seed;
    int         _optimizationsteps;
    double      _optimizationfactor;


};


void EImport::Initialize(Property *options) {
    _stochastic     = false;

    
    // update options with the VOTCASHARE defaults   
    UpdateWithDefaults( options );
    string key = "options." + Identify();

    if (options->exists(key + ".energies")) {
         _energiesFile = options->get(key + ".energies").as< string >();
    }
    
    if (options->exists(key + ".probabilityfile_h")) {
        _probabilityfile_h = options->get(key+".probabilityfile_h").as< string >();
        _stochastic = true;
    }
    else{
        _probabilityfile_h = "";
    }
    if (options->exists(key + ".probabilityfile_e")) {
        _probabilityfile_e = options->get(key+".probabilityfile_e").as< string >();
        _stochastic = true;
    }
    else{
        _probabilityfile_e = "";
    }
     if (_stochastic == true) {
        cout << endl << "... ... Creating stochastic site energies based on provided spatial correlation file(s)." << flush;
    }
    if (options->exists(key + ".sigma_h")) {
        _sigma_h = options->get(key+".sigma_h").as< double >();
    }
    else{
        _sigma_h = 1;
    }
    if (options->exists(key + ".sigma_e")) {
        _sigma_e = options->get(key+".sigma_e").as< double >();
    }
    else{
        _sigma_e = 1;
    }
    if (options->exists(key + ".avgestatic_h")) {
        _avgestatic_h = options->get(key+".avgestatic_h").as< double >();
    }
    else{
        _avgestatic_h = 0;
    }
    if (options->exists(key + ".avgestatic_e")) {
        _avgestatic_e = options->get(key+".avgestatic_e").as< double >();
    }
    else{
        _avgestatic_e = 0;
    }
    if (options->exists(key + ".cutoff")) {
        _cutoff = options->get(key+".cutoff").as< double >();
    }
    else{
        _cutoff = 0;
    }
    _seed = 1;
    if (options->exists(key + ".seed")) {
         _seed = options->get(key + ".seed").as< int >();
    }
    _optimizationsteps = 10000;
    if (options->exists(key + ".optimizationsteps")) {
         _optimizationsteps = options->get(key + ".optimizationsteps").as< int >();
    }
    _optimizationfactor = 0.97;
    if (options->exists(key + ".optimizationfactor")) {
         _optimizationfactor = options->get(key + ".optimizationfactor").as< double >();
    }
    


    

    if (options->exists(key+".reset")) {
        int reset = options->get(key+".reset").as<int>();
        _reset = (reset == 1) ? true : false;
        if (_reset) {
            cout << endl
                 << "... ... NOTE: Resetting site energies to zero."
                 << flush;
        }

    }
    else {
        _reset = false;
    }
}

double EImport::SphereIntersection(double radius, double distance){
    double overlapratio;
    if(distance <= 2*radius){
        overlapratio = 1./16. * (4.*radius+distance) * pow((2.*radius-distance),2.) / pow(radius,3.);
        //cout << "radius = " << radius << endl;
        //cout << "distance = " << distance << endl;
        //cout << "intersectvolume = " << intersectvolume << endl;
    }
    else{
        //cout << "doing nothing " << endl;
        overlapratio = 0;
    }
    return overlapratio;
}

double EImport::Correlation(double distance, vector<double> distances, vector<double>bcoeff, double sigma){
        double correlation = 0;
        for(unsigned j=0; j<bcoeff.size(); j++){
            correlation += SphereIntersection(distances[j+1]/2.,distance) * bcoeff[j] ;
        }
        correlation /=  (sigma*sigma);
        // cout << "correlation "<< correlation << endl;
        return correlation;
}

void EImport::StochasticEnergies(Topology *top, string &filename, int state) {
    double sigma = 1;
    double avgestatic = 0;
    if(state == 1){
        cout << endl << "... ... calculating stochastic hole-type energies." << endl;
        sigma = _sigma_h;
        avgestatic = _avgestatic_h;
        cout << endl << "... ... sigma = " << sigma << " eV" << endl;
    }
    else if(state == -1){
        cout << endl << "... ... calculating stochastic electron-type energies." << endl;
        sigma = _sigma_e;
        avgestatic = _avgestatic_e;
        cout << endl << "... ... sigma = " << sigma << " eV" << endl;
    }
    // read in probability function
    cout << "... ... reading in spatial correlation function from file "+filename ;
    if (_cutoff != 0) {cout << " up to cutoff " << _cutoff << " nm";}
    cout << "." << endl;
    vector<double> distances;
    vector<double> kappas;

    std::string line;
    std::ifstream intt;
    intt.open(filename.c_str());
    int linenumber = 0;
    int Nkappa = 0;
    if (intt.is_open() ) {
        while ( intt.good() ) {

            std::getline(intt, line);
            
            if (linenumber > 1){
                vector<string> split;
                Tokenizer toker(line, " \t");
                toker.ToVector(split);

                if ( !split.size()      ||
                      split[0] == "!"   ||
                      split[0].substr(0,1) == "!" ) { continue; }
             
                double distance          = boost::lexical_cast<double>(split[0]);
                double kappa             = boost::lexical_cast<double>(split[1]);
                // double sigma             = boost::lexical_cast<double>(split[2]);
                if((distance<=_cutoff || _cutoff == 0) && kappa >= 0){
                    distances.push_back(distance);
                    kappas.push_back(kappa);
                    cout << "        " << distance << " nm:   " << kappa << endl;
                    Nkappa ++;
                }
            }
            linenumber++;
        }
        Nkappa -= 1;
    }
    else { cout << endl << "ERROR: No such file " << filename << endl;
           throw std::runtime_error("Supply input probability file."); }  
    
    
     // Initialise random number generator
    if(votca::tools::globals::verbose) { cout << endl << "Initialising random number generator with seed: " << _seed << endl; }
    cout << endl << "Initialising random number generator with seed: " << _seed << endl;
    srand(_seed); 
    votca::tools::Random2 *RandomVariable = new votca::tools::Random2();
    RandomVariable->init(rand(), rand(), rand(), rand());
    
     // calculate b_i coefficients    
    cout << "... ... calculating b coefficients" << endl;
    cout << "... ... If the program freezes please set 'optimizationfactor' to a smaller value in your options file." << endl;

    vector<double> bcoeff;
    bcoeff.resize(Nkappa+1);
    distances.push_back(2*distances[Nkappa]-distances[Nkappa-1]);
    // cout << "distances[0]=" << distances[0] << ", distances[1]=" << distances[1] << endl;
            
    bcoeff[Nkappa] = sigma*sigma / SphereIntersection(distances[Nkappa+1]/2,distances[Nkappa]) * kappas[Nkappa];
    // cout << "        b[" << Nkappa+1 << "] = " << bcoeff[Nkappa] << endl;
    
    bool interpolatable = true;
    for(int nit = Nkappa-1; nit >=0; nit--){
        cout << "\r... ...       b[" << nit+1 << "] ..." << flush;
        // cout << endl << endl << "b[" << nit+1 << "] (distance " << distances[nit] << ")"<< endl; 
        double sum = 0;

        for (int j=nit+1; j<= Nkappa; j++){
            sum += SphereIntersection(distances[j+1]/2.,distances[nit]) * bcoeff[j];
            // cout << "  " << j+1 << ": " << bcoeff[j] << endl;
        }
        bcoeff[nit] = 1/SphereIntersection(distances[nit+1]/2.,distances[nit]) * (sigma*sigma*kappas[nit] - sum);
        if(bcoeff[nit] < 0){interpolatable = false;}
        // cout << "        b[" << nit+1 << "] = " << bcoeff[nit] << endl;

        while(bcoeff[nit] < 0){
            for(int k = nit+1; k <= Nkappa; k++){
                if(std::abs(bcoeff[k]) < 1E-5) {bcoeff[k] = 0;}
                if(bcoeff[k] == 0){break;}
                else{
                    bcoeff[k] *= _optimizationfactor;
                    if(std::abs(bcoeff[k]) < 1E-5) {bcoeff[k] = 0;}
                    // cout << "changing b[" << k+1 << "] = " << bcoeff[k] << endl;
                    double sum = 0;
                    for (int j=nit+1; j<= Nkappa; j++){
                        sum += SphereIntersection(distances[j+1]/2.,distances[nit]) * bcoeff[j];
                    }
                    bcoeff[nit] = 1/SphereIntersection(distances[nit+1]/2.,distances[nit]) * (sigma*sigma*kappas[nit] - sum);
                    // cout << "        b[" << nit+1 << "] = " << bcoeff[nit] << endl;
                    if (bcoeff[nit] >= 0) {break;}
                }
            }
        }
    }

    
    if(interpolatable == false){
        // try to improve solution by slightly modifying the parameters
        double mse = 0;
        for(unsigned j = 0; j<distances.size(); j++){
            mse += pow(Correlation(distances[j], distances, bcoeff, sigma) - kappas[j], 2);
        }
        cout << "\n        NOTE: The given correlation data could not exactly be reproduced with the intersecting sphere model." << endl;
        cout << "              An approximative solution will be calculated. Check the expected correlation function" << endl;
        cout << "              against the result. If the approximation is bad it might help to choose less interpolation" << endl;
        cout << "              points, i.e., choose a larger spatial resolution in eanalyze." << endl;
        cout << "              If the long range part is not reproduced well it might help to remove the first point," << endl;
        cout << "              since this algorithm focuses on a good reproduction of the short range part." << endl << endl;
        cout << "... ...       initial mean square error: " << mse  << endl;
        cout << "... ...       optimizing result..." << endl;
        
        
        for(int iteration = 0; iteration < _optimizationsteps; iteration ++){
            cout << "\r... ...       " << int(double(iteration+1)/double(_optimizationsteps)*100.) << " % done" << flush;
            // which coefficient to change
            int nit = rand() % Nkappa;
            vector<double>btry = bcoeff;
            
            for(double i = -100; i<=100; i++){
                double diff = double(i)*std::abs(bcoeff[nit])/1000.;
                btry[nit] = bcoeff[nit] + diff*bcoeff[nit];
                double thismse = 0;
                for(unsigned j = 0; j<distances.size(); j++){
                    thismse += pow(Correlation(distances[j], distances, btry, sigma) - kappas[j], 2);
                }
                if(thismse<mse && btry[nit] >= 0){
                    mse = thismse;
                    bcoeff = btry;
                }
            }
            
        }
        cout << endl << "... ...       final mean square error: " << mse << endl << endl;
    }
    
    
    double sum =0;
    for(unsigned i=0; i<bcoeff.size(); i++){sum += bcoeff[i];}
    double acoeff = sigma*sigma-sum;
    if(acoeff<0){
        acoeff = 0;
        cout << "        WARNING: a was smaller than 0. It has now been set to 0, the energetic disorder will probably be wrong." << endl;
    }
    cout << "        a = " << acoeff << endl << endl;
    
    for(int nit = 0; nit <=Nkappa; nit++){
         cout << "        b[" << nit+1 << "] = " << bcoeff[nit] << endl;
    }

   
    FILE *out;
    string tag = boost::lexical_cast<string>("eimport.sitecorr_") + ( (state == -1) ? "e" : "h" ) + "_compare.out";
    out = fopen(tag.c_str(), "w");

    fprintf(out, "# EIMPORT: COMPARISION OF IMPORTED AND REALISED CORRELATION.\n");
    fprintf(out, "# STATE %1d\n", state);

    for (unsigned i = 0; i < distances.size()-1; ++i) {
        fprintf(out, "%4.7f %4.7f %4.7f\n", distances[i], kappas[i], Correlation(distances[i], distances, bcoeff, sigma));
    }
    fclose(out);

    tag = boost::lexical_cast<string>("eimport.sitecorr_") + ( (state == -1) ? "e" : "h" ) + "_interpolation.out";
    out = fopen(tag.c_str(), "w");
    cout << "        An interpolated expected correlation function will be written out to " << tag.c_str() << endl;
    // cout << "        Specify borders: " << endl;
    // cout << "        minimum (data: " << distances[0] << " nm): ";
    double min = 0.1;
    // std::cin >> min;
    double max = 8;
    // cout << "        maximum (data: " << distances[Nkappa] << " nm): ";
    // std::cin >> max;
    

    fprintf(out, "# INTERPOLATING STOCHASTIC CORRELATION FUNCTION.\n");
    fprintf(out, "# STATE %1d\n", state);

    for (int i = 0; i < 1000; ++i) {
        double thisr = min+ double(i)/1000.*(max-min);
        fprintf(out, "%4.7f %4.7f\n", thisr, Correlation(thisr, distances, bcoeff, sigma));
    }
    fclose(out);

    
    // generate site energies
    vector <Segment*> segments = top->Segments();
    vector <Segment*>::iterator seg1;
    vector <Segment*>::iterator seg2;
    
    int Nmolecules = segments.size();
    // calculate random numbers
    cout << endl << "... ... generating "  << (Nkappa+1+1)*Nmolecules << " random numbers ...";
    vector< vector<double> > X;
    for(int A = 0; A<Nmolecules; A++){
        vector<double> dummy;
        X.push_back(dummy);
        for(unsigned i = 0; i<bcoeff.size()+1; i++){
            X[A].push_back( RandomVariable->rand_gaussian(1.0) );
        }
    }
    cout << " done." << endl;

    
    
    double AVG = 0.0;
    int countE = 0;
    vector<double> energies;
    for (seg1 = segments.begin(); seg1!= segments.end(); seg1++){
        cout << "\r... ... ..." << " calculating energy for segment ID = "
             << (*seg1)->getId() << flush;
        int molA = (*seg1)->getId()-1;
        
        vector<double> molsinrange(bcoeff.size(), 0.0);
        vector<double> randompart(bcoeff.size()+1,0.0);

        for (seg2 = segments.begin(); seg2!= segments.end(); seg2++){
            int molB = (*seg2)->getId()-1;
            vec r1 = (*seg1)->getPos();
            vec r2 = (*seg2)->getPos();
            double distance = abs( top->PbShortestConnect(r1, r2));
            for(unsigned nit=0; nit<bcoeff.size(); nit++){
                if(distance<=distances[nit+1]/2 ) { // r=0 included ?????
                    molsinrange[nit] += 1;
                    randompart[nit]  += X[molB][nit];
                    
                }
            }
            
        }
        
        double energy = sqrt(acoeff) * X[molA][bcoeff.size()];
        for(unsigned i =0; i<bcoeff.size(); i++){
            if(molsinrange[i] > 0){
                energy += sqrt(bcoeff[i]/molsinrange[i]) * randompart[i];
            }
            else {
                cout << "NOTHING IN RANGE FOR RADIUS " << distances[i] << " for molecule " << molA+1 << endl;
            }
        }
        if(molA <= 5 || molA >= Nmolecules-5){ cout << "   energy[" << molA+1 << "] = " << energy << endl;}
        if(molA == 5 ){ cout << "......" << endl;}
        AVG += energy;
        energies.push_back(energy);
        countE += 1;
    }
    if(countE > 0){
    AVG /= countE;
    }

    double VAR = 0.0;
    double STD = 0.0;
    // Calculate variance
    vector< double > ::iterator eit;
    for (eit = energies.begin(); eit < energies.end(); ++eit) {
        VAR += ((*eit) - AVG)*((*eit) - AVG) / countE;
    }
   
    STD = sqrt(VAR);
    
    cout << "\r... ... ..." << " uncorrected mean:  " << AVG << endl;
    cout << "\r... ... ..." << " uncorrected sigma: " << STD << endl;
    
    // Write corrected energies
    cout << "\r... ... ..." << " Setting energies to adjusted standard deviation." << endl;
    countE = 0;
    for (seg1 = segments.begin(); seg1!= segments.end(); seg1++){
        (*seg1)->setEMpoles(state, (energies[countE]*sqrt(sigma/STD)-AVG+avgestatic));
        countE ++;
    }


    cout << "done." << endl;
}

bool EImport::EvaluateFrame(Topology *top) {

    if (_reset == true) {

        vector<Segment*> ::iterator sit;
        for (sit = top->Segments().begin();
             sit < top->Segments().end();
             ++sit) {

            (*sit)->setEMpoles( 0, 0.0);
            (*sit)->setEMpoles(-1, 0.0);
            (*sit)->setEMpoles(+1, 0.0);
        }
    }
    
    else if(_stochastic == true)  {
      if(_probabilityfile_h != ""){this->StochasticEnergies(top, _probabilityfile_h, 1);}
      if(_probabilityfile_e != ""){this->StochasticEnergies(top, _probabilityfile_e, -1);}
    }
        
    

    else {

        string filename = _energiesFile;
        std::string line;
        std::ifstream intt;
        intt.open(filename.c_str());

        if (intt.is_open() ) {

            while ( intt.good() ) {

                std::getline(intt, line);
                vector<string> split;
                Tokenizer toker(line, " \t");
                toker.ToVector(split);

                if ( !split.size()      ||
                      split[0] == "!"   ||
                      split[0].substr(0,1) == "!" ) { continue; }


              // Sample line
              // 10 Alq3 0 0.00000000  +1 -0.49966334  0 20  +1 10  SPH 2366 ...

                int    id   = boost::lexical_cast<int>(split[0]);
                string name = boost::lexical_cast<string>(split[1]);

                Segment *seg = top->getSegment(id);

                if (seg->getName() != name) {
                    cout << endl
                         << "... ... ERROR: Structure of input file does "
                            "not match state file. "
                         << "SEG ID " << id
                         << flush;
                    continue;
                }

                int    tell = boost::lexical_cast<int>(split[6]);
                if (tell*tell == 0) {

                    int state_N = boost::lexical_cast<int>(split[2]);
                    int state_C = boost::lexical_cast<int>(split[4]);

                    assert(state_N == 0);
                    assert(state_C*state_C == 1);

                    double e_N = boost::lexical_cast<double>(split[3]);
                    double e_C = boost::lexical_cast<double>(split[5]);

                    seg->setEMpoles(state_N, e_N);
                    seg->setEMpoles(state_C, e_C);
                }

                else if (tell*tell == 1) {

                    int state_N = boost::lexical_cast<int>(split[2]);
                    int state_A = boost::lexical_cast<int>(split[4]);
                    int state_C = boost::lexical_cast<int>(split[6]);

                    assert(state_N == 0);
                    assert(state_A == -1);
                    assert(state_C == +1);

                    double e_N = boost::lexical_cast<double>(split[3]);
                    double e_A = boost::lexical_cast<double>(split[5]);
                    double e_C = boost::lexical_cast<double>(split[7]);

                    seg->setEMpoles(state_N, e_N);
                    seg->setEMpoles(state_A, e_A);
                    seg->setEMpoles(state_C, e_C);
                }

                else {
                    cout << endl << "... ... ERROR: Wrong input format. "
                         << "SEG ID " << id << flush;
                    continue;
                }
            }
        }
        else { cout << endl << "ERROR: No such file " << filename << endl; }
    }

    return 1;

}




}}


#endif
