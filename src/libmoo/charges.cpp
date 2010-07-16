#include "charges.h"

int multipoles::read_crg_eps(const char * namefile ){
        ifstream in(namefile);
        int i=0;
        string word;
	double mpl;

        int count = 0;
        while( in  ) {
                count ++;
                in >> word;
                if (count %3 == 0){
                    sscanf(word.c_str(), "%lf", &mpl);
                    mpls.push_back(mpl);
                }
        }
        in.close();
        return 0;
}



