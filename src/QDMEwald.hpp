#ifndef QDMEWALD_INCLUDED
#define QDMEWALD_INCLUDED

#include "KSpace.hpp"
#include "RSpace.hpp"

template <class T> QDMEwald
{
    public:
        /**
         *
         * Constructor for the class (base-version)
         *
         */
        QDMEwald(){}
        /**
         *
         * Constructor for the class
         *
         * @param alpha     Ewald sum splitting parameter
         * @param k_max     floating point cut-off radius for the k-space
         *                  part
         * @param r_max     floating point cut-off radius for the r-space
         *                  part
         */
        QDMEwald( T alpha, T k_max, T r_max ) : kspace(alpha,k_max), rspace(alpha,r_max){}

    private:
        KSpace kspace;
        RSpace rspace;
};

#endif
