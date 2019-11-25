#ifndef KSPACE_CLASS_INCLUDED
#define KSPACE_CLASS_INCLUDED

#include <vector>
#include <cmath>
#include "Kokkos_Core.hpp"

/**
 *
 * Class to compute the k-space part of an Ewald sum for
 * quadrupole, dipole and/or monopole systems
 *
 */

template <class T> class KSpace
{
    public:
        /**
         *
         * Constructor for the class (base-version)
         *
         */
        KSpace(){}
        /**
         *
         * Constructor for the class
         *
         * @param alpha     Ewald sum splitting parameter
         * @param k_max     floating point cut-off radius for the k-space
         *                  part
         */
        KSpace( const T _alpha, const T _k_max ) { init_params(_alpha, _k_max); }
        /*
         *
         * Initialization and computation of internal parameters
         *
         * @param alpha     Ewald sum splitting parameter
         * @param k_max     floating point cut-off radius for the k-space
         *                  part
         */
        void init_params( const T, const T );
        /**
         *
         * Computation of the kspace part of the Ewald sum, fills the
         * private variables of the class which can be accessed by the
         * corresponding getters after the computation
         *
         * @param l     system length (cubic system)
         * @param xyz   particle positions (xyz(0), xyz(1), ...)
         * @param q     particle charges(q,q,q,...)
         * @param d     particle dipole moments(d(0,0),d(0,1),
         *              d(1,0),d(1,1),d(2,0),d(2,1),...)
         * @param q     particle quadrupole moments( q(0,0),...,q(0,8),q(1,0),...,q(1,8),...)
         */
        void compute( const T,
                      const std::vector<T>&, 
                      const std::vector<T>&,
                      const std::vector<T>&,
                      const std::vector<T>& );

    private:
        T alpha;                        ///< splitting parameter of the Ewald sum
        T gamma;                        ///< transformed splitting parameter of the Ewald sum
        int k_sq_int;                   ///< square value of the k-space limit
        int k_max_int;                  ///< k-space image limit

        T pot_energy;                   ///< potential energy computed for the k-space part
        T virial;                       ///< virial computed for the k-space part
        Kokkos::View<T*[3]> f;          ///< forces computed for the k-space part
        Kokkos::View<T*[3]> tqe;        ///< torque computed for the k-space part

        bool is_monopole;               ///< compute the monopole contributions
        bool is_dipole;                 ///< compute the dipole contributions
        bool is_quadrupole;             ///< compute the quadrupole contributions

        Kokkos::View<T***> cos_fac;     ///< cosine based exponential factors
        Kokkos::View<T***> sin_fac;     ///< sine based exponential factors

        Kokkos::View<T*> ak;            ///< AK coefficients 

        /*
         *
         * Compute the exponential factors for each particle
         *
         * @param xyz   std::vector of type T containing the positions of the
         *              particles to be used in the computation of the exponential
         *              factors
         * @param l     system length
         */
        void compute_exponentials(Kokkos::View<T*[3]>, const T);        

        /*
         *
         * Computation of the AK coefficients for the matrix multiplications 
         *
         * @param l     system length
         */
        void compute_ak(const T);         
};

/*
 *
 * Initialization and computation of internal parameters
 *
 * @param alpha     Ewald sum splitting parameter
 * @param k_max     floating point cut-off radius for the k-space
 *                  part
 */
template <class T> void KSpace<T>::init_params( const T _alpha, const T _k_max )
{
    alpha = _alpha;
    gamma = -0.25 / ( alpha * alpha );
    k_max_int = (int)_k_max;
    k_sq_int = (int)std::floor(sqrt(_k_max));
}

/*
 *
 * Computation of the AK coefficients for the matrix multiplications 
 * (currently only implemented for host-space)
 *
 * @param l         system size
 */
template <class T> void KSpace<T>::compute_ak( 
                                                const T l
                                             )
{
    T expf = 0.0;

    // transformed length
    T rcl = 2.0 * M_PI / l;

    if (alpha > (T)1e-12)
    {
        expf = std::exp( gamma * rcl * rcl ); 
    }

    ak = Kokkos::View<T*>(
                            "AK coefficients",
                            k_sq_int 
                         );
    
    Kokkos::parallel_for(k_sq_int, KOKKOS_LAMBDA(const int k)
        {
            T rksq = (T)k * rcl * rcl;
            T eksq = std::pow(expf,(T)k);
            ak(k) = eksq / rksq; 
        }
    );
}

/*
 *
 * Computation of the exponential factors for the k-space part
 * (currently only implemented for host-space)
 *
 * @param xyz       particles positons
 * @param l         system size
 */
template <class T> void KSpace<T>::compute_exponentials( 
                                                            Kokkos::View<T*[3]> xyz,
                                                            const T l
                                                       )
{
    // get number of particles
    size_t N = xyz.size() / 3.0;
    // create Kokkos views of sufficient size to store the factors
    cos_fac = Kokkos::View<T***>(
                                    "cosine exponential factors",
                                    N,
                                    2*(k_max_int)+1,
                                    3
                                 );
    sin_fac = Kokkos::View<T***>(
                                    "sine exponential factors",
                                    N,
                                    2*(k_max_int)+1,
                                    3
                                 );

    // to deal with negative k-values
    int offset = k_max_int;

    // initialize the first factors (k == 0)
    Kokkos::parallel_for(N, KOKKOS_LAMBDA(const int n)
        {
            for (int d = 0; d < 3; ++d)
            {
                cos_fac(n,0+offset,d) = 1.0;
                sin_fac(n,0+offset,d) = 0.0;
            }
        }
    );

    // transformed length
    T rcl = 2.0 * M_PI / l;

    // compute the exponential factors (k == 1 / k == -1)
    Kokkos::parallel_for(N, KOKKOS_LAMBDA(const int n)
        {
            for(int d = 0; d < 3; ++d)
            {
                cos_fac(n,1+offset,d) = std::cos(rcl * xyz(3*n,d));
                sin_fac(n,1+offset,d) = std::sin(rcl * xyz(3*n,d));
                cos_fac(n,-1+offset,d) = cos_fac(n,1+offset,d);
                sin_fac(n,-1+offset,d) = -sin_fac(n,1+offset,d);
            }
        }
    );

    Kokkos::parallel_for(N, KOKKOS_LAMBDA(const int n)
        {
            for(int k = 2; k <= k_max_int; ++k)
            {
                for(int d = 0; d < 3; ++d)
                {
                    cos_fac(n,k+offset,d) = cos_fac(n,k-1+offset,d) * cos_fac(n,1+offset,d) -
                                            sin_fac(n,k-1+offset,d) * sin_fac(n,1+offset,d);
                    sin_fac(n,k+offset,d) = sin_fac(n,k-1+offset,d) * cos_fac(n,1+offset,d) -
                                            cos_fac(n,k-1+offset,d) * sin_fac(n,1+offset,d);
                    cos_fac(n,-k+offset,d) = cos_fac(n,k+offset,d);
                    sin_fac(n,-k+offset,d) = -sin_fac(n,k+offset,d);
                }
            }
        }
    );

}

/**
 *
 * Computation of the kspace part of the Ewald sum, fills the
 * private variables of the class which can be accessed by the
 * corresponding getters after the computation
 *
 * @param _l     system length (cubic system)
 * @param _xyz   particle positions (xyz(0), xyz(1), ...)
 * @param _q     particle charges(q,q,q,...)
 * @param _d     particle dipole moments(d(0,0),d(0,1),
 *              d(1,0),d(1,1),d(2,0),d(2,1),...)
 * @param _Q     particle quadrupole moments( Q(0,0),...,Q(0,8),Q(1,0),...,Q(1,8),...)
 */
template <class T> void KSpace<T>::compute( 
                                            const T _l,
                                            const std::vector<T>& _xyz,
                                            const std::vector<T>& _q,
                                            const std::vector<T>& _d,
                                            const std::vector<T>& _Q
                                          )
{
    // get number of particles
    size_t N = _xyz.size() / 3;

    // store C++ vectors to Kokkos Views
    Kokkos::View<T*[3]> xyz("positions",N);
    Kokkos::View<T*> q("charges",N);
    Kokkos::View<T*[3]> d("dipole moments",N);
    Kokkos::View<T*[3][3]> Q("quadrupole moments",N);

    Kokkos::parallel_for(N, KOKKOS_LAMBDA( const int n )
        {
            q(n) = _q.at(n);
            for (int i = 0; i < 3; ++i)
            {
                xyz(n,i) = _xyz.at(3*n+i);
                d(n,i) = _d.at(3*n+i);
                for (int j = 0; j < 3; ++j)
                {
                    Q(n,i,j) = _Q.at(9*n + 3*i + j);
                }
            }
        }
    );

    // compute exponential factors
    compute_exponentials(xyz,_l);

    // compute AK coefficients
    compute_ak(_l);
}

#endif
