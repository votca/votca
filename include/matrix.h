/* 
 * File:   matrix.h
 * Author: ruehle
 *
 * Created on July 11, 2007, 3:24 PM
 */

#ifndef _mat_H
#define	_mat_H

#include "types.h"
#include "ostream"
#include "vec.h"

class matrix
{
public:
        
    matrix() {};
    matrix(const real &v) { *this=v; }
    matrix(const matrix &m) { *this=m; }
    matrix(real  arr[9]) {*this=arr; }
    matrix(const vec& a, const vec& b, const vec& c){
        _m[0]=a.getX(); _m[1]=b.getX(); _m[2]=c.getX();
        _m[3]=a.getY(); _m[4]=b.getY(); _m[5]=c.getY();
        _m[6]=a.getZ(); _m[7]=b.getZ(); _m[8]=c.getZ();
    } // takes three vectors and creates a matrix with them as columns
    

    matrix &operator=(const real &v);
    matrix &operator=(const matrix &v);
    matrix &operator=(real [9]);
    //vec &operator+=(const vec &v);
    //vec &operator-=(const vec &v);
    matrix &operator*=(const double &d){
        for(size_t i=0; i<9; ++i) _m[i] *=d;
        return *this;
    }
    //matrix &operator*(const double &d){
    	
    //}
    matrix &operator/=(const double &d){
        for(size_t i=0; i<9; ++i) _m[i] /=d;
        return *this;
    }  
    matrix &operator-=(const matrix &v){
        for(size_t i=0; i<9; ++i) _m[i] -= v._m[i];
        return *this;
    }
    matrix &operator+=(const matrix &v){
        for(size_t i=0; i<9; ++i) _m[i] += v._m[i];
        return *this;
    }  
    
    void ZeroMatrix();
    void UnitMatrix();
    
    void set(const byte_t &i, const byte_t &j, const real &v) { _m[i*3+j] = v; }
    const real &get(const byte_t &i, const byte_t &j) const { return _m[i*3+j]; }

    vec getRow(const byte_t &i) const { return vec(&_m[i*3]); }
    vec getCol(const byte_t &i) const { return vec(_m[i], _m[i+3], _m[i+6]); }
    
    real *operator[](size_t i) { return &_m[i*3]; }       
    
    struct eigensystem_t {
        real eigenvalues[3];
        vec eigenvecs[3];
        
        eigensystem_t operator+=(const eigensystem_t &e) {
            eigenvalues[0]+=e.eigenvalues[0];
            eigenvalues[1]+=e.eigenvalues[1];
            eigenvalues[2]+=e.eigenvalues[2];
            eigenvecs[0]+=e.eigenvecs[0];
            eigenvecs[1]+=e.eigenvecs[1];
            eigenvecs[2]+=e.eigenvecs[2];
        }        

        eigensystem_t operator*=(const double &f) {
            eigenvalues[0]*=f;
            eigenvalues[1]*=f;
            eigenvalues[2]*=f;
            eigenvecs[0]*=f;
            eigenvecs[1]*=f;
            eigenvecs[2]*=f;
        }        
        
        void zero() {
            eigenvalues[0]=eigenvalues[1]=eigenvalues[2]=0;
            eigenvecs[0]=eigenvecs[1]=eigenvecs[2]=vec(0.,0.,0.);
        }
    };
    void RandomRotation();
    void SolveEigensystem(eigensystem_t &out);
    
    matrix &Transpose(){
        std::swap( _m[1], _m[3]);
        std::swap( _m[2], _m[6]);
        std::swap( _m[5], _m[7]);
        return *this;
    }
    ///matrix multiplication
    matrix  operator * (const matrix & a){
        matrix r;
        r._m[0] = _m[0] * a._m[0] + _m[1] * a._m[3] + _m[2] * a._m[6];
        r._m[1] = _m[0] * a._m[1] + _m[1] * a._m[4] + _m[2] * a._m[7];
        r._m[2] = _m[0] * a._m[2] + _m[1] * a._m[5] + _m[2] * a._m[8];
        
        r._m[3] = _m[3] * a._m[0] + _m[4] * a._m[3] + _m[5] * a._m[6];
        r._m[4] = _m[3] * a._m[1] + _m[4] * a._m[4] + _m[5] * a._m[7];
        r._m[5] = _m[3] * a._m[2] + _m[4] * a._m[5] + _m[5] * a._m[8];
        
        r._m[6] = _m[6] * a._m[0] + _m[7] * a._m[3] + _m[8] * a._m[6];
        r._m[7] = _m[6] * a._m[1] + _m[7] * a._m[4] + _m[8] * a._m[7];
        r._m[8] = _m[6] * a._m[2] + _m[7] * a._m[5] + _m[8] * a._m[8];
        
        return r;
    }
       
    vec operator * ( const vec & a){
       return vec( _m[0] * a.getX() + _m[1] * a.getY() + _m[2] * a.getZ(), 
                _m[3] * a.getX() + _m[4] * a.getY() + _m[5] * a.getZ(),
                _m[6] * a.getX() + _m[7] * a.getY() + _m[8] * a.getZ() ); 
    }
    friend matrix operator*(const double &, const matrix &);
    friend vec operator*(const vec &, const matrix &);
	  private:
    real _m[9];
};

inline matrix &matrix::operator=(const real &v)
{
    for(size_t i=0; i<9; ++i)
        _m[i] = v;
    return *this;
}

inline matrix &matrix::operator=(const matrix &m)
{
    for(size_t i=0; i<9; ++i)
        _m[i] = m._m[i];
    return *this;
}

inline matrix &matrix::operator=(real arr[9])
{
    for(size_t i=0; i<9; ++i)
        _m[i] = arr[i];
    return *this;
}

inline void matrix::UnitMatrix()
{
    ZeroMatrix();
    _m[0] = _m[4] = _m[8] = 1.0;
}

inline void matrix::ZeroMatrix()
{
    for(size_t i=0; i<9; ++i)
        _m[i] = 0.;//(*this) = 0.;
    _m[0] = _m[4] = _m[8] = 1.0;
}

inline std::ostream &operator<<(std::ostream &out, matrix& m)
{
      out << '|' << m[0][0] << ',' << m[0][1] << ',' << m[0][2] << '|' << std::endl;
      out << '|' << m[1][0] << ',' << m[1][1] << ',' << m[1][2] << '|' << std::endl;
      out << '|' << m[2][0] << ',' << m[2][1] << ',' << m[2][2] << '|' << std::endl;
      return out;
}

inline matrix operator|( const vec & a, const vec & b){
    matrix res;
    res.set(0,0, a.getX() * b.getX());
    res.set(0,1, a.getX() * b.getY());
    res.set(0,2, a.getX() * b.getZ());
    res.set(1,0, a.getY() * b.getX());
    res.set(1,1, a.getY() * b.getY());
    res.set(1,2, a.getY() * b.getZ());
    res.set(2,0, a.getZ() * b.getX());
    res.set(2,1, a.getZ() * b.getY());
    res.set(2,2, a.getZ() * b.getZ());
    return res;
}


inline matrix operator*(const matrix & r, const double &d){
    return ( matrix(r) *= d);
}

inline matrix operator*(const double &d, const matrix &m)
{
        matrix mm;
	for(size_t i=0; i<9; ++i)
	    mm._m[i] = d*m._m[i];
	return mm;
}

inline vec operator*(const matrix & r, const vec &a){
        return matrix(r) * a;
}


inline matrix operator/(const matrix & r,const double &d){
    return ( matrix(r) /= d);
}  
inline matrix operator+(const matrix & r, const matrix & v){
    return ( matrix(r) += v);
}
inline matrix operator-(const matrix & r, const matrix & v){
    return ( matrix(r) -= v);
}
    

/* provided the matrix a diagonalizes it and returns the eigenvalues 
   lambda0 = a[0][0], lambda1 = a[1][1], lambda2= a[2][2], ...
   as well as the corresponding eigenvectors v[][0], v[][1], v[][2] 
*/
int cjcbi(matrix &a, matrix &v, double eps=1e-10, int jt=100);



#endif	/* _matrix_H */

