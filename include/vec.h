/* 
 * Copyright 2009 The VOTCA Development Team (http://www.votca.org)
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

#ifndef _vec_H
#define	_vec_H

#include <iostream>
#include <math.h>

/**
    \brief Vector class for a 3 component vector

    This class represents a 3 component vector to store e.g. postitions, velocities, forces, ...
    Operators for basic vector-vector and vector-scalar operations are defined.
    you can access the elements with the functions x(), y(), z(), both reading and writing is possible;
    x + v.x();
    v.x() = 5.;
*/

class vec {
public:
    
    vec();
    vec(const vec &v);
    vec(const double r[3]);
    vec(const double &x, const double &y, const double &z);
    
    vec &operator=(const vec &v);
    vec &operator+=(const vec &v);
    vec &operator-=(const vec &v);
    vec &operator*=(const double &d);
    vec &operator/=(const double &d);
    
    double &x() { return _x; }
    double &y() { return _y; }
    double &z() { return _z; }
    
    void setX(const double &x) { _x = x; }
    void setY(const double &y) { _y = y; }
    void setZ(const double &z) { _z = z; }
    
    const double &getX() const { return _x; }
    const double &getY() const { return _y; }
    const double &getZ() const { return _z; }
    
    vec &normalize();
    
    private:
        double _x, _y, _z;
};

inline vec::vec() {}

inline vec::vec(const vec &v)
    : _x(v._x), _y(v._y), _z(v._z) {}
        
inline vec::vec(const double r[3])
    : _x(r[0]), _y(r[1]), _z(r[2]) {}
    
inline vec::vec(const double &x, const double &y, const double &z)
        : _x(x), _y(y), _z(z) {}
    
inline bool operator==(const vec &v1, const vec &v2)
{
    return ((v1.getX()==v2.getX()) && (v1.getY()==v2.getY()) && (v1.getZ()==v2.getZ()));
}

inline bool operator!=(const vec &v1, const vec &v2)
{
    return ((v1.getX()!=v2.getX()) || (v1.getY()!=v2.getY()) || (v1.getZ()==v2.getZ()));
}

inline vec &vec::operator=(const vec &v)
{ 
        _x=v._x; _y=v._y; _z=v._z;
        return *this;
}    

inline vec &vec::operator+=(const vec &v)
{ 
        _x+=v._x; _y+=v._y; _z+=v._z;
        return *this;
}    
        
inline vec &vec::operator-=(const vec &v)
{ 
        _x-=v._x; _y-=v._y; _z-=v._z;
        return *this;
}    

inline vec &vec::operator*=(const double &d)
{ 
        _x*=d; _y*=d; _z*=d;
        return *this;
}    

inline vec &vec::operator/=(const double &d)
{ 
        _x/=d; _y/=d; _z/=d;
        return *this;
}    

inline vec operator+(const vec &v1, const vec &v2)
{
    return (vec(v1)+=v2);
}

inline vec operator-(const vec &v1, const vec &v2)
{
    return (vec(v1)-=v2);
}

inline vec operator-(const vec &v1){
    return vec (-v1.getX(), -v1.getY(), -v1.getZ());
}

inline vec operator*(const vec &v1, const double &d)
{
    return (vec(v1)*=d);
}

inline vec operator*(const double &d, const vec &v1)
{
    return (vec(v1)*=d);
}

inline vec operator/(const vec &v1, const double &d)
{
    return (vec(v1)/=d);
}

inline std::ostream &operator<<(std::ostream &out, const vec& v)
{
      out << '<' << v.getX() << ',' << v.getY() << ',' << v.getZ() << '>';
      return out;
}
    
/// dot product
inline double operator*(const vec &v1, const vec &v2)
{
    return v1.getX()*v2.getX() + v1.getY()*v2.getY() + v1.getZ()*v2.getZ();
}

/// cross product
inline vec operator^(const vec &v1, const vec &v2)
{
    return vec(
        v1.getY()*v2.getZ() - v1.getZ()*v2.getY(),
        v1.getZ()*v2.getX() - v1.getX()*v2.getZ(),
        v1.getX()*v2.getY() - v1.getY()*v2.getX()
    );
}

inline double abs(const vec &v)
{
    return sqrt(v*v);
}

inline vec &vec::normalize()
{ 
    return ((*this)*=1./abs(*this));
}



#endif	/* _vec_H */

