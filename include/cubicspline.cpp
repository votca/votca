/* 
 * File:   cubicspline.cpp
 * Author: ruehle
 *
 * Created on August 22, 2008, 4:44 PM
 */

#include <cubicspline.h>

void CubicSpline::Fit(ub::vector<double> x, ub::vector<double> y)
{
    ub::matrix<double> A;
    A.resize(_r.size*4, x.size());
}

