/*
 * EMElasticityFunctions.hpp
 *
 *  Created on: 28/apr/2014
 *      Author: srossi
 */

#ifndef EMELASTICITYFUNCTIONS_HPP_
#define EMELASTICITYFUNCTIONS_HPP_


#include <lifev/core/LifeV.hpp>
#include <lifev/core/array/MatrixSmall.hpp>
#include <lifev/core/array/VectorSmall.hpp>


namespace LifeV
{

namespace Elasticity
{

Real J (const LifeV::MatrixSmall<3, 3>& F)
{
    return F.determinant();
}

LifeV::MatrixSmall<3, 3> H (const LifeV::MatrixSmall<3, 3>& F)
{
    return F.cofactor();
}

Real Jm23 (const LifeV::MatrixSmall<3, 3>& F)
{
    return std::pow (F.determinant(), -2.0 / 3.0);
}

Real Jm43 (const LifeV::MatrixSmall<3, 3>& F)
{
    return std::pow (F.determinant(), -4.0 / 3.0);
}

Real Jm53 (const LifeV::MatrixSmall<3, 3>& F)
{
    return std::pow (F.determinant(), -5.0 / 3.0);
}

Real Jm83 (const LifeV::MatrixSmall<3, 3>& F)
{
    return std::pow (F.determinant(), -8.0 / 3.0);
}

Real I1 (const LifeV::MatrixSmall<3, 3>& F)
{
    return F.dot (F);
}

Real I1bar (const LifeV::MatrixSmall<3, 3>& F)
{
    return std::pow (J (F), -2.0 / 3.0) * I1 (F);
}

Real I2 (const LifeV::MatrixSmall<3, 3>& F)
{
    return 0.5 * ( I1(F) * I1(F) - I1(F.transpose() * F ) );
}

Real I2bar (const LifeV::MatrixSmall<3, 3>& F)
{
    return  std::pow (J (F), 4 / -3.0) * I2 (F);
}

Real I4 (const LifeV::MatrixSmall<3, 3>& F, const LifeV::VectorSmall<3>& f0)
{
    auto f = F * f0;
    return f.dot (f);
}

Real I4 (const LifeV::VectorSmall<3>& f)
{
    return f.dot (f);
}

Real I8 (const LifeV::MatrixSmall<3, 3>& F, const LifeV::VectorSmall<3>& f0, const LifeV::VectorSmall<3>& s0)
{
    auto f = F * f0;
    auto s = F * s0;
    return f.dot (s);
}


Real RegularizedHeaviside (const LifeV::Real x)
{
    return 0.5 * ( 1 + std::tanh (1000.0 * x) );
}

Real dRegularizedHeaviside (const LifeV::Real x)
{
    Real sech = 1.0 / std::cosh (1000.0 * x);
    return 0.5 * 1000.0 *  sech * sech;
}

Real Heaviside (const LifeV::Real x)
{
    if (x > 0)
    {
        return 1.0;
    }
    else
    {
        return 0.0;
    }
}


VectorSmall<3> crossProduct (const VectorSmall<3>& v1, const VectorSmall<3>& v2)
{
    VectorSmall<3> v;
    v[0] = v1[1] * v2[2] - v1[2] * v2[1];
    v[1] = v1[2] * v2[0] - v1[0] * v2[2];
    v[2] = v1[0] * v2[1] - v1[1] * v2[0];
    return v;
}



}// Elasticity

}//LifeV

#endif /* EMELASTICITYFUNCTIONS_HPP_ */
