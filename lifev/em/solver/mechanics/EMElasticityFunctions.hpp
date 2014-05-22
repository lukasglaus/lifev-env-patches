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


namespace LifeV
{

namespace Elasticity
{

Real J(const LifeV::MatrixSmall<3, 3>& F)
{
	return F.determinant();
}

Real Jm23(const LifeV::MatrixSmall<3, 3>& F)
{
	return std::pow(F.determinant(), -2.0/3.0);
}

Real I1(const LifeV::MatrixSmall<3, 3>& F)
{
	return F.dot(F);
}

Real I1bar(const LifeV::MatrixSmall<3, 3>& F)
{
	return std::pow(J(F), -2.0/3.0) * I1(F);
}

Real I4(const LifeV::MatrixSmall<3, 3>& F, const LifeV::VectorSmall<3>& f0)
{
	auto f = F * f0;
	return f.dot(f);
}

Real I4(const LifeV::VectorSmall<3>& f)
{
	return f.dot(f);
}

Real RegularizedHeaviside(const LifeV::Real x)
{
	return 0.5 * ( 1 + std::tanh(100.0 * x) );
}

}// Elasticity

}//LifeV

#endif /* EMELASTICITYFUNCTIONS_HPP_ */
