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

}

}//LifeV

#endif /* EMELASTICITYFUNCTIONS_HPP_ */
