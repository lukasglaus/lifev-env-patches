/*
 * EMMaterialFunctions.hpp
 *
 *  Created on: 28/apr/2014
 *      Author: srossi
 */

#ifndef FUNCTIONSNEOHOOKEANLINEARIZED_HPP_
#define FUNCTIONSNEOHOOKEANLINEARIZED_HPP_

//#include <lifev/em/solver/mechanics/EMETAAssembler.hpp>
#include <lifev/em/solver/mechanics/materials/functions/EMMaterialFunctions.hpp>

//using namespace LifeV;

//namespace MaterialFunctions
namespace LifeV
{

namespace MaterialFunctions
{


template <class Mesh> using ETFESpacePtr_Type = boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > >;
template <class Mesh> using FESpacePtr_Type = boost::shared_ptr< FESpace< Mesh, MapEpetra >  >;

typedef VectorEpetra           vector_Type;
typedef boost::shared_ptr<vector_Type>         vectorPtr_Type;

typedef MatrixEpetra<Real>           matrix_Type;
typedef boost::shared_ptr<matrix_Type>         matrixPtr_Type;



////////////////////////////////////////////////////////////////////////
//  LINEARIZED NEO HOOKEAN FUNCTIONS
////////////////////////////////////////////////////////////////////////
template <class Mesh>
class LinearizedNeoHookeanVolumetric : public virtual EMMaterialFunctions<Mesh> {
public:
	typedef typename MaterialFunctions::EMMaterialFunctions<Mesh>::return_Type return_Type;

    virtual return_Type operator() (const MatrixSmall<3, 3>& F)
    {
    	return 2 * M_mu;
    }

//    NeoHookeanLinearized() : M_mu(4960), M {} // 0.496 KPa
    LinearizedNeoHookeanVolumetric(Real mu = 4960.) : M_mu(mu) {} // 0.496 KPa
    LinearizedNeoHookeanVolumetric (const LinearizedNeoHookeanVolumetric& neoHookeanLinearized)
    {
    	M_mu = neoHookeanLinearized.M_mu;
    }
    virtual ~LinearizedNeoHookeanVolumetric() {}

	inline virtual void computeJacobian( const vector_Type& disp,
								  ETFESpacePtr_Type<Mesh>  dispETFESpace,
								  FESpacePtr_Type<Mesh>    dispFESpace,
								  matrixPtr_Type           jacobianPtr)
	{
		EMAssembler::computeLinearizedVolumetricJacobianTerms(disp,dispETFESpace,dispFESpace, jacobianPtr, this->getMe());
	}

	inline virtual void computeResidual( const vector_Type& disp,
								  ETFESpacePtr_Type<Mesh>  dispETFESpace,
								  FESpacePtr_Type<Mesh>    dispFESpace,
								  vectorPtr_Type           residualVectorPtr)
	{
		EMAssembler::computeLinearizedVolumetricResidualTerms(disp,dispETFESpace,dispFESpace, residualVectorPtr, this->getMe());
	}

private:
    Real M_mu;
};

template <class Mesh>
class LinearizedNeoHookeanDeviatoric : public virtual EMMaterialFunctions<Mesh> {
public:
	typedef typename MaterialFunctions::EMMaterialFunctions<Mesh>::return_Type return_Type;

    virtual return_Type operator() (const MatrixSmall<3, 3>& F)
    {
    	//using the volumetric function U = (J-1)^2 + (ln J)^2 ->  U''(1) = 4;
    	return 4 * M_bulk - 2.0 / 3.0 * M_mu;
    }

//    NeoHookeanLinearized() : M_mu(4960), M {} // 0.496 KPa
    LinearizedNeoHookeanDeviatoric(Real mu = 4960., Real bulk = 2000.) : M_mu(mu), M_bulk(bulk) {} // 0.496 KPa
    LinearizedNeoHookeanDeviatoric (const LinearizedNeoHookeanDeviatoric& neoHookeanLinearized)
    {
    	M_mu = neoHookeanLinearized.M_mu;
    	M_bulk = neoHookeanLinearized.M_bulk;
    }
    virtual ~LinearizedNeoHookeanDeviatoric() {}

	inline virtual void computeJacobian( const vector_Type& disp,
								  ETFESpacePtr_Type<Mesh>  dispETFESpace,
								  FESpacePtr_Type<Mesh>    dispFESpace,
								  matrixPtr_Type           jacobianPtr)
	{
		EMAssembler::computeLinearizedDeviatoricJacobianTerms(disp,dispETFESpace,dispFESpace, jacobianPtr, this->getMe());
	}

	inline virtual void computeResidual( const vector_Type& disp,
								  ETFESpacePtr_Type<Mesh>  dispETFESpace,
								  FESpacePtr_Type<Mesh>    dispFESpace,
								  vectorPtr_Type           residualVectorPtr)
	{
		EMAssembler::computeLinearizedDeviatoricResidualTerms(disp,dispETFESpace,dispFESpace, residualVectorPtr, this->getMe());
	}

private:
    Real M_mu;
    Real M_bulk;
};

} //EMMaterialFunctions

} //LifeV
#endif /* EMMATERIALFUNCTIONS_HPP_ */
