/*
 * EMMaterialFunctions.hpp
 *
 *  Created on: 28/apr/2014
 *      Author: srossi
 */

#ifndef FUNCTIONSNEOHOOKEAN_HPP_
#define FUNCTIONSNEOHOOKEAN_HPP_

#include <lifev/em/solver/mechanics/EMElasticityFunctions.hpp>

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
//  NEO HOOKEAN FUNCTIONS
////////////////////////////////////////////////////////////////////////
template <class Mesh>
class NeoHookean : public virtual EMMaterialFunctions<Mesh> {
public:
	typedef typename MaterialFunctions::EMMaterialFunctions<Mesh>::return_Type return_Type;

    virtual return_Type operator() (const MatrixSmall<3, 3>& F)
    {
    	return M_mu;
    }

    NeoHookean() : M_mu(4960) {} // 0.496 KPa
    NeoHookean(Real mu) : M_mu(mu) {} // 0.496 KPa
    NeoHookean (const NeoHookean& neoHookean)
    {
    	M_mu = neoHookean.M_mu;
    }
    virtual ~NeoHookean() {}

	inline virtual void computeJacobian( const vector_Type& disp,
								  ETFESpacePtr_Type<Mesh>  dispETFESpace,
								  FESpacePtr_Type<Mesh>    dispFESpace,
								  matrixPtr_Type           jacobianPtr)
	{
		EMAssembler::computeI1JacobianTerms(disp,dispETFESpace,dispFESpace, jacobianPtr, this->getMe());
	}

	inline virtual void computeResidual( const vector_Type& disp,
								  ETFESpacePtr_Type<Mesh>  dispETFESpace,
								  FESpacePtr_Type<Mesh>    dispFESpace,
								  vectorPtr_Type           residualVectorPtr)
	{
		EMAssembler::computeI1ResidualTerms(disp,dispETFESpace,dispFESpace, residualVectorPtr, this->getMe());
	}

private:
    Real M_mu;
};


} //EMMaterialFunctions

} //LifeV
#endif /* EMMATERIALFUNCTIONS_HPP_ */
