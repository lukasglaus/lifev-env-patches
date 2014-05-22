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
    	return 0.5 * M_mu;
    }

//    NeoHookean() : M_mu(4960) {} // 0.496 KPa
    NeoHookean(Real mu = 4960) : M_mu(mu) {} // 0.496 KPa
    NeoHookean (const NeoHookean& neoHookean)
    {
    	M_mu = neoHookean.M_mu;
    }
    virtual ~NeoHookean() {}

	inline virtual void computeJacobian( const vector_Type& disp,
								  boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > >  dispETFESpace,
							         const vector_Type& fibers,
							         const vector_Type& sheets,
								  matrixPtr_Type           jacobianPtr)
	{
		EMAssembler::computeI1JacobianTerms(disp,dispETFESpace, jacobianPtr, this->getMe());
	}

	inline virtual void computeResidual( const vector_Type& disp,
								  boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > >  dispETFESpace,
							         const vector_Type& fibers,
							         const vector_Type& sheets,
								  vectorPtr_Type           residualVectorPtr)
	{
		EMAssembler::computeI1ResidualTerms(disp,dispETFESpace, residualVectorPtr, this->getMe());
	}

private:
    Real M_mu;
};


} //EMMaterialFunctions

} //LifeV
#endif /* EMMATERIALFUNCTIONS_HPP_ */
