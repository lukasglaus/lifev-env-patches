/*
 * EMMaterialFunctions.hpp
 *
 *  Created on: 28/apr/2014
 *      Author: srossi
 */

#ifndef FUNCTIONSVOLUMETRIC_HPP_
#define FUNCTIONSVOLUMETRIC_HPP_

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
//  VOLUMETRIC FUNCTIONS
////////////////////////////////////////////////////////////////////////
template <class Mesh>
class Volumetric : public virtual MaterialFunctions::EMMaterialFunctions<Mesh> {
public:
	typedef typename MaterialFunctions::EMMaterialFunctions<Mesh>::return_Type return_Type;

    virtual return_Type operator() (const MatrixSmall<3, 3>& F)
    {
        Real J = Elasticity::J(F);
        return M_bulk * ( (J-1) + std::log(J) / J );
    }

    Volumetric(Real bulk = 2000.0) : M_bulk(bulk) {}
    Volumetric (const Volumetric& v) { M_bulk = v.M_bulk; }
    virtual ~Volumetric() {}

	inline virtual void computeJacobian( const vector_Type& disp,
								  ETFESpacePtr_Type<Mesh>  dispETFESpace,
								  FESpacePtr_Type<Mesh>    dispFESpace,
								  matrixPtr_Type           jacobianPtr)
	{
		EMAssembler::computeVolumetricJacobianTerms(disp,dispETFESpace,dispFESpace, jacobianPtr, this->getMe());
	}

	inline virtual void computeResidual( const vector_Type& disp,
								  ETFESpacePtr_Type<Mesh>  dispETFESpace,
								  FESpacePtr_Type<Mesh>    dispFESpace,
								  vectorPtr_Type           residualVectorPtr)
	{
		EMAssembler::computeVolumetricResidualTerms(disp,dispETFESpace,dispFESpace, residualVectorPtr, this->getMe());
	}

private:
    Real M_bulk;
};

template <class Mesh>
class dVolumetric : public virtual EMMaterialFunctions<Mesh> {
public:
	typedef typename MaterialFunctions::EMMaterialFunctions<Mesh>::return_Type return_Type;

    virtual return_Type operator() (const MatrixSmall<3, 3>& F)
    {
        Real J = Elasticity::J(F);
    	return M_bulk  / J / J * ( 1 + J * J + std::log(J)) ;
    }

    dVolumetric(Real bulk = 2000.0) : M_bulk(bulk) {}
    dVolumetric (const dVolumetric& v) { M_bulk = v.M_bulk; }
    virtual ~dVolumetric() {}

	inline virtual void computeJacobian( const vector_Type& disp,
								  ETFESpacePtr_Type<Mesh>  dispETFESpace,
								  FESpacePtr_Type<Mesh>    dispFESpace,
								  matrixPtr_Type           jacobianPtr)
	{
		EMAssembler::computeVolumetricJacobianTermsSecondDerivative(disp,dispETFESpace,dispFESpace, jacobianPtr, this->getMe());
	}

private:
    Real M_bulk;
};



} //EMMaterialFunctions

} //LifeV
#endif /* EMMATERIALFUNCTIONS_HPP_ */
