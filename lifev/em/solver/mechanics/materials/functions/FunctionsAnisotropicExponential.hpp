/*
 * EMMaterialFunctions.hpp
 *
 *  Created on: 28/apr/2014
 *      Author: srossi
 */

#ifndef FUNCTIONSAnisotropicExponential_HPP_
#define FUNCTIONSAnisotropicExponential_HPP_

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
//  ISOTROPIC EXPONENTIAL FUNCTIONS
////////////////////////////////////////////////////////////////////////

template <class Mesh>
class AnisotropicExponential : public virtual EMMaterialFunctions<Mesh> {
public:
	typedef typename MaterialFunctions::EMMaterialFunctions<Mesh>::return_Type return_Type;

    virtual return_Type operator() (const VectorSmall<3>& f)
    {
    	using namespace ExpressionAssembly;
    	LifeV::Real I4 = f.dot(f);
    	return M_a * ( I4 - 1.0 ) * Elasticity::Heaviside(I4-1);
//    	return M_a * ( I4 - 1.0 ) * std::exp( M_b * ( I4 - 1 ) * ( I4 - 1 ) ); // * Elasticity::RegularizedHeaviside(I4-1);
    }

    AnisotropicExponential(Real a = 1853.50, Real b = 15.972, bool useFibers = true) : M_a(a),
    		                                                                          M_b(b),
    		                                                                          M_useFibers(useFibers) {} // 0.33 KPa
    AnisotropicExponential(const AnisotropicExponential& AnisotropicExponential)
    {
    	M_a = AnisotropicExponential.M_a;
    	M_b = AnisotropicExponential.M_b;
    	M_useFibers = AnisotropicExponential.M_useFibers;
    }
    virtual ~AnisotropicExponential() {}

	inline virtual void computeJacobian( const vector_Type& disp,
			boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > >  dispETFESpace,
								         const vector_Type& fibers,
								         const vector_Type& sheets,
								         matrixPtr_Type     jacobianPtr)
	{
		EMAssembler::computeI4JacobianTerms(disp,dispETFESpace, fibers, jacobianPtr, this->getMe());
	}

	inline virtual void computeResidual( const vector_Type& disp,
			boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > >  dispETFESpace,
							         const vector_Type& fibers,
							         const vector_Type& sheets,
							         vectorPtr_Type     residualVectorPtr)
	{
		EMAssembler::computeI4ResidualTerms(disp,dispETFESpace, fibers, residualVectorPtr, this->getMe());
	}

private:
    Real M_a;
    Real M_b;
    bool M_useFibers; //use fiber vector if true,  else uses sheets

};

template <class Mesh>
class dAnisotropicExponential : public virtual EMMaterialFunctions<Mesh> {
public:
	typedef typename MaterialFunctions::EMMaterialFunctions<Mesh>::return_Type return_Type;

    virtual return_Type operator() (const VectorSmall<3>& f)
    {
    	using namespace ExpressionAssembly;
    	LifeV::Real I4 = f.dot(f);
    	return M_a * Elasticity::Heaviside(I4-1);
  //  	return 2.0 * M_a * M_b * ( I4 - 1.0 ) * ( I4 - 1.0 ) * std::exp( M_b * ( I4 - 1 ) * ( I4 - 1 ) ); // * Elasticity::RegularizedHeaviside(I4-1);
    }

//    dAnisotropicExponential() : M_a(3330), M_b(9.242) {} // 0.33 KPa
    dAnisotropicExponential(Real a = 1853.50, Real b = 15.972, bool useFibers = true) : M_a(a),
    		                                                                           M_b(b),
    		                                                                           M_useFibers(useFibers) {} // 0.33 KPa
    dAnisotropicExponential(const dAnisotropicExponential& dAnisotropicExponential)
    {
    	M_a = dAnisotropicExponential.M_a;
    	M_b = dAnisotropicExponential.M_b;
    	M_useFibers = dAnisotropicExponential.M_useFibers;
    }
    virtual ~dAnisotropicExponential() {}


	inline virtual void computeJacobian( const vector_Type& disp,
			boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > > dispETFESpace,
		     const vector_Type& fibers,
		     const vector_Type& sheets,
			 matrixPtr_Type           jacobianPtr)
	{
		EMAssembler::computeI4JacobianTermsSecondDerivative(disp,dispETFESpace, fibers, jacobianPtr, this->getMe());
	}

private:
    Real M_a;
    Real M_b;
    bool M_useFibers;
};


} //EMMaterialFunctions

} //LifeV
#endif /* EMMATERIALFUNCTIONS_HPP_ */
