/*
 * EMMaterialFunctions.hpp
 *
 *  Created on: 28/apr/2014
 *      Author: srossi
 */

#ifndef FUNCTIONSANAnisotropicExponential_HPP_
#define FUNCTIONSANAnisotropicExponential_HPP_

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
    	auto I4 = Elasticity::I4(f);

    	return M_a / 2.0 * ( I4 - 1 ) * std::exp( M_b * ( I4 - 1 ) * ( I4 - 1 ) ) * Elasticity::RegularizedHeaviside(I4-1);
    }

    AnisotropicExponential(Real a = 185350, Real b = 15.972, bool useFibers = true) : M_a(a),
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
	//	EMAssembler::computeI1JacobianTerms(disp,dispETFESpace, jacobianPtr, this->getMe());
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

    virtual return_Type operator() (const MatrixSmall<3, 3>& F)
    {
    	auto I1bar = Elasticity::I1bar(F);
    	return M_a * M_b / 2.0 * std::exp( M_b * ( I1bar - 3 ) );
    }

    dAnisotropicExponential() : M_a(3330), M_b(9.242) {} // 0.33 KPa
    dAnisotropicExponential(Real a, Real b) : M_a(a), M_b(b) {} // 0.33 KPa
    dAnisotropicExponential(const dAnisotropicExponential& dAnisotropicExponential)
    {
    	M_a = dAnisotropicExponential.M_a;
    	M_b = dAnisotropicExponential.M_b;
    }
    virtual ~dAnisotropicExponential() {}


	inline virtual void computeJacobian( const vector_Type& disp,
			boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > > dispETFESpace,
		     const vector_Type& fibers,
		     const vector_Type& sheets,
								  matrixPtr_Type           jacobianPtr)
	{
		EMAssembler::computeI1JacobianTermsSecondDerivative(disp,dispETFESpace, jacobianPtr, this->getMe());
	}
private:
    Real M_a;
    Real M_b;
};


} //EMMaterialFunctions

} //LifeV
#endif /* EMMATERIALFUNCTIONS_HPP_ */
