/*
 * EMETAAssmebler.hpp
 *
 *  Created on: 29/apr/2014
 *      Author: srossi
 */

#ifndef EMETAJACOBIANASSMEBLER_HPP_
#define EMETAJACOBIANASSMEBLER_HPP_


#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/array/MatrixEpetra.hpp>

//ET include for assemblings
#include <lifev/eta/fem/ETFESpace.hpp>
#include <lifev/eta/expression/Integrate.hpp>

#include <lifev/em/solver/mechanics/EMMechanicalExpressions.hpp>
//#include <lifev/em/solver/mechanics/materials/EMMaterialFunctions.hpp>

#include <lifev/em/util/EMUtility.hpp>
//#include <boost/typeof/typeof.hpp>

namespace LifeV
{

typedef VectorEpetra           vector_Type;
typedef boost::shared_ptr<vector_Type>         vectorPtr_Type;

typedef MatrixEpetra<Real>           matrix_Type;
typedef boost::shared_ptr<matrix_Type>         matrixPtr_Type;

//typedef EMMaterialFunctions  materialFunctions_Type;
//typedef boost::shared_ptr<materialFunctions_Type> materialFunctionsPtr_Type;

namespace EMAssembler
{

template <class Mesh> using ETFESpacePtr_Type = boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > >;
template <class Mesh> using scalarETFESpacePtr_Type = boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 1 > >;
//template <class Mesh> using FESpacePtr_Type = boost::shared_ptr< FESpace< Mesh, MapEpetra >  >;


template< typename Mesh, typename FunctorPtr >
void
computeLinearizedVolumetricJacobianTerms( const vector_Type& disp,
								ETFESpacePtr_Type<Mesh>  dispETFESpace,
								matrixPtr_Type           jacobianPtr,
								FunctorPtr               W)
{
	using namespace ExpressionAssembly;
	//
	auto I = value(EMUtility::identity());
	auto dF = _dF;
	auto dFT = transpose(dF);
	auto dP = eval(W, I ) * (dF + dFT) * value(1./2.);

	std::cout << "Computing linear volumetric Jacobian terms ... \n";
	integrate ( elements ( dispETFESpace->mesh() ) ,
			quadRuleTetra4pt,
				dispETFESpace,
				dispETFESpace,
				dot ( dP  , grad (phi_i) )
				) >> jacobianPtr;
}

template< typename Mesh, typename FunctorPtr >
void
computeLinearizedDeviatoricJacobianTerms( const vector_Type& disp,
								ETFESpacePtr_Type<Mesh>  dispETFESpace,
								matrixPtr_Type           jacobianPtr,
								FunctorPtr               W)
{
	using namespace ExpressionAssembly;
	//
	auto I = value(EMUtility::identity());
	auto dF = _dF;
	auto dFT = transpose(dF);
	auto dP = eval(W, I) * trace(dF + dFT) * value(1./2.) * I;

	std::cout << "Computing linear deviatoric Jacobian terms ... \n";
	integrate ( elements ( dispETFESpace->mesh() ) ,
			quadRuleTetra4pt,
				dispETFESpace,
				dispETFESpace,
				dot ( dP  , grad (phi_i) )
				) >> jacobianPtr;
}

template <typename Mesh, typename FunctorPtr>
void
computeI1JacobianTermsSecondDerivative( const vector_Type& disp,
										ETFESpacePtr_Type<Mesh>  dispETFESpace,
										matrixPtr_Type           jacobianPtr,
										FunctorPtr                 dW1)
{
	using namespace ExpressionAssembly;
	//
	std::cout << "Computing I1 jacobian terms with second derivative of the energy ... \n";

	auto dP = eval(dW1, _F(dispETFESpace, disp, 0)) * (_dI1bardF(dispETFESpace, disp, 0)) * (_dI1bar(dispETFESpace, disp, 0));
	integrate ( elements ( dispETFESpace->mesh() ) ,
			quadRuleTetra4pt,
				dispETFESpace,
				dispETFESpace,
				dot ( dP , grad (phi_i) )
			  ) >> jacobianPtr;
}

template <typename Mesh, typename FunctorPtr >
void
computeI1JacobianTerms( const vector_Type& disp,
                        ETFESpacePtr_Type<Mesh>  dispETFESpace,
                        matrixPtr_Type           jacobianPtr,
                        FunctorPtr                 W1)
{
	using namespace ExpressionAssembly;

	//
	std::cout << "Computing I1 jacobian terms  ... \n";

	auto dP = eval(W1, _F(dispETFESpace, disp, 0)) * (_d2I1bardF(dispETFESpace, disp, 0));
	integrate ( elements ( dispETFESpace->mesh() ) ,
			quadRuleTetra4pt,
				dispETFESpace,
				dispETFESpace,
				dot ( dP , grad (phi_i) )
			  ) >> jacobianPtr;
}



template <typename Mesh, typename FunctorPtr >
void
computeVolumetricJacobianTermsSecondDerivative( const vector_Type& disp,
								ETFESpacePtr_Type<Mesh>  dispETFESpace,
								matrixPtr_Type           jacobianPtr,
								FunctorPtr                 dWvol)
{
	{
		using namespace ExpressionAssembly;

		//
		std::cout << "Computing Volumetric jacobian terms with second derivative of the energy ... \n";

		auto dP = eval(dWvol, _F(dispETFESpace, disp, 0) ) * (_dJdF(dispETFESpace, disp, 0)) * (_dJ(dispETFESpace, disp, 0));
		integrate ( elements ( dispETFESpace->mesh() ) ,
				quadRuleTetra4pt,
					dispETFESpace,
					dispETFESpace,
					dot ( dP , grad (phi_i) )
				  ) >> jacobianPtr;

	}
}


template <typename Mesh, typename FunctorPtr >
void
computeVolumetricJacobianTerms( const vector_Type& disp,
								ETFESpacePtr_Type<Mesh>  dispETFESpace,
								matrixPtr_Type           jacobianPtr,
								FunctorPtr                Wvol)
{
	{
		using namespace ExpressionAssembly;
		//
		std::cout << "Computing Volumetric jacobian terms  ... \n";

//		BOOST_AUTO_TPL (dP, eval(Wvol, _F) * (_d2JdF) );
		integrate ( elements ( dispETFESpace->mesh() ) ,
				quadRuleTetra4pt,
					dispETFESpace,
					dispETFESpace,
					dot ( eval(Wvol, _F(dispETFESpace, disp, 0)) * (_d2JdF(dispETFESpace, disp, 0)), grad (phi_i) )
				  ) >> jacobianPtr;

	}
}


template <typename Mesh, typename FunctorPtr >
void
computeFiberActiveStressJacobianTerms( const vector_Type& disp,
								ETFESpacePtr_Type<Mesh>  dispETFESpace,
						        const vector_Type& fibers,
						        const vector_Type& sheets,
                          const vector_Type& activation,
								scalarETFESpacePtr_Type<Mesh>  activationETFESpace,
								matrixPtr_Type           jacobianPtr,
								FunctorPtr               W)
{
	//
	std::cout << "Computing Fibers Active Stress Jacobian terms ... \n";

	using namespace ExpressionAssembly;

	auto dF = _dF;
	auto I = value(EMUtility::identity());

    auto F = _F(dispETFESpace, disp, 0);
    auto f0 = value(dispETFESpace, fibers);
    auto df = grad(phi_i) * f0;//dot( grad(phi_j) * outerProduct( fiber0, fiber0 ), grad(phi_i) ) )
    auto dfxf0 = outerProduct (df, f0);
    auto H = value(activationETFESpace, activation);
    auto Wa = eval(W, H);
    auto Wm = eval(W, I);
//    auto dP = Wa * Wm * dfxf0;
    auto dP = Wa * Wm * dfxf0;
	integrate ( elements ( dispETFESpace->mesh() ) ,
			quadRuleTetra4pt,
			dispETFESpace,
			dispETFESpace,
			dot ( dP, grad (phi_i) )
		  ) >> jacobianPtr;

}



}//EMAssembler

}//LifeV

#endif /* EMETAASSMEBLER_HPP_ */
