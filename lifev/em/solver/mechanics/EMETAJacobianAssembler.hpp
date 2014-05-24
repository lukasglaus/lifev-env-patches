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

	auto dP = eval(W1, _F(dispETFESpace, disp, 0)) * (_d2I1bardF(dispETFESpace, disp, 0) );
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
    auto df = dF * f0;//dot( grad(phi_j) * outerProduct( fiber0, fiber0 ), grad(phi_i) ) )
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


template< typename Mesh, typename FunctorPtr >
void
computeI4JacobianTerms( const vector_Type& disp,
						ETFESpacePtr_Type<Mesh>  dispETFESpace,
					    const vector_Type& fibers,
					    matrixPtr_Type     jacobianPtr,
						FunctorPtr         W4)
{
	using namespace ExpressionAssembly;
	//
	std::cout << "EMETA - Computing I4 jacobian terms ... \n";

    auto F = _F(dispETFESpace, disp, 0);

    auto f0 = value(dispETFESpace, fibers);
	auto dF = _dF;
    auto f = F * f0;
//    auto fxf0 = outerProduct (f, f0);
    auto df = dF * f0;
    auto dfxf0 = dF * outerProduct (f0, f0);
//    auto H = value(activationETFESpace, activation);
//    auto Wa = eval(W, H);
    auto W4f = eval(W4, f);
    auto dP = value(2.0) * W4f *  dfxf0;
	integrate ( elements ( dispETFESpace->mesh() ) ,
			    quadRuleTetra4pt,
				dispETFESpace,
				dispETFESpace,
				dot ( dP , grad (phi_i) )
				) >> jacobianPtr;
}


template< typename Mesh, typename FunctorPtr >
void
computeI4JacobianTermsSecondDerivative( const vector_Type& disp,
						ETFESpacePtr_Type<Mesh>  dispETFESpace,
					    const vector_Type& fibers,
					    matrixPtr_Type     jacobianPtr,
						FunctorPtr         W4)
{
	using namespace ExpressionAssembly;
	//
	std::cout << "EMETA - Computing I4 jacobian terms second derivative ... \n";

    auto F = _F(dispETFESpace, disp, 0);

    auto f0 = value(dispETFESpace, fibers);
	auto dF = _dF;
    auto f = F * f0;
    auto fxf0 = outerProduct (f, f0);
    auto df = dF * f0;
    auto dfxf0 = outerProduct (df, f0);
//    auto H = value(activationETFESpace, activation);
//    auto Wa = eval(W, H);
    auto W4f = eval(W4, f);
    auto dP = value(2.0) * W4f * dot( fxf0, dF ) * fxf0;
	integrate ( elements ( dispETFESpace->mesh() ) ,
			    quadRuleTetra4pt,
				dispETFESpace,
				dispETFESpace,
				dot ( dP , grad (phi_i) )
				) >> jacobianPtr;
}


template< typename Mesh, typename FunctorPtr >
void
computeI8JacobianTerms( const vector_Type& disp,
		                boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > >  dispETFESpace,
					    const vector_Type& fibers,
					    const vector_Type& sheets,
					    matrixPtr_Type     jacobianPtr,
					    FunctorPtr         W8)
{
	using namespace ExpressionAssembly;
//    auto F = _F(dispETFESpace, disp, 0);
//    auto dF = _dF;
//    auto f0 = value(dispETFESpace, fibers);
//    auto f = F * f0;
//    auto s0 = value(dispETFESpace, sheets);
//    auto s = F * s0;
//    auto df = dF * f0;
//    auto ds = dF * s0;
    auto I8bar = _I8fsbar(dispETFESpace, disp, 0, fibers, sheets);

    auto P = eval (W8, I8bar )
    	   * _dI8fsbardF(dispETFESpace, disp, 0, fibers, sheets);

//    auto fxs = _dI8fsbardF(dispETFESpace, disp, 0, fibers, sheets);
//    auto fxs = _dI8fsdF(dispETFESpace, disp, 0, fibers, sheets);
    auto W = eval (W8, I8bar );

	std::cout << "EMETA - Computing I8 jacobian terms ... \n";
	integrate ( elements ( dispETFESpace->mesh() ) ,
			    quadRuleTetra4pt,
				dispETFESpace,
				dispETFESpace,
				dot ( W * _d2I8fsbardF(dispETFESpace, disp, 0, fibers, sheets) , grad (phi_i) )
				) >> jacobianPtr;
}


template< typename Mesh, typename FunctorPtr >
void
computeI8JacobianTermsSecondDerivative( const vector_Type& disp,
										boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > >  dispETFESpace,
										const vector_Type& fibers,
										const vector_Type& sheets,
										matrixPtr_Type     jacobianPtr,
										FunctorPtr         dW8)
{
	using namespace ExpressionAssembly;
//    auto F = _F(dispETFESpace, disp, 0);
//    auto dF = _dF;
//    auto f0 = value(dispETFESpace, fibers);
//    auto f = F * f0;
//    auto s0 = value(dispETFESpace, sheets);
//    auto s = F * s0;
//    auto df = dF * f0;
//    auto ds = dF * s0;
//    auto dfxs = dot( dF, ( outerProduct(f, s0) + outerProduct(s, f0) ) );
//    auto W = eval(W8, f, s);
//
//    auto FmT = minusT(F);
//    auto dFmT = value(-1.0) * FmT * transpose(dF) * FmT;
//    auto I8 = dot(f, s);
//    auto dFdev8fs = (outerProduct (dF * s0, f0) + outerProduct (dF * f0, s0) )
//                    - value (2. / 3) * ( (dot (F * s0, dF * f0) + dot (F * f0, dF * s0) ) *FmT + I8 * dFmT);

//    auto d2I8fseiso = value (2.0) * Jm23 * (dFdev8fs - value (2. / 3) *dot (FmT, dF) *Fdev8fs) );

    auto I8bar = _I8fsbar(dispETFESpace, disp, 0, fibers, sheets);

	auto dP = eval (dW8, I8bar)
			* (_dI8fsbardF(dispETFESpace, disp, 0, fibers, sheets) ) * (_dI8fsbar(dispETFESpace, disp, 0, fibers, sheets) );

    std::cout << "EMETA - Computing I8 jacobian terms second derivative ... \n";
	integrate ( elements ( dispETFESpace->mesh() ) ,
			    quadRuleTetra4pt,
				dispETFESpace,
				dispETFESpace,
				dot (  dP, grad (phi_i) )
				) >> jacobianPtr;
}


}//EMAssembler

}//LifeV

#endif /* EMETAASSMEBLER_HPP_ */
