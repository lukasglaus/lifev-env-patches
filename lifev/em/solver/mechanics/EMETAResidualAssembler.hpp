/*
 * EMETAAssmebler.hpp
 *
 *  Created on: 29/apr/2014
 *      Author: srossi
 */

#ifndef EMETARESIDUALASSMEBLER_HPP_
#define EMETARESIDUALASSMEBLER_HPP_


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

//static MatrixSmall<3, 3> IdMatrix = EMUtility::identity();

//typedef EMMaterialFunctions  materialFunctions_Type;
//typedef boost::shared_ptr<materialFunctions_Type> materialFunctionsPtr_Type;

namespace EMAssembler
{

template <class Mesh> using ETFESpacePtr_Type = boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > >;
template <class Mesh> using scalarETFESpacePtr_Type = boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 1 > >;
//template <class Mesh> using FESpacePtr_Type = boost::shared_ptr< FESpace< Mesh, MapEpetra >  >;


template< typename Mesh, typename FunctorPtr >
void
computeLinearizedVolumetricResidualTerms( const vector_Type& disp,
								ETFESpacePtr_Type<Mesh>  dispETFESpace,
								vectorPtr_Type           residualVectorPtr,
								FunctorPtr               W)
{
	using namespace ExpressionAssembly;
	//
	auto I = value(EMUtility::identity());
	auto GradU = _Grad_u(dispETFESpace, disp, 0);
	auto GradUT = transpose(GradU);
	auto P = eval(W, I) * (GradU + GradUT) * value(1./2.);

	std::cout << "Computing linear volumetric residual terms ... \n";
	integrate ( elements ( dispETFESpace->mesh() ) ,
			quadRuleTetra4pt,
				dispETFESpace,
				dot ( P  , grad (phi_i) )
				) >> residualVectorPtr;
}

template< typename Mesh, typename FunctorPtr >
void
computeLinearizedDeviatoricResidualTerms( const vector_Type& disp,
								ETFESpacePtr_Type<Mesh>  dispETFESpace,
								vectorPtr_Type           residualVectorPtr,
								FunctorPtr               W)
{
	using namespace ExpressionAssembly;
	//
	auto I = value(EMUtility::identity());
	auto GradU = _Grad_u(dispETFESpace, disp, 0);
	auto GradUT = transpose(GradU);
	auto P = eval(W, I) * trace(GradU + GradUT) * value(1./2.) * I;

	std::cout << "Computing linear deviatoric residual terms ... \n";
	integrate ( elements ( dispETFESpace->mesh() ) ,
			    quadRuleTetra4pt,
				dispETFESpace,
				dot ( P  , grad (phi_i) )
				) >> residualVectorPtr;
}

template< typename Mesh, typename FunctorPtr >
void
computeI1ResidualTerms( const vector_Type& disp,
		                ETFESpacePtr_Type<Mesh>  dispETFESpace,
		                vectorPtr_Type           residualVectorPtr,
                        FunctorPtr                  W1)
{
	using namespace ExpressionAssembly;
	//
//	auto I = value(EMUtility::identity());
//	//auto u = value(dispETFESpace, disp, 0);
//	auto F = I + grad(dispETFESpace, disp, 0);
//	auto J = det(F);
//	auto Jm23 = pow(J, -2.0/3);
//	auto I1 = dot(F, F);
//	auto FinvT = minusT(F);
//	auto P = 4960 * Jm23 * (F - I1 * value(1.0/3.0) * FinvT);
//	std::cout << "Computing I1 residual terms ... \n";
	integrate ( elements ( dispETFESpace->mesh() ) ,
			quadRuleTetra4pt,
				dispETFESpace,
		//		dot(P, grad (phi_i) )
				dot ( eval (W1, _F(dispETFESpace, disp, 0)) * _dI1bar(dispETFESpace, disp, 0), grad (phi_i) )
				) >> residualVectorPtr;
}


template <typename Mesh, typename FunctorPtr >
void
computeVolumetricResidualTerms( const vector_Type& disp,
								ETFESpacePtr_Type<Mesh>  dispETFESpace,
								vectorPtr_Type           residualVectorPtr,
								FunctorPtr                  Wvol)
{
	//
	std::cout << "Computing Volumetric residual terms ... \n";

	using namespace ExpressionAssembly;


	integrate ( elements ( dispETFESpace->mesh() ) ,
			quadRuleTetra4pt,
			dispETFESpace,
			dot ( eval(Wvol, _F(dispETFESpace, disp, 0)) * _dJ(dispETFESpace, disp, 0), grad (phi_i) )
		  ) >> residualVectorPtr;

}


template <typename Mesh, typename FunctorPtr >
void
computeFiberActiveStressResidualTerms( const vector_Type& disp,
								ETFESpacePtr_Type<Mesh>  dispETFESpace,
							  const vector_Type& fibers,
						        const vector_Type& sheets,
							  const vector_Type& activation,
								scalarETFESpacePtr_Type<Mesh>  activationETFESpace,
								vectorPtr_Type           residualVectorPtr,
								FunctorPtr               W)
{
	//
	std::cout << "Computing Fibers Active Stress residual terms ... \n";

	using namespace ExpressionAssembly;

    auto F = _F(dispETFESpace, disp, 0);
	auto I = value(EMUtility::identity());

    auto f0 = value(dispETFESpace, fibers);
    auto f = F * f0;
    auto fxf0 = outerProduct (f, f0);
    auto H = value(activationETFESpace, activation);
    auto Wa = eval(W, H);
    auto Wm = eval(W, I);
    auto P = Wa /* Wm */ *  fxf0;

    integrate ( elements ( dispETFESpace->mesh() ) ,
			    quadRuleTetra4pt,
			    dispETFESpace,
				dot ( Wa* F * outerProduct(f0, f0), grad (phi_i) )
			    //			    dot ( P, grad (phi_i) )
		      ) >> residualVectorPtr;

}


template< typename Mesh, typename FunctorPtr >
void
computeI4ResidualTerms( const vector_Type& disp,
						ETFESpacePtr_Type<Mesh>  dispETFESpace,
					    const vector_Type& fibers,
						vectorPtr_Type     residualVectorPtr,
						FunctorPtr         W4)
{
	using namespace ExpressionAssembly;
	//
    auto F = _F(dispETFESpace, disp, 0);

    auto f0 = value(dispETFESpace, fibers);
    auto f = F * f0;
    auto fxf0 = value(2.0) * outerProduct (f, f0);
//    auto H = value(activationETFESpace, activation);
//    auto Wa = eval(W, H);
//    auto W4 = eval(W4, I);
  //  auto P = Wa /* Wm */ *  fxf0;
	integrate ( elements ( dispETFESpace->mesh() ) ,
			quadRuleTetra4pt,
				dispETFESpace,
		//		dot(P, grad (phi_i) )
				dot ( eval (W4, F) * fxf0 , grad (phi_i) )
				) >> residualVectorPtr;
}


}//EMAssembler

}//LifeV

#endif /* EMETAASSMEBLER_HPP_ */
