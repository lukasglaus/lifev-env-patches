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
template <class Mesh> using FESpacePtr_Type = boost::shared_ptr< FESpace< Mesh, MapEpetra >  >;


template< typename Mesh, typename FunctorPtr >
void
computeLinearizedVolumetricResidualTerms( const vector_Type& disp,
								ETFESpacePtr_Type<Mesh>  dispETFESpace,
								FESpacePtr_Type<Mesh>    dispFESpace,
								vectorPtr_Type           residualVectorPtr,
								FunctorPtr               W)
{
	using namespace ExpressionAssembly;
	//
	auto GradU = _Grad_u(dispETFESpace, disp, 0);
	auto GradUT = transpose(GradU);
	auto P = eval(W, value(IdMatrix)) * (GradU + GradUT) * value(1./2.);

	std::cout << "Computing linear volumetric residual terms ... \n";
	integrate ( elements ( dispETFESpace->mesh() ) ,
				dispFESpace->qr(),
				dispETFESpace,
				dot ( P  , grad (phi_i) )
				) >> residualVectorPtr;
}

template< typename Mesh, typename FunctorPtr >
void
computeLinearizedDeviatoricResidualTerms( const vector_Type& disp,
								ETFESpacePtr_Type<Mesh>  dispETFESpace,
								FESpacePtr_Type<Mesh>    dispFESpace,
								vectorPtr_Type           residualVectorPtr,
								FunctorPtr               W)
{
	using namespace ExpressionAssembly;
	//
	auto GradU = _Grad_u(dispETFESpace, disp, 0);
	auto GradUT = transpose(GradU);
	auto P = eval(W, value(IdMatrix)) * trace(GradU + GradUT) * value(1./2.) * IdMatrix;

	std::cout << "Computing linear deviatoric residual terms ... \n";
	integrate ( elements ( dispETFESpace->mesh() ) ,
				dispFESpace->qr(),
				dispETFESpace,
				dot ( P  , grad (phi_i) )
				) >> residualVectorPtr;
}

template< typename Mesh, typename FunctorPtr >
void
computeI1ResidualTerms( const vector_Type& disp,
		                ETFESpacePtr_Type<Mesh>  dispETFESpace,
		                FESpacePtr_Type<Mesh>    dispFESpace,
		                vectorPtr_Type           residualVectorPtr,
                        FunctorPtr                  W1)
{
	using namespace ExpressionAssembly;
	//
	std::cout << "Computing I1 residual terms ... \n";
	integrate ( elements ( dispETFESpace->mesh() ) ,
				dispFESpace->qr(),
				dispETFESpace,
				dot ( eval (W1, _F(IdMatrix, dispETFESpace, disp, 0)) * _dI1bar(IdMatrix, dispETFESpace, disp, 0), grad (phi_i) )
				) >> residualVectorPtr;
}


template <typename Mesh, typename FunctorPtr >
void
computeVolumetricResidualTerms( const vector_Type& disp,
								ETFESpacePtr_Type<Mesh>  dispETFESpace,
								FESpacePtr_Type<Mesh>    dispFESpace,
								vectorPtr_Type           residualVectorPtr,
								FunctorPtr                  Wvol)
{
	//
	std::cout << "Computing Volumetric residual terms ... \n";

	using namespace ExpressionAssembly;


	integrate ( elements ( dispETFESpace->mesh() ) ,
			dispFESpace->qr(),
			dispETFESpace,
			dot ( eval(Wvol, _F(IdMatrix, dispETFESpace, disp, 0)) * _dJ(IdMatrix, dispETFESpace, disp, 0), grad (phi_i) )
		  ) >> residualVectorPtr;

}



}//EMAssembler

}//LifeV

#endif /* EMETAASSMEBLER_HPP_ */
