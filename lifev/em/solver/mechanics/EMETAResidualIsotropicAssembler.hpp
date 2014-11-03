/*
 * EMETAAssmebler.hpp
 *
 *  Created on: 29/apr/2014
 *      Author: srossi
 */

#ifndef EMETARESIDUALISOTROPICASSMEBLER_HPP_
#define EMETARESIDUALISOTROPICASSMEBLER_HPP_

#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/array/MatrixEpetra.hpp>

//ET include for assemblings
#include <lifev/eta/fem/ETFESpace.hpp>
#include <lifev/eta/expression/Integrate.hpp>

#include <lifev/em/solver/mechanics/EMMechanicalExpressions.hpp>

#include <lifev/em/util/EMUtility.hpp>
#include <lifev/em/solver/EMETAFunctors.hpp>

namespace LifeV
{

typedef VectorEpetra           vector_Type;
typedef boost::shared_ptr<vector_Type>         vectorPtr_Type;

typedef MatrixEpetra<Real>           matrix_Type;
typedef boost::shared_ptr<matrix_Type>         matrixPtr_Type;


namespace EMAssembler
{



template< typename Mesh, typename FunctorPtr >
void
computeLinearizedVolumetricResidualTerms ( const vector_Type& disp,
                                           boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > >  dispETFESpace,
                                           vectorPtr_Type           residualVectorPtr,
                                           FunctorPtr               W)
{
    using namespace ExpressionAssembly;
    //
    auto I = value (EMUtility::identity() );
    auto GradU = _Grad_u (dispETFESpace, disp, 0);
    auto GradUT = transpose (GradU);
    auto P = eval (W, I) * (GradU + GradUT) * value (1. / 2.);

	if(disp.comm().MyPID() == 0)
    std::cout << "EMETA - Computing linear volumetric residual terms ... \n";
    integrate ( elements ( dispETFESpace->mesh() ) ,
                quadRuleTetra4pt,
                dispETFESpace,
                dot ( P  , grad (phi_i) )
              ) >> residualVectorPtr;
}

template< typename Mesh, typename FunctorPtr >
void
computeLinearizedDeviatoricResidualTerms ( const vector_Type& disp,
                                           boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > >  dispETFESpace,
                                           vectorPtr_Type           residualVectorPtr,
                                           FunctorPtr               W)
{
    using namespace ExpressionAssembly;
    //
    auto I = value (EMUtility::identity() );
    auto GradU = _Grad_u (dispETFESpace, disp, 0);
    auto GradUT = transpose (GradU);
    auto P = eval (W, I) * trace (GradU + GradUT) * value (1. / 2.) * I;

	if(disp.comm().MyPID() == 0)
    std::cout << "EMETA - Computing linear deviatoric residual terms ... \n";
    integrate ( elements ( dispETFESpace->mesh() ) ,
                quadRuleTetra4pt,
                dispETFESpace,
                dot ( P  , grad (phi_i) )
              ) >> residualVectorPtr;
}

template< typename Mesh, typename FunctorPtr >
void
computeI1ResidualTerms ( const vector_Type& disp,
                         boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > >  dispETFESpace,
                         vectorPtr_Type           residualVectorPtr,
                         FunctorPtr                  W1)
{
    using namespace ExpressionAssembly;

	if(disp.comm().MyPID() == 0)
    std::cout << "EMETA - Computing I1 residual terms ... \n";
    integrate ( elements ( dispETFESpace->mesh() ) ,
                quadRuleTetra4pt,
                dispETFESpace,
                dot ( eval (W1, _F (dispETFESpace, disp, 0) ) * _dI1bar (dispETFESpace, disp, 0), grad (phi_i) )
              ) >> residualVectorPtr;
}


template< typename Mesh, typename FunctorPtr >
void
computeI2ResidualTerms ( const vector_Type& disp,
                         boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > >  dispETFESpace,
                         vectorPtr_Type           residualVectorPtr,
                         FunctorPtr                  W2)
{
    using namespace ExpressionAssembly;

	if(disp.comm().MyPID() == 0)
    std::cout << "EMETA - Computing I2 residual terms ... \n";
    integrate ( elements ( dispETFESpace->mesh() ) ,
                quadRuleTetra4pt,
                dispETFESpace,
                dot ( eval (W2, _F (dispETFESpace, disp, 0) ) * _dI2bar (dispETFESpace, disp, 0), grad (phi_i) )
              ) >> residualVectorPtr;
}


template <typename Mesh, typename FunctorPtr >
void
computeVolumetricResidualTerms ( const vector_Type& disp,
                                 boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > >  dispETFESpace,
                                 vectorPtr_Type           residualVectorPtr,
                                 FunctorPtr                  Wvol)
{
    //
	if(disp.comm().MyPID() == 0)
		std::cout << "EMETA - Computing Volumetric residual terms ... \n";

    using namespace ExpressionAssembly;


    integrate ( elements ( dispETFESpace->mesh() ) ,
                quadRuleTetra4pt,
                dispETFESpace,
                dot ( eval (Wvol, _F (dispETFESpace, disp, 0) ) * _dJ (dispETFESpace, disp, 0), grad (phi_i) )
              ) >> residualVectorPtr;

}

template< typename Mesh, typename FunctorPtr >
void
computeI1ResidualTerms ( const vector_Type& disp,
                         boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > >  dispETFESpace,
                         const vector_Type& fibers,
                         const vector_Type& sheets,
                         vectorPtr_Type     residualVectorPtr,
                         FunctorPtr         W1)
{
    using namespace ExpressionAssembly;

    auto f_0 = _v0 (dispETFESpace, fibers);
    auto s_0 = _v0 (dispETFESpace, sheets);

    boost::shared_ptr<orthonormalizeFibers> normalize0 (new orthonormalizeFibers);
    boost::shared_ptr<orthonormalizeFibers> normalize1 (new orthonormalizeFibers (1) );
    auto f0 = eval (normalize0, f_0);

    auto s_00 = s_0 - dot (f0, s_0) * s_0;

    auto s0 = eval (normalize1, s_00);

	if(disp.comm().MyPID() == 0)
    std::cout << "EMETA - Computing I1 residual terms ... \n";
    integrate ( elements ( dispETFESpace->mesh() ) ,
                quadRuleTetra4pt,
                dispETFESpace,
                dot ( eval (W1, _F (dispETFESpace, disp, 0), f0, s0) * _dI1bar (dispETFESpace, disp, 0), grad (phi_i) )
              ) >> residualVectorPtr;
}


}//EMAssembler

}//LifeV

#endif /* EMETAASSMEBLER_HPP_ */
