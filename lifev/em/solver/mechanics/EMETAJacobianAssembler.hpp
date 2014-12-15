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
#include <lifev/em/solver/EMETAFunctors.hpp>

//#include <boost/typeof/typeof.hpp>

namespace LifeV
{

typedef VectorEpetra           vector_Type;
typedef boost::shared_ptr<vector_Type>         vectorPtr_Type;

typedef MatrixEpetra<Real>           matrix_Type;
typedef boost::shared_ptr<matrix_Type>         matrixPtr_Type;


namespace EMAssembler
{

template <typename Mesh, typename FunctorPtr >
void
computeFiberActiveStressJacobianTerms ( const vector_Type& disp,
                                        boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > > dispETFESpace,
                                        const vector_Type& fibers,
                                        const vector_Type& sheets,
                                        const vector_Type& activation,
                                        boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 1 > > activationETFESpace,
                                        matrixPtr_Type           jacobianPtr,
                                        FunctorPtr               W)
{
    //
	if(disp.comm().MyPID() == 0)
    std::cout << "Computing Fibers Active Stress Jacobian terms ... \n";

    using namespace ExpressionAssembly;

    auto dF = _dF;
    auto I = value (EMUtility::identity() );

    auto F = _F (dispETFESpace, disp, 0);
    auto f0 = value (dispETFESpace, fibers);
    auto df = dF * f0;//dot( grad(phi_j) * outerProduct( fiber0, fiber0 ), grad(phi_i) ) )
    auto dfxf0 = outerProduct (df, f0);
    auto H = value (activationETFESpace, activation);
    auto Wa = eval (W, H);
    auto Wm = eval (W, I);
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
computeI4JacobianTerms ( const vector_Type& disp,
                         boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > >  dispETFESpace,
                         const vector_Type& fibers,
                         matrixPtr_Type     jacobianPtr,
                         FunctorPtr         W4)
{
    using namespace ExpressionAssembly;

    auto f_0 = _v0 (dispETFESpace, fibers);

    boost::shared_ptr<orthonormalizeFibers> normalize0 (new orthonormalizeFibers);
    auto f0 = eval (normalize0, f_0);

    auto dP = eval (W4, _I4 ( dispETFESpace, disp, 0, f0 ) )
              *  _d2I4dF (dispETFESpace, f0 );

	if(disp.comm().MyPID() == 0)
    std::cout << "EMETA - Computing I4 f jacobian terms ... \n";

	integrate ( elements ( dispETFESpace->mesh() ) ,
                quadRuleTetra4pt,
                dispETFESpace,
                dispETFESpace,
                dot (  dP, grad (phi_i) )
              ) >> jacobianPtr;
}


template< typename Mesh, typename FunctorPtr >
void
computeI4JacobianTermsSecondDerivative ( const vector_Type& disp,
                                         boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > >  dispETFESpace,
                                         const vector_Type& fibers,
                                         matrixPtr_Type     jacobianPtr,
                                         FunctorPtr         dW4)
{
    using namespace ExpressionAssembly;
    //
	if(disp.comm().MyPID() == 0)
    std::cout << "EMETA - Computing I4 f jacobian terms second derivative ... \n";

    auto f_0 = _v0 (dispETFESpace, fibers);

    boost::shared_ptr<orthonormalizeFibers> normalize0 (new orthonormalizeFibers);
    auto f0 = eval (normalize0, f_0);

    auto dP = eval (dW4, _I4 ( dispETFESpace, disp, 0, f0 ) )
              *  _dI4dF ( dispETFESpace, disp, 0, f0 )
              *  _dI4 ( dispETFESpace, disp, 0, f0 );

    integrate ( elements ( dispETFESpace->mesh() ) ,
                quadRuleTetra4pt,
                dispETFESpace,
                dispETFESpace,
                dot ( dP , grad (phi_i) )
              ) >> jacobianPtr;
}


template< typename Mesh, typename FunctorPtr >
void
computeI4JacobianTerms ( const vector_Type& disp,
                         boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > >  dispETFESpace,
                         const vector_Type& fibers,
                         const vector_Type& sheets,
                         matrixPtr_Type     jacobianPtr,
                         FunctorPtr         W4)
{
    using namespace ExpressionAssembly;

    auto f_0 = _v0 (dispETFESpace, fibers);
    auto s_0 = _v0 (dispETFESpace, sheets);

    boost::shared_ptr<orthonormalizeFibers> normalize0 (new orthonormalizeFibers);
    boost::shared_ptr<orthonormalizeFibers> normalize1 (new orthonormalizeFibers (1) );
    auto f0 = eval (normalize0, f_0);

    auto s_00 = s_0 - dot (f0, s_0) * f0;

    auto s0 = eval (normalize1, s_00);

    auto dP = eval (W4, _I4 ( dispETFESpace, disp, 0, s0 ) )
              *  _d2I4dF (dispETFESpace, s0 );

	if(disp.comm().MyPID() == 0)
    std::cout << "EMETA - Computing I4 f jacobian terms ... \n";
    integrate ( elements ( dispETFESpace->mesh() ) ,
                quadRuleTetra4pt,
                dispETFESpace,
                dispETFESpace,
                dot (  dP, grad (phi_i) )
              ) >> jacobianPtr;
}


template< typename Mesh, typename FunctorPtr >
void
computeI4JacobianTermsSecondDerivative ( const vector_Type& disp,
                                         boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > >  dispETFESpace,
                                         const vector_Type& fibers,
                                         const vector_Type& sheets,
                                         matrixPtr_Type     jacobianPtr,
                                         FunctorPtr         dW4)
{
    using namespace ExpressionAssembly;
    //
	if(disp.comm().MyPID() == 0)
    std::cout << "EMETA - Computing I4 s jacobian terms second derivative ... \n";

    auto f_0 = _v0 (dispETFESpace, fibers);
    auto s_0 = _v0 (dispETFESpace, sheets);

    boost::shared_ptr<orthonormalizeFibers> normalize0 (new orthonormalizeFibers);
    boost::shared_ptr<orthonormalizeFibers> normalize1 (new orthonormalizeFibers (1) );
    auto f0 = eval (normalize0, f_0);

    auto s_00 = s_0 - dot (f0, s_0) * s_0;

    auto s0 = eval (normalize1, s_00);


    auto dP = eval (dW4, _I4 ( dispETFESpace, disp, 0, s0 ) )
              *  _dI4dF ( dispETFESpace, disp, 0, s0 )
              *  _dI4 ( dispETFESpace, disp, 0, s0 );

    integrate ( elements ( dispETFESpace->mesh() ) ,
                quadRuleTetra4pt,
                dispETFESpace,
                dispETFESpace,
                dot ( dP , grad (phi_i) )
              ) >> jacobianPtr;
}


//template< typename Mesh, typename FunctorPtr >
//void
//computeI4JacobianTermsFung( const vector_Type& disp,
//      boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > >  dispETFESpace,
//                      const vector_Type& fibers,
//                      const vector_Type& sheets,
//                      matrixPtr_Type     jacobianPtr,
//                      FunctorPtr         W4,
//                      Int Case = 0)
//{
//  using namespace ExpressionAssembly;
//
//  auto f_0 = _v0(dispETFESpace, fibers);
//    auto s_0 = _v0(dispETFESpace, sheets);
//
//    boost::shared_ptr<orthonormalizeFibers> normalize0(new orthonormalizeFibers);
//    boost::shared_ptr<orthonormalizeFibers> normalize1(new orthonormalizeFibers(1));
//    auto f0 = eval(normalize0, f_0);
//
//    auto s_00 = s_0 - dot(f0, s_0) * s_0;
//
//    auto s0 = eval(normalize1, s_00);
//
//    if(Case = 0)
//    {
//  auto dP = eval (W4, _F( dispETFESpace, disp, 0), f0, s0 )
//         *  _d2I4dF(dispETFESpace, f0 );
//
//  std::cout << "EMETA - Computing I4 f jacobian terms (Fung)... \n";
//  integrate ( elements ( dispETFESpace->mesh() ) ,
//              quadRuleTetra4pt,
//              dispETFESpace,
//              dispETFESpace,
//              dot (  dP, grad (phi_i) )
//              ) >> jacobianPtr;
//    }
//    else
//    {
//  auto dP = eval (W4, _F( dispETFESpace, disp, 0), f0, s0 )
//         *  _d2I4dF(dispETFESpace, s0 );
//
//  std::cout << "EMETA - Computing I4 s jacobian terms (Fung)... \n";
//  integrate ( elements ( dispETFESpace->mesh() ) ,
//              quadRuleTetra4pt,
//              dispETFESpace,
//              dispETFESpace,
//              dot (  dP, grad (phi_i) )
//              ) >> jacobianPtr;
//    }
//}

//template< typename Mesh, typename FunctorPtr >
//void
//computeI4JacobianTermsFungSecondDerivative( const vector_Type& disp,
//      boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > >  dispETFESpace,
//                      const vector_Type& fibers,
//                      const vector_Type& sheets,
//                      matrixPtr_Type     jacobianPtr,
//                      FunctorPtr         dW4,
//                      Int Case = 0)
//{
//  using namespace ExpressionAssembly;
//  //
//
//  auto f_0 = _v0(dispETFESpace, fibers);
//    auto s_0 = _v0(dispETFESpace, sheets);
//
//    boost::shared_ptr<orthonormalizeFibers> normalize0(new orthonormalizeFibers);
//    boost::shared_ptr<orthonormalizeFibers> normalize1(new orthonormalizeFibers(1));
//    auto f0 = eval(normalize0, f_0);
//
//    auto s_00 = s_0 - dot(f0, s_0) * s_0;
//
//    auto s0 = eval(normalize1, s_00);
//
//    if(Case == 0)
//    {
//      std::cout << "EMETA - Computing I4 f jacobian terms second derivative (Fung) ... \n";
//
//      auto dP = eval (dW4, _F( dispETFESpace, disp, 0 ), f0, s0 )
//             *  _dI4dF( dispETFESpace, disp, 0, f0 )
//             *  _dI4( dispETFESpace, disp, 0, f0 );
//
//      integrate ( elements ( dispETFESpace->mesh() ) ,
//                  quadRuleTetra4pt,
//                  dispETFESpace,
//                  dispETFESpace,
//                  dot ( dP , grad (phi_i) )
//                  ) >> jacobianPtr;
//    }
//    else
//    {
//      std::cout << "EMETA - Computing I4 s jacobian terms second derivative (Fung) ... \n";
//
//      auto dP = eval (dW4, _F( dispETFESpace, disp, 0), f0, s0 )
//             *  _dI4dF( dispETFESpace, disp, 0, s0 )
//             *  _dI4( dispETFESpace, disp, 0, s0 );
//
//      integrate ( elements ( dispETFESpace->mesh() ) ,
//                  quadRuleTetra4pt,
//                  dispETFESpace,
//                  dispETFESpace,
//                  dot ( dP , grad (phi_i) )
//                  ) >> jacobianPtr;
//    }
//
//}


//template< typename Mesh, typename FunctorPtr >
//void
//computeI4JacobianTermsMixedFungSecondDerivative( const vector_Type& disp,
//      boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > >  dispETFESpace,
//                      const vector_Type& fibers,
//                      const vector_Type& sheets,
//                      matrixPtr_Type     jacobianPtr,
//                      FunctorPtr         dW4,
//                      Int Case = 0)
//{
//  using namespace ExpressionAssembly;
//  //
//
//  auto f_0 = _v0(dispETFESpace, fibers);
//    auto s_0 = _v0(dispETFESpace, sheets);
//
//    boost::shared_ptr<orthonormalizeFibers> normalize0(new orthonormalizeFibers);
//    boost::shared_ptr<orthonormalizeFibers> normalize1(new orthonormalizeFibers(1));
//    auto f0 = eval(normalize0, f_0);
//
//    auto s_00 = s_0 - dot(f0, s_0) * s_0;
//
//    auto s0 = eval(normalize1, s_00);
//
//    if(Case == 0)
//    {
//      std::cout << "EMETA - Computing I4 f jacobian mixed terms second derivative dI1 (Fung) ... \n";
//
//      auto dP = eval (dW4, _F( dispETFESpace, disp, 0 ), f0, s0 )
//             *  _dI1bardF( dispETFESpace, disp, 0 )
//             *  _dI4( dispETFESpace, disp, 0, f0 );
//
//      integrate ( elements ( dispETFESpace->mesh() ) ,
//                  quadRuleTetra4pt,
//                  dispETFESpace,
//                  dispETFESpace,
//                  dot ( dP , grad (phi_i) )
//                  ) >> jacobianPtr;
//    }
//    else
//    {
//      std::cout << "EMETA - Computing I4 s jacobian mixed terms second derivative dI1 (Fung) ... \n";
//
//      auto dP = eval (dW4, _F( dispETFESpace, disp, 0), f0, s0 )
//                     *  _dI1bardF( dispETFESpace, disp, 0 )
//             *  _dI4( dispETFESpace, disp, 0, s0 );
//
//      integrate ( elements ( dispETFESpace->mesh() ) ,
//                  quadRuleTetra4pt,
//                  dispETFESpace,
//                  dispETFESpace,
//                  dot ( dP , grad (phi_i) )
//                  ) >> jacobianPtr;
//    }
//
//}



//template< typename Mesh, typename FunctorPtr >
//void
//computeI4JacobianTermsMixedFungSecondDerivative( const vector_Type& disp,
//      boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > >  dispETFESpace,
//                      const vector_Type& fibers,
//                      const vector_Type& sheets,
//                      matrixPtr_Type     jacobianPtr,
//                      FunctorPtr         dW4,
//                      Int Case,
//                      Int ShearCase)
//{
//  using namespace ExpressionAssembly;
//  //
//
//  auto f_0 = _v0(dispETFESpace, fibers);
//    auto s_0 = _v0(dispETFESpace, sheets);
//
//    boost::shared_ptr<orthonormalizeFibers> normalize0(new orthonormalizeFibers);
//    boost::shared_ptr<orthonormalizeFibers> normalize1(new orthonormalizeFibers(1));
//    auto f0 = eval(normalize0, f_0);
//
//    auto s_00 = s_0 - dot(f0, s_0) * s_0;
//
//    auto s0 = eval(normalize1, s_00);
//  boost::shared_ptr<CrossProduct> wedge(new CrossProduct);
//  auto n0 = eval( wedge, f0, s0);
//
//    if(Case == 0)
//    {
//        if(ShearCase == 0)
//        {
//      std::cout << "EMETA - Computing I4 f jacobian mixed terms second derivative dI8 fs (Fung) ... \n";
//
//      auto dP = eval (dW4, _F( dispETFESpace, disp, 0 ), f0, s0 )
//             *  _dI8dF( dispETFESpace, disp, 0, f0, s0 )
//             *  _dI4( dispETFESpace, disp, 0, f0 );
//
//      integrate ( elements ( dispETFESpace->mesh() ) ,
//                  quadRuleTetra4pt,
//                  dispETFESpace,
//                  dispETFESpace,
//                  dot ( dP , grad (phi_i) )
//                  ) >> jacobianPtr;
//        }
//        else if(ShearCase == 1)
//        {
//      std::cout << "EMETA - Computing I4 f jacobian mixed terms second derivative dI8 fn (Fung) ... \n";
//
//      auto dP = eval (dW4, _F( dispETFESpace, disp, 0 ), f0, s0 )
//             *  _dI8dF( dispETFESpace, disp, 0, f0, n0 )
//             *  _dI4( dispETFESpace, disp, 0, f0 );
//
//      integrate ( elements ( dispETFESpace->mesh() ) ,
//                  quadRuleTetra4pt,
//                  dispETFESpace,
//                  dispETFESpace,
//                  dot ( dP , grad (phi_i) )
//                  ) >> jacobianPtr;
//        }
//        else
//        {
//      std::cout << "EMETA - Computing I4 f jacobian mixed terms second derivative dI8 sn (Fung) ... \n";
//
//      auto dP = eval (dW4, _F( dispETFESpace, disp, 0 ), f0, s0 )
//             *  _dI8dF( dispETFESpace, disp, 0, s0, n0 )
//             *  _dI4( dispETFESpace, disp, 0, f0 );
//
//      integrate ( elements ( dispETFESpace->mesh() ) ,
//                  quadRuleTetra4pt,
//                  dispETFESpace,
//                  dispETFESpace,
//                  dot ( dP , grad (phi_i) )
//                  ) >> jacobianPtr;
//        }
//    }
//    else
//    {
//        if(ShearCase == 0)
//        {
//      std::cout << "EMETA - Computing I4 s jacobian mixed terms second derivative dI8 fs (Fung) ... \n";
//
//      auto dP = eval (dW4, _F( dispETFESpace, disp, 0 ), f0, s0 )
//             *  _dI8dF( dispETFESpace, disp, 0, f0, s0 )
//             *  _dI4( dispETFESpace, disp, 0, s0 );
//
//      integrate ( elements ( dispETFESpace->mesh() ) ,
//                  quadRuleTetra4pt,
//                  dispETFESpace,
//                  dispETFESpace,
//                  dot ( dP , grad (phi_i) )
//                  ) >> jacobianPtr;
//        }
//        else if(ShearCase == 1)
//        {
//      std::cout << "EMETA - Computing I4 s jacobian mixed terms second derivative dI8 fn (Fung) ... \n";
//
//      auto dP = eval (dW4, _F( dispETFESpace, disp, 0 ), f0, s0 )
//             *  _dI8dF( dispETFESpace, disp, 0, f0, n0 )
//             *  _dI4( dispETFESpace, disp, 0, s0 );
//
//      integrate ( elements ( dispETFESpace->mesh() ) ,
//                  quadRuleTetra4pt,
//                  dispETFESpace,
//                  dispETFESpace,
//                  dot ( dP , grad (phi_i) )
//                  ) >> jacobianPtr;
//        }
//        else
//        {
//      std::cout << "EMETA - Computing I4 s jacobian mixed terms second derivative dI8 sn (Fung) ... \n";
//
//      auto dP = eval (dW4, _F( dispETFESpace, disp, 0 ), f0, s0 )
//             *  _dI8dF( dispETFESpace, disp, 0, s0, n0 )
//             *  _dI4( dispETFESpace, disp, 0, s0 );
//
//      integrate ( elements ( dispETFESpace->mesh() ) ,
//                  quadRuleTetra4pt,
//                  dispETFESpace,
//                  dispETFESpace,
//                  dot ( dP , grad (phi_i) )
//                  ) >> jacobianPtr;
//        }
//    }
//
//}



template< typename Mesh, typename FunctorPtr >
void
computeI8JacobianTerms ( const vector_Type& disp,
                         boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > >  dispETFESpace,
                         const vector_Type& fibers,
                         const vector_Type& sheets,
                         matrixPtr_Type     jacobianPtr,
                         FunctorPtr         W8,
                         bool orthonormalize = true)
{
    using namespace ExpressionAssembly;

    if (orthonormalize)
    {
        auto f_0 = _v0 (dispETFESpace, fibers);
        auto s_0 = _v0 (dispETFESpace, sheets);

        boost::shared_ptr<orthonormalizeFibers> normalize0 (new orthonormalizeFibers);
        boost::shared_ptr<orthonormalizeFibers> normalize1 (new orthonormalizeFibers (1) );
        auto f0 = eval (normalize0, f_0);

        auto s_00 = s_0 - dot (f0, s_0) * f0;

        auto s0 = eval (normalize1, s_00);

        auto dP = eval (W8, _I8 ( dispETFESpace, disp, 0, f0, s0 ) )
                  *  _d2I8dF (dispETFESpace, f0, s0 );

    	if(disp.comm().MyPID() == 0)
        std::cout << "EMETA - Computing I8 jacobian terms orthonormalizing ... \n";
        integrate ( elements ( dispETFESpace->mesh() ) ,
                    quadRuleTetra4pt,
                    dispETFESpace,
                    dispETFESpace,
                    dot ( dP, grad (phi_i) )
                  ) >> jacobianPtr;
    }
    else
    {
        auto f_0 = _v0 (dispETFESpace, fibers);
        auto s_0 = _v0 (dispETFESpace, sheets);

        boost::shared_ptr<orthonormalizeFibers> normalize0 (new orthonormalizeFibers);
        boost::shared_ptr<orthonormalizeFibers> normalize1 (new orthonormalizeFibers (1) );

        auto f0 = eval (normalize0, f_0);
        auto s0 = eval (normalize1, s_0);


        auto dP = eval (W8, _I8 ( dispETFESpace, disp, 0, f0, s0 ) )
                  *  _d2I8dF (dispETFESpace, f0, s0 );

    	if(disp.comm().MyPID() == 0)
        std::cout << "EMETA - Computing I8 jacobian terms ... \n";
        integrate ( elements ( dispETFESpace->mesh() ) ,
                    quadRuleTetra4pt,
                    dispETFESpace,
                    dispETFESpace,
                    dot ( dP, grad (phi_i) )
                  ) >> jacobianPtr;
    }

}


template< typename Mesh, typename FunctorPtr >
void
computeI8JacobianTermsSecondDerivative ( const vector_Type& disp,
                                         boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > >  dispETFESpace,
                                         const vector_Type& fibers,
                                         const vector_Type& sheets,
                                         matrixPtr_Type     jacobianPtr,
                                         FunctorPtr         dW8,
                                         bool orthonormalize = true)
{
    using namespace ExpressionAssembly;

    if (orthonormalize)
    {
        auto f_0 = _v0 (dispETFESpace, fibers);
        auto s_0 = _v0 (dispETFESpace, sheets);

        boost::shared_ptr<orthonormalizeFibers> normalize0 (new orthonormalizeFibers);
        boost::shared_ptr<orthonormalizeFibers> normalize1 (new orthonormalizeFibers (1) );
        auto f0 = eval (normalize0, f_0);

        auto s_00 = s_0 - dot (f0, s_0) * f0;

        auto s0 = eval (normalize1, s_00);

        auto dP = eval (dW8, _I8 ( dispETFESpace, disp, 0, f0, s0 ) )
                  *  _dI8dF ( dispETFESpace, disp, 0, f0, s0 )
                  *  _dI8 ( dispETFESpace, disp, 0, f0, s0 );

    	if(disp.comm().MyPID() == 0)
    		std::cout << "EMETA - Computing I8 jacobian terms second derivative orthonormalizing ... \n";
        integrate ( elements ( dispETFESpace->mesh() ) ,
                    quadRuleTetra4pt,
                    dispETFESpace,
                    dispETFESpace,
                    dot (  dP, grad (phi_i) )
                  ) >> jacobianPtr;
    }
    else
    {
        auto f_0 = _v0 (dispETFESpace, fibers);
        auto s_0 = _v0 (dispETFESpace, sheets);

        boost::shared_ptr<orthonormalizeFibers> normalize0 (new orthonormalizeFibers);
        boost::shared_ptr<orthonormalizeFibers> normalize1 (new orthonormalizeFibers (1) );

        auto f0 = eval (normalize0, f_0);
        auto s0 = eval (normalize1, s_0);


        auto dP = eval (dW8, _I8 ( dispETFESpace, disp, 0, f0, s0 ) )
                  *  _dI8dF ( dispETFESpace, disp, 0, f0, s0 )
                  *  _dI8 ( dispETFESpace, disp, 0, f0, s0 );

    	if(disp.comm().MyPID() == 0)
        std::cout << "EMETA - Computing I8 jacobian terms second derivative ... \n";
        integrate ( elements ( dispETFESpace->mesh() ) ,
                    quadRuleTetra4pt,
                    dispETFESpace,
                    dispETFESpace,
                    dot (  dP, grad (phi_i) )
                  ) >> jacobianPtr;
    }


}


//
//template< typename Mesh, typename FunctorPtr >
//void
//computeI8JacobianTermsFung( const vector_Type& disp,
//                      boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > >  dispETFESpace,
//                      const vector_Type& fibers,
//                      const vector_Type& sheets,
//                      matrixPtr_Type     jacobianPtr,
//                      FunctorPtr         W8,
//                      Int Case = 0,
//                      bool orthonormalize = true)
//{
//  using namespace ExpressionAssembly;
//
//  if(orthonormalize)
//  {
//      auto f_0 = _v0(dispETFESpace, fibers);
//      auto s_0 = _v0(dispETFESpace, sheets);
//
//      boost::shared_ptr<orthonormalizeFibers> normalize0(new orthonormalizeFibers);
//      boost::shared_ptr<orthonormalizeFibers> normalize1(new orthonormalizeFibers(1));
//      auto f0 = eval(normalize0, f_0);
//
//      auto s_00 = s_0 - dot(f0, s_0) * s_0;
//
//      auto s0 = eval(normalize1, s_00);
//
//      if(Case == 0)
//      {
//          auto dP = eval (W8, _F( dispETFESpace, disp, 0), f0, s0 )
//                 *  _d2I8dF(dispETFESpace, f0, s0 );
//
//          std::cout << "EMETA - Computing I8 fs jacobian terms Fung orthonormalizing ... \n";
//
//          integrate ( elements ( dispETFESpace->mesh() ) ,
//                      quadRuleTetra4pt,
//                      dispETFESpace,
//                      dispETFESpace,
//                      dot ( dP, grad (phi_i) )
//                      ) >> jacobianPtr;
//      }
//      else
//      {
//          boost::shared_ptr<CrossProduct> wedge(new CrossProduct);
//          auto n0 = eval( wedge, f0, s0);
//
//          if(Case == 1)
//          {
//              auto dP = eval (W8, _F( dispETFESpace, disp, 0), f0, s0 )
//                     *  _d2I8dF(dispETFESpace, f0, n0 );
//
//              std::cout << "EMETA - Computing I8 fn jacobian terms Fung orthonormalizing ... \n";
//
//              integrate ( elements ( dispETFESpace->mesh() ) ,
//                          quadRuleTetra4pt,
//                          dispETFESpace,
//                          dispETFESpace,
//                          dot ( dP, grad (phi_i) )
//                          ) >> jacobianPtr;
//          }
//          else
//          {
//              auto dP = eval (W8, _F( dispETFESpace, disp, 0), f0, s0 )
//                     *  _d2I8dF(dispETFESpace, s0, n0 );
//
//              std::cout << "EMETA - Computing I8 fn jacobian terms Fung orthonormalizing ... \n";
//
//              integrate ( elements ( dispETFESpace->mesh() ) ,
//                          quadRuleTetra4pt,
//                          dispETFESpace,
//                          dispETFESpace,
//                          dot ( dP, grad (phi_i) )
//                          ) >> jacobianPtr;
//          }
//      }
//  }
//  else
//  {
//      auto f_0 = _v0(dispETFESpace, fibers);
//      auto s_0 = _v0(dispETFESpace, sheets);
//
//      boost::shared_ptr<orthonormalizeFibers> normalize0(new orthonormalizeFibers);
//      boost::shared_ptr<orthonormalizeFibers> normalize1(new orthonormalizeFibers(1));
//      boost::shared_ptr<orthonormalizeFibers> normalize2(new orthonormalizeFibers(2));
//
//
//      auto f0 = eval(normalize0, f_0);
//      auto s0 = eval(normalize1, s_0);
//
//      if(Case == 0)
//      {
//          auto dP = eval (W8, _I8( dispETFESpace, disp, 0, f0, s0 ) )
//                 *  _d2I8dF(dispETFESpace, f0, s0 );
//
//
//          std::cout << "EMETA - Computing I8 jacobian terms ... \n";
//          integrate ( elements ( dispETFESpace->mesh() ) ,
//                      quadRuleTetra4pt,
//                      dispETFESpace,
//                      dispETFESpace,
//                      dot ( dP, grad (phi_i) )
//                      ) >> jacobianPtr;
//      }
//      else
//      {
//          boost::shared_ptr<CrossProduct> wedge(new CrossProduct);
//          auto n_0 = eval( wedge, f0, s0);
//          auto n0 = eval(normalize2, n_0);
//
//          if(Case == 1)
//          {
//              auto dP = eval (W8, _F( dispETFESpace, disp, 0), f0, s0 )
//                     *  _d2I8dF(dispETFESpace, f0, n0 );
//
//              std::cout << "EMETA - Computing I8 fn jacobian terms Fung orthonormalizing ... \n";
//
//              integrate ( elements ( dispETFESpace->mesh() ) ,
//                          quadRuleTetra4pt,
//                          dispETFESpace,
//                          dispETFESpace,
//                          dot ( dP, grad (phi_i) )
//                          ) >> jacobianPtr;
//          }
//          else
//          {
//              auto dP = eval (W8, _F( dispETFESpace, disp, 0), f0, s0 )
//                     *  _d2I8dF(dispETFESpace, s0, n0 );
//
//              std::cout << "EMETA - Computing I8 fn jacobian terms Fung orthonormalizing ... \n";
//
//              integrate ( elements ( dispETFESpace->mesh() ) ,
//                          quadRuleTetra4pt,
//                          dispETFESpace,
//                          dispETFESpace,
//                          dot ( dP, grad (phi_i) )
//                          ) >> jacobianPtr;
//          }
//      }
//
//  }
//
//}


//template< typename Mesh, typename FunctorPtr >
//void
//computeI8JacobianTermsFungSecondDerivative( const vector_Type& disp,
//                                      boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > >  dispETFESpace,
//                                      const vector_Type& fibers,
//                                      const vector_Type& sheets,
//                                      matrixPtr_Type     jacobianPtr,
//                                      FunctorPtr         dW8,
//                                      Int Case = 0,
//                                      bool orthonormalize = true)
//{
//  using namespace ExpressionAssembly;
//
//  if(orthonormalize)
//  {
//      auto f_0 = _v0(dispETFESpace, fibers);
//      auto s_0 = _v0(dispETFESpace, sheets);
//
//      boost::shared_ptr<orthonormalizeFibers> normalize0(new orthonormalizeFibers);
//      boost::shared_ptr<orthonormalizeFibers> normalize1(new orthonormalizeFibers(1));
//
//      auto f0 = eval(normalize0, f_0);
//
//      auto s_00 = s_0 - dot(f0, s_0) * s_0;
//
//      auto s0 = eval(normalize1, s_00);
//
//
//      if(Case == 0)
//      {
//          auto dP = eval (dW8, _F( dispETFESpace, disp, 0), f0, s0 )
//                 *  _dI8dF( dispETFESpace, disp, 0, f0, s0 )
//                 *  _dI8( dispETFESpace, disp, 0, f0, s0 );
//
//          std::cout << "EMETA - Computing I8 fs jacobian terms Fung second derivative orthonormalizing ... \n";
//
//          integrate ( elements ( dispETFESpace->mesh() ) ,
//                      quadRuleTetra4pt,
//                      dispETFESpace,
//                      dispETFESpace,
//                      dot (  dP, grad (phi_i) )
//                      ) >> jacobianPtr;
//      }
//      else
//      {
//          boost::shared_ptr<CrossProduct> wedge(new CrossProduct);
//          auto n0 = eval( wedge, f0, s0);
//
//          if(Case == 1)
//          {
//              auto dP = eval (dW8, _F( dispETFESpace, disp, 0), f0, s0 )
//                     *  _dI8dF( dispETFESpace, disp, 0, f0, n0 )
//                     *  _dI8( dispETFESpace, disp, 0, f0, n0 );
//
//              std::cout << "EMETA - Computing I8 fn jacobian terms Fung second derivative orthonormalizing ... \n";
//
//              integrate ( elements ( dispETFESpace->mesh() ) ,
//                          quadRuleTetra4pt,
//                          dispETFESpace,
//                          dispETFESpace,
//                          dot (  dP, grad (phi_i) )
//                          ) >> jacobianPtr;
//          }
//          else
//          {
//              auto dP = eval (dW8, _F( dispETFESpace, disp, 0), f0, s0 )
//                     *  _dI8dF( dispETFESpace, disp, 0, s0, n0 )
//                     *  _dI8( dispETFESpace, disp, 0, s0, n0 );
//
//              std::cout << "EMETA - Computing I8 sn jacobian terms second derivative orthonormalizing ... \n";
//
//              integrate ( elements ( dispETFESpace->mesh() ) ,
//                          quadRuleTetra4pt,
//                          dispETFESpace,
//                          dispETFESpace,
//                          dot (  dP, grad (phi_i) )
//                          ) >> jacobianPtr;
//          }
//      }
//
//  }
//  else
//  {
//      auto f_0 = _v0(dispETFESpace, fibers);
//      auto s_0 = _v0(dispETFESpace, sheets);
//
//      boost::shared_ptr<orthonormalizeFibers> normalize0(new orthonormalizeFibers);
//      boost::shared_ptr<orthonormalizeFibers> normalize1(new orthonormalizeFibers(1));
//      boost::shared_ptr<orthonormalizeFibers> normalize2(new orthonormalizeFibers(2));
//
//      auto f0 = eval(normalize0, f_0);
//      auto s0 = eval(normalize1, s_0);
//
//
//      if(Case == 0)
//      {
//          auto dP = eval (dW8, _F( dispETFESpace, disp, 0), f0, s0 )
//                 *  _dI8dF( dispETFESpace, disp, 0, f0, s0 )
//                 *  _dI8( dispETFESpace, disp, 0, f0, s0 );
//
//          std::cout << "EMETA - Computing I8 jacobian terms Fung second derivative  ... \n";
//
//          integrate ( elements ( dispETFESpace->mesh() ) ,
//                      quadRuleTetra4pt,
//                      dispETFESpace,
//                      dispETFESpace,
//                      dot (  dP, grad (phi_i) )
//                      ) >> jacobianPtr;
//      }
//      else
//      {
//          boost::shared_ptr<CrossProduct> wedge(new CrossProduct);
//          auto n_0 = eval( wedge, f0, s0);
//          auto n0 = eval(normalize2, n_0);
//
//          if(Case == 1)
//          {
//              auto dP = eval (dW8, _F( dispETFESpace, disp, 0), f0, s0 )
//                     *  _dI8dF( dispETFESpace, disp, 0, f0, n0 )
//                     *  _dI8( dispETFESpace, disp, 0, f0, n0 );
//
//              std::cout << "EMETA - Computing I8 fn jacobian terms Fung second derivative ... \n";
//
//              integrate ( elements ( dispETFESpace->mesh() ) ,
//                          quadRuleTetra4pt,
//                          dispETFESpace,
//                          dispETFESpace,
//                          dot (  dP, grad (phi_i) )
//                          ) >> jacobianPtr;
//          }
//          else
//          {
//              auto dP = eval (dW8, _F( dispETFESpace, disp, 0), f0, s0 )
//                     *  _dI8dF( dispETFESpace, disp, 0, s0, n0 )
//                     *  _dI8( dispETFESpace, disp, 0, s0, n0 );
//
//              std::cout << "EMETA - Computing I8 sn jacobian terms second derivative... \n";
//
//              integrate ( elements ( dispETFESpace->mesh() ) ,
//                          quadRuleTetra4pt,
//                          dispETFESpace,
//                          dispETFESpace,
//                          dot (  dP, grad (phi_i) )
//                          ) >> jacobianPtr;
//          }
//
//      }
//  }
//
//
//}


//
//template< typename Mesh, typename FunctorPtr >
//void
//computeI8JacobianMixedTermsFungSecondDerivative( const vector_Type& disp,
//                                      boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > >  dispETFESpace,
//                                      const vector_Type& fibers,
//                                      const vector_Type& sheets,
//                                      matrixPtr_Type     jacobianPtr,
//                                      FunctorPtr         dW8,
//                                      Int Case = 0,
//                                      bool orthonormalize = true)
//{
//  using namespace ExpressionAssembly;
//
//  if(orthonormalize)
//  {
//      auto f_0 = _v0(dispETFESpace, fibers);
//      auto s_0 = _v0(dispETFESpace, sheets);
//
//      boost::shared_ptr<orthonormalizeFibers> normalize0(new orthonormalizeFibers);
//      boost::shared_ptr<orthonormalizeFibers> normalize1(new orthonormalizeFibers(1));
//
//      auto f0 = eval(normalize0, f_0);
//
//      auto s_00 = s_0 - dot(f0, s_0) * s_0;
//
//      auto s0 = eval(normalize1, s_00);
//
//
//      if(Case == 0)
//      {
//          auto dP = eval (dW8, _F( dispETFESpace, disp, 0), f0, s0 )
//                 *  _dI1bardF( dispETFESpace, disp, 0 )
//                 *  _dI8( dispETFESpace, disp, 0, f0, s0 );
//
//          std::cout << "EMETA - Computing I8 fs jacobian terms Fung second derivative dI1 orthonormalizing ... \n";
//
//          integrate ( elements ( dispETFESpace->mesh() ) ,
//                      quadRuleTetra4pt,
//                      dispETFESpace,
//                      dispETFESpace,
//                      dot (  dP, grad (phi_i) )
//                      ) >> jacobianPtr;
//      }
//      else
//      {
//          boost::shared_ptr<CrossProduct> wedge(new CrossProduct);
//          auto n0 = eval( wedge, f0, s0);
//
//          if(Case == 1)
//          {
//              auto dP = eval (dW8, _F( dispETFESpace, disp, 0), f0, s0 )
//                         *  _dI1bardF( dispETFESpace, disp, 0 )
//                     *  _dI8( dispETFESpace, disp, 0, f0, n0 );
//
//              std::cout << "EMETA - Computing I8 fn jacobian terms Fung second derivative dI1 orthonormalizing ... \n";
//
//              integrate ( elements ( dispETFESpace->mesh() ) ,
//                          quadRuleTetra4pt,
//                          dispETFESpace,
//                          dispETFESpace,
//                          dot (  dP, grad (phi_i) )
//                          ) >> jacobianPtr;
//          }
//          else
//          {
//              auto dP = eval (dW8, _F( dispETFESpace, disp, 0), f0, s0 )
//                         *  _dI1bardF( dispETFESpace, disp, 0 )
//                     *  _dI8( dispETFESpace, disp, 0, s0, n0 );
//
//              std::cout << "EMETA - Computing I8 sn jacobian terms second derivative dI1 orthonormalizing ... \n";
//
//              integrate ( elements ( dispETFESpace->mesh() ) ,
//                          quadRuleTetra4pt,
//                          dispETFESpace,
//                          dispETFESpace,
//                          dot (  dP, grad (phi_i) )
//                          ) >> jacobianPtr;
//          }
//      }
//
//  }
//
//}




//template< typename Mesh, typename FunctorPtr >
//void
//computeI8JacobianTermsMixedFungSecondDerivative( const vector_Type& disp,
//      boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > >  dispETFESpace,
//                      const vector_Type& fibers,
//                      const vector_Type& sheets,
//                      matrixPtr_Type     jacobianPtr,
//                      FunctorPtr         dW4,
//                      Int Case,
//                      Int ShearCase)
//{
//  using namespace ExpressionAssembly;
//  //
//
//  auto f_0 = _v0(dispETFESpace, fibers);
//    auto s_0 = _v0(dispETFESpace, sheets);
//
//    boost::shared_ptr<orthonormalizeFibers> normalize0(new orthonormalizeFibers);
//    boost::shared_ptr<orthonormalizeFibers> normalize1(new orthonormalizeFibers(1));
//    auto f0 = eval(normalize0, f_0);
//
//    auto s_00 = s_0 - dot(f0, s_0) * s_0;
//
//    auto s0 = eval(normalize1, s_00);
//  boost::shared_ptr<CrossProduct> wedge(new CrossProduct);
//  auto n0 = eval( wedge, f0, s0);
//
//    if(Case == 0)
//    {
//        if(ShearCase == 0)
//        {
//      std::cout << "EMETA - Computing I8 fs jacobian mixed terms second derivative dI4 f (Fung) ... \n";
//
//      auto dP = eval (dW4, _F( dispETFESpace, disp, 0 ), f0, s0 )
//             *  _dI4dF( dispETFESpace, disp, 0, f0 )
//             *  _dI8( dispETFESpace, disp, 0, f0, s0 );
//
//      integrate ( elements ( dispETFESpace->mesh() ) ,
//                  quadRuleTetra4pt,
//                  dispETFESpace,
//                  dispETFESpace,
//                  dot ( dP , grad (phi_i) )
//                  ) >> jacobianPtr;
//        }
//        else if(ShearCase == 1)
//        {
//          std::cout << "EMETA - Computing I8 fn jacobian mixed terms second derivative dI4 f (Fung) ... \n";
//
//          auto dP = eval (dW4, _F( dispETFESpace, disp, 0 ), f0, s0 )
//                 *  _dI4dF( dispETFESpace, disp, 0, f0 )
//                 *  _dI8( dispETFESpace, disp, 0, f0, n0 );
//
//          integrate ( elements ( dispETFESpace->mesh() ) ,
//                      quadRuleTetra4pt,
//                      dispETFESpace,
//                      dispETFESpace,
//                      dot ( dP , grad (phi_i) )
//                      ) >> jacobianPtr;
//        }
//        else
//        {
//          std::cout << "EMETA - Computing I8 sn jacobian mixed terms second derivative dI4 f (Fung) ... \n";
//
//          auto dP = eval (dW4, _F( dispETFESpace, disp, 0 ), f0, s0 )
//                 *  _dI4dF( dispETFESpace, disp, 0, f0 )
//                 *  _dI8( dispETFESpace, disp, 0, s0, n0 );
//
//          integrate ( elements ( dispETFESpace->mesh() ) ,
//                      quadRuleTetra4pt,
//                      dispETFESpace,
//                      dispETFESpace,
//                      dot ( dP , grad (phi_i) )
//                      ) >> jacobianPtr;
//        }
//    }
//    else
//    {
//        if(ShearCase == 0)
//        {
//      std::cout << "EMETA - Computing I8 fs jacobian mixed terms second derivative dI4 s (Fung) ... \n";
//
//      auto dP = eval (dW4, _F( dispETFESpace, disp, 0 ), f0, s0 )
//             *  _dI4dF( dispETFESpace, disp, 0, s0 )
//             *  _dI8( dispETFESpace, disp, 0, f0, s0 );
//
//      integrate ( elements ( dispETFESpace->mesh() ) ,
//                  quadRuleTetra4pt,
//                  dispETFESpace,
//                  dispETFESpace,
//                  dot ( dP , grad (phi_i) )
//                  ) >> jacobianPtr;
//        }
//        else if(ShearCase == 1)
//        {
//          std::cout << "EMETA - Computing I8 fn jacobian mixed terms second derivative dI4 s (Fung) ... \n";
//
//          auto dP = eval (dW4, _F( dispETFESpace, disp, 0 ), f0, s0 )
//                 *  _dI4dF( dispETFESpace, disp, 0, s0 )
//                 *  _dI8( dispETFESpace, disp, 0, f0, n0 );
//
//          integrate ( elements ( dispETFESpace->mesh() ) ,
//                      quadRuleTetra4pt,
//                      dispETFESpace,
//                      dispETFESpace,
//                      dot ( dP , grad (phi_i) )
//                      ) >> jacobianPtr;
//        }
//        else
//        {
//          std::cout << "EMETA - Computing I8 sn jacobian mixed terms second derivative dI4 s (Fung) ... \n";
//
//          auto dP = eval (dW4, _F( dispETFESpace, disp, 0 ), f0, s0 )
//                 *  _dI4dF( dispETFESpace, disp, 0, s0 )
//                 *  _dI8( dispETFESpace, disp, 0, s0, n0 );
//
//          integrate ( elements ( dispETFESpace->mesh() ) ,
//                      quadRuleTetra4pt,
//                      dispETFESpace,
//                      dispETFESpace,
//                      dot ( dP , grad (phi_i) )
//                      ) >> jacobianPtr;
//        }
//    }
//}



}//EMAssembler

}//LifeV

#endif /* EMETAASSMEBLER_HPP_ */
