/*
 * EMETAAssmebler.hpp
 *
 *  Created on: 29/apr/2014
 *      Author: srossi
 */

#ifndef EMETAJACOBIANISOTROPICASSMEBLER_HPP_
#define EMETAJACOBIANISOTROPICASSMEBLER_HPP_


#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/array/MatrixEpetra.hpp>

//ET include for assemblings
#include <lifev/eta/fem/ETFESpace.hpp>
#include <lifev/eta/expression/Integrate.hpp>

#include <lifev/em/solver/mechanics/EMMechanicalExpressions.hpp>

#include <lifev/em/solver/EMETAFunctors.hpp>


namespace LifeV
{

typedef VectorEpetra           vector_Type;
typedef boost::shared_ptr<vector_Type>         vectorPtr_Type;

typedef MatrixEpetra<Real>           matrix_Type;
typedef boost::shared_ptr<matrix_Type>         matrixPtr_Type;

namespace EMAssembler
{

    
template <typename Mesh, typename FunctorPtr1, typename FunctorPtr2, typename FunctorPtr3, typename FunctorPtr4 >
void
computePMRCJacobianTerms ( const vector_Type& disp,
                           boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > >  dispETFESpace,
                           matrixPtr_Type jacobianPtr,
                           FunctorPtr1 Wvol, FunctorPtr2 dWvol, FunctorPtr3 W1, FunctorPtr4 W2)
{
    {
        using namespace ExpressionAssembly;
        
        auto F = _F (dispETFESpace, disp, 0);
        
        //
        if(disp.comm().MyPID() == 0) std::cout << "EMETA - Computing Volumetric jacobian terms  ... \n";
        
        auto dPvol = eval (Wvol, F ) * _d2JdF (F, _dF);
//        integrate ( elements ( dispETFESpace->mesh() ) ,
//                   quadRuleTetra4pt,
//                   dispETFESpace,
//                   dispETFESpace,
//                   dot ( dPvol, grad (phi_i) )
//                   ) >> jacobianPtr;
        
        //
        if(disp.comm().MyPID() == 0) std::cout << "EMETA - Computing Volumetric jacobian terms with second derivative of the energy ... \n";
        
        auto dPdvol = eval (dWvol, F ) * _dJdF (F, _dF) * _dJ (F);
        integrate ( elements ( dispETFESpace->mesh() ) ,
                   quadRuleTetra4pt,
                   dispETFESpace,
                   dispETFESpace,
                   dot ( dPdvol + dPvol , grad (phi_i) )
                   ) >> jacobianPtr;
        
        //
        if(disp.comm().MyPID() == 0) std::cout << "Computing I1 jacobian terms  ... \n";
        
        auto dPdF = eval (W1, F ) * _d2I1bardF(F, _dF);
//        integrate ( elements ( dispETFESpace->mesh() ) ,
//                   quadRuleTetra4pt,
//                   dispETFESpace,
//                   dispETFESpace,
//                   dot ( dPdF , grad (phi_i) )
//                   ) >> jacobianPtr;
        
        //
        if(disp.comm().MyPID() == 0) std::cout << "Computing I2 jacobian terms  ... \n";
        
        auto dP2 = eval (W2, F ) * _d2I2bardF (F, _dF);
        integrate ( elements ( dispETFESpace->mesh() ) ,
                   quadRuleTetra4pt,
                   dispETFESpace,
                   dispETFESpace,
                   dot ( dP2 + dPdF , grad (phi_i) )
                   ) >> jacobianPtr;
        
    }
}


template <typename Mesh, typename FunctorPtr >
void
computeVolumetricJacobianTerms ( const vector_Type& disp,
                                 boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > >  dispETFESpace,
                                 matrixPtr_Type           jacobianPtr,
                                 FunctorPtr                Wvol)
{
    {
        using namespace ExpressionAssembly;
        //
    	if(disp.comm().MyPID() == 0)
    		std::cout << "EMETA - Computing Volumetric jacobian terms  ... \n";

    	auto F = _F (dispETFESpace, disp, 0);


    	auto dP = eval (Wvol, F ) * _d2JdF (F, _dF);
        integrate ( elements ( dispETFESpace->mesh() ) ,
                    quadRuleTetra4pt,
                    dispETFESpace,
                    dispETFESpace,
                    dot ( dP, grad (phi_i) )
                  ) >> jacobianPtr;

    }
}

template <typename Mesh, typename FunctorPtr >
void
computeVolumetricJacobianTermsSecondDerivative ( const vector_Type& disp,
                                                 boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > > dispETFESpace,
                                                 matrixPtr_Type           jacobianPtr,
                                                 FunctorPtr                 dWvol)
{
    {
        using namespace ExpressionAssembly;

        //
    	if(disp.comm().MyPID() == 0)
    		std::cout << "EMETA - Computing Volumetric jacobian terms with second derivative of the energy ... \n";

    	auto F = _F (dispETFESpace, disp, 0);

        auto dP = eval (dWvol, F ) * _dJdF (F, _dF) * _dJ (F);
        integrate ( elements ( dispETFESpace->mesh() ) ,
                    quadRuleTetra4pt,
                    dispETFESpace,
                    dispETFESpace,
                    dot ( dP , grad (phi_i) )
                  ) >> jacobianPtr;

    }
}

template< typename Mesh, typename FunctorPtr >
void
computeLinearizedVolumetricJacobianTerms ( const vector_Type& disp,
                                           boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > > dispETFESpace,
                                           matrixPtr_Type           jacobianPtr,
                                           FunctorPtr               W)
{
    using namespace ExpressionAssembly;
    //
    auto I = _I;
    auto dF = _dF;
    auto dFT = transpose (dF);
    auto dP = eval (W, I ) * (dF + dFT) * value (1. / 2.);

	if(disp.comm().MyPID() == 0)
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
computeLinearizedDeviatoricJacobianTerms ( const vector_Type& disp,
                                           boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > >  dispETFESpace,
                                           matrixPtr_Type           jacobianPtr,
                                           FunctorPtr               W)
{
    using namespace ExpressionAssembly;
    //
    auto I = _I;
    auto dF = _dF;
    auto dFT = transpose (dF);
    auto dP = eval (W, I) * trace (dF + dFT) * value (1. / 2.) * I;

	if(disp.comm().MyPID() == 0)
    std::cout << "Computing linear deviatoric Jacobian terms ... \n";
    integrate ( elements ( dispETFESpace->mesh() ) ,
                quadRuleTetra4pt,
                dispETFESpace,
                dispETFESpace,
                dot ( dP  , grad (phi_i) )
              ) >> jacobianPtr;
}

template <typename Mesh, typename FunctorPtr >
void
computeI1JacobianTerms ( const vector_Type& disp,
                         boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > >  dispETFESpace,
                         matrixPtr_Type           jacobianPtr,
                         FunctorPtr                 W1)
{
    using namespace ExpressionAssembly;

    //
	if(disp.comm().MyPID() == 0)
		std::cout << "Computing I1 jacobian terms  ... \n";
//	auto I = _I;
//	auto dF = grad(phi_j);
//	auto GradU = _Grad_u(dispETFESpace, disp, 0);
//	auto F = I + GradU;
//	auto FmT = minusT(F);
//	auto dFmTdF = value(-1.0) * FmT * transpose(dF) * FmT;
//	auto J = det(F);
//	auto Jm23 = pow(J, - 2. / 3.);
//	auto dJm23 = value(- 2. / 3. ) * Jm23 * FmT;
//	auto d2Jm23dF = value( -2. / 3. ) * ( dot( dJm23, dF ) * FmT + Jm23 * dFmTdF );
//	auto I1 = dot(F, F);
//	auto dI1 = value(2.0) * F;
//	auto d2I1 = value(2.0) * dF;
//	auto I1bar = Jm23 * I1;
//	auto dI1bar = dJm23 * I1 + Jm23 * dI1;
//	auto d2I1bardF = dot(dJm23, dF) * dI1 + Jm23 * d2I1 + dJm23 * dot(dI1, dF) + d2Jm23dF * I1;
//	auto P = eval (W1, _F (dispETFESpace, disp, 0) ) * dI1bar ;
//
//	auto dPdF = eval (W1, _F (dispETFESpace, disp, 0) ) * d2I1bardF;

	auto F = _F (dispETFESpace, disp, 0);
	auto dPdF = eval (W1, F ) * _d2I1bardF(F, _dF);

    integrate ( elements ( dispETFESpace->mesh() ) ,
                quadRuleTetra4pt,
                dispETFESpace,
                dispETFESpace,
                dot ( dPdF , grad (phi_i) )
              ) >> jacobianPtr;
}


template <typename Mesh, typename FunctorPtr>
void
computeI1JacobianTermsSecondDerivative ( const vector_Type& disp,
                                         boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > >  dispETFESpace,
                                         matrixPtr_Type           jacobianPtr,
                                         FunctorPtr                 dW1)
{
    using namespace ExpressionAssembly;
    //
	if(disp.comm().MyPID() == 0)
		std::cout << "Computing I1 jacobian terms with second derivative of the energy ... \n";

//	auto I = _I;
//	auto dF = grad(phi_j);
//	auto GradU = _Grad_u(dispETFESpace, disp, 0);
//	auto F = I + GradU;
//	auto FmT = minusT(F);
//	auto dFmTdF = value(-1.0) * FmT * transpose(dF) * FmT;
//	auto J = det(F);
//	auto Jm23 = pow(J, 2 / (-3.) );
//	auto dJm23 = value(2/(-3.)) * pow(J, 2 / (-3.) ) * FmT;
//	auto d2Jm23dF = value( 2/(-3.) ) * ( dot( dJm23, dF ) * FmT + Jm23 * dFmTdF );
//	auto I1 = dot(F, F);
//	auto dI1 = value(2.0) * F;
//	auto I1bar = Jm23 * I1;
//	auto dI1bar = dJm23 * I1 + Jm23 * dI1;
//	auto dI1bardF = dot( dI1bar, dF);
////	auto d2I1bardF = dot(dJm23, dF) * dI1 + Jm23 * dF + dJm23 * dot(dI1, dF) + d2Jm23dF * I1;
//	auto P = eval (dW1, _F (dispETFESpace, disp, 0) ) * dI1bar ;
//
//    auto dPdF = eval (dW1, _F (dispETFESpace, disp, 0) ) * ( dI1bardF ) * ( dI1bar );

	auto F = _F (dispETFESpace, disp, 0);
    auto dPdF = eval (dW1, F ) * _dI1bardF (F, _dF) * _dI1bar (F) ;

    integrate ( elements ( dispETFESpace->mesh() ) ,
                quadRuleTetra4pt,
                dispETFESpace,
                dispETFESpace,
                dot ( dPdF , grad (phi_i) )
              ) >> jacobianPtr;
}


template <typename Mesh, typename FunctorPtr>
void
computeI1JacobianMixedTermsSecondDerivative ( const vector_Type& disp,
                                              boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > >  dispETFESpace,
                                              matrixPtr_Type           jacobianPtr,
                                              FunctorPtr                 dW1dI2)
{
    using namespace ExpressionAssembly;
    //
	if(disp.comm().MyPID() == 0)
    std::cout << "Computing I1 jacobian mixed terms with second derivative of the energy  (derivative on I2)... \n";

	auto F = _F (dispETFESpace, disp, 0);

    auto dP = eval (dW1dI2, F ) * _dI2bardF (F, _dF)  * dI1bar (F);
    integrate ( elements ( dispETFESpace->mesh() ) ,
                quadRuleTetra4pt,
                dispETFESpace,
                dispETFESpace,
                dot ( dP , grad (phi_i) )
              ) >> jacobianPtr;
}


template <typename Mesh, typename FunctorPtr >
void
computeI2JacobianTerms ( const vector_Type& disp,
                         boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > >  dispETFESpace,
                         matrixPtr_Type           jacobianPtr,
                         FunctorPtr                 W2)
{
    using namespace ExpressionAssembly;

    //
	if(disp.comm().MyPID() == 0)
    std::cout << "Computing I2 jacobian terms  ... \n";
	auto F = _F (dispETFESpace, disp, 0);

    auto dP = eval (W2, F ) * _d2I2bardF (F, _dF);
    integrate ( elements ( dispETFESpace->mesh() ) ,
                quadRuleTetra4pt,
                dispETFESpace,
                dispETFESpace,
                dot ( dP , grad (phi_i) )
              ) >> jacobianPtr;
}

template <typename Mesh, typename FunctorPtr>
void
computeI2JacobianTermsSecondDerivative ( const vector_Type& disp,
                                         boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > >  dispETFESpace,
                                         matrixPtr_Type           jacobianPtr,
                                         FunctorPtr                 dW2)
{
    using namespace ExpressionAssembly;
    //
	if(disp.comm().MyPID() == 0)
    std::cout << "Computing I2 jacobian terms with second derivative of the energy ... \n";
	auto F = _F (dispETFESpace, disp, 0);

    auto dP = eval (dW2, F ) * _dI2bardF (F, _dF)  * _dI2bar (F);
    integrate ( elements ( dispETFESpace->mesh() ) ,
                quadRuleTetra4pt,
                dispETFESpace,
                dispETFESpace,
                dot ( dP , grad (phi_i) )
              ) >> jacobianPtr;
}



template <typename Mesh, typename FunctorPtr>
void
computeI2JacobianMixedTermsSecondDerivative ( const vector_Type& disp,
                                              boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > >  dispETFESpace,
                                              matrixPtr_Type           jacobianPtr,
                                              FunctorPtr                 dW2dI1)
{
    using namespace ExpressionAssembly;
    //
	if(disp.comm().MyPID() == 0)
    std::cout << "Computing I2 jacobian mixed terms with second derivative of the energy (derivative on I1 ) ... \n";
	auto F = _F (dispETFESpace, disp, 0);

    auto dP = eval (dW2dI1, F ) * (_dI1bardF (F, _dF ) ) * (_dI2bar (F) );
    integrate ( elements ( dispETFESpace->mesh() ) ,
                quadRuleTetra4pt,
                dispETFESpace,
                dispETFESpace,
                dot ( dP , grad (phi_i) )
              ) >> jacobianPtr;
}





template <typename Mesh, typename FunctorPtr>
void
computeI1JacobianTermsSecondDerivative ( const vector_Type& disp,
                                         boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > >  dispETFESpace,
                                         const vector_Type& fibers,
                                         const vector_Type& sheets,
                                         matrixPtr_Type           jacobianPtr,
                                         FunctorPtr                 dW1)
{
    using namespace ExpressionAssembly;

    auto f_0 = _v0 (dispETFESpace, fibers);
    auto s_0 = _v0 (dispETFESpace, sheets);

    boost::shared_ptr<orthonormalizeFibers> normalize0 (new orthonormalizeFibers);
    boost::shared_ptr<orthonormalizeFibers> normalize1 (new orthonormalizeFibers (1) );
    auto f0 = eval (normalize0, f_0);

    auto s_00 = s_0 - dot (f0, s_0) * s_0;

    auto s0 = eval (normalize1, s_00);
    //
	if(disp.comm().MyPID() == 0)
    std::cout << "Computing I1 jacobian terms with second derivative of the energy  (Fung)... \n";

	auto F = _F (dispETFESpace, disp, 0);

    auto dP = eval (dW1, F, f0, s0) * ( _dI1bardF (F, _dF) ) * ( _dI1bar (F) );
    integrate ( elements ( dispETFESpace->mesh() ) ,
                quadRuleTetra4pt,
                dispETFESpace,
                dispETFESpace,
                dot ( dP , grad (phi_i) )
              ) >> jacobianPtr;
}

template <typename Mesh, typename FunctorPtr>
void
computeI1JacobianTermsSecondDerivative ( const vector_Type& disp,
                                         boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > >  dispETFESpace,
                                         const vector_Type& fibers,
                                         const vector_Type& sheets,
                                         matrixPtr_Type           jacobianPtr,
                                         FunctorPtr                 dW1,
                                         Int Case)
{
    using namespace ExpressionAssembly;

    auto f_0 = _v0 (dispETFESpace, fibers);
    boost::shared_ptr<orthonormalizeFibers> normalize0 (new orthonormalizeFibers);
    auto f0 = eval (normalize0, f_0);

    auto s_0 = _v0 (dispETFESpace, sheets);
    boost::shared_ptr<orthonormalizeFibers> normalize1 (new orthonormalizeFibers (1) );
    auto s_00 = s_0 - dot (f0, s_0) * s_0;
    auto s0 = eval (normalize1, s_00);

    boost::shared_ptr<CrossProduct> wedge (new CrossProduct);
    auto n0 = eval ( wedge, f0, s0);

	auto F = _F (dispETFESpace, disp, 0);

    if (Case == 0)
    {
    	if(disp.comm().MyPID() == 0)
        std::cout << "Computing I1 jacobian terms with second derivative of the energy dI4 f (Fung)... \n";

        auto dP = eval (dW1, F, f0, s0) * ( _dI4dF ( F, f0, _dF ) ) * (_dI1bar (F) );
        integrate ( elements ( dispETFESpace->mesh() ) ,
                    quadRuleTetra4pt,
                    dispETFESpace,
                    dispETFESpace,
                    dot ( dP , grad (phi_i) )
                  ) >> jacobianPtr;
    }
    else if (Case == 1)
    {
    	if(disp.comm().MyPID() == 0)
        std::cout << "Computing I1 jacobian terms with second derivative of the energy dI4 s (Fung)... \n";

        auto dP = eval (dW1, F, f0, s0) * ( _dI4dF ( F, s0, _dF ) ) * (_dI1bar (F) );
        integrate ( elements ( dispETFESpace->mesh() ) ,
                    quadRuleTetra4pt,
                    dispETFESpace,
                    dispETFESpace,
                    dot ( dP , grad (phi_i) )
                  ) >> jacobianPtr;
    }
    else if (Case == 2)
    {
    	if(disp.comm().MyPID() == 0)
        std::cout << "Computing I1 jacobian terms with second derivative of the energy dI8 fs (Fung)... \n";

        auto dP = eval (dW1, F, f0, s0)
                  * ( _dI8dF ( F, f0, s0, _dF ) )
                  * (_dI1bar (F) );
        integrate ( elements ( dispETFESpace->mesh() ) ,
                    quadRuleTetra4pt,
                    dispETFESpace,
                    dispETFESpace,
                    dot ( dP , grad (phi_i) )
                  ) >> jacobianPtr;
    }
    else if (Case == 3)
    {
    	if(disp.comm().MyPID() == 0)
        std::cout << "Computing I1 jacobian terms with second derivative of the energy dI8 fn (Fung)... \n";

        auto dP = eval (dW1, F, f0, s0)
                  * ( _dI8dF ( F, f0, n0, _dF ) )
                  * (_dI1bar (F) );
        integrate ( elements ( dispETFESpace->mesh() ) ,
                    quadRuleTetra4pt,
                    dispETFESpace,
                    dispETFESpace,
                    dot ( dP , grad (phi_i) )
                  ) >> jacobianPtr;
    }
    else if (Case == 4)
    {
    	if(disp.comm().MyPID() == 0)
        std::cout << "Computing I1 jacobian terms with second derivative of the energy dI8 sn (Fung)... \n";

        auto dP = eval (dW1, F, f0, s0)
                  * ( _dI8dF ( F, s0, n0, _dF ) )
                  * (_dI1bar (F) );
        integrate ( elements ( dispETFESpace->mesh() ) ,
                    quadRuleTetra4pt,
                    dispETFESpace,
                    dispETFESpace,
                    dot ( dP , grad (phi_i) )
                  ) >> jacobianPtr;
    }
}





}//EMAssembler

}//LifeV

#endif /* EMETAASSMEBLER_HPP_ */
