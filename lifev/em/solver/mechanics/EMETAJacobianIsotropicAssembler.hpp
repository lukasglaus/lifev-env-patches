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
        std::cout << "Computing Volumetric jacobian terms  ... \n";

        //      BOOST_AUTO_TPL (dP, eval(Wvol, _F) * (_d2JdF) );
        integrate ( elements ( dispETFESpace->mesh() ) ,
                    quadRuleTetra4pt,
                    dispETFESpace,
                    dispETFESpace,
                    dot ( eval (Wvol, _F (dispETFESpace, disp, 0) ) * (_d2JdF (dispETFESpace, disp, 0) ), grad (phi_i) )
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
        std::cout << "Computing Volumetric jacobian terms with second derivative of the energy ... \n";

        auto dP = eval (dWvol, _F (dispETFESpace, disp, 0) ) * (_dJdF (dispETFESpace, disp, 0) ) * (_dJ (dispETFESpace, disp, 0) );
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
    auto I = value (EMUtility::identity() );
    auto dF = _dF;
    auto dFT = transpose (dF);
    auto dP = eval (W, I ) * (dF + dFT) * value (1. / 2.);

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
    auto I = value (EMUtility::identity() );
    auto dF = _dF;
    auto dFT = transpose (dF);
    auto dP = eval (W, I) * trace (dF + dFT) * value (1. / 2.) * I;

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
computeI1JacobianTermsSecondDerivative ( const vector_Type& disp,
                                         boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > >  dispETFESpace,
                                         matrixPtr_Type           jacobianPtr,
                                         FunctorPtr                 dW1)
{
    using namespace ExpressionAssembly;
    //
    std::cout << "Computing I1 jacobian terms with second derivative of the energy ... \n";

    auto dP = eval (dW1, _F (dispETFESpace, disp, 0) ) * (_dI1bardF (dispETFESpace, disp, 0) ) * (_dI1bar (dispETFESpace, disp, 0) );
    integrate ( elements ( dispETFESpace->mesh() ) ,
                quadRuleTetra4pt,
                dispETFESpace,
                dispETFESpace,
                dot ( dP , grad (phi_i) )
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
    std::cout << "Computing I1 jacobian terms  ... \n";

    auto dP = eval (W1, _F (dispETFESpace, disp, 0) ) * (_d2I1bardF (dispETFESpace, disp, 0) );
    integrate ( elements ( dispETFESpace->mesh() ) ,
                quadRuleTetra4pt,
                dispETFESpace,
                dispETFESpace,
                dot ( dP , grad (phi_i) )
              ) >> jacobianPtr;
}

template <typename Mesh, typename FunctorPtr >
void
computeI1JacobianTerms ( const vector_Type& disp,
                         boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > >  dispETFESpace,
                         const vector_Type& fibers,
                         const vector_Type& sheets,
                         matrixPtr_Type           jacobianPtr,
                         FunctorPtr                 W1)
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
    std::cout << "Computing I1 jacobian terms (Fung) ... \n";

    auto dP = eval (W1, _F (dispETFESpace, disp, 0), f0, s0) * (_d2I1bardF (dispETFESpace, disp, 0) );
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
    std::cout << "Computing I1 jacobian terms with second derivative of the energy  (Fung)... \n";

    auto dP = eval (dW1, _F (dispETFESpace, disp, 0), f0, s0) * (_dI1bardF (dispETFESpace, disp, 0) ) * (_dI1bar (dispETFESpace, disp, 0) );
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


    if (Case == 0)
    {
        std::cout << "Computing I1 jacobian terms with second derivative of the energy dI4 f (Fung)... \n";

        auto dP = eval (dW1, _F (dispETFESpace, disp, 0), f0, s0) * ( _dI4dF ( dispETFESpace, disp, 0, f0 ) ) * (_dI1bar (dispETFESpace, disp, 0) );
        integrate ( elements ( dispETFESpace->mesh() ) ,
                    quadRuleTetra4pt,
                    dispETFESpace,
                    dispETFESpace,
                    dot ( dP , grad (phi_i) )
                  ) >> jacobianPtr;
    }
    else if (Case == 1)
    {
        std::cout << "Computing I1 jacobian terms with second derivative of the energy dI4 s (Fung)... \n";

        auto dP = eval (dW1, _F (dispETFESpace, disp, 0), f0, s0) * ( _dI4dF ( dispETFESpace, disp, 0, s0 ) ) * (_dI1bar (dispETFESpace, disp, 0) );
        integrate ( elements ( dispETFESpace->mesh() ) ,
                    quadRuleTetra4pt,
                    dispETFESpace,
                    dispETFESpace,
                    dot ( dP , grad (phi_i) )
                  ) >> jacobianPtr;
    }
    else if (Case == 2)
    {
        std::cout << "Computing I1 jacobian terms with second derivative of the energy dI8 fs (Fung)... \n";

        auto dP = eval (dW1, _F (dispETFESpace, disp, 0), f0, s0)
                  * ( _dI8dF ( dispETFESpace, disp, 0, f0, s0 ) )
                  * (_dI1bar (dispETFESpace, disp, 0) );
        integrate ( elements ( dispETFESpace->mesh() ) ,
                    quadRuleTetra4pt,
                    dispETFESpace,
                    dispETFESpace,
                    dot ( dP , grad (phi_i) )
                  ) >> jacobianPtr;
    }
    else if (Case == 3)
    {
        std::cout << "Computing I1 jacobian terms with second derivative of the energy dI8 fn (Fung)... \n";

        auto dP = eval (dW1, _F (dispETFESpace, disp, 0), f0, s0)
                  * ( _dI8dF ( dispETFESpace, disp, 0, f0, n0 ) )
                  * (_dI1bar (dispETFESpace, disp, 0) );
        integrate ( elements ( dispETFESpace->mesh() ) ,
                    quadRuleTetra4pt,
                    dispETFESpace,
                    dispETFESpace,
                    dot ( dP , grad (phi_i) )
                  ) >> jacobianPtr;
    }
    else if (Case == 4)
    {
        std::cout << "Computing I1 jacobian terms with second derivative of the energy dI8 sn (Fung)... \n";

        auto dP = eval (dW1, _F (dispETFESpace, disp, 0), f0, s0)
                  * ( _dI8dF ( dispETFESpace, disp, 0, s0, n0 ) )
                  * (_dI1bar (dispETFESpace, disp, 0) );
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
