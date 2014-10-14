/*
 * EMMaterialFunctions.hpp
 *
 *  Created on: 28/apr/2014
 *      Author: srossi
 */

#ifndef FUNCTIONSISOTROPICFUNG_HPP_
#define FUNCTIONSISOTROPICFUNG_HPP_

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

namespace IsotropicFungFunctions
{

template<class Mesh>
class W1
{
public:
    typedef typename MaterialFunctions::EMMaterialFunctions<Mesh>::return_Type return_Type;

    W1(Real C = 3330., Real b = 9.242) : M_C (C), M_b (b)  {}
    W1 (const W1& w1)
    {
        M_C = w1.M_C;
        M_b = w1.M_b;
    }
    virtual return_Type operator() (const MatrixSmall<3, 3>& F)
    {
        auto I1bar = Elasticity::I1bar (F);
        auto I2bar = Elasticity::I2bar (F);
        auto Q = 0.25 * M_b * ( I1bar * I1bar - 2.0 * I1bar - 2.0 * I2bar + 3.0  );
//        auto Q = 0.25 * M_b * ( I1bar * I1bar - 2.0 * I1bar + 2.0 * I2bar - 9.0  );
//        auto Q = M_b * ( 0.75 - 0.5 * I1bar + 0.25 * I1bar * I1bar - 0.5 * I2bar );
        auto dQdI1bar = 0.5 * M_b * (  I1bar - 1.0);
        return 0.5 * M_C * std::exp(Q) * dQdI1bar;
    }


    Real M_C;
    Real M_b;
};


template<class Mesh>
class dW1dI1
{
public:
    typedef typename MaterialFunctions::EMMaterialFunctions<Mesh>::return_Type return_Type;

    dW1dI1(Real C = 3330., Real b = 9.242) : M_C (C), M_b (b)  {}
    dW1dI1 (const dW1dI1& w1)
    {
        M_C = w1.M_C;
        M_b = w1.M_b;
    }
    virtual return_Type operator() (const MatrixSmall<3, 3>& F)
    {
        auto I1bar = Elasticity::I1bar (F);
        auto I2bar = Elasticity::I2bar (F);
//        auto Q = M_b * ( 0.75 - 0.5 * I1bar + 0.25 * I1bar * I1bar - 0.5 * I2bar );
        auto Q = 0.25 * M_b * ( I1bar * I1bar - 2.0 * I1bar - 2.0 * I2bar + 3.0  );
//        auto Q = 0.25 * M_b * ( I1bar * I1bar - 2.0 * I1bar + 2.0 * I2bar - 9.0  );
        auto dQdI1bar = 0.5 * M_b * ( I1bar - 1.0 );
        auto d2QdI1bar = 0.5 * M_b;

        return 0.5 * M_C * std::exp(Q) * d2QdI1bar + 0.5 * M_C * std::exp(Q) * dQdI1bar * dQdI1bar;
    }


    Real M_C;
    Real M_b;
};


template<class Mesh>
class dW1dI2
{
public:
    typedef typename MaterialFunctions::EMMaterialFunctions<Mesh>::return_Type return_Type;

    dW1dI2(Real C = 3330., Real b = 9.242) : M_C (C), M_b (b)  {}
    dW1dI2 (const dW1dI2& w1)
    {
        M_C = w1.M_C;
        M_b = w1.M_b;
    }
    virtual return_Type operator() (const MatrixSmall<3, 3>& F)
    {
        auto I1bar = Elasticity::I1bar (F);
        auto I2bar = Elasticity::I2bar (F);
//        auto Q = M_b * ( 0.75 - 0.5 * I1bar + 0.25 * I1bar * I1bar - 0.5 * I2bar );
        auto Q = 0.25 * M_b * ( I1bar * I1bar - 2.0 * I1bar - 2.0 * I2bar + 3.0  );
//        auto Q = 0.25 * M_b * ( I1bar * I1bar - 2.0 * I1bar + 2.0 * I2bar - 9.0  );
        auto dQdI1bar = 0.5 * M_b * ( I1bar - 1.0 );
        auto dQdI2bar =   - 0.5 * M_b;

//        return 0.0;
        return 0.5 * M_C * std::exp(Q) * dQdI2bar * dQdI1bar;
    }


    Real M_C;
    Real M_b;
};

template<class Mesh>
class W2
{
public:
    typedef typename MaterialFunctions::EMMaterialFunctions<Mesh>::return_Type return_Type;

    W2(Real C = 3330., Real b = 9.242) : M_C (C), M_b (b)  {}
    W2 (const W2& w2)
    {
        M_C = w2.M_C;
        M_b = w2.M_b;
    }
    virtual return_Type operator() (const MatrixSmall<3, 3>& F)
    {
        auto I1bar = Elasticity::I1bar (F);
        auto I2bar = Elasticity::I2bar (F);
//        auto Q = M_b * ( 0.75 - 0.5 * I1bar + 0.25 * I1bar * I1bar - 0.5 * I2bar );
        auto Q = 0.25 * M_b * ( I1bar * I1bar - 2.0 * I1bar - 2.0 * I2bar + 3.0  );
//        auto Q = 0.25 * M_b * ( I1bar * I1bar - 2.0 * I1bar + 2.0 * I2bar - 9.0  );
        auto dQdI2bar =  - 0.5 * M_b;
        return 0.5 * M_C * std::exp(Q) * dQdI2bar;
    }


    Real M_C;
    Real M_b;
};



template<class Mesh>
class dW2dI2
{
public:
    typedef typename MaterialFunctions::EMMaterialFunctions<Mesh>::return_Type return_Type;

    dW2dI2(Real C = 3330., Real b = 9.242) : M_C (C), M_b (b)  {}
    dW2dI2 (const dW2dI2& w1)
    {
        M_C = w1.M_C;
        M_b = w1.M_b;
    }
    virtual return_Type operator() (const MatrixSmall<3, 3>& F)
    {
        auto I1bar = Elasticity::I1bar (F);
        auto I2bar = Elasticity::I2bar (F);
//        auto Q = M_b * ( 0.75 - 0.5 * I1bar + 0.25 * I1bar * I1bar - 0.5 * I2bar );
        auto Q = 0.25 * M_b * ( I1bar * I1bar - 2.0 * I1bar - 2.0 * I2bar + 3.0  );
//        auto Q = 0.25 * M_b * ( I1bar * I1bar - 2.0 * I1bar + 2.0 * I2bar - 9.0  );

        auto dQdI2bar =  - 0.5 * M_b;
        return 0.5 * M_C * std::exp(Q) * dQdI2bar * dQdI2bar;
    }


    Real M_C;
    Real M_b;
};


template<class Mesh>
class dW2dI1
{
public:
    typedef typename MaterialFunctions::EMMaterialFunctions<Mesh>::return_Type return_Type;

    dW2dI1(Real C = 3330., Real b = 9.242) : M_C (C), M_b (b)  {}
    dW2dI1 (const dW2dI1& w1)
    {
        M_C = w1.M_C;
        M_b = w1.M_b;
    }
    virtual return_Type operator() (const MatrixSmall<3, 3>& F)
    {
        auto I1bar = Elasticity::I1bar (F);
        auto I2bar = Elasticity::I2bar (F);
//        auto Q = M_b * ( 0.75 - 0.5 * I1bar + 0.25 * I1bar * I1bar - 0.5 * I2bar );
        auto Q = 0.25 * M_b * ( I1bar * I1bar - 2.0 * I1bar - 2.0 * I2bar + 3.0  );
//        auto Q = 0.25 * M_b * ( I1bar * I1bar - 2.0 * I1bar + 2.0 * I2bar - 9.0  );
        auto dQdI1bar = 0.5 * M_b * (I1bar - 1.0);
        auto dQdI2bar =  - 0.5 * M_b;
//        Q = 0.0;

//        return 0.0;
        return 0.5 * M_C * std::exp(Q) * dQdI2bar * dQdI1bar;
    }


    Real M_C;
    Real M_b;
};


}

template <class Mesh>
class IsotropicFung : public virtual EMMaterialFunctions<Mesh>
{
public:
    typedef typename MaterialFunctions::EMMaterialFunctions<Mesh>::return_Type return_Type;



    //    IsotropicFung() : M_a(3330), M_b(9.242) {} // 0.33 KPa
    IsotropicFung (Real C = 3330., Real b = 9.242) : M_W1    ( new IsotropicFungFunctions::W1<Mesh>    (C,b) ),
    		                                         M_W2    ( new IsotropicFungFunctions::W2<Mesh>    (C,b) ),
    		                                         M_dW1dI1( new IsotropicFungFunctions::dW1dI1<Mesh>(C,b) ),
    		                                         M_dW2dI2( new IsotropicFungFunctions::dW2dI2<Mesh>(C,b) ),
    		                                         M_dW1dI2( new IsotropicFungFunctions::dW1dI2<Mesh>(C,b) ),
    		                                         M_dW2dI1( new IsotropicFungFunctions::dW2dI1<Mesh>(C,b) ) {} // 0.33 KPa
    IsotropicFung (const IsotropicFung& isotropicFung)
    {
        *M_W1 = *(isotropicFung.M_W1);
        *M_W2 = *(isotropicFung.M_W2);
    }
    virtual ~IsotropicFung() {}

    inline virtual void computeJacobian ( const vector_Type& disp,
                                          boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > >  dispETFESpace,
                                          const vector_Type& fibers,
                                          const vector_Type& sheets,
                                          matrixPtr_Type           jacobianPtr)
    {
        EMAssembler::computeI1JacobianTerms (disp, dispETFESpace, jacobianPtr, M_W1 );
        EMAssembler::computeI2JacobianTerms (disp, dispETFESpace, jacobianPtr, M_W2 );

        EMAssembler::computeI1JacobianTermsSecondDerivative(disp, dispETFESpace, jacobianPtr, M_dW1dI1 );
        EMAssembler::computeI2JacobianTermsSecondDerivative(disp, dispETFESpace, jacobianPtr, M_dW2dI2 );

        EMAssembler::computeI1JacobianMixedTermsSecondDerivative(disp, dispETFESpace, jacobianPtr, M_dW1dI2 );
        EMAssembler::computeI2JacobianMixedTermsSecondDerivative(disp, dispETFESpace, jacobianPtr, M_dW2dI1 );
    }
    inline virtual void computeResidual ( const vector_Type& disp,
                                          boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > >  dispETFESpace,
                                          const vector_Type& fibers,
                                          const vector_Type& sheets,
                                          vectorPtr_Type           residualVectorPtr)
    {
        EMAssembler::computeI1ResidualTerms (disp, dispETFESpace, residualVectorPtr, M_W1 );
        EMAssembler::computeI2ResidualTerms (disp, dispETFESpace, residualVectorPtr, M_W2 );
    }

    void showMe()
    {
        std::cout << "Isotropic Fung Function\n";
        std::cout << "Coefficient C: " << M_W1->M_C;
        std::cout << ", coefficient b: " << M_W1->M_b << "\n";
    }


    boost::shared_ptr< IsotropicFungFunctions::W1<Mesh> > M_W1;
    boost::shared_ptr< IsotropicFungFunctions::W2<Mesh> > M_W2;
    boost::shared_ptr< IsotropicFungFunctions::dW1dI1<Mesh> > M_dW1dI1;
    boost::shared_ptr< IsotropicFungFunctions::dW2dI2<Mesh> > M_dW2dI2;
    boost::shared_ptr< IsotropicFungFunctions::dW1dI2<Mesh> > M_dW1dI2;
    //this is equivalente to the M_dW1dI2, for performances use just one of them
    boost::shared_ptr< IsotropicFungFunctions::dW2dI1<Mesh> > M_dW2dI1;
};


} //EMMaterialFunctions

} //LifeV
#endif /* EMMATERIALFUNCTIONS_HPP_ */
