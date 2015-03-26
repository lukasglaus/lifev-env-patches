/*
 * EMMaterialFunctions.hpp
 *
 *  Created on: 28/apr/2014
 *      Author: srossi
 */

#ifndef FUNCTIONSISOTROPICEXPONENTIAL_HPP_
#define FUNCTIONSISOTROPICEXPONENTIAL_HPP_

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
//  FIBER EXPONENTIAL FUNCTIONS
////////////////////////////////////////////////////////////////////////

template <class Mesh>
class IsotropicExponential : public virtual EMMaterialFunctions<Mesh>
{
public:
    typedef typename MaterialFunctions::EMMaterialFunctions<Mesh>::return_Type return_Type;
    typedef EMData          data_Type;


    virtual return_Type operator() (const MatrixSmall<3, 3>& F)
    {
        auto I1bar = Elasticity::I1bar (F);
        return 0.5 * M_a * std::exp ( M_b * ( I1bar - 3 ) );
    }

    //    IsotropicExponential() : M_a(3330), M_b(9.242) {} // 0.33 KPa
    IsotropicExponential (Real a = 3330., Real b = 9.242) : M_a (a), M_b (b) {} // 0.33 KPa
    IsotropicExponential (const IsotropicExponential& isotropicExponential)
    {
        M_a = isotropicExponential.M_a;
        M_b = isotropicExponential.M_b;
    }
    virtual ~IsotropicExponential() {}

    virtual void computeJacobian ( const vector_Type& disp,
                                          boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > >  dispETFESpace,
                                          const vector_Type& fibers,
                                          const vector_Type& sheets,
                                          matrixPtr_Type           jacobianPtr)
    {
        EMAssembler::computeI1JacobianTerms (disp, dispETFESpace, jacobianPtr, this->getMe() );
    }

    virtual void computeResidual ( const vector_Type& disp,
                                          boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > >  dispETFESpace,
                                          const vector_Type& fibers,
                                          const vector_Type& sheets,
                                          vectorPtr_Type           residualVectorPtr)
    {
        EMAssembler::computeI1ResidualTerms (disp, dispETFESpace, residualVectorPtr, this->getMe() );
    }

    void showMe()
    {
        std::cout << "Isotropic Exponential Function\n";
        std::cout << "Coefficient a: " << M_a;
        std::cout << ", coefficient b: " << M_b << "\n";
    }

    void setParameters (data_Type& data)
    {
    	M_a = data.solidParameter<Real>("a");
    	M_b = data.solidParameter<Real>("b");
    }


protected:
    Real M_a;
    Real M_b;
};

template <class Mesh>
class dIsotropicExponential : public virtual EMMaterialFunctions<Mesh>
{
public:
    typedef typename MaterialFunctions::EMMaterialFunctions<Mesh>::return_Type return_Type;

    virtual return_Type operator() (const MatrixSmall<3, 3>& F)
    {
        auto I1bar = Elasticity::I1bar (F);
        return 0.5 * M_a * M_b * std::exp ( M_b * ( I1bar - 3 ) );
    }

    //    dIsotropicExponential() : M_a(3330), M_b(9.242) {} // 0.33 KPa
    dIsotropicExponential (Real a = 3330., Real b = 9.242) : M_a (a), M_b (b) {} // 0.33 KPa
    dIsotropicExponential (const dIsotropicExponential& dIsotropicExponential)
    {
        M_a = dIsotropicExponential.M_a;
        M_b = dIsotropicExponential.M_b;
    }
    virtual ~dIsotropicExponential() {}

    inline virtual void computeJacobian ( const vector_Type& disp,
                                          boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > >  dispETFESpace,
                                          const vector_Type& fibers,
                                          const vector_Type& sheets,
                                          matrixPtr_Type           jacobianPtr)
    {
        EMAssembler::computeI1JacobianTermsSecondDerivative (disp, dispETFESpace, jacobianPtr, this->getMe() );
    }

    void showMe()
    {
        std::cout << "Derivative Isotropic Exponential Function\n";
        std::cout << "Coefficient a: " << M_a;
        std::cout << ", coefficient b: " << M_b << "\n";
    }

    typedef EMData          data_Type;

    void setParameters (data_Type& data)
    {
    	M_a = data.solidParameter<Real>("a");
    	M_b = data.solidParameter<Real>("b");
    }

protected:
    Real M_a;
    Real M_b;
};


} //EMMaterialFunctions

} //LifeV
#endif /* EMMATERIALFUNCTIONS_HPP_ */
