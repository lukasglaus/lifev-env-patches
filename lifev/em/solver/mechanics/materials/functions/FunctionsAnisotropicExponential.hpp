/*
 * EMMaterialFunctions.hpp
 *
 *  Created on: 28/apr/2014
 *      Author: srossi
 */

#ifndef FUNCTIONSAnisotropicExponential_HPP_
#define FUNCTIONSAnisotropicExponential_HPP_

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
//  ISOTROPIC EXPONENTIAL FUNCTIONS
////////////////////////////////////////////////////////////////////////

template <class Mesh>
class AnisotropicExponential : public virtual EMMaterialFunctions<Mesh>
{
public:
    typedef typename MaterialFunctions::EMMaterialFunctions<Mesh>::return_Type return_Type;

    virtual return_Type operator() (const VectorSmall<3>& f)
    {
        LifeV::Real I4 = f.dot (f);
        LifeV::Real I4m1 = I4 - 1.0;
        return M_a * I4m1 * std::exp (M_b * I4m1 * I4m1 ) * Elasticity::RegularizedHeaviside (I4m1);
    }

    virtual return_Type operator() (const LifeV::Real& I4)
    {
        LifeV::Real I4m1 = I4 - 1.0;
        return M_a * I4m1 * std::exp (M_b * I4m1 * I4m1 ) * Elasticity::RegularizedHeaviside (I4m1);
    }


    AnisotropicExponential (Real a = 185350, Real b = 15.972, bool useFibers = true) : M_a (a),
        M_b (b),
        M_useFibers (useFibers) {} // 0.33 KPa
    AnisotropicExponential (const AnisotropicExponential& AnisotropicExponential)
    {
        M_a = AnisotropicExponential.M_a;
        M_b = AnisotropicExponential.M_b;
        M_useFibers = AnisotropicExponential.M_useFibers;
    }
    virtual ~AnisotropicExponential() {}

    inline virtual void computeJacobian ( const vector_Type& disp,
                                          boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > >  dispETFESpace,
                                          const vector_Type& fibers,
                                          const vector_Type& sheets,
                                          matrixPtr_Type     jacobianPtr)
    {
        if (M_useFibers)
        {
            EMAssembler::computeI4JacobianTerms (disp, dispETFESpace, fibers, jacobianPtr, this->getMe() );
        }
        else
        {
            EMAssembler::computeI4JacobianTerms (disp, dispETFESpace, fibers,  sheets, jacobianPtr, this->getMe() );
        }
    }

    inline virtual void computeResidual ( const vector_Type& disp,
                                          boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > >  dispETFESpace,
                                          const vector_Type& fibers,
                                          const vector_Type& sheets,
                                          vectorPtr_Type     residualVectorPtr)
    {
        if (M_useFibers)
        {
            EMAssembler::computeI4ResidualTerms (disp, dispETFESpace, fibers, residualVectorPtr, this->getMe() );
        }
        else
        {
            EMAssembler::computeI4ResidualTerms (disp, dispETFESpace, fibers, sheets, residualVectorPtr, this->getMe() );
        }
    }

    void showMe()
    {
        std::cout << "Anisotropic Exponential Function\n";
        std::cout << "Coefficient a: " << M_a;
        std::cout << ", coefficient b: " << M_b;
        std::cout << ". Use Fibers? "  << M_useFibers << "\n";
    }

    void setParametersFromGetPot (GetPot& data)
    {
        M_a = data ( "solid/physics/af", 185350);
        M_b = data ( "solid/physics/bf", 15.972);
        M_useFibers = data ( "solid/physics/fibers", true);
    }

private:
    Real M_a;
    Real M_b;
    bool M_useFibers; //use fiber vector if true,  else uses sheets

};

template <class Mesh>
class dAnisotropicExponential : public virtual EMMaterialFunctions<Mesh>
{
public:
    typedef typename MaterialFunctions::EMMaterialFunctions<Mesh>::return_Type return_Type;

    virtual return_Type operator() (const VectorSmall<3>& f)
    {
        LifeV::Real I4 = f.dot (f);
        LifeV::Real I4m1 = I4 - 1.0;
        return M_a * std::exp ( M_b * I4m1 * I4m1 )
               * ( 1.0 + 2.0 * M_b * I4m1 * I4m1 ) * Elasticity::RegularizedHeaviside (I4m1)
               + M_a * I4m1 * std::exp (M_b * I4m1 * I4m1 ) * Elasticity::dRegularizedHeaviside (I4m1);
    }

    virtual return_Type operator() (const LifeV::Real& I4)
    {
        LifeV::Real I4m1 = I4 - 1.0;
        return M_a * std::exp ( M_b * I4m1 * I4m1 )
               * ( 1.0 + 2.0 * M_b * I4m1 * I4m1 ) * Elasticity::RegularizedHeaviside (I4m1) // * Elasticity::RegularizedHeaviside(I4-1);
               + M_a * I4m1 * std::exp (M_b * I4m1 * I4m1 ) * Elasticity::dRegularizedHeaviside (I4m1);
    }

    //    dAnisotropicExponential() : M_a(3330), M_b(9.242) {} // 0.33 KPa
    dAnisotropicExponential (Real a = 185350, Real b = 15.972, bool useFibers = true) : M_a (a),
        M_b (b),
        M_useFibers (useFibers) {} // 0.33 KPa
    dAnisotropicExponential (const dAnisotropicExponential& dAnisotropicExponential)
    {
        M_a = dAnisotropicExponential.M_a;
        M_b = dAnisotropicExponential.M_b;
        M_useFibers = dAnisotropicExponential.M_useFibers;
    }
    virtual ~dAnisotropicExponential() {}


    inline virtual void computeJacobian ( const vector_Type& disp,
                                          boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > > dispETFESpace,
                                          const vector_Type& fibers,
                                          const vector_Type& sheets,
                                          matrixPtr_Type           jacobianPtr)
    {
        if (M_useFibers)
        {
            EMAssembler::computeI4JacobianTermsSecondDerivative (disp, dispETFESpace, fibers, jacobianPtr, this->getMe() );
        }
        else
        {
            EMAssembler::computeI4JacobianTermsSecondDerivative (disp, dispETFESpace, fibers, sheets, jacobianPtr, this->getMe() );
        }
    }

    void showMe()
    {
        std::cout << "Derivative Anisotropic Exponential Function\n";
        std::cout << "Coefficient a: " << M_a;
        std::cout << ", coefficient b: " << M_b;
        std::cout << ". Use Fibers? "  << M_useFibers << "\n";
    }

private:
    Real M_a;
    Real M_b;
    bool M_useFibers;
};


} //EMMaterialFunctions

} //LifeV
#endif /* EMMATERIALFUNCTIONS_HPP_ */
