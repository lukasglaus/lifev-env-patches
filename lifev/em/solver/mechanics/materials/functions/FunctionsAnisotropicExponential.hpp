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

    enum Field { Fibers, Sheets };

    typedef typename MaterialFunctions::EMMaterialFunctions<Mesh>::return_Type return_Type;
    typedef EMData          data_Type;

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

//    class W4
//    {
//    public:
//        typedef typename MaterialFunctions::EMMaterialFunctions<Mesh>::return_Type return_Type;
//
//        virtual return_Type operator() (const VectorSmall<3>& f)
//        {
//
//            LifeV::Real I4 = f.dot (f);
//            LifeV::Real I4m1 = I4 - 1.0;
//            return M_a * I4m1 * std::exp (M_b * I4m1 * I4m1 ) * Elasticity::RegularizedHeaviside (I4m1);
//        }
//
//        virtual return_Type operator() (const LifeV::Real& I4)
//        {
//            LifeV::Real I4m1 = I4 - 1.0;
//            return M_a * I4m1 * std::exp (M_b * I4m1 * I4m1 ) * Elasticity::RegularizedHeaviside (I4m1);
//        }
//
//        W4(Real a = 185350, Real b = 15.972) : M_a(a), M_b(b) {}
//        W4(const W4& W4)
//        {
//        	M_a = W4.M_a;
//        	M_b = W4.M_b;
//        }
//
//        Real M_a;
//        Real M_b;
//
//    };


//    class dW4 : public W4
//    {
//    public:
//        typedef typename MaterialFunctions::EMMaterialFunctions<Mesh>::return_Type return_Type;
//        typedef W4 super;
//
//        return_Type operator() (const VectorSmall<3>& f)
//        {
//            LifeV::Real I4 = f.dot (f);
//            LifeV::Real I4m1 = I4 - 1.0;
//            return this->M_a * std::exp ( this->M_b * I4m1 * I4m1 )
//                   * ( 1.0 + 2.0 * this->M_b * I4m1 * I4m1 ) * Elasticity::RegularizedHeaviside (I4m1);
////                   + this->M_a * I4m1 * std::exp (this->M_b * I4m1 * I4m1 ) * Elasticity::dRegularizedHeaviside (I4m1);
//        }
//
//        return_Type operator() (const LifeV::Real& I4)
//        {
//            LifeV::Real I4m1 = I4 - 1.0;
//            return this->M_a * std::exp ( this->M_b * I4m1 * I4m1 )
//                   * ( 1.0 + 2.0 * this->M_b * I4m1 * I4m1 ) * Elasticity::RegularizedHeaviside (I4m1);
////                   + this->M_a * I4m1 * std::exp (this->M_b * I4m1 * I4m1 ) * Elasticity::dRegularizedHeaviside (I4m1);
//        }
//
//        dW4(Real a = 185350, Real b = 15.972) : super(a,b) {}
//        dW4(const dW4& dW4) : super(dW4) {}
//
//
//
//    };

    AnisotropicExponential (Real a = 185350, Real b = 15.972, Field fibers = Fibers ) :
        M_anisotropyField (fibers),
        M_a(a), M_b(b) {}
//        M_W4 ( new W4(a, b)),
//        M_dW4 ( new dW4(a, b)){} // 0.33 KPa
    AnisotropicExponential (const AnisotropicExponential& AnisotropicExponential)
    {
        M_a = AnisotropicExponential.M_a;
        M_b = AnisotropicExponential.M_b;
        M_anisotropyField = AnisotropicExponential.M_anisotropyField;
    }
    virtual ~AnisotropicExponential() {}

    virtual void computeJacobian ( const vector_Type& disp,
                                          boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > >  dispETFESpace,
                                          const vector_Type& fibers,
                                          const vector_Type& sheets,
                                          matrixPtr_Type     jacobianPtr)
    {
        if (M_anisotropyField == Fibers)
        {
            EMAssembler::computeI4JacobianTerms (disp, dispETFESpace, fibers, jacobianPtr, this->getMe() );
//            EMAssembler::computeI4JacobianTermsSecondDerivative (disp, dispETFESpace, fibers, jacobianPtr, M_dW4 );
        }
        else
        {
            EMAssembler::computeI4JacobianTerms (disp, dispETFESpace, fibers,  sheets, jacobianPtr, this->getMe() );
//            EMAssembler::computeI4JacobianTermsSecondDerivative (disp, dispETFESpace, fibers, sheets, jacobianPtr, M_dW4 );
        }
    }

    virtual void computeResidual ( const vector_Type& disp,
                                          boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > >  dispETFESpace,
                                          const vector_Type& fibers,
                                          const vector_Type& sheets,
                                          vectorPtr_Type     residualVectorPtr)
    {
        if (M_anisotropyField == Fibers)
        {
            EMAssembler::computeI4ResidualTerms (disp, dispETFESpace, fibers, residualVectorPtr, this->getMe());
        }
        else
        {
            EMAssembler::computeI4ResidualTerms (disp, dispETFESpace, fibers, sheets, residualVectorPtr, this->getMe() );
        }
    }

    void showMe()
    {
        std::cout << "Anisotropic Exponential Function\n";
//        std::cout << "Coefficient a: " << M_W4  -> M_a;
//        std::cout << ", coefficient b: " << M_W4  -> M_b;
        std::cout << "Coefficient a: " << M_a;
        std::cout << ", coefficient b: " << M_b;
        if(M_anisotropyField == Fibers)
        	std::cout << ". Using fibers.\n";
        else
        	std::cout << ". Using sheets.\n";
    }

    void setParametersFromGetPot (GetPot& data)
    {
//        M_a = data ( "solid/physics/af", 185350);
//        M_b = data ( "solid/physics/bf", 15.972);
    }

    void setParameters (data_Type& data)
    {
        if(M_anisotropyField == Fibers)
		{
//        	M_W4  -> M_a = data.solidParameter<Real>("af");
//        	M_W4  -> M_b = data.solidParameter<Real>("bf");
//        	M_dW4 -> M_a = data.solidParameter<Real>("af");
//        	M_dW4 -> M_b = data.solidParameter<Real>("bf");
        	M_a = data.solidParameter<Real>("af");
        	M_b = data.solidParameter<Real>("bf");
		}
        else
		{
//        	M_W4  -> M_a = data.solidParameter<Real>("as");
//        	M_W4  -> M_b = data.solidParameter<Real>("bs");
//        	M_dW4 -> M_a = data.solidParameter<Real>("as");
//        	M_dW4 -> M_b = data.solidParameter<Real>("bs");
        	M_a = data.solidParameter<Real>("as");
        	M_b = data.solidParameter<Real>("bs");
		}

    }

protected:
    Field M_anisotropyField; //use fiber vector if true,  else uses sheets
    Real M_a;
    Real M_b;
//    boost::shared_ptr< W4 >   M_W4;
//    boost::shared_ptr< dW4 >   M_dW4;


};






template <class Mesh>
class dAnisotropicExponential : public virtual AnisotropicExponential<Mesh>
{
public:
    typedef typename MaterialFunctions::EMMaterialFunctions<Mesh>::return_Type return_Type;
    typedef AnisotropicExponential<Mesh>  super;

    virtual return_Type operator() (const VectorSmall<3>& f)
    {
        LifeV::Real I4 = f.dot (f);
        LifeV::Real I4m1 = I4 - 1.0;
        return this->M_a * std::exp ( this->M_b * I4m1 * I4m1 )
               * ( 1.0 + 2.0 * this->M_b * I4m1 * I4m1 ) * Elasticity::RegularizedHeaviside (I4m1);
//               + this->M_a * I4m1 * std::exp (this->M_b * I4m1 * I4m1 ) * Elasticity::dRegularizedHeaviside (I4m1);
    }

    virtual return_Type operator() (const LifeV::Real& I4)
    {
        LifeV::Real I4m1 = I4 - 1.0;
        return this->M_a * std::exp ( this->M_b * I4m1 * I4m1 )
               * ( 1.0 + 2.0 * this->M_b * I4m1 * I4m1 ) * Elasticity::RegularizedHeaviside (I4m1);
//               + this->M_a * I4m1 * std::exp (this->M_b * I4m1 * I4m1 ) * Elasticity::dRegularizedHeaviside (I4m1);
    }

    //    dAnisotropicExponential() : M_a(3330), M_b(9.242) {} // 0.33 KPa
    dAnisotropicExponential (Real a = 185350, Real b = 15.972, typename super::Field fibers = super::Fibers ) : super(a, b, fibers)
    {} // 0.33 KPa
    dAnisotropicExponential (const dAnisotropicExponential& dAnisotropicExponential)
    {
    	this->M_a = dAnisotropicExponential.super::M_a;
    	this->M_b = dAnisotropicExponential.super::M_b;
    	this->M_anisotropyField = dAnisotropicExponential.super::M_anisotropyField;
    }
    virtual ~dAnisotropicExponential() {}

    virtual void computeResidual ( const vector_Type& disp,
                                          boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > >  dispETFESpace,
                                          const vector_Type& fibers,
                                          const vector_Type& sheets,
                                          vectorPtr_Type     residualVectorPtr)
    {

    }

    virtual void computeJacobian ( const vector_Type& disp,
                                          boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > > dispETFESpace,
                                          const vector_Type& fibers,
                                          const vector_Type& sheets,
                                          matrixPtr_Type           jacobianPtr)
    {
        if (this->M_anisotropyField == super::Fibers)
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
        std::cout << "Coefficient a: " << this->M_a;
        std::cout << ", coefficient b: " << this->M_b;
        if(this->M_anisotropyField == super::Fibers)
        	std::cout << ". Using fibers.\n";
        else
        	std::cout << ". Using sheets.\n";
    }


};


} //EMMaterialFunctions

} //LifeV
#endif /* EMMATERIALFUNCTIONS_HPP_ */
