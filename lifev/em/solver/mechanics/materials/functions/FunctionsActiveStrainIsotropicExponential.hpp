/*
 * FunctionActiveStrainIsotropicExponential.hpp
 *
 *  Created on: 28/apr/2014
 *      Author: srossi, Luca Barbarotta
 */

#ifndef FUNCTIONSACTIVESTRAINISOTROPICEXPONENTIAL_HPP_
#define FUNCTIONSACTIVESTRAINISOTROPICEXPONENTIAL_HPP_

#include <lifev/em/solver/mechanics/EMElasticityFunctions.hpp>

//#include <lifev/em/solver/mechanics/EMETAAssembler.hpp>
#include <lifev/em/solver/mechanics/materials/functions/FunctionsIsotropicExponential.hpp>
#include <lifev/em/solver/mechanics/materials/functions/EMActiveStrainMaterialFunctions.hpp>

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
class ActiveStrainIsotropicExponential : public virtual IsotropicExponential<Mesh>,
                                         public virtual EMActiveStrainMaterialFunctions<Mesh>
{
public:
    typedef typename MaterialFunctions::EMMaterialFunctions<Mesh>::return_Type return_Type;
    typedef IsotropicExponential<Mesh> super;
    typedef typename super::data_Type data_Type;
    typedef EMActiveStrainMaterialFunctions<Mesh> activeStrain;


    virtual return_Type operator() (const MatrixSmall<3, 3>& F)
    {
        auto I1bar = Elasticity::I1bar (F);
        return this->M_a / 2.0 * std::exp ( this->M_b * ( I1bar - 3 ) );
    }

    //    IsotropicExponential() : M_a(3330), M_b(9.242) {} // 0.33 KPa
    ActiveStrainIsotropicExponential (Real a = 3330., Real b = 9.242) : super(a,b) {} // 0.33 KPa
    ActiveStrainIsotropicExponential (const ActiveStrainIsotropicExponential& ActiveStrainIsotropicExponential)
    {
        this->M_a = ActiveStrainIsotropicExponential.M_a;
        this->M_b = ActiveStrainIsotropicExponential.M_b;
        this->M_activeStrainType = ActiveStrainIsotropicExponential.M_activeStrainType;
        this->M_activeStrainOrthotropicParameter = ActiveStrainIsotropicExponential.M_activeStrainOrthotropicParameter;
    }
    virtual ~ActiveStrainIsotropicExponential() {}
    
    virtual  void computeJacobian ( const vector_Type& disp,
                                    boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > >  dispETFESpace,
                                    const vector_Type& fibers,
                                    const vector_Type& sheets,
                                    const vectorPtr_Type& fiberActivation,
                                    const vectorPtr_Type& sheetActivation,
                                    const vectorPtr_Type& normalActivation,
                                    boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 1 > >  activationETFESpace,
                                    matrixPtr_Type           jacobianPtr)
    {
        // EMAssembler::computeI1JacobianTerms (disp, dispETFESpace, jacobianPtr, this->getMe() );
        EMAssembler::ActiveStrainNearlyIncompressible::computeActiveStrainI1barJacobianTerms (disp,
                                                                                              dispETFESpace,
                                                                                              fibers,
                                                                                              sheets,
                                                                                              fiberActivation,
                                                                                              sheetActivation,
                                                                                              normalActivation,
                                                                                              activationETFESpace,
                                                                                              jacobianPtr,
                                                                                              this->getMe(),
                                                                                              this->M_activeStrainOrthotropicParameter);

    }

    virtual void computeResidual ( const vector_Type& disp,
                                    boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > >  dispETFESpace,
                                    const vector_Type& fibers,
                                    const vector_Type& sheets,
                                    const vectorPtr_Type& fiberActivation,
                                    const vectorPtr_Type& sheetActivation,
                                    const vectorPtr_Type& normalActivation,
                                    boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 1 > >  activationETFESpace,
                                   vectorPtr_Type           residualVectorPtr)
    {
        EMAssembler::ActiveStrainNearlyIncompressible::computeActiveStrainI1barResidualTerms( disp,
                                                                                              dispETFESpace,
                                                                                              fibers,
                                                                                              sheets,
                                                                                              fiberActivation,
                                                                                              sheetActivation,
                                                                                              normalActivation,
                                                                                              activationETFESpace,
                                                                                              residualVectorPtr,
                                                                                              this->getMe(),
                                                                                              this->M_activeStrainOrthotropicParameter);
        
        // EMAssembler::computeI1ResidualTerms (disp, dispETFESpace, residualVectorPtr, this->getMe() );
    }

    void showMe()
    {
        std::cout << "Active Strain Isotropic Exponential Function\n";
        std::cout << "Coefficient a: " << this->M_a;
        std::cout << ", coefficient b: " << this->M_b << "\n";
    }

    virtual void setParameters (data_Type& data)
    {
    	activeStrain::setParameters(data);
    	super::setParameters(data);
    }

private:

};

template <class Mesh>
class dActiveStrainIsotropicExponential : public virtual IsotropicExponential<Mesh>,
                                          public virtual EMActiveStrainMaterialFunctions<Mesh>
{
public:
    typedef typename MaterialFunctions::EMMaterialFunctions<Mesh>::return_Type return_Type;
    typedef IsotropicExponential<Mesh> super;
    typedef typename super::data_Type data_Type;
    typedef EMActiveStrainMaterialFunctions<Mesh> activeStrain;

    dActiveStrainIsotropicExponential (Real a = 3330., Real b = 9.242) : super(a,b) {} // 0.33 KPa
    dActiveStrainIsotropicExponential (const dActiveStrainIsotropicExponential& dActiveStrainIsotropicExponential)
    {
        this->M_a = dActiveStrainIsotropicExponential.M_a;
        this->M_b = dActiveStrainIsotropicExponential.M_b;
        this->M_activeStrainType = dActiveStrainIsotropicExponential.M_activeStrainType;
        this->M_activeStrainOrthotropicParameter = dActiveStrainIsotropicExponential.M_activeStrainOrthotropicParameter;
    }

    virtual return_Type operator() (const MatrixSmall<3, 3>& F)
    {
        auto I1bar = Elasticity::I1bar (F);
        return this->M_a * this->M_b / 2.0 * std::exp ( this->M_b * ( I1bar - 3 ) );
    }

    virtual ~dActiveStrainIsotropicExponential() {}


    virtual  void computeJacobian ( const vector_Type& disp,
                                    boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > >  dispETFESpace,
                                    const vector_Type& fibers,
                                    const vector_Type& sheets,
                                    const vectorPtr_Type& fiberActivation,
                                    const vectorPtr_Type& sheetActivation,
                                    const vectorPtr_Type& normalActivation,
                                    boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 1 > >  activationETFESpace,
                                    matrixPtr_Type           jacobianPtr)
    
    {
        EMAssembler::ActiveStrainNearlyIncompressible::computeActiveStrainI1barJacobianTermsSecondDerivative (disp,
                                                                                                              dispETFESpace,
                                                                                                              fibers,
                                                                                                              sheets,
                                                                                                              fiberActivation,
                                                                                                              sheetActivation,
                                                                                                              normalActivation,
                                                                                                              activationETFESpace,
                                                                                                              jacobianPtr,
                                                                                                              this->getMe(),
                                                                                                              this->M_activeStrainOrthotropicParameter);
        
    }

    void showMe()
    {
        std::cout << "Derivative Active Strain Isotropic Exponential Function\n";
        std::cout << "Coefficient a: " << this->M_a;
        std::cout << ", coefficient b: " << this->M_b << "\n";
    }

    virtual void setParameters (data_Type& data)
    {
    	activeStrain::setParameters(data);
    	super::setParameters(data);
    }

private:

};


} //MaterialFunctions

} //LifeV
#endif /* FUNCTIONSACTIVESTRAINISOTROPICEXPONENTIAL_HPP_ */
