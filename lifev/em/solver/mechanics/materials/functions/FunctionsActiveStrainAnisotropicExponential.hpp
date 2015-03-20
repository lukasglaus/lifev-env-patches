/*
 * FunctionActiveStrainAnisotropicExponential.hpp
 *
 *  Created on: 3/mar/2015
 *      Author: srossi, Luca Barbarotta
 */

#ifndef FUNCTIONS_ACTIVE_STRAIN_ANISOTROPIC_EXPONENTIAL_HPP_
#define FUNCTIONS_ACTIVE_STRAIN_ANISOTROPIC_EXPONENTIAL_HPP_

#include <lifev/em/solver/mechanics/EMElasticityFunctions.hpp>

//#include <lifev/em/solver/mechanics/EMETAAssembler.hpp>
#include <lifev/em/solver/mechanics/materials/functions/FunctionsAnisotropicExponential.hpp>
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
class ActiveStrainAnisotropicExponential : public virtual AnisotropicExponential<Mesh>,
                                           public virtual EMActiveStrainMaterialFunctions<Mesh>
{
public:
    typedef typename MaterialFunctions::EMMaterialFunctions<Mesh>::return_Type return_Type;
    typedef AnisotropicExponential<Mesh> super;
    typedef typename super::data_Type data_Type;
    typedef EMActiveStrainMaterialFunctions<Mesh> activeStrain;
    typedef typename super::Field Field;
    

    // virtual return_Type operator() (const VectorSmall<3>& f)
    // {
    //     return this->super::operator(f);
    // }

    // virtual return_Type operator() (const LifeV::Real& I4)
    // {
    //     return this->super::operator(I4);
    // }

    //    AnisotropicExponential() : M_a(3330), M_b(9.242) {} // 0.33 KPa
    ActiveStrainAnisotropicExponential (Real a = 185350., Real b = 15.972, Field fibers = Field::Fibers, bool NISplitting=true) : super(a,b,fibers), M_NISplitting(NISplitting) {} // 0.33 KPa
    ActiveStrainAnisotropicExponential (const ActiveStrainAnisotropicExponential& ActiveStrainAnisotropicExponential)
    {
        this->M_a = ActiveStrainAnisotropicExponential.M_a;
        this->M_b = ActiveStrainAnisotropicExponential.M_b;
        M_NISplitting = ActiveStrainAnisotropicExponential.M_NISplitting;
        this->M_anisotropyField = ActiveStrainAnisotropicExponential.M_anisotropyField;
        this->M_activeStrainType = ActiveStrainAnisotropicExponential.M_activeStrainType;
        this->M_activeStrainOrthotropicParameter = ActiveStrainAnisotropicExponential.M_activeStrainOrthotropicParameter;
    }
    virtual ~ActiveStrainAnisotropicExponential() {}
    
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
        if (this->M_anisotropyField == Field::Fibers)
        {
            if ( M_NISplitting )
            {
                EMAssembler::ActiveStrainNearlyIncompressible::computeActiveStrainI4barFibersJacobianTerms
                    (disp,
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
            else
            {
                EMAssembler::ActiveStrain::computeActiveStrainI4FibersJacobianTerms (disp,
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
        }
        else
        {
            if ( M_NISplitting )
            {
                EMAssembler::ActiveStrainNearlyIncompressible::computeActiveStrainI4barSheetsJacobianTerms
                    (disp,
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
            else
            {
                EMAssembler::ActiveStrain::computeActiveStrainI4SheetsJacobianTerms (disp,
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
        }

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
        if (this->M_anisotropyField == Field::Fibers)
        {
            if ( M_NISplitting )
            {
                EMAssembler::ActiveStrainNearlyIncompressible::computeActiveStrainI4barFibersResidualTerms
                    (disp,
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
            }
            else
            {
                EMAssembler::ActiveStrain::computeActiveStrainI4FibersResidualTerms (disp,
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
            }
        }
        else
        {
            if ( M_NISplitting )
            {
                EMAssembler::ActiveStrainNearlyIncompressible::computeActiveStrainI4barSheetsResidualTerms
                    (disp,
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
            }
            else
            {
                EMAssembler::ActiveStrain::computeActiveStrainI4SheetsResidualTerms (disp,
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
            }
        }
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
    // Used to include the Nearly-incompressible splitting inside the
    // active strain formulation
    bool const M_NISplitting;
    
};

template <class Mesh>
class dActiveStrainAnisotropicExponential : public virtual AnisotropicExponential<Mesh>,
                                            public virtual EMActiveStrainMaterialFunctions<Mesh>
{
public:
    typedef typename MaterialFunctions::EMMaterialFunctions<Mesh>::return_Type return_Type;
    typedef AnisotropicExponential<Mesh> super;
    typedef typename super::data_Type data_Type;
    typedef EMActiveStrainMaterialFunctions<Mesh> activeStrain;
    typedef typename super::Field Field;
    
    dActiveStrainAnisotropicExponential (Real a = 185350., Real b = 15.972, Field fibers = Field::Fibers, bool NISplitting=true) : super(a,b,fibers), M_NISplitting(NISplitting) {} // 0.33 KPa
    dActiveStrainAnisotropicExponential (const dActiveStrainAnisotropicExponential& dActiveStrainAnisotropicExponential)
    {
        this->M_a = dActiveStrainAnisotropicExponential.M_a;
        this->M_b = dActiveStrainAnisotropicExponential.M_b;
        M_NISplitting = dActiveStrainAnisotropicExponential.M_NISplitting;
        this->M_anisotropyField = dActiveStrainAnisotropicExponential.M_anisotropyField;
        this->M_activeStrainType = dActiveStrainAnisotropicExponential.M_activeStrainType;
        this->M_activeStrainOrthotropicParameter = dActiveStrainAnisotropicExponential.M_activeStrainOrthotropicParameter;
    }
    
    virtual ~dActiveStrainAnisotropicExponential() {}

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
        if (this->M_anisotropyField == Field::Fibers)
        {
            if ( M_NISplitting )
            {
                EMAssembler::ActiveStrainNearlyIncompressible::computeActiveStrainI4barFibersJacobianTermsSecondDerivative
                    (disp,
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
            else
            {
                EMAssembler::ActiveStrain::computeActiveStrainI4FibersJacobianTermsSecondDerivative (disp,
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
            
        }
        else
        {
            if ( M_NISplitting )
            {
                EMAssembler::ActiveStrainNearlyIncompressible::computeActiveStrainI4barSheetsJacobianTermsSecondDerivative
                    (disp,
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
            else
            {
                EMAssembler::ActiveStrain::computeActiveStrainI4SheetsJacobianTermsSecondDerivative (disp,
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
        }
        
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
    // Used to include the Nearly-incompressible splitting inside the
    // active strain formulation
    bool const M_NISplitting;
};


} //MaterialFunctions

} //LifeV
#endif /* FUNCTIONS_ACTIVE_STRAIN_ANISOTROPIC_EXPONENTIAL_HPP_ */
