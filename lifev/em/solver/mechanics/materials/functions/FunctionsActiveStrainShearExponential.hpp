/*
 * FunctionActiveStrainShearExponential.hpp
 *
 *  Created on: 03/mar/2015
 *      Author: Luca Barbarotta
 */

#ifndef FUNCTION_ACTIVE_STRAIN_SHEAR_EXPONENTIAL_HPP_
#define FUNCTION_ACTIVE_STRAIN_SHEAR_EXPONENTIAL_HPP_

#include <lifev/em/solver/mechanics/materials/functions/FunctionsShearExponential.hpp>
#include <lifev/em/solver/mechanics/materials/functions/EMActiveStrainMaterialFunctions.hpp>

namespace LifeV
{
    namespace MaterialFunctions
    {

        template <class Mesh>
        class ActiveStrainShearExponential : public virtual ShearExponential<Mesh>,
                                             public virtual EMActiveStrainMaterialFunctions<Mesh>
        {
        public:
            typedef typename MaterialFunctions::EMMaterialFunctions<Mesh>::return_Type return_Type;
            typedef ShearExponential<Mesh> super;
            typedef typename super::data_Type data_Type;
            typedef EMActiveStrainMaterialFunctions<Mesh> activeStrain;

            ActiveStrainShearExponential (Real a = 4170., Real b = 11.602, bool NISplitting=true) : super(a,b), M_NISplitting(NISplitting) {} // 0.33 KPa
            ActiveStrainShearExponential (const ActiveStrainShearExponential& ActiveStrainShearExponential)
            {
                this->M_a = ActiveStrainShearExponential.M_a;
                this->M_b = ActiveStrainShearExponential.M_b;
                M_NISplitting = ActiveStrainShearExponential.M_NISplitting;
                this->M_activeStrainType = ActiveStrainShearExponential.M_activeStrainType;
                this->M_activeStrainOrthotropicParameter = ActiveStrainShearExponential.M_activeStrainOrthotropicParameter;
            }
            virtual ~ActiveStrainShearExponential() {}

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
                if ( M_NISplitting )
                {
                    EMAssembler::ActiveStrainNearlyIncompressible::computeActiveStrainI8barJacobianTerms (disp,
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
                    EMAssembler::ActiveStrain::computeActiveStrainI8JacobianTerms (disp,
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
                if ( M_NISplitting )
                {
                    EMAssembler::ActiveStrainNearlyIncompressible::computeActiveStrainI8barResidualTerms (disp,
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
                    EMAssembler::ActiveStrain::computeActiveStrainI8ResidualTerms (disp,
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
            
            void showMe()
            {
                std::cout << "Active Strain Shear Exponential Function\n";
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

        }; // End Class ActiveStrainShearExponential

        template <class Mesh>
        class dActiveStrainShearExponential : public virtual dShearExponential<Mesh>,
                                              public virtual EMActiveStrainMaterialFunctions<Mesh>
        {
        public:
            typedef typename MaterialFunctions::EMMaterialFunctions<Mesh>::return_Type return_Type;
            typedef dShearExponential<Mesh> super;
            typedef typename super::data_Type data_Type;
            typedef EMActiveStrainMaterialFunctions<Mesh> activeStrain;

            dActiveStrainShearExponential (Real a = 4170., Real b = 11.602, bool NISplitting = true) : super(a,b), M_NISplitting(NISplitting) {} // 0.33 KPa
            dActiveStrainShearExponential (const dActiveStrainShearExponential& dActiveStrainShearExponential)
            {
                this->M_a = dActiveStrainShearExponential.M_a;
                this->M_b = dActiveStrainShearExponential.M_b;
                M_NISplitting = dActiveStrainShearExponential.M_NISplitting;
                this->M_activeStrainType = dActiveStrainShearExponential.M_activeStrainType;
                this->M_activeStrainOrthotropicParameter = dActiveStrainShearExponential.M_activeStrainOrthotropicParameter;
            }
            virtual ~dActiveStrainShearExponential() {}

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
                if ( M_NISplitting )
                {
                    EMAssembler::ActiveStrainNearlyIncompressible::computeActiveStrainI8barJacobianTermsSecondDerivative
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
                    EMAssembler::ActiveStrain::computeActiveStrainI8JacobianTermsSecondDerivative (disp,
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

            void showMe()
            {
                std::cout << "Derivative Active Strain Shear Exponential Function\n";
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
            
        }; // End Class ActiveStrainShearExponential

    } // End Namespace MaterialFunctions
    
} // End Namespace LifeV

#endif 
