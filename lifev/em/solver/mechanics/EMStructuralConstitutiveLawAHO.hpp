//@HEADER
/*
*******************************************************************************

    Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
    Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University

    This file is part of LifeV.

    LifeV is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    LifeV is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with LifeV.  If not, see <http://www.gnu.org/licenses/>.

*******************************************************************************
*/
//@HEADER

/*!
 *  @file
 *  @brief This file contains an abstract class to implement different kinds of materials for structural dynamic problems (St. Venant-Kirchhoff, Neo-Hookean and Exponential materials right now )
 *
 *  @version 1.0
 *  @date 04-2014
 *  @author Simone Rossi <simone.rossi@epfl.ch>
 */

#ifndef _EMSTRUCTURALCONSTITUTIVELAW_H_
#define _EMSTRUCTURALCONSTITUTIVELAW_H_ 1

#include <boost/typeof/typeof.hpp>
#include <lifev/core/array/MatrixSmall.hpp>

#include <lifev/core/filter/ExporterHDF5.hpp>
#include <lifev/structure/solver/StructuralConstitutiveLaw.hpp>
//#include <lifev/em/solver/mechanics/materials/EMMaterial.hpp>
//#include <lifev/em/solver/mechanics/EMStructuralConstitutiveLaw.hpp>
#include <lifev/em/solver/mechanics/materials/EMMaterialType.hpp>
//#include <lifev/em/solver/mechanics/materials/EMActiveStressMaterial.hpp>
//#include <lifev/em/solver/mechanics/materials/EMActiveStrainMaterial.hpp>

//#include <lifev/em/util/EMUtility.hpp>
#include <lifev/electrophysiology/util/HeartUtility.hpp>
#include <lifev/em/solver/mechanics/materials/MaterialsList.hpp>
//#include <lifev/em/solver/mechanics/materials/functions/FunctionsList.hpp>

#include <lifev/em/solver/mechanics/EMMechanicalExpressions.hpp>
#include <lifev/em/solver/EMETAFunctors.hpp>



namespace LifeV
{
/*!
  \class EMStructuralConstitutiveLaw
  \brief
  This class is an abstract class to define different type of models for the arterial wall.
  This class has just pure virtual methods. They are implemented in the specific class for one material
*/

//forward declaration
class EMData;

template <typename MeshType>
class EMStructuralConstitutiveLaw  : public StructuralConstitutiveLaw<MeshType>
{
public:

    //!@name Type definitions
    //@{
    typedef typename boost::shared_ptr<const Displayer>   displayerPtr_Type;
    typedef StructuralConstitutiveLaw<MeshType> super;
    typedef MatrixEpetra<Real>            matrix_Type;
    typedef boost::shared_ptr<matrix_Type>         matrixPtr_Type;
    typedef VectorEpetra           vector_Type;
    typedef boost::shared_ptr<vector_Type>         vectorPtr_Type;

    typedef ETFESpace<MeshType, MapEpetra, 3, 3 >         ETFESpace_Type;
    typedef boost::shared_ptr<ETFESpace_Type>             ETFESpacePtr_Type;

    typedef FESpace< MeshType, MapEpetra >                FESpace_Type;
    typedef boost::shared_ptr<FESpace_Type>               FESpacePtr_Type;

    typedef MeshType                                        mesh_Type;
    typedef ETFESpace< mesh_Type, MapEpetra, 3, 1 >                        scalarETFESpace_Type;
    typedef boost::shared_ptr<ETFESpace< mesh_Type, MapEpetra, 3, 1 > >    scalarETFESpacePtr_Type;

    typedef MapEpetra map_Type;
    typedef boost::shared_ptr<map_Type>             mapPtr_Type;

    //    typedef EMMaterial<MeshType>                              material_Type;

    typedef StructuralConstitutiveLawData          data_Type;
    typedef typename boost::shared_ptr<data_Type>  dataPtr_Type;

    //    typedef FactorySingleton<Factory<StructuralConstitutiveLaw<MeshType>, std::string> >  StructureMaterialFactory;

    typedef std::vector< typename MeshType::element_Type* > vectorVolumes_Type;

    typedef std::map< UInt, vectorVolumes_Type>           mapMarkerVolumes_Type;
    typedef boost::shared_ptr<mapMarkerVolumes_Type>      mapMarkerVolumesPtr_Type;

    typedef std::vector<UInt>                             vectorIndexes_Type;
    typedef std::map< UInt, vectorIndexes_Type>           mapMarkerIndexes_Type;
    typedef boost::shared_ptr<mapMarkerIndexes_Type>      mapMarkerIndexesPtr_Type;

    typedef EMMaterialType<MeshType>                      material_Type;
    typedef boost::shared_ptr<material_Type>              materialPtr_Type;

    typedef EMPassiveMaterialType<MeshType>               passiveMaterial_Type;
    typedef boost::shared_ptr<passiveMaterial_Type>       passiveMaterialPtr_Type;

    typedef EMActiveMaterialType<MeshType>                activeMaterial_Type;
    typedef boost::shared_ptr<activeMaterial_Type>        activeMaterialPtr_Type;

    //@}



    //! @name Constructor &  Deconstructor
    //@{

    //! Empty constructor
    EMStructuralConstitutiveLaw();
    //! Destructor
    virtual ~EMStructuralConstitutiveLaw() {}

    //@}

    //getters


    //SOME METHODS FOR ACTIVATED MATERIALS
    inline void setFiberVector ( const vector_Type& fiberVector)
    {
        M_fiberVectorPtr.reset ( new vector_Type ( fiberVector ) );
        ElectrophysiologyUtility::normalize (*M_fiberVectorPtr);
    }

    inline void setSheetVector ( const vector_Type& sheetVector)
    {
        M_sheetVectorPtr.reset ( new vector_Type ( sheetVector ) );
        ElectrophysiologyUtility::normalize (*M_sheetVectorPtr);
    }

    inline void setFiberVectorPtr ( const vectorPtr_Type fiberVectorPtr)
    {
        M_fiberVectorPtr = fiberVectorPtr;
    }

    inline void setSheetVectorPtr ( const vectorPtr_Type sheetVectorPtr)
    {
        M_sheetVectorPtr = sheetVectorPtr;
    }



    inline  vectorPtr_Type const fiberVectorPtr() const
    {
        return M_fiberVectorPtr;
    }

    inline  vectorPtr_Type fiberVectorPtr()
    {
        return M_fiberVectorPtr;
    }

    inline  vectorPtr_Type const sheetVectorPtr() const
    {
        return M_sheetVectorPtr;
    }

    inline  vectorPtr_Type sheetVectorPtr()
    {
        return M_sheetVectorPtr;
    }

    inline  scalarETFESpacePtr_Type scalarETFESpacePtr()
    {
        return M_scalarETFESpacePtr;
    }


    inline void setupFiberVector ( std::string& name, boost::shared_ptr<mesh_Type> mesh )
    {
        ElectrophysiologyUtility::importFibers ( M_fiberVectorPtr, name, mesh  );
        ElectrophysiologyUtility::normalize (*M_fiberVectorPtr);
    }

    inline void setupFiberVector ( std::string& name, std::string& path )
    {
        /* These two lines are giving compilation errors that I have not been able to solve.
           The function importFibers is not defined in the ElectroPhysiologyUtility or the
           Heartutility. Paolo Tricerri, June, 10th, 2014
        */
       //ElectrophysiologyUtility::importFibers ( M_fiberVectorPtr, name, path  );
        //ElectrophysiologyUtility::normalize (*M_fiberVectorPtr);
    }

    void setupFiberVector ( Real fx, Real fy, Real fz )
    {
        ElectrophysiologyUtility::setupFibers ( *M_fiberVectorPtr, fx, fy, fz  );
        ElectrophysiologyUtility::normalize (*M_fiberVectorPtr);
    }

    inline void setupSheetVector ( std::string& name, boost::shared_ptr<mesh_Type> mesh )
    {
        ElectrophysiologyUtility::importFibers ( M_sheetVectorPtr, name, mesh  );
        ElectrophysiologyUtility::normalize (*M_sheetVectorPtr);
    }

    inline void setupSheetVector ( std::string& name, std::string& path )
    {
        /* These two lines are giving compilation errors that I have not been able to solve.
           The function importFibers is not defined in the ElectroPhysiologyUtility or the
           Heartutility. Paolo Tricerri, June, 10th, 2014
        */
        //ElectrophysiologyUtility::importFibers ( M_sheetVectorPtr, name, path  );
        //ElectrophysiologyUtility::normalize (*M_sheetVectorPtr);
    }

    void setupSheetVector ( Real sx, Real sy, Real sz )
    {
        ElectrophysiologyUtility::setupFibers ( *M_sheetVectorPtr, sx, sy, sz);
        ElectrophysiologyUtility::normalize (*M_sheetVectorPtr);
    }

    inline virtual vectorPtr_Type& fiberActivationPtr()
    {
    	return M_fiberActivationPtr;
    }
    inline virtual vectorPtr_Type& sheetActivationPtr()
    {
    	return M_sheetActivationPtr;
    }
    inline virtual vectorPtr_Type& normalActivationPtr()
    {
    	return M_normalActivationPtr;
    }


    //! Setup the created object of the class StructuralConstitutiveLaw
    /*!
      \param dFespace: the FiniteElement Space
      \param monolithicMap: the MapEpetra
      \param offset: the offset parameter used assembling the matrices
    */
    virtual void setup ( const FESpacePtr_Type& dFESpace,
                         const ETFESpacePtr_Type& ETFESpace,
                         const boost::shared_ptr<const MapEpetra>&   monolithicMap,
                         const UInt offset, const dataPtr_Type& dataMaterial,
                         const displayerPtr_Type& displayer  );


    //! Computes the linear part of the stiffness matrix StructuralSolver::buildSystem
    /*!
      \param dataMaterial the class with Material properties data
    */
    virtual  void computeLinearStiff ( dataPtr_Type& dataMaterial,
                                       const mapMarkerVolumesPtr_Type /*mapsMarkerVolumes*/,
                                       const mapMarkerIndexesPtr_Type /*mapsMarkerIndexes*/ ) {}

    //! Updates the Jacobian matrix in StructuralSolver::updateJacobian
    /*!
      \param disp: solution at the k-th iteration of NonLinearRichardson Method
      \param dataMaterial: a pointer to the dataType member in StructuralSolver class to get the
                           material coefficients (e.g. Young modulus, Poisson ratio..)
      \param displayer: a pointer to the Dysplaier member in the StructuralSolver class
    */
    virtual  void updateJacobianMatrix ( const vector_Type& disp, const dataPtr_Type& dataMaterial,
                                                const mapMarkerVolumesPtr_Type mapsMarkerVolumes,
                                                const mapMarkerIndexesPtr_Type mapsMarkerIndexes,
                                                const displayerPtr_Type& displayer );

    //! Computes the new Stiffness matrix in StructuralSolver given a certain displacement field.
    //! This function is used both in StructuralSolver::evalResidual and in
    //! StructuralSolver::updateSystem since the matrix is the expression of the matrix is the same.
    //!This is virtual and not pure virtual since in the linear St. Venant-Kirchhoff law it is not needed.
    /*!
      \param sol:  the solution vector
      \param factor: scaling factor used in FSI
      \param dataMaterial: a pointer to the dataType member in StructuralSolver class to get the
                           material coefficients (e.g. Young modulus, Poisson ratio..)
      \param displayer: a pointer to the Dysplaier member in the StructuralSolver class
    */
    inline virtual  void computeStiffness ( const vector_Type& disp, Real factor, const dataPtr_Type& dataMaterial,
                                            const mapMarkerVolumesPtr_Type mapsMarkerVolumes,
                                            const mapMarkerIndexesPtr_Type mapsMarkerIndexes,
                                            const displayerPtr_Type& displayer );

    //! Computes the deformation Gradient F, the cofactor of F Cof(F),
    //! the determinant of F J = det(F), the trace of C Tr(C).
    /*!
      \param dk_loc: local displacement vector
    */
    inline virtual  void computeKinematicsVariables ( const VectorElemental& dk_loc ) {}


    //! Output of the class
    /*!
       \param fileNamelinearStiff the filename where to apply the spy method for the linear part of the Stiffness matrix
       \param fileNameStiff the filename where to apply the spy method for the Stiffness matrix
    */
    virtual void showMe ( std::string const& fileNameStiff, std::string const& fileNameJacobian ) {}


    //! Compute the First Piola Kirchhoff Tensor
    /*!
       \param firstPiola Epetra_SerialDenseMatrix that has to be filled
       \param tensorF Epetra_SerialDenseMatrix the deformation gradient
       \param cofactorF Epetra_SerialDenseMatrix cofactor of F
       \param invariants std::vector with the invariants of C and the detF
       \param material UInt number to get the material parameteres form the VenantElasticData class
    */
    inline virtual void computeLocalFirstPiolaKirchhoffTensor ( Epetra_SerialDenseMatrix& firstPiola,
                                                                const Epetra_SerialDenseMatrix& tensorF,
                                                                const Epetra_SerialDenseMatrix& cofactorF,
                                                                const std::vector<Real>& invariants,
                                                                const UInt material)
    {
        // Volumetric part
        auto Pvol = tensorF;
        Pvol.Scale( invariants[3] * (3500000 / 2.0) * (invariants[3] - 1 + (1.0 / invariants[3]) * std::log (invariants[3]) ) );

        // Isotropic part
        auto Piso = tensorF;
        Piso.Scale( 3300 * std::exp( 9.242*(invariants[0] - 3) ) );

        // Anisotropic fiber part
        auto Pi4 = tensorF;
        if ( invariants[1] > 1. ) Pi4.Scale( 2 * invariants[1] * 185350 * (invariants[1] - 1) * std::exp( 15.972*std::pow(invariants[1] - 1, 2.0 ) ) );
        else Pi4.Scale( 0. );
        
        // Active stress part
        auto Pact = tensorF;
        Pact.Scale( invariants[1] * 0.5 * std::pow(invariants[2], 2.0) * 1300000 );

        // Assemble first piola kirchhoff tensor
        firstPiola.Scale(0.0);
        if ( std::fabs(invariants[4]) < 0.1 || std::fabs(invariants[4] - 1.0) < 0.1 )
        {
            firstPiola += Pvol;
            firstPiola += Piso;
            firstPiola += Pi4;
        }
        if ( std::fabs(invariants[4]) < 0.1 || std::fabs(invariants[4] - 2.0) < 0.1 )
        {
            firstPiola += Pact;
        }
    }
    
    //! Get the Stiffness matrix
    matrixPtr_Type const stiffMatrix() const
    {
        return super::M_jacobian;
    }

    //! Get the stiffness vector
    vectorPtr_Type const stiffVector() const
    {
        return M_residualVectorPtr;
    }

    virtual void apply ( const vector_Type& sol, vector_Type& res,
                         const mapMarkerVolumesPtr_Type mapsMarkerVolumes,
                         const mapMarkerIndexesPtr_Type mapsMarkerIndexes) {};



//    inline vectorPtr_Type activationPtr()
//    {
//        return M_fiberActivationPtr;
//    }

    inline void setFiberActivationPtr (vectorPtr_Type activationPtr)
    {
        M_fiberActivationPtr = activationPtr;
    }
    inline void setFiberActivation (vector_Type& activation)
    {
        *M_fiberActivationPtr = activation;
    }
    inline void setSheetActivationPtr (vectorPtr_Type activationPtr)
    {
        M_sheetActivationPtr = activationPtr;
    }
    inline void setSheetActivation (vector_Type& activation)
    {
        *M_sheetActivationPtr = activation;
    }
    inline void setNormalActivationPtr (vectorPtr_Type activationPtr)
    {
        M_normalActivationPtr = activationPtr;
    }
    inline void setNormalActivation (vector_Type& activation)
    {
        *M_normalActivationPtr = activation;
    }



    void setParameters(EMData& data);


    void showMaterialParameters()
    {
    	if(M_passiveMaterialPtr) M_passiveMaterialPtr-> showMe();
    	if(M_activeStressMaterialPtr) M_activeStressMaterialPtr-> showMe();
    }

    //@}



protected:
    virtual void setupVectorsParameters ( void ) {}
    //ET finite element space for scalar variables
    scalarETFESpacePtr_Type                        M_scalarETFESpacePtr;

    vectorPtr_Type                                 M_fiberVectorPtr;

    vectorPtr_Type                                 M_sheetVectorPtr;

    vectorPtr_Type                                 M_residualVectorPtr;

    //    MatrixSmall<3,3>                               M_identity;

    //    std::vector<materialPtr_Type>                M_materialPtrs;
    //materialPtr_Type                               M_materialPtr;
    //    Real                                         M_bulkModulus;

    passiveMaterialPtr_Type                               M_passiveMaterialPtr;
    activeMaterialPtr_Type                               M_activeStressMaterialPtr;

    vectorPtr_Type                                 M_fiberActivationPtr;
    vectorPtr_Type                                 M_sheetActivationPtr;
    vectorPtr_Type                                 M_normalActivationPtr;



};

template<typename MeshType>
EMStructuralConstitutiveLaw<MeshType>::EMStructuralConstitutiveLaw() :
    super                           ( ),
    M_scalarETFESpacePtr            ( ),
    M_fiberVectorPtr                ( ),
    M_sheetVectorPtr                ( ),
    M_passiveMaterialPtr            ( ),
    M_activeStressMaterialPtr       ( ),
    M_residualVectorPtr             ( ),
    M_fiberActivationPtr                 ( )
{}

template <typename MeshType>
void
EMStructuralConstitutiveLaw<MeshType>::setup ( const FESpacePtr_Type&                      dFESpace,
                                               const ETFESpacePtr_Type&                    dETFESpace,
                                               const boost::shared_ptr<const MapEpetra>&   monolithicMap,
                                               const UInt                                  offset,
                                               const dataPtr_Type&                         dataMaterial,
                                               const displayerPtr_Type&                    displayer)
{
    this->M_displayer = displayer;
    this->M_dataMaterial  = dataMaterial;

    //    std::cout<<"I am setting up the Material"<<std::endl;

    this->M_dispFESpace                 = dFESpace;
    this->M_dispETFESpace               = dETFESpace;
    this->M_localMap                    = monolithicMap;
    this->M_offset                      = offset;
    this->M_dataMaterial                = dataMaterial;
    this->M_displayer                   = displayer;

    M_residualVectorPtr.reset ( new vector_Type (*this->M_localMap, Repeated) );
    //   M_identity = EMUtility::identity();

//    M_fiberVectorPtr.reset             ( new vector_Type (*this->M_localMap, Repeated) );
//    M_sheetVectorPtr.reset             ( new vector_Type (*this->M_localMap, Repeated) );
    M_fiberVectorPtr.reset             ( new vector_Type (*this->M_localMap, Unique) );
    M_sheetVectorPtr.reset             ( new vector_Type (*this->M_localMap, Unique) );
    M_scalarETFESpacePtr.reset         ( new scalarETFESpace_Type ( dETFESpace -> mesh(),
    																&( dETFESpace -> refFE() ),
                                                                    dFESpace->map().commPtr() ) );

    M_fiberActivationPtr.reset (new vector_Type (M_scalarETFESpacePtr -> map() ) );

    std::string passiveMaterialType ( dataMaterial -> passiveType() );
    std::string activeStressMaterialType (dataMaterial -> activeStressType() );
    std::cout << "\n===========================";
    std::cout << "\nPassive Type: " << passiveMaterialType;
    std::cout << "\nActive Stress Type: " << activeStressMaterialType;
    std::cout << "\n===========================\n";

    if (activeStressMaterialType != "NO_DEFAULT_ACTIVESTRESS_TYPE")
    {
        M_activeStressMaterialPtr.reset (activeMaterial_Type::EMActiveMaterialFactory::instance().createObject ( activeStressMaterialType ) );
        if(dFESpace->map().commPtr() ->MyPID() == 0)
        {
			std::cout << "\nCreated Active Stress Material!\n";
			M_activeStressMaterialPtr-> showMe();
			std::cout << "\nShowed Active Stress Material!\n";
        }

    }
    if (passiveMaterialType != "NO_DEFAULT_PASSIVE_TYPE")
    {
        //    M_materialPtr.reset(new EMMaterial<MeshType>(materialType));
        M_passiveMaterialPtr.reset (passiveMaterial_Type::EMPassiveMaterialFactory::instance().createObject ( passiveMaterialType ) );
        if(dFESpace->map().commPtr() ->MyPID() == 0)
        {
			std::cout << "\nCreated Passive Material!\n";
			M_passiveMaterialPtr -> showMe();
			std::cout << "\nCreated Passive Material!\n";
        }
    }

    //    // The 2 is because the law uses three parameters (mu, bulk).
    //    // another way would be to set up the number of constitutive parameters of the law
    //    // in the data file to get the right size. Note the comment below.
    //    this->M_vectorsParameters.reset ( new vectorsParameters_Type ( 2 ) );
    //
    //    this->setupVectorsParameters();
}


template <typename MeshType>
void EMStructuralConstitutiveLaw<MeshType>::updateJacobianMatrix ( const vector_Type&       disp,
                                                                   const dataPtr_Type&      dataMaterial,
                                                                   const mapMarkerVolumesPtr_Type mapsMarkerVolumes,
                                                                   const mapMarkerIndexesPtr_Type mapsMarkerIndexes,
                                                                   const displayerPtr_Type& displayer )
{
    this->M_jacobian.reset (new matrix_Type (*this->M_localMap) );
    //    matrixPtr_Type jac(new matrix_Type(*this->M_localMap));

    //displayer->leaderPrint (" \n*********************************\n  ");
    //displayer->leaderPrint (" Non-Linear S-  Computing the EM material  Jacobian"     );
    //displayer->leaderPrint (" \n*********************************\n  ");
    * (this->M_jacobian) *= 0.0;
    
    
    
    
    
    
    
    
    
    if (M_passiveMaterialPtr)
    {

        M_passiveMaterialPtr -> computeJacobian (disp,
                                                 this->M_dispETFESpace,
                                                 *M_fiberVectorPtr,
                                                 *M_sheetVectorPtr,
                                                 this->M_jacobian);
    }
    if (M_activeStressMaterialPtr)
        M_activeStressMaterialPtr -> computeJacobian ( disp,
                                                       this->M_dispETFESpace,
                                                       *M_fiberVectorPtr,
                                                       *M_sheetVectorPtr,
                                                       M_fiberActivationPtr,
                                                       M_sheetActivationPtr,
                                                       M_normalActivationPtr,
                                                       M_scalarETFESpacePtr,
                                                       this->M_jacobian);


    //  computeJacobian(disp);

    this->M_jacobian->globalAssemble();
    //displayer->leaderPrint (" \n*********************************\n\n  ");
    //std::cout << std::endl;
}

template <typename MeshType>
void EMStructuralConstitutiveLaw<MeshType>::computeStiffness ( const vector_Type&       disp,
                                                               Real                     /*factor*/,
                                                               const dataPtr_Type&      dataMaterial,
                                                               const mapMarkerVolumesPtr_Type mapsMarkerVolumes,
                                                               const mapMarkerIndexesPtr_Type mapsMarkerIndexes,
                                                               const displayerPtr_Type& displayer )
{
    //displayer->leaderPrint (" \n******************************************************************\n  ");
    //displayer->leaderPrint (" Non-Linear S-  Computing the EM material residual vector"     );
    //displayer->leaderPrint (" \n******************************************************************\n  ");
    * (M_residualVectorPtr) *= 0.0;
    vectorPtr_Type vec (new vector_Type ( M_residualVectorPtr -> map() ) );
    Real passive_residual(0.0);
    
    
    class HeavisideFct
    {
    public:
        typedef Real return_Type;
        return_Type operator() (const Real& I4f)
        {
            return (I4f > 0. ? 1. : 0.);
        }
        
        HeavisideFct() {}
        HeavisideFct (const HeavisideFct&) {}
        ~HeavisideFct() {}
    };

    boost::shared_ptr<HeavisideFct> heaviside (new HeavisideFct);
    
    if (M_passiveMaterialPtr)
    {
        
        using namespace ExpressionAssembly;
        
        auto F = _F (super::M_dispETFESpace, disp, 0);
        auto I = _I;
        auto J = det(F);
        auto Jm23 = pow(J, 2 / (-3.) );
        auto FmT = minusT(F);
        auto I1 = dot(F, F);
        auto dI1bar = value(2.0) * Jm23 * (  F + value(1/(-3.)) * I1 * FmT );

        
        
        
        auto f_0 = _v0 (super::M_dispETFESpace, *M_fiberVectorPtr);
        auto s_0 = _v0 (super::M_dispETFESpace, *M_sheetVectorPtr);
        
        boost::shared_ptr<orthonormalizeFibers> normalize0 (new orthonormalizeFibers);
        auto f0 = eval (normalize0, f_0);
        auto f = F * f0;
        
        boost::shared_ptr<orthonormalizeFibers> normalize1 (new orthonormalizeFibers (1) );
        auto s0 = eval (normalize1, f0, s_0);
        auto s = F * s0;

        boost::shared_ptr<CrossProduct> wedge (new CrossProduct);
        auto n0 = eval ( wedge, f0, s0);
        
        
        
        // Volumetric
        auto dWvol = ( 3500000 * ( J + J * log(J) - 1. ) ) / ( 2 * J );
        auto Pvol = dWvol * _dJ (F);
        integrate ( elements ( super::M_dispETFESpace->mesh() ) ,
                   quadRuleTetra4pt,
                   super::M_dispETFESpace,
                   dot ( Pvol , grad (phi_i) )
                   ) >> M_residualVectorPtr;
        
        
        // I1
        auto I1bar = _I1bar (F);
        auto dW1 = 0.5 * 3300 * exp ( 9.242 * ( I1bar - 3 ) );
        auto Pw1 = dW1 * dI1bar ;
        //	auto P = eval (W1, F ) * _dI1bar(F) ;
        integrate ( elements ( super::M_dispETFESpace->mesh() ) ,
                   quadRuleTetra4pt,
                   super::M_dispETFESpace,
                   dot ( Pw1, grad (phi_i) )
                   ) >> M_residualVectorPtr;


        // I4f
        auto I4f = dot (f,f);
        auto I4m1f = I4f - 1.0;
        auto dW4f = 185350 * I4m1f * exp (15.972 * I4m1f * I4m1f ) * eval(heaviside, I4m1f);
        auto Pw4f = dW4f * _dI4 (F, f0);
        //
        //    auto P = eval (W4, _I4bar ( dispETFESpace, disp, 0, f0 ) )
        //             * _dI4bar (dispETFESpace, disp, 0, f0);
        integrate ( elements ( super::M_dispETFESpace->mesh() ) ,
                   quadRuleTetra4pt,
                   super::M_dispETFESpace,
                   dot ( Pw4f , grad (phi_i) )
                   ) >> M_residualVectorPtr;
        
        
        // I4s
        auto I4s = dot (s,s);
        auto I4m1s = I4s - 1.0;
        auto dW4s = 25640 * I4m1s * exp (10.446 * I4m1s * I4m1s ) *  eval(heaviside, I4m1s);
        auto Pw4s = dW4s * _dI4 ( F, s0);
        integrate ( elements ( super::M_dispETFESpace->mesh() ) ,
                   quadRuleTetra4pt,
                   super::M_dispETFESpace,
                   dot ( Pw4s , grad (phi_i) )
                   ) >> M_residualVectorPtr;
        
        
        // I8fs
        auto I8fs = dot (f,s);
        auto dW8fs = 4170 * I8fs * exp ( 11.602 * I8fs * I8fs );
        auto Pw8fs = dW8fs * _dI8 ( F, f0, s0);
        integrate ( elements ( super::M_dispETFESpace->mesh() ),
                   quadRuleTetra4pt,
                   super::M_dispETFESpace,
                   dot (  Pw8fs, grad (phi_i) )
                   ) >> M_residualVectorPtr;


        
    // Active strain I1E
    auto gf = value (M_scalarETFESpacePtr, *M_fiberActivationPtr);

    auto k = 4.0;
    
    //auto gm = value(-1.0) * ( gf ) / ( ( gf ) + 1.0 );
    //auto go = gf * ( k + gf * k + value(1.0) );
    //auto gmn = value(-1.0) * ( k*gf ) / ( ( k*gf ) + 1.0 ) ;
        
    auto FAinv = _FAinvO( gf, k, f0, s0, n0 ); // I + gm * outerProduct(f0, f0) + go * outerProduct(s0, s0) + gmn * outerProduct(n0, n0);
    auto FE =  F * FAinv;
    
    auto I1barE = _I1bar (FE);
    auto dWI1E = 3300 / 2.0 * exp ( 9.242 * ( I1barE - 3 ) );

    auto PI1barE = dWI1E * _dI1bar (FE) * FAinv;
    
    integrate ( elements ( super::M_dispETFESpace->mesh() ) ,
               quadRuleTetra4pt,
               super::M_dispETFESpace,
               dot ( PI1barE, grad (phi_i) )
               ) >> M_residualVectorPtr;

    
    // Active strain I4fE
    auto I4fbarE =  _I4bar (FE, f0);
    auto I4fbarEm1 = I4fbarE - 1.0;
    auto dW4fE = 185350 * I4fbarEm1 * exp (15.972 * I4fbarEm1 * I4fbarEm1 ) * eval(heaviside, I4fbarEm1);
    auto Pw4fE = dW4fE * _dI4bar (FE, f0) * FAinv;
    integrate ( elements ( super::M_dispETFESpace->mesh() ) ,
               quadRuleTetra4pt,
               super::M_dispETFESpace,
               dot ( Pw4fE, grad (phi_i) )
               ) >> M_residualVectorPtr;
    
    
    // Active strain I4sE
    auto I4sbarE =  _I4bar (FE, s0);
    auto I4sbarEm1 = I4sbarE - 1.0;
    auto dW4sE = 25640 * I4sbarEm1 * exp (10.446 * I4sbarEm1 * I4sbarEm1 ) * eval(heaviside, I4sbarEm1);
    auto Pw4sE = dW4sE * _dI4bar (FE, s0) * FAinv;
    integrate ( elements ( super::M_dispETFESpace->mesh() ) ,
               quadRuleTetra4pt,
               super::M_dispETFESpace,
               dot ( Pw4sE, grad (phi_i) )
               ) >> M_residualVectorPtr;
    
    
    // Active strain I8fsE
    auto I8barfsE = _I8bar(FE, f0, s0);
    auto dW8fsE = 4170 * I8barfsE * exp ( 11.602 * I8barfsE * I8barfsE );
        
    auto Pw8fsE = dW8fsE * _dI8bar (FE, f0, s0) * FAinv;
    
    integrate ( elements ( super::M_dispETFESpace->mesh() ) ,
               quadRuleTetra4pt,
               super::M_dispETFESpace,
               dot ( Pw8fsE, grad (phi_i) )
               ) >> M_residualVectorPtr;
    
 
    }
    
//    using namespace ExpressionAssembly;
//    
//    ///this->M_stiff.reset (new vector_Type (*this->M_localMap) );
//    
//    displayer->leaderPrint (" \n******************************************************************\n  ");
//    displayer->leaderPrint (" Non-Linear S-  Computing the Holzapfel-Ogden nonlinear stiffness vector"     );
//    displayer->leaderPrint (" \n******************************************************************\n  ");
//    
//    ///M_stiff.reset (new vector_Type (*this->M_localMap) );
//    ///* (M_stiff) *= 0.0;
//    
//    // Activation is Fa = gamma (fo x fo) + 1/sqrt(gamma) * (I - fo x fo)
//    BOOST_AUTO_TPL (gamma,  value (1.0) + value (this->M_activationSpace, *M_Gammaf) );
//    
//    MatrixSmall<3, 3> Id;
//    Id (0, 0) = 1.;
//    Id (0, 1) = 0., Id (0, 2) = 0.;
//    Id (1, 0) = 0.;
//    Id (1, 1) = 1., Id (1, 2) = 0.;
//    Id (2, 0) = 0.;
//    Id (2, 1) = 0., Id (2, 2) = 1.;
//    
//    // Strain-energy terms
//    boost::shared_ptr<StrainEnergyHO::dW1> dW1fun (new StrainEnergyHO::dW1 (M_aiso, M_biso) );
//    boost::shared_ptr<StrainEnergyHO::d2W1> d2W1fun (new StrainEnergyHO::d2W1 (M_aiso, M_biso) );
//    
//    boost::shared_ptr<StrainEnergyHO::dW4> dW4ffun (new StrainEnergyHO::dW4 (M_af, M_bf) );
//    boost::shared_ptr<StrainEnergyHO::d2W4> d2W4ffun (new StrainEnergyHO::d2W4 (M_af, M_bf) );
//    
//    boost::shared_ptr<StrainEnergyHO::dW4> dW4sfun (new StrainEnergyHO::dW4 (M_as, M_bs) );
//    boost::shared_ptr<StrainEnergyHO::d2W4> d2W4sfun (new StrainEnergyHO::d2W4 (M_as, M_bs) );
//    
//    boost::shared_ptr<StrainEnergyHO::dW8> dW8fsfun (new StrainEnergyHO::dW8 (M_afs, M_bfs) );
//    boost::shared_ptr<StrainEnergyHO::d2W8> d2W8fsfun (new StrainEnergyHO::d2W8 (M_afs, M_bfs) );
//    
//    boost::shared_ptr<StrainEnergyHO::dWvol> dWvolfun (new StrainEnergyHO::dWvol (M_kappa) );
//    boost::shared_ptr<StrainEnergyHO::d2Wvol> d2Wvolfun (new StrainEnergyHO::d2Wvol (M_kappa) );
//    
//    // Kinematics
//    BOOST_AUTO_TPL (I,      value (Id) );
//    BOOST_AUTO_TPL (Grad_u, grad (this->M_dispETFESpace, disp, this->M_offset) );
//    BOOST_AUTO_TPL (F,      Grad_u + I);
//    BOOST_AUTO_TPL (FmT,    minusT (F) );
//    BOOST_AUTO_TPL (J,      det (F) );
//    // Fibres
//    BOOST_AUTO_TPL (f0,     value (this->M_dispETFESpace, *M_fiberVector) );
//    BOOST_AUTO_TPL (s0,     value (this->M_dispETFESpace, *M_sheetVector) );
//    // Invariants
//    BOOST_AUTO_TPL (I1,     dot (F, F) );
//    BOOST_AUTO_TPL (I4f,    dot (F * f0, F * f0) );
//    BOOST_AUTO_TPL (I4s,    dot (F * s0, F * s0) );
//    BOOST_AUTO_TPL (I8fs,   dot (F * f0, F * s0) );
//    // Reduced invariants
//    BOOST_AUTO_TPL (Jm23,    pow (J, -2. / 3) );
//    BOOST_AUTO_TPL (I1iso,   Jm23 * I1);
//    BOOST_AUTO_TPL (I4fiso,  Jm23 * I4f);
//    BOOST_AUTO_TPL (I4siso,  Jm23 * I4s);
//    BOOST_AUTO_TPL (I8fsiso, Jm23 * I8fs);
//    // Generalised invariants
//    BOOST_AUTO_TPL (gammac, pow (gamma, -2) - gamma);
//    
//    BOOST_AUTO_TPL (I1eiso,   gamma * I1iso + gammac * I4fiso);
//    BOOST_AUTO_TPL (I4feiso,  pow (gamma, -2) * I4fiso);
//    BOOST_AUTO_TPL (I4seiso,  gamma * I4siso);
//    BOOST_AUTO_TPL (I8fseiso, pow (gamma, -1. / 2) * I8fsiso);
//    
//    // Strain-energy derivatives
//    BOOST_AUTO_TPL (dW1,    eval (dW1fun, I1eiso) );
//    BOOST_AUTO_TPL (d2W1,   eval (d2W1fun, I1eiso) );
//    BOOST_AUTO_TPL (dW4f,   eval (dW4ffun, I4feiso) );
//    BOOST_AUTO_TPL (d2W4f,  eval (d2W4ffun, I4feiso) );
//    BOOST_AUTO_TPL (dW4s,   eval (dW4sfun, I4seiso) );
//    BOOST_AUTO_TPL (d2W4s,  eval (d2W4sfun, I4seiso) );
//    BOOST_AUTO_TPL (dW8fs,  eval (dW8fsfun, I8fseiso) );
//    BOOST_AUTO_TPL (d2W8fs, eval (d2W8fsfun, I8fseiso) );
//    
//    BOOST_AUTO_TPL (dWvol,  eval (dWvolfun, J) );
//    BOOST_AUTO_TPL (d2Wvol, eval (d2Wvolfun, J) );
//    
//    // Deviatorics
//    BOOST_AUTO_TPL (Fdev1,   F - value (1. / 3) *I1 * FmT);
//    BOOST_AUTO_TPL (Fdev4f,  outerProduct (F * f0, f0) - value (1. / 3) *I4f * FmT);
//    BOOST_AUTO_TPL (Fdev4s,  outerProduct (F * s0, s0) - value (1. / 3) *I4s * FmT);
//    BOOST_AUTO_TPL (Fdev8fs, value (0.5) * (outerProduct (F * s0, f0) + outerProduct (F * f0, s0) ) - value (1. / 3) *I8fs * FmT);
//    
//    // Derivatives
//    BOOST_AUTO_TPL (dF,   grad (phi_j) );
//    BOOST_AUTO_TPL (dJ,   J * FmT);
//    BOOST_AUTO_TPL (dFmT, value (-1.0) *FmT * transpose (dF) *FmT);
//    
//    BOOST_AUTO_TPL (dI1eiso,   value (2.) *Jm23 * (gamma * Fdev1 + gammac * Fdev4f) );
//    BOOST_AUTO_TPL (dI4feiso,  value (2.) *Jm23 * pow (gamma, -2) *Fdev4f);
//    BOOST_AUTO_TPL (dI4seiso,  value (2.) *Jm23 * gamma * Fdev4s);
//    BOOST_AUTO_TPL (dI8fseiso, value (2.) *Jm23 * pow (gamma, -1. / 2) *Fdev8fs);
//    
//    BOOST_AUTO_TPL (dFdev1,   dF + value (-1. / 3) * (value (2.0) *dot (F, dF) *FmT + I1 * dFmT) );
//    BOOST_AUTO_TPL (dFdev4f,  outerProduct (dF * f0, f0) - value (1. / 3) * (value (2.0) *dot (F * f0, dF * f0) *FmT + I4f * dFmT) );
//    BOOST_AUTO_TPL (dFdev4s,  outerProduct (dF * s0, s0) - value (1. / 3) * (value (2.0) *dot (F * s0, dF * s0) *FmT + I4s * dFmT) );
//    BOOST_AUTO_TPL (dFdev8fs,   value (0.5) * (outerProduct (dF * s0, f0) + outerProduct (dF * f0, s0) )
//                    - value (1. / 3) * ( (dot (F * s0, dF * f0) + dot (F * f0, dF * s0) ) *FmT + I8fs * dFmT) );
//    
//    BOOST_AUTO_TPL (d2I1eiso,    value (2.) *Jm23 * (gamma * dFdev1 + gammac * dFdev4f
//                                                     - value (2. / 3) *dot (FmT, dF) * (gamma * Fdev1 + gammac * Fdev4f) ) );
//    BOOST_AUTO_TPL (d2I4feiso,   value (2.) *Jm23 * pow (gamma, -2) * (dFdev4f
//                                                                       - value (2. / 3) *dot (FmT, dF) *Fdev4f) );
//    BOOST_AUTO_TPL (d2I4seiso,   value (2.) *Jm23 * gamma * (dFdev4s
//                                                             - value (2. / 3) *dot (FmT, dF) *Fdev4s) );
//    BOOST_AUTO_TPL (d2I8fseiso,   value (2.) *Jm23 * pow (gamma, -1. / 2) * (dFdev8fs
//                                                                             - value (2. / 3) *dot (FmT, dF) *Fdev8fs) );
//    
//    // Isochoric part
//    BOOST_AUTO_TPL (Piso_1,   dW1 * dI1eiso);
//    BOOST_AUTO_TPL (Piso_4f,  dW4f * dI4feiso);
//    BOOST_AUTO_TPL (Piso_4s,  dW4s * dI4seiso);
//    BOOST_AUTO_TPL (Piso_8fs, dW8fs * dI8fseiso);
//    
//    BOOST_AUTO_TPL (dPiso_1,   d2W1 * dot (dI1eiso, dF) *dI1eiso + dW1 * d2I1eiso );
//    BOOST_AUTO_TPL (dPiso_4f,  d2W4f * dot (dI4feiso, dF) *dI4feiso + dW4f * d2I4feiso );
//    BOOST_AUTO_TPL (dPiso_4s,  d2W4s * dot (dI4seiso, dF) *dI4seiso + dW4s * d2I4seiso );
//    BOOST_AUTO_TPL (dPiso_8fs, d2W8fs * dot (dI8fseiso, dF) *dI8fseiso + dW8fs * d2I8fseiso );
//    
//    // Volumetric part
//    BOOST_AUTO_TPL (Pvol,  dWvol * dJ);
//    BOOST_AUTO_TPL (dPvol, (d2Wvol + J * dWvol) *dot (FmT, dF) *dJ + dWvol * J * dFmT);
//    
//    // Assemble residual
//    integrate ( elements ( this->M_dispETFESpace->mesh() ),
//               this->M_dispFESpace->qr(),
//               this->M_dispETFESpace,
//               dot (Piso_1 + Pvol, grad (phi_i) )
//               ) >> M_stiff;
//    // Assemble residual
//    if (M_af > 0)
//    {
//        integrate ( elements ( this->M_dispETFESpace->mesh() ),
//                   this->M_dispFESpace->qr(),
//                   this->M_dispETFESpace,
//                   dot ( Piso_4f , grad (phi_i) )
//                   ) >> M_stiff;
//    }
//    // Assemble residual
//    if (M_as > 0)
//    {
//        
//        integrate ( elements ( this->M_dispETFESpace->mesh() ),
//                   this->M_dispFESpace->qr(),
//                   this->M_dispETFESpace,
//                   dot ( Piso_4s , grad (phi_i) )
//                   ) >> M_stiff;
//    }
//    // Assemble residual
//    if (M_afs > 0)
//    {
//        integrate ( elements ( this->M_dispETFESpace->mesh() ),
//                   this->M_dispFESpace->qr(),
//                   this->M_dispETFESpace,
//                   dot ( Piso_8fs, grad (phi_i) )
//                   ) >> M_stiff;
//    }
//    this->M_stiff->globalAssemble();
//    
//    
//    
//    
//    
//    
//    
//    
//    if (M_passiveMaterialPtr)
//    {
//        M_passiveMaterialPtr -> computeResidual (disp,
//                                                 this->M_dispETFESpace,
//                                                 *M_fiberVectorPtr,
//                                                 *M_sheetVectorPtr,
//                                                 M_residualVectorPtr);
//        passive_residual = M_residualVectorPtr -> norm2();
//        if(vec->map().commPtr() ->MyPID() == 0)
//        {
//        	//std::cout << "\nPassive Residual: " << passive_residual << "\n";
//        }
//    }

//    if (M_activeStressMaterialPtr)
//	{
//    	M_activeStressMaterialPtr -> computeResidual ( disp,
//                                                      super::M_dispETFESpace,
//                                                       *M_fiberVectorPtr,
//                                                       *M_sheetVectorPtr,
//                                                       M_fiberActivationPtr,
//                                                       M_sheetActivationPtr,
//                                                       M_normalActivationPtr,
//                                                       M_scalarETFESpacePtr,
//                                                       M_residualVectorPtr);
//
//        if(vec->map().commPtr() ->MyPID() == 0)
//        {
//            //std::cout << "\nActive Residual: " << M_residualVectorPtr -> norm2() - passive_residual << "\n";
//        }
//	}
    //  computeResidual(disp);
    this->M_residualVectorPtr->globalAssemble();
}

template <typename MeshType>
void EMStructuralConstitutiveLaw<MeshType>::setParameters(EMData& data)
{
    if (M_activeStressMaterialPtr)
    {
        M_activeStressMaterialPtr-> setParameters(data);
    }

    if (M_passiveMaterialPtr)
    {
        M_passiveMaterialPtr-> setParameters(data);
    }
}




template <typename MeshType>
inline StructuralConstitutiveLaw<MeshType>* createEMMaterial()
{
    return new EMStructuralConstitutiveLaw<MeshType>();
}
namespace
{
static bool registerEM = StructuralConstitutiveLaw<LifeV::RegionMesh<LinearTetra> >::StructureMaterialFactory::instance().registerProduct ( "EMMaterial", &createEMMaterial<LifeV::RegionMesh<LinearTetra> > );
}


}
#endif /*_STRUCTURALMATERIAL_H*/
