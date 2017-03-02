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

//#include <lifev/em/solver/mechanics/EMMechanicalExpressions.hpp>
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

    
    
//    typedef StructuralConstitutiveLaw<Mesh>                 super;
//    
//    typedef typename super::data_Type                data_Type;
//    
//    typedef typename super::vector_Type              vector_Type;
//    typedef typename super::matrix_Type              matrix_Type;
//    
//    typedef typename super::matrixPtr_Type           matrixPtr_Type;
//    typedef typename super::vectorPtr_Type           vectorPtr_Type;
//    typedef typename super::dataPtr_Type             dataPtr_Type;
//    typedef typename super::displayerPtr_Type        displayerPtr_Type;
    
//    typedef typename super::mapMarkerVolumesPtr_Type mapMarkerVolumesPtr_Type;
//    typedef typename super::mapMarkerVolumes_Type mapMarkerVolumes_Type;
    typedef typename mapMarkerVolumes_Type::const_iterator mapIterator_Type;
    
    
    
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
    
    void computeStiffness2 (    const vector_Type& sol,
                                Real /*factor*/,
                                const dataPtr_Type& dataMaterial,
                                const mapMarkerVolumesPtr_Type mapsMarkerVolumes,
                                const displayerPtr_Type& displayer );


    //! Computes the deformation Gradient F, the cofactor of F Cof(F),
    //! the determinant of F J = det(F), the trace of C Tr(C).
    /*!
      \param dk_loc: local displacement vector
    */
    inline virtual  void computeKinematicsVariables ( const VectorElemental& dk_loc ); // {}

    inline virtual  void computeAnisotropyVariables ( const VectorElemental& fk_loc, const VectorElemental& sk_loc, const VectorElemental& fAk_loc );
    
    
    
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
        
    }
    
    
    
    Epetra_SerialDenseVector matrixTimesVector( const Epetra_SerialDenseMatrix& A, const Epetra_SerialDenseVector& x ) const
    {
        Epetra_SerialDenseVector y (3);
        for (UInt i (0); i < 3; ++i)
        {
            for (UInt j (0); j < 3; ++j)
            {
                y(i) = A(i,j) * x(j);
            }
        }
        return y;
    }
    
    
    Epetra_SerialDenseMatrix tensorProduct( const Epetra_SerialDenseVector& v, const Epetra_SerialDenseVector& w ) const
    {
        Epetra_SerialDenseMatrix m (3,3);
        for (UInt i (0); i < 3; ++i)
        {
            for (UInt j (0); j < 3; ++j)
            {
                m(i,j) = v(i) * w(j);
            }
        }
        return m;
    }
    
    
    inline virtual void computeLocalFirstPiolaKirchhoffTensor_ ( Epetra_SerialDenseMatrix& firstPiola,
                                                               const Epetra_SerialDenseMatrix& tensorF,
                                                               const Epetra_SerialDenseMatrix& cofactorF,
                                                               const Epetra_SerialDenseVector& fiber,
                                                               const Epetra_SerialDenseVector& sheet,
                                                               const std::vector<Real>& invariants,
                                                               const UInt material)
    {
        auto I1 = invariants[0];
        auto J = invariants[1];
        auto I4f = invariants[2];
        auto I4s = invariants[3];
        auto I8fs = invariants[4];
        
        auto gammaf = invariants[5];
        auto gamman = 4.0 * gammaf;
        auto gammas = 1.0 / ( (1.0+gammaf)*(1.0+gamman) ) - 1.0;
        
        auto g1 = 1 - ( gamman * ( gamman + 2 ) / std::pow( gamman+1 , 2.0 ) );
        auto g4f = ( gamman * ( gamman + 2 ) / std::pow( gamman+1 , 2.0 ) ) - ( gammaf * ( gammaf + 2 ) / std::pow( gammaf+1 , 2.0 ) );
        auto g4s = ( gamman * ( gamman + 2 ) / std::pow( gamman+1 , 2.0 ) ) - ( gammas * ( gammas + 2 ) / std::pow( gammas+1 , 2.0 ) );
        
        auto I1E = g1 * I1 + g4f * I4f + g4s * I4s;
        auto JE = J;
        auto I1barE = std::pow(JE, -2.0/3.0 ) * I1E;
        auto I4fE = I4f / std::pow( gammaf + 1 , 2.0 );
        auto I4sE = I4s / std::pow( gammas + 1 , 2.0 );
        auto I8fsE = I8fs / ( (gammaf + 1) * (gammas + 1) );
        
        auto W1E = 0.5 * 3300 * std::exp( 9.242 * (I1barE - 3) );
        auto W4fE = 185350 * (I4fE - 1) * std::exp( 15.972 * std::pow(I4fE - 1, 2.0) ) * (I4fE > 1.0);
        auto W4sE = 25640 * (I4sE - 1) * std::exp (10.446 * std::pow(I4sE - 1, 2.0) ) * (I4sE > 1.0);
        auto W8fsE = 4170 * I8fsE * std::exp ( 11.602 * I8fsE * I8fsE );

        auto f = matrixTimesVector(tensorF, fiber);
        auto s = matrixTimesVector(tensorF, sheet);
        auto f_f0 = tensorProduct(f, fiber);
        auto s_s0 = tensorProduct(s, sheet);
        auto f_s0 = tensorProduct(f, sheet);
        auto s_f0 = tensorProduct(s, fiber);

        
        // Pvol
        Epetra_SerialDenseMatrix Pvol (3,3);
        Pvol.Scale(0.0);
        Pvol += cofactorF;
        Pvol.Scale( J * (3500000 / 2.0) * (J - 1.0 + (1.0 / J) * std::log(J) ) );

        // P1
        Epetra_SerialDenseMatrix P1 (3,3);
        P1.Scale(0.0);
        P1 += cofactorF;
        P1.Scale(-I1/3);
        P1 += tensorF;
        P1.Scale( 2.0 * g1 * W1E * std::pow(J, -2.0/3.0 ) );

        // P4f
        Epetra_SerialDenseMatrix P4f (3,3);
        P4f.Scale(0.0);
        P4f += f_f0;
        P4f.Scale ( 2.0 * ( g4f * W1E + W4fE / std::pow( gammaf + 1.0 , 2.0 ) ) );

        // P4s
        Epetra_SerialDenseMatrix P4s (3,3);
        P4s.Scale(0.0);
        P4s += s_s0;
        P4s.Scale ( 2 * ( g4s * W1E + W4sE / std::pow( gammas + 1.0 , 2.0 ) ) );

        // P8fs
        Epetra_SerialDenseMatrix P8fs (3,3);
        P8fs.Scale(0.0);
        P8fs += f_s0;
        P8fs += s_f0;
        P8fs.Scale( W8fsE / ( (gammaf + 1.0) * (gammas + 1.0) ) );

        // Assemble first piola kirchhoff tensor
        firstPiola.Scale(0.0);
        firstPiola += Pvol;
        firstPiola += P1;
        firstPiola += P4f;
        firstPiola += P4s;
        firstPiola += P8fs;
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

    
    
    
    

    
    

    void orthonormalize (const boost::shared_ptr<boost::multi_array<Real, 2> > & v, const boost::shared_ptr<boost::multi_array<Real, 2> > & v1, const UInt & ig) const
    {
        normalize(v, ig);
        
        Real dotProduct = (*v)[0][ig] * (*v1)[0][ig] + (*v)[1][ig] * (*v1)[1][ig] + (*v)[2][ig] * (*v1)[2][ig];
        (*v)[0][ig] = (*v)[0][ig] - dotProduct * (*v1)[0][ig];
        (*v)[1][ig] = (*v)[1][ig] - dotProduct * (*v1)[1][ig];
        (*v)[2][ig] = (*v)[2][ig] - dotProduct * (*v1)[2][ig];
        
        normalize(v, ig);
    }
        
        
        
    void normalize (const boost::shared_ptr<boost::multi_array<Real, 2> > & v, const UInt & ig) const
    {
        Real norm = std::sqrt ( (*v)[0][ig] * (*v)[0][ig] + (*v)[1][ig] * (*v)[1][ig] + (*v)[2][ig] * (*v)[2][ig]);
        if ( norm >= 1e-13 )
        {
            (*v)[0][ig] = (*v)[0][ig] / norm;
            (*v)[1][ig] = (*v)[1][ig] / norm;
            (*v)[2][ig] = (*v)[2][ig] / norm;
        }
        else
        {
            (*v)[0][ig] = 1.0;
            (*v)[1][ig] = 0.0;
            (*v)[2][ig] = 0.0;
        }
    }
        
    
    void crossProduct (boost::shared_ptr<boost::multi_array<Real, 2> > & v, const boost::shared_ptr<boost::multi_array<Real, 2> > & v1, const boost::shared_ptr<boost::multi_array<Real, 2> > & v2, const UInt & ig) const
    {
        (*v)[0][ig] = (*v1)[1][ig] * (*v2)[2][ig] - (*v1)[2][ig] * (*v2)[1][ig];
        (*v)[1][ig] = (*v1)[2][ig] * (*v2)[0][ig] - (*v1)[0][ig] * (*v2)[2][ig];
        (*v)[2][ig] = (*v1)[0][ig] * (*v2)[1][ig] - (*v1)[1][ig] * (*v2)[0][ig];
    }
        
    
    
    // Source term : Int { coef * exp(coefExp *(  Ic_iso -3 )) * ( J^(-2/3)* (F : \nabla v) - 1/3 * (Ic_iso / J) * (CofF : \nabla v) ) }
    void  source_P1isoE_Exp (   Real                                coef,
                                Real                                coefExp,
                                const boost::multi_array<Real, 3 >& CofFk,
                                const boost::multi_array<Real, 3 >& Fk,
                                const std::vector<Real>&            Jk,
                                VectorElemental&                    elvec,
                                const CurrentFE&                    fe )
    {
        
        Real s;
        
        for ( UInt icoor = 0; icoor < nDimensions; ++icoor )
        {
            VectorElemental::vector_view vec =  elvec.block ( icoor );
            for ( UInt i = 0; i < fe.nbFEDof(); ++i )
            {
                s = 0.0;
                for ( UInt k = 0; k < nDimensions; ++k )
                {
                    for ( UInt ig = 0; ig < fe.nbQuadPt(); ++ig )
                    {
                        s += exp ( coefExp * ( (*M_I1Ebar)[ ig ] - 3.0 ) ) * (
                        
                            ( 1 - (*M_nAk)[ig] * ((*M_nAk)[ig] + 2) / std::pow( (*M_nAk)[ig] + 1 , 2.0 ) ) *
                        
                            ( pow ( Jk[ig], (-2.0/3.0) ) * Fk[icoor][k][ig] - 1.0/3.0 * ( 1/Jk[ig] ) * (*M_trCisok)[ ig ] * CofFk[icoor][k][ig] ) +
                                                                              
                            ( (*M_nAk)[ig] * ((*M_nAk)[ig] + 2) / std::pow( (*M_nAk)[ig] + 1 , 2.0 ) - (*M_fAk)[ig] * ((*M_fAk)[ig] + 2) / std::pow( (*M_fAk)[ig] + 1 , 2.0 ) ) *
                                                                              
                            (*M_fk)[icoor][ig] * (*M_f0k)[k][ig] +
                                                                              
                            ( (*M_nAk)[ig] * ((*M_nAk)[ig] + 2) / std::pow( (*M_nAk)[ig] + 1 , 2.0 ) - (*M_sAk)[ig] * ((*M_sAk)[ig] + 2) / std::pow( (*M_sAk)[ig] + 1 , 2.0 ) ) *
                          
                            (*M_sk)[icoor][ig] * (*M_s0k)[k][ig]
                                                                              
                            ) *
                        
                            fe.phiDer ( i, k, ig ) * fe.weightDet ( ig );
                        
                    }
                }
                vec ( i ) += s * coef;
            }
        }
    }
    
    
    
    // Source term : Int { coef * exp(coefExp *(  Ic_iso -3 )) * ( J^(-2/3)* (F : \nabla v) - 1/3 * (Ic_iso / J) * (CofF : \nabla v) ) }
    void  source_P4fE_Exp (   Real                                coef,
                             Real                                coefExp,
                             const boost::multi_array<Real, 3 >& CofFk,
                             const boost::multi_array<Real, 3 >& Fk,
                             const std::vector<Real>&            Jk,
                             VectorElemental&                    elvec,
                             const CurrentFE&                    fe )
    {
        
        Real s;
        
        for ( UInt icoor = 0; icoor < nDimensions; ++icoor )
        {
            VectorElemental::vector_view vec =  elvec.block ( icoor );
            for ( UInt i = 0; i < fe.nbFEDof(); ++i )
            {
                s = 0.0;
                for ( UInt k = 0; k < nDimensions; ++k )
                {
                    for ( UInt ig = 0; ig < fe.nbQuadPt(); ++ig )
                    {
                        s += ( (*M_I4fE)[ ig ] - 1.0 ) * exp ( coefExp * ( (*M_I4fE)[ ig ] - 1.0 ) * ( (*M_I4fE)[ ig ] - 1.0 ) ) * (
                                                                              
                            ( 2 / std::pow( (*M_fAk)[ig] + 1 , 2.0 ) ) *
                        
                            (*M_fk)[icoor][ig] * (*M_f0k)[k][ig] *
                            
                            ( (*M_I4fE)[ ig ] > 0. )
                          
                            ) *
                        
                            fe.phiDer ( i, k, ig ) * fe.weightDet ( ig );
                        
                    }
                }
                vec ( i ) += s * coef;
            }
        }
    }

    
    void  source_P4sE_Exp (   Real                                coef,
                           Real                                coefExp,
                           const boost::multi_array<Real, 3 >& CofFk,
                           const boost::multi_array<Real, 3 >& Fk,
                           const std::vector<Real>&            Jk,
                           VectorElemental&                    elvec,
                           const CurrentFE&                    fe )
    {
        
        Real s;
        
        for ( UInt icoor = 0; icoor < nDimensions; ++icoor )
        {
            VectorElemental::vector_view vec =  elvec.block ( icoor );
            for ( UInt i = 0; i < fe.nbFEDof(); ++i )
            {
                s = 0.0;
                for ( UInt k = 0; k < nDimensions; ++k )
                {
                    for ( UInt ig = 0; ig < fe.nbQuadPt(); ++ig )
                    {
                        s += ( (*M_I4sE)[ ig ] - 1.0 ) * exp ( coefExp * ( (*M_I4sE)[ ig ] - 1.0 ) * ( (*M_I4sE)[ ig ] - 1.0 ) ) * (
                                                                                                                                    
                            2 / std::pow( (*M_sAk)[ig] + 1 , 2.0 ) *
                            
                            (*M_sk)[icoor][ig] * (*M_s0k)[k][ig]
                            
                            ) *
                        
                            fe.phiDer ( i, k, ig ) * fe.weightDet ( ig );
                        
                    }
                }
                vec ( i ) += s * coef;
            }
        }
    }
    
    
    void  source_P8fsE_Exp (   Real                                coef,
                           Real                                coefExp,
                           const boost::multi_array<Real, 3 >& CofFk,
                           const boost::multi_array<Real, 3 >& Fk,
                           const std::vector<Real>&            Jk,
                           VectorElemental&                    elvec,
                           const CurrentFE&                    fe )
    {
        
        Real s;
        
        for ( UInt icoor = 0; icoor < nDimensions; ++icoor )
        {
            VectorElemental::vector_view vec =  elvec.block ( icoor );
            for ( UInt i = 0; i < fe.nbFEDof(); ++i )
            {
                s = 0.0;
                for ( UInt k = 0; k < nDimensions; ++k )
                {
                    for ( UInt ig = 0; ig < fe.nbQuadPt(); ++ig )
                    {
                        s += (*M_I8fsE)[ ig ] * exp ( coefExp * (*M_I8fsE)[ ig ] * (*M_I8fsE)[ ig ] ) * (
                                                                                                                                    
                            1 / ( ((*M_fAk)[ig] + 1) * ((*M_sAk)[ig] + 1) ) *
                            
                            ( (*M_fk)[icoor][ig] * (*M_f0k)[k][ig] + (*M_sk)[icoor][ig] * (*M_s0k)[k][ig] )
                            
                            ) *
                        
                            fe.phiDer ( i, k, ig ) * fe.weightDet ( ig );
                        
                    }
                }
                vec ( i ) += s * coef;
            }
        }
    }
    
    
    

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

    
    
    
    
    //! Local stress vector
    boost::scoped_ptr<VectorElemental>                 M_elvecK;
    
    //! Elementary matrices
    boost::scoped_ptr<MatrixElemental>                 M_elmatK;
    
    //! Vector: stiffness non-linear
    boost::shared_ptr<boost::multi_array<Real, 3> > M_CofFk;

    //vectorPtr_Type                         M_stiff;
    
    //! First Piola-Kirchhoff stress tensor
    vectorPtr_Type                                M_FirstPiolaKStress;
    
    //! Local tensors initialization
    boost::shared_ptr<boost::multi_array<Real, 3> > M_Fk;    
    boost::shared_ptr<std::vector<Real> > M_Jack; // J
    boost::shared_ptr<std::vector<Real> > M_trCisok; // I1bar
    boost::shared_ptr<std::vector<Real> > M_trCk; // I1

    
    boost::shared_ptr<boost::multi_array<Real, 2> > M_fk;
    boost::shared_ptr<boost::multi_array<Real, 2> > M_sk;
    boost::shared_ptr<boost::multi_array<Real, 2> > M_nk;

    boost::shared_ptr<boost::multi_array<Real, 2> > M_f0k;
    boost::shared_ptr<boost::multi_array<Real, 2> > M_s0k;
    boost::shared_ptr<boost::multi_array<Real, 2> > M_n0k;
    
    boost::shared_ptr<std::vector<Real> > M_fAk;
    boost::shared_ptr<std::vector<Real> > M_sAk;
    boost::shared_ptr<std::vector<Real> > M_nAk;

    boost::shared_ptr<std::vector<Real> > M_I1Ebar;
    boost::shared_ptr<std::vector<Real> > M_I4f;
    boost::shared_ptr<std::vector<Real> > M_I4fE;
    boost::shared_ptr<std::vector<Real> > M_I4s;
    boost::shared_ptr<std::vector<Real> > M_I4sE;
    boost::shared_ptr<std::vector<Real> > M_I8fs;
    boost::shared_ptr<std::vector<Real> > M_I8fsE;
    

    
    
    class OrthonormalizeVector
    {
    public:
        typedef LifeV::VectorSmall<3> return_Type;
        
        return_Type operator() (const VectorSmall<3>& v)
        {
            return normalize(v, 0);
        }
        
        return_Type operator() (const VectorSmall<3>& v, const VectorSmall<3>& w)
        {
            auto f (v);
            auto s (w);
            
            s = normalize(s, 1);
            s = s - s.dot (f) * f;
            s = normalize(s, 2);
            
            return s;
        }
        
        return_Type normalize(const VectorSmall<3>& v, const UInt& comp)
        {
            auto V (v);
            Real norm = std::sqrt (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
//            if ( norm >= 1e-13 )
            {
                V[0] = v[0] / norm;
                V[1] = v[1] / norm;
                V[2] = v[2] / norm;
            }
//            else
//            {
//                V *= 0.0;
//                V[comp] = 1.0;
//            }
            return V;
        }
        
        OrthonormalizeVector () {}
        ~OrthonormalizeVector () {}
    };
    
    
    class CrossProduct
    {
    public:
        typedef LifeV::VectorSmall<3> return_Type;
        
        return_Type operator() (const LifeV::VectorSmall<3>& v1, const LifeV::VectorSmall<3>& v2)
        {
            VectorSmall<3> v;
            v[0] = v1[1] * v2[2] - v1[2] * v2[1];
            v[1] = v1[2] * v2[0] - v1[0] * v2[2];
            v[2] = v1[0] * v2[1] - v1[1] * v2[0];
            return v;
        }
        
        CrossProduct() {}
        ~CrossProduct() {}
    };
    
    
    class FAInverse
    {
    public:
        typedef LifeV::MatrixSmall<3,3> return_Type;
        
        return_Type operator() (const LifeV::VectorSmall<3>& f0, const LifeV::VectorSmall<3>& s0, const LifeV::VectorSmall<3>& n0, const Real& gf)
        {
            MatrixSmall<3,3> I;
            I(0,0) = 1.; I(0,1) = 0., I(0,2) = 0.;
            I(1,0) = 0.; I(1,1) = 1., I(1,2) = 0.;
            I(2,0) = 0.; I(2,1) = 0., I(2,2) = 1.;
            
            auto gn = 4 * gf;
            auto gs = 1 / ( (gf + 1) * (gn + 1) ) - 1;
            
            MatrixSmall<3,3> FA;
            FA = I - gf/(gf+1) * outerProduct(f0) - gs/(gs+1) * outerProduct(s0) - gn/(gn+1) * outerProduct(n0);
            
            return FA;
        }
        
        return_Type outerProduct(const LifeV::VectorSmall<3>& i0)
        {
            MatrixSmall<3,3> M;
            for (UInt i (0); i < 3; ++i)
            {
                for (UInt j (0); j < 3; ++j)
                {
                    M(i,j) = i0(i) * i0(j);
                }
            }
            return M;
        }
        
        FAInverse() {}
        ~FAInverse() {}
    };
    
    
    class HeavisideFct
    {
    public:
        typedef Real return_Type;
        
        return_Type operator() (const Real& I4f)
        {
            return (I4f > 0. ? 1. : 0.);
        }
        
        HeavisideFct() {}
        ~HeavisideFct() {}
    };
    
    class DeformationGradientReAssembler
    {
    public:
        typedef LifeV::MatrixSmall<3,3> return_Type;
        
        return_Type operator() (const LifeV::VectorSmall<3>& Fx, const LifeV::VectorSmall<3>& Fy, const LifeV::VectorSmall<3>& Fz)
        {
            MatrixSmall<3,3> v;
            v(0,0) = Fx[0];
            v(1,0) = Fx[1];
            v(2,0) = Fx[2];
            v(0,1) = Fy[0];
            v(1,1) = Fy[1];
            v(2,1) = Fy[2];
            v(0,2) = Fz[0];
            v(1,2) = Fz[1];
            v(2,2) = Fz[2];
            return v;
        }
        
        DeformationGradientReAssembler() {}
        ~DeformationGradientReAssembler() {}
    };

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

    this->M_dispFESpace                 = dFESpace;
    this->M_dispETFESpace               = dETFESpace;
    this->M_localMap                    = monolithicMap;
    this->M_offset                      = offset;
    this->M_dataMaterial                = dataMaterial;
    this->M_displayer                   = displayer;

    M_residualVectorPtr.reset ( new vector_Type (*this->M_localMap, Repeated) );
    //   M_identity = EMUtility::identity();

    M_fiberVectorPtr.reset             ( new vector_Type (*this->M_localMap, Repeated) );
    M_sheetVectorPtr.reset             ( new vector_Type (*this->M_localMap, Repeated) );
//    M_fiberVectorPtr.reset             ( new vector_Type (*this->M_localMap, Unique) );
//    M_sheetVectorPtr.reset             ( new vector_Type (*this->M_localMap, Unique) );
    M_scalarETFESpacePtr.reset         ( new scalarETFESpace_Type ( dETFESpace -> mesh(),
    																&( dETFESpace -> refFE() ),
                                                                    dFESpace->map().commPtr() ) );

    M_fiberActivationPtr.reset (new vector_Type (M_scalarETFESpacePtr -> map() ) );

    std::string passiveMaterialType ( dataMaterial -> passiveType() );
    std::string activeStressMaterialType (dataMaterial -> activeStressType() );
    displayer->leaderPrint ("\n===========================");
    displayer->leaderPrint ("\nActive strain Holzapfel-Ogden material created");
    //std::cout << "\nPassive Type: " << passiveMaterialType;
    //std::cout << "\nActive Stress Type: " << activeStressMaterialType;
    displayer->leaderPrint ("\n===========================");

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
    
    
    
    
    //M_stiff.reset                     ( new vector_Type (*this->M_localMap) );
    
    M_FirstPiolaKStress.reset        ( new vector_Type (*this->M_localMap) );
    M_elvecK.reset            ( new VectorElemental (this->M_dispFESpace->fe().nbFEDof(), nDimensions) );
    this->M_elmatK.reset                ( new MatrixElemental ( this->M_dispFESpace->fe().nbFEDof(), nDimensions, nDimensions ) );
    
    //! Local tensors initilization
    M_Fk.reset ( new boost::multi_array<Real, 3> (boost::extents[nDimensions][nDimensions][dFESpace->fe().nbQuadPt()]) );
    M_CofFk.reset ( new boost::multi_array<Real, 3> (boost::extents[nDimensions][nDimensions][dFESpace->fe().nbQuadPt()]) );
    
    M_Jack.reset ( new std::vector<Real> (dFESpace->fe().nbQuadPt(), 0.0) ); // J
    M_trCisok.reset ( new std::vector<Real> (dFESpace->fe().nbQuadPt(), 0.0) ); // I1bar
    M_trCk.reset ( new std::vector<Real> (dFESpace->fe().nbQuadPt(), 0.0) ); // I1

    M_I1Ebar.reset ( new std::vector<Real> (dFESpace->fe().nbQuadPt(), 0.0) );
    M_I4f.reset ( new std::vector<Real> (dFESpace->fe().nbQuadPt(), 0.0) );
    M_I4fE.reset ( new std::vector<Real> (dFESpace->fe().nbQuadPt(), 0.0) );
    M_I4s.reset ( new std::vector<Real> (dFESpace->fe().nbQuadPt(), 0.0) );
    M_I4sE.reset ( new std::vector<Real> (dFESpace->fe().nbQuadPt(), 0.0) );
    M_I8fs.reset ( new std::vector<Real> (dFESpace->fe().nbQuadPt(), 0.0) );
    M_I8fsE.reset ( new std::vector<Real> (dFESpace->fe().nbQuadPt(), 0.0) );

    M_fk.reset ( new boost::multi_array<Real, 2> (boost::extents[nDimensions][dFESpace->fe().nbQuadPt()]) );
    M_sk.reset ( new boost::multi_array<Real, 2> (boost::extents[nDimensions][dFESpace->fe().nbQuadPt()]) );
    M_nk.reset ( new boost::multi_array<Real, 2> (boost::extents[nDimensions][dFESpace->fe().nbQuadPt()]) );
    
    M_f0k.reset ( new boost::multi_array<Real, 2> (boost::extents[nDimensions][dFESpace->fe().nbQuadPt()]) );
    M_s0k.reset ( new boost::multi_array<Real, 2> (boost::extents[nDimensions][dFESpace->fe().nbQuadPt()]) );
    M_n0k.reset ( new boost::multi_array<Real, 2> (boost::extents[nDimensions][dFESpace->fe().nbQuadPt()]) );

    M_fAk.reset ( new std::vector<Real> (dFESpace->fe().nbQuadPt(), 0.0) );
    M_sAk.reset ( new std::vector<Real> (dFESpace->fe().nbQuadPt(), 0.0) );
    M_nAk.reset ( new std::vector<Real> (dFESpace->fe().nbQuadPt(), 0.0) );
    

    //    // The 2 is because the law uses three parameters (mu, bulk).
    //    // another way would be to set up the number of constitutive parameters of the law
    //    // in the data file to get the right size. Note the comment below.
    //    this->M_vectorsParameters.reset ( new vectorsParameters_Type ( 2 ) );
    //
    //    this->setupVectorsParameters();
}

//const std::vector<Vector3D> currentPosition(const VectorEpetra& disp) const
//{
//    Int nLocalDof = disp.blockMap().NumGlobalElements(); //disp.epetraVector().MyLength();
//    Int nComponentLocalDof = nLocalDof / 3;
//    
//    std::vector<Vector3D> boundaryCoordinates(M_boundaryPoints.size());
//    
//    for ( unsigned int i (0) ; i < M_boundaryPoints.size() ; ++i )
//    {
//        int root; int LID;
//        disp.blockMap().RemoteIDList(1, &M_boundaryPoints[i], &root, &LID);
//        
//        if ( disp.blockMap().MyGID( M_boundaryPoints[i] ) )
//        {
//            Vector3D pointCoordinates;
//            
//            UInt iGID = M_boundaryPoints[i];
//            UInt jGID = M_boundaryPoints[i] + nComponentLocalDof;
//            UInt kGID = M_boundaryPoints[i] + 2 * nComponentLocalDof;
//            
//            pointCoordinates[0] = M_fullMesh.point (iGID).x() + disp[iGID];
//            pointCoordinates[1] = M_fullMesh.point (iGID).y() + disp[jGID];
//            pointCoordinates[2] = M_fullMesh.point (iGID).z() + disp[kGID];
//            
//            boundaryCoordinates[i] = pointCoordinates;
//        }
//        
//        MPI_Bcast(&boundaryCoordinates[i], 3, MPI_DOUBLE, root, MPI_COMM_WORLD);
//    }
//    
//    MPI_Barrier(MPI_COMM_WORLD);
//    return boundaryCoordinates;
//}
//    
//    Int p1nCompLocalDof = p1PositionVector.epetraVector().MyLength() / 3;
//    for (int j (0); j < p1nCompLocalDof; j++)
//    {
//        UInt iGID = p1PositionVector.blockMap().GID (j);
//        UInt jGID = p1PositionVector.blockMap().GID (j + p1nCompLocalDof);
//        UInt kGID = p1PositionVector.blockMap().GID (j + 2 * p1nCompLocalDof);
//        
//        p1PositionVector[iGID] = M_fullMesh.point (iGID).x();
//        p1PositionVector[jGID] = M_fullMesh.point (iGID).y();
//        p1PositionVector[kGID] = M_fullMesh.point (iGID).z();
//    }
//
//    
//VectorEpetra pathologic activation ( VectorEpetra& vec, boost::shared_ptr<  RegionMesh<LinearTetra> > fullMesh, Real value, std::vector<UInt> flags)
//{
//    VectorEpetra fiberActivation ( M_fiberActivationPtr );
//
//    for ( int j (0); j < fiberActivation.epetraVector().MyLength() ; ++j )
//    {
//        for ( UInt k (0); k < flags.size(); k++ )
//        {
//            if ( fullMesh -> point ( vec.blockMap().GID (j) ).markerID() == flags.at (k) )
//            {
//                if ( vec.blockMap().LID ( vec.blockMap().GID (j) ) != -1 )
//                {
//                    (vec) ( vec.blockMap().GID (j) ) = value;
//                }
//            }
//        }
//    }
//}

    
    

    
    
    template <typename MeshType>
    void EMStructuralConstitutiveLaw<MeshType>::computeKinematicsVariables ( const VectorElemental& dk_loc )
    {
        
        Real s;
        
        //! loop on quadrature points (ig)
        for ( UInt ig = 0; ig < this->M_dispFESpace->fe().nbQuadPt(); ig++ )
        {
            //! loop on space coordinates (icoor)
            for ( UInt icoor = 0; icoor < nDimensions; icoor++ )
            {
                //! loop  on space coordinates (jcoor)
                for ( UInt jcoor = 0; jcoor < nDimensions; jcoor++ )
                {
                    s = 0.0;
                    for ( UInt i = 0; i < this->M_dispFESpace->fe().nbFEDof(); i++ )
                    {
                        //! \grad u^k at a quadrature point
                        s += this->M_dispFESpace->fe().phiDer ( i, jcoor, ig ) *
                        dk_loc[ i + icoor * this->M_dispFESpace->fe().nbFEDof() ];
                    }
                    //! gradient of displacement
                    (*M_Fk) [ icoor ][ jcoor ][ig ] = s;
                }
            }
        }
        
        //! loop on quadrature points (ig)
        for ( UInt ig = 0; ig < this->M_dispFESpace->fe().nbQuadPt(); ig++ )
        {
            //! loop on space coordinates (icoor)
            for ( UInt  icoor = 0; icoor < nDimensions; icoor++ )
            {
                //! deformation gradient Fk
                (*M_Fk) [ icoor ][ icoor ][ ig ] +=  1.0;
            }
        }
        
        Real a, b, c, d, e, f, g, h, i;
        
        for ( UInt ig = 0; ig < this->M_dispFESpace->fe().nbQuadPt(); ig++ )
        {
            a = (*M_Fk) [ 0 ][ 0 ][ ig ];
            b = (*M_Fk) [ 0 ][ 1 ][ ig ];
            c = (*M_Fk) [ 0 ][ 2 ][ ig ];
            d = (*M_Fk) [ 1 ][ 0 ][ ig ];
            e = (*M_Fk) [ 1 ][ 1 ][ ig ];
            f = (*M_Fk) [ 1 ][ 2 ][ ig ];
            g = (*M_Fk) [ 2 ][ 0 ][ ig ];
            h = (*M_Fk) [ 2 ][ 1 ][ ig ];
            i = (*M_Fk) [ 2 ][ 2 ][ ig ];
            
            //! determinant of deformation gradient Fk
            (*M_Jack) [ig] = a * ( e * i - f * h ) - b * ( d * i - f * g ) + c * ( d * h - e * g );
            
            ASSERT_PRE ( (*M_Jack) [ig] > 0, "Negative Jacobian. Error!" );
            
            (*M_CofFk) [ 0 ][ 0 ][ ig ] =   ( e * i - f * h );
            (*M_CofFk) [ 0 ][ 1 ][ ig ] = - ( d * i - g * f );
            (*M_CofFk) [ 0 ][ 2 ][ ig ] =   ( d * h - e * g );
            (*M_CofFk) [ 1 ][ 0 ][ ig ] = - ( b * i - c * h );
            (*M_CofFk) [ 1 ][ 1 ][ ig ] =   ( a * i - c * g );
            (*M_CofFk) [ 1 ][ 2 ][ ig ] = - ( a * h - g * b );
            (*M_CofFk) [ 2 ][ 0 ][ ig ] =   ( b * f - c * e );
            (*M_CofFk) [ 2 ][ 1 ][ ig ] = - ( a * f - c * d );
            (*M_CofFk) [ 2 ][ 2 ][ ig ] =   ( a * e - d * b );
        }
        
        //! loop on quadrature points
        for ( UInt ig = 0; ig < this->M_dispFESpace->fe().nbQuadPt(); ig++ )
        {
            s = 0.0;
            for ( UInt i = 0; i < nDimensions; i++)
            {
                for ( UInt j = 0; j < nDimensions; j++)
                {
                    //! trace of  C1 = (F1k^t F1k)
                    s +=  (*M_Fk) [ i ][ j ][ ig ] * (*M_Fk) [ i ][ j ][ ig ];
                }
            }
            (*M_trCk) [ ig ] = s;
        }
        
        for ( UInt ig = 0; ig <  this->M_dispFESpace->fe().nbQuadPt(); ig++ )
        {
            //! trace of deviatoric C
            (*M_trCisok) [ ig ] =  pow ( (*M_Jack) [ ig ], -2. / 3.) * (*M_trCk) [ ig ];
        }
        
    }
    
    
    template <typename MeshType>
    void EMStructuralConstitutiveLaw<MeshType>::computeAnisotropyVariables ( const VectorElemental& fk_loc, const VectorElemental& sk_loc, const VectorElemental& fAk_loc )
    {
        
        //=========================================//
        //  f0, s0, n0, f, s, n, fA
        //=========================================//
        
        Real sf, ss;
        
        //! loop on quadrature points (ig)
        for ( UInt ig = 0; ig < this->M_dispFESpace->fe().nbQuadPt(); ig++ )
        {
            //! loop on space coordinates (icoor)
            for ( UInt icoor = 0; icoor < nDimensions; icoor++ )
            {
                sf = 0.0; ss = 0.0;
                for ( UInt i = 0; i < this->M_dispFESpace->fe().nbFEDof(); i++ )
                {
                    //! \grad u^k at a quadrature point
                    sf += this->M_dispFESpace->fe().phi ( i, ig ) * fk_loc[ i + icoor * this->M_dispFESpace->fe().nbFEDof() ];
                    ss += this->M_dispFESpace->fe().phi ( i, ig ) * sk_loc[ i + icoor * this->M_dispFESpace->fe().nbFEDof() ];
                }
                
                (*M_f0k) [ icoor ][ ig ] = sf;
                (*M_s0k) [ icoor ][ ig ] = ss;
                
            }
            
            normalize(M_f0k, ig);
            normalize(M_s0k, ig);
            
            crossProduct(M_n0k, M_f0k, M_s0k, ig);
            
            for ( UInt mcoor = 0; mcoor < nDimensions; mcoor++ )
            {
                (*M_fk) [ mcoor ][ ig ] = 0.0;
                (*M_sk) [ mcoor ][ ig ] = 0.0;
                (*M_nk) [ mcoor ][ ig ] = 0.0;

                for ( UInt ncoor = 0; ncoor < nDimensions; ncoor++ )
                {
                    (*M_fk) [ mcoor ][ ig ] += (*M_Fk) [mcoor] [ncoor] [ig] * (*M_f0k) [ ncoor ][ ig ];
                    (*M_sk) [ mcoor ][ ig ] += (*M_Fk) [mcoor] [ncoor] [ig] * (*M_s0k) [ ncoor ][ ig ];
                    (*M_nk) [ mcoor ][ ig ] += (*M_Fk) [mcoor] [ncoor] [ig] * (*M_n0k) [ ncoor ][ ig ];
                }
            }
            
        }
        
        
        //=========================================//
        //  Orthotropic activation
        //=========================================//

        Real k = 4.0;
        Real sfA;
        
        for ( UInt ig = 0; ig < this->M_dispFESpace->fe().nbQuadPt(); ig++ )
        {
            sfA = 0.0;
            for ( UInt i = 0; i < this->M_dispFESpace->fe().nbFEDof(); i++ )
            {
                sfA += this->M_dispFESpace->fe().phi ( i, ig ) * fAk_loc[ i ];
            }
            
            (*M_fAk) [ ig ] = sfA;
            (*M_nAk) [ ig ] = k * (*M_fAk) [ ig ];
            (*M_sAk) [ ig ] = 1 / ( ((*M_fAk) [ ig ] + 1) * ((*M_nAk) [ ig ] + 1) ) - 1;
        }
        

//        // Coefficients for FAinv
//        auto gm = value(-1.0) * ( gf ) / ( ( gf ) + 1.0 );
//        auto go = gf * ( k + gf * k + value(1.0) );
//        auto gmn = value(-1.0) * ( k*gf ) / ( ( k*gf ) + 1.0 ) ;
        
        
        //=========================================//
        //  Invariants
        //=========================================//
     
        for ( UInt ig = 0; ig < this->M_dispFESpace->fe().nbQuadPt(); ig++ )
        {
            // I4f
            (*M_I4f)[ig] = (*M_fk)[0][ig] * (*M_fk)[0][ig] + (*M_fk)[1][ig] * (*M_fk)[1][ig] + (*M_fk)[2][ig] * (*M_fk)[2][ig];
            (*M_I4fE)[ig] = (*M_I4f)[ig] / std::pow( (*M_fAk)[ig] + 1 , 2.0 );
            
            // I4s
            (*M_I4s)[ig] = (*M_sk)[0][ig] * (*M_sk)[0][ig] + (*M_sk)[1][ig] * (*M_sk)[1][ig] + (*M_sk)[2][ig] * (*M_sk)[2][ig];
            (*M_I4sE)[ig] = (*M_I4s)[ig] / std::pow( (*M_sAk)[ig] + 1 , 2.0 );
            
            // I8fs
            (*M_I8fs)[ig] = (*M_fk)[0][ig] * (*M_sk)[0][ig] + (*M_fk)[1][ig] * (*M_sk)[1][ig] + (*M_fk)[2][ig] * (*M_sk)[2][ig];
            (*M_I8fsE)[ig] = (*M_I4f)[ig] / ( ((*M_fAk)[ig] + 1) * ((*M_sAk)[ig] + 1) );
            
            // I1Ebar
            Real coeffI1bar = 1 - (*M_nAk)[ig] * ((*M_nAk)[ig] + 2) / std::pow( (*M_nAk)[ig] + 1 , 2.0 );
            Real coeffI4f = (*M_nAk)[ig] * ((*M_nAk)[ig] + 2) / std::pow( (*M_nAk)[ig] + 1 , 2.0 ) - (*M_fAk)[ig] * ((*M_fAk)[ig] + 2) / std::pow( (*M_fAk)[ig] + 1 , 2.0 );
            Real coeffI4s = (*M_nAk)[ig] * ((*M_nAk)[ig] + 2) / std::pow( (*M_nAk)[ig] + 1 , 2.0 ) - (*M_sAk)[ig] * ((*M_sAk)[ig] + 2) / std::pow( (*M_sAk)[ig] + 1 , 2.0 );
            (*M_I1Ebar)[ig] = coeffI1bar * (*M_trCisok)[ig] + coeffI4f * (*M_I4f)[ig] + coeffI4s * (*M_I4s)[ig];
        }
        
    }
    

template <typename MeshType>
void EMStructuralConstitutiveLaw<MeshType>::computeStiffness2 ( const vector_Type& sol,
                                                           Real /*factor*/,
                                                           const dataPtr_Type& dataMaterial,
                                                           const mapMarkerVolumesPtr_Type mapsMarkerVolumes,
                                                           const displayerPtr_Type& displayer )
{
    //this->M_stiff.reset (new vector_Type (*this->M_localMap) );
    
    displayer->leaderPrint (" \n*********************************\n  ");
    displayer->leaderPrint (" Non-Linear S- computeStiffness2 ");
    displayer->leaderPrint (" \n*********************************\n  ");
    
    displayer->leaderPrint ( this->M_dispFESpace->qr() );
    displayer->leaderPrint ( this->M_dispFESpace->fe().nbQuadPt() );
    displayer->leaderPrint ( this->M_dispFESpace->fe().nbFEDof() );
    displayer->leaderPrint ( this->M_dispFESpace->fe().quadRule() );

    
    
    UInt totalDof   = this->M_dispFESpace->dof().numTotalDof();
    UInt dim = this->M_dispFESpace->dim();
    
    VectorElemental dk_loc ( this->M_dispFESpace->fe().nbFEDof(), nDimensions );
    VectorElemental fk_loc ( this->M_dispFESpace->fe().nbFEDof(), nDimensions );
    VectorElemental sk_loc ( this->M_dispFESpace->fe().nbFEDof(), nDimensions );
    VectorElemental fAk_loc ( this->M_dispFESpace->fe().nbFEDof(), 1 );

    vector_Type dRep (sol, Repeated);
    vector_Type fRep (*M_fiberVectorPtr, Repeated);
    vector_Type sRep (*M_sheetVectorPtr, Repeated);
    vector_Type fARep (*M_fiberActivationPtr, Repeated);
    
    const UInt nbElements (this->M_dispFESpace->mesh()->numElements() );
//    const UInt fieldDim (this->M_dispFESpace->fieldDim() );
//    const UInt nbTotalDof (this->M_dispFESpace->dof().numTotalDof() );
    
    
//    mapIterator_Type it;
    
//    for ( it = (*mapsMarkerVolumes).begin(); it != (*mapsMarkerVolumes).end(); it++ )
//    {
//        
//        //Given the marker pointed by the iterator, let's extract the material parameters
//        UInt marker = it->first;
    
    for (UInt iterElement (0); iterElement < nbElements; ++iterElement)
    {
        
//        M_diffCFE->update ( M_fespace->mesh()->element (iterElement), UPDATE_DPHI | UPDATE_WDET );
        
        Real bulk = 3500000; //dataMaterial->bulk (marker);
        Real alpha = 3300; //dataMaterial->alpha (marker);

        Real gamma = 9.242; //dataMaterial->gamma (marker);
//        std::cout << "iterator: " << it->first << "\t" << it->second.size() << std::endl;
//        for ( UInt j (0); j < it->second.size(); j++ )
//        {
//            this->M_dispFESpace->fe().updateFirstDerivQuadPt ( * (it->second[j]) );
        
            //this->M_dispFESpace->fe().update ( this->M_dispFESpace->mesh()->element (iterElement), UPDATE_DPHI | UPDATE_WDET | UPDATE_PHI_VECT);
            this->M_dispFESpace->fe().updateFirstDerivQuadPt ( this->M_dispFESpace->mesh()->element (iterElement) );

            UInt eleID = this->M_dispFESpace->fe().currentLocalId();
        
            //std::cout << "elem: " << eleID << "\t" << iterElement << std::endl;

            for ( UInt iNode = 0 ; iNode < ( UInt ) this->M_dispFESpace->fe().nbFEDof() ; iNode++ )
            {
                UInt  iloc = this->M_dispFESpace->fe().patternFirst ( iNode );
                
                for ( UInt iComp = 0; iComp < nDimensions; ++iComp )
                {
                    UInt ig = this->M_dispFESpace->dof().localToGlobalMap ( eleID, iloc ) + iComp * dim + this->M_offset;
                    dk_loc[ iloc + iComp * this->M_dispFESpace->fe().nbFEDof() ] = dRep[ig];
                    fk_loc[ iloc + iComp * this->M_dispFESpace->fe().nbFEDof() ] = fRep[ig];
                    sk_loc[ iloc + iComp * this->M_dispFESpace->fe().nbFEDof() ] = sRep[ig];
                }
                
                fAk_loc[iloc] = fARep[this->M_dispFESpace->dof().localToGlobalMap ( eleID, iloc ) + this->M_offset];
            }
            
            this->M_elvecK->zero();
            
            this->computeKinematicsVariables ( dk_loc );

            this->computeAnisotropyVariables ( fk_loc , sk_loc , fAk_loc );

        
            //! Stiffness for non-linear terms of the Neo-Hookean model
            /*!
             The results of the integrals are stored at each step into elvecK, until to build K matrix of the bilinear form
             */
            //! Volumetric part
            /*!
             Source term Pvol: int { bulk /2* (J1^2 - J1  + log(J1) ) * 1/J1 * (CofF1 : \nabla v) }
             */
            //AssemblyElementalStructure::source_Pvol ( 0.5 * bulk, (*M_CofFk), (*M_Jack), *this->M_elvecK,  this->M_dispFESpace->fe() );
            
            //! Isochoric part
            /*!
             Source term P1iso_Exp: int { alpha * exp(gamma *(  Ic1_iso -3 )) *
             ( J1^(-2/3)* (F1 : \nabla v) - 1/3 * (Ic1_iso / J1) * (CofF1 : \nabla v) ) }
             */
        
        
        
            // M_fk, M_sk, M_I1Ebar, gamman, gammaf, gammas
        
            //source_P1isoE_Exp ( 3300 , 9.242 , (*M_CofFk), (*M_Fk), (*M_Jack), *this->M_elvecK, this->M_dispFESpace->fe() );
            
            //source_P4fE_Exp ( 185350, 15.972, (*M_CofFk), (*M_Fk), (*M_Jack), *this->M_elvecK, this->M_dispFESpace->fe() );

            //source_P4sE_Exp ( 25640, 10.446, (*M_CofFk), (*M_Fk), (*M_Jack), *this->M_elvecK, this->M_dispFESpace->fe() );

            //source_P8fsE_Exp ( 4170.0, 11.602, (*M_CofFk), (*M_Fk), (*M_Jack), *this->M_elvecK, this->M_dispFESpace->fe() );

        
            for ( UInt ic = 0; ic < nDimensions; ++ic )
            {
                /*!
                 M_elvecK is assemble into *vec_stiff vector that is recall
                 from updateSystem(matrix_ptrtype& mat_stiff, vector_ptr_type& vec_stiff)
                 */
                assembleVector ( *this->M_residualVectorPtr,
                                *this->M_elvecK,
                                this->M_dispFESpace->fe(),
                                this->M_dispFESpace->dof(), ic, this->M_offset +  ic * totalDof );
            }
        
    }
    
    this->M_residualVectorPtr->globalAssemble();
}

    

    
    
template <typename MeshType>
void EMStructuralConstitutiveLaw<MeshType>::updateJacobianMatrix ( const vector_Type&       disp,
                                                                   const dataPtr_Type&      dataMaterial,
                                                                   const mapMarkerVolumesPtr_Type mapsMarkerVolumes,
                                                                   const mapMarkerIndexesPtr_Type mapsMarkerIndexes,
                                                                   const displayerPtr_Type& displayer )
{
    this->M_jacobian.reset (new matrix_Type (*this->M_localMap) );

    * (this->M_jacobian) *= 0.0;
    
    MatrixSmall<3,3> I;
    I(0,0) = 1.; I(0,1) = 0., I(0,2) = 0.;
    I(1,0) = 0.; I(1,1) = 1., I(1,2) = 0.;
    I(2,0) = 0.; I(2,1) = 0., I(2,2) = 1.;
    

    //std::vector<VectorEpetra> defF = computeGlobalDeformationGradientVector(this->M_dispFESpace, disp);

    
    boost::shared_ptr<HeavisideFct> heaviside (new HeavisideFct);
    boost::shared_ptr<CrossProduct> crossProduct (new CrossProduct);
    boost::shared_ptr<OrthonormalizeVector> orthonormalizeVector (new OrthonormalizeVector);
    //boost::shared_ptr<DeformationGradientReAssembler> defGReAssembler (new DeformationGradientReAssembler);

    LifeChrono chrono;
    chrono.start();

    {
        using namespace ExpressionAssembly;
        
//#define deformationGradientTensor ( value(I) + grad(super::M_dispETFESpace, disp, 0) )

        
        //todo: FE ganz ausgeschrieben, nur einmal integrate
        
//        auto Fx = value(super::M_dispETFESpace, defF[0]);
//        auto Fy = value(super::M_dispETFESpace, defF[1]);
//        auto Fz = value(super::M_dispETFESpace, defF[2]);
//        auto F = eval(defGReAssembler, Fx, Fy, Fz);
        
        auto F = value(I) + grad(super::M_dispETFESpace, disp, 0);
        
        auto dF = grad(phi_j);
        auto FmT = minusT(F);
        auto J = det(F);
        auto dJ = J * FmT;
        auto Jm23 = pow(J, 2 / (-3.) );
        auto I1 = dot(F, F);
        auto dI1bar = value(2.0) * Jm23 * ( F + value(1/(-3.)) * I1 * FmT );
        
        
        // Anisotropy
        auto f_0 = value (super::M_dispETFESpace, *M_fiberVectorPtr);
        auto s_0 = value (super::M_dispETFESpace, *M_sheetVectorPtr);
        auto f0 = eval (orthonormalizeVector, f_0);
        auto s0 = eval (orthonormalizeVector, f0, s_0);
        auto n0 = eval (crossProduct, f0, s0);
        auto f = F * f0;
        auto s = F * s0;

        
        // Orthotropic activation
        auto k = 4.0;
        auto gf = value (M_scalarETFESpacePtr, *M_fiberActivationPtr);
        auto gn = k * gf;
        auto gs = 1 / ( (gf + 1) * (gn + 1) ) - 1;
        auto gm = value(-1.0) * ( gf ) / ( ( gf ) + 1.0 );
        auto go = gf * ( k + gf * k + value(1.0) );
        auto gmn = value(-1.0) * ( k*gf ) / ( ( k*gf ) + 1.0 ) ;
        
        
        // Active strain
        auto FAinv = I + (value(-1.0) * ( value (M_scalarETFESpacePtr, *M_fiberActivationPtr) ) / ( ( value (M_scalarETFESpacePtr, *M_fiberActivationPtr) ) + 1.0 )) * outerProduct(eval (orthonormalizeVector, value (super::M_dispETFESpace, *M_fiberVectorPtr)), eval (orthonormalizeVector, value (super::M_dispETFESpace, *M_fiberVectorPtr))) + (value (M_scalarETFESpacePtr, *M_fiberActivationPtr) * ( 4.0 + value (M_scalarETFESpacePtr, *M_fiberActivationPtr) * 4.0 + value(1.0) )) * outerProduct(eval (orthonormalizeVector, f0,  value (super::M_dispETFESpace, *M_sheetVectorPtr)), eval (orthonormalizeVector, f0,  value (super::M_dispETFESpace, *M_sheetVectorPtr))) + (value(-1.0) * ( 4.0*value (M_scalarETFESpacePtr, *M_fiberActivationPtr) ) / ( ( 4.0*value (M_scalarETFESpacePtr, *M_fiberActivationPtr) ) + 1.0 )) * outerProduct(eval (crossProduct, f0, s0), eval (crossProduct,  eval (orthonormalizeVector, value (super::M_dispETFESpace, *M_fiberVectorPtr)), eval (orthonormalizeVector, f0,  value (super::M_dispETFESpace, *M_sheetVectorPtr))));
        //auto FAinv = I + gm * outerProduct(f0, f0) + go * outerProduct(s0, s0) + gmn * outerProduct(n0, n0);
        auto FAinv = eval(FAInverse, f0, s0, n0, gf, gs, gn);
        auto FE =  F * FAinv;
        auto dFE = dF * FAinv;
        auto FEmT = minusT(FE);
        auto JE = det(FE);
        auto dJE = JE * FEmT;
        auto JEm23 = pow(JE, 2 / (-3.) );
        auto I1E = dot(FE, FE);
        auto dI1E = 2 * FE;
        auto I1barE = pow ( det(FE), 2 / -3.0 ) *  dot( FE, FE );
        auto dI1barE = value(2.0) * JEm23 * ( FE + value(1/(-3.)) * I1E * FEmT );

        
        // Pvol
        auto dJdF = dot(dJ, dF);
        auto dFT = transpose( dF );
        auto dFmTdF = value (-1.0) * FmT * dFT * FmT;
        auto d2JdF = dJdF * FmT + J * dFmTdF;
        auto dWvol = ( 3500000 * ( J + J * log(J) - 1. ) ) / ( 2 * J );
        auto dPvol = dWvol * d2JdF;
        
        auto ddWvol = ( 3500000 * ( J + 1. ) ) / ( 2. * J * J );
        auto ddPvol = ddWvol * dJdF * dJ;

        
        // P1E
        auto dJEm23 = value(-2.0/3.0) * JEm23 * FEmT;
        auto dJEm23dFE = dot(dJEm23, dFE);
        auto dI1EdFE = dot(dI1E, dFE);
        auto d2I1EdFE = 2 * dFE;
        auto dFEmTdFE = value (-1.0) * FEmT * transpose(dFE) * FEmT;
        auto d2JEm23dFE = value(-2.0/3.0) * ( JEm23 *  dFEmTdFE + dJEm23dFE * FEmT );
        auto d2I1barEdFE = dJEm23dFE * dI1E + JEm23 * d2I1EdFE + I1E * d2JEm23dFE + dI1EdFE * dJEm23;
        auto dW1E = 3300 / 2.0 * exp ( 9.242 * ( I1barE - 3 ) );
        auto dP1E = dW1E * d2I1barEdFE * FAinv;
        
        auto dI1barEdFE = dot( dI1barE, dFE );
        auto ddW1E = 3300 * 9.242 / 2.0 * exp ( 9.242 * ( I1barE - 3 ) );
        auto ddP1E = ddW1E * dI1barEdFE * dI1barE * FAinv;
        this->M_displayer->leaderPrint ("\nIntegrate P1E in \n");
        integrate ( elements ( super::M_dispETFESpace->mesh() ) ,
                   quadRuleTetra4pt,
                   super::M_dispETFESpace,
                   super::M_dispETFESpace,
                   dot ( (
                         3300 / 2.0 * exp ( 9.242 * ( pow ( det(FE), 2 / -3.0 ) *  dot( FE, FE ) - 3 ) ) *
                         ( dot(value(-2.0/3.0) * pow(det(FE), 2 / (-3.) ) * minusT(FE), grad(phi_j)*FAinv) * 2*FE + pow(det(FE), 2 / (-3.) ) * 2*grad(phi_j)*FAinv + dot(FE, FE) * value(-2.0/3.0) * ( pow(det(FE), 2 / (-3.) ) *  value (-1.0) * minusT(FE) * transpose(grad(phi_j)*FAinv) * minusT(FE) + dot(value(-2.0/3.0) * pow(det(FE), 2 / (-3.) ) * minusT(FE), grad(phi_j)*FAinv) * minusT(FE) ) + dot(2*FE, grad(phi_j)*FAinv) *  value(-2.0/3.0) * pow(det(FE), 2 / (-3.) ) * minusT(FE) )
                         +
                         3300 * 9.242 / 2.0 * exp ( 9.242 * ( pow ( det(FE), 2 / -3.0 ) *  dot( FE, FE ) - 3 ) ) *
                         dot(  value(2.0) * pow(det(FE), 2 / (-3.) ) * ( FE + value(1/(-3.)) * dot(FE, FE) *  minusT(FE) ), grad(phi_j)*FAinv ) *
                         value(2.0) * pow(det(FE), 2 / (-3.) ) * ( FE + value(1/(-3.)) * dot(FE, FE) * minusT(FE) )
                         ) *
                         FAinv
                        
                        , grad (phi_i) )
                   ) >> this->M_jacobian;
        
        this->M_displayer->leaderPrint ("\ndone in ", chrono.diff(),"\n");
        
        // P4fE
        auto I4fE = dot (f,f) / pow (gf + 1, 2.0);
        auto I4m1fE = I4fE - 1.0;
        auto dW4fE = 185350 * I4m1fE * exp (15.972 * I4m1fE * I4m1fE ) * eval(heaviside, I4m1fE);
        auto d2I4fEdFE = value(2.0) * outerProduct( dFE * f0, f0 );
        auto dP4fE = dW4fE * d2I4fEdFE * FAinv;
        
        auto dI4fE = value(2.0) * outerProduct( FE*f0, f0 );
        auto dI4fEdFE =  dot ( dI4fE, dFE );
        auto ddW4fE = 185350 * exp ( 15.972 * I4m1fE * I4m1fE ) * ( 1.0 + 2.0 * 15.972 * I4m1fE * I4m1fE ) * eval(heaviside, I4m1fE);
        auto ddP4fE = ddW4fE * dI4fEdFE * dI4fE * FAinv;

        
//        integrate ( elements ( super::M_dispETFESpace->mesh() ) ,
//                   quadRuleTetra4pt,
//                   super::M_dispETFESpace,
//                   super::M_dispETFESpace,
//                   dot ( ( 185350 * (dot (f,f) / pow (gf + 1, 2.0) - 1.0) * exp (15.972 * (dot (f,f) / pow (gf + 1, 2.0) - 1.0) * (dot (f,f) / pow (gf + 1, 2.0) - 1.0) ) * eval(heaviside, (dot (f,f) / pow (gf + 1, 2.0) - 1.0)) *  value(2.0) * outerProduct( grad(phi_j)*FAinv * f0, f0 ) +
//                       185350 * exp ( 15.972 * (dot (f,f) / pow (gf + 1, 2.0) - 1.0) * (dot (f,f) / pow (gf + 1, 2.0) - 1.0) ) * ( 1.0 + 2.0 * 15.972 * (dot (f,f) / pow (gf + 1, 2.0) - 1.0) * (dot (f,f) / pow (gf + 1, 2.0) - 1.0) ) * eval(heaviside, (dot (f,f) / pow (gf + 1, 2.0) - 1.0)) * dot ( value(2.0) * outerProduct( (F*FAinv)*f0, f0 ), grad(phi_j)*FAinv ) * value(2.0) * outerProduct( (F*FAinv)*f0, f0 )
//                       ) * FAinv
//                        
//                        , grad (phi_i) )
//                   ) >> this->M_jacobian;
        
        
        // P4sE
        auto I4sE = dot (s,s) / pow (gs + 1, 2.0);
        auto I4m1sE = I4sE - 1.0;
        auto dW4sE = 25640 * I4m1sE * exp (10.446 * I4m1sE * I4m1sE ) * eval(heaviside, I4m1sE);
        auto d2I4sEdFE = value(2.0) * outerProduct( dFE * s0, s0 );
        auto dP4sE = dW4sE * d2I4sEdFE * FAinv;

        auto dI4sE = value(2.0) * outerProduct( FE*s0, s0 );
        auto dI4sEdFE =  dot ( dI4sE, dFE );
        auto ddW4sE = 25640 * exp ( 10.446 * I4m1sE * I4m1sE ) * ( 1.0 + 2.0 * 10.446 * I4m1sE * I4m1sE ) * eval(heaviside, I4m1sE);
        auto ddP4sE = ddW4sE * dI4sEdFE * dI4sE * FAinv;
        
//        integrate ( elements ( super::M_dispETFESpace->mesh() ) ,
//                   quadRuleTetra4pt,
//                   super::M_dispETFESpace,
//                   super::M_dispETFESpace,
//                   dot ( ( 25640 *  (dot (s,s) / pow (gs + 1, 2.0)-1) * exp (10.446 *  (dot (s,s) / pow (gs + 1, 2.0)-1) *  (dot (s,s) / pow (gs + 1, 2.0)-1) ) * eval(heaviside,  (dot (s,s) / pow (gs + 1, 2.0)-1)) * value(2.0) * outerProduct( (grad(phi_j)*FAinv) * s0, s0 )
//                        +
//                         25640 * exp ( 10.446 *  (dot (s,s) / pow (gs + 1, 2.0)-1) *  (dot (s,s) / pow (gs + 1, 2.0)-1) ) * ( 1.0 + 2.0 * 10.446 *  (dot (s,s) / pow (gs + 1, 2.0)-1) *  (dot (s,s) / pow (gs + 1, 2.0)-1) ) * eval(heaviside,  (dot (s,s) / pow (gs + 1, 2.0)-1)) *  dot (  value(2.0) * outerProduct( (F*FAinv)*s0, s0 ) , grad(phi_j)*FAinv ) * value(2.0) * outerProduct( (F*FAinv)*s0, s0 )
//                        
//                        ) * FAinv
//                        
//                        , grad (phi_i) )
//                   ) >> this->M_jacobian;


        
        // P8fsE
        auto I8fsE = dot (f,s) / ( (gf + 1) * (gs + 1) );
        auto dW8fsE = 4170 * I8fsE * exp ( 11.602 * I8fsE * I8fsE );
        auto d2I8EdFE = dFE * ( outerProduct( f0, s0 ) + outerProduct( s0, f0 ) );
        auto dP8fsE = dW8fsE * d2I8EdFE * FAinv;

        auto dI8E = FE * ( outerProduct( f0, s0 ) + outerProduct( s0, f0 ) );
        auto dI8EdFE = dot ( dI8E , dFE );
        auto ddW8fsE = 4170.0 * exp ( 11.602 * I8fsE * I8fsE ) * ( 2.0 * 11.602 * I8fsE * I8fsE + 1.0 );
        auto ddP8fsE = ddW8fsE * dI8EdFE * dI8E * FAinv;
        
        
//        integrate ( elements ( super::M_dispETFESpace->mesh() ) ,
//                   quadRuleTetra4pt,
//                   super::M_dispETFESpace,
//                   super::M_dispETFESpace,
//                   dot ( ( 4170 * dot (f,s) / ( (gf + 1) * (gs + 1) ) * exp ( 11.602 * dot (f,s) / ( (gf + 1) * (gs + 1) ) * dot (f,s) / ( (gf + 1) * (gs + 1) ) ) *  grad(phi_j)*FAinv * ( outerProduct( f0, s0 ) + outerProduct( s0, f0 ) )
//                        +
//                         4170.0 * exp ( 11.602 * dot (f,s) / ( (gf + 1) * (gs + 1) ) * dot (f,s) / ( (gf + 1) * (gs + 1) ) ) * ( 2.0 * 11.602 * dot (f,s) / ( (gf + 1) * (gs + 1) ) * dot (f,s) / ( (gf + 1) * (gs + 1) ) + 1.0 ) * dot (  (F*FAinv) * ( outerProduct( f0, s0 ) + outerProduct( s0, f0 ) ) , grad(phi_j)*FAinv ) * (F*FAinv) * ( outerProduct( f0, s0 ) + outerProduct( s0, f0 ) )
//                        
//                        ) * FAinv
//                        
//                        , grad (phi_i) )
//                   ) >> this->M_jacobian;

        this->M_displayer->leaderPrint ("\nIntegrate total in \n");
        // Sum up contributions and integrate
        auto dP = dPvol + ddPvol /*+ dP1E + ddP1E*/ + dP4fE + ddP4fE + dP4sE + ddP4sE + dP8fsE + ddP8fsE;
        integrate ( elements ( super::M_dispETFESpace->mesh() ) ,
                   quadRuleTetra4pt,
                   super::M_dispETFESpace,
                   super::M_dispETFESpace,
                   dot ( dP, grad (phi_i) )
                   ) >> this->M_jacobian;
        this->M_displayer->leaderPrint ("\ndone in ", chrono.diff(),"\n");
        
    }
    
    this->M_jacobian->globalAssemble();
}

template <typename MeshType>
void EMStructuralConstitutiveLaw<MeshType>::computeStiffness ( const vector_Type&       disp,
                                                               Real                     factor,
                                                               const dataPtr_Type&      dataMaterial,
                                                               const mapMarkerVolumesPtr_Type mapsMarkerVolumes,
                                                               const mapMarkerIndexesPtr_Type mapsMarkerIndexes,
                                                               const displayerPtr_Type& displayer )
{
    
    * (M_residualVectorPtr) *= 0.0;
    
    LifeChrono chrono;
    chrono.start();
    //computeStiffness2(disp, factor, dataMaterial, mapsMarkerVolumes, displayer);
    chrono.stop();
    this->M_displayer->leaderPrint ("computeStiffness2 function computed in ", chrono.globalDiff ( *this->M_displayer->comm() ), " s\n" );
    
    MatrixSmall<3,3> I;
    I(0,0) = 1.; I(0,1) = 0., I(0,2) = 0.;
    I(1,0) = 0.; I(1,1) = 1., I(1,2) = 0.;
    I(2,0) = 0.; I(2,1) = 0., I(2,2) = 1.;
    
    
    boost::shared_ptr<HeavisideFct> heaviside (new HeavisideFct);
    boost::shared_ptr<CrossProduct> crossProduct (new CrossProduct);
    boost::shared_ptr<OrthonormalizeVector> orthonormalizeVector (new OrthonormalizeVector);

    
    {
        using namespace ExpressionAssembly;
        
        
        auto F = I + grad(super::M_dispETFESpace, disp, 0);
        auto J = det(F);
        auto Jm23 = pow(J, 2 / (-3.) );
        auto FmT = minusT(F);
        auto I1 = dot(F, F);
        auto dI1bar = value(2.0) * Jm23 * ( F + value(1/(-3.)) * I1 * FmT );

        
        // Anisotropy
//        auto f_0 = value (super::M_dispETFESpace, *M_fiberVectorPtr);
//        auto s_0 = value (super::M_dispETFESpace, *M_sheetVectorPtr);
//        auto f0 = eval (orthonormalizeVector, f_0);
//        auto s0 = eval (orthonormalizeVector, f0, s_0);
//        auto n0 = eval (crossProduct, f0, s0);
//        auto f = F * f0;
//        auto s = F * s0;
        
        
        // Orthotropic activation
        auto k = 4.0;
        auto gf = value (M_scalarETFESpacePtr, *M_fiberActivationPtr);
        auto gn = k * gf;
        auto gs = 1 / ( (gf + 1) * (gn + 1) ) - 1;
        auto gm = value(-1.0) * ( gf ) / ( ( gf ) + 1.0 );
        auto go = gf * ( k + gf * k + value(1.0) );
        auto gmn = value(-1.0) * ( k*gf ) / ( ( k*gf ) + 1.0 ) ;
        
        
        // Active strain
        auto f0 = eval (orthonormalizeVector, value (super::M_dispETFESpace, *M_fiberVectorPtr));
        auto s0 = eval (orthonormalizeVector, f0,  value (super::M_dispETFESpace, *M_sheetVectorPtr));
        auto n0 = eval (crossProduct, f0, s0);
        auto f = F * f0;
        auto s = F * s0;

        auto FAinv = I + (value(-1.0) * ( value (M_scalarETFESpacePtr, *M_fiberActivationPtr) ) / ( ( value (M_scalarETFESpacePtr, *M_fiberActivationPtr) ) + 1.0 )) * outerProduct(eval (orthonormalizeVector, value (super::M_dispETFESpace, *M_fiberVectorPtr)), eval (orthonormalizeVector, value (super::M_dispETFESpace, *M_fiberVectorPtr))) + (value (M_scalarETFESpacePtr, *M_fiberActivationPtr) * ( 4.0 + value (M_scalarETFESpacePtr, *M_fiberActivationPtr) * 4.0 + value(1.0) )) * outerProduct(eval (orthonormalizeVector, f0,  value (super::M_dispETFESpace, *M_sheetVectorPtr)), eval (orthonormalizeVector, f0,  value (super::M_dispETFESpace, *M_sheetVectorPtr))) + (value(-1.0) * ( 4.0*value (M_scalarETFESpacePtr, *M_fiberActivationPtr) ) / ( ( 4.0*value (M_scalarETFESpacePtr, *M_fiberActivationPtr) ) + 1.0 )) * outerProduct(eval (crossProduct, f0, s0), eval (crossProduct,  eval (orthonormalizeVector, value (super::M_dispETFESpace, *M_fiberVectorPtr)), eval (orthonormalizeVector, f0,  value (super::M_dispETFESpace, *M_sheetVectorPtr))));
        //auto FAinv = I + gm * outerProduct(f0, f0) + go * outerProduct(s0, s0) + gmn * outerProduct(n0, n0);
        auto FE =  F * FAinv;
        
        
        // Pvol
        auto dWvol = ( 3500000 * ( J + J * log(J) - 1. ) ) / ( 2 * J );
        auto dJ = det(F) * minusT(F);
        auto Pvol = dWvol * dJ;
        

        // P1E
        auto I1barE = pow ( det(FE), 2 / -3.0 ) *  dot( FE, FE );
        auto dI1barE = pow ( det(FE), 2 / -3.0 ) * ( value(2.0) * FE + dot( FE, FE ) * value(-2.0/3.0) * minusT(FE) );
        auto dWI1E = 3300 / 2.0 * exp ( 9.242 * ( I1barE - 3 ) );
        auto P1E = dWI1E * dI1barE * FAinv;
        
        
        // P4fE
        auto I4fE = dot (f,f) / pow (gf + 1, 2.0);
        auto I4m1fE = I4fE - 1.0;
        auto dW4fE = 185350 * I4m1fE * exp (15.972 * I4m1fE * I4m1fE ) * eval(heaviside, I4m1fE);
        auto dI4fE = pow(gf + 1, -2.0);
        auto dI4f = value(2.0) * outerProduct( f, f0 );
        auto P4fE = dW4fE * dI4fE * dI4f;
 
        
        // P4sE
        auto I4sE = dot (s,s) / pow (gs + 1, 2.0);
        auto I4m1sE = I4sE - 1.0;
        auto dW4sE = 25640 * I4m1sE * exp (10.446 * I4m1sE * I4m1sE ) * eval(heaviside, I4m1sE);
        auto dI4sE = pow(gs + 1, -2.0);
        auto dI4s = value(2.0) * outerProduct( s, s0 );
        auto P4sE = dW4sE * dI4sE * dI4s;

        
        // P8fsE
        auto I8fsE = dot (f,s) / ( (gf + 1) * (gs + 1) );
        auto dW8fsE = 4170 * I8fsE * exp ( 11.602 * I8fsE * I8fsE );
        auto dI8fsE = 1 / ( (gf + 1) * (gs + 1) );
        auto dI8fs = F * ( outerProduct( f0, s0 ) + outerProduct( s0, f0 ) );
        auto P8fsE = dW8fsE * dI8fsE * dI8fs;
        
        
        // Sum up contributions and integrate
        auto P = Pvol + P1E + P4fE + P4sE + P8fsE;
        integrate ( elements ( super::M_dispETFESpace->mesh() ) ,
                   quadRuleTetra4pt,
                   super::M_dispETFESpace,
                   dot ( P, grad (phi_i) )
                   ) >> M_residualVectorPtr;
    
    }

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
