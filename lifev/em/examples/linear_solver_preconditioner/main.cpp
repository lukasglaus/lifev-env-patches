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
    @file
    @brief

    @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @date 30-03-2011
 */


#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_RCP.hpp>


#include <lifev/core/LifeV.hpp>
#include <lifev/core/mesh/RegionMesh3DStructured.hpp>
#include <lifev/core/mesh/MeshData.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>
#include <lifev/core/mesh/MeshUtility.hpp>
#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/fem/BCManage.hpp>
#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/solver/ADRAssembler.hpp>
#include <lifev/core/algorithm/PreconditionerLinearSolver.hpp>
#include <lifev/core/algorithm/SolverAztecOO.hpp>
#include <lifev/core/algorithm/LinearSolver.hpp>
#include <lifev/core/function/Laplacian.hpp>

#include <lifev/core/filter/ExporterHDF5.hpp>

#include <lifev/structure/fem/AssemblyElementalStructure.hpp>



#define TEST_TOLERANCE 1e-13

using namespace LifeV;

namespace
{
typedef RegionMesh<LinearTetra>           mesh_Type;
typedef boost::shared_ptr<mesh_Type>      meshPtr_Type;
typedef MatrixEpetra<Real>                matrix_Type;
typedef boost::shared_ptr<matrix_Type>    matrix_ptrType;
typedef VectorEpetra                      vector_Type;
typedef boost::shared_ptr<VectorEpetra>   vectorPtr_Type;
typedef FESpace< mesh_Type, MapEpetra >   fespace_Type;
typedef boost::shared_ptr< fespace_Type > fespacePtr_Type;

typedef LifeV::Preconditioner             basePrec_Type;
typedef boost::shared_ptr<basePrec_Type>  basePrecPtr_Type;
typedef LifeV::PreconditionerLinearSolver prec_Type;
typedef boost::shared_ptr<prec_Type>      precPtr_Type;
typedef boost::function < Real ( Real const&,
                                 Real const&,
                                 Real const&,
                                 Real const&,
                                 UInt const& ) > function_Type;
}



//void jacobian ( vectorPtr_Type displacement, boost::shared_ptr<matrix_Type> jacobian , fespacePtr_Type M_fespace)
//{
////    //Initialization of the deformationF matrix
////    M_deformationF.reset  ( new matrixSerialDense_Type ( 3,3 ) );
////    
////    matrix_Type Jacobian ( jacobian );
////    Jacobian *= 0.0;
////    
////
////
////    //Loop over the volumes to compute J = det(F)
////    //Inside the loop, the determinant is store in the appropriate positions
////    vector_Type dRep (*displacement, Repeated);
////    UInt totalDof = M_dispFESpace->dof().numTotalDof();
////    VectorElemental dk_loc ( M_dispFESpace->fe().nbFEDof(), this->M_dispFESpace->fieldDim() );
////    VectorElemental elVecDet ( M_dispFESpace->fe().nbFEDof(), this->M_dispFESpace->fieldDim() );
////    
////    //Building fake quadrature rules to compute the deformation gradient F at the nodes
////    QuadratureRule fakeQuadratureRule;
////    
////    Real refElemArea (0); //area of reference element
////    //compute the area of reference element
////    for (UInt iq = 0; iq < M_dispFESpace->qr().nbQuadPt(); iq++)
////    {
////        refElemArea += M_dispFESpace->qr().weight (iq);
////    }
////    
////    Real wQuad (refElemArea / M_dispFESpace->refFE().nbDof() );
////    
////    //Setting the quadrature Points = DOFs of the element and weight = 1
////    std::vector<GeoVector> coords = M_dispFESpace->refFE().refCoor();
////    std::vector<Real> weights (M_dispFESpace->fe().nbFEDof(), wQuad);
////    fakeQuadratureRule.setDimensionShape ( shapeDimension (M_dispFESpace->refFE().shape() ), M_dispFESpace->refFE().shape() );
////    fakeQuadratureRule.setPoints (coords, weights);
////    
////    //Set the new quadrature rule
////    M_dispFESpace->setQuadRule (fakeQuadratureRule);
////    //Loop over the volumes. No markerIDs are necessary on this loop for J = det(F)
////    for ( UInt i = 0; i < M_dispFESpace->mesh()->numVolumes(); ++i )
////    {
////        
////        //Vectors for the deformation tensor
////        std::vector<matrixSerialDense_Type> vectorDeformationF (M_dispFESpace->fe().nbFEDof(), *M_deformationF);
////        M_invariants.resize   ( M_dispFESpace->fieldDim() + 1 );
////        
////        
////        M_dispFESpace->fe().updateFirstDerivQuadPt ( M_dispFESpace->mesh()->volumeList ( i ) );
////        
////        UInt eleID = M_dispFESpace->fe().currentLocalId();
////        
////        //Extracting the local displacement
////        for ( UInt iNode = 0; iNode < ( UInt ) M_dispFESpace->fe().nbFEDof(); iNode++ )
////        {
////            UInt  iloc = M_dispFESpace->fe().patternFirst ( iNode );
////            
////            for ( UInt iComp = 0; iComp < this->M_dispFESpace->fieldDim(); ++iComp )
////            {
////                UInt ig = M_dispFESpace->dof().localToGlobalMap ( eleID, iloc ) + iComp * M_dispFESpace->dim() + this->M_offset;
////                dk_loc[iloc + iComp * M_dispFESpace->fe().nbFEDof()] = dRep[ig];
////            }
////        }
////        
////        //computing the tensor F
////        AssemblyElementalStructure::computeLocalDeformationGradient ( dk_loc, vectorDeformationF, M_dispFESpace->fe() );
////        
////        
////        //Cycle over the nDofs/element
////        for ( UInt nDOF = 0; nDOF < ( UInt ) M_dispFESpace->fe().nbFEDof(); nDOF++ )
////        {
////            UInt  iloc = M_dispFESpace->fe().patternFirst ( nDOF );
////            
////            //computing the determinant
////            AssemblyElementalStructure::computeInvariantsRightCauchyGreenTensor ( M_invariants, vectorDeformationF[nDOF] );
////            
////            //Assembling elemental vector
////            (elVecDet) [ iloc ] = M_invariants[3];
////            (elVecDet) [ iloc + M_dispFESpace->fe().nbFEDof() ] = 0.0;
////            (elVecDet) [ iloc + 2 * M_dispFESpace->fe().nbFEDof() ] = 0.0;
////            
////        }
////        
////        //multiplying it for the patch area
////        reconstructElementaryVector ( elVecDet, patchAreaR, i );
////        
////        //assembling it into the global vector
////        for ( UInt ic = 0; ic < this->M_dispFESpace->fieldDim(); ++ic )
////        {
////            assembleVector (vectorJacobian, elVecDet, M_dispFESpace->fe(), M_dispFESpace->dof(), ic, this->M_offset +  ic * totalDof );
////        }
////        
////        vectorDeformationF.clear();
////        M_invariants.clear();
////        
////    }
////    
////    vectorJacobian.globalAssemble();
////    
////    chrono.stop();
////    M_Displayer->leaderPrintMax ("done in ", chrono.diff() );
////    
////    jacobianDistribution = vectorJacobian;
//}

class HeartSolver
{
public:
void computeKinematicsVariables( const VectorElemental& dk_loc , boost::scoped_ptr<CurrentFE>& M_diffCFE )
{
    
    Real s;
    
    //! loop on quadrature points (ig)
    for ( UInt ig = 0; ig < M_diffCFE->nbQuadPt(); ig++ )
    {
        //! loop on space coordinates (icoor)
        for ( UInt icoor = 0; icoor < nDimensions; icoor++ )
        {
            //! loop  on space coordinates (jcoor)
            for ( UInt jcoor = 0; jcoor < nDimensions; jcoor++ )
            {
                s = 0.0;
                for ( UInt i = 0; i < M_diffCFE->nbFEDof(); i++ )
                {
                    //! \grad u^k at a quadrature point
                    s += M_diffCFE->phiDer( i, jcoor, ig ) * dk_loc[ i + icoor * M_diffCFE->nbFEDof() ];
                }
                //! gradient of displacement
                (*M_Fk)[ icoor ][ jcoor ][ig ] = s;
            }
        }
    }
    
    //! loop on quadrature points (ig)
    for ( UInt ig = 0; ig < M_diffCFE->nbQuadPt(); ig++ )
    {
        //! loop on space coordinates (icoor)
        for ( UInt  icoor = 0;icoor < nDimensions; icoor++ )
        {
            //! deformation gradient Fk
            (*M_Fk)[ icoor ][ icoor ][ ig ] +=  1.0;
        }
    }
    
    Real a,b,c,d,e,f,g,h,i;
    
    for( UInt ig=0; ig< M_diffCFE->nbQuadPt(); ig++ )
    {
        a = (*M_Fk)[ 0 ][ 0 ][ ig ];
        b = (*M_Fk)[ 0 ][ 1 ][ ig ];
        c = (*M_Fk)[ 0 ][ 2 ][ ig ];
        d = (*M_Fk)[ 1 ][ 0 ][ ig ];
        e = (*M_Fk)[ 1 ][ 1 ][ ig ];
        f = (*M_Fk)[ 1 ][ 2 ][ ig ];
        g = (*M_Fk)[ 2 ][ 0 ][ ig ];
        h = (*M_Fk)[ 2 ][ 1 ][ ig ];
        i = (*M_Fk)[ 2 ][ 2 ][ ig ];
        
        //! determinant of deformation gradient Fk
        (*M_Jack)[ig] = a*( e*i - f*h ) - b*( d*i - f*g ) + c*( d*h - e*g );
        
        ASSERT_PRE((*M_Jack)[ig] > 0, "Negative Jacobian. Error!" );
        
        (*M_CofFk)[ 0 ][ 0 ][ ig ] =   ( e*i - f*h );
        (*M_CofFk)[ 0 ][ 1 ][ ig ] = - ( d*i - g*f );
        (*M_CofFk)[ 0 ][ 2 ][ ig ] =   ( d*h - e*g );
        (*M_CofFk)[ 1 ][ 0 ][ ig ] = - ( b*i - c*h );
        (*M_CofFk)[ 1 ][ 1 ][ ig ] =   ( a*i - c*g );
        (*M_CofFk)[ 1 ][ 2 ][ ig ] = - ( a*h - g*b );
        (*M_CofFk)[ 2 ][ 0 ][ ig ] =   ( b*f - c*e );
        (*M_CofFk)[ 2 ][ 1 ][ ig ] = - ( a*f - c*d );
        (*M_CofFk)[ 2 ][ 2 ][ ig ] =   ( a*e - d*b );
    }
    
    //! loop on quadrature points
    for ( UInt ig = 0;ig < M_diffCFE->nbQuadPt(); ig++ )
    {
        s = 0.0;
        for ( UInt i = 0; i < nDimensions; i++)
        {
            for ( UInt j = 0; j < nDimensions; j++)
            {
                //! trace of  C1 = (F1k^t F1k)
                s +=  (*M_Fk)[ i ][ j ][ ig ] * (*M_Fk)[ i ][ j ][ ig ];
            }
        }
        (*M_trCk)[ ig ] = s;
    }
    
    for ( UInt ig = 0; ig <  M_diffCFE->nbQuadPt(); ig++ )
    {
        //! trace of deviatoric C
        (*M_trCisok)[ ig ] =  pow((*M_Jack)[ ig ], -2./3.) * (*M_trCk)[ ig ];
    }
}

    
    
    //! Local tensors initialization
    boost::shared_ptr<boost::multi_array<Real, 3> > M_Fk;
    boost::shared_ptr<boost::multi_array<Real, 3> > M_CofFk;
    
    boost::shared_ptr<std::vector<Real> > M_Jack;
    boost::shared_ptr<std::vector<Real> > M_trCisok;
    boost::shared_ptr<std::vector<Real> > M_trCk;
    
void
computeJacobian (matrix_ptrType jacobian, vectorPtr_Type disp, fespacePtr_Type M_fespace)
{
    typedef CurrentFE                                    currentFE_type;
    typedef boost::scoped_ptr<currentFE_type>            currentFE_ptrType;
    typedef MatrixElemental                              localMatrix_type;
    typedef boost::scoped_ptr<localMatrix_type>          localMatrix_ptrType;

    currentFE_ptrType M_diffCFE;
    localMatrix_ptrType M_localDiff;
    
    M_diffCFE.reset (new currentFE_type (M_fespace->refFE(), M_fespace->fe().geoMap(), M_fespace->qr() ) );
    M_localDiff.reset (new localMatrix_type (M_fespace->fe().nbFEDof(),
                                             M_fespace->fieldDim(),
                                             M_fespace->fieldDim() ) );

    M_Fk.reset ( new boost::multi_array<Real, 3>(boost::extents[3][3][M_fespace->fe().nbQuadPt()]) );
    M_CofFk.reset ( new boost::multi_array<Real, 3>(boost::extents[3][3][M_fespace->fe().nbQuadPt()]) );
    
    M_Jack.reset ( new std::vector<Real>(M_fespace->fe().nbQuadPt(),0.0) );
    M_trCisok.reset ( new std::vector<Real>(M_fespace->fe().nbQuadPt(),0.0) );
    M_trCk.reset ( new std::vector<Real>(M_fespace->fe().nbQuadPt(),0.0) );

    
    const UInt offsetLeft = 0;
    const UInt offsetUp = 0;
    const Real coefficient = 1.0;
    
    // Some constants
    const UInt nbElements (M_fespace->mesh()->numElements() );
    const UInt fieldDim (M_fespace->fieldDim() );
    const UInt nbTotalDof (M_fespace->dof().numTotalDof() );
    
    
    //this->M_jacobian.reset(new matrix_Type(*this->M_localMap));
    
    VectorElemental dk_loc(M_fespace->fe().nbFEDof(), 3);
    
    vector_Type dRep(*disp, Repeated);
    
    //! Number of displacement components
    UInt nc = nDimensions;
    

        Real bulk = 350000;
        Real alpha = 5;
        Real gamma = 5;
    //std::cout << "c" << std::endl;

    // Loop over the elements
    for (UInt iterElement (0); iterElement < nbElements; ++iterElement)
    {
        // Update the diffusion current FE
        M_diffCFE->update ( M_fespace->mesh()->element (iterElement), UPDATE_DPHI | UPDATE_WDET );
        
        //M_fespace->fe().updateFirstDerivQuadPt ( M_dispFESpace->mesh()->volumeList ( i ) );

            UInt eleID = M_diffCFE->currentLocalId();
            
            for ( UInt iNode = 0; iNode < ( UInt ) M_diffCFE->nbFEDof(); iNode++ )
            {
                UInt  iloc = M_diffCFE->patternFirst( iNode );
                
                for ( UInt iComp = 0; iComp < nDimensions; ++iComp )
                {
                    UInt ig = M_fespace->dof().localToGlobalMap( eleID, iloc ) + iComp*M_fespace->dim() + 0;
                    dk_loc[iloc + iComp*M_fespace->fe().nbFEDof()] = dRep[ig];
                }
            }
        //std::cout << iterElement << " / " << nbElements << std::endl;

            M_localDiff->zero();
            
            //! Computes F, Cof(F), J = det(F), Tr(C)
            computeKinematicsVariables( dk_loc, M_diffCFE );

            //! VOLUMETRIC PART
            //! 1. Stiffness matrix: int { 1/2 * bulk * ( 2 - 1/J + 1/J^2 ) * ( CofF : \nabla \delta ) (CofF : \nabla v) }
            AssemblyElementalStructure::stiff_Jac_Pvol_1term( 0.5 *  bulk, (*M_CofFk), (*M_Jack), *M_localDiff, *M_diffCFE );
            
            //! 2. Stiffness matrix: int { 1/2 * bulk * ( 1/J- 1 - log(J)/J^2 ) * ( CofF1 [\nabla \delta]^t CofF ) : \nabla v }
            AssemblyElementalStructure::stiff_Jac_Pvol_2term( 0.5 *  bulk, (*M_CofFk), (*M_Jack), *M_localDiff, *M_diffCFE );
            
            //! ISOCHORIC PART
            //! 1. Stiffness matrix : int { - 2/3 * alpha * J^(-5/3) * exp( gamma*( Ic_iso - 3) )* ( 1. + coefExp * Ic_iso )
            //!                      *( CofF : \nabla \delta ) ( F : \nabla \v ) }
            AssemblyElementalStructure::stiff_Jac_P1iso_Exp_1term( (-2.0/3.0) * alpha, gamma, (*M_CofFk), (*M_Fk), (*M_Jack), (*M_trCisok), *M_localDiff, *M_diffCFE );
            
            //! 2. Stiffness matrix : int { 2 * alpha * gamma * J^(-4/3) * exp( gamma*( Ic_iso - 3) ) *
            //!             ( F : \nabla \delta ) ( F : \nabla \v ) }
            AssemblyElementalStructure::stiff_Jac_P1iso_Exp_2term( 2.0 * alpha * gamma, gamma, (*M_Fk), (*M_Jack), (*M_trCisok), *M_localDiff, *M_diffCFE );
            
            //! 3. Stiffness matrix : int { 2.0/9.0 *  alpha * J^-2 * Ic_iso * exp( gamma*( Ic_iso - 3) )*
            //!                       ( 1. + gamma * Ic_iso )( CofF : \nabla \delta ) ( CofF : \nabla \v )}
            AssemblyElementalStructure::stiff_Jac_P1iso_Exp_3term( (2.0/9.0) * alpha, gamma, (*M_CofFk), (*M_Jack), (*M_trCisok), *M_localDiff, *M_diffCFE );
            
            //! 4. Stiffness matrix : int { -2.0/3.0 *  alpha * J^(-5/3) * exp( gamma*( Ic_iso - 3) )
            //!                    * ( 1. + gamma * Ic_iso )( F : \nabla \delta ) ( CofF : \nabla \v ) }
            AssemblyElementalStructure::stiff_Jac_P1iso_Exp_4term( (-2.0/3.0) * alpha, gamma, (*M_CofFk), (*M_Fk), (*M_Jack), (*M_trCisok), *M_localDiff, *M_diffCFE );
            
            //! 5. Stiffness matrix : int {  alpha * J^(-2/3) * exp( gamma*( Ic_iso - 3)) (\nabla \delta: \nabla \v)}
            AssemblyElementalStructure::stiff_Jac_P1iso_Exp_5term( alpha, gamma, (*M_Jack), (*M_trCisok), *M_localDiff, *M_diffCFE );
            
            //! 6. Stiffness matrix : int { 1.0/3.0 * alpha * J^(-2) * Ic_iso *  exp(gamma*( Ic_iso - 3)) *
            //!                       (CofF [\nabla \delta]^t CofF ) : \nabla \v  }
            AssemblyElementalStructure::stiff_Jac_P1iso_Exp_6term( (1.0/3.0) * alpha, gamma, (*M_CofFk), (*M_Jack), (*M_trCisok), *M_localDiff, *M_diffCFE );
            
            //! assembling
            for ( UInt ic = 0; ic < nc; ++ic )
            {
                for ( UInt jc = 0; jc < nc; jc++ )
                {
                    assembleMatrix( *jacobian,
                                    *M_localDiff,
                                    *M_diffCFE,
                                    M_fespace->dof(),
                                    ic,
                                    jc,
                                    0 +  ic*nbTotalDof,
                                    0 +  jc*nbTotalDof );
                }
            }
        }
    }
    
    




void computeStiffness( vectorPtr_Type M_stiff, vectorPtr_Type sol , fespacePtr_Type M_fespace )
{
    typedef CurrentFE                                    currentFE_type;
    typedef boost::scoped_ptr<currentFE_type>            currentFE_ptrType;
    typedef VectorElemental                              localVector_type;
    typedef boost::scoped_ptr<localVector_type>          localVector_ptrType;
    
    currentFE_ptrType M_diffCFE;
    localVector_ptrType M_elvecK;
    
    M_diffCFE.reset (new currentFE_type (M_fespace->refFE(), M_fespace->fe().geoMap(), M_fespace->qr() ) );
    M_elvecK.reset (new localVector_type (M_fespace->fe().nbFEDof(),
                                             M_fespace->fieldDim() ) );


    //std::cout << "c" << std::endl;

    
    VectorElemental dk_loc(M_fespace->fe().nbFEDof(), 3);
    
    vector_Type dRep(*sol, Repeated);
    //std::cout << "c" << std::endl;

    // Some constants
    const UInt nbElements (M_fespace->mesh()->numElements() );
    const UInt fieldDim (M_fespace->fieldDim() );
    const UInt nbTotalDof (M_fespace->dof().numTotalDof() );
    //std::cout << "c" << std::endl;

    // Loop over the elements
    for (UInt iterElement (0); iterElement < nbElements; ++iterElement)
    {
        
        Real bulk =350000;
        Real alpha = 5;
        Real gamma = 5;
        
        std::cout << "d" << std::endl;

            M_diffCFE->update ( M_fespace->mesh()->element (iterElement), UPDATE_DPHI | UPDATE_WDET );
            
        std::cout << "d" << std::endl;

            UInt eleID = M_diffCFE->currentLocalId();
            
            for ( UInt iNode = 0; iNode < ( UInt ) M_diffCFE->nbFEDof(); iNode++ )
            {
                UInt  iloc = M_diffCFE->patternFirst( iNode );
                
                for ( UInt iComp = 0; iComp < nDimensions; ++iComp )
                {
                    UInt ig = M_fespace->dof().localToGlobalMap( eleID, iloc ) + iComp*M_fespace->dim() + 0;
                    dk_loc[iloc + iComp*M_fespace->fe().nbFEDof()] = dRep[ig];
                }
            }

        std::cout << iterElement << " / " << nbElements << std::endl;

            M_elvecK->zero();
            
            computeKinematicsVariables( dk_loc , M_diffCFE );
            
            //! Stiffness for non-linear terms of the Neo-Hookean model
            /*!
             The results of the integrals are stored at each step into elvecK, until to build K matrix of the bilinear form
             */
            //! Volumetric part
            /*!
             Source term Pvol: int { bulk /2* (J1^2 - J1  + log(J1) ) * 1/J1 * (CofF1 : \nabla v) }
             */
            AssemblyElementalStructure::source_Pvol( 0.5 * bulk, (*M_CofFk), (*M_Jack), *M_elvecK, M_fespace->fe() );
            
            //! Isochoric part
            /*!
             Source term P1iso_Exp: int { alpha * exp(gamma *(  Ic1_iso -3 )) *
             ( J1^(-2/3)* (F1 : \nabla v) - 1/3 * (Ic1_iso / J1) * (CofF1 : \nabla v) ) }
             */
            AssemblyElementalStructure::source_P1iso_Exp( alpha, gamma, (*M_CofFk), (*M_Fk), (*M_Jack), (*M_trCisok), *M_elvecK, M_fespace->fe() );
            
            for ( UInt ic = 0; ic < nDimensions; ++ic )
            {
                /*!
                 M_elvecK is assemble into *vec_stiff vector that is recall
                 from updateSystem(matrix_ptrtype& mat_stiff, vector_ptr_type& vec_stiff)
                 */
                assembleVector( *M_stiff,
                                *M_elvecK,
                                M_fespace->fe(),
                                M_fespace->dof(),
                                ic,
                                0 +  ic*nbTotalDof );
            }
        
    }
    
    M_stiff->globalAssemble();
}
    
    
};





void
addDiffusion (matrix_ptrType matrix, fespacePtr_Type M_fespace)
{
    typedef CurrentFE                                    currentFE_type;
    typedef boost::scoped_ptr<currentFE_type>            currentFE_ptrType;
    typedef MatrixElemental                              localMatrix_type;
    typedef boost::scoped_ptr<localMatrix_type>          localMatrix_ptrType;

    currentFE_ptrType M_diffCFE;
    localMatrix_ptrType M_localDiff;
    
    M_diffCFE.reset (new currentFE_type (M_fespace->refFE(), M_fespace->fe().geoMap(), M_fespace->qr() ) );
    M_localDiff.reset (new localMatrix_type (M_fespace->fe().nbFEDof(),
                                             M_fespace->fieldDim(),
                                             M_fespace->fieldDim() ) );
    
    const UInt offsetLeft = 0;
    const UInt offsetUp = 0;
    const Real coefficient = 1.0;

//    // Check that the fespace is set
//    ASSERT (M_fespace != 0, "No FE space for assembling the diffusion!");
//    
//    M_diffusionAssemblyChrono.start();
    
    // Some constants
    const UInt nbElements (M_fespace->mesh()->numElements() );
    const UInt fieldDim (M_fespace->fieldDim() );
    const UInt nbTotalDof (M_fespace->dof().numTotalDof() );
    
    // Loop over the elements
    for (UInt iterElement (0); iterElement < nbElements; ++iterElement)
    {
        // Update the diffusion current FE
        M_diffCFE->update ( M_fespace->mesh()->element (iterElement), UPDATE_DPHI | UPDATE_WDET );
        
        // Clean the local matrix
        M_localDiff->zero();
        
        // local stiffness
        AssemblyElemental::stiffness (*M_localDiff, *M_diffCFE, coefficient, fieldDim);
        
        // Assembly
        for (UInt iFieldDim (0); iFieldDim < fieldDim; ++iFieldDim)
        {
            assembleMatrix ( *matrix,
                            *M_localDiff,
                            *M_diffCFE,
                            *M_diffCFE,
                            M_fespace->dof(),
                            M_fespace->dof(),
                            iFieldDim, iFieldDim,
                            iFieldDim * nbTotalDof + offsetLeft, iFieldDim * nbTotalDof + offsetUp );
        }
    }
    
//    M_diffusionAssemblyChrono.stop();
}


static Real bcZero (const Real& /*t*/, const Real&  X, const Real& /*Y*/, const Real& /*Z*/, const ID& i)
{
//    switch (i)
//    {
//        case 0:
//            return  0.0;
//            break;
//        case 1:
//            return  300000.0;
//            break;
//        case 2:
//            return 0.0;
//            break;
//            
//    }
    return 0;//(0.01*X*X);
}

static Real bcForce (const Real& /*t*/, const Real&  X, const Real& Y, const Real& /*Z*/, const ID& i)
{
    if ( Y > 0.5 ) return 0.001;
    return 0;
}



int main ( int argc, char** argv )
{
    // +-----------------------------------------------+
    // |            Initialization of MPI              |
    // +-----------------------------------------------+
#ifdef HAVE_MPI
    MPI_Init (&argc, &argv);
#endif
    bool verbose ( true );

    {

#ifdef HAVE_MPI
        boost::shared_ptr<Epetra_Comm> Comm ( new Epetra_MpiComm ( MPI_COMM_WORLD ) );
#else
        boost::shared_ptr<Epetra_Comm> Comm ( new Epetra_SerialComm );
#endif

        verbose = (Comm->MyPID() == 0);

        // +-----------------------------------------------+
        // |               Loading the data                |
        // +-----------------------------------------------+
        if ( verbose )
        {
            std::cout << std::endl << "[Loading the data]" << std::endl;
        }
        LifeChrono globalChrono;

        globalChrono.start();

        // ********** GetPot **********
        GetPot command_line ( argc, argv );
        const std::string dataFileName = command_line.follow ( "data", 2, "-f", "--file" );
        GetPot dataFile ( dataFileName );
        // ****************************

        // Space discretization
        const UInt numMeshElem    = dataFile ( "mesh/num_elements", 10);
        const std::string uOrder  = dataFile ( "finite_element/velocity", "P2" );

        // +-----------------------------------------------+
        // |               Loading the mesh                |
        // +-----------------------------------------------+
        meshPtr_Type fullMeshPtr ( new mesh_Type );
        regularMesh3D ( *fullMeshPtr,
                        1,
                        numMeshElem, numMeshElem, numMeshElem,
                        false,
                        1.0, 1.0, 1.0,
                        0.0, 0.0, 0.0 );

        meshPtr_Type meshPtr;
        MeshPartitioner< mesh_Type >   meshPart ( fullMeshPtr, Comm );
        meshPtr = meshPart.meshPartition();
        
        fullMeshPtr.reset(); //Freeing the global mesh to save memory

        // +-----------------------------------------------+
        // |            Creating the FE spaces             |
        // +-----------------------------------------------+
        fespacePtr_Type uFESpace ( new fespace_Type ( meshPtr, "P2", 3, Comm ) );

        // +-----------------------------------------------+
        // |             Boundary conditions               |
        // +-----------------------------------------------+
        BCHandler bcHandler;
        BCFunctionBase bcZeroFunc ( bcZero );
        BCFunctionBase bcForceFunc ( bcForce );

        //BCFunctionBase fRHS ( Laplacian::f );

        //for ( UInt iDirichlet ( 1 ); iDirichlet <= 26; ++iDirichlet )
        {
            bcHandler.addBC ( "Wall", 1, Essential, Full, bcZeroFunc, 3 );
        }
        
        //for ( UInt iDirichlet ( 1 ); iDirichlet <= 1000; ++iDirichlet )
        {
            bcHandler.addBC ( "Force", 3, Essential, Component, bcForceFunc , 1 );
        }


        // Update the BCHandler (internal data related to FE)
        bcHandler.bcUpdate ( *uFESpace->mesh(), uFESpace->feBd(), uFESpace->dof() );

        // +-----------------------------------------------+
        // |              Matrices Assembly                |
        // +-----------------------------------------------+
        ADRAssembler<mesh_Type, matrix_Type, vector_Type> adrAssembler;
        adrAssembler.setup ( uFESpace, uFESpace );
        
        boost::shared_ptr<matrix_Type> systemMatrix ( new matrix_Type ( uFESpace->map() ) );
        
        //ADRAssembler<mesh_Type, matrix_Type, vector_Type> adrAssembler;
        //adrAssembler.setup ( uFESpace, uFESpace );
        //adrAssembler.addDiffusion ( systemMatrix, 1.0 );

        //addDiffusion(systemMatrix, uFESpace);
        
        boost::shared_ptr<vector_Type> solution;
        solution.reset ( new vector_Type ( uFESpace->map(), Unique ) );
        
        std::cout << "set up heartsolver" << std::endl;
        HeartSolver heartSolver;
        heartSolver.computeJacobian(systemMatrix, solution, uFESpace);
        std::cout << "set up heartsolver done" << std::endl;

        // +-----------------------------------------------+
        // |            Solver initialization              |
        // +-----------------------------------------------+
        if ( verbose )
        {
            std::cout << std::endl << "[Solvers initialization]" << std::endl;
        }
        prec_Type* precRawPtr;
        basePrecPtr_Type precPtr;
        precRawPtr = new prec_Type;
        precRawPtr->setDataFromGetPot ( dataFile, "prec/LinearSolver" );
        precPtr.reset ( precRawPtr );

        // Linear solver
        Teuchos::RCP< Teuchos::ParameterList > belosList = Teuchos::rcp ( new Teuchos::ParameterList );
        belosList = Teuchos::getParametersFromXmlFile ( "SolverParamList.xml" );

        LinearSolver linearSolver;
        linearSolver.setCommunicator ( Comm );
        linearSolver.setParameters ( *belosList );
        linearSolver.setPreconditioner ( precPtr );

        
        linearSolver.showMe();

        // +-----------------------------------------------+
        // |                   Simulation                  |
        // +-----------------------------------------------+
        if ( verbose )
        {
            std::cout << std::endl << "[Initialization of the simulation]" << std::endl;
        }
        if ( verbose )
        {
            std::cout << "Creation of vectors... " << std::flush;
        }

        boost::shared_ptr<vector_Type> rhs;
        rhs.reset ( new vector_Type ( uFESpace->map(), Repeated ) );

        //adrAssembler.addMassRhs ( *rhs, fRHS, 0.0 );
        
        heartSolver.computeStiffness(rhs, solution, uFESpace);

        rhs->globalAssemble();
        
        if ( verbose )
        {
            std::cout << "done" << std::endl;
        }


        if ( verbose )
        {
            std::cout << "Applying BC... " << std::flush;
        }
        systemMatrix->globalAssemble();
        boost::shared_ptr<vector_Type> rhsBC;
        rhsBC.reset ( new vector_Type ( *rhs, Unique ) );
        bcManage ( *systemMatrix, *rhsBC, *uFESpace->mesh(), uFESpace->dof(), bcHandler, uFESpace->feBd(), 1.0, 0.0 );
        if ( verbose )
        {
            std::cout << "done" << std::endl;
        }

        if ( verbose )
        {
            std::cout << std::endl << "Solving the system with LinearSolver (Belos)... " << std::endl;
        }
//        boost::shared_ptr<vector_Type> solution;
//        solution.reset ( new vector_Type ( uFESpace->map(), Unique ) );
        linearSolver.setOperator ( systemMatrix );
        linearSolver.setRightHandSide ( rhsBC );
        linearSolver.solve ( solution );



        // +-----------------------------------------------+
        // |            Ending the simulation              |
        // +-----------------------------------------------+
        globalChrono.stop();
        if ( verbose )
        {
            std::cout << std::endl << "Total simulation time: " << globalChrono.diff() << " s." << std::endl;
        }


        typedef LifeV::RegionMesh<LinearTetra>                              mesh_Type;
        typedef LifeV::ExporterHDF5<mesh_Type >                             hdf5Filter_Type;
        hdf5Filter_Type exporter ( dataFile, "jacDist" );
        exporter.setMeshProcId (uFESpace->mesh(), uFESpace->map().comm().MyPID() );
        exporter.addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, "displacementField", uFESpace, solution, UInt (0) );
        exporter.postProcess ( 0.0 );

        
    }
    
    
    
#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    return 0;
}
