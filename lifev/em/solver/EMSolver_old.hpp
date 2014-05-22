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
 @brief Class for solving the Monodomain equations in electrophysiology.

 @date 02-2013
 @author Simone Rossi <simone.rossi@epfl.ch>

 @last update 02-2013

 This class provides interfaces to solve the monodomain equation
 ( reaction diffusion equation ) using the ETA framework.
 The solution can be performed using three different methods:
 -operator splitting method (at this point available only with forward Euler
 for the reaction step and backward Euler for the diffusion step. );
 -Ionic Currents Interpolation (at this point only forward Euler);
 -State Variable interpolation (at this point only forward Euler).
 */

#ifndef _EMSOLVERf_H_
#define _EMSOLVERf_H_


#include <lifev/em/solver/electrophysiology/EMMonodomainSolver.hpp>
//
//#include <lifev/em/solver/EMStructuralOperator.hpp>
//#include <lifev/em/solver/EMGeneralizedActiveHolzapfelOgdenMaterial.hpp>
//#include <lifev/em/solver/EMActiveStrainSolver.hpp>
//#include <lifev/core/interpolation/RBFlocallyRescaledVectorial.hpp>
//#include <lifev/core/interpolation/RBFlocallyRescaledScalar.hpp>
//#include <lifev/core/interpolation/RBFrescaledVectorial.hpp>
//#include <lifev/core/interpolation/RBFrescaledScalar.hpp>
////#include <lifev/core/interpolation/RBFscalar.hpp>
//#include <lifev/core/interpolation/RBFvectorial.hpp>
//
//#include <lifev/bc_interface/3D/bc/BCInterface3D.hpp>
//
//
//#include <lifev/em/solver/EMEvaluate.hpp>

namespace LifeV
{

//! EMSolver - Class featuring the solution of the electromechanical problem with monodomain equation

template<typename Mesh , typename IonicModel>
class EMSolver
{

//    //!Monodomain Solver
//    /*!
//     The monodomain equation reads
//     \f \Chi
//
//     */
//
//public:
//
//    //! @name Type definitions
//    //@{
//
//    typedef Mesh mesh_Type;
//    typedef boost::shared_ptr<mesh_Type> meshPtr_Type;
//
//    typedef VectorEpetra vector_Type;
//    typedef boost::shared_ptr<VectorEpetra> vectorPtr_Type;
//
//    typedef std::vector<vectorPtr_Type> vectorOfPtr_Type;
//
//    typedef MatrixEpetra<Real> matrix_Type;
//    typedef boost::shared_ptr<matrix_Type> matrixPtr_Type;
//
//    typedef Epetra_Comm comm_Type;
//    typedef boost::shared_ptr<comm_Type> commPtr_Type;
//
//    typedef ETFESpace<mesh_Type, MapEpetra, 3, 1> ETFESpace_Type;
//    typedef boost::shared_ptr<ETFESpace<mesh_Type, MapEpetra, 3, 1> > ETFESpacePtr_Type;
//
//    typedef ETFESpace<mesh_Type, MapEpetra, 3, 3> ETFESpaceVectorial_Type;
//    typedef boost::shared_ptr< ETFESpaceVectorial_Type > ETFESpaceVectorialPtr_Type;
//
//    typedef FESpace<mesh_Type, MapEpetra> feSpace_Type;
//    typedef boost::shared_ptr<feSpace_Type> feSpacePtr_Type;
//
//    typedef LinearSolver linearSolver_Type;
//    typedef boost::shared_ptr<LinearSolver> linearSolverPtr_Type;
//
//    typedef ExporterHDF5< mesh_Type >          exporter_Type;
//    typedef boost::shared_ptr<exporter_Type>                       exporterPtr_Type;
//    //  typedef Exporter<mesh_Type> exporter_Type;    //                IOFile_Type;
//    //  typedef boost::shared_ptr<exporter_Type> exporterPtr_Type; //                IOFilePtr_Type;
//
//    typedef LifeV::Preconditioner basePrec_Type;
//    typedef boost::shared_ptr<basePrec_Type> basePrecPtr_Type;
//    typedef LifeV::PreconditionerIfpack prec_Type;
//    typedef boost::shared_ptr<prec_Type> precPtr_Type;
//
//    typedef IonicModel ionicModel_Type;
//    typedef ElectroIonicModel superIonicModel;
//    typedef boost::shared_ptr<ionicModel_Type> ionicModelPtr_Type;
//
//    typedef Teuchos::ParameterList list_Type;
//
//    typedef boost::function <
//    Real (const Real& t, const Real& x, const Real& y, const Real& z,
//          const ID& i) > function_Type;
//
//    typedef MatrixSmall<3, 3>                          matrixSmall_Type;
//
//    typedef EMMonodomainSolver<mesh_Type, ionicModel_Type>                  monodomainSolver_Type;
//    typedef boost::shared_ptr<monodomainSolver_Type>                    monodomainSolverPtr_Type;
//
//    typedef FESpace< RegionMesh<LinearTetra>, MapEpetra >               FESpace_Type;
//    typedef boost::shared_ptr<FESpace_Type>                        FESpacePtr_Type;
//
//    typedef ETFESpace< RegionMesh<LinearTetra>, MapEpetra, 3, 1 >       scalarETFESpace_Type;
//    typedef boost::shared_ptr<scalarETFESpace_Type>                      scalarETFESpacePtr_Type;
//    typedef ETFESpace< RegionMesh<LinearTetra>, MapEpetra, 3, 3 >       solidETFESpace_Type;
//    typedef boost::shared_ptr<solidETFESpace_Type>                      solidETFESpacePtr_Type;
//
//    typedef StructuralConstitutiveLawData           structureData_Type;
//    typedef boost::shared_ptr<structureData_Type>           structureDataPtr_Type;
//
//    typedef EMStructuralOperator< RegionMesh<LinearTetra> > structuralOperator_Type;
//    typedef boost::shared_ptr< structuralOperator_Type > structuralOperatorPtr_Type;
//
//
//    typedef BCHandler                                          bc_Type;
//    typedef boost::shared_ptr< bc_Type >                       bcPtr_Type;
//    typedef StructuralOperator< RegionMesh<LinearTetra> >       physicalSolver_Type;
//    typedef BCInterface3D< bc_Type, physicalSolver_Type >      bcInterface_Type;
//    typedef boost::shared_ptr< bcInterface_Type >              bcInterfacePtr_Type;
//
//    typedef EMActiveStrainSolver<mesh_Type> activeStrain_Type;
//    typedef boost::shared_ptr< activeStrain_Type >  activeStrainPtr_Type;
//
//    typedef RBFInterpolation<mesh_Type>           interpolation_Type;
//    typedef boost::shared_ptr<interpolation_Type> interpolationPtr_Type;
//
//    typedef BCVector                                bcVector_Type;
//    typedef boost::shared_ptr<bcVector_Type>        bcVectorPtr_Type;
//    typedef std::vector<bcVectorPtr_Type>           bcVectorPtrs_Type;
//
//    typedef MapEpetra                               map_Type;
//    typedef boost::shared_ptr<map_Type>             mapPtr_Type;
//    ///////////////////////////////////////////////////////////////////////////
//
//    inline monodomainSolverPtr_Type    monodomainPtr()
//    {
//        return  M_monodomainPtr;
//    }
//    inline bool                     usingDifferentMeshes()
//    {
//        return M_usingDifferentMeshes;
//    }
//    inline structureDataPtr_Type       solidDataPtr()
//    {
//        return M_solidDataPtr;
//    }
//    inline structuralOperatorPtr_Type  solidPtr()
//    {
//        return M_solidPtr;
//    }
//    inline bcInterfacePtr_Type         solidBCPtr()
//    {
//        return M_solidBCPtr;
//    }
//    inline activeStrainPtr_Type         activationPtr()
//    {
//        return M_activationPtr;
//    }
//    inline Real                         monodomainTimeStep()
//    {
//        return M_monodomainTimeStep;
//    }
//    inline Real                         solidTimeStep()
//    {
//        return M_solidTimeStep;
//    }
//    inline Real                         lvPressure()
//    {
//        return M_lvPressure;
//    }
//    inline Real                         rvPressure()
//    {
//        return M_rvPressure;
//    }
//    inline Real                         lvVolume()
//    {
//        return M_lvVolume;
//    }
//    inline Real                         rvVolume()
//    {
//        return M_rvVolume;
//    }
//    inline std::vector<UInt>            lvFlags()
//    {
//        return M_lvFlags;
//    }
//    inline std::vector<UInt>            rvFlags()
//    {
//        return M_rvFlags;
//    }
//    inline bcVectorPtrs_Type            bcVectorPtrs()
//    {
//        return M_bcVectorPtrs;
//    }
//    inline bool                     oneWayCoupling()
//    {
//        return M_oneWayCoupling;
//    }
//    inline vectorPtr_Type               activationSolidPtr()
//    {
//        return M_activationSolidPtr;
//    }
//    inline FESpacePtr_Type              activationSolidFESpacePtr()
//    {
//        return M_solidActivationFESpacePtr;
//    }
//    inline FESpace_Type&                activationSolidFESpace()
//    {
//        return *M_solidActivationFESpacePtr;
//    }
//    inline meshPtr_Type                 fullSolidMesh()
//    {
//        return M_fullSolidMesh;
//    }
//
//    inline void setMonodomainPtr (monodomainSolverPtr_Type p)
//    {
//        M_monodomainPtr = p;
//    }
//    inline void setMonodomainPtr (monodomainSolver_Type& p)
//    {
//        *M_monodomainPtr = p;
//    }
//    inline void setUsingDifferentMeshes (bool p)
//    {
//        M_usingDifferentMeshes = p;
//    }
//    inline void setSolidDataPtr (structureDataPtr_Type p )
//    {
//        M_solidDataPtr = p;
//    }
//    inline void setSolidDataPtr (structureData_Type& p )
//    {
//        *M_solidDataPtr = p;
//    }
//    inline void setSolidPtr (structuralOperatorPtr_Type p)
//    {
//        M_solidPtr = p;
//    }
//    inline void setSolidPtr (structuralOperator_Type& p)
//    {
//        *M_solidPtr = p;
//    }
//    inline void setSolidBCPtr (bcInterfacePtr_Type p)
//    {
//        M_solidBCPtr = p;
//    }
//    inline void setSolidBCPtr (bcInterface_Type& p)
//    {
//        *M_solidBCPtr = p;
//    }
//    inline void setActivationPtr (activeStrainPtr_Type p)
//    {
//        M_activationPtr = p;
//    }
//    inline void setActivationPtr (activeStrain_Type& p)
//    {
//        *M_activationPtr = p;
//    }
//    inline void setMonodomainTimeStep (Real p)
//    {
//        M_monodomainTimeStep = p;
//    }
//    inline void setSolidTimeStep (Real p)
//    {
//        M_solidTimeStep = p;
//    }
//    inline void setLVPressure (Real p)
//    {
//        M_lvPressure = p;
//    }
//    inline void setRVPressure (Real p)
//    {
//        M_rvPressure = p;
//    }
//    inline void setLVVolume (Real p)
//    {
//        M_lvVolume = p;
//    }
//    inline void setRVVolume (Real p)
//    {
//        M_rvVolume = p;
//    }
//    inline void setLVFlags (std::vector<UInt> p)
//    {
//        M_lvFlags = p;
//    }
//    inline void setRVFlags (std::vector<UInt> p)
//    {
//        M_rvFlags = p;
//    }
//    inline void setBCVectorPtrs (bcVectorPtrs_Type p)
//    {
//        M_bcVectorPtrs = p;
//    }
//    inline void setOneWayCoupling ( bool p)
//    {
//        M_oneWayCoupling = p;
//    }
//
//
//
//    //@}
//
//    //! @name Constructors & Destructor
//    //@{
//    EMSolver( );
//
//    EMSolver (  Teuchos::ParameterList& parameterList,
//                const std::string data_file_name,
//                commPtr_Type comm );    //!Empty Constructor
//
//    EMSolver (  structuralOperatorPtr_Type solidPtr,
//                Teuchos::ParameterList& parameterList,
//                const std::string data_file_name,
//                commPtr_Type comm );
//
//    virtual ~EMSolver() {};
//
//    void setup (Teuchos::ParameterList& parameterList,
//                const std::string data_file_name,
//                std::string parameterListName = "ParamList.xml");
//
//    void setupMonodomainMatrix (Teuchos::ParameterList& parameterList);
//
//    void setupProjectionSolver (GetPot& dataFile, commPtr_Type comm);
//
//    void setupProjectionRhs();
//
//    void projection();
//
//    void updateMonodomainMatrix();
//
//    void updateMonodomain();
//
//    void updateSolid();
//
//    void solveOneMonodomainStep (int subiter = 1);
//
//    void solveOneReactionStepFE (int subiter = 1);
//
//    void solveOneICIStep();
//
//    void solveOneSVIStep();
//
//    void solveOneHLStep();
//
//    void solveOneOSStep (int subiter = 1);
//
//    void solveOneDiffusionStep();
//
//    void solveOneActivationStep();
//
//    inline void solveSolid()
//    {
//        M_solidPtr -> iterate ( M_solidBCPtr -> handler() );
//    }
//
//    void setupExporters (std::string dir = "./");
//
//    void exportSolution (Real time = 0.0);
//
//    void closeExporters();
//
//    void setupPreloadBC ( GetPot& dataFile, map_Type map );
//
//    void preloadRamp ( Real dt = 0.1);
//
//    void createReferencePositionVector();
//
//    void setupBC (const std::string data_file_name, map_Type map  );
//
//    void computeLVVolume (Real nx, Real ny, Real nz);
//    void computeRVVolume (Real nx, Real ny, Real nz);
//
//    void exportSolidFibersDirection (std::string dir = "./" );
//    void exportMonodomainFibersDirection (std::string dir = "./");
//    void exportSolidSheetsDirection (std::string dir = "./" );
//    void exportFibersAndSheetsFields (std::string dir = "./" );
//    void exportActivationTime (std::string dir = "./" );
//
//    inline void importSolidFibers (Teuchos::ParameterList& parameterList);
//    inline void importSolidSheets (Teuchos::ParameterList& parameterList);
//    inline void importMonodomainFibers (Teuchos::ParameterList& parameterList);
//
//    void setFibersAndSheets (Teuchos::ParameterList& parameterList);
//
//    void setupInterpolants (std::string parameterListName, Teuchos::ParameterList& parameterList, GetPot& dataFile);
//
//
//    inline void createLVPositionVector (Real nx, Real ny, Real nz)
//    {
//        M_lvPositionVectorPtr.reset (new vector_Type ( M_solidPtr -> activeMaterial() -> activationSpace() -> map() ) );
//        createPositionVector (nx, ny, nz, M_lvFlags, *M_lvPositionVectorPtr);
//    }
//    inline void createRVPositionVector (Real nx, Real ny, Real nz)
//    {
//        M_rvPositionVectorPtr.reset (new vector_Type ( M_solidPtr -> activeMaterial() -> activationSpace() -> map() ) );
//        createPositionVector (nx, ny, nz, M_rvFlags, *M_rvPositionVectorPtr );
//    }
//
//    void createPositionVector (Real nx, Real ny, Real nz, std::vector<UInt> flags, vector_Type& vec);
//
//    inline void registerActivationTime ( Real time, Real threshold = 0.0)
//    {
//        M_monodomainPtr -> registerActivationTime (*M_activationTimePtr, time, threshold);
//    }
//
//    inline void setSolidFibers (vector_Type& fibers)
//    {
//        M_solidPtr -> activeMaterial() -> setFiberVector (fibers);
//    }
//    inline void setSolidFibers (vectorPtr_Type fibers)
//    {
//        M_solidPtr -> activeMaterial() -> setFiberVector (*fibers);
//    }
//    inline void setSolidSheets (vector_Type& sheets)
//    {
//        M_solidPtr -> activeMaterial() -> setSheetVector (sheets);
//    }
//    inline void setSolidSheets (vectorPtr_Type sheets)
//    {
//        M_solidPtr -> activeMaterial() -> setSheetVector (*sheets);
//    }
//    //monodomain -> setFiberPtr( electroFibers )
//    inline void setMonodomainFibers (vectorPtr_Type fibers)
//    {
//        M_monodomainPtr -> setFiberPtr (fibers);
//    }
//    inline void setMonodomainFibers (vector_Type& fibers)
//    {
//        M_monodomainPtr -> setFiber (fibers);
//    }
//
//
//    //  void update();
//
//    inline void  setPotentialOnBoundary (Real value, UInt flag)
//    {
//        ElectrophysiologyUtility::setValueOnBoundary ( * (M_monodomainPtr -> potentialPtr() ), M_monodomainPtr -> fullMeshPtr(), value, flag);
//    }
//    inline void  setPotentialFromFunction (function_Type f, Real time = 0.0)
//    {
//        M_monodomainPtr -> setPotentialFromFunction ( f, time );
//    }
//
//    inline Real isochoricPressure (Real Cp, Real pn, Real dV)
//    {
//        return (pn + dV / Cp);
//    }
//    inline void isochoricLVPressure (Real Cp)
//    {
//        M_lvPressure = isochoricPressure (Cp, M_lvPressure, M_lvdV);
//    }
//    inline void isochoricRVPressure (Real Cp)
//    {
//        M_rvPressure = isochoricPressure (Cp, M_rvPressure, M_rvdV);
//    }
//    inline Real windkesselPressure (Real R, Real C, Real pn, Real dV, Real dt = 1.0)
//    {
//        return  pn - ( pn * dt / C / R + dV / C );
//    }
//
//    inline  void setSolutionMethod ( std::string p)
//    {
//        std::map< std::string, solutionMethod_Type > solutionMap;
//        solutionMap["OS"] = OS;
//        solutionMap["ICI"] = ICI;
//        solutionMap["HL"] = HL;
//        solutionMap["SVI"] = SVI;
//        M_solutionMethod   = solutionMap[p];
//    }
//    enum solutionMethod_Type { OS, ICI, HL, SVI };
//    solutionMethod_Type                 M_solutionMethod;
//    /*!
//    */
//    monodomainSolverPtr_Type    M_monodomainPtr;
//    bool                        M_usingDifferentMeshes;
//    structureDataPtr_Type       M_solidDataPtr;
//    structuralOperatorPtr_Type  M_solidPtr;
//    bcInterfacePtr_Type         M_solidBCPtr;
//    activeStrainPtr_Type        M_activationPtr;
//    Real                        M_monodomainTimeStep;
//    Real                        M_solidTimeStep;
//    exporterPtr_Type            M_monodomainExporterPtr;
//    exporterPtr_Type            M_solidExporterPtr;
//    exporterPtr_Type            M_activationExporterPtr;
//
//    //Coarse To Fine ( C2F )
//    vectorPtr_Type              M_monodomainDisplacementPtr;
//    interpolationPtr_Type       M_C2FPtr;
//    //Fine To Coarse ( F2C )
//    vectorPtr_Type              M_activationSolidPtr;
//    interpolationPtr_Type       M_F2CPtr;
//    meshPtr_Type                M_fullSolidMesh;
//
//    vectorPtr_Type              M_activationTimePtr;
//
//    Real                        M_lvPressure;
//    Real                        M_rvPressure;
//    Real                        M_lvVolume;
//    Real                        M_rvVolume;
//    Real                        M_lvdV;
//    Real                        M_rvdV;
//    std::ofstream               M_lvPVexporter;
//    std::ofstream               M_rvPVexporter;
//    bool                        M_lvPV;
//    bool                        M_rvPV;
//
//    std::vector<UInt>           M_lvFlags;
//    std::vector<UInt>           M_rvFlags;
//
//    bcVectorPtrs_Type           M_bcVectorPtrs;
//
//    vectorPtr_Type              M_referencePositionPtr;
//
//    Real                        M_lvPreloadPressure;
//    Real                        M_rvPreloadPressure;
//
//    vectorPtr_Type              M_lvPositionVectorPtr;//one of the orthogonal vectors with respect  to the normal of the lid
//    vectorPtr_Type              M_rvPositionVectorPtr;//one of the orthogonal vectors with respect  to the normal of the lid
//
//    bool                        M_oneWayCoupling;
//
//    commPtr_Type                M_comm;
//
//    matrixPtr_Type              M_projectionMassMatrix;
//    vectorPtr_Type              M_projectionRhs;
//    basePrecPtr_Type            M_precPtr;
//    linearSolver_Type           M_projectionSolver;
//
//    FESpacePtr_Type             M_solidActivationFESpacePtr;
//
//
//private:
    //    void initSolid();
    //    void initMonodomain();
    //    void initActivation();

};
//
//// ===================================================
////! Constructors
//// ===================================================
//template<typename Mesh, typename IonicModel>
//EMSolver<Mesh, IonicModel>::EMSolver() :
//    M_solutionMethod (SVI),
//    M_monodomainPtr(),
//    M_usingDifferentMeshes (false),
//    M_solidDataPtr(),
//    M_solidPtr(),
//    M_solidBCPtr(),
//    M_activationPtr(),
//    M_monodomainTimeStep (0.01),
//    M_solidTimeStep (1.0),
//    M_monodomainExporterPtr(),
//    M_solidExporterPtr(),
//    M_activationExporterPtr(),
//    M_monodomainDisplacementPtr(),
//    M_C2FPtr(),
//    M_activationSolidPtr(),
//    M_F2CPtr(),
//    M_fullSolidMesh(),
//    M_activationTimePtr(),
//    M_lvPressure (0.0),
//    M_rvPressure (0.0),
//    M_lvVolume (0.0),
//    M_rvVolume (0.0),
//    M_lvPVexporter(),
//    M_rvPVexporter(),
//    M_lvPV (false),
//    M_rvPV (false),
//    M_bcVectorPtrs(),
//    M_referencePositionPtr(),
//    M_lvPreloadPressure (0.0),
//    M_rvPreloadPressure (0.0),
//    M_lvPositionVectorPtr(),
//    M_rvPositionVectorPtr(),
//    M_oneWayCoupling (false),
//    M_comm(),
//    M_projectionMassMatrix(),
//    M_projectionRhs(),
//    M_precPtr(),
//    M_projectionSolver(),
//    M_solidActivationFESpacePtr()
//{}
//
//template<typename Mesh, typename IonicModel>
//EMSolver<Mesh, IonicModel>::EMSolver (   Teuchos::ParameterList& parameterList,
//                                         const std::string data_file_name, commPtr_Type comm )
//{
//    M_comm = comm;
//    if (M_comm->MyPID() == 0)
//    {
//        std::cout << "\n==========================================";
//        std::cout << "\n\t EM SOLVER: 'YOU ROCK!!!!' ";
//        std::cout << "\n==========================================";
//    }
//
//
//    M_solidTimeStep = parameterList.get ("emdt", 1.0);
//    M_monodomainTimeStep = parameterList.get ("dt", 0.01);
//    M_lvPressure = parameterList.get ("lv_pressure", 0.0);
//    M_rvPressure  = parameterList.get ("rv_pressure", 0.0);
//    M_lvVolume  = 0.0;
//    M_rvVolume  = 0.0;
//    M_lvPV = parameterList.get ("lv_pv", false);
//    M_rvPV = parameterList.get ("rv_pv", false);
//    M_oneWayCoupling = parameterList.get ("one_way_coupling", false);
//    std::string solutionMethod = parameterList.get ( "solutionMethod", "SVI" );
//    setSolutionMethod (solutionMethod);
//    GetPot dataFile (data_file_name);
//    //Initializing monodomain solver
//    if (M_comm->MyPID() == 0)
//    {
//        std::cout << "\n==========================================";
//        std::cout << "\n\t Initializing Monodomain Solver";
//        std::cout << "\n==========================================";
//    }
//
//    std::string meshName = parameterList.get ("mesh_name", "lid16.mesh");
//    std::string meshPath = parameterList.get ("mesh_path", "./");
//    ionicModelPtr_Type ionicPtr ( new IonicModel() );
//    M_monodomainPtr.reset ( new monodomainSolver_Type ( meshName, meshPath, dataFile, ionicPtr ) );
//    M_monodomainPtr -> setInitialConditions();
//    M_monodomainPtr-> setParameters ( parameterList );
//
//    M_usingDifferentMeshes = false;
//
//    //Initializing structural solver
//
//    if (M_comm->MyPID() == 0)
//    {
//        std::cout << "\n==========================================";
//        std::cout << "\n\t Initializing Solid Data and Solid mesh";
//        std::cout << "\n==========================================";
//    }
//
//    //Material data
//    M_solidDataPtr.reset ( new structureData_Type() );
//    M_solidDataPtr -> setup (dataFile);
//
//    //mesh
//    std::string solidMeshName = parameterList.get ("solid_mesh_name", "no_solid_mesh");
//    if (solidMeshName == "no_solid_mesh")
//    {
//        solidMeshName   = dataFile ( "solid/pace_discretization/mesh_file",   "no_solid_mesh" );
//    }
//
//    if (solidMeshName == "no_solid_mesh" || solidMeshName == meshName )
//    {
//        M_usingDifferentMeshes = false;
//    }
//    else
//    {
//        M_usingDifferentMeshes = true;
//    }
//
//    std::string solidMeshPath = parameterList.get ("solid_mesh_path", "");
//    if (solidMeshPath == "")
//    {
//        solidMeshPath   = dataFile ( "solid/pace_discretization/mesh_dir",  "" );
//    }
//
//
//    M_fullSolidMesh.reset (new mesh_Type ( comm ) );
//    meshPtr_Type localSolidMesh (new mesh_Type ( comm ) );
//    if ( M_usingDifferentMeshes  )
//    {
//        localSolidMesh.reset (new mesh_Type ( comm ) );
//        MeshUtility::loadMesh (localSolidMesh, M_fullSolidMesh,  solidMeshName,  solidMeshPath );
//    }
//    else
//    {
//        M_fullSolidMesh = M_monodomainPtr -> fullMeshPtr();
//        localSolidMesh = M_monodomainPtr -> localMeshPtr();
//    }
//
//    //FESPACEs
//    std::string dOrder =  dataFile ( "solid/space_discretization/order", "P1");
//    FESpacePtr_Type dFESpace ( new FESpace_Type ( localSolidMesh,
//                                                  dOrder,
//                                                  3,
//                                                  localSolidMesh -> comm() ) );
//
//    solidETFESpacePtr_Type dETFESpace ( new solidETFESpace_Type (localSolidMesh, &feTetraP1, comm) );
//
//    //boundary conditions
//    if (M_comm->MyPID() == 0)
//    {
//        std::cout << "\n==========================================";
//        std::cout << "\n\t Creating BC handler";
//        std::cout << "\n==========================================";
//        std::cout << "\n";
//    }
//    setupBC ( data_file_name, dFESpace -> map() );
//
//
//    //setup structural operator
//    if (M_comm->MyPID() == 0)
//    {
//        std::cout << "\n==========================================";
//        std::cout << "\n\t Initializing Structure Solver";
//        std::cout << "\n==========================================";
//        std::cout << "\n";
//    }
//    M_solidPtr.reset (new structuralOperator_Type() );
//    M_solidPtr -> setup (M_solidDataPtr,
//                         dFESpace,
//                         dETFESpace,
//                         M_solidBCPtr -> handler(),
//                         comm);
//    M_solidPtr -> setDataFromGetPot (dataFile);
//
//    //activation
//    if (M_comm->MyPID() == 0)
//    {
//        std::cout << "\n==========================================";
//        std::cout << "\n\t Initializing Activation Solver";
//        std::cout << "\n==========================================";
//        std::cout << "\n";
//    }
//    M_activationPtr.reset ( new activeStrain_Type (parameterList, dataFile, M_monodomainPtr -> localMeshPtr(), comm) );
//
//    M_activationTimePtr.reset (new vector_Type ( M_monodomainPtr -> feSpacePtr() -> map() ) );
//    *M_activationTimePtr = -1.0;
//
//
//
//}
//
//
//
//template<typename Mesh, typename IonicModel>
//EMSolver<Mesh, IonicModel>::EMSolver (   structuralOperatorPtr_Type solidPtr, Teuchos::ParameterList& parameterList,
//                                         const std::string data_file_name, commPtr_Type comm )
//{
//    if (M_comm->MyPID() == 0)
//    {
//        std::cout << "\n==========================================";
//        std::cout << "\n\t EM SOLVER: 'YOU ROCK!!!!' ";
//        std::cout << "\n==========================================";
//    }
//
//    M_solidTimeStep = parameterList.get ("emdt", 1.0);
//    M_monodomainTimeStep = parameterList.get ("dt", 0.01);
//    M_lvPressure = parameterList.get ("lv_pressure", 0.0);
//    M_rvPressure  = parameterList.get ("rv_pressure", 0.0);
//    M_lvVolume  = 0.0;
//    M_rvVolume  = 0.0;
//    M_lvPV = parameterList.get ("lv_pv", false);
//    M_rvPV = parameterList.get ("rv_pv", false);
//    M_oneWayCoupling = parameterList.get ("one_way_coupling", false);
//
//    GetPot dataFile (data_file_name);
//    //Initializing monodomain solver
//    if (M_comm->MyPID() == 0)
//    {
//        std::cout << "\n==========================================";
//        std::cout << "\n\t Initializing Monodomain Solver";
//        std::cout << "\n==========================================";
//    }
//
//    std::string meshName = parameterList.get ("mesh_name", "lid16.mesh");
//    std::string meshPath = parameterList.get ("mesh_path", "./");
//    ionicModelPtr_Type ionicPtr ( new IonicModel() );
//    M_monodomainPtr.reset ( new monodomainSolver_Type ( meshName, meshPath, dataFile, ionicPtr ) );
//    M_monodomainPtr -> setInitialConditions();
//    M_monodomainPtr-> setParameters ( parameterList );
//
//    M_usingDifferentMeshes = false;
//
//    //Initializing structural solver
//
//    if (M_comm->MyPID() == 0)
//    {
//        std::cout << "\n==========================================";
//        std::cout << "\n\t Initializing Solid Data and Solid mesh";
//        std::cout << "\n==========================================";
//    }
//
//    //Material data
//    M_solidDataPtr = solidPtr -> data();
//
//
//    //mesh
//    std::string solidMeshName = parameterList.get ("solid_mesh_name", "no_solid_mesh");
//    if (solidMeshName == "no_solid_mesh")
//    {
//        solidMeshName   = dataFile ( "solid/space_discretization/mesh_file",   "no_solid_mesh" );
//    }
//
//    if (solidMeshName == "no_solid_mesh" || solidMeshName == meshName )
//    {
//        M_usingDifferentMeshes = false;
//    }
//    else
//    {
//        M_usingDifferentMeshes = true;
//    }
//
//    std::string solidMeshPath = parameterList.get ("solid_mesh_path", "");
//    if (solidMeshPath == "")
//    {
//        solidMeshPath   = dataFile ( "solid/space_discretization/mesh_dir",  "" );
//    }
//
//
//    M_fullSolidMesh.reset (new mesh_Type ( M_comm ) );
//    meshPtr_Type localSolidMesh (new mesh_Type ( M_comm ) );
//    if ( M_usingDifferentMeshes  )
//    {
//        localSolidMesh.reset (new mesh_Type ( M_comm ) );
//        MeshUtility::loadMesh (localSolidMesh, M_fullSolidMesh,  solidMeshName,  solidMeshPath );
//    }
//    else
//    {
//        M_fullSolidMesh = M_monodomainPtr -> fullMeshPtr();
//        localSolidMesh = M_monodomainPtr -> localMeshPtr();
//    }
//
//    //FESPACEs
//    //    std::string dOrder =  dataFile ( "solid/space_discretization/order", "P1");
//    //    FESpacePtr_Type dFESpace ( new FESpace_Type ( localSolidMesh,
//    //                                                 dOrder,
//    //                                                 3,
//    //                                                 localSolidMesh -> M_comm() ) );
//    //
//    //    solidETFESpacePtr_Type dETFESpace ( new solidETFESpace_Type (localSolidMesh, &feTetraP1, M_comm) );
//
//    //boundary conditions
//    if (M_comm->MyPID() == 0)
//    {
//        std::cout << "\n==========================================";
//        std::cout << "\n\t Creating BC handler";
//        std::cout << "\n==========================================";
//        std::cout << "\n";
//    }
//
//    M_solidBCPtr.reset ( new bcInterface_Type() );
//    M_solidBCPtr-> setHandler ( solidPtr -> bcHandler() );
//    M_solidBCPtr -> setPhysicalSolver (solidPtr);
//
//    //setup structural operator
//    if (M_comm->MyPID() == 0)
//    {
//        std::cout << "\n==========================================";
//        std::cout << "\n\t Initializing Structure Solver";
//        std::cout << "\n==========================================";
//        std::cout << "\n";
//    }
//    M_solidPtr = solidPtr;
//    //    M_solidPtr -> setup (M_solidDataPtr,
//    //                      dFESpace,
//    //                      dETFESpace,
//    //                      M_solidBCPtr -> handler(),
//    //                      M_comm);
//    //    M_solidPtr -> setDataFromGetPot (dataFile);
//
//    //activation
//    if (M_comm->MyPID() == 0)
//    {
//        std::cout << "\n==========================================";
//        std::cout << "\n\t Initializing Activation Solver";
//        std::cout << "\n==========================================";
//        std::cout << "\n";
//    }
//    M_activationPtr.reset ( new activeStrain_Type (parameterList, dataFile, M_monodomainPtr -> localMeshPtr(), M_comm) );
//
//
//    M_activationTimePtr.reset (new vector_Type ( M_monodomainPtr -> feSpacePtr() -> map() ) );
//    *M_activationTimePtr = -1.0;
//}
//
//
//template<typename Mesh, typename IonicModel>
//void EMSolver<Mesh, IonicModel>::setup (Teuchos::ParameterList& parameterList,
//                                        const std::string data_file_name,
//                                        std::string parameterListName )
//{
//    if (M_usingDifferentMeshes)
//    {
//        M_monodomainDisplacementPtr.reset ( new vector_Type ( M_monodomainPtr -> displacementETFESpacePtr() -> map() ) );
//        M_activationSolidPtr.reset ( new vector_Type ( M_solidPtr -> activeMaterial() -> activationSpace() -> map() ) );
//        GetPot dataFile (data_file_name);
//        setupInterpolants (parameterListName, parameterList, dataFile);
//    }
//    else
//    {
//        M_monodomainDisplacementPtr = M_solidPtr -> displacementPtr();
//        M_activationSolidPtr = M_activationPtr -> gammafPtr();
//    }
//
//    if (M_oneWayCoupling == false)
//    {
//        M_monodomainPtr -> setDisplacementPtr ( M_monodomainDisplacementPtr );
//    }
//    M_solidPtr -> activeMaterial() -> setGammaf (*M_activationSolidPtr);
//    setupMonodomainMatrix (parameterList);
//
//    createReferencePositionVector();
//
//}
//
//template<typename Mesh, typename IonicModel>
//void EMSolver<Mesh, IonicModel>::setupProjectionSolver (GetPot& dataFile, commPtr_Type comm)
//{
//    if (M_comm -> MyPID() == 0)
//    {
//        std::cout << "\nEM Solver: Setting L2 projection Linear Solver";
//    }
//    prec_Type* precRawPtr;
//    precRawPtr = new prec_Type;
//    precRawPtr->setDataFromGetPot (dataFile, "prec");
//    M_precPtr.reset (precRawPtr);
//
//    Teuchos::RCP < Teuchos::ParameterList > solverParamList = Teuchos::rcp (
//                                                                  new Teuchos::ParameterList);
//
//    std::string xmlpath = dataFile ("activation/activation_xml_path",
//                                    "./");
//    std::string xmlfile = dataFile ("activation/activation_xml_file",
//                                    "ParamList.xml");
//
//    solverParamList = Teuchos::getParametersFromXmlFile (xmlpath + xmlfile);
//    //linearSolver_Type linearSolver;
//    M_projectionSolver.setCommunicator ( comm );
//    M_projectionSolver.setParameters ( *solverParamList );
//    M_projectionSolver.setPreconditioner ( M_precPtr );
//    M_projectionSolver.setOperator ( M_projectionMassMatrix );
//}
//
//template<typename Mesh, typename IonicModel>
//void EMSolver<Mesh, IonicModel>::setupMonodomainMatrix (Teuchos::ParameterList& parameterList)
//{
//    M_monodomainPtr -> setupMassMatrix();
//    M_monodomainPtr -> setupStiffnessMatrix();
//    M_monodomainPtr -> setupGlobalMatrix();
//}
//
//
//template<typename Mesh, typename IonicModel>
//void EMSolver<Mesh, IonicModel>::updateMonodomainMatrix()
//{
//    if (M_oneWayCoupling == false)
//    {
//        M_monodomainPtr -> setDisplacementPtr ( M_monodomainDisplacementPtr );
//        M_monodomainPtr -> setupMassMatrix();
//        M_monodomainPtr -> setupStiffnessMatrix();
//        M_monodomainPtr -> setupGlobalMatrix();
//    }
//}
//
//template<typename Mesh, typename IonicModel>
//void EMSolver<Mesh, IonicModel>::updateMonodomain()
//{
//    if (M_usingDifferentMeshes)
//    {
//        M_C2FPtr -> updateRhs ( M_solidPtr -> displacementPtr() );
//        M_C2FPtr -> interpolate();
//        M_C2FPtr -> solution (  M_monodomainDisplacementPtr );
//    }
//    updateMonodomainMatrix();
//}
//
//
//template<typename Mesh, typename IonicModel>
//void EMSolver<Mesh, IonicModel>::updateSolid()
//{
//    if (M_usingDifferentMeshes)
//    {
//        //      if(M_F2CPtr)
//        //      {
//        M_F2CPtr -> updateRhs ( M_activationPtr -> gammafPtr() );
//        M_F2CPtr -> interpolate();
//        M_F2CPtr -> solution (  M_activationSolidPtr );
//
//        //  setupProjectionRhs();
//        //  projection();
//        //      }
//        //      else
//        //      {
//        //  //  feToFEInterpolate
//        //      EMEvaluate utility;
//        //      utility.evaluate(M_activationPtr -> gammaf(), *M_activationSolidPtr, *(M_activationPtr ->FESpacePtr()), *M_solidActivationFESpacePtr, *M_fullSolidMesh);
//        //      }
//    }
//
//    M_solidPtr -> activeMaterial() -> setGammaf ( *M_activationSolidPtr );
//
//    M_activationPtr -> computeGammasAndGamman (M_solidPtr -> activeMaterial() -> gammaf(),
//                                               M_solidPtr -> activeMaterial() -> gammas(),
//                                               M_solidPtr -> activeMaterial() -> gamman() );
//}
//
//
//template<typename Mesh, typename IonicModel>
//void EMSolver<Mesh, IonicModel>::solveOneMonodomainStep (int subiter)
//{
//    if ( M_comm->MyPID() == 0 )
//    {
//        cout << "\n------------------";
//        cout << "\nMonodomain Solver ";
//        cout << "\n------------------";
//    }
//    LifeChrono timer;
//    timer.start();
//    switch (M_solutionMethod)
//    {
//        case    OS:
//
//            if ( M_comm->MyPID() == 0 )
//            {
//                std::cout << "\nEM SOLVER: OS: Solving reactions";
//            }
//            for (int j (0); j < subiter; j++)
//            {
//                solveOneReactionStepFE (subiter);
//            }
//            if ( M_comm->MyPID() == 0 )
//            {
//                std::cout << "\t ... done in " << timer.diff() << "s\n";
//                std::cout << "\nEM SOLVER: OS: Solving diffusion";
//            }
//            solveOneDiffusionStep();
//            break;
//
//        case ICI:
//
//            if (M_comm -> MyPID() == 0 )
//            {
//                std::cout << "\nEM SOLVER: ICI: solving gating variables";
//            }
//            M_monodomainPtr -> solveOneStepGatingVariablesFE();
//            if (M_comm->MyPID() == 0 )
//            {
//                std::cout << "\t ... done in " << timer.diff() << "s\n";
//                std::cout << "\nEM SOLVER: ICI: Solving potential";
//            }
//            M_monodomainPtr ->solveOneICIStep();
//            break;
//
//        case HL:
//
//            * (M_monodomainPtr -> appliedCurrentPtr() ) *= 10.0;
//            if (M_comm -> MyPID() == 0 )
//            {
//                std::cout << "\nEM SOLVER: HL: solving gating variables";
//            }
//            M_monodomainPtr -> solveOneStepGatingVariablesFE();
//            if (M_comm->MyPID() == 0 )
//            {
//                std::cout << "\t ... done in " << timer.diff() << "s\n";
//                std::cout << "\nEM SOLVER: HL: Solving potential";
//            }
//            M_monodomainPtr ->solveOneICIStep ( * (M_activationPtr -> massMatrixPtr() ) );
//            break;
//
//        default:
//
//            if (M_comm -> MyPID() == 0 )
//            {
//                std::cout << "\nEM SOLVER: SVI: solving gating variables";
//            }
//            M_monodomainPtr -> solveOneStepGatingVariablesFE();
//            if (M_comm->MyPID() == 0 )
//            {
//                std::cout << "\t ... done in " << timer.diff() << "s\n";
//                std::cout << "\nEM SOLVER: SVI: Solving potential";
//            }
//            M_monodomainPtr ->solveOneSVIStep();
//            break;
//    }
//
//    if ( M_comm->MyPID() == 0 )
//    {
//        std::cout << "\t ... done in " << timer.diffCumul() << "s\n";
//        std::cout << "\nEM SOLVER: Total time " << timer.diff() << "s\n";
//    }
//}
//
//template<typename Mesh, typename IonicModel>
//void EMSolver<Mesh, IonicModel>::solveOneReactionStepFE (int subiter)
//{
//    if (M_comm -> MyPID() == 0 )
//    {
//        std::cout << "\nEM SOLVER: solving monodomain with OS: reaction step\n";
//    }
//    for (int j (0); j < subiter; j++)
//    {
//        M_monodomainPtr -> solveOneReactionStepFE (subiter);
//    }
//}
//
//template<typename Mesh, typename IonicModel>
//void EMSolver<Mesh, IonicModel>::solveOneICIStep()
//{
//    if (M_comm -> MyPID() == 0 )
//    {
//        std::cout << "\nEM SOLVER: ICI: solving monodomain with ICI\n";
//    }
//    M_monodomainPtr -> solveOneStepGatingVariablesFE();
//    M_monodomainPtr ->solveOneICIStep();
//}
//
//template<typename Mesh, typename IonicModel>
//void EMSolver<Mesh, IonicModel>::solveOneSVIStep()
//{
//    if (M_comm -> MyPID() == 0 )
//    {
//        std::cout << "\nEM SOLVER: solving monodomain with SVI\n";
//    }
//    M_monodomainPtr -> solveOneStepGatingVariablesFE();
//    M_monodomainPtr ->solveOneSVIStep();
//}
//
//template<typename Mesh, typename IonicModel>
//void EMSolver<Mesh, IonicModel>::solveOneHLStep()
//{
//    if (M_comm -> MyPID() == 0 )
//    {
//        std::cout << "\nEM SOLVER: solving monodomain with ICI half lumping\n";
//    }
//    * (M_monodomainPtr -> appliedCurrentPtr() ) *= 10.0;
//    M_monodomainPtr -> solveOneStepGatingVariablesFE();
//    M_monodomainPtr ->solveOneICIStep ( * (M_activationPtr -> massMatrixPtr() ) );
//}
//
//template<typename Mesh, typename IonicModel>
//void EMSolver<Mesh, IonicModel>::solveOneOSStep (int subiter)
//{
//    solveOneReactionStepFE (subiter);
//    solveOneDiffusionStep();
//}
//
//
//template<typename Mesh, typename IonicModel>
//void EMSolver<Mesh, IonicModel>::solveOneDiffusionStep()
//{
//    if (M_comm -> MyPID() == 0 )
//    {
//        std::cout << "\nEM SOLVER: solving monodomain with OS: reaction step\n";
//    }
//    * (M_monodomainPtr -> rhsPtrUnique() ) *= 0.0;
//    M_monodomainPtr -> updateRhs();
//    M_monodomainPtr -> solveOneDiffusionStepBE();
//}
//
//template<typename Mesh, typename IonicModel>
//void EMSolver<Mesh, IonicModel>::solveOneActivationStep()
//{
//    M_activationPtr -> setupRhs (*M_monodomainDisplacementPtr,
//                                 M_monodomainPtr -> displacementETFESpacePtr(),
//                                 * (M_monodomainPtr -> globalSolution().at (3) ),
//                                 M_monodomainPtr -> ETFESpacePtr(),
//                                 M_monodomainPtr -> timeStep() );
//    M_activationPtr -> solve();
//}
//
//
//
//
//template<typename Mesh, typename IonicModel>
//void EMSolver<Mesh, IonicModel>::setupExporters (std::string dir)
//{
//    if (M_comm)
//    {
//        std::cout << "\n==========================================";
//        std::cout << "\n\t Setting up the exporters";
//        std::cout << "\n==========================================";
//    }
//    M_monodomainExporterPtr.reset (new exporter_Type() );
//    M_monodomainPtr -> setupExporter (*M_monodomainExporterPtr, "ElectroOutput", dir);
//    M_activationExporterPtr.reset (new exporter_Type() );
//    M_activationPtr -> setupExporter (*M_activationExporterPtr, M_comm, dir);
//
//    M_solidExporterPtr.reset (new exporter_Type() );
//    M_solidExporterPtr -> setMeshProcId ( M_solidPtr -> mesh(), M_comm->MyPID() );
//    M_solidExporterPtr -> setPrefix ( "StructureOutput" );
//    M_solidExporterPtr -> setPostDir ( dir );
//    M_solidExporterPtr -> addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, "displacement", M_solidPtr -> dispFESpacePtr(), M_solidPtr -> displacementPtr(), UInt (0) );
//
//
//    if (M_usingDifferentMeshes)
//    {
//        FESpacePtr_Type gfSolidFESpace ( new FESpace_Type ( M_solidPtr -> dispFESpace().mesh(),
//                                                            "P1",   1,  M_comm) );
//        M_solidExporterPtr -> addVariable ( ExporterData<RegionMesh<LinearTetra> >::ScalarField,
//                                            "gammaf",
//                                            gfSolidFESpace,
//                                            M_activationSolidPtr,
//                                            UInt (0) );
//
//        FESpacePtr_Type displacementMonodomainFESpace ( new FESpace_Type ( M_monodomainPtr -> localMeshPtr(),
//                                                                           "P1",   3,   M_comm ) );
//        M_monodomainExporterPtr -> addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField,
//                                                 "interpolated_displacement",
//                                                 displacementMonodomainFESpace,
//                                                 M_monodomainDisplacementPtr,
//                                                 UInt (0) );
//    }
//
//    if (M_lvPV)
//    {
//        if ( M_comm->MyPID() == 0 )
//        {
//            std::string outputFile = dir + "/LV_PV.txt";
//            M_lvPVexporter.open ( outputFile.c_str() );
//            M_lvPVexporter << "Time" << ", " << "LV_pressure" << ", " << "LV_volume" << "\n";
//        }
//    }
//    if (M_rvPV)
//    {
//        if ( M_comm->MyPID() == 0 )
//        {
//            std::string outputFile = dir + "/RV_PV.txt";
//            M_rvPVexporter.open ( outputFile.c_str() );
//            M_rvPVexporter << "Time" << ", " << "RV_pressure" << ", " << "RV_volume" << "\n";
//        }
//    }
//
//
//}
//
//
//template<typename Mesh, typename IonicModel>
//void EMSolver<Mesh, IonicModel>::exportSolution (Real time)
//{
//    M_monodomainExporterPtr -> postProcess (time);
//    M_activationExporterPtr -> postProcess (time);
//    M_solidExporterPtr -> postProcess (time);
//
//    if (M_lvPV)
//    {
//        if ( M_comm->MyPID() == 0 )
//        {
//            M_lvPVexporter << time << ", " << - M_lvPressure  * 3 / 4 << ", " << M_lvVolume << "\n";
//        }
//    }
//    if (M_rvPV)
//    {
//        if ( M_comm->MyPID() == 0 )
//        {
//            M_rvPVexporter << time << ", " << - M_rvPressure * 3 / 4 << ", " << M_rvVolume << "\n";
//        }
//    }
//
//}
//
//
//template<typename Mesh, typename IonicModel>
//void EMSolver<Mesh, IonicModel>::closeExporters()
//{
//    M_monodomainExporterPtr -> closeFile();
//    M_activationExporterPtr -> closeFile();
//    M_solidExporterPtr -> closeFile();
//    if (M_lvPV)
//    {
//        if ( M_comm->MyPID() == 0 )
//        {
//            M_lvPVexporter.close();
//        }
//    }
//    if (M_rvPV)
//    {
//        if ( M_comm->MyPID() == 0 )
//        {
//            M_rvPVexporter.close();
//        }
//    }
//}
//
//
//template<typename Mesh, typename IonicModel>
//void EMSolver<Mesh, IonicModel>::setupPreloadBC ( GetPot& dataFile, map_Type map  )
//{
//
//
//    M_lvPreloadPressure = dataFile ( "solid/boundary_conditions/lv_preload_pressure", 0.0 );
//    vectorPtr_Type  lvPressureVector (new vector_Type ( map, Repeated ) );
//    UInt numTotalDof = 3 * map.mapSize();
//    *lvPressureVector = -M_lvPreloadPressure;
//    if (M_lvFlags.empty() == false)
//    {
//        int numLVFlags = M_lvFlags.size();
//        M_bcVectorPtrs.resize (numLVFlags);
//
//        for (int i (0); i < numLVFlags; i++)
//        {
//            M_bcVectorPtrs.at (i).reset (new bcVector_Type ( *lvPressureVector, numTotalDof, 1) );
//            M_solidBCPtr -> handler() -> addBC ("LV_Endocardium", M_lvFlags[i], Natural, Full, *M_bcVectorPtrs[i], 3);
//        }
//    }
//
//    M_rvPreloadPressure = dataFile ( "solid/boundary_conditions/rv_preload_pressure", 0.0 );
//    vectorPtr_Type  rvPressureVector (new vector_Type ( map, Repeated ) );
//    *rvPressureVector = -M_rvPreloadPressure;
//    if (M_rvFlags.empty() == false)
//    {
//        int numRVFlags = M_rvFlags.size();
//        for (int i (0); i < numRVFlags; i++)
//        {
//            M_bcVectorPtrs.push_back ( * (new bcVectorPtr_Type ( new bcVector_Type ( *rvPressureVector, numTotalDof, 1) ) ) ); // );
//            M_solidBCPtr -> handler() -> addBC ("RV_Endocardium", M_rvFlags[i], Natural, Full, * (M_bcVectorPtrs.back() ), 3);
//        }
//    }
//}
//template<typename Mesh, typename IonicModel>
//void EMSolver<Mesh, IonicModel>::preloadRamp ( Real dt  )
//{
//    M_solidPtr -> displayer().leaderPrint ("\n==========================================");
//    M_solidPtr -> displayer().leaderPrint ("\n\t Starting Preload Ramp ");
//    M_solidPtr -> displayer().leaderPrint ("\n==========================================");
//    M_solidPtr -> displayer().leaderPrint ("\n");
//
//    if (M_lvFlags.empty() == false || M_rvFlags.empty() == false)
//    {
//        vectorPtr_Type  aux (new vector_Type ( M_solidPtr -> displacementPtr() -> map(), Repeated ) );
//        for (Real t (0.0); t < 1; )
//        {
//            t += dt;
//
//
//            M_solidPtr -> displayer().leaderPrint ("\n==========================================");
//            M_solidPtr -> displayer().leaderPrint ("\nRamp Time: ", t);
//            M_solidPtr -> displayer().leaderPrint ("\ntimestep: ", dt);
//            M_solidPtr -> displayer().leaderPrint ("\n==========================================");
//            M_solidPtr -> displayer().leaderPrint ("\n");
//
//
//            Real lvp = t * M_lvPreloadPressure;
//            Real rvp = t * M_rvPreloadPressure;
//            int numLVFlags = M_lvFlags.size();
//            for (int i (0); i < numLVFlags; i++)
//            {
//                *aux = -lvp;
//                M_bcVectorPtrs[i].reset ( new bcVector_Type (*aux, M_solidPtr -> dispFESpacePtr() -> dof().numTotalDof(), 1) );
//                M_solidBCPtr -> handler() -> modifyBC (M_lvFlags[i], * (M_bcVectorPtrs[i]) );
//            }
//            for (int i (numLVFlags); i < numLVFlags + M_rvFlags.size(); i++)
//            {
//                *aux = -rvp;
//                M_bcVectorPtrs[i].reset ( new bcVector_Type (*aux, M_solidPtr -> dispFESpacePtr() -> dof().numTotalDof(), 1) );
//                M_solidBCPtr -> handler() -> modifyBC (M_rvFlags[i - numLVFlags], * (M_bcVectorPtrs[i]) );
//            }
//            solveSolid();
//        }
//    }
//}
//
//template<typename Mesh, typename IonicModel>
//void EMSolver<Mesh, IonicModel>::createReferencePositionVector()
//{
//    M_referencePositionPtr.reset ( new vector_Type ( M_solidPtr -> displacementPtr() -> map() ) );
//    Int nLocalDof = M_referencePositionPtr -> epetraVector().MyLength();
//    Int nComponentLocalDof = nLocalDof / 3;
//    for (int k (0); k < nComponentLocalDof; k++)
//    {
//        UInt iGID = M_referencePositionPtr -> blockMap().GID (k);
//        UInt jGID = M_referencePositionPtr -> blockMap().GID (k + nComponentLocalDof);
//        UInt kGID = M_referencePositionPtr -> blockMap().GID (k + 2 * nComponentLocalDof);
//
//        (* (M_referencePositionPtr) ) [iGID] = M_fullSolidMesh -> point (iGID).x();
//        (* (M_referencePositionPtr) ) [jGID] = M_fullSolidMesh -> point (iGID).y();
//        (* (M_referencePositionPtr) ) [kGID] = M_fullSolidMesh -> point (iGID).z();
//    }
//}
//
//template<typename Mesh, typename IonicModel>
//void EMSolver<Mesh, IonicModel>::setupBC ( const std::string data_file_name, map_Type map   )
//{
//    GetPot dataFile (data_file_name);
//    M_solidBCPtr.reset ( new bcInterface_Type() );
//    M_solidBCPtr->createHandler();
//    M_solidBCPtr->fillHandler ( data_file_name, "solid" );
//
//    UInt lvFlagsNumber = dataFile.vector_variable_size ( "solid/boundary_conditions/lv_flags"  );
//    //  UInt pippo = dataFile.vector_variable_size ( ( "/" + "/" + subSection + "/list" ).data() )
//    if (lvFlagsNumber == 0)
//    {
//        if (M_comm->MyPID() == 0)
//        {
//            std::cout << "\nWARNING: You have not set the LV Flags!\n";
//        }
//    }
//    else
//    {
//        M_lvFlags.resize (lvFlagsNumber);
//        for (UInt i (0); i < lvFlagsNumber; i++)
//        {
//            M_lvFlags[i] =  dataFile ( ( "solid/boundary_conditions/lv_flags" ), -1, i );
//        }
//    }
//
//    UInt rvFlagsNumber = dataFile.vector_variable_size ( "solid/boundary_conditions/rv_flags" );
//    if (rvFlagsNumber == 0)
//    {
//        if (M_comm->MyPID() == 0)
//        {
//            std::cout << "\nWARNING: You have not set the RV Flags!\n";
//        }
//    }
//    else
//    {
//        M_rvFlags.resize (rvFlagsNumber);
//        for (UInt i (0); i < rvFlagsNumber; i++)
//        {
//            M_rvFlags[i] =  dataFile ( "solid/boundary_conditions/rv_flags", -1, i );
//        }
//    }
//    setupPreloadBC ( dataFile, map );
//}
//
//
//template<typename Mesh, typename IonicModel>
//void EMSolver<Mesh, IonicModel>::exportSolidFibersDirection (std::string dir )
//{
//    exporter_Type exp;
//    exp.setMeshProcId ( M_solidPtr -> mesh(), M_comm->MyPID() );
//    exp.setPostDir ( dir );
//    exp.setPrefix ("SolidFibesrDirection");
//    exp.addVariable (ExporterData<mesh_Type>::VectorField,
//                     "solid_fibers",
//                     M_solidPtr -> dispFESpacePtr(),
//                     M_solidPtr -> activeMaterial() -> fiberVectorPtr(), UInt (0) );
//    exp.postProcess (0);
//    exp.closeFile();
//}
//
//
//template<typename Mesh, typename IonicModel>
//void EMSolver<Mesh, IonicModel>::exportMonodomainFibersDirection (std::string dir)
//{
//    M_monodomainPtr -> exportFiberDirection (dir);
//}
//
//
//template<typename Mesh, typename IonicModel>
//void EMSolver<Mesh, IonicModel>::exportSolidSheetsDirection (std::string dir )
//{
//    exporter_Type exp;
//    exp.setMeshProcId ( M_solidPtr -> mesh(), M_comm->MyPID() );
//    exp.setPostDir ( dir );
//    exp.setPrefix ("SolidSheetsDirection");
//    exp.addVariable (ExporterData<mesh_Type>::VectorField,
//                     "solid_sheets",
//                     M_solidPtr -> dispFESpacePtr(),
//                     M_solidPtr -> activeMaterial() -> sheetVectorPtr(), UInt (0) );
//    exp.postProcess (0);
//    exp.closeFile();
//}
//
//template<typename Mesh, typename IonicModel>
//void EMSolver<Mesh, IonicModel>::exportFibersAndSheetsFields (std::string dir )
//{
//    exportSolidSheetsDirection ( dir);
//    exportSolidFibersDirection (dir);
//    exportMonodomainFibersDirection (dir);
//}
//
//
//template<typename Mesh, typename IonicModel>
//void EMSolver<Mesh, IonicModel>::exportActivationTime (std::string dir )
//{
//    exporter_Type exp;
//    exp.setMeshProcId ( M_monodomainPtr -> localMeshPtr(), M_comm->MyPID() );
//    exp.setPostDir ( dir );
//    exp.setPrefix ("ActivationTime");
//    exp.addVariable (ExporterData<mesh_Type>::ScalarField,
//                     "activation_time",
//                     M_monodomainPtr -> feSpacePtr(),
//                     M_activationTimePtr, UInt (0) );
//    exp.postProcess (0);
//    exp.closeFile();
//}
//
//template<typename Mesh, typename IonicModel>
//void EMSolver<Mesh, IonicModel>::importSolidFibers (Teuchos::ParameterList& parameterList)
//{
//    std::string solidFibersFile = parameterList.get ("solid_fibers_file", "");
//    std::string solidFibersField = parameterList.get ("solid_fibers_field", "");
//    ElectrophysiologyUtility::importVectorField ( M_solidPtr -> activeMaterial() -> fiberVectorPtr(),
//                                      solidFibersFile,
//                                      solidFibersField,
//                                      M_solidPtr -> mesh() );
//}
//template<typename Mesh, typename IonicModel>
//void EMSolver<Mesh, IonicModel>::importSolidSheets (Teuchos::ParameterList& parameterList)
//{
//    std::string solidSheetsFile = parameterList.get ("solid_sheets_file", "");
//    std::string solidSheetsField = parameterList.get ("solid_sheets_field", "");
//    if ( solidSheetsFile != "" )
//    {
//        ElectrophysiologyUtility::importVectorField ( M_solidPtr -> activeMaterial() -> sheetVectorPtr(),
//                                          solidSheetsFile,
//                                          solidSheetsField,
//                                          M_solidPtr -> mesh() );
//    }
//}
//
//template<typename Mesh, typename IonicModel>
//void EMSolver<Mesh, IonicModel>::importMonodomainFibers (Teuchos::ParameterList& parameterList)
//{
//    std::string fibersFile = parameterList.get ("fibers_file", "");
//    std::string fibersField = parameterList.get ("fibers_field", "");
//    vectorPtr_Type monodomainFibers (new vector_Type ( M_monodomainPtr -> displacementETFESpacePtr() -> map() ) );
//    ElectrophysiologyUtility::importFibers (monodomainFibers, fibersFile, M_monodomainPtr -> localMeshPtr() );
//    //    HeartUtility::importVectorField( monodomainFibers,
//    //                                   fibersFile,
//    //                                   fibersField,
//    //                                   M_monodomainPtr -> localMeshPtr() );
//    setMonodomainFibers (monodomainFibers);
//}
//
//template<typename Mesh, typename IonicModel>
//void EMSolver<Mesh, IonicModel>::setFibersAndSheets (Teuchos::ParameterList& parameterList)
//{
//    std::cout << "\nImporting Solid Fibers\n\n";
//    this->importSolidFibers (parameterList);
//
//    std::cout << "\nImporting Solid Sheets\n\n";
//    importSolidSheets (parameterList);
//
//    std::cout << "\nImporting Monodomain Fibers\n\n";
//    importMonodomainFibers (parameterList);
//
//
//    M_activationPtr -> setFiberPtr ( M_monodomainPtr -> fiberPtr() );
//}
//
//template<typename Mesh, typename IonicModel>
//void EMSolver<Mesh, IonicModel>::setupInterpolants (std::string parameterListName, Teuchos::ParameterList& parameterList, GetPot& dataFile)
//{
//
//    //  if(M_monodomainPtr -> commPtr() -> MyPID() == 0)
//    //  {
//    //        cout << "\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXX";
//    //      std::cout << "\n Activation Solid NORM2: " << M_activationSolidPtr -> norm2();
//    //        cout << "\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXX";
//    //  }
//
//
//    Teuchos::RCP< Teuchos::ParameterList > belosList = Teuchos::rcp ( new Teuchos::ParameterList );
//    belosList = Teuchos::getParametersFromXmlFile ( parameterListName );
//
//    int nFlags = 1;
//    std::vector<int> flags (nFlags);
//    flags[0] = -1;
//
//    std::string c2f = parameterList.get ("c2f", "RBFrescaledVectorial");
//    if (M_comm -> MyPID() == 0)
//    {
//        std::cout << "\nINTERPOLATION C2F: from coarse to fine using " << c2f;
//    }
//    M_C2FPtr.reset ( interpolation_Type::InterpolationFactory::instance().createObject ( c2f ) );
//    M_C2FPtr->setup ( M_fullSolidMesh,
//                      M_solidPtr -> mesh(),
//                      M_monodomainPtr -> fullMeshPtr(),
//                      M_monodomainPtr -> localMeshPtr(),
//                      flags);
//    M_C2FPtr -> setRadius ( 2.0 * (double) MeshUtility::MeshStatistics::computeSize (* (M_fullSolidMesh) ).maxH );
//    M_C2FPtr -> setupRBFData ( M_solidPtr -> displacementPtr(), M_monodomainDisplacementPtr, dataFile, belosList);
//    if (c2f == "RBFvectorial")
//    {
//        M_C2FPtr->setBasis ("TPS");
//    }
//    M_C2FPtr->buildOperators();
//    M_C2FPtr->interpolate();
//    M_C2FPtr->solution (M_monodomainDisplacementPtr);
//
//
//    ////////////////////////////////////////////////////////////////////////////////
//    std::string f2c = parameterList.get ("f2c", "RBFrescaledScalar");
//    if (M_monodomainPtr -> commPtr() -> MyPID() == 0)
//    {
//        std::cout << "\nINTERPOLATION F2C: from fine to coarse using " << f2c;
//    }
//
//    M_F2CPtr.reset ( interpolation_Type::InterpolationFactory::instance().createObject ( f2c ) );
//    M_F2CPtr->setup ( M_monodomainPtr -> fullMeshPtr(),
//                      M_monodomainPtr -> localMeshPtr(),
//                      M_fullSolidMesh,
//                      M_solidPtr -> mesh(),
//                      flags);
//    //WARNING
//    M_solidPtr -> displayer().leaderPrint ("\nWARNING!!! Setting the Radius of interpolation using the full monodomain mesh.");
//    M_solidPtr -> displayer().leaderPrint ("\nWARNING!!! You shoul use the full activation mesh, but it's not coded yet...\n");
//
//
//    M_F2CPtr -> setRadius ( (double) MeshUtility::MeshStatistics::computeSize (* ( M_monodomainPtr -> fullMeshPtr() ) ).maxH );
//    M_F2CPtr -> setupRBFData ( M_activationPtr -> gammafPtr(), M_activationSolidPtr , dataFile, belosList);
//    M_F2CPtr -> buildOperators();
//    M_F2CPtr->interpolate();
//    M_F2CPtr->solution (M_activationSolidPtr);
//
//    ////////////////////////////////////////////////////////////////////////////////
//
//
//
//    //  M_projectionMassMatrix.reset(new matrix_Type(  M_solidPtr -> material() -> activationSpace() -> map() ) ) ;
//    //      {
//    //          using namespace ExpressionAssembly;
//    //
//    //          integrate( elements( M_solidPtr -> mesh() ),
//    //                      M_activationPtr -> FESpacePtr() -> qr(),
//    //                      M_solidPtr -> material() -> activationSpace(),
//    //                      M_solidPtr -> material() -> activationSpace(),
//    //                      phi_i * phi_j ) >> M_projectionMassMatrix;
//    //
//    //      }
//    //      std::cout << "\nProjection Mass Matrix Assembled!!!\n";
//    //      M_projectionMassMatrix -> globalAssemble();
//
//
//    //      setupProjectionSolver(dataFile,M_comm);
//    //      setupProjectionRhs();
//    //   cout << "\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXX";
//    //  M_projectionMassMatrix ->  spy("l2proj");
//    //  std::cout << "\n Projection Rhs SIZE: " << M_projectionRhs ->  size();
//    //  std::cout << "\n Projection Rhs NORM2: " << M_projectionRhs -> norm2();
//    //   cout << "\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXX";
//
//    //      projection();
//    //      M_solidActivationFESpacePtr.reset(new FESpace_Type(M_solidPtr -> mesh(), "P1" , 1, M_comm) );
//
//}
//
//template<typename Mesh, typename IonicModel>
//void EMSolver<Mesh, IonicModel>::setupProjectionRhs()
//{
//
//    vectorPtr_Type tmpRhs ( new vector_Type ( M_solidPtr -> activeMaterial() -> activationSpace() -> map(), Repeated ) );
//
//    if (M_activationPtr -> gammafPtr() )
//    {
//        {
//            using namespace ExpressionAssembly;
//
//            integrate ( elements ( M_solidPtr -> mesh() ),
//                        M_activationPtr -> FESpacePtr() -> qr(),
//                        M_solidPtr -> activeMaterial() -> activationSpace(),
//                        value (M_activationPtr -> ETFESpacePtr(), M_activationPtr -> gammaf() ) * phi_i  ) >> tmpRhs;
//
//        }
//    }
//    else
//    {
//        *tmpRhs *= 0.0;
//    }
//
//
//    M_projectionRhs.reset ( new vector_Type ( *tmpRhs, Unique ) );
//
//    //std::cout << "\nProjection rhs vector Assembled!!!\n";
//    M_projectionRhs -> globalAssemble();
//    M_projectionSolver.setRightHandSide (M_projectionRhs);
//}
//
//template<typename Mesh, typename IonicModel>
//void EMSolver<Mesh, IonicModel>::projection()
//{
//    if (M_comm -> MyPID() == 0)
//    {
//        std::cout << "\n L2 projection: Solving";
//    }
//
//    M_projectionSolver.solve (M_activationSolidPtr);
//}
//
//template<typename Mesh, typename IonicModel>
//void EMSolver<Mesh, IonicModel>::createPositionVector (Real nx, Real ny, Real nz, std::vector<UInt> flags, vector_Type& vec)
//{
//    Real e1x = 1.0;
//    Real e1y = 0.0;
//    Real e1z = 0.0;
//    Real e2x = 0.0;
//    Real e2y = 1.0;
//    Real e2z = 0.0;
//    Real e3x = 0.0;
//    Real e3y = 0.0;
//    Real e3z = 1.0;
//
//    Real sp = nx * e1x + ny * e1y + nz * e1z;
//    Real t1x, t1y, t1z;
//    if (sp != 1)
//    {
//        t1x = e1x - sp * nx;
//        t1y = e1y - sp * ny;
//        t1z = e1z - sp * nz;
//        Real  norm = std::sqrt ( t1x * t1x + t1y * t1y + t1z * t1z );
//        t1x = t1x / norm;
//        t1y = t1y / norm;
//        t1z = t1z / norm;
//    }
//    else
//    {
//        t1x = e2x;
//        t1y = e2y;
//        t1z = e2z;
//    }
//    Real t1e1 = t1x * e1x + t1y * e1y + t1z * e1z;
//    Real t1e2 = t1x * e2x + t1y * e2y + t1z * e2z;
//    Real t1e3 = t1x * e3x + t1y * e3y + t1z * e3z;
//
//    Int nLocalDof = M_referencePositionPtr -> epetraVector().MyLength();
//    Int nComponentLocalDof = nLocalDof / 3;
//
//    for (int k (0); k < nComponentLocalDof; k++)
//    {
//        UInt iGID = M_referencePositionPtr-> blockMap().GID (k);
//        UInt jGID = M_referencePositionPtr-> blockMap().GID (k + nComponentLocalDof);
//        UInt kGID = M_referencePositionPtr-> blockMap().GID (k + 2 * nComponentLocalDof);
//
//        vec[iGID] = (* (M_referencePositionPtr) ) [iGID] * t1e1
//                    + (* (M_referencePositionPtr) ) [jGID] * t1e2
//                    + (* (M_referencePositionPtr) ) [kGID] * t1e3;
//    }
//}
//
//template<typename Mesh, typename IonicModel>
//void EMSolver<Mesh, IonicModel>::computeLVVolume (Real nx, Real ny, Real nz)
//{
//
//    Real Vn = M_lvVolume;
//    MatrixSmall<3, 3> Id;
//    Id (0, 0) = 1.;
//    Id (0, 1) = 0., Id (0, 2) = 0.;
//    Id (1, 0) = 0.;
//    Id (1, 1) = 1., Id (1, 2) = 0.;
//    Id (2, 0) = 0.;
//    Id (2, 1) = 0., Id (2, 2) = 1.;
//
//
//    Real e1x = 1.0;
//    Real e1y = 0.0;
//    Real e1z = 0.0;
//    Real e2x = 0.0;
//    Real e2y = 1.0;
//    Real e2z = 0.0;
//    Real e3x = 0.0;
//    Real e3y = 0.0;
//    Real e3z = 1.0;
//
//    Real sp = nx * e1x + ny * e1y + nz * e1z;
//    Real t1x, t1y, t1z;
//    if (sp != 1)
//    {
//        t1x = e1x - sp * nx;
//        t1y = e1y - sp * ny;
//        t1z = e1z - sp * nz;
//        Real  norm = std::sqrt ( t1x * t1x + t1y * t1y + t1z * t1z );
//        t1x = t1x / norm;
//        t1y = t1y / norm;
//        t1z = t1z / norm;
//    }
//    else
//    {
//        t1x = e2x;
//        t1y = e2y;
//        t1z = e2z;
//    }
//
//    Real t1e1 = t1x * e1x + t1y * e1y + t1z * e1z;
//    Real t1e2 = t1x * e2x + t1y * e2y + t1z * e2z;
//    Real t1e3 = t1x * e3x + t1y * e3y + t1z * e3z;
//
//    Int nLocalDof = M_referencePositionPtr -> epetraVector().MyLength();
//    Int nComponentLocalDof = nLocalDof / 3;
//
//    vectorPtr_Type tmp ( new vector_Type ( *M_lvPositionVectorPtr ) );
//    for (int k (0); k < nComponentLocalDof; k++)
//    {
//        UInt iGID = M_referencePositionPtr-> blockMap().GID (k);
//        UInt jGID = M_referencePositionPtr-> blockMap().GID (k + nComponentLocalDof);
//        UInt kGID = M_referencePositionPtr-> blockMap().GID (k + 2 * nComponentLocalDof);
//
//        (*tmp) [iGID] = (*M_lvPositionVectorPtr) [iGID]
//                        + (* (M_solidPtr -> displacementPtr() ) ) [iGID] * t1e1
//                        + (* (M_solidPtr -> displacementPtr() ) ) [jGID] * t1e2
//                        + (* (M_solidPtr -> displacementPtr() ) ) [kGID] * t1e3;
//    }
//
//
//    VectorSmall<3> E1;
//    E1 (0) = t1x;
//    E1 (1) = t1y;
//    E1 (2) = t1z;
//
//    vectorPtr_Type intergral (new vector_Type ( M_lvPositionVectorPtr -> map() ) );
//
//    {
//        using namespace ExpressionAssembly;
//
//        BOOST_AUTO_TPL (I, value (Id) );
//        BOOST_AUTO_TPL (vE1, value (E1) );
//        BOOST_AUTO_TPL (Grad_u, grad (M_solidPtr -> dispETFESpacePtr(), M_solidPtr -> displacement(), 0) );
//        BOOST_AUTO_TPL (x, value (M_solidPtr -> activeMaterial() -> activationSpace(), *M_lvPositionVectorPtr) );
//        BOOST_AUTO_TPL (F, (Grad_u + I) );
//        BOOST_AUTO_TPL (FmT, minusT (F) );
//        BOOST_AUTO_TPL (J, det (F) );
//        BOOST_AUTO_TPL (x1, dot (x, vE1) );
//
//        QuadratureBoundary myBDQR (buildTetraBDQR (quadRuleTria4pt) );
//
//        *intergral *= 0.0;
//        for (int i (0); i < M_lvFlags.size(); i++)
//        {
//            integrate (boundary (M_solidPtr->mesh(), M_lvFlags[i]), myBDQR, M_solidPtr -> activeMaterial() -> activationSpace(),
//                       value (-1.0) * J * dot (vE1, FmT * Nface) * phi_i) >> intergral;
//        }
//        intergral->globalAssemble();
//
//        M_lvVolume = M_lvPositionVectorPtr->dot (*intergral);
//
//    }
//    M_lvdV = M_lvVolume - Vn;
//}
//
//
//template<typename Mesh, typename IonicModel>
//void EMSolver<Mesh, IonicModel>::computeRVVolume (Real nx, Real ny, Real nz)
//{
//
//    Real Vn = M_rvVolume;
//    MatrixSmall<3, 3> Id;
//    Id (0, 0) = 1.;
//    Id (0, 1) = 0., Id (0, 2) = 0.;
//    Id (1, 0) = 0.;
//    Id (1, 1) = 1., Id (1, 2) = 0.;
//    Id (2, 0) = 0.;
//    Id (2, 1) = 0., Id (2, 2) = 1.;
//
//
//    Real e1x = 1.0;
//    Real e1y = 0.0;
//    Real e1z = 0.0;
//    Real e2x = 0.0;
//    Real e2y = 1.0;
//    Real e2z = 0.0;
//    Real e3x = 0.0;
//    Real e3y = 0.0;
//    Real e3z = 1.0;
//
//    Real sp = nx * e1x + ny * e1y + nz * e1z;
//    Real t1x, t1y, t1z;
//    if (sp != 1)
//    {
//        t1x = e1x - sp * nx;
//        t1y = e1y - sp * ny;
//        t1z = e1z - sp * nz;
//        Real  norm = std::sqrt ( t1x * t1x + t1y * t1y + t1z * t1z );
//        t1x = t1x / norm;
//        t1y = t1y / norm;
//        t1z = t1z / norm;
//    }
//    else
//    {
//        t1x = e2x;
//        t1y = e2y;
//        t1z = e2z;
//    }
//
//    Real t1e1 = t1x * e1x + t1y * e1y + t1z * e1z;
//    Real t1e2 = t1x * e2x + t1y * e2y + t1z * e2z;
//    Real t1e3 = t1x * e3x + t1y * e3y + t1z * e3z;
//
//    Int nLocalDof = M_referencePositionPtr -> epetraVector().MyLength();
//    Int nComponentLocalDof = nLocalDof / 3;
//
//    vectorPtr_Type tmp ( new vector_Type ( *M_rvPositionVectorPtr ) );
//    for (int k (0); k < nComponentLocalDof; k++)
//    {
//        UInt iGID = M_referencePositionPtr-> blockMap().GID (k);
//        UInt jGID = M_referencePositionPtr-> blockMap().GID (k + nComponentLocalDof);
//        UInt kGID = M_referencePositionPtr-> blockMap().GID (k + 2 * nComponentLocalDof);
//
//        (*tmp) [iGID] = (*M_rvPositionVectorPtr) [iGID]
//                        + (* (M_solidPtr -> displacementPtr() ) ) [iGID] * t1e1
//                        + (* (M_solidPtr -> displacementPtr() ) ) [jGID] * t1e2
//                        + (* (M_solidPtr -> displacementPtr() ) ) [kGID] * t1e3;
//    }
//
//
//    VectorSmall<3> E1;
//    E1 (0) = t1x;
//    E1 (1) = t1y;
//    E1 (2) = t1z;
//
//    vectorPtr_Type intergral (new vector_Type ( M_rvPositionVectorPtr -> map() ) );
//
//    {
//        using namespace ExpressionAssembly;
//
//        BOOST_AUTO_TPL (I, value (Id) );
//        BOOST_AUTO_TPL (vE1, value (E1) );
//        BOOST_AUTO_TPL (Grad_u, grad (M_solidPtr -> dispETFESpacePtr(), M_solidPtr -> displacement(), 0) );
//        BOOST_AUTO_TPL (x, value (M_solidPtr -> activeMaterial() -> activationSpace(), *M_rvPositionVectorPtr) );
//        BOOST_AUTO_TPL (F, (Grad_u + I) );
//        BOOST_AUTO_TPL (FmT, minusT (F) );
//        BOOST_AUTO_TPL (J, det (F) );
//        BOOST_AUTO_TPL (x1, dot (x, vE1) );
//
//        QuadratureBoundary myBDQR (buildTetraBDQR (quadRuleTria4pt) );
//
//        *intergral *= 0.0;
//        for (int i (0); i < M_rvFlags.size(); i++)
//        {
//            integrate (boundary (M_solidPtr->mesh(), M_rvFlags[i]), myBDQR, M_solidPtr -> activeMaterial() -> activationSpace(),
//                       value (-1.0) * J * dot (vE1, FmT * Nface) * phi_i) >> intergral;
//        }
//        intergral->globalAssemble();
//
//        M_rvVolume = M_rvPositionVectorPtr->dot (*intergral);
//
//    }
//    M_rvdV = M_rvVolume - Vn;
//}
//
//
//
//
////  }
//

} // namespace LifeV

#endif //_MONODOMAINSOLVER_H_
