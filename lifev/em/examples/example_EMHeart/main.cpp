//============================================//
// Includes
//============================================//

#include <lifev/core/LifeV.hpp>

// Passive material
#include <lifev/electrophysiology/solver/ElectroETAMonodomainSolver.hpp>
#include <lifev/electrophysiology/util/HeartUtility.hpp>
#include <lifev/structure/solver/StructuralConstitutiveLawData.hpp>
#include <lifev/em/solver/mechanics/EMStructuralOperator.hpp>
#include <lifev/em/solver/mechanics/EMStructuralConstitutiveLaw.hpp>
#include <lifev/bc_interface/3D/bc/BCInterface3D.hpp>

// Solver
#include <lifev/em/solver/EMSolver.hpp>

// Exporter
#include <lifev/core/filter/ExporterEnsight.hpp>
#include <lifev/core/filter/ExporterVTK.hpp>
#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5.hpp>
#endif
#include <lifev/core/filter/ExporterEmpty.hpp>

// Resize mesh
#include <lifev/core/mesh/MeshUtility.hpp>

// B.C. modification
#include <lifev/core/fem/BCVector.hpp>

// Circulation
#include <lifev/em/solver/circulation/Circulation.hpp>

// Volume computation
#include <lifev/em/solver/circulation/CirculationVolumeIntegrator.hpp>

// Heart solver
#include <lifev/em/solver/HeartSolver.hpp>

// PatchBC
#include <lifev/em/examples/example_EMHeart/PatchBC.hpp>

// Track nan
// #include <fenv.h>



using namespace LifeV;


int main (int argc, char** argv)
{

//    feenableexcept(FE_INVALID | FE_OVERFLOW);
    
    //============================================//
    // Typedefs
    //============================================//
    
    typedef RegionMesh<LinearTetra>                         mesh_Type;
    typedef boost::shared_ptr<mesh_Type>                    meshPtr_Type;
    
    typedef boost::function < Real (const Real & t,
                                    const Real &   x,
                                    const Real &   y,
                                    const Real & z,
                                    const ID&   /*i*/ ) >   function_Type;
    
    typedef VectorEpetra                                    vector_Type;
    typedef boost::shared_ptr<vector_Type>                  vectorPtr_Type;
    
    typedef BCVector                                        bcVector_Type;
    typedef boost::shared_ptr<bcVector_Type>                bcVectorPtr_Type;
    
    typedef EMMonodomainSolver<mesh_Type>                   monodomain_Type;

    typedef FESpace< mesh_Type, MapEpetra >                 solidFESpace_Type;
    typedef boost::shared_ptr<solidFESpace_Type>            solidFESpacePtr_Type;
    
    typedef ETFESpace< mesh_Type, MapEpetra, 3, 3 >         solidETFESpace_Type;
    typedef boost::shared_ptr<solidETFESpace_Type>          solidETFESpacePtr_Type;
    
    typedef StructuralConstitutiveLawData                   constitutiveLaw_Type;
    typedef boost::shared_ptr<constitutiveLaw_Type>         constitutiveLawPtr_Type;
    
    typedef BCHandler                                       bc_Type;
    typedef StructuralOperator< mesh_Type >                 physicalSolver_Type;
    typedef BCInterface3D< bc_Type, physicalSolver_Type >   bcInterface_Type;
    typedef boost::shared_ptr< bcInterface_Type >           bcInterfacePtr_Type;
    

    //============================================//
    // Communicator and displayer
    //============================================//
#ifdef HAVE_MPI
    MPI_Init ( &argc, &argv );
#endif

    boost::shared_ptr<Epetra_Comm>  comm ( new Epetra_MpiComm (MPI_COMM_WORLD) );
    Displayer displayer ( comm );
    
    displayer.leaderPrint("\n\nEMHeart running ...\n\n");
    
    //============================================//
    // Read data file and create output folder
    //============================================//
    GetPot command_line (argc, argv);
    const std::string data_file_name = command_line.follow ("data", 2, "-f", "--file");
    GetPot dataFile (data_file_name);
    std::string problemFolder = EMUtility::createOutputFolder (command_line, *comm);

    
    //============================================//
    // Electromechanic solver
    //============================================//
    EMSolver<mesh_Type, monodomain_Type> solver(comm);
    
    
    //============================================//
    // Body circulation
    //============================================//
    const std::string circulationInputFile = command_line.follow ("circulation", 2, "-cif", "--cifile");
    const std::string circulationOutputFile = command_line.follow ( (problemFolder + "solution.dat").c_str(), 2, "-cof", "--cofile");
    
    Circulation circulationSolver( circulationInputFile );
    
    // Flow rate between two vertices
    auto Q = [&circulationSolver] (const std::string& N1, const std::string& N2) { return circulationSolver.solution ( std::vector<std::string> {N1, N2} ); };
    auto p = [&circulationSolver] (const std::string& N1) { return circulationSolver.solution ( N1 ); };
    
    
    //============================================//
    // Heart solver
    //============================================//
    HeartSolver<EMSolver<mesh_Type, monodomain_Type> > heartSolver (solver, circulationSolver);
    
    
    //============================================//
    // Setup material data
    //============================================//
    EMData emdata;
    emdata.setup (dataFile);
    
    
    //============================================//
    // Load mesh
    //============================================//
    std::string meshName = dataFile("solid/space_discretization/mesh_file", "cube4.mesh");
    std::string meshPath = dataFile("solid/space_discretization/mesh_dir", "./");
    
    solver.loadMesh (meshName, meshPath);
    
    
    //============================================//
    // Resize mesh
    //============================================//
    if ( 0 == comm->MyPID() ) std::cout << "\nResizing mesh ... " << '\r' << std::flush;

    std::vector<Real> scale (3, dataFile("solid/space_discretization/mesh_scaling", 1.0));
    std::vector<Real> rotate { dataFile("solid/space_discretization/mesh_rotation_0", 0.0) , dataFile("solid/space_discretization/mesh_rotation_1", 0.0) , dataFile("solid/space_discretization/mesh_rotation_2", 0.0) };
    std::vector<Real> translate { dataFile("solid/space_discretization/mesh_translation_0", 0.0) , dataFile("solid/space_discretization/mesh_translation_1", 0.0) , dataFile("solid/space_discretization/mesh_translation_2", 0.0) };
    
    MeshUtility::MeshTransformer<mesh_Type> transformerFull (* (solver.fullMeshPtr() ) );
    MeshUtility::MeshTransformer<mesh_Type> transformerLocal (* (solver.localMeshPtr() ) );
    
    transformerFull.transformMesh (scale, rotate, translate);
    transformerLocal.transformMesh (scale, rotate, translate);
    
    if ( 0 == comm->MyPID() ) std::cout << "Resizing mesh done" << '\r' << std::flush;
    if ( 0 == comm->MyPID() ) solver.fullMeshPtr()->showMe();

    
    //============================================//
    // Setup solver (including fe-spaces & b.c.)
    //============================================//
    EMAssembler::quadRule.setQuadRule( dataFile ( "solid/space_discretization/quad_rule", "4pt") );
    solver.setup (dataFile);
    
    auto& disp = solver.structuralOperatorPtr() -> displacement();
    auto FESpace = solver.structuralOperatorPtr() -> dispFESpacePtr();
    auto dETFESpace = solver.electroSolverPtr() -> displacementETFESpacePtr();
    auto ETFESpace = solver.electroSolverPtr() -> ETFESpacePtr();
    
    
    //============================================//
    // Setup anisotropy vectors
    //============================================//
    bool anisotropy = dataFile ( "solid/space_discretization/anisotropic", false );

    if ( anisotropy )
    {
        std::string fiberFileName  =  dataFile ( "solid/space_discretization/fiber_name", "FiberDirection");
        std::string sheetFileName  =  dataFile ( "solid/space_discretization/sheet_name", "SheetsDirection");
        std::string fiberFieldName =  dataFile ( "solid/space_discretization/fiber_fieldname", "fibers");
        std::string sheetFieldName =  dataFile ( "solid/space_discretization/sheet_fieldname", "sheets");
        std::string fiberDir       =  meshPath; //dataFile ( "solid/space_discretization/fiber_dir", "./");
        std::string sheetDir       =  meshPath; //dataFile ( "solid/space_discretization/sheet_dir", "./");
        std::string elementOrder   =  dataFile ( "solid/space_discretization/order", "P1");

        solver.setupFiberVector ( fiberFileName, fiberFieldName, fiberDir, elementOrder );
        solver.setupMechanicalSheetVector ( sheetFileName, sheetFieldName, sheetDir, elementOrder );
    }
    else
    {
        solver.setupFiberVector (1., 0., 0.);
        solver.setupSheetVector (0., 1., 0.);
    }
    
    
    //============================================//
    // Initialize electrophysiology
    //============================================//
    solver.initialize();
    
    
    //============================================//
    // Building Matrices
    //============================================//
    solver.oneWayCoupling();
    solver.structuralOperatorPtr()->setNewtonParameters(dataFile);
    solver.buildSystem();
    
    if ( 0 == comm->MyPID() )
    {
        std::cout << "\nNode number: " << disp.epetraVector().MyLength() << std::endl;
    }

    
    //============================================//
    // Setup exporters for EMSolver
    //============================================//
    solver.setupExporters (problemFolder);
        
    
    //============================================//
    // Electric stimulus function
    //============================================//
    function_Type stim = &HeartSolver<EMSolver<mesh_Type, monodomain_Type> >::Iapp;

    
    //============================================//
    // Create patches for essential patch b.c.
    //============================================//

    UInt nDispPatchBC = dataFile.vector_variable_size ( ( "solid/boundary_conditions/listEssentialPatchBC" ) );
    
    for ( UInt i (0) ; i < nDispPatchBC ; ++i )
    {
        std::string patchName = dataFile ( ( "solid/boundary_conditions/listEssentialPatchBC" ), " ", i );
        Real patchFlag = dataFile ( ("solid/boundary_conditions/" + patchName + "/flag").c_str(), 0 );
        Real patchRadius = dataFile ( ("solid/boundary_conditions/" + patchName + "/radius").c_str(), 1.0 );
        
        Vector3D patchCenter;
        for ( UInt j (0); j < 3; ++j )
        {
            patchCenter[j] = dataFile ( ("solid/boundary_conditions/" + patchName + "/center").c_str(), 0, j );
        }

        heartSolver.createPatch(solver, patchCenter, patchRadius, patchFlag, (900+i));
    }
    
    
    //============================================//
    // Create patches for natural patch b.c.
    //============================================//
    
    UInt nForcePatchBC = dataFile.vector_variable_size ( ( "solid/boundary_conditions/listNaturalPatchBC" ) );
    
    for ( UInt i (0) ; i < nForcePatchBC ; ++i )
    {
        std::string patchName = dataFile ( ( "solid/boundary_conditions/listNaturalPatchBC" ), " ", i );
        Real patchFlag = dataFile ( ("solid/boundary_conditions/" + patchName + "/flag").c_str(), 0 );
        Real patchRadius = dataFile ( ("solid/boundary_conditions/" + patchName + "/radius").c_str(), 1.0 );
        
        Vector3D patchCenter;
        for ( UInt j (0); j < 3; ++j )
        {
            patchCenter[j] = dataFile ( ("solid/boundary_conditions/" + patchName + "/center").c_str(), 0, j );
        }
        
        heartSolver.createPatch(solver, patchCenter, patchRadius, patchFlag, (800+i));
    }
    
    
    //============================================//
    // Create displacement patch b.c.
    //============================================//
    std::vector<Real> patchDisplacement;
    std::vector<Vector3D> patchDirection;
    
    std::vector<vectorPtr_Type> patchDispVecPtr;
    std::vector<bcVectorPtr_Type> patchDispBCVecPtr;
    
    const VectorEpetra dispPreload (disp);
    
    if ( 0 == comm->MyPID() )
    {
        std::cout << "\n*****************************************************************";
        std::cout << "\nCreate Essential Patch BC";
        std::cout << "\n*****************************************************************\n";
    }
    
    for ( UInt i (0) ; i < nDispPatchBC ; ++i )
    {
        std::string patchName = dataFile ( ( "solid/boundary_conditions/listPatchBC" ), " ", i );
        Real patchFlag = dataFile ( ("solid/boundary_conditions/" + patchName + "/flag").c_str(), 0 );
        Real patchRadius = dataFile ( ("solid/boundary_conditions/" + patchName + "/radius").c_str(), 1.0 );
        patchDisplacement.push_back( dataFile ( ("solid/boundary_conditions/" + patchName + "/displacement").c_str(), 1.0 ) );
        
        Vector3D patchCenter;
        for ( UInt j (0); j < 3; ++j )
        {
            patchCenter[j] = dataFile ( ("solid/boundary_conditions/" + patchName + "/center").c_str(), 0, j );
        }
        
        Vector3D pd;
        for ( UInt j (0); j < 3; ++j )
        {
            pd[j] = dataFile ( ("solid/boundary_conditions/" + patchName + "/direction").c_str(), 0, j );
        }
        patchDirection.push_back(pd);
        
        UInt componentSize = dataFile.vector_variable_size ( ("solid/boundary_conditions/" + patchName + "/component").c_str() );
        std::vector<ID> patchComponent (componentSize);
        for ( UInt j (0); j < componentSize; ++j )
        {
            patchComponent[j] = dataFile ( ("solid/boundary_conditions/" + patchName + "/component").c_str(), 0, j );
        }
        
        patchDispVecPtr.push_back ( heartSolver.directionalVectorField(FESpace, patchDirection[i], 1e-10) );
        *patchDispVecPtr[i] += dispPreload;
        patchDispBCVecPtr.push_back ( bcVectorPtr_Type( new bcVector_Type( *patchDispVecPtr[i], solver.structuralOperatorPtr() -> dispFESpacePtr() -> dof().numTotalDof(), 1 ) ) );
        solver.bcInterfacePtr() -> handler()->addBC (patchName, (900+i),  Essential, Component, *patchDispBCVecPtr[i], patchComponent);
    }
    
    if ( 0 == comm->MyPID() )
    {
        std::cout << "\n*****************************************************************";
        std::cout << "\nCreate Essential Patch BC done";
        std::cout << "\n*****************************************************************\n";
    }
    
    Real tmax = dataFile ( "solid/patches/tmax", 0. );
    Real tduration = dataFile ( "solid/patches/tduration", 0. );
    
    auto modifyEssentialPatchBC = [&] (const Real& time)
    {
        for ( UInt i (0) ; i < nDispPatchBC ; ++i )
        {
            Real currentPatchDisp = heartSolver.sinSquared(time, patchDisplacement[i], tmax, tduration);
            
            patchDispVecPtr[i] = heartSolver.directionalVectorField(FESpace, patchDirection[i], currentPatchDisp);
            *patchDispVecPtr[i] += dispPreload;
            
            patchDispBCVecPtr[i].reset( new bcVector_Type( *patchDispVecPtr[i], FESpace->dof().numTotalDof(), 1 ) );
            solver.bcInterfacePtr()->handler()->modifyBC((900+i), *patchDispBCVecPtr[i]);
        }
    };
    
    if ( 0 == comm->MyPID() )
    {
        std::cout << "\n*****************************************************************";
        std::cout << "\nCreate Essential Patch BC update done";
        std::cout << "\n*****************************************************************\n";
    }
    
    //============================================//
    // Create force patch b.c.
    //============================================//
    std::vector<Real> patchForce;
    std::vector<Vector3D> patchForceDirection;
    
    std::vector<vectorPtr_Type> patchForceVecPtr;
    std::vector<bcVectorPtr_Type> patchForceBCVecPtr;
    
    if ( 0 == comm->MyPID() )
    {
        std::cout << "\n*****************************************************************";
        std::cout << "\nCreate Natural Patch BC";
        std::cout << "\n*****************************************************************\n";
    }
    
    for ( UInt i (0) ; i < nForcePatchBC ; ++i )
    {
        std::string patchName = dataFile ( ( "solid/boundary_conditions/listNaturalPatchBC" ), " ", i );
        Real patchFlag = dataFile ( ("solid/boundary_conditions/" + patchName + "/flag").c_str(), 0 );
        Real patchRadius = dataFile ( ("solid/boundary_conditions/" + patchName + "/radius").c_str(), 1.0 );
        patchForce.push_back( dataFile ( ("solid/boundary_conditions/" + patchName + "/force").c_str(), 1.0 ) );
        
        Vector3D patchCenter;
        for ( UInt j (0); j < 3; ++j )
        {
            patchCenter[j] = dataFile ( ("solid/boundary_conditions/" + patchName + "/center").c_str(), 0, j );
        }
        
        Vector3D pd;
        for ( UInt j (0); j < 3; ++j )
        {
            pd[j] = dataFile ( ("solid/boundary_conditions/" + patchName + "/direction").c_str(), 0, j );
        }
        patchForceDirection.push_back(pd);
        
        UInt componentSize = dataFile.vector_variable_size ( ("solid/boundary_conditions/" + patchName + "/component").c_str() );
        std::vector<ID> patchComponent (componentSize);
        for ( UInt j (0); j < componentSize; ++j )
        {
            patchComponent[j] = dataFile ( ("solid/boundary_conditions/" + patchName + "/component").c_str(), 0, j );
        }
        
        patchForceVecPtr.push_back ( heartSolver.directionalVectorField(FESpace, patchForceDirection[i], 1e-10) );
        patchForceBCVecPtr.push_back ( bcVectorPtr_Type( new bcVector_Type( *patchForceVecPtr[i], solver.structuralOperatorPtr() -> dispFESpacePtr() -> dof().numTotalDof(), 1 ) ) );
        solver.bcInterfacePtr() -> handler()->addBC (patchName, (800+i), Natural, Component, *patchForceBCVecPtr[i], patchComponent);
    }
    
    if ( 0 == comm->MyPID() )
    {
        std::cout << "\n*****************************************************************";
        std::cout << "\nCreate Natural Patch BC done";
        std::cout << "\n*****************************************************************\n";
    }
    
    auto modifyNaturalPatchBC = [&] (const Real& time)
    {
        for ( UInt i (0) ; i < nForcePatchBC ; ++i )
        {
            Real currentPatchForce = heartSolver.sinSquared(time, patchForce[i], tmax, tduration);
            patchForceVecPtr[i] = heartSolver.directionalVectorField(FESpace, patchForceDirection[i], currentPatchForce);
            
            patchForceBCVecPtr[i].reset( new bcVector_Type( *patchForceVecPtr[i], FESpace->dof().numTotalDof(), 1 ) );
            solver.bcInterfacePtr()->handler()->modifyBC((800+i), *patchForceBCVecPtr[i]);
        }
    };
    
    if ( 0 == comm->MyPID() )
    {
        std::cout << "\n*****************************************************************";
        std::cout << "\nCreate Natural Patch BC update done";
        std::cout << "\n*****************************************************************\n";
    }
    
    solver.bcInterfacePtr() -> handler() -> bcUpdate( *solver.structuralOperatorPtr() -> dispFESpacePtr() -> mesh(), solver.structuralOperatorPtr() -> dispFESpacePtr() -> feBd(), solver.structuralOperatorPtr() -> dispFESpacePtr() -> dof() );
    
    if ( 0 == comm->MyPID() ) solver.bcInterfacePtr() -> handler() -> showMe();
    
    
    //============================================//
    // Pressure b.c. on endocardia
    //============================================//
    std::vector<vectorPtr_Type> pVecPtrs;
    std::vector<bcVectorPtr_Type> pBCVecPtrs;

    std::vector<ID> flagsBC;
    std::vector<UInt> ventIdx;

    // Endocardia
    UInt nVarBC = dataFile.vector_variable_size ( ( "solid/boundary_conditions/listVariableBC" ) );
    for ( UInt i (0) ; i < nVarBC ; ++i )
    {
        std::string varBCSection = dataFile ( ( "solid/boundary_conditions/listVariableBC" ), " ", i );
        flagsBC.push_back( dataFile ( ("solid/boundary_conditions/" + varBCSection + "/flag").c_str(), 0 ) );
        ventIdx.push_back ( dataFile ( ("solid/boundary_conditions/" + varBCSection + "/index").c_str(), 0 ) );
        
        pVecPtrs.push_back ( vectorPtr_Type ( new vector_Type ( solver.structuralOperatorPtr() -> displacement().map(), Repeated ) ) );
        pBCVecPtrs.push_back ( bcVectorPtr_Type( new bcVector_Type( *pVecPtrs[i], solver.structuralOperatorPtr() -> dispFESpacePtr() -> dof().numTotalDof(), 1 ) ) );
        solver.bcInterfacePtr() -> handler() -> addBC(varBCSection, flagsBC[i], Natural, Full, *pBCVecPtrs[i], 3);
    }
    

    // Functions to modify b.c.
    auto modifyPressureBC = [&] (const std::vector<Real>& bcValues)
    {
        for ( UInt i (0) ; i < nVarBC ; ++i )
        {
            *pVecPtrs[i] = - bcValues[ ventIdx[i] ] * 1333.224;
            pBCVecPtrs[i].reset ( ( new bcVector_Type (*pVecPtrs[i], solver.structuralOperatorPtr() -> dispFESpacePtr() -> dof().numTotalDof(), 1) ) );
            solver.bcInterfacePtr() -> handler() -> modifyBC(flagsBC[i], *pBCVecPtrs[i]);
        }
    };

    
    //============================================//
    // Update and print bcHandler
    //============================================//
    solver.bcInterfacePtr() -> handler() -> bcUpdate( *solver.structuralOperatorPtr() -> dispFESpacePtr() -> mesh(), solver.structuralOperatorPtr() -> dispFESpacePtr() -> feBd(), solver.structuralOperatorPtr() -> dispFESpacePtr() -> dof() );
    
    if ( 0 == comm->MyPID() ) solver.bcInterfacePtr() -> handler() -> showMe();
    
    
    //============================================//
    // Volume integrators
    //============================================//
    std::vector<int> LVFlags;
    std::vector<int> RVFlags;

    for ( UInt i (0) ; i < nVarBC ; ++i )
    {
        std::string varBCSection = dataFile ( ( "solid/boundary_conditions/listVariableBC" ), " ", i );
        ID flag  =  dataFile ( ("solid/boundary_conditions/" + varBCSection + "/flag").c_str(), 0 );
        ID index =  dataFile ( ("solid/boundary_conditions/" + varBCSection + "/index").c_str(), 0 );
        switch ( index )
        {
            case 0: LVFlags.push_back ( flag );
                break;
            case 1: RVFlags.push_back ( flag );
                break;
            default:
                break;
        }
    }
    
    VolumeIntegrator LV (LVFlags, "Left Ventricle", solver.fullMeshPtr(), solver.localMeshPtr(), ETFESpace, FESpace);
    VolumeIntegrator RV (RVFlags, "Right Ventricle", solver.fullMeshPtr(), solver.localMeshPtr(), ETFESpace, FESpace);

    
    //============================================//
    // Set variables and functions
    //============================================//
    const Real dt_activation = solver.data().electroParameter<Real>("timestep");
    const Real dt_loadstep =  dataFile ( "solid/time_discretization/dt_loadstep", 1.0 );
    const Real activationLimit_loadstep =  dataFile ( "solid/time_discretization/activation_limit_loadstep", 0.0 );
    const Real dt_mechanics = solver.data().solidParameter<Real>("timestep");
    const Real dt_save = dataFile ( "exporter/save", 10. );
    const Real endtime = solver.data().electroParameter<Real>("endtime");
    const UInt mechanicsLoadstepIter = static_cast<UInt>( dt_loadstep / dt_activation );
    const UInt mechanicsCouplingIter = static_cast<UInt>( dt_mechanics / dt_activation );
    const UInt maxiter = static_cast<UInt>( endtime / dt_activation ) ;
    
    const Real pPerturbationFe = dataFile ( "solid/coupling/pPerturbationFe", 1e-2 );
    const Real pPerturbationCirc = dataFile ( "solid/coupling/pPerturbationCirc", 1e-3 );
    const Real couplingError = dataFile ( "solid/coupling/couplingError", 1e-6 );
    const UInt couplingJFeSubIter = dataFile ( "solid/coupling/couplingJFeSubIter", 1 );
    const UInt couplingJFeSubStart = dataFile ( "solid/coupling/couplingJFeSubStart", 1 );
    const UInt couplingJFeIter = dataFile ( "solid/coupling/couplingJFeIter", 1 );
    
    const Real dpMax = dataFile ( "solid/coupling/dpMax", 0.1 );
    
    std::vector<std::vector<std::string> > bcNames { { "lv" , "p" } , { "rv" , "p" } };
    std::vector<double> bcValues { p ( "lv" ) , p ( "rv") };
    std::vector<double> bcValuesPre ( bcValues );

    VectorSmall<2> VCirc, VCircNew, VCircPert, VFe, VFeNew, VFePert, R, dp;
    MatrixSmall<2,2> JFe, JCirc, JR;

    UInt iter (0);
    Real t (0);
    
    auto printCoupling = [&] ( std::string label ) { if ( 0 == comm->MyPID() )
    {
        std::cout << "\n===============================================================";
        std::cout << "\nCoupling: " << label;
        std::cout << "\nNewton iteration nr. " << iter << " at time " << t;
        std::cout << "\nLV - Pressure: \t\t\t" << bcValues[0];
        std::cout << "\nLV - FE-Volume: \t\t" << VFeNew[0];
        std::cout << "\nLV - Circulation-Volume: \t" << VCircNew[0];
        std::cout << "\nLV - Residual: \t\t\t" << std::abs(VFeNew[0] - VCircNew[0]);
        std::cout << "\nRV - Pressure: \t\t\t" << bcValues[1];
        std::cout << "\nRV - FE-Volume : \t\t" << VFeNew[1];
        std::cout << "\nRV - Circulation-Volume: \t" << VCircNew[1];
        std::cout << "\nRV - Residual: \t\t\t" << std::abs(VFeNew[1] - VCircNew[1]);
        //std::cout << "\nJFe   = " << JFe;
        //std::cout << "\nJCirc = " << JCirc;
        //std::cout << "\nJR    = " << JR;
        std::cout << "\n===============================================================\n"; }
    };
    
    
    //============================================//
    // Load restart file
    //============================================//
    
    std::string restartInput = command_line.follow ("noRestart", 2, "-r", "--restart");
    const bool restart ( restartInput != "noRestart" );

    if ( restart )
    {
        const std::string restartDir = command_line.follow (problemFolder.c_str(), 2, "-rd", "--restartDir");
        
        Real dtExport = 10.;
        
        // Set time variable
        const unsigned int restartInputStr = std::stoi(restartInput);
        const unsigned int nIter = (restartInputStr - 1) * dtExport / dt_mechanics;
        t = nIter * dt_mechanics;

        // Set time exporter time index
        solver.setTimeIndex(restartInputStr + 1);

        // Load restart solutions from output files
        std::string polynomialDegree = dataFile ( "solid/space_discretization/order", "P2");

        ElectrophysiologyUtility::importVectorField ( solver.structuralOperatorPtr() -> displacementPtr(), "MechanicalSolution" , "displacement", solver.localMeshPtr(), restartDir, polynomialDegree, restartInput );

        for ( unsigned int i = 0; i < solver.electroSolverPtr()->globalSolution().size() ; ++i )
        {
            ElectrophysiologyUtility::importScalarField (solver.electroSolverPtr()->globalSolution().at(i), "ElectroSolution" , ("Variable" + std::to_string(i)), solver.localMeshPtr(), restartDir, polynomialDegree, restartInput );
        }

        ElectrophysiologyUtility::importScalarField (solver.activationModelPtr() -> fiberActivationPtr(), "ActivationSolution" , "Activation", solver.localMeshPtr(), restartDir, polynomialDegree, restartInput );
        //ElectrophysiologyUtility::importScalarField (solver.activationTimePtr(), "ActivationTimeSolution" , "Activation Time", solver.localMeshPtr(), restartDir, polynomialDegree, restartInput );

        circulationSolver.restartFromFile ( restartDir + "solution.dat" , nIter );
        
        // Set boundary mechanics conditions
        bcValues = { p ( "lv" ) , p ( "rv" ) };
        bcValuesPre = { p ( "lv" ) , p ( "rv" ) };
        modifyPressureBC(bcValues);
    }

    
    //============================================//
    // Preload
    //============================================//
    
    if ( ! restart )
    {
        solver.structuralOperatorPtr() -> data() -> dataTime() -> setTime(0.0);
        
        const int preloadSteps = dataFile ( "solid/boundary_conditions/numPreloadSteps", 0);
        
        auto preloadPressure = [] (std::vector<double> p, const int& step, const int& steps)
        {
            for (auto& i : p) {i *= double(step) / double(steps);}
            return p;
        };
        
        LifeChrono chronoSave;
        chronoSave.start();

        solver.saveSolution (-1.0);
        
        if ( 0 == comm->MyPID() )
        {
            std::cout << "\n*****************************************************************";
            std::cout << "\nData stored in " << chronoSave.diff() << " s";
            std::cout << "\n*****************************************************************\n";
        }
        
        LifeChrono chronoPreload;
        chronoPreload.start();
        
        for (int i (1); i <= preloadSteps; i++)
        {
            if ( 0 == comm->MyPID() )
            {
                std::cout << "\n*****************************************************************";
                std::cout << "\nPreload step: " << i << " / " << preloadSteps;
                std::cout << "\n*****************************************************************\n";
            }

            // Update pressure b.c.
            modifyPressureBC(preloadPressure(bcValues, i, preloadSteps));

            // Solve mechanics
            solver.bcInterfacePtr() -> updatePhysicalSolverVariables();
            solver.solveMechanics();
            // solver.saveSolution (i-1);
        }

        auto maxI4fValue ( solver.activationModelPtr()->I4f().maxValue() );
        auto minI4fValue ( solver.activationModelPtr()->I4f().minValue() );
        
        if ( 0 == comm->MyPID() ) std::cout << "\nI4fmax/I4fmin = " << maxI4fValue << "/" << minI4fValue << std::endl;

        if ( 0 == comm->MyPID() )
        {
            std::cout << "\n*****************************************************************";
            std::cout << "\nPreload done in: " << chronoPreload.diff();
            std::cout << "\n*****************************************************************\n";
        }

    }


    //============================================//
    // Time loop
    //============================================//
    
    VFe[0] = LV.volume(disp, dETFESpace, - 1);
    VFe[1] = RV.volume(disp, dETFESpace, 1);
    VCirc = VFe;
    
    VectorEpetra dispCurrent ( disp );
    ID bdPowerFlag  =  dataFile ( ("solid/boundary_conditions/LVEndo/flag") , 0 );
    
    printCoupling("Initial values");
    
    auto perturbedPressure = [] (std::vector<double> p, const double& dp)
    {
        for (auto& i : p) {i += dp;}
        return p;
    };
    
    auto perturbedPressureComp = [] (std::vector<double> p, const double& dp, int comp)
    {
        p[comp] += dp;
        return p;
    };

    if ( ! restart )
    {
        solver.saveSolution (t);
        circulationSolver.exportSolution( circulationOutputFile );
    }

    for (int k (1); k <= maxiter; k++)
    {
        if ( 0 == comm->MyPID() )
        {
            std::cout << "\n*****************************************************************";
            std::cout << "\nTIME = " << t+dt_activation;
            std::cout << "\n*****************************************************************\n";
        }

        t = t + dt_activation;

        //============================================//
        // Solve electrophysiology and activation
        //============================================//

        auto maxI4fValue ( solver.activationModelPtr()->I4f().maxValue() );
        auto minI4fValue ( solver.activationModelPtr()->I4f().minValue() );
        
        if ( 0 == comm->MyPID() ) std::cout << "\nI4fmax/I4fmin = " << maxI4fValue << "/" << minI4fValue << std::endl;
        
        solver.solveElectrophysiology (stim, t);
        solver.solveActivation (dt_activation);

        //============================================//
        // Load steps mechanics (activation & b.c.)
        //============================================//

        auto minActivationValue ( solver.activationModelPtr() -> fiberActivationPtr() -> minValue() );

        if ( k % mechanicsLoadstepIter == 0 && k % mechanicsCouplingIter != 0 && minActivationValue < activationLimit_loadstep )
        {
            if ( 0 == comm->MyPID() )
            {
                std::cout << "\n*****************************************************************";
                std::cout << "\nLoad step at time = " << t;
                std::cout << "\nMinimal activation value = " << minActivationValue;
                std::cout << "\n*****************************************************************\n";
            }

            // Linear b.c. extrapolation
            auto bcValuesLoadstep ( bcValues );
            bcValuesLoadstep[0] = bcValues[0] + ( bcValues[0] - bcValuesPre[0] ) * ( k % mechanicsCouplingIter ) / mechanicsCouplingIter;
            bcValuesLoadstep[1] = bcValues[1] + ( bcValues[1] - bcValuesPre[1] ) * ( k % mechanicsCouplingIter ) / mechanicsCouplingIter;

            if ( 0 == comm->MyPID() )
            {
                std::cout << "\n***************************************************************";
                std::cout << "\nLV-Pressure extrapolation from " <<  bcValues[0] << " to " <<  bcValuesLoadstep[0];
                std::cout << "\nRV-Pressure extrapolation from " <<  bcValues[1] << " to " <<  bcValuesLoadstep[1];
                std::cout << "\n***************************************************************\n\n";
            }

            // Load step mechanics
            solver.structuralOperatorPtr() -> data() -> dataTime() -> setTime(t);
            modifyPressureBC(bcValuesLoadstep);
            modifyEssentialPatchBC(t);
            modifyNaturalPatchBC(t);
            solver.bcInterfacePtr() -> updatePhysicalSolverVariables();
            solver.solveMechanics();
        }
        
        
        //============================================//
        // Iterate mechanics / circulation
        //============================================//
        
        if ( k % mechanicsCouplingIter == 0 )
        {
            iter = 0;
            const double dt_circulation ( dt_mechanics / 1000 );
            solver.structuralOperatorPtr() -> data() -> dataTime() -> setTime(t);
            
            modifyEssentialPatchBC(t);
            modifyNaturalPatchBC(t);
            
            //============================================//
            // 4th order Adam-Bashforth pressure extrapol.
            //============================================//
            
            heartSolver.extrapolate4thOrderAdamBashforth(bcValues, bcValuesPre, dpMax);
            
            if ( 0 == comm->MyPID() )
            {
                std::cout << "\n***************************************************************";
                std::cout << "\nMinimal activation value = " << minActivationValue;
                std::cout << "\nLV-Pressure extrapolation from " <<  bcValuesPre[0] << " to " <<  bcValues[0];
                std::cout << "\nRV-Pressure extrapolation from " <<  bcValuesPre[1] << " to " <<  bcValues[1];
                std::cout << "\n***************************************************************\n\n";
            }

            
            //============================================//
            // Solve mechanics
            //============================================//
            modifyPressureBC(bcValues);
            solver.bcInterfacePtr() -> updatePhysicalSolverVariables();
            solver.solveMechanics();
            
            VFeNew[0] = LV.volume(disp, dETFESpace, - 1);
            VFeNew[1] = RV.volume(disp, dETFESpace, 1);

            //============================================//
            // Solve circlation
            //============================================//
            circulationSolver.iterate(dt_circulation, bcNames, bcValues, iter);
            VCircNew[0] = VCirc[0] + dt_circulation * ( Q("la", "lv") - Q("lv", "sa") );
            VCircNew[1] = VCirc[1] + dt_circulation * ( Q("ra", "rv") - Q("rv", "pa") );

            //============================================//
            // Residual computation
            //============================================//
            R = VFeNew - VCircNew;
            printCoupling("Residual Computation");

            //============================================//
            // Newton iterations
            //============================================//
            while ( R.norm() > couplingError )
            {
                ++iter;

                //============================================//
                // Jacobian circulation
                //============================================//

                // Left ventricle
                circulationSolver.iterate(dt_circulation, bcNames, perturbedPressureComp(bcValues, pPerturbationCirc, 0), iter);
                VCircPert[0] = VCirc[0] + dt_circulation * ( Q("la", "lv") - Q("lv", "sa") );
                VCircPert[1] = VCirc[1] + dt_circulation * ( Q("ra", "rv") - Q("rv", "pa") );

                JCirc(0,0) = ( VCircPert[0] - VCircNew[0] ) / pPerturbationCirc;
                JCirc(1,0) = ( VCircPert[1] - VCircNew[1] ) / pPerturbationCirc;

                // Right ventricle
                circulationSolver.iterate(dt_circulation, bcNames, perturbedPressureComp(bcValues, pPerturbationCirc, 1), iter);
                VCircPert[0] = VCirc[0] + dt_circulation * ( Q("la", "lv") - Q("lv", "sa") );
                VCircPert[1] = VCirc[1] + dt_circulation * ( Q("ra", "rv") - Q("rv", "pa") );

                JCirc(0,1) = ( VCircPert[0] - VCircNew[0] ) / pPerturbationCirc;
                JCirc(1,1) = ( VCircPert[1] - VCircNew[1] ) / pPerturbationCirc;


                //============================================//
                // Jacobian fe
                //============================================//

                const bool jFeIter ( ! ( k % (couplingJFeIter * mechanicsCouplingIter) ) );
                const bool jFeSubIter ( ! ( (iter - couplingJFeSubStart) % couplingJFeSubIter) && iter >= couplingJFeSubStart );
                const bool jFeEmpty ( JFe.norm() == 0 );

                if ( jFeIter || jFeSubIter || jFeEmpty )
                {
                    JFe *= 0.0;
                    dispCurrent = disp;

                    // Left ventricle
                    modifyPressureBC(perturbedPressureComp(bcValues, pPerturbationFe, 0));
                    solver.bcInterfacePtr() -> updatePhysicalSolverVariables();
                    solver.solveMechanicsLin();

                    VFePert[0] = LV.volume(disp, dETFESpace, - 1);
                    VFePert[1] = RV.volume(disp, dETFESpace, 1);

                    JFe(0,0) = ( VFePert[0] - VFeNew[0] ) / pPerturbationFe;
                    JFe(1,0) = ( VFePert[1] - VFeNew[1] ) / pPerturbationFe;

                    disp = dispCurrent;

                    // Right ventricle
                    modifyPressureBC(perturbedPressureComp(bcValues, pPerturbationFe, 1));
                    solver.bcInterfacePtr() -> updatePhysicalSolverVariables();
                    solver.solveMechanicsLin();

                    VFePert[0] = LV.volume(disp, dETFESpace, - 1);
                    VFePert[1] = RV.volume(disp, dETFESpace, 1);

                    JFe(0,1) = ( VFePert[0] - VFeNew[0] ) / pPerturbationFe;
                    JFe(1,1) = ( VFePert[1] - VFeNew[1] ) / pPerturbationFe;

                    disp = dispCurrent;
                }

                //============================================//
                // Update pressure b.c.
                //============================================//
                JR = JFe - JCirc;

                if ( JR.determinant() != 0 )
                {
                    dp = ( JR | R );
                    if ( iter > 5 ) dp *= 0.7;
                    if ( iter > 20 ) dp *= 0.5;
                    bcValues[0] -= std::min( std::max( dp(0) , - dpMax ) , dpMax );
                    bcValues[1] -= std::min( std::max( dp(1) , - dpMax ) , dpMax );
                }

                printCoupling("Pressure Update");

                //============================================//
                // Solve circulation
                //============================================//
                circulationSolver.iterate(dt_circulation, bcNames, bcValues, iter);
                VCircNew[0] = VCirc[0] + dt_circulation * ( Q("la", "lv") - Q("lv", "sa") );
                VCircNew[1] = VCirc[1] + dt_circulation * ( Q("ra", "rv") - Q("rv", "pa") );

                //============================================//
                // Solve mechanics
                //============================================//
                modifyPressureBC(bcValues);
                solver.bcInterfacePtr() -> updatePhysicalSolverVariables();
                solver.solveMechanics();

                VFeNew[0] = LV.volume(disp, dETFESpace, - 1);
                VFeNew[1] = RV.volume(disp, dETFESpace, 1);

                //============================================//
                // Residual update
                //============================================//
                R = VFeNew - VCircNew;
                printCoupling("Residual Update");
            }
 
            if ( 0 == comm->MyPID() )
            {
                std::cout << "\n*****************************************************************";
                std::cout << "\nCoupling converged after " << iter << " iteration" << ( iter > 1 ? "s" : "" );
                std::cout << "\n*****************************************************************\n\n";
            }
            
            //============================================//
            // Update volume variables
            //============================================//
            VCirc = VCircNew;
            VFe = VFeNew;
            
            //============================================//
            // Export circulation solution
            //============================================//
            if ( 0 == comm->MyPID() ) circulationSolver.exportSolution( circulationOutputFile );
            
        }
        
        //============================================//
        // Export FE-solution
        //============================================//
        bool save ( std::abs(std::remainder(t, dt_save)) < 0.01 );
        if ( save ) solver.saveSolution(t);

    }

    
    //============================================//
    // Close all exporters
    //============================================//
    solver.closeExporters();
    

#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return 0;
}
