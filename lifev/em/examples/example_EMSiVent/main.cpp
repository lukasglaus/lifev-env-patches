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

// Track nan
// #include <fenv.h>



// Namespaces
using namespace LifeV;


//============================================//
// Functions
//============================================//
               
Real Iapp (const Real& t, const Real&  X, const Real& Y, const Real& Z, const ID& /*i*/)
{
    bool coords ( Y > 1.5 && Y < 3 );
    bool time ( fmod(t, 800.) < 4 && fmod(t, 800.) > 2);
    return ( coords && time ? 30 : 0 );
}

Real potentialMultiplyerFcn (const Real& t, const Real&  X, const Real& Y, const Real& Z, const ID& /*i*/)
{
    bool time ( fmod(t, 800.) < 4 && fmod(t, 800.) > 2);
    return 1.4 * time; // ( Y < 2.5 && Y > 0.5 ? 1.0 : 0.0 );
}


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
    // Declare communicator and solver
    //============================================//
    
#ifdef HAVE_MPI
    MPI_Init ( &argc, &argv );
#endif

    boost::shared_ptr<Epetra_Comm>  comm ( new Epetra_MpiComm (MPI_COMM_WORLD) );
    if ( comm->MyPID() == 0 )
    {
        cout << "% using MPI" << std::endl;
    }
    //displayer->leaderPrint ("% using MPI");

    EMSolver<mesh_Type, monodomain_Type> solver(comm);

    
    //============================================//
    // Read data file and create output folder
    //============================================//

    GetPot command_line (argc, argv);
    const std::string data_file_name = command_line.follow ("data", 2, "-f", "--file");
    GetPot dataFile (data_file_name);
    std::string problemFolder = EMUtility::createOutputFolder (command_line, *comm);

    
    //============================================//
    // Setup material data
    //============================================//
    
    EMData emdata;
    emdata.setup (dataFile);
    
    
    //============================================//
    // Load mesh
    //============================================//
    
    if ( comm->MyPID() == 0 )
    {
        std::cout << "Load mesh...\n";
    }
    
    const std::string meshFile = command_line.follow ("biVent", 2, "-m", "--mesh");
    std::string elementOrder   =  dataFile ( "solid/space_discretization/order", "P1");
    std::string meshName;
    std::string meshPath;
    
    if ( meshFile == "fine" && elementOrder == "P1" )
    {
        meshName = "biVentFine.mesh";
        meshPath = "biVentFine/";
    }
    else if ( meshFile == "med" && elementOrder == "P1" )
    {
        meshName = "biVentMedium.mesh";
        meshPath = "biVentMedium/";
    }
    else if ( meshFile == "coa" && elementOrder == "P1" )
    {
        meshName = "biVentCoarse.mesh";
        meshPath = "biVentCoarse/";
    }
    else if ( meshFile == "med" && elementOrder == "P2" )
    {
        meshName = "biVentMedium.mesh";
        meshPath = "biVentMediumP2/";
    }
    else if ( meshFile == "coa" && elementOrder == "P2" )
    {
        meshName = "biVentCoarse.mesh";
        meshPath = "biVentCoarseP2/";
    }
    else
    {
        meshName = dataFile("solid/space_discretization/mesh_file", "cube4.mesh");
        meshPath = dataFile("solid/space_discretization/mesh_dir", "./");
    }
    
    solver.loadMesh (meshName, meshPath);
    
    if ( comm->MyPID() == 0 )
    {
        std::cout << " Done!" << endl;
    }
    
    
    //============================================//
    // Resize mesh
    //============================================//
    
    if ( comm->MyPID() == 0 )
    {
        std::cout << "Resizing mesh..." << endl;
    }
    
    std::vector<Real> scale (3, dataFile("solid/space_discretization/mesh_scaling", 1.0));
    std::vector<Real> rotate { dataFile("solid/space_discretization/mesh_rotation_0", 0.0) , dataFile("solid/space_discretization/mesh_rotation_1", 0.0) , dataFile("solid/space_discretization/mesh_rotation_2", 0.0) };
    std::vector<Real> translate { dataFile("solid/space_discretization/mesh_translation_0", 0.0) , dataFile("solid/space_discretization/mesh_translation_1", 0.0) , dataFile("solid/space_discretization/mesh_translation_2", 0.0) };
    
    MeshUtility::MeshTransformer<mesh_Type> transformerFull (* (solver.fullMeshPtr() ) );
    MeshUtility::MeshTransformer<mesh_Type> transformerLocal (* (solver.localMeshPtr() ) );
    
    transformerFull.transformMesh (scale, rotate, translate);
    transformerLocal.transformMesh (scale, rotate, translate);
    
    if ( comm->MyPID() == 0 )
    {
        std::cout << " Done!" << endl;
    }
    
    
    //============================================//
    // Setup solver (including fe-spaces & b.c.)
    //============================================//
    
    if( 0 == comm->MyPID() )
    {
    	std::cout << "Setting up EM solver ... ";
    }
    
    solver.setup (dataFile);
    
    if( 0 == comm->MyPID() )
    {
    	std::cout << " done!" << std::endl;
    }

    
    //============================================//
    // Setup anisotropy vectors
    //============================================//
    
    if( 0 == comm->MyPID() )
    {
    	std::cout << "Setting up anisotropy vectors ...";
    }
 
    bool anisotropy = dataFile ( "solid/space_discretization/anisotropic", false );

    if ( anisotropy )
    {
        std::string fiberFileName  =  dataFile ( "solid/space_discretization/fiber_name", "FiberDirection");
        std::string sheetFileName  =  dataFile ( "solid/space_discretization/sheet_name", "SheetsDirection");
        std::string fiberFieldName =  dataFile ( "solid/space_discretization/fiber_fieldname", "fibers");
        std::string sheetFieldName =  dataFile ( "solid/space_discretization/sheet_fieldname", "sheets");
        std::string fiberDir       =  meshPath; //dataFile ( "solid/space_discretization/fiber_dir", "./");
        std::string sheetDir       =  meshPath; //dataFile ( "solid/space_discretization/sheet_dir", "./");
        //std::string elementOrder   =  dataFile ( "solid/space_discretization/order", "P1");

        solver.setupFiberVector ( fiberFileName, fiberFieldName, fiberDir, elementOrder );
        solver.setupMechanicalSheetVector ( sheetFileName, sheetFieldName, sheetDir, elementOrder );
    }
    else
    {
        solver.setupFiberVector (1., 0., 0.);
        solver.setupSheetVector (0., 1., 0.);
    }
    
    if( 0 == comm->MyPID() )
    {
    	std::cout << " done!" << std::endl;
    }

    
    //============================================//
    // Initialize electrophysiology
    //============================================//
    
    if( 0 == comm->MyPID() )
    {
    	std::cout << "Initialize electrophysiology ... ";
    }
    
    solver.initialize();
    
    // Set potential on certain flags
    UInt lvendo = dataFile( "electrophysiology/flags/lvendo", 36 );
    //UInt rvendo = dataFile( "electrophysiology/flags/rvendo", 37 );
    //UInt rvseptum = dataFile( "electrophysiology/flags/rvseptum", 38 );
    ElectrophysiologyUtility::setValueOnBoundary ( * (solver.electroSolverPtr()->potentialPtr() ), solver.fullMeshPtr(), 1.0, lvendo );
    //ElectrophysiologyUtility::setValueOnBoundary ( * (solver.electroSolverPtr()->potentialPtr() ), solver.fullMeshPtr(), 1.0, rvendo );
    //ElectrophysiologyUtility::setValueOnBoundary ( * (solver.electroSolverPtr()->potentialPtr() ), solver.fullMeshPtr(), 1.0, rvseptum);
    
    // Restrict the potential set by a function
    vectorPtr_Type potentialMultiplyer ( new vector_Type ( solver.electroSolverPtr()->potentialPtr()->map() ) ); // or: vectorPtr_Type potentialMultiplyer ( new vector_Type ( *solver.electroSolverPtr()->potentialPtr() ) );
    function_Type potMult = &potentialMultiplyerFcn;
    solver.electroSolverPtr()->feSpacePtr()->interpolate( potMult, *potentialMultiplyer, 0 );
    *solver.electroSolverPtr()->potentialPtr() *= *potentialMultiplyer;
    
    if( 0 == comm->MyPID() )
    {
    	std::cout << " done!" << std::endl;
    }

    
    //============================================//
    // Building Matrices
    //============================================//

    if( 0 == comm->MyPID() )
    {
    	std::cout << "Building matrices ... ";
    }
    
    solver.oneWayCoupling();
    solver.structuralOperatorPtr()->setNewtonParameters(dataFile);
    solver.buildSystem();
    
    if( 0 == comm->MyPID() )
    {
    	std::cout << " done!" << std::endl;
    }

    
    //============================================//
    // Setup exporter for EMSolver
    //============================================//
    
    if ( 0 == comm->MyPID() )
    {
        std::cout << "Setting up exporters .. " << std::endl;
    }

    solver.setupExporters (problemFolder);
    
    if ( 0 == comm->MyPID() )
    {
        std::cout << " done!" << std::endl;
    }
    
    
    //============================================//
    // Setup vector & exporter for activation time
    //============================================//
    
    vectorPtr_Type activationTimeVector ( new vector_Type ( solver.electroSolverPtr()->potentialPtr() -> map() ) );
    *activationTimeVector = -1.0;
    
    ExporterHDF5< RegionMesh <LinearTetra> > activationTimeExporter;
    activationTimeExporter.setMeshProcId (solver.localMeshPtr(), solver.comm()->MyPID() );
    activationTimeExporter.addVariable (ExporterData<mesh_Type>::ScalarField, "Activation Time", solver.electroSolverPtr()->feSpacePtr(), activationTimeVector, UInt (0) );
    activationTimeExporter.setPrefix ("ActivationTime");
    activationTimeExporter.setPostDir (problemFolder);
    
    
    //============================================//
    // Electric stimulus function
    //============================================//
    
    function_Type stim = &Iapp;
    
    
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
    // Kept-normal boundary conditions
    //============================================//

    // Get b.c. flags
    // ID LvFlag =  dataFile ( "solid/boundary_conditions/LvFlag", 0);
    
    // Boundary vector normal in deformed configuration
    // solver.structuralOperatorPtr() -> setBCFlag( LvFlag );
    
    //solver.bcInterfacePtr() -> handler() -> addBC("LvPressure", LVFlag, Natural, Full, *pLvBCVectorPtr, 3); // BC for using function which keeps bc normal
    // Todo: Normal boundary condition!!
    
    
    //============================================//
    // Modifiable-value boundary condition
    //============================================//

    std::vector<vectorPtr_Type> pVecPtrs;
    std::vector<bcVectorPtr_Type> pBCVecPtrs;
    
    UInt nVarBC = dataFile.vector_variable_size ( ( "solid/boundary_conditions/listVariableBC" ) );
    for ( UInt i (0) ; i < nVarBC ; ++i )
    {
        std::string varBCSection = dataFile ( ( "solid/boundary_conditions/listVariableBC" ), " ", i );
        ID flag  =  dataFile ( ("solid/boundary_conditions/" + varBCSection + "/flag").c_str(), 0 );
        ID index =  dataFile ( ("solid/boundary_conditions/" + varBCSection + "/index").c_str(), 0 );

        pVecPtrs.push_back ( vectorPtr_Type ( new vector_Type ( solver.structuralOperatorPtr() -> displacement().map(), Repeated ) ) );
        pBCVecPtrs.push_back ( bcVectorPtr_Type( new bcVector_Type( *pVecPtrs[i], solver.structuralOperatorPtr() -> dispFESpacePtr() -> dof().numTotalDof(), 1 ) ) );
        //solver.bcInterfacePtr() -> handler() -> addBC(varBCSection, flag, Natural, Normal, *pBCVecPtrs[i]);
        solver.bcInterfacePtr() -> handler() -> addBC(varBCSection, flag, Natural, Full, *pBCVecPtrs[i], 3);
    }

    solver.bcInterfacePtr() -> handler() -> bcUpdate( *solver.structuralOperatorPtr() -> dispFESpacePtr() -> mesh(), solver.structuralOperatorPtr() -> dispFESpacePtr() -> feBd(), solver.structuralOperatorPtr() -> dispFESpacePtr() -> dof() );
    
    if ( 0 == comm->MyPID() ) solver.bcInterfacePtr() -> handler() -> showMe();

    
    auto modifyFeBC = [&] (const std::vector<Real>& bcValues)
    {
        for ( UInt i (0) ; i < nVarBC ; ++i )
        {
            std::string varBCSection = dataFile ( ( "solid/boundary_conditions/listVariableBC" ), " ", i );
            ID flag  =  dataFile ( ("solid/boundary_conditions/" + varBCSection + "/flag").c_str(), 0 );
            ID index =  dataFile ( ("solid/boundary_conditions/" + varBCSection + "/index").c_str(), 0 );
            *pVecPtrs[i] = - bcValues[index] * 1333.224; // 0.001333224
            pBCVecPtrs[i].reset ( ( new bcVector_Type (*pVecPtrs[i], solver.structuralOperatorPtr() -> dispFESpacePtr() -> dof().numTotalDof(), 1) ) );
            solver.bcInterfacePtr() -> handler() -> modifyBC(flag, *pBCVecPtrs[i]);
        }
    };


    //============================================//
    // Volume integrators
    //============================================//
    
    auto& disp = solver.structuralOperatorPtr() -> displacement();
    auto FESpace = solver.structuralOperatorPtr() -> dispFESpacePtr();
    auto dETFESpace = solver.electroSolverPtr() -> displacementETFESpacePtr();
    auto ETFESpace = solver.electroSolverPtr() -> ETFESpacePtr();
    
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
    //VolumeIntegrator RV (RVFlags, "Right Ventricle", solver.fullMeshPtr(), solver.localMeshPtr(), ETFESpace, FESpace);

    
    //============================================//
    // Set variables and functions
    //============================================//
    
    Real dt_activation = solver.data().electroParameter<Real>("timestep");
    Real dt_mechanics = solver.data().solidParameter<Real>("timestep");
    Real endtime = solver.data().electroParameter<Real>("endtime");
    UInt saveIter = static_cast<UInt>( dt_mechanics / dt_activation );
    UInt maxiter = static_cast<UInt>( endtime / dt_activation ) ;
    
    Real pPerturbationFe = dataFile ( "solid/coupling/pPerturbationFe", 1e-2 );
    Real pPerturbationCirc = dataFile ( "solid/coupling/pPerturbationCirc", 1e-3 );
    Real couplingError = dataFile ( "solid/coupling/couplingError", 1e-6 );
    UInt couplingJFeSubIter = dataFile ( "solid/coupling/couplingJFeSubIter", 1 );
    UInt couplingJFeSubStart = dataFile ( "solid/coupling/couplingJFeSubStart", 1 );
    UInt couplingJFeIter = dataFile ( "solid/coupling/couplingJFeIter", 1 );

    Real dpMax = dataFile ( "solid/coupling/dpMax", 0.1 );
    
    std::vector<std::vector<std::string> > bcNames { { "lv" , "p" } /*, { "rv" , "p" } */};
    std::vector<double> bcValues { p ( "lv" ) /*, p ( "rv") */};
    
    VectorSmall<1> VCirc, VCircNew, VCircPert, VFe, VFeNew, VFePert, R, dp;
    MatrixSmall<1,1> JFe, JCirc, JR;

    UInt iter (0);
    Real t (0);
    
    auto printCoupling = [&] ( std::string label ) { if ( 0 == comm->MyPID() )
    {
        std::cout << "\n****************************** Coupling: " << label << " *******************************";
        std::cout << "\nNewton iteration nr. " << iter << " at time " << t;
        std::cout << "\nLV - Pressure: " << bcValues[0];
        std::cout << "\nLV - FE-Volume (Current - Pert - New - J): \t\t" << VFe[0] << "\t" << VFePert[0] << "\t" << VFeNew[0];
        std::cout << "\nLV - Circulation-Volume (Current - Pert - New - J): \t" << VCirc[0] << "\t" << VCircPert[0] << "\t" << VCircNew[0];
        std::cout << "\nLV - Residual = " << std::abs(VFeNew[0] - VCircNew[0]);
//        std::cout << "\nRV - Pressure: " << bcValues[1];
//        std::cout << "\nRV - FE-Volume (Current - Pert - New - J): \t\t" << VFe[1] << "\t" << VFePert[1] << "\t" << VFeNew[1];
//        std::cout << "\nRV - Circulation-Volume (Current - Pert - New - J): \t" << VCirc[1] << "\t" << VCircPert[1] << "\t" << VCircNew[1];
//        std::cout << "\nRV - Residual = " << std::abs(VFeNew[1] - VCircNew[1]);
        std::cout << "\nJFe   = " << JFe;
        std::cout << "\nJCirc = " << JCirc;
        std::cout << "\nJR    = " << JR;
        std::cout << "\n****************************** Coupling: " << label << " *******************************\n\n"; }
    };
    
    auto pipeToString = [] ( const char* command )
    {
        FILE* file = popen( command, "r" ) ;
        std::ostringstream stm ;
        char line[6] ;
        fgets( line, 6, file );
        stm << line;
        pclose(file) ;
        return stm.str() ;
    };
    
    
    //============================================//
    // Load restart file
    //============================================//
    
    std::string restartInput = command_line.follow ("noRestart", 2, "-r", "--restart");
    const bool restart ( restartInput != "noRestart" );

//    if ( restart )
//    {
//        const std::string restartDir = command_line.follow (problemFolder.c_str(), 2, "-rd", "--restartDir");
//        
//        // Get most recent restart index
//        if ( restartInput == "." )
//        {
//            restartInput = pipeToString( ("tail -n 1 " + restartDir + "solution.dat | awk -F '[. ]' '{print $1 \".\" $2}' | awk '{printf \"%05g\", $1*1000/" + std::to_string(dt_activation) + " + 1}'").c_str() );
//        }
//        
//        std::cout << comm->MyPID() << " ----------------------- " << restartInput << std::endl;
//
//        
//        // Set time variable
//        const unsigned int restartInputStr = std::stoi(restartInput);
//        const unsigned int nIter = (restartInputStr - 1) / saveIter;
//        t = nIter * dt_mechanics;
//        
//        activationTimeExporter.setTimeIndex(nIter * saveIter + 2);
//        solver.setTimeIndex(nIter * saveIter + 2);
//
//        std::string polynomialDegree = dataFile ( "solid/space_discretization/order", "P1");
//        
//        // Todo: add other electrophysiology variables! solver.electroSolverPtr() -> importSolution ("ElectroSolution", problemFolder, t);
//        
//        ElectrophysiologyUtility::importVectorField ( solver.structuralOperatorPtr() -> displacementPtr(), "MechanicalSolution" , "displacement", solver.localMeshPtr(), restartDir, polynomialDegree, restartInput );
//
//        ElectrophysiologyUtility::importScalarField (solver.electroSolverPtr()->globalSolution().at(0), "ElectroSolution" , "Variable0", solver.localMeshPtr(), restartDir, polynomialDegree, restartInput );
//        ElectrophysiologyUtility::importScalarField (solver.electroSolverPtr()->globalSolution().at(1), "ElectroSolution" , "Variable1", solver.localMeshPtr(), restartDir, polynomialDegree, restartInput );
//        ElectrophysiologyUtility::importScalarField (solver.electroSolverPtr()->globalSolution().at(2), "ElectroSolution" , "Variable2", solver.localMeshPtr(), restartDir, polynomialDegree, restartInput );
//        ElectrophysiologyUtility::importScalarField (solver.electroSolverPtr()->globalSolution().at(3), "ElectroSolution" , "Variable3", solver.localMeshPtr(), restartDir, polynomialDegree, restartInput );
//        
//        ElectrophysiologyUtility::importScalarField (solver.activationModelPtr() -> fiberActivationPtr(), "ActivationSolution" , "Activation", solver.localMeshPtr(), restartDir, polynomialDegree, restartInput );
//
//        
//        ElectrophysiologyUtility::importScalarField (activationTimeVector, "ActivationTime" , "Activation Time", solver.localMeshPtr(), restartDir, polynomialDegree, restartInput );
//        
//        
//        std::cout << restartDir + "solution.dat" << "    " << nIter << std::endl;
//        
//        // Circulation
//        circulationSolver.restartFromFile ( restartDir + "solution.dat" , nIter );
//        //circulationSolver.exportSolution( circulationOutputFile );
//
//        // Coupling boundary conditions
//        bcValues = { p ( "lv" ) , p ( "rv" ) };
//        
//        std::cout << "Volume 0: " << LV.volume(disp, dETFESpace, - 1) << " " << RV.volume(disp, dETFESpace, 1) << std::endl;
//        
//        modifyFeBC(bcValues);
//        solver.bcInterfacePtr() -> updatePhysicalSolverVariables();
//        solver.solveMechanics();
//        
//        std::cout << "Volume 1: " << LV.volume(disp, dETFESpace, - 1) << " " << RV.volume(disp, dETFESpace, 1) << std::endl;
//        
//        // Adjust time step
//        Real timestepFactor = dataFile ( "solid/time_discretization/timestepRestartFactor", 1. );
//        dt_mechanics *= timestepFactor;
//        saveIter = static_cast<UInt>( dt_mechanics / dt_activation );
//
//    }

    
    //============================================//
    // Preload
    //============================================//
    
    if ( ! restart )
    {
        const int preloadSteps = dataFile ( "solid/boundary_conditions/numPreloadSteps", 0);
        
        auto preloadPressure = [] (std::vector<double> p, const int& step, const int& steps)
        {
            for (auto& i : p) {i *= double(step) / double(steps);}
            return p;
        };
        
        solver.saveSolution (-1.0);
        activationTimeExporter.postProcess(-1.0);
        
        for (int i (1); i <= preloadSteps; i++)
        {
            if ( 0 == comm->MyPID() )
            {
                std::cout << "\n*********************";
                std::cout << "\nPreload step: " << i << " / " << preloadSteps;
                std::cout << "\n*********************\n";
            }
            
            // Update pressure b.c.
            modifyFeBC(preloadPressure(bcValues, i, preloadSteps));

            // Solve mechanics
            solver.bcInterfacePtr() -> updatePhysicalSolverVariables();
            solver.solveMechanics();
            //solver.saveSolution ( i );
            //activationTimeExporter.postProcess( i );

        }
    }
    
    
    //============================================//
    // Time loop
    //============================================//
    
    VFe[0] = LV.volume(disp, dETFESpace, - 1);
    //VFe[1] = RV.volume(disp, dETFESpace, 1);
    VCirc = VFe;
    
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
        activationTimeExporter.postProcess(t);
        circulationSolver.exportSolution( circulationOutputFile );
    }
    
    for (int k (1); k <= maxiter; k++)
    {
        if ( 0 == comm->MyPID() )
        {
            std::cout << "\n*********************";
            std::cout << "\nTIME = " << t+dt_activation;
            std::cout << "\n*********************\n";
        }

        t = t + dt_activation;

        //============================================//
        // Solve electrophysiology and activation
        //============================================//
        
        solver.electroSolverPtr() -> registerActivationTime (*activationTimeVector, t, 0.9);
        solver.solveElectrophysiology (stim, t);
        solver.solveActivation (dt_activation);

        if ( k % saveIter == 0 )
        {
            iter = 0;
            const double dt_circulation ( dt_mechanics / 1000 );
            solver.structuralOperatorPtr() -> data() -> dataTime() -> setTime(t);
            
            //============================================//
            // Solve mechanics
            //============================================//
            modifyFeBC(bcValues);
            solver.bcInterfacePtr() -> updatePhysicalSolverVariables();
            solver.solveMechanics();
            
            VFeNew[0] = LV.volume(disp, dETFESpace, - 1);

            //============================================//
            // Solve circlation
            //============================================//
            circulationSolver.iterate(dt_circulation, bcNames, bcValues, iter);
            VCircNew[0] = VCirc[0] + dt_circulation * ( Q("la", "lv") - Q("lv", "sa") );

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
                JCirc(0,0) = ( VCircPert[0] - VCircNew[0] ) / pPerturbationCirc;

                //============================================//
                // Jacobian fe
                //============================================//

                const bool jFeIter ( ! ( k % (couplingJFeIter * saveIter) ) );
                const bool jFeSubIter ( ! ( (iter - couplingJFeSubStart) % couplingJFeSubIter) && iter >= couplingJFeSubStart );
                const bool jFeEmpty ( JFe.norm() == 0 );
                
                if ( jFeIter || jFeSubIter || jFeEmpty )
                {
                    JFe *= 0.0;

                    // Left ventricle
                    modifyFeBC(perturbedPressureComp(bcValues, pPerturbationFe, 0));
                    solver.bcInterfacePtr() -> updatePhysicalSolverVariables();
                    solver.solveMechanics();
                    
                    VFePert[0] = LV.volume(disp, dETFESpace, - 1);
                    JFe(0,0) = ( VFePert[0] - VFeNew[0] ) / pPerturbationFe;
                }

                //============================================//
                // Update pressure b.c.
                //============================================//
                JR = JFe - JCirc;

                if ( JR.determinant() != 0 )
                {
                    dp = ( JR | R );
                    if ( iter > 5 ) dp *= 0.7;
                    if ( iter > 20 ) dp *= 0.7;
                    bcValues[0] -= std::min( std::max( dp(0) , - dpMax ) , dpMax );
                    bcValues[1] = bcValues[0] / 6;
                }
                
                printCoupling("Pressure Update");

                //============================================//
                // Solve circulation
                //============================================//
                circulationSolver.iterate(dt_circulation, bcNames, bcValues, iter);
                VCircNew[0] = VCirc[0] + dt_circulation * ( Q("la", "lv") - Q("lv", "sa") );

                //============================================//
                // Solve mechanics
                //============================================//
                modifyFeBC(bcValues);
                solver.bcInterfacePtr() -> updatePhysicalSolverVariables();
                solver.solveMechanics();
                
                VFeNew[0] = LV.volume(disp, dETFESpace, - 1);

                //============================================//
                // Residual update
                //============================================//
                R = VFeNew - VCircNew;
                printCoupling("Residual Update");
            }
 
            if ( 0 == comm->MyPID() )
            {
                std::cout << "\n******************************************";
                std::cout << "\nCoupling converged after " << iter << " iteration" << ( iter > 1 ? "s" : "" );
                std::cout << "\n******************************************\n\n";
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
        solver.saveSolution(t);
        activationTimeExporter.postProcess(t);
    }

    
    //============================================//
    // Close all exporters
    //============================================//
    
    solver.closeExporters();
    activationTimeExporter.closeFile();
    

#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return 0;
}
