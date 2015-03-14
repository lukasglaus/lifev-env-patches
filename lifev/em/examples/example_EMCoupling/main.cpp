//********************************************//
// Includes
//********************************************//

#include <lifev/core/LifeV.hpp>

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

// Namespaces
using namespace LifeV;


//********************************************//
// Functions
//********************************************//

Real Iapp (const Real& t, const Real&  X, const Real& Y, const Real& Z, const ID& /*i*/)
{
    return 0;
}

Real potentialMultiplyerFcn (const Real& t, const Real&  X, const Real& Y, const Real& Z, const ID& /*i*/)
{
    return ( Z > 0 && Y < 2.5 && Y > 0.5 ? 1.0 : 0.0 );
}


int main (int argc, char** argv)
{
    //********************************************//
    // Typedefs
    //********************************************//
    
    typedef RegionMesh<LinearTetra>                         mesh_Type;
    typedef boost::shared_ptr<mesh_Type>                    meshPtr_Type;
    
    typedef boost::function < Real (const Real & t,
                                    const Real &   x,
                                    const Real &   y,
                                    const Real & z,
                                    const ID&   /*i*/ ) >   function_Type;
    typedef VectorEpetra                                    vector_Type;
    typedef boost::shared_ptr<vector_Type>                  vectorPtr_Type;
    
    typedef EMMonodomainSolver<mesh_Type>                   monodomain_Type;


    //********************************************//
    // Declare comm and solver
    //********************************************//
    
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


    //********************************************//
    // Read data file and create output folder
    //********************************************//

    GetPot command_line (argc, argv);
    const std::string data_file_name = command_line.follow ("data", 2, "-f", "--file");
    GetPot dataFile (data_file_name);
    std::string problemFolder = EMUtility::createOutputFolder (command_line, *comm);

    
    //********************************************//
    // Load mesh
    //********************************************//
    
    if ( comm->MyPID() == 0 )
    {
        std::cout << "Load mesh...\n";
    }
    
    std::string meshName = dataFile("solid/space_discretization/mesh_file", "lid16.mesh");
    std::string meshPath = dataFile("solid/space_discretization/mesh_dir", "lid16.mesh");

    solver.loadMesh (meshName, meshPath);
    
    if ( comm->MyPID() == 0 )
    {
        std::cout << " Done!" << endl;
    }
    
    
    //********************************************//
    // Resize mesh
    //********************************************//
    
    if ( comm->MyPID() == 0 )
    {
        std::cout << "Resizing mesh..." << endl;
    }
    
    std::vector<Real> scale (3, 0.1);
    std::vector<Real> rotate (3, 0.0);
    std::vector<Real> translate (3, 0.0);
    MeshUtility::MeshTransformer<mesh_Type> transformerFull (* (solver.fullMeshPtr() ) );
    MeshUtility::MeshTransformer<mesh_Type> transformerLocal (* (solver.localMeshPtr() ) );
    transformerFull.transformMesh (scale, rotate, translate);
    transformerLocal.transformMesh (scale, rotate, translate);
    
    if ( comm->MyPID() == 0 )
    {
        std::cout << " Done!" << endl;
    }

    
    //********************************************//
    // Setup solver
    //********************************************//
    
    if( 0 == comm->MyPID() )
    {
    	std::cout << "Setting up EM solver ... ";
    }
    
    solver.setup (dataFile);
    
    if( 0 == comm->MyPID() )
    {
    	std::cout << " done!" << std::endl;
    }

    
    //********************************************//
    // Setup anisotropy vectors
    //********************************************//
    
    if( 0 == comm->MyPID() )
    {
    	std::cout << "Setting up anisotropy vectors ...";
    }

    solver.setupFiberVector (1., 0., 0.);
    string fiber_dir = dataFile( "electrophysiology/fiber/fiber", "./" );
    string fiber_name = dataFile( "electrophysiology/fiber/filename", "FiberDirection" );
    string fiber_fieldname = dataFile( "electrophysiology/fiber/fieldname", "fibers" );
    solver.setupElectroFiberVector ( fiber_name, fiber_fieldname, fiber_dir );
    // Or: VectorSmall<3> fibers; fibers[0] = 1; fibers[1] = 0; fibers[2] = 0; solver.setupElectroFiberVector(fibers);

    if( 0 == comm->MyPID() )
    {
    	std::cout << " done!" << std::endl;
    }

    
    //********************************************//
    // Initialize electrophysiology
    //********************************************//
    
    if( 0 == comm->MyPID() )
    {
    	std::cout << "Initialize electrophysiology ... ";
    }
    
    solver.initialize();
    
    // Set potential on certain flags
    UInt lvendo = dataFile( "electrophysiology/flags/lvendo", 36 );
    UInt rvendo = dataFile( "electrophysiology/flags/rvendo", 37 );
    UInt rvseptum = dataFile( "electrophysiology/flags/rvseptum", 38 );
    ElectrophysiologyUtility::setValueOnBoundary ( * (solver.electroSolverPtr()->potentialPtr() ), solver.fullMeshPtr(), 1.0, lvendo );
    ElectrophysiologyUtility::setValueOnBoundary ( * (solver.electroSolverPtr()->potentialPtr() ), solver.fullMeshPtr(), 1.0, rvendo );
    ElectrophysiologyUtility::setValueOnBoundary ( * (solver.electroSolverPtr()->potentialPtr() ), solver.fullMeshPtr(), 1.0, rvseptum);
    
    // Restrict the potential set by a function
    vectorPtr_Type potentialMultiplyer ( new vector_Type ( solver.electroSolverPtr()->potentialPtr()->map() ) ); // or: vectorPtr_Type potentialMultiplyer ( new vector_Type ( *solver.electroSolverPtr()->potentialPtr() ) );
    function_Type potMult = &potentialMultiplyerFcn;
    solver.electroSolverPtr()->feSpacePtr()->interpolate( potMult, *potentialMultiplyer, 0 );
    *solver.electroSolverPtr()->potentialPtr() *= *potentialMultiplyer;
    
    if( 0 == comm->MyPID() )
    {
    	std::cout << " done!" << std::endl;
    }

    
    //********************************************//
    // Building Matrices
    //********************************************//

    if( 0 == comm->MyPID() )
    {
    	std::cout << "Building matrices ... ";
    }

    solver.buildElectroSystem();
    
    if( 0 == comm->MyPID() )
    {
    	std::cout << " done!" << std::endl;
    }

    function_Type stim = &Iapp;

    
    //********************************************//
    // Setup exporter for EMSolver
    //********************************************//
    
    if ( 0 == comm->MyPID() )
    {
        std::cout << "Setting up exporters .. " << std::endl;
    }

    solver.setupExporters (problemFolder);
    
    if ( 0 == comm->MyPID() )
    {
        std::cout << " done!" << std::endl;
    }
    
    
    //********************************************//
    // Setup vector and exporter for activation t.
    //********************************************//
    
    vectorPtr_Type activationTimeVector ( new vector_Type ( solver.electroSolverPtr()->potentialPtr() -> map() ) );
    *activationTimeVector = -1.0;
    
    ExporterHDF5< RegionMesh <LinearTetra> > activationTimeExporter;
    activationTimeExporter.setMeshProcId (solver.localMeshPtr(), solver.comm()->MyPID() );
    activationTimeExporter.addVariable (ExporterData<mesh_Type>::ScalarField, "Activation Time",
                                        solver.electroSolverPtr()->feSpacePtr(), activationTimeVector, UInt (0) );
    activationTimeExporter.setPrefix ("ActivationTime");
    activationTimeExporter.setPostDir (problemFolder);
    
    
    //********************************************//
    // Time loop
    //********************************************//
    
    Real dt_activation = solver.data().electroParameter<Real>("timestep");
    Real endtime = solver.data().electroParameter<Real>("endtime");
    
    UInt maxiter = static_cast<UInt>( endtime / dt_activation ) ;
    Real t = 0;
    solver.saveSolution (0.0);

    for (int k (1); k <= maxiter; k++)
    {
        solver.electroSolverPtr() -> registerActivationTime (*activationTimeVector, t, 0.9);
        
        std::cout << "\n*********************";
        std::cout << "\nTIME = " << t+dt_activation;
        std::cout << "\n*********************\n";

        t = t + dt_activation;
        solver.solveElectrophysiology (stim, t);

        solver.saveSolution(t);
        activationTimeExporter.postProcess(t);// usually stored after time loop
    }

    
    //********************************************//
    // Close all exporters
    //********************************************//
    
    solver.closeExporters();
    activationTimeExporter.closeFile();



#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return 0;
}
