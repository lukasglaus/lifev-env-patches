#include <lifev/core/LifeV.hpp>

#include <lifev/em/solver/EMSolver.hpp>

#include <lifev/core/filter/ExporterEnsight.hpp>
#include <lifev/core/filter/ExporterVTK.hpp>
#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5.hpp>
#endif
#include <lifev/core/filter/ExporterEmpty.hpp>


// Resize mesh
#include <lifev/core/mesh/MeshUtility.hpp>


using namespace LifeV;

Real Iapp (const Real& t, const Real&  X, const Real& Y, const Real& Z, const ID& /*i*/)
{
    return 0;
}

Real potentialMultiplyerFcn (const Real& t, const Real&  X, const Real& Y, const Real& Z, const ID& /*i*/)
{
    return ( Y < 5 && Y > 4 ? 1.0 : 0.0 );
}



int main (int argc, char** argv)
{

    typedef RegionMesh<LinearTetra>                         mesh_Type;
    typedef boost::shared_ptr<mesh_Type>                    meshPtr_Type;
    typedef boost::function < Real (const Real & t,
                                    const Real &   x,
                                    const Real &   y,
                                    const Real & z,
                                    const ID&   /*i*/ ) >   function_Type;
    typedef VectorEpetra                                      vector_Type;
    typedef boost::shared_ptr<vector_Type>                    vectorPtr_Type;
    


#ifdef HAVE_MPI
    MPI_Init ( &argc, &argv );
#endif


    //===========================================================
    //===========================================================
    //              ELECTROPHYSIOLOGY
    //===========================================================
    //===========================================================


    boost::shared_ptr<Epetra_Comm>  comm ( new Epetra_MpiComm (MPI_COMM_WORLD) );
    if ( comm->MyPID() == 0 )
    {
        cout << "% using MPI" << endl;
    }

    typedef EMMonodomainSolver<mesh_Type>  monodomain_Type;
    EMSolver<mesh_Type, monodomain_Type> solver(comm);


    //********************************************//
    // In the parameter list we need to specify   //
    // the mesh name and the mesh path.           //
    //********************************************//
    if ( comm->MyPID() == 0 )
    {
        std::cout << "Reading Mesh Name and Path...\n";
    }

    GetPot command_line (argc, argv);
    const std::string data_file_name = command_line.follow ("data", 2, "-f", "--file");
    GetPot dataFile (data_file_name);
    
    std::string meshName = dataFile("solid/space_discretization/mesh_file", "lid16.mesh");
    std::string meshPath = dataFile("solid/space_discretization/mesh_dir", "lid16.mesh");

    solver.loadMesh (meshName, meshPath);
    
    // Resize mesh if necessary
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

 
    

    std::string problemFolder = EMUtility::createOutputFolder (command_line, *comm);




    if( 0 == comm->MyPID() )
    {
    	std::cout << "Setting up EM solver ... ";
    }
    
    solver.setup (dataFile);
    
    if( 0 == comm->MyPID() )
    {
    	std::cout << " done!" << std::endl;
    }

    if( 0 == comm->MyPID() )
    {
    	std::cout << "Setting up anisotropy vectors ...";
    }

//    VectorSmall<3> fibers;
//    fibers[0] = 1;
//    fibers[1] = 0;
//    fibers[2] = 0;
//    
//    solver.setupElectroFiberVector(fibers);
    solver.setupFiberVector (1., 0., 0.);
    string foldername = dataFile( "electrophysiology/fiber/foldername", "./" );
    string filename = dataFile( "electrophysiology/fiber/filename", "FiberDirection" );
    string fieldname = dataFile( "electrophysiology/fiber/fieldname", "fibers" );
    solver.setupElectroFiberVector ( filename, fieldname, foldername );
//    solver.setupSheetVector (0., 1., 0.);

    if( 0 == comm->MyPID() )
    {
    	std::cout << " done!" << std::endl;
    }

    if( 0 == comm->MyPID() )
    {
    	std::cout << "Setting solver properties ... ";
    }

    //solver.oneWayCoupling();
    
    if( 0 == comm->MyPID() )
    {
    	std::cout << " done!" << std::endl;
    }

    //Here we initialize the electrophysiology
    //with given intial conditions
    if( 0 == comm->MyPID() )
    {
    	std::cout << "Initialize electrophysiology ... ";
    }
    
    UInt endoFlag = dataFile( "electrophysiology/flags/lvendo", 36 );
    
    solver.initialize();
    ElectrophysiologyUtility::setValueOnBoundary ( * (solver.electroSolverPtr()->potentialPtr() ), solver.fullMeshPtr(), 1.0, endoFlag );
    
    // Create band
    vectorPtr_Type potentialMultiplyer ( new vector_Type ( solver.electroSolverPtr()->potentialPtr()->map() ) );
    // or: vectorPtr_Type potentialMultiplyer ( new vector_Type ( *solver.electroSolverPtr()->potentialPtr() ) );
    
    function_Type potMult = &potentialMultiplyerFcn;
    solver.electroSolverPtr()->feSpacePtr()->interpolate( potMult, *potentialMultiplyer, 0 );
    
    *solver.electroSolverPtr()->potentialPtr() *= *potentialMultiplyer;
    
    if( 0 == comm->MyPID() )
    {
    	std::cout << " done!" << std::endl;
    }

    //Here we are building the matrices
    //mass matrix for mechanic and the others for electrophysiology

    if( 0 == comm->MyPID() )
    {
    	std::cout << "Buildin matrices ... ";
    }

    solver.buildElectroSystem();
    
    if( 0 == comm->MyPID() )
    {
    	std::cout << " done!" << std::endl;
    }


    
    
    function_Type stim = &Iapp;

    Real dt_activation = solver.data().electroParameter<Real>("timestep");
    Real endtime = solver.data().electroParameter<Real>("endtime");

    UInt maxiter = static_cast<UInt>( endtime / dt_activation ) ;
    Real t = 0;

    if( 0 == comm->MyPID() )
    {
    	std::cout << "Setting up exporters .. ";
    }

    solver.setupExporters (problemFolder);

    if( 0 == comm->MyPID() )
    {
    	std::cout << " done!" << std::endl;
    }
    
    
// Create vector to check activation times
    // ---------------------------------------------------------------
    // We want to save the activation times in the domains.
    // Therefore, we create a vector which is initialized with the value -1.
    // At every timestep, we will check if the nodes in the mesh have
    // been activated, that is we check if the value of the potential
    // is bigger than a given threshold (which was defined at the beninning
    // when choosing the ionic model).
    // Moreover, we want to export the activation time. We therefore create
    // another HDF5 exporter to save the activation times on a separate file.
    // ---------------------------------------------------------------
    
    vectorPtr_Type activationTimeVector ( new vector_Type ( solver.electroSolverPtr()->potentialPtr() -> map() ) );
    *activationTimeVector = -1.0;
    
    ExporterHDF5< RegionMesh <LinearTetra> > activationTimeExporter;
    activationTimeExporter.setMeshProcId (solver.localMeshPtr(), solver.comm()->MyPID() );
    activationTimeExporter.addVariable (ExporterData<mesh_Type>::ScalarField, "Activation Time",
                                        solver.electroSolverPtr()->feSpacePtr(), activationTimeVector, UInt (0) );
    activationTimeExporter.setPrefix ("ActivationTime");
    activationTimeExporter.setPostDir (problemFolder);
    
    
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

    solver.closeExporters();
    activationTimeExporter.closeFile();



#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return 0;
}
