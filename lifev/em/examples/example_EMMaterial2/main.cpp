#include <lifev/core/LifeV.hpp>

#include <lifev/em/solver/EMSolver.hpp>

#include <lifev/core/filter/ExporterEnsight.hpp>
#include <lifev/core/filter/ExporterVTK.hpp>
#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5.hpp>
#endif
#include <lifev/core/filter/ExporterEmpty.hpp>


using namespace LifeV;

Real Iapp (const Real& t, const Real&  X, const Real& Y, const Real& Z, const ID& /*i*/)
{
    return (X < 0.25 && t < 0.5 ? 10 : 0);
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
    // Import parameters from an xml list. Use    //
    // Teuchos to create a list from a given file //
    // in the execution directory.                //
    //********************************************//



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
    
    if ( comm->MyPID() == 0 )
    {
        std::cout << " Done!" << endl;
    }

    //********************************************//
    // We need the GetPot datafile for to setup   //
    // the preconditioner.                        //
    //********************************************//
    
    

    std::string problemFolder = EMUtility::createOutputFolder (command_line, *comm);


    //===========================================================
    //===========================================================
    //              SOLID MECHANICS
    //===========================================================
    //===========================================================


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

    solver.initialize();
    
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

    solver.saveSolution (0.0);

    for (int k (1); k <= maxiter; k++)
    {
        std::cout << "\n*********************";
        std::cout << "\nTIME = " << t+dt_activation;
        std::cout << "\n*********************\n";

        t = t + dt_activation;
        solver.solveElectrophysiology (stim, t);

        solver.saveSolution(t);




    }

      solver.closeExporters();




#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return 0;
}
