#include <lifev/core/LifeV.hpp>
#include <lifev/electrophysiology/solver/IonicModels/IonicMinimalModel.hpp>
#include <lifev/electrophysiology/solver/IonicModels/IonicAlievPanfilov.hpp>

#include <lifev/em/solver/mechanics/EMStructuralOperator.hpp>
#include <lifev/em/solver/electrophysiology/EMMonodomainSolver.hpp>

#include <lifev/core/filter/ExporterEnsight.hpp>
#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5.hpp>
#endif
#include <lifev/core/filter/ExporterEmpty.hpp>

#include <lifev/eta/fem/ETFESpace.hpp>
#include <lifev/eta/expression/Integrate.hpp>


#include <lifev/bc_interface/3D/bc/BCInterface3D.hpp>

#include <lifev/em/solver/activation/activeStressModels/ActiveStressRossiModel14.hpp>

#include <sys/stat.h>

using namespace LifeV;


Real bcZero (const Real& /*t*/, const Real&  /*X*/, const Real& /*Y*/, const Real& /*Z*/, const ID& /*i*/)
{
    return  0.;
}

Real bcTraction (const Real& /*t*/, const Real&  /*X*/, const Real& /*Y*/, const Real& /*Z*/, const ID& /*i*/)
{
    return  50.;
}

Real d0 (const Real& /*t*/, const Real&  /*X*/, const Real& /*Y*/, const Real& /*Z*/, const ID& /*i*/)
{
    return  0.;
}

Real Iapp (const Real& t, const Real&  X, const Real& Y, const Real& Z, const ID& /*i*/)
{
    Real r = 0.1;
    Real t0 = 2;
    if (X < r && Y < r && Z < r && t < t0)
    {
        return 5.0;
    }
    else
    {
        return 0.0;
    }
}

Real initialV0left (const Real& /*t*/, const Real&  X, const Real& /*Y*/, const Real& /*Z*/, const ID& /*i*/)
{
    if ( X == 0 )
    {
        return 1.0;
    }
    else
    {
        return  0.;
    }
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
    //    typedef IonicMinimalModel                   ionicModel_Type;
    //    typedef boost::shared_ptr< ionicModel_Type >  ionicModelPtr_Type;
    //
    //    typedef EMMonodomainSolver< mesh_Type, ionicModel_Type >        monodomainSolver_Type;
    //    typedef boost::shared_ptr< monodomainSolver_Type >  monodomainSolverPtr_Type;
    typedef VectorEpetra                vector_Type;
    typedef boost::shared_ptr<vector_Type> vectorPtr_Type;

    typedef MatrixEpetra<Real> matrix_Type;
    typedef boost::shared_ptr<matrix_Type> matrixPtr_Type;

    typedef BCHandler                                          bc_Type;
    typedef boost::shared_ptr< bc_Type >                       bcPtr_Type;
    typedef StructuralOperator< RegionMesh<LinearTetra> >      physicalSolver_Type;
    typedef BCInterface3D< bc_Type, physicalSolver_Type >              bcInterface_Type;
    typedef boost::shared_ptr< bcInterface_Type >              bcInterfacePtr_Type;



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

    //********************************************//
    // Import parameters from an xml list. Use    //
    // Teuchos to create a list from a given file //
    // in the execution directory.                //
    //********************************************//

    if ( comm->MyPID() == 0 )
    {
        std::cout << "Importing parameters list...";
    }
    Teuchos::ParameterList parameterList = * ( Teuchos::getParametersFromXmlFile ( "ParamList.xml" ) );
    if ( comm->MyPID() == 0 )
    {
        std::cout << " Done!" << endl;
    }


    //********************************************//
    // In the parameter list we need to specify   //
    // the mesh name and the mesh path.           //
    //********************************************//
    if ( comm->MyPID() == 0 )
    {
        std::cout << "Reading Mesh Name and Path...\n";
    }

    std::string meshName = parameterList.get ("mesh_name", "lid16.mesh");
    std::string meshPath = parameterList.get ("mesh_path", "./");

    meshPtr_Type localSolidMesh ( new mesh_Type ( comm ) );
    meshPtr_Type fullSolidMesh ( new mesh_Type ( comm ) );
    MeshUtility::loadMesh (localSolidMesh, fullSolidMesh, meshName, meshPath);
    if ( comm->MyPID() == 0 )
    {
        std::cout << " Done!" << endl;
    }

    //********************************************//
    // We need the GetPot datafile for to setup   //
    // the preconditioner.                        //
    //********************************************//
    GetPot command_line (argc, argv);
    std::string problemFolder = command_line.follow ( "Output", 2, "-o", "--output" );
    // Create the problem folder
    if ( problemFolder.compare ("./") )
    {
        problemFolder += "/";

        if ( comm->MyPID() == 0 )
        {
            mkdir ( problemFolder.c_str(), 0777 );
        }
    }

    const string data_file_name = command_line.follow ("data", 2, "-f", "--file");
    GetPot dataFile (data_file_name);


    //===========================================================
    //===========================================================
    //              SOLID MECHANICS
    //===========================================================
    //===========================================================



    if ( comm->MyPID() == 0 )
    {
        std::cout << "monodomain: passed!" << std::endl;
    }

    typedef FESpace< RegionMesh<LinearTetra>, MapEpetra >               solidFESpace_Type;
    typedef boost::shared_ptr<solidFESpace_Type>                        solidFESpacePtr_Type;

    typedef ETFESpace< RegionMesh<LinearTetra>, MapEpetra, 3, 1 >       scalarETFESpace_Type;
    typedef boost::shared_ptr<scalarETFESpace_Type>                      scalarETFESpacePtr_Type;
    typedef ETFESpace< RegionMesh<LinearTetra>, MapEpetra, 3, 3 >       solidETFESpace_Type;
    typedef boost::shared_ptr<solidETFESpace_Type>                      solidETFESpacePtr_Type;
    if ( comm->MyPID() == 0 )
    {
        std::cout << "\n\ninitialization bc handler" << std::endl;
    }



    if ( comm->MyPID() == 0 )
    {
        std::cout << "\nparameters" << std::endl;
    }

    //    Real rho, poisson, young, bulk, alpha, gamma, mu;
    //    rho     = dataFile ( "solid/physics/density", 1. );
    //    young   = dataFile ( "solid/physics/young",   1. );
    //    poisson = dataFile ( "solid/physics/poisson", 1. );
    //    bulk    = dataFile ( "solid/physics/bulk",    1. );
    //    alpha   = dataFile ( "solid/physics/alpha",   1. );
    //    gamma   = dataFile ( "solid/physics/gamma",   1. );
    //    mu   = dataFile ( "solid/physics/mu",   1. );
    //    //  M_gammaf  = dataFile ( "solid/physics/gammaf",  0. );
    //
    //    if ( comm->MyPID() == 0 )
    //    {
    //        std::cout << "density = " << rho     << std::endl
    //                  << "young   = " << young   << std::endl
    //                  << "poisson = " << poisson << std::endl
    //                  << "bulk    = " << bulk    << std::endl
    //                  << "alpha   = " << alpha   << std::endl
    //                  << "gamma   = " << gamma   << std::endl;
    //    }


    if ( comm->MyPID() == 0 )
    {
        std::cout << "\ninitialization constitutive law" << std::endl;
    }

    boost::shared_ptr<StructuralConstitutiveLawData> dataStructure (new StructuralConstitutiveLawData( ) );
    dataStructure->setup (dataFile);

    if ( comm->MyPID() == 0 )
    {
        std::cout << "setup spaces" << std::endl;
    }

    std::string dOrder =  dataFile ( "solid/space_discretization/order", "P1");
    solidFESpacePtr_Type dFESpace ( new solidFESpace_Type (localSolidMesh, dOrder, 3, comm) );
    solidETFESpacePtr_Type dETFESpace ( new solidETFESpace_Type (localSolidMesh, & (dFESpace->refFE() ), & (dFESpace->fe().geoMap() ), comm) );


    if ( comm->MyPID() == 0 )
    {
        std::cout << "\nsetup boundary conditions" << std::endl;
    }

    //! #################################################################################
    //! BOUNDARY CONDITIONS
    //! #################################################################################
    bcInterfacePtr_Type                     solidBC ( new bcInterface_Type() );
    solidBC->createHandler();
    solidBC->fillHandler ( data_file_name, "solid" );
    solidBC->handler()->bcUpdate ( *dFESpace->mesh(), dFESpace->feBd(), dFESpace->dof() );


    if ( comm->MyPID() == 0 )
    {
        std::cout << "\nsetup structural operator" << std::endl;
    }
    //! 1. Constructor of the structuralSolver
    typedef EMStructuralOperator< RegionMesh<LinearTetra> > solid_Type;
    typedef boost::shared_ptr<solid_Type>                 solidPtr_Type;
    solid_Type solid;
    //    solid.setup(dataStructure,dFESpace, dETFESpace, solidBC -> handler(), comm);
    boost::shared_ptr<BCHandler> BCh (  solidBC -> handler() );
    solidBC -> handler() -> showMe();
    solid.setup (dataStructure, dFESpace, dETFESpace, solidBC -> handler(),  comm);
    BCh -> showMe();
    solid.EMMaterial() -> setupFiberVector (1, 0, 0);



    solid.setDataFromGetPot (dataFile);


    solid.buildSystem (1.0);
    vectorPtr_Type rhs (new vector_Type (solid.displacement(), Unique) );
    vectorPtr_Type disp (new vector_Type (solid.displacement(), Unique) );
    vectorPtr_Type initialDisplacement (new vector_Type (solid.displacement(), Unique) );
    solid.initialize ( initialDisplacement );


    boost::shared_ptr< Exporter<RegionMesh<LinearTetra> > > exporter;
    exporter.reset ( new ExporterHDF5<RegionMesh<LinearTetra> > ( dataFile, "structure" ) );

    exporter->setPostDir ( problemFolder );
    exporter->setMeshProcId ( localSolidMesh, comm->MyPID() );

    //    vectorPtr_Type solidDisp ( new vector_Type (solid.displacement(), exporter->mapType() ) );
    exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, "displacement", dFESpace, solid.displacementPtr(), UInt (0) );

    exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, "Fibers", dFESpace, solid.EMMaterial() -> fiberVectorPtr(), UInt (0) );









    typedef ElectroIonicModel ionicModel_Type;
    typedef boost::shared_ptr<ionicModel_Type> ionicModelPtr_Type;

    ionicModelPtr_Type ionicModel (new IonicAlievPanfilov() );

    //    typedef ElectroETAMonodomainSolver<mesh_Type, ionicModel_Type>  monodomain_Type;
    typedef EMMonodomainSolver<mesh_Type>  monodomain_Type;
    monodomain_Type monodomain ( meshName, meshPath, dataFile , ionicModel );
    std::cout << "Monodomain 1 Done!\n";

    monodomain.setDisplacementPtr ( solid.displacementPtr() );
    std::cout << "Displacement set!\n";

    monodomain.setMechanicsModifiesConductivity (false);
    std::cout << "Conductivity set!\n";

    monodomain.setInitialConditions();
    std::cout << "Initial Conditions set!\n";
    monodomain.setParameters ( parameterList );
    std::cout << "Parameters set!\n";

    monodomain.setFiberPtr ( solid.EMMaterial()->fiberVectorPtr() );
    std::cout << "Fiber vector set!\n";

    monodomain.setupMatrices();
    std::cout << "Matrices set!\n";

    std::cout << "\nSetting up the exporter ... " ;
    ExporterHDF5< RegionMesh <LinearTetra> > exporter2;
    monodomain.setupExporter ( exporter2, "ElectroSolution" , problemFolder);
    std::cout << " exporting initial solution ... \n" ;
    //   monodomain.exportSolution ( exporter2, 0);

    ExporterHDF5< RegionMesh <LinearTetra> > exporter3;
    monodomain.setupExporter ( exporter3, "ElectroSolution2" , problemFolder);
    std::cout << " exporting initial solution ... \n" ;
    //  monodomain.exportSolution ( exporter3, 0);

    function_Type stim = &Iapp;


    Real dt = 0.02;
    Real t = 0;

    EMData emdata;
    emdata.setup(dataFile);
    ActiveStressRossiModel14 activationModel;
    activationModel.setup(emdata, monodomain.potentialPtr()->map() );
    activationModel.setVariablesPtr(monodomain);
    std::cout << "Activation Model set!\n";

    vectorPtr_Type activationPtr ( activationModel.fiberActivationPtr() );


   solid.EMMaterial()->setFiberActivationPtr ( activationPtr);

    exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::ScalarField, "Activation", monodomain.feSpacePtr(), activationPtr, UInt (0) );
    exporter2.addVariable ( ExporterData<RegionMesh<LinearTetra> >::ScalarField, "Activation", monodomain.feSpacePtr(), activationPtr, UInt (0) );
    exporter3.addVariable ( ExporterData<RegionMesh<LinearTetra> >::ScalarField, "Activation", monodomain.feSpacePtr(), activationPtr, UInt (0) );

    for (int k (0); k <= 2000; k++)
    {
        std::cout << "\n*********************";
        std::cout << "\nTIME = " << t;
        std::cout << "\n*********************\n";
        if (k % 25 == 0)
        {
            solid.iterate ( solidBC -> handler() );
            exporter->postProcess ( t );
            monodomain.exportSolution ( exporter2, t);
        }

        monodomain.setAppliedCurrentFromFunction ( stim, t );
        monodomain.solveOneStepGatingVariablesFE();
        monodomain.solveOneICIStep();

        activationModel.solveModel ( dt );
        t += dt;

    }

    exporter2.closeFile();
    monodomain.setupMatrices();
    monodomain.setInitialConditions();
//    activationModel.activation() *= 0.0;
    monodomain.setMechanicsModifiesConductivity (true);

    t = 0;
    for (int k (0); k <= 2000; k++)
    {
        if (k % 50 == 0)
        {
            monodomain.exportSolution ( exporter3, t);
        }

        monodomain.setAppliedCurrentFromFunction ( stim, t );
        monodomain.solveOneStepGatingVariablesFE();
        monodomain.solveOneICIStep();

        activationModel.solveModel ( dt );
        t += dt;

    }


    exporter3.closeFile();

    exporter -> closeFile();
    ionicModel.reset();

#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return 0;
}
