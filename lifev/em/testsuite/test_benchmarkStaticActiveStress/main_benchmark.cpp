#include <lifev/core/LifeV.hpp>
#include <lifev/electrophysiology/solver/ElectroETAMonodomainSolver.hpp>
#include <lifev/electrophysiology/solver/IonicModels/IonicMinimalModel.hpp>

#include <lifev/structure/solver/StructuralConstitutiveLawData.hpp>

#include <lifev/structure/solver/StructuralConstitutiveLaw.hpp>
#include <lifev/structure/solver/StructuralOperator.hpp>

#include <lifev/em/solver/mechanics/EMStructuralOperator.hpp>
#include <lifev/em/solver/mechanics/EMStructuralConstitutiveLaw.hpp>
#include <lifev/em/solver/EMETAFunctors.hpp>
#include <lifev/em/solver/EMSolver.hpp>


#include <lifev/core/filter/ExporterEnsight.hpp>
#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5.hpp>
#endif
#include <lifev/core/filter/ExporterEmpty.hpp>

#include <lifev/eta/fem/ETFESpace.hpp>
#include <lifev/eta/expression/Integrate.hpp>

//#include <lifev/em/solver/mechanics/materials/EMMaterial.hpp>


#include <lifev/bc_interface/3D/bc/BCInterface3D.hpp>
#include <lifev/structure/solver/NeoHookeanMaterialNonLinear.hpp>

using namespace LifeV;





template<typename space> void ComputeBC( const boost::shared_ptr<RegionMesh<LinearTetra> > localMesh,
						const VectorEpetra& disp,
						boost::shared_ptr<VectorEpetra> bcVectorPtr,
						const boost::shared_ptr< space > dETFESpace,
						int bdFlag);




int main (int argc, char** argv)
{

#ifdef HAVE_MPI
    MPI_Init ( &argc, &argv );
#endif

    //start with the communicator
    boost::shared_ptr<Epetra_Comm>  comm ( new Epetra_MpiComm (MPI_COMM_WORLD) );
    if ( comm->MyPID() == 0 )
    {
        cout << "% using MPI" << std::endl;
    }

    //===========================================================
    //===========================================================
    //              READ DATAFILE and CREATE OUTPUT FOLDER
    //===========================================================
    //===========================================================
    GetPot command_line (argc, argv);
    const string data_file_name = command_line.follow ("data", 2, "-f", "--file");
    GetPot dataFile (data_file_name);

    //When launching the executable use the flag:  -o OutputFolderName
    std::string problemFolder = EMUtility::createOutputFolder (command_line, *comm);


    //===========================================================
    //===========================================================
    //              LOAD MESH
    //===========================================================
    //===========================================================
    if ( comm->MyPID() == 0 )
    {
        std::cout << "Reading Mesh Name and Path...\n";
    }

    std::string meshName = dataFile ( "solid/space_discretization/mesh_file", "" );
    std::string meshPath = dataFile ( "solid/space_discretization/mesh_dir", "./" );

    typedef RegionMesh<LinearTetra>                         mesh_Type;
    typedef boost::shared_ptr<mesh_Type>                    meshPtr_Type;

    meshPtr_Type localSolidMesh ( new mesh_Type ( comm ) );
    meshPtr_Type fullSolidMesh ( new mesh_Type ( comm ) );

    MeshUtility::loadMesh (localSolidMesh, fullSolidMesh, meshName, meshPath);

    if ( comm->MyPID() == 0 )
    {
        std::cout << " Done!" << std::endl;
    }

    //===========================================================
    //===========================================================
    //              FINITE ELEMENT SPACES
    //===========================================================
    //===========================================================

    //Define the finite element space for bc and exporter
    // and the ET finite element space for assembly
    if ( comm->MyPID() == 0 )
    {
        std::cout << "setup spaces ... ";
    }

    typedef FESpace< RegionMesh<LinearTetra>, MapEpetra >               solidFESpace_Type;
    typedef boost::shared_ptr<solidFESpace_Type>                        solidFESpacePtr_Type;

    typedef ETFESpace< RegionMesh<LinearTetra>, MapEpetra, 3, 3 >       solidETFESpace_Type;
    typedef boost::shared_ptr<solidETFESpace_Type>                      solidETFESpacePtr_Type;



    std::string dOrder =  dataFile ( "solid/space_discretization/order", "P1");
    solidFESpacePtr_Type dFESpace ( new solidFESpace_Type (localSolidMesh, dOrder, 3, comm) );
    solidETFESpacePtr_Type dETFESpace ( new solidETFESpace_Type (localSolidMesh, & (dFESpace->refFE() ), & (dFESpace->fe().geoMap() ), comm) );

    if ( comm->MyPID() == 0 )
    {
        std::cout << " Done!" << std::endl;
    }

    //===========================================================
    //===========================================================
    //              BOUNDARY CONDITIONS Part I
    //===========================================================
    //===========================================================
    if ( comm->MyPID() == 0 )
    {
        std::cout << "setup bc ... ";
    }

    typedef BCHandler                                          bc_Type;
    typedef boost::shared_ptr< bc_Type >                       bcPtr_Type;
    typedef StructuralOperator< RegionMesh<LinearTetra> >      physicalSolver_Type;
    typedef BCInterface3D< bc_Type, physicalSolver_Type >      bcInterface_Type;
    typedef boost::shared_ptr< bcInterface_Type >              bcInterfacePtr_Type;

    bcInterfacePtr_Type                     solidBC ( new bcInterface_Type() );
    solidBC->createHandler();
    solidBC->fillHandler ( data_file_name, "solid" );

    if ( comm->MyPID() == 0 )
    {
        std::cout << " Done!" << std::endl;
    }

    //===========================================================
    //===========================================================
    //              SOLID MECHANICS
    //===========================================================
    //===========================================================

    if ( comm->MyPID() == 0 )
    {
        std::cout << "setup constitutive law data ... ";
    }

    //Setup the data of the constitutive law
    boost::shared_ptr<StructuralConstitutiveLawData> dataStructure (new StructuralConstitutiveLawData( ) );
    dataStructure->setup (dataFile);

    if ( comm->MyPID() == 0 )
    {
        std::cout << " Done!" << std::endl;
    }

    //Setup the structural operator
    if ( comm->MyPID() == 0 )
    {
        std::cout << "setup structural operator ... " << std::endl;
    }
    //! 1. Constructor of the structuralSolver
    EMStructuralOperator< RegionMesh<LinearTetra> > solid;

    solid.setup ( dataStructure, dFESpace, dETFESpace, solidBC -> handler(), comm);
    solid.setDataFromGetPot (dataFile);
    solid.EMMaterial() -> setupFiberVector( 1.0, 0.0, 0.0);
    solid.EMMaterial() -> setupSheetVector( 0.0, 1.0, 0.0);
    if(solid.EMMaterial()->sheetVectorPtr()) std::cout << "I have sheets!\n";
    else
    {
        std::cout << "I don't have sheets!\n";
        return 0.0;
    }

    solid.buildSystem (1.0);

    if ( comm->MyPID() == 0 )
    {
        std::cout << " Done!" << std::endl;
    }


    //===========================================================
    //===========================================================
    //              IMPORTER
    //===========================================================
    //===========================================================
    vectorPtr_Type importedSolutionPtr;
    bool import = dataFile("importer/import",false);
    Real importTime = 0.0;

    if(import)
    {
		importedSolutionPtr.reset( new vector_Type ( solid.displacement().map() ) );

		boost::shared_ptr< Exporter<RegionMesh<LinearTetra> > > importer;
		importer.reset(new ExporterHDF5<RegionMesh<LinearTetra> > () );

		std::string prefix = dataFile("importer/prefix","");
		std::cout <<  "Prefix: " << prefix <<"\n";
		importer->setPrefix (prefix);

		std::string postDir = dataFile("importer/prefix","");
		std::cout <<  "PostDir: " << postDir <<"\n";
		importer->setPostDir (postDir);

		importer->setMeshProcId (dFESpace->mesh(),
								 dFESpace->map().comm().MyPID() );
		//! Exporter data
		typedef ExporterData<mesh_Type>                                     IOData_Type;
		importer->addVariable (IOData_Type::VectorField, "displacement", dFESpace,
				importedSolutionPtr, static_cast<UInt> (0) );

		importTime = dataFile("importer/import_time", 0.0);
		importer->importFromTime (importTime);
		importer->closeFile();

		solid.initialize(importedSolutionPtr);
    }
    //===========================================================
    //===========================================================
    //              BOUNDARY CONDITIONS Part II
    //===========================================================
    //===========================================================

    //ADD THE NATURAL CONDITION ON THE DEFORMED CONFIGURATION
    Int flag =  dataFile ( "solid/boundary_conditions/flag", 500);
    Real FinalPressure =  dataFile ( "solid/boundary_conditions/pressure", 1000.0);

    vectorPtr_Type boundaryVectorPtr(new vector_Type ( solid.displacement().map(), Repeated ) );
    ComputeBC<solidETFESpace_Type>(localSolidMesh, solid.displacement(), boundaryVectorPtr, dETFESpace, flag);
    boost::shared_ptr<BCVector> bcVectorPtr( new BCVector (*boundaryVectorPtr, dFESpace -> dof().numTotalDof(), 0 ) );
    solidBC -> handler() -> addBC("DeformedSide", flag, Natural, Full, *bcVectorPtr, 3);
    solidBC -> handler() -> bcUpdate( *dFESpace->mesh(), dFESpace->feBd(), dFESpace->dof() );


    //===========================================================
    //===========================================================
    //              CREATE EXPORTER and SAVE INITIAL SOLUTION
    //===========================================================
    //===========================================================

    boost::shared_ptr< Exporter<RegionMesh<LinearTetra> > > exporter;
    exporter.reset ( new ExporterHDF5<RegionMesh<LinearTetra> > ( dataFile, "structure" ) );
    exporter->setPostDir ( problemFolder );
    exporter->setMeshProcId ( localSolidMesh, comm->MyPID() );
    exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, "displacement", dFESpace, solid.displacementPtr(), UInt (0) );
    exporter->postProcess ( 0 );

    //===========================================================
    //===========================================================
    //         SOLVE
    //===========================================================
    //===========================================================
    Real dt =  dataFile ( "solid/time_discretization/timestep", 0.1);
    Real endTime = dataFile ( "solid/time_discretization/endtime", 1.0);
    bool updatedTimeStep1 = false;
    bool updatedTimeStep2 = false;


    //===========================================================
    //         SOLVE PART I: Increase value of the pressure
    //===========================================================
    Real time(importTime);
    std::cout << "Starting from time: " << importTime << ", and finishing at: " << endTime << "\n";
    for ( ; time < endTime ; )
    {
    	time += dt;

        if ( comm->MyPID() == 0 )
        {
			std::cout << "\n=====================================================\n";
			std::cout << "============= TIME: " << time ;
			std::cout << "\n=====================================================\n";
        }
        ComputeBC<solidETFESpace_Type>(localSolidMesh, solid.displacement(), boundaryVectorPtr, dETFESpace, flag);
        *boundaryVectorPtr *= (FinalPressure * time);
        Real norm = boundaryVectorPtr -> norm2();
        std::cout << "Norm BC: " << norm << ", Pressure: " << FinalPressure * time << "\n";
        bcVectorPtr.reset( new BCVector (*boundaryVectorPtr, dFESpace -> dof().numTotalDof(), 0 ) );
	    solidBC -> handler() -> modifyBC(flag, *bcVectorPtr);

    	solid.iterate ( solidBC -> handler() );

        exporter->postProcess ( time );
    }

    //===========================================================
    //         SOLVE PART II: Iterate till convergence
    //===========================================================
    Real tol = 1e-8;
    Real res = 2*tol;
    int iteration(0);
    while( res > tol )
    {
    	iteration++;
    	time += dt;

        if ( comm->MyPID() == 0 )
        {
			std::cout << "\n=====================================================\n";
			std::cout << "============= TIME: " << time << ", Iteration: " << iteration;
			std::cout << "\n=====================================================\n";
        }
    	VectorEpetra disp_tn(solid.displacement());

        ComputeBC<solidETFESpace_Type>(localSolidMesh, solid.displacement(), boundaryVectorPtr, dETFESpace, flag);
        *boundaryVectorPtr *= (FinalPressure);
        Real norm = boundaryVectorPtr -> norm2();
        std::cout << "Norm BC: " << norm << ", Pressure: " << FinalPressure << "\n";
        bcVectorPtr.reset( new BCVector (*boundaryVectorPtr, dFESpace -> dof().numTotalDof(), 0 ) );
	    solidBC -> handler() -> modifyBC(flag, *bcVectorPtr);

    	solid.iterate ( solidBC -> handler() );

        exporter->postProcess ( time );

        disp_tn -= solid.displacement();

        res = disp_tn.norm2();

        if ( comm->MyPID() == 0 )
        {
			std::cout << "\n============= Residual: " << res ;
			std::cout << "\n=====================================================\n";
        }
    }

    //===========================================================
    //===========================================================
    //              CLOSE EXPORTER
    //===========================================================
    //===========================================================
    exporter -> closeFile();


#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return 0;
}



template<typename space> void ComputeBC( const boost::shared_ptr<RegionMesh<LinearTetra> > localMesh,
						const VectorEpetra& disp,
						boost::shared_ptr<VectorEpetra> bcVectorPtr,
						const boost::shared_ptr< space > dETFESpace,
						int bdFlag)
{

	*bcVectorPtr *= 0.0;

    MatrixSmall<3,3> Id;
    Id(0,0) = 1.; Id(0,1) = 0., Id(0,2) = 0.;
    Id(1,0) = 0.; Id(1,1) = 1., Id(1,2) = 0.;
    Id(2,0) = 0.; Id(2,1) = 0., Id(2,2) = 1.;

	boost::shared_ptr<VectorEpetra> intergral( new VectorEpetra( disp.map() ) );

	{
	  	using namespace ExpressionAssembly;

	  	auto I = value(Id);
	  	auto Grad_u = grad( dETFESpace, disp, 0);
	  	auto F =  Grad_u + I;
	  	auto FmT = minusT(F);
	  	auto J = det(F);

    	QuadratureBoundary myBDQR (buildTetraBDQR (quadRuleTria7pt) );

        integrate ( boundary ( localMesh, bdFlag),
        		    myBDQR,
        		    dETFESpace,
        		    J * dot( FmT * Nface,  phi_i)
                  ) >> bcVectorPtr;

        bcVectorPtr -> globalAssemble();


	}
}
