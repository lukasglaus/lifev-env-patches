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
#include <lifev/core/filter/ExporterVTK.hpp>

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

    EMData emdata;
    emdata.setup (dataFile);

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
    solid.EMMaterial()->setParameters(emdata);

    solid.EMMaterial()->showMaterialParameters();
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


    std::cout << "\nNumber of DOF: " << solid.displacement().size() << std::endl;

    //===========================================================
    //===========================================================
    //              IMPORTER
    //===========================================================
    //===========================================================
    vectorPtr_Type importedSolutionPtr;
    bool import = dataFile("importer/import",false);
    Real importTime = 0.0;


	importedSolutionPtr.reset( new vector_Type ( solid.displacement().map() ) );

	boost::shared_ptr< Exporter<RegionMesh<LinearTetra> > > importer;
	importer.reset(new ExporterHDF5<RegionMesh<LinearTetra> > () );

	std::string prefix = dataFile("importer/prefix","");
	std::cout <<  "Prefix: " << prefix <<"\n";
	importer->setPrefix (prefix);

	std::string postDir = dataFile("importer/postDir","");
	std::cout <<  "PostDir: " << postDir <<"\n";
//	importer->setPostDir (postDir);
//	std::cout <<  "PostDir: " << postDir <<"\n";
//	importer->setMeshProcId (dFESpace->mesh(),
//							 dFESpace->map().comm().MyPID() );
	//! Exporter data

	std::cout <<  "\nAdding variable\n";
	typedef ExporterData<mesh_Type>                                     IOData_Type;
	importer->addVariable (IOData_Type::VectorField, "displacement", dFESpace,
			importedSolutionPtr, static_cast<UInt> (0) );

	std::cout <<  "\nImporting from time ";
	importTime = dataFile("importer/import_time", 0.0);
	std::cout <<  importTime << "\n";
	std::cout <<  "\nReading: ";

	importer->importFromTime (importTime);
	std::cout <<  "\nDone!";

	importer->closeFile();

	std::cout <<  "\nDone!";

	std::cout << "\nNorm of the solution: " << importedSolutionPtr->norm2() << std::endl;

    //*************************************************************//
    // We save the potential field just computed on a file
    // using VTK.
    //*************************************************************//
	std::cout <<  "\nVTK Export: ";
    ExporterVTK< mesh_Type > exporterVTK;
	std::cout <<  "\nSetting MeshProcId ";
    exporterVTK.setMeshProcId ( localSolidMesh, comm -> MyPID() );
	std::cout <<  "\nSetting output folder";
    exporterVTK.setPostDir (problemFolder);
	std::cout <<  "\nSetting output filename";
    exporterVTK.setPrefix ("structure");
	std::cout <<  "\nAdding variable to exporter";
    exporterVTK.addVariable ( ExporterData<mesh_Type>::VectorField,  "displacement", solid.dispFESpacePtr(), importedSolutionPtr, UInt (0) );
	std::cout <<  "\nSaving.";
    exporterVTK.postProcess (0);
	std::cout <<  " Done!\n";
    exporterVTK.closeFile();



#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return 0;
}



