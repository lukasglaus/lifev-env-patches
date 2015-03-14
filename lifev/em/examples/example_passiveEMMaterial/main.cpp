#include <lifev/core/LifeV.hpp>
//include <lifev/electrophysiology/solver/ElectroETAMonodomainSolver.hpp>
//#include <lifev/electrophysiology/solver/IonicModels/IonicMinimalModel.hpp>
#include <lifev/electrophysiology/util/HeartUtility.hpp>
#include <lifev/core/mesh/MeshLoadingUtility.hpp>


#include <lifev/structure/solver/StructuralConstitutiveLawData.hpp>

//#include <lifev/structure/solver/StructuralConstitutiveLaw.hpp>
//#include <lifev/structure/solver/StructuralOperator.hpp>

#include <lifev/em/solver/mechanics/EMStructuralOperator.hpp>
#include <lifev/em/solver/mechanics/EMStructuralConstitutiveLaw.hpp>
//#include <lifev/em/solver/EMETAFunctors.hpp>

#include <lifev/core/filter/ExporterEnsight.hpp>
#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5.hpp>
#endif
#include <lifev/core/filter/ExporterEmpty.hpp>

//#include <lifev/eta/fem/ETFESpace.hpp>
//#include <lifev/eta/expression/Integrate.hpp>

//#include <lifev/em/solver/mechanics/materials/EMMaterial.hpp>


#include <lifev/bc_interface/3D/bc/BCInterface3D.hpp>
//#include <lifev/structure/solver/NeoHookeanMaterialNonLinear.hpp>

using namespace LifeV;


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
    
    
    //********************************************//
    // Declare comm
    //********************************************//
    
#ifdef HAVE_MPI
    MPI_Init ( &argc, &argv );
#endif
    
    boost::shared_ptr<Epetra_Comm>  comm ( new Epetra_MpiComm (MPI_COMM_WORLD) );
    if ( comm->MyPID() == 0 )
    {
        cout << "% using MPI" << endl;
    }
    
    
    //********************************************//
    // Read data file and create output folder
    //********************************************//
    
    GetPot command_line (argc, argv);
    const std::string data_file_name = command_line.follow ("data", 2, "-f", "--file");
    GetPot dataFile (data_file_name);
    std::string problemFolder = EMUtility::createOutputFolder (command_line, *comm);

        
    //********************************************//
    // Setup material data
    //********************************************//
    
    EMData emdata;
    emdata.setup (dataFile);

        
    //********************************************//
    // Load mesh
    //********************************************//
    
    if ( comm->MyPID() == 0 )
    {
        std::cout << "Reading Mesh Name and Path...\n";
    }

    std::string meshName = dataFile ( "solid/space_discretization/mesh_file", "" );
    std::string meshPath = dataFile ( "solid/space_discretization/mesh_dir", "./" );

    meshPtr_Type localSolidMesh ( new mesh_Type ( comm ) );
    meshPtr_Type fullSolidMesh ( new mesh_Type ( comm ) );

    MeshUtility::loadMesh (localSolidMesh, fullSolidMesh, meshName, meshPath);

    if ( comm->MyPID() == 0 )
    {
        std::cout << " Done!" << std::endl;
    }
        

    //********************************************//
    // FE-Spaces
    //********************************************//

    // Define the finite element space for bc and exporter
    // and the ET finite element space for assembly
    if ( comm->MyPID() == 0 )
    {
        std::cout << "setup spaces ... ";
    }

    std::string dOrder =  dataFile ( "solid/space_discretization/order", "P1");
    solidFESpacePtr_Type dFESpace ( new solidFESpace_Type (localSolidMesh, dOrder, 3, comm) );
    solidETFESpacePtr_Type dETFESpace ( new solidETFESpace_Type (localSolidMesh, & (dFESpace->refFE() ), & (dFESpace->fe().geoMap() ), comm) );

    if ( comm->MyPID() == 0 )
    {
        std::cout << " Done!" << std::endl;
    }


    //********************************************//
    // Boundary Conditions
    //********************************************//
    
    if ( comm->MyPID() == 0 )
    {
        std::cout << "setup bc ... ";
    }

    bcInterfacePtr_Type solidBC ( new bcInterface_Type() );
    solidBC->createHandler();
    solidBC->fillHandler ( data_file_name, "solid" );
    solidBC->handler()->bcUpdate ( *dFESpace->mesh(), dFESpace->feBd(), dFESpace->dof() );

    if ( comm->MyPID() == 0 )
    {
        std::cout << " Done!" << std::endl;
    }
    
    
    //********************************************//
    // Solid Mechanics
    //********************************************//
    
    //Setup data of constitutive law
    if ( comm->MyPID() == 0 )
    {
        std::cout << "setup constitutive law data ... ";
    }

    constitutiveLawPtr_Type dataStructure ( new constitutiveLaw_Type() );
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

    EMStructuralOperator< mesh_Type > solid;

    // solidBC is copied inside the StructuralOperator
    // any changes to it don't affect the M_BCh inside the solver
    // after a change in solidBC add it in the solver!
    solid.setup ( dataStructure, dFESpace, dETFESpace, solidBC -> handler(), comm);
    solid.setDataFromGetPot (dataFile);
    solid.EMMaterial()->setParameters(emdata);

    // load fibers and sheets fields
    std::string fiberFileName  =  dataFile ( "solid/space_discretization/fiber_name", "fiber");
    std::string sheetFileName  =  dataFile ( "solid/space_discretization/sheet_name", "sheet");
    std::string fiberFieldName =  dataFile ( "solid/space_discretization/fiber_fieldname", "fiber");
    std::string sheetFieldName =  dataFile ( "solid/space_discretization/sheet_fieldname", "sheet");
    std::string fiberDir       =  dataFile ( "solid/space_discretization/fiber_dir", "./");
    std::string sheetDir       =  dataFile ( "solid/space_discretization/sheet_dir", "./");

    if ( comm->MyPID() == 0 )
    {
        std::cout << "Importing fibers field\n";
        std::cout << "Fibers file name: " << fiberFileName  << "\n";
        std::cout << "Fibers field name: " << fiberFieldName << "\n";
        std::cout << "Fibers dir name: " << fiberDir       << "\n";
    }
    ElectrophysiologyUtility::importVectorField (  solid.EMMaterial()->fiberVectorPtr(),
                                                   fiberFileName,
                                                   fiberFieldName,
                                                   localSolidMesh,
                                                   fiberDir,
                                                   dOrder );
    if ( comm->MyPID() == 0 )
    {
        std::cout << "Importing sheets field\n";
        std::cout << "Sheets file name: " << fiberFileName  << "\n";
        std::cout << "Sheets field name: " << fiberFieldName << "\n";
        std::cout << "Sheets dir name: " << fiberDir       << "\n";
    }
    ElectrophysiologyUtility::importVectorField (  solid.EMMaterial()->sheetVectorPtr(),
                                                   sheetFileName,
                                                   sheetFieldName,
                                                   localSolidMesh,
                                                   sheetDir,
                                                   dOrder );

//    solid.EMMaterial() -> setupFiberVector( 1.0, 0.0, 0.0);
//    solid.EMMaterial() -> setupSheetVector( 0.0, 1.0, 0.0);

    solid.setNewtonParameters(dataFile);
    solid.buildSystem (1.0);

    if ( comm->MyPID() == 0 )
    {
        std::cout << " Done!" << std::endl;
    }

    
    //********************************************//
    // Create exporter and initial solution
    //********************************************//
    
    // Exporter displacement
    boost::shared_ptr< Exporter<RegionMesh<LinearTetra> > > exporter;
    exporter.reset ( new ExporterHDF5<RegionMesh<LinearTetra> > ( dataFile, "structure" ) );
    exporter->setPostDir ( problemFolder );
    exporter->setMeshProcId ( localSolidMesh, comm->MyPID() );
    exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, "displacement", dFESpace, solid.displacementPtr(), UInt (0) );
    exporter->postProcess ( 0 );

    // Exporter fibers and sheets
    boost::shared_ptr< Exporter<RegionMesh<LinearTetra> > > fibers_exporter;
    fibers_exporter.reset ( new ExporterHDF5<RegionMesh<LinearTetra> > ( dataFile, "fibers" ) );
    fibers_exporter->setPostDir ( problemFolder );
    fibers_exporter->setMeshProcId ( localSolidMesh, comm->MyPID() );
    fibers_exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, "fibers", dFESpace, solid.EMMaterial()->fiberVectorPtr(), UInt (0) );
    fibers_exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, "sheets", dFESpace, solid.EMMaterial()->sheetVectorPtr(), UInt (0) );
    fibers_exporter->postProcess ( 0 );

    std::cout << "Displacement vector size: " << solid.displacementPtr()->size() << std::endl;
    std::cout << "Fiber vector size: " << solid.EMMaterial()->fiberVectorPtr()->size() << std::endl;

    
    //********************************************//
    // Solve
    //********************************************//
    
    Real dt =  dataFile ( "solid/time_discretization/timestep", 0.1);
    Real endTime =  dataFile ( "solid/time_discretization/endtime", 1.0);
    ID LVFlag =  dataFile ( "solid/boundary_conditions/LV_flag", 0);
    Real LVPreloadPressure =  dataFile ( "solid/boundary_conditions/LV_preload_pressure", 0.0);
    bool deformedPressure =  dataFile ( "solid/boundary_conditions/deformed_pressure", 1 );

    if ( deformedPressure)
    {
        std::cout << "Setting pressure in the deformed configuration\n";
    }
    else
    {
        std::cout << "Setting pressure in the reference configuration\n";
    }
    
    solid.setBCFlag( LVFlag );

    for (Real time (0.0); time < endTime;)
    {
        time += dt;
        std::cout << "\nTime: " << time << std::endl;
        solid.data() -> dataTime() -> updateTime();
        std::cout << "----- Preload Time: " << time << std::endl;

        solidBC -> updatePhysicalSolverVariables();

//        solid.bcH() = solidBC -> handler();

        solid.setLVPressureBC( -time*LVPreloadPressure ); // RV ?

        solid.iterate ( solidBC -> handler() , deformedPressure );

        exporter->postProcess ( time );
    }

    
    //********************************************//
    // Close exporters
    //********************************************//
    
    exporter -> closeFile();
    fibers_exporter->closeFile();

#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return 0;
}
