#include <lifev/core/LifeV.hpp>
#include <lifev/electrophysiology/solver/ElectroETAMonodomainSolver.hpp>
#include <lifev/electrophysiology/solver/IonicModels/IonicGoldbeter.hpp>
#include <lifev/structure/solver/StructuralConstitutiveLawData.hpp>
#include <lifev/structure/solver/StructuralConstitutiveLaw.hpp>
#include <lifev/structure/solver/StructuralOperator.hpp>
#include <lifev/structure/solver/NeoHookeanActivatedMaterial.hpp>
//#include <lifev/structure/solver/GeneralizedActiveHolzapfelOgdenMaterial.hpp>
#include <lifev/em/solver/EMETAFunctors.hpp>

#include <lifev/core/filter/ExporterEnsight.hpp>
#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5.hpp>
#endif
#include <lifev/core/filter/ExporterEmpty.hpp>

#include <lifev/eta/fem/ETFESpace.hpp>
#include <lifev/eta/expression/Integrate.hpp>
#include <lifev/bc_interface/3D/bc/BCInterface3D.hpp>

using namespace LifeV;


Real initialSphereOnCylinder (const Real& /*t*/, const Real&  X, const Real& Y, const Real& Z, const ID& /*i*/)
{

    double r = std::sqrt (pow (X, 2) + pow (Y - 60, 2) + pow (Z - 12, 2) );
    double auxexp = 1.0 - 1.0 / (1.0 + exp (-50.0 * (r - 10) ) );

    return 0.1 + 3.5 * auxexp;
}

Real initialSphereOnCell (const Real& /*t*/, const Real&  X, const Real& Y, const Real& Z, const ID& /*i*/)
{

    double r = std::sqrt (pow (X - 42, 2) + pow (Y - 46, 2) + pow (Z - 7, 2) );
    double auxexp = 1.0 - 1.0 / (1.0 + exp (-50.0 * (r - 6.5) ) );

    return 0.1 + 3.5 * auxexp;
}


Real initialStimulus (const Real& /*t*/, const Real&  /*X*/, const Real& Y, const Real& /*Z*/, const ID& /*i*/)
{
    if ( Y == 0 )
    {
        return 3.5;
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
    typedef boost::function < Real (const Real& /*t*/,
                                    const Real &   x,
                                    const Real &   y,
                                    const Real& /*z*/,
                                    const ID&   /*i*/ ) >   function_Type;
    typedef IonicGoldbeter ionicModel_Type;
    typedef boost::shared_ptr< ionicModel_Type >  ionicModelPtr_Type;

    typedef ElectroETAMonodomainSolver< mesh_Type, ionicModel_Type >        monodomainSolver_Type;
    typedef boost::shared_ptr< monodomainSolver_Type >  monodomainSolverPtr_Type;
    typedef VectorEpetra              vector_Type;
    typedef boost::shared_ptr<vector_Type> vectorPtr_Type;

    typedef MatrixEpetra<Real> matrix_Type;
    typedef boost::shared_ptr<matrix_Type> matrixPtr_Type;
    typedef BCHandler                                          bc_Type;
    typedef boost::shared_ptr< bc_Type >                       bcPtr_Type;
    typedef EMStructuralOperator< RegionMesh<LinearTetra> >     physicalSolver_Type;
    typedef BCInterface3D< bc_Type, physicalSolver_Type >              bcInterface_Type;
    typedef boost::shared_ptr< bcInterface_Type >              bcInterfacePtr_Type;

#ifdef HAVE_MPI
    MPI_Init ( &argc, &argv );
#endif



    //===========================================================
    //                ELECTROPHYSIOLOGY
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
    // mesh name and the mesh path               //
    //********************************************//
    if ( comm->MyPID() == 0 )
    {
        std::cout << "Reading Mesh Name and Path...";
    }

    std::string meshName = parameterList.get ("mesh_name", "lid16.mesh");
    std::string meshPath = parameterList.get ("mesh_path", "./");

    meshPtr_Type mesh ( new mesh_Type ( comm ) );
    meshPtr_Type fullMesh ( new mesh_Type ( comm ) );
    MeshUtility::fillWithFullMesh (mesh, fullMesh, meshName, meshPath);
    if ( comm->MyPID() == 0 )
    {
        std::cout << " Done!" << endl;
    }

    //********************************************//
    // GetPot datafile to setup the preconditioner//
    //********************************************//
    GetPot command_line (argc, argv);
    const string data_file_name = command_line.follow ("data", 2, "-f", "--file");
    GetPot dataFile (data_file_name);

    //********************************************//
    // Creates a new model object.  The           //
    // model input are the parameters. Pass  the  //
    // parameter list in the constructor          //
    //********************************************//
    if ( comm->MyPID() == 0 )
    {
        std::cout << "Building Constructor for Intracellular Calcium with parameters ... ";
    }
    ionicModelPtr_Type  ionicModel ( new ionicModel_Type() );
    if ( comm->MyPID() == 0 )
    {
        std::cout << " Done!" << endl;
    }

    //********************************************//
    // set up the monodomain solver               //
    //********************************************//
    if ( comm->MyPID() == 0 )
    {
        std::cout << "Building Monodomain Solvers... ";
    }

    monodomainSolverPtr_Type monodomain ( new monodomainSolver_Type ( meshName, meshPath, dataFile, ionicModel ) );
    if ( comm->MyPID() == 0 )
    {
        std::cout << " Splitting solver done... ";
    }

    monodomain -> setInitialConditions();
    function_Type f = &initialSphereOnCell;
    monodomain -> setPotentialFromFunction (f);

    for (int i (0); i < ionicModel -> Size(); i++ )
    {
        std::cout << "Norm Inf variable " << i  << " = " <<  (  * ( monodomain -> globalSolution().at (i) ) ).normInf() << std::endl;
    }
    monodomain -> setParameters ( parameterList );

    //********************************************//
    // Creating exporters to save the solution    //
    //********************************************//
    ExporterHDF5< RegionMesh <LinearTetra> > exp;

    if ( comm->MyPID() == 0 )
    {
        cout << "\nExporter setup:  " ;
    }

    monodomain -> setupExporter ( exp, parameterList.get ("OutputFile", "output") );

    if ( comm->MyPID() == 0 )
    {
        cout << "\nExport at 0:  " ;
    }

    monodomain -> exportSolution ( exp, 0.0 );

    if ( comm->MyPID() == 0 )
    {
        cout << "\nsolve system:  " ;
    }

    //===========================================================
    //                SOLID MECHANICS
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
        std::cout << "\nparameters" << std::endl;
        std::cout << "\ninitialization constitutive law" << std::endl;
    }

    boost::shared_ptr<StructuralConstitutiveLawData> dataStructure (new StructuralConstitutiveLawData( ) );
    dataStructure->setup (dataFile);

    if ( comm->MyPID() == 0 )
    {
        std::cout << "setup spaces" << std::endl;
    }
    std::string dOrder =  dataFile ( "solid/space_discretization/order", "P1");
    solidFESpacePtr_Type dFESpace ( new solidFESpace_Type (monodomain -> localMeshPtr(), dOrder, 3, comm) );
    solidFESpacePtr_Type aFESpace ( new solidFESpace_Type (monodomain -> localMeshPtr(), dOrder, 1, comm) );
    solidETFESpacePtr_Type dETFESpace ( new solidETFESpace_Type (monodomain -> localMeshPtr(), & (dFESpace->refFE() ), & (dFESpace->fe().geoMap() ), comm) );
    scalarETFESpacePtr_Type aETFESpace ( new scalarETFESpace_Type (monodomain -> localMeshPtr(), & (aFESpace->refFE() ), & (aFESpace->fe().geoMap() ), comm) );

    if ( comm->MyPID() == 0 )
    {
        std::cout << "\nsetup boundary conditions" << std::endl;
    }

    //===========================================================
    //            BOUNDARY CONDITIONS
    //===========================================================

    //! Setting solid BCs from data file
    bcInterfacePtr_Type solidBC ( new bcInterface_Type() );
    solidBC->createHandler();
    solidBC->fillHandler ( data_file_name, "solid" );

    if ( comm->MyPID() == 0 )
    {
        std::cout << "\nsetup structural operator" << std::endl;
    }

    //! 1. Constructor of the structuralSolver
    EMStructuralOperator< RegionMesh<LinearTetra> > solid;
    solid.setup (dataStructure,
                 dFESpace,
                 dETFESpace,
                 solidBC -> handler(),
                 comm);
    if ( comm->MyPID() == 0 )
    {
        std::cout << "\ninitial guess" << std::endl;
    }

    solid.setDataFromGetPot (dataFile);


    //===========================================================
    //                FIBERS
    //===========================================================
    vectorPtr_Type fibers ( new vector_Type ( dFESpace -> map() ) );
    vectorPtr_Type emDisp;

    if ( comm->MyPID() == 0 )
    {
        std::cout << "\nread fibers" << std::endl;
    }

    ElectrophysiologyUtility::importFibers (fibers, parameterList.get ("fibers_file", ""), monodomain -> localMeshPtr() );


    if ( comm->MyPID() == 0 )
    {
        std::cout << "\nset fibers" << std::endl;
    }

    monodomain -> setFiberPtr ( fibers );
    monodomain -> exportFiberDirection();
    emDisp = solid.displacementPtr();
    solid.material() -> setFiberVector ( *fibers );

    //********************************************//
    // Global Electrophys matrix: mass + stiffness//
    //********************************************//
    if ( comm->MyPID() == 0 )
    {
        cout << "\nSetup operators:  dt = " << monodomain -> timeStep() << "\n" ;
    }
    // monodomain -> setDisplacementPtr( emDisp );
    monodomain -> setupLumpedMassMatrix();
    monodomain -> setupStiffnessMatrix();
    monodomain -> setupGlobalMatrix();

    if ( comm->MyPID() == 0 )
    {
        cout << "Done! \n" ;
    }




    if ( comm->MyPID() == 0 )
    {
        std::cout << "\nset gammaf and fibers" << std::endl;
    }

    vectorPtr_Type gammaf ( new vector_Type ( monodomain -> globalSolution().at (0) -> map() ) );

    *gammaf *= 0.0;

    //obs
    solid.material() -> setGammaf ( *gammaf );

    if ( comm->MyPID() == 0 )
    {
        std::cout << "\nbuild solid system" << std::endl;
    }

    solid.buildSystem (1.0);
    vectorPtr_Type rhs (new vector_Type (solid.displacement(), Unique) );
    vectorPtr_Type initialDisplacement (new vector_Type (solid.displacement(), Unique) );
    solid.initialize ( initialDisplacement );


    //MPI_Barrier (MPI_COMM_WORLD);

    if ( comm->MyPID() == 0 )
    {
        std::cout << "\nsetup solid exporter" << std::endl;
    }

    boost::shared_ptr< Exporter<RegionMesh<LinearTetra> > > exporter;
    exporter.reset ( new ExporterHDF5<RegionMesh<LinearTetra> > ( dataFile, "structure" ) );

    exporter->setPostDir ( "./" );
    exporter->setMeshProcId ( monodomain -> localMeshPtr(), comm->MyPID() );

    vectorPtr_Type solidDisp ( new vector_Type (solid.displacement(), exporter->mapType() ) );
    exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField,
                            "displacement", dFESpace, solidDisp, UInt (0) );
    exporter->postProcess ( 0 );


    //================================================================//
    //     COUPLING SOLVER                    //
    //================================================================//
    ExporterHDF5< RegionMesh <LinearTetra> > expGammaf;
    expGammaf.setMeshProcId (monodomain -> localMeshPtr(), comm->MyPID() );
    expGammaf.setPrefix ("gammaf");

    //obs
    matrixPtr_Type mass (new matrix_Type ( monodomain -> massMatrixPtr() -> map() ) ) ;
    {
        using namespace ExpressionAssembly;

        integrate (elements (monodomain -> localMeshPtr() ), monodomain -> feSpacePtr() -> qr(), monodomain -> ETFESpacePtr(),
                   monodomain -> ETFESpacePtr(), phi_i * phi_j) >> mass;

    }
    mass -> globalAssemble();





    vectorPtr_Type rhsActivation ( new vector_Type ( *gammaf ) );
    *rhsActivation *= 0;

    if ( comm->MyPID() == 0 )
    {
        std::cout << "\nSolve system" << std::endl;
    }

    //==================================================================//
    //    SETUP LINEAR SOLVER & ACTIVATION                      //
    //==================================================================//

    if ( comm->MyPID() == 0 )
    {
        std::cout << "\nset up linear solver... Does it work???" << std::endl;
    }
    typedef LinearSolver linearSolver_Type;
    typedef boost::shared_ptr<LinearSolver> linearSolverPtr_Type;
    typedef LifeV::Preconditioner basePrec_Type;
    typedef boost::shared_ptr<basePrec_Type> basePrecPtr_Type;
    typedef LifeV::PreconditionerIfpack prec_Type;
    typedef boost::shared_ptr<prec_Type> precPtr_Type;


    prec_Type* precRawPtr;
    basePrecPtr_Type precPtr;
    precRawPtr = new prec_Type;
    precRawPtr->setDataFromGetPot (dataFile, "prec");
    precPtr.reset (precRawPtr);

    if ( comm->MyPID() == 0 )
    {
        std::cout << "\nprec done!!!!" << std::endl;
    }

    Teuchos::RCP < Teuchos::ParameterList > solverParamList = Teuchos::rcp (
                                                                  new Teuchos::ParameterList);

    std::string xmlpath = dataFile ("electrophysiology/monodomain_xml_path",
                                    "./");
    std::string xmlfile = dataFile ("electrophysiology/monodomain_xml_file",
                                    "MonodomainSolverParamList.xml");

    solverParamList = Teuchos::getParametersFromXmlFile (xmlpath + xmlfile);

    if ( comm->MyPID() == 0 )
    {
        std::cout << "\nreading file done!!!!" << std::endl;
    }

    linearSolver_Type linearSolver;
    linearSolver.setCommunicator ( comm );
    linearSolver.setParameters ( *solverParamList );
    linearSolver.setPreconditioner ( precPtr );
    linearSolver.setOperator ( mass );
    //  linearSolver.setOperator( monodomain -> massMatrixPtr() );

    if ( comm->MyPID() == 0 )
    {
        std::cout << "\nIt does!!!!" << std::endl;
    }




    //Initial condition for gammaf
    *gammaf = -0.015;

    vectorPtr_Type tmpRhsActivation ( new vector_Type ( rhsActivation -> map(), Repeated ) );
    expGammaf.addVariable (ExporterData<mesh_Type>::ScalarField, "gammaf",
                           monodomain -> feSpacePtr(), gammaf, UInt (0) );
    expGammaf.addVariable (ExporterData<mesh_Type>::ScalarField, "rhs",
                           monodomain -> feSpacePtr(), rhsActivation, UInt (0) );
    expGammaf.postProcess (0.0);



    //===========================================================
    //===========================================================
    //                TIME LOOP
    //===========================================================
    //===========================================================
    int k (0);
    Real emdt = parameterList.get ("emdt", 1.0);
    int iter ( (emdt / monodomain -> timeStep() ) );
    Real saveStep = parameterList.get ("save_step", 1.0);
    int saveIter ( (saveStep / monodomain -> timeStep() ) );

    for ( Real t (0.0); t < monodomain -> endTime(); )
    {
        t = t + monodomain -> timeStep();
        k++;
        monodomain -> solveOneSplittingStep ( exp, t );
        *tmpRhsActivation *= 0;
        {
            using namespace ExpressionAssembly;

            BOOST_AUTO_TPL (Ca, value ( aETFESpace, * (monodomain -> globalSolution().at (0)  ) ) );
            BOOST_AUTO_TPL (Gammaf, value ( aETFESpace, *gammaf ) );
            BOOST_AUTO_TPL (activationEquation, value (-0.5) *Ca + value (-2.5) *Gammaf); //Laadhari et al. IJNME 2013

            integrate ( elements ( monodomain -> localMeshPtr() ),
                        monodomain -> feSpacePtr() -> qr() ,
                        monodomain -> ETFESpacePtr(),
                        activationEquation * phi_i
                      ) >> tmpRhsActivation;

        }
        *rhsActivation *= 0;
        *rhsActivation = ( * (mass) * ( *gammaf ) );
        *rhsActivation += ( monodomain -> timeStep() * *tmpRhsActivation );

        linearSolver.setRightHandSide (rhsActivation);
        linearSolver.solve (gammaf);

        if ( k % iter == 0)
        {
            solid.material() -> setGammaf ( *gammaf );
            solid.iterate ( solidBC -> handler() );
            *solidDisp = solid.displacement();
        }

        //if ( k % saveIter == 0){
        monodomain -> exportSolution (exp, t);
        expGammaf.postProcess (t);
        exporter->postProcess ( t );
        //}
    }
    exp.closeFile();
    expGammaf.closeFile();

    exporter -> closeFile();
    if ( comm->MyPID() == 0 )
    {
        std::cout << "Active strain example: Passed!" << std::endl;
    }

#undef Gammaf
#undef Ca
#undef activationEquation

#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return 0;
}
