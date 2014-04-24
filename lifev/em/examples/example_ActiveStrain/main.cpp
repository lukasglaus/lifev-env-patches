#include <lifev/core/LifeV.hpp>
#include <lifev/electrophysiology/solver/ElectroETAMonodomainSolver.hpp>
#include <lifev/electrophysiology/solver/IonicModels/IonicMinimalModel.hpp>

#include <lifev/structure/solver/StructuralConstitutiveLawData.hpp>

#include <lifev/structure/solver/StructuralConstitutiveLaw.hpp>
#include <lifev/structure/solver/StructuralOperator.hpp>
#include <lifev/structure/solver/NeoHookeanActivatedMaterial.hpp>
#include <lifev/em/solver/EMETAFunctors.hpp>

#include <lifev/core/filter/ExporterEnsight.hpp>
#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5.hpp>
#endif
#include <lifev/core/filter/ExporterEmpty.hpp>

#include <lifev/eta/fem/ETFESpace.hpp>
#include <lifev/eta/expression/Integrate.hpp>

#include <lifev/core/interpolation/RBFhtp.hpp>
#include <lifev/core/interpolation/RBFhtpVectorial.hpp>

#include <lifev/core/interpolation/RBFlocallyRescaledVectorial.hpp>
#include <lifev/core/interpolation/RBFlocallyRescaledScalar.hpp>
#include <lifev/core/interpolation/RBFrescaledVectorial.hpp>
#include <lifev/core/interpolation/RBFrescaledScalar.hpp>


#include <lifev/bc_interface/3D/bc/BCInterface3D.hpp>

using namespace LifeV;


Real bcZero (const Real& /*t*/, const Real&  /*X*/, const Real& /*Y*/, const Real& /*Z*/, const ID& /*i*/)
{
    return  0.;
}
Real d0 (const Real& /*t*/, const Real&  /*X*/, const Real& /*Y*/, const Real& /*Z*/, const ID& /*i*/)
{
    return  0.;
}

Real initialVlid (const Real& /*t*/, const Real&  X, const Real& /*Y*/, const Real& /*Z*/, const ID& /*i*/)
{
    if ( X < 0.1 )
    {
        return 1.0;
    }
    else
    {
        return  0.;
    }
}

Real fiberRotation (const Real& /*t*/, const Real&  X, const Real& Y, const Real& Z, const ID& i)
{
    Real R = std::sqrt ( X * X + Y * Y);
    //Real teta = std::atan( Y / X );
    Real fz = 0.0;
    Real fx =  Y / R;
    Real fy = - X / R;
    Real sx = X / R;
    Real sy = Y / R;
    Real m = -1.9040;
    Real q = 3.5224;
    Real theta = m * R + q;

    //  f01a f001*cos(teta)+f001*s01^2*(1-cos(teta))+s01*s02*f002*(1-cos(teta))
    //  f02a s01*s02*f001*(1-cos(teta))+f002*cos(teta)+f002*s02^2*(1-cos(teta))
    //  f03a s01*f002*sin(teta)-s02*f001*sin(teta)

    switch (i)
    {
        case 0:
            return  fx * std::cos (theta) + fx * sx * sx * ( 1.0  - std::cos (theta) ) + sx * sy * fy * ( 1.0  - std::cos (theta) );
            break;
        case 1:
            return sx * sy * fy *  ( 1.0  - std::cos (theta) ) + fy * std::cos (theta) + fy * sy * sy * ( 1.0  - std::cos (theta) ) ;
            break;
        case 2:
            return sx * fy * std::sin (theta) - sy * fx * std::sin (theta);
            break;
        default:
            ERROR_MSG ("This entry is not allowed: ud_functions.hpp");
            return 0.;
            break;
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
    typedef IonicMinimalModel                   ionicModel_Type;
    typedef boost::shared_ptr< ionicModel_Type >  ionicModelPtr_Type;

    typedef ElectroETAMonodomainSolver< mesh_Type, ionicModel_Type >        monodomainSolver_Type;
    typedef boost::shared_ptr< monodomainSolver_Type >  monodomainSolverPtr_Type;
    typedef VectorEpetra                vector_Type;
    typedef boost::shared_ptr<vector_Type> vectorPtr_Type;

    typedef MatrixEpetra<Real> matrix_Type;
    typedef boost::shared_ptr<matrix_Type> matrixPtr_Type;

    typedef BCHandler                                          bc_Type;
    typedef boost::shared_ptr< bc_Type >                       bcPtr_Type;
    typedef  EMStructuralOperator< RegionMesh<LinearTetra> >      physicalSolver_Type;
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

    //    meshPtr_Type mesh ( new mesh_Type ( comm ) );
    //    meshPtr_Type fullMesh ( new mesh_Type ( comm ) );
    //    MeshUtility::fillWithFullMesh (mesh, fullMesh, meshName, meshPath);
    if ( comm->MyPID() == 0 )
    {
        std::cout << " Done!" << endl;
    }

    //********************************************//
    // We need the GetPot datafile for to setup   //
    // the preconditioner.                        //
    //********************************************//
    GetPot command_line (argc, argv);
    const string data_file_name = command_line.follow ("data", 2, "-f", "--file");
    GetPot dataFile (data_file_name);

    //********************************************//
    // Creates a new model object representing the//
    // model from Aliev and Panfilov 1996.  The   //
    // model input are the parameters. Pass  the  //
    // parameter list in the constructor          //
    //********************************************//
    if ( comm->MyPID() == 0 )
    {
        std::cout << "Building Constructor for Minimal Model with parameters ... ";
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

    //    ionicModel -> initialize( monodomain -> globalSolution() );
    monodomain -> setInitialConditions();

    for (int i (0); i < ionicModel -> Size(); i++ )
    {
        std::cout << "Norm Inf variable " << i  << " = " <<  (  * ( monodomain -> globalSolution().at (i) ) ).normInf() << std::endl;
    }

    //HeartUtility::setValueOnBoundary( *(monodomain -> potentialPtr() ), monodomain -> fullMeshPtr(), 1.0, 30 );
    //    function_Type Vlid = &initialVlid;
    //    monodomain -> setPotentialFromFunction( Vlid );
    ElectrophysiologyUtility::setValueOnBoundary ( * (monodomain -> potentialPtr() ), monodomain -> fullMeshPtr(), 1.0, 21 );
    ElectrophysiologyUtility::setValueOnBoundary ( * (monodomain -> potentialPtr() ), monodomain -> fullMeshPtr(), 1.0, 22 );

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

    Real rho, poisson, young, bulk, alpha, gamma, mu;
    rho     = dataFile ( "solid/physics/density", 1. );
    young   = dataFile ( "solid/physics/young",   1. );
    poisson = dataFile ( "solid/physics/poisson", 1. );
    bulk    = dataFile ( "solid/physics/bulk",    1. );
    alpha   = dataFile ( "solid/physics/alpha",   1. );
    gamma   = dataFile ( "solid/physics/gamma",   1. );
    mu   = dataFile ( "solid/physics/mu",   1. );
    //  M_gammaf  = dataFile ( "solid/physics/gammaf",  0. );

    if ( comm->MyPID() == 0 )
    {
        std::cout << "density = " << rho     << std::endl
                  << "young   = " << young   << std::endl
                  << "poisson = " << poisson << std::endl
                  << "bulk    = " << bulk    << std::endl
                  << "alpha   = " << alpha   << std::endl
                  << "gamma   = " << gamma   << std::endl;
    }


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

    meshPtr_Type fullSolidMesh;
    meshPtr_Type localSolidMesh;

    std::string solidMeshName = parameterList.get ("solid_mesh_name", "no_solid_mesh");
    bool usingDifferentMeshes = true;
    if (solidMeshName == "no_solid_mesh" || solidMeshName == meshName )
    {
        usingDifferentMeshes = false;
    }
    else
    {
        if ( comm->MyPID() == 0 )
        {
            std::cout << "\nI am using different meshes: " << meshName << " for the electro and " << solidMeshName << " for the solid\n\n" << std::endl;
        }
    }

    if ( usingDifferentMeshes  )
    {
        fullSolidMesh.reset (new mesh_Type ( comm ) );
        localSolidMesh.reset (new mesh_Type ( comm ) );
        MeshUtility::fillWithFullMesh (localSolidMesh, fullSolidMesh,  solidMeshName,  parameterList.get ("solid_mesh_path", "") );
    }
    else
    {
        fullSolidMesh = monodomain -> fullMeshPtr();
        localSolidMesh = monodomain -> localMeshPtr();
    }

    std::string dOrder =  dataFile ( "solid/space_discretization/order", "P1");
    solidFESpacePtr_Type dFESpace ( new solidFESpace_Type (localSolidMesh, dOrder, 3, comm) );
    solidFESpacePtr_Type aFESpace ( new solidFESpace_Type (monodomain -> localMeshPtr(), dOrder, 1, comm) );
    solidETFESpacePtr_Type dETFESpace ( new solidETFESpace_Type (localSolidMesh, & (dFESpace->refFE() ), & (dFESpace->fe().geoMap() ), comm) );
    scalarETFESpacePtr_Type aETFESpace ( new scalarETFESpace_Type (monodomain -> localMeshPtr(), & (aFESpace->refFE() ), & (aFESpace->fe().geoMap() ), comm) );
    solidFESpacePtr_Type solidaFESpace ( new solidFESpace_Type (localSolidMesh, "P1", 1, comm) );


    if ( comm->MyPID() == 0 )
    {
        std::cout << "\nsetup boundary conditions" << std::endl;
    }

    //! #################################################################################
    //! BOUNDARY CONDITIONS
    //! #################################################################################
    vector <ID> compx (1), compy (1), compz (1), compxy (2), compxz (2), compyz (2);
    compx[0] = 0;
    compy[0] = 1, compz[0] = 2;
    compxy[0] = 0;
    compxy[1] = 1;
    compxz[0] = 0;
    compxz[1] = 2;
    compyz[0] = 1;
    compyz[1] = 2;

    //   boost::shared_ptr<BCHandler> BCh ( new BCHandler() );
    //
    //    BCFunctionBase zero (bcZero);
    //    BCFunctionBase load (Private::boundaryLoad);
    //
    //
    //    //! =================================================================================
    //    //! BC for quarter ring
    //    //! =================================================================================
    ////    BCh->addBC ("EdgesIn",      29,  Essential, Component, zero,    compz);
    ////    BCh->addBC ("EdgesIn",      31,  Essential, Component, zero,    compy);
    ////    BCh->addBC ("EdgesIn",      32,  Essential, Component, zero,    compx);
    //    //! =================================================================================
    //    //! BC for idealHeart
    //    //! =================================================================================
    //    BCh->addBC ("EdgesIn",      40,  Essential, Full, zero,    3);
    //! =================================================================================
    bcInterfacePtr_Type                     solidBC ( new bcInterface_Type() );
    solidBC->createHandler();
    solidBC->fillHandler ( data_file_name, "solid" );

    //M_FSIoperator->setSolidBC ( M_solidBC->handler() );

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

    //    function_Type fibersDirection = &fiberRotation;
    //  vectorPtr_Type fibersRotated( new vector_Type( dFESpace -> map() ) );
    // dFESpace -> interpolate ( static_cast< FESpace< RegionMesh<LinearTetra>, MapEpetra >::function_Type > ( fibersDirection ), *fibersRotated , 0);

    //===========================================================
    //===========================================================
    //              FIBERS
    //===========================================================
    //===========================================================

    vectorPtr_Type solidFibers ( new vector_Type ( dFESpace -> map() ) );

    if ( comm->MyPID() == 0 )
    {
        std::cout << "\nread fibers" << std::endl;
    }

    ElectrophysiologyUtility::importFibers (solidFibers, parameterList.get ("solid_fiber_file", ""), localSolidMesh );

    if ( comm->MyPID() == 0 )
    {
        std::cout << "\nset fibers" << std::endl;
    }

    solid.material() -> setFiberVector ( *solidFibers );

    //     monodomain -> setupFibers();

    vectorPtr_Type gammaf ( new vector_Type ( ( monodomain -> globalSolution().at (3) ) -> map() ) );
    vectorPtr_Type solidGammaf;
    vectorPtr_Type emDisp;
    solidFESpacePtr_Type electroFiberFESpace;
    solidETFESpacePtr_Type electrodETFESpace;
    if (usingDifferentMeshes)
    {



        electroFiberFESpace.reset ( new solidFESpace_Type (monodomain -> localMeshPtr(), "P1", 3, comm) );
        electrodETFESpace.reset ( new solidETFESpace_Type (monodomain -> localMeshPtr(), & (dFESpace->refFE() ), & (dFESpace->fe().geoMap() ), comm) );

        vectorPtr_Type electroFibers ( new vector_Type ( electroFiberFESpace -> map() ) );
        ElectrophysiologyUtility::importFibers (electroFibers, parameterList.get ("fiber_file", ""), monodomain -> localMeshPtr() );
        monodomain -> setFiberPtr ( electroFibers );
        emDisp.reset (  new vector_Type ( electroFibers -> map() ) );
        solidGammaf.reset ( new vector_Type ( solidaFESpace -> map() ) );

    }
    else
    {
        solidGammaf = gammaf;
        monodomain -> setFiberPtr ( solidFibers );
        emDisp = solid.displacementPtr();
        electroFiberFESpace = dFESpace;
        electrodETFESpace = dETFESpace;
    }


    monodomain -> exportFiberDirection();
    //********************************************//
    // Create the global matrix: mass + stiffness in ELECTROPHYSIOLOGY //
    //********************************************//
    if ( comm->MyPID() == 0 )
    {
        cout << "\nSetup operators:  dt = " << monodomain -> timeStep() << "\n" ;
    }

    monodomain -> setDisplacementPtr ( emDisp );
    monodomain -> setupLumpedMassMatrix();
    monodomain -> setupStiffnessMatrix();
    monodomain -> setupGlobalMatrix();

    if ( comm->MyPID() == 0 )
    {
        cout << "Done! \n" ;
    }

    //==================================================================//
    //==================================================================//
    //                 SETUP Activation                                //
    //==================================================================//
    //==================================================================//

    //   vectorPtr_Type gammaf( new vector_Type( monodomain -> globalSolution().at(3) -> map() ) );
    *gammaf *= 0;
    *solidGammaf *= 0;

    if ( comm->MyPID() == 0 )
    {
        std::cout << "\nset gammaf and fibers" << std::endl;
    }
    solid.material() -> setGammaf ( *solidGammaf );


    //==================================================================//
    //==================================================================//
    //                 SETUP INTERPOLATION                             //
    //==================================================================//
    //==================================================================//

    typedef RBFInterpolation<mesh_Type>           interpolation_Type;
    typedef boost::shared_ptr<interpolation_Type> interpolationPtr_Type;

    //Coarse To Fine ( C2F )
    interpolationPtr_Type C2F;
    //Fine To Coarse ( F2C )
    interpolationPtr_Type F2C;
    if (usingDifferentMeshes)
    {
        Teuchos::RCP< Teuchos::ParameterList > belosList = Teuchos::rcp ( new Teuchos::ParameterList );
        belosList = Teuchos::getParametersFromXmlFile ( "ParamList.xml" );

        if ( comm->MyPID() == 0 )
        {
            std::cout << "\nresetting" << std::endl;
        }
        std::string c2f = parameterList.get ("c2f", "RBFrescaledVectorial");
        std::string f2c = parameterList.get ("f2c", "RBFrescaledScalar");

        C2F.reset ( interpolation_Type::InterpolationFactory::instance().createObject ( c2f ) );
        F2C.reset ( interpolation_Type::InterpolationFactory::instance().createObject ( f2c ) );
        //       C2F.reset ( interpolation_Type::InterpolationFactory::instance().createObject ( "RBFlocallyRescaledVectorial" ) );
        //       F2C.reset ( interpolation_Type::InterpolationFactory::instance().createObject ( "RBFlocallyRescaledScalar" ) );

        int nFlags = 1;
        std::vector<int> flags (nFlags);
        flags[0] = -1;

        if ( comm->MyPID() == 0 )
        {
            std::cout << "\nintepolation setup" << std::endl;
        }

        C2F->setup ( fullSolidMesh , localSolidMesh, monodomain -> fullMeshPtr(), monodomain -> localMeshPtr(), flags);
        F2C->setup ( monodomain -> fullMeshPtr(), monodomain -> localMeshPtr(), fullSolidMesh , localSolidMesh, flags);


        if ( comm->MyPID() == 0 )
        {
            std::cout << "\nC2F: set Radius, ..." << std::endl;
        }

        C2F->setRadius ( (double) MeshUtility::MeshStatistics::computeSize (* (fullSolidMesh) ).maxH );

        if ( comm->MyPID() == 0 )
        {
            std::cout << "\nC2F: set data..." << std::endl;
        }
        C2F->setupRBFData ( solid.displacementPtr(), emDisp , dataFile, belosList);

        if ( comm->MyPID() == 0 )
        {
            std::cout << "\nC2F: Build operator..." << std::endl;
        }
        C2F->buildOperators();

        if ( comm->MyPID() == 0 )
        {
            std::cout << "\nF2C: set Radius, data, and build operator..." << std::endl;
        }

        F2C->setRadius ( (double) MeshUtility::MeshStatistics::computeSize (* (monodomain -> fullMeshPtr() ) ).maxH );
        F2C->setupRBFData ( gammaf, solidGammaf , dataFile, belosList);
        F2C->buildOperators();

        if ( comm->MyPID() == 0 )
        {
            std::cout << "\nC2F: interpolate..." << std::endl;
        }

        C2F->interpolate();
        C2F->solution (emDisp);

        if ( comm->MyPID() == 0 )
        {
            std::cout << "\nF2C: interpolate..." << std::endl;
        }
        F2C->interpolate();
        F2C->solution (solidGammaf);

    }

    //     function_Type initialGuess = &d0;
    //     vectorPtr_Type initd( new vector_Type( dFESpace -> map() ) );
    //     dFESpace -> interpolate ( static_cast< FESpace< RegionMesh<LinearTetra>, MapEpetra >::function_Type > ( initialGuess ), *initd , 0);





    //     if ( comm->MyPID() == 0 )
    //    {
    //        std::cout << "\nnorm inf gammaf: " << solid.material() -> gammaf() -> normInf() << std::endl;
    //        std::cout << "\nnorm inf fiber: " << solid.material() -> fiberVector() -> normInf() << std::endl;
    //    }
    if ( comm->MyPID() == 0 )
    {
        std::cout << "\nbuild solid system" << std::endl;
    }

    solid.buildSystem (1.0);
    vectorPtr_Type rhs (new vector_Type (solid.displacement(), Unique) );
    vectorPtr_Type disp (new vector_Type (solid.displacement(), Unique) );
    vectorPtr_Type initialDisplacement (new vector_Type (solid.displacement(), Unique) );
    solid.initialize ( initialDisplacement );


    MPI_Barrier (MPI_COMM_WORLD);

    if ( comm->MyPID() == 0 )
    {
        std::cout << "\nsetup solid exporter" << std::endl;
    }

    boost::shared_ptr< Exporter<RegionMesh<LinearTetra> > > exporter;
    exporter.reset ( new ExporterHDF5<RegionMesh<LinearTetra> > ( dataFile, "structure" ) );

    exporter->setPostDir ( "./" );
    exporter->setMeshProcId ( localSolidMesh, comm->MyPID() );

    vectorPtr_Type solidDisp ( new vector_Type (solid.displacement(), exporter->mapType() ) );
    exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, "displacement", dFESpace, solidDisp, UInt (0) );
    exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::ScalarField, "solid_gammaf", solidaFESpace, solidGammaf, UInt (0) );
    exporter->postProcess ( 0 );




    //================================================================//
    //================================================================//
    //                    SETUP COUPLING SOLVER                       //
    //                                                                //
    //================================================================//
    //================================================================//
    ExporterHDF5< RegionMesh <LinearTetra> > expGammaf;
    expGammaf.setMeshProcId (monodomain -> localMeshPtr(), comm->MyPID() );
    expGammaf.setPrefix ("gammaf");


    //      expGammaf.addVariable(ExporterData<mesh_Type>::ScalarField, "gammaf",
    //              monodomain -> feSpacePtr(), gammaf, UInt(0));
    //    expGammaf.postProcess(0.0);
    //    Real min =  0.2;
    //    Real max =  0.85;
    //
    //    Real beta = -0.3;
    //
    //    HeartUtility::rescaleVector(*gammaf, min, max, beta);


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
    //==================================================================//
    //                    SETUP LINEAR SOLVER                             //
    //                        ACTIVATION                                  //
    //==================================================================//
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










    boost::shared_ptr<FLRelationshipGamma> flg (new FLRelationshipGamma);
    boost::shared_ptr<FLRelationship> fl (new FLRelationship);

    boost::shared_ptr<HeavisideFct> H (new HeavisideFct);
#define deformationGradientTensor ( grad( electrodETFESpace, *emDisp, 0) + value( solid.material()-> identity() ) )
#define RIGHTCAUCHYGREEN transpose(deformationGradientTensor) * deformationGradientTensor
#define firstInvariantC trace( RIGHTCAUCHYGREEN )
#define fiber0       ( value( electrodETFESpace, *( monodomain -> fiberPtr() ) ) )
#define fiber        ( deformationGradientTensor * fiber0 )
#define I4f             dot( fiber, fiber)
#define Ca    ( value( aETFESpace, *( monodomain -> globalSolution().at(3)  ) ) )
#define Ca2         (  Ca * Ca )
#define dCa         ( Ca + value(-0.02155) )
#define Gammaf          ( value( aETFESpace, *gammaf ) )
#define Pa         ( value(-2.5) * eval(H, dCa ) * eval(H, dCa )/*eval( fl,  I4f) */ * eval( flg,  Gammaf) + value(1.0) / ( ( GammaPlusOne ) * ( GammaPlusOne ) * ( GammaPlusOne ) * ( GammaPlusOne ) ) )
    //      #define Pa            beta * eval( fl,  I4f)

#define GammaPlusOne ( Gammaf + value(1.0) )
    //#define dW1       ( ( I4f - value(1.0) ) )
#define dW1     ( ( Gammaf * Gammaf + 2.0 * Gammaf ) )
#define dW      ( ( ( dW1 ) + value(1.0) / ( ( GammaPlusOne ) * ( GammaPlusOne ) * ( GammaPlusOne ) ) ) )

    //   #define dgGammaf ( value(-1.0) + value(-2.0) / ( GammaPlusOne ) + value(2.0) * Gammaf * ( Gammaf + value(2.0)  )* pow( GammaPlusOne, -3 ) )
    //   #define activationEquation value(-1.0) * ( Pa  -  ( value(2.0) * GammaPlusOne * firstInvariantC + dgGammaf * I4f )  * value( mu / 2.0 ) ) / beta
#define activationEquation value(0.0005) * (Pa - dW) / ( Ca2 )

    vectorPtr_Type tmpRhsActivation ( new vector_Type ( rhsActivation -> map(), Repeated ) );

    expGammaf.addVariable (ExporterData<mesh_Type>::ScalarField, "gammaf",
                           monodomain -> feSpacePtr(), gammaf, UInt (0) );
    expGammaf.addVariable (ExporterData<mesh_Type>::VectorField, "interpolated displacement",
                           monodomain -> feSpacePtr(), emDisp, UInt (0) );
    expGammaf.addVariable (ExporterData<mesh_Type>::ScalarField, "rhs",
                           monodomain -> feSpacePtr(), rhsActivation, UInt (0) );


    expGammaf.postProcess (0.0);


    //===========================================================
    //===========================================================
    //              TIME LOOP
    //===========================================================
    //===========================================================
    Real emdt = parameterList.get ("emdt", 1.0);
    int iter ( (emdt / monodomain -> timeStep() ) );
    int k (0);
    Real saveStep = parameterList.get ("save_step", 1.0);
    int saveIter ( (saveStep / monodomain -> timeStep() ) );



    bool twoWayCoupling = parameterList.get ("two_way", false);
    for ( Real t (0.0); t < monodomain -> endTime(); )
    {
        t = t + monodomain -> timeStep();
        k++;

        monodomain -> solveOneSplittingStep();


        //        *gammaf = *( monodomain -> globalSolution().at(3) );
        //        Real min =  0.2;
        //        Real max =  0.85;
        //
        //        Real beta = -0.3;

        //        HeartUtility::rescaleVector(*gammaf, min, max, beta);

        //    if(finiteElement)
        //    {

        //    }
        //        else
        //        {
        //
        //            int size = gammaf -> epetraVector().MyLength();
        //            int i = 0;
        //            for( int j(0); j< size; j++)
        //            {
        //                i = gammaf -> blockMap().GID(j);
        //                gammaf(i) = gammaf(i) + monodomain -> timeStep() *
        //            }
        //
        //        }
        *tmpRhsActivation *= 0;
        {
            using namespace ExpressionAssembly;




            integrate ( elements ( monodomain -> localMeshPtr() ),
                        monodomain -> feSpacePtr() -> qr() ,
                        monodomain -> ETFESpacePtr(),
                        activationEquation  * phi_i
                      ) >> tmpRhsActivation;

        }
        *rhsActivation *= 0;
        *rhsActivation = ( * (mass) * ( *gammaf ) );
        *rhsActivation += ( ( monodomain -> timeStep() * *tmpRhsActivation ) );

        linearSolver.setRightHandSide (rhsActivation);
        linearSolver.solve (gammaf);

        if ( k % iter == 0)
        {




            if (usingDifferentMeshes)
            {
                F2C -> updateRhs ( gammaf );
                F2C -> interpolate();
                F2C -> solution ( solidGammaf );
            }

            solid.material() -> setGammaf ( *solidGammaf );
            solid.iterate ( solidBC -> handler() );

            //        timeAdvance->shiftRight ( solid.displacement() );

            *solidDisp = solid.displacement();


            if (usingDifferentMeshes)
            {
                C2F -> updateRhs ( solid.displacementPtr() );
                C2F -> interpolate();
                C2F -> solution ( emDisp );
            }



            if (twoWayCoupling)
            {
                monodomain -> setupStiffnessMatrix();
                monodomain -> setupGlobalMatrix();

            }
        }
        //*solidVel  = timeAdvance->firstDerivative();
        //*solidAcc  = timeAdvance->secondDerivative();
        if ( k % saveIter == 0)
        {

            monodomain -> exportSolution (exp, t);
            expGammaf.postProcess (t);

            exporter->postProcess ( t );
        }
    }
    exp.closeFile();
    expGammaf.closeFile();

    exporter -> closeFile();
    if ( comm->MyPID() == 0 )
    {
        std::cout << "Active strain example: Passed!" << std::endl;
    }

#undef Gammaf
#undef deformationGradientTensor
#undef RIGHTCAUCHYGREEN
#undef firstInvariantC
#undef fiber0
#undef fiber
#undef Pa
#undef I4f
#undef beta
#undef GammaPlusOne
#undef dgGammaf
#undef activationEquation
#undef Ca
#undef Ca2
#undef dCa
#undef dW1
#undef dW



#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return 0;
}
