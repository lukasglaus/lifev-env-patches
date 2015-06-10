//********************************************//
// Includes
//********************************************//

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



// Namespaces
using namespace LifeV;


//********************************************//
// Volume computation
//********************************************//

void createPositionVector (const RegionMesh<LinearTetra>& mesh,
                           VectorEpetra& positionVector)
{
    Int nLocalDof = positionVector.epetraVector().MyLength();
    Int nComponentLocalDof = nLocalDof / 3;
    for (int k (0); k < nComponentLocalDof; k++)
    {
        UInt iGID = positionVector.blockMap().GID (k);
        UInt jGID = positionVector.blockMap().GID (k + nComponentLocalDof);
        UInt kGID = positionVector.blockMap().GID (k + 2 * nComponentLocalDof);

        positionVector[iGID] = mesh.point (iGID).x();
        positionVector[jGID] = mesh.point (iGID).y();
        positionVector[kGID] = mesh.point (iGID).z();
    }
    
}

Real openEndVolume (const RegionMesh<LinearTetra>& mesh,
                           VectorEpetra& positionVector)
{

    int bcflag (36);
    int component(1);
    int direction(1);
    
    
    // Determine points at boundary
    std::set< unsigned int> vertexIds;
    for (UInt iBFaceIn = 0; iBFaceIn < mesh.numBFaces(); ++iBFaceIn)
    {
        UInt markerIdIn = mesh.boundaryFace(iBFaceIn).markerID();
        if ( markerIdIn == bcflag )
        {
            for (UInt iBFaceOut = 0; iBFaceOut < mesh.numBFaces(); ++iBFaceOut)
            {
                UInt markerIdOut = mesh.boundaryFace(iBFaceOut).markerID();
                if ( markerIdOut != bcflag )
                {
                    for (UInt iBPointIn = 0; iBPointIn < mesh.boundaryFace(iBFaceIn).S_numPoints; ++ iBPointIn)
                    {
                        for (UInt iBPointOut = 0; iBPointOut < mesh.boundaryFace(iBFaceOut).S_numPoints; ++ iBPointOut)
                        {
                            UInt pointIdIn = mesh.boundaryFace(iBFaceIn).point(iBPointIn).id();
                            UInt pointIdOut = mesh.boundaryFace(iBFaceOut).point(iBPointOut).id();

                            if ( pointIdIn == pointIdOut )
                            {
                                vertexIds.insert(pointIdIn);
                            }
                        }
                    }
                }
            }
        }
    }
    
    
    // Determine center of point cloud
    Vector3D center;
    for (auto it = vertexIds.begin(); it != vertexIds.end(); ++it)
    {
        center += mesh.point(*it).coordinates() / vertexIds.size();
    }
    
    
    // Order points by angles
    std::map<double, unsigned int> vertexIdsOrdered;
    unsigned int compOne = (1 + component) % 3;
    unsigned int compTwo = (2 + component) % 3;
    for (auto it = vertexIds.begin(); it != vertexIds.end(); ++it)
    {
        double angle = std::atan2( mesh.point(*it).coordinate(compOne) - center(compOne) , mesh.point(*it).coordinate(compTwo) - center(compTwo) );
        vertexIdsOrdered.insert( std::pair<double, unsigned int> (angle, *it) );
    }

    
    // Compute projected area
    double area (0.0);
    auto itNext = vertexIdsOrdered.begin();
    unsigned int i (0);
    for (auto it = vertexIdsOrdered.begin(); it != vertexIdsOrdered.end(); ++it)
    {
        if ( i++ < vertexIdsOrdered.size() - 1 ) std::advance(itNext, 1);
        else std::advance(itNext, - (vertexIdsOrdered.size() -1));

        area += 0.5 * ( mesh.point(itNext->second).coordinate(compTwo) + mesh.point(it->second).coordinate(compTwo) ) * ( mesh.point(itNext->second).coordinate(compOne) - mesh.point(it->second).coordinate(compOne) );
    }
    
    std::cout << "------------------------------------ nP: " << vertexIds.size() << std::endl << area << std::endl;
    
    return ( direction * center(component) * area );
    
}

template<typename space> Real ComputeVolume (const boost::shared_ptr<RegionMesh<LinearTetra> > localMesh,
                                             VectorEpetra positionVector,
                                             const VectorEpetra& disp,
                                             const boost::shared_ptr <ETFESpace<RegionMesh<LinearTetra>, MapEpetra, 3, 1> > ETFESpace,
                                             const boost::shared_ptr<space> dETFESpace,
                                             int bdFlag,
                                             boost::shared_ptr<Epetra_Comm> comm)
{
    
    Real fluidVolume;
    
    MatrixSmall<3, 3> Id;
    Id (0, 0) = 1., Id (0, 1) = 0., Id (0, 2) = 0.;
    Id (1, 0) = 0., Id (1, 1) = 1., Id (1, 2) = 0.;
    Id (2, 0) = 0., Id (2, 1) = 0., Id (2, 2) = 1.;
    VectorSmall<3> E1;
    E1 (0) = 1., E1 (1) = 0., E1 (2) = 0.; ////
    
    Int nLocalDof = positionVector.epetraVector().MyLength();
    for (int k (0); k < nLocalDof; k++)
    {
        UInt iGID = positionVector.blockMap().GID ( k ); ////
        
        positionVector[iGID] += disp[iGID];
    }
    
    boost::shared_ptr<VectorEpetra> intergral ( new VectorEpetra ( positionVector.map() ) );
    
    {
        using namespace ExpressionAssembly;
        
        BOOST_AUTO_TPL (I, value (Id) );
        BOOST_AUTO_TPL (vE1, value (E1) );
        BOOST_AUTO_TPL (Grad_u, grad (dETFESpace, disp, 0) ); ////
        BOOST_AUTO_TPL (x, value (ETFESpace, positionVector) );
        BOOST_AUTO_TPL (F, (Grad_u + I) );
        BOOST_AUTO_TPL (FmT, minusT (F) );
        BOOST_AUTO_TPL (J, det (F) );
        BOOST_AUTO_TPL (x1, dot (x, vE1) );
        
        QuadratureBoundary myBDQR (buildTetraBDQR (quadRuleTria4pt) );
        
        *intergral *= 0.0;
        integrate (boundary (localMesh, bdFlag), myBDQR, ETFESpace,
                   value (-1.0) * J * dot (vE1, FmT * Nface) * phi_i) >> intergral;
        
        intergral->globalAssemble();
        //        *position = *positionR;
        //        for
        //        fluidVolume = position ->
        
        fluidVolume = positionVector.dot (*intergral);
        
        if (comm->MyPID() == 0)
        {
            std::cout << "\nFluid volume: " << fluidVolume << " in processor "
            << comm->MyPID() << std::endl;
        }
        return fluidVolume;
    }
}


//********************************************//
// Functions
//********************************************//

Real Iapp (const Real& t, const Real&  X, const Real& Y, const Real& Z, const ID& /*i*/)
{
    return ( Y > 1.5 && Y < 3 /*std::abs(X) < 1 && std::abs(Z-3) < 1 && Y < 0*/ /*&& Y < 0.25 && Z < 0.25 */ && t < 2 ? 30 : 0 );
    // setAppliedCurrent in electrophys. module.
}

Real potentialMultiplyerFcn (const Real& t, const Real&  X, const Real& Y, const Real& Z, const ID& /*i*/)
{
    return 0; // ( Y < 2.5 && Y > 0.5 ? 1.0 : 0.0 );
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
    

    //********************************************//
    // Declare communicator and solver
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
    // Setup material data
    //********************************************//
    
    EMData emdata;
    emdata.setup (dataFile);
    
    
    //********************************************//
    // Load mesh
    //********************************************//
    
    if ( comm->MyPID() == 0 )
    {
        std::cout << "Load mesh...\n";
    }
    
    std::string meshName = dataFile("solid/space_discretization/mesh_file", "cube4.mesh");
    std::string meshPath = dataFile("solid/space_discretization/mesh_dir", "cube4.mesh");

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
    
    Real meshScaling = dataFile("solid/space_discretization/mesh_scaling", 1.0);
    std::vector<Real> scale (3, meshScaling);
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
    // Setup solver (including fe-spaces & b.c.)
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
 
    bool anisotropy = dataFile ( "solid/space_discretization/anisotropic", false );

    if ( anisotropy )
    {
        std::string fiberFileName  =  dataFile ( "solid/space_discretization/fiber_name", "FiberDirection");
        std::string sheetFileName  =  dataFile ( "solid/space_discretization/sheet_name", "SheetsDirection");
        std::string fiberFieldName =  dataFile ( "solid/space_discretization/fiber_fieldname", "fibers");
        std::string sheetFieldName =  dataFile ( "solid/space_discretization/sheet_fieldname", "sheets");
        std::string fiberDir       =  dataFile ( "solid/space_discretization/fiber_dir", "./");
        std::string sheetDir       =  dataFile ( "solid/space_discretization/sheet_dir", "./");
        
        solver.setupFiberVector ( fiberFileName, fiberFieldName, fiberDir );
        solver.setupMechanicalSheetVector ( sheetFileName, sheetFieldName, sheetDir );
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
    
    solver.twoWayCoupling();
    solver.structuralOperatorPtr()->setNewtonParameters(dataFile);
    solver.buildSystem();
    
    if( 0 == comm->MyPID() )
    {
    	std::cout << " done!" << std::endl;
    }

    
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
    // Setup vector & exporter for activation time
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
    // Electric stimulus function
    //********************************************//
    
    function_Type stim = &Iapp;
    
    
    //********************************************//
    // Body circulation
    //********************************************//
    
    Circulation circulationSolver("inputfile");

    std::vector<std::vector<std::string> > bcNames { { "lv" , "p" } , { "rv" , "p" } };
    std::vector<double> bcValues(2, 5);
    
    circulationSolver.exportSolution( "solution.txt" );
    circulationSolver.iterate(0.01, bcNames, bcValues, 0);
    circulationSolver.exportSolution( "solution.txt" );

    
    //********************************************//
    // Boundary conditions
    //********************************************//

    // Get b.c. flags
    ID LvFlag =  dataFile ( "solid/boundary_conditions/LvFlag", 0);
    
    // Get preload pressure and steps
    Real pPreloadLvBC = dataFile ( "solid/boundary_conditions/LvPreloadPressure", 0.0);
    int preloadSteps = dataFile ( "solid/boundary_conditions/numPreloadSteps", 0);
    
    // Declare pressure b.c. obtained from circulation
    Real pLvBC;
    
    // Boundary vector normal in deformed configuration
    solver.structuralOperatorPtr() -> setBCFlag( LvFlag );

    // Create b.c. for endocardium which can be modified during time loop
    vectorPtr_Type endoLvVectorPtr( new vector_Type ( solver.structuralOperatorPtr() -> displacement().map(), Repeated ) );
    bcVectorPtr_Type pLvBCVectorPtr( new bcVector_Type( *endoLvVectorPtr, solver.structuralOperatorPtr() -> dispFESpacePtr() -> dof().numTotalDof(), 1 ) );
    solver.bcInterfacePtr() -> handler() -> addBC("LvPressure", LvFlag, Natural, Full, *pLvBCVectorPtr, 3);
    solver.bcInterfacePtr() -> handler() -> bcUpdate( *solver.structuralOperatorPtr() -> dispFESpacePtr() -> mesh(), solver.structuralOperatorPtr() -> dispFESpacePtr() -> feBd(), solver.structuralOperatorPtr() -> dispFESpacePtr() -> dof() );

    
    //********************************************//
    // Volume computation
    //********************************************//
    
    vectorPtr_Type referencePosition( new vector_Type( solver.structuralOperatorPtr() -> displacement().map() ) );
    createPositionVector( *solver.fullMeshPtr(), *referencePosition );

    Real fluidVolume =  ComputeVolume( solver.localMeshPtr(), *referencePosition, solver.structuralOperatorPtr() -> displacement(), solver.electroSolverPtr() -> ETFESpacePtr(), solver.electroSolverPtr() -> displacementETFESpacePtr() , 36, comm );

    auto a = openEndVolume(*solver.fullMeshPtr(), *referencePosition);

    //********************************************//
    // Preload
    //********************************************//
    
    for (int i (1); i <= preloadSteps; i++)
    {
        std::cout << "\n*********************";
        std::cout << "\nPreload step: " << i << " / " << preloadSteps;
        std::cout << "\n*********************\n";
        
        // Update pressure b.c.
        *endoLvVectorPtr = - pPreloadLvBC * i / preloadSteps;
        pLvBCVectorPtr.reset ( ( new bcVector_Type (*endoLvVectorPtr, solver.structuralOperatorPtr() -> dispFESpacePtr() -> dof().numTotalDof(), 1) ) );
        solver.bcInterfacePtr() -> handler() -> modifyBC(LvFlag, *pLvBCVectorPtr);
        
        // Solve mechanics
        solver.bcInterfacePtr() -> updatePhysicalSolverVariables();
        solver.solveMechanics();
    }
    

    //********************************************//
    // Time loop
    //********************************************//
    
    Real dt_activation = solver.data().electroParameter<Real>("timestep");
    Real dt_mechanics = solver.data().solidParameter<Real>("timestep");
    Real endtime = solver.data().electroParameter<Real>("endtime");
    UInt saveIter = static_cast<UInt>( dt_mechanics / dt_activation );
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
        solver.solveActivation (dt_activation);
        
        if ( k % saveIter == 0 )
        {
            // Update pressure b.c.
            pLvBC = pPreloadLvBC; //Circulation::computePressure(1/pPreloadLvBC);
            *endoLvVectorPtr = - pLvBC;
            pLvBCVectorPtr.reset ( ( new bcVector_Type (*endoLvVectorPtr, solver.structuralOperatorPtr() -> dispFESpacePtr() -> dof().numTotalDof(), 1) ) );
            solver.bcInterfacePtr() -> handler() -> modifyBC(LvFlag, *pLvBCVectorPtr);
            
            // Solve mechanics
            solver.structuralOperatorPtr() -> data() -> dataTime() -> setTime(t);
            solver.bcInterfacePtr() -> updatePhysicalSolverVariables();
            solver.solveMechanics();
        }
        
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
