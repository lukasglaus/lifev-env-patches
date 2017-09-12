//============================================//
// Includes
//============================================//

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

// Volume computation
#include <lifev/em/solver/circulation/CirculationVolumeIntegrator.hpp>

// Heart solver
#include <lifev/em/solver/HeartSolver.hpp>

// Track nan
// #include <fenv.h>





// Namespaces
using namespace LifeV;


//============================================//
// Functions
//============================================//

const VectorEpetra
undeformedPositionVector (const boost::shared_ptr<RegionMesh<LinearTetra> > fullMeshPtr, const boost::shared_ptr<FESpace<RegionMesh<LinearTetra>, MapEpetra >> dFeSpace)
{
    // New P1 Space
    FESpace<RegionMesh<LinearTetra> , MapEpetra > p1FESpace ( dFeSpace->mesh(), "P1", 3, dFeSpace->map().commPtr() );
    
    // Create P1 VectorEpetra
    VectorEpetra p1PositionVector (p1FESpace.map());
    
    // Fill P1 vector with mesh values
    Int p1nCompLocalDof = p1PositionVector.epetraVector().MyLength() / 3;
    for (int j (0); j < p1nCompLocalDof; j++)
    {
        UInt iGID = p1PositionVector.blockMap().GID (j);
        UInt jGID = p1PositionVector.blockMap().GID (j + p1nCompLocalDof);
        UInt kGID = p1PositionVector.blockMap().GID (j + 2 * p1nCompLocalDof);
        
        p1PositionVector[iGID] = fullMeshPtr->point (iGID).x();
        p1PositionVector[jGID] = fullMeshPtr->point (iGID).y();
        p1PositionVector[kGID] = fullMeshPtr->point (iGID).z();
    }
    
    // Interpolate position vector from P1-space to current space
    VectorEpetra positionVector ( dFeSpace->map() );
    positionVector = dFeSpace -> feToFEInterpolate(p1FESpace, p1PositionVector);
    
    return positionVector;
}


Real patchForce (const Real& t, const Real& Tmax, const Real& tmax, const Real& tduration)
{
    bool time ( fmod(t-tmax+0.5*tduration, 800.) < tduration && fmod(t-tmax+0.5*tduration, 800.) > 0);
    Real force = std::pow( std::sin(fmod(t-tmax+0.5*tduration, 800.)*3.14159265359/tduration) , 2 ) * Tmax;
    return ( time ? force : 0 );
}

Real patchFunction (const Real& t, const Real&  X, const Real& Y, const Real& Z, const ID& /*i*/)
{
    Real disp = std::pow( std::sin(fmod(t, 800.) * 3.14159265359/300) , 2 )*15;
    return disp;
}

Real Iapp (const Real& t, const Real&  X, const Real& Y, const Real& Z, const ID& /*i*/)
{
    bool coords ( Y < -7. );
    //bool coords ( Y > 4. ); //( Y > 1.5 && Y < 3 );
    bool time ( fmod(t, 800.) < 4 && fmod(t, 800.) > 2);
    return ( coords && time ? 30 : 0 );
}

Real potentialMultiplyerFcn (const Real& t, const Real&  X, const Real& Y, const Real& Z, const ID& /*i*/)
{
    bool time ( fmod(t, 800.) < 4 && fmod(t, 800.) > 2);
    return 1.4 * time; // ( Y < 2.5 && Y > 0.5 ? 1.0 : 0.0 );
}


//class Stimulus
//{
//public:
//    Stimulus(Real ta, Real tb, Real y, Real f) :
//        M_ta    (ta),
//        M_tb    (tb),
//        M_y     (y),
//        M_f     (f)
//    {}
//    
//    Real Iapp (const Real& t, const Real&  X, const Real& Y, const Real& Z, const ID& i)
//    {
//        bool coords ( Y < M_y );
//        //bool coords ( Y > 4. ); //( Y > 1.5 && Y < 3 );
//        bool time ( fmod(t, 800.) < M_tb && fmod(t, 800.) > M_ta);
//        return ( coords && time ? M_f : 0 );
//    }
//    
//private:
//    Real M_ta;
//    Real M_tb;
//    Real M_y;
//    Real M_f;
//    
//};


int main (int argc, char** argv)
{

//    feenableexcept(FE_INVALID | FE_OVERFLOW);
    
    
    //============================================//
    // Typedefs
    //============================================//
    
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
    

    //============================================//
    // Communicator and displayer
    //============================================//
#ifdef HAVE_MPI
    MPI_Init ( &argc, &argv );
#endif

    boost::shared_ptr<Epetra_Comm>  comm ( new Epetra_MpiComm (MPI_COMM_WORLD) );
    
    Displayer displayer ( comm );
    displayer.leaderPrint ("\nUsing MPI\n");

    
    //============================================//
    // Read data file and create output folder
    //============================================//
    GetPot command_line (argc, argv);
    const std::string data_file_name = command_line.follow ("data", 2, "-f", "--file");
    GetPot dataFile (data_file_name);
    std::string problemFolder = EMUtility::createOutputFolder (command_line, *comm);

    
    //============================================//
    // Electromechanic solver
    //============================================//
    EMSolver<mesh_Type, monodomain_Type> solver(comm);
    
    
    //============================================//
    // Body circulation
    //============================================//
    const std::string circulationInputFile = command_line.follow ("circulation", 2, "-cif", "--cifile");
    const std::string circulationOutputFile = command_line.follow ( (problemFolder + "solution.dat").c_str(), 2, "-cof", "--cofile");
    
    Circulation circulationSolver( circulationInputFile );
    
    // Flow rate between two vertices
    auto Q = [&circulationSolver] (const std::string& N1, const std::string& N2) { return circulationSolver.solution ( std::vector<std::string> {N1, N2} ); };
    auto p = [&circulationSolver] (const std::string& N1) { return circulationSolver.solution ( N1 ); };
    
    
    //============================================//
    // Heart solver
    //============================================//
    HeartSolver<EMSolver<mesh_Type, monodomain_Type> > heartSolver (solver, circulationSolver);
    
    
    //============================================//
    // Setup material data
    //============================================//
    EMData emdata;
    emdata.setup (dataFile);
    
    
    //============================================//
    // Load mesh
    //============================================//
    displayer.leaderPrint ("\nLoading mesh ... ");

    std::string meshName = dataFile("solid/space_discretization/mesh_file", "cube4.mesh");
    std::string meshPath = dataFile("solid/space_discretization/mesh_dir", "./");
    
    solver.loadMesh (meshName, meshPath);
    
    displayer.leaderPrint ("\ndone!");

    
    //============================================//
    // Resize mesh
    //============================================//
    displayer.leaderPrint ("\nResizing mesh ... ");

    std::vector<Real> scale (3, dataFile("solid/space_discretization/mesh_scaling", 1.0));
    std::vector<Real> rotate { dataFile("solid/space_discretization/mesh_rotation_0", 0.0) , dataFile("solid/space_discretization/mesh_rotation_1", 0.0) , dataFile("solid/space_discretization/mesh_rotation_2", 0.0) };
    std::vector<Real> translate { dataFile("solid/space_discretization/mesh_translation_0", 0.0) , dataFile("solid/space_discretization/mesh_translation_1", 0.0) , dataFile("solid/space_discretization/mesh_translation_2", 0.0) };
    
    MeshUtility::MeshTransformer<mesh_Type> transformerFull (* (solver.fullMeshPtr() ) );
    MeshUtility::MeshTransformer<mesh_Type> transformerLocal (* (solver.localMeshPtr() ) );
    
    transformerFull.transformMesh (scale, rotate, translate);
    transformerLocal.transformMesh (scale, rotate, translate);
    
    displayer.leaderPrint ("\ndone!");
    
    
    //============================================//
    // Setup solver (including fe-spaces & b.c.)
    //============================================//
    displayer.leaderPrint ("\nSetting up EM solver ... ");
    
    EMAssembler::quadRule.setQuadRule( dataFile ( "solid/space_discretization/quad_rule", "4pt") );
    solver.setup (dataFile);
    
    displayer.leaderPrint ("\ndone!");

    
    //============================================//
    // Setup anisotropy vectors
    //============================================//
    displayer.leaderPrint ("\nSetting up anisotropy vectors ... ");

    bool anisotropy = dataFile ( "solid/space_discretization/anisotropic", false );

    if ( anisotropy )
    {
        std::string fiberFileName  =  dataFile ( "solid/space_discretization/fiber_name", "FiberDirection");
        std::string sheetFileName  =  dataFile ( "solid/space_discretization/sheet_name", "SheetsDirection");
        std::string fiberFieldName =  dataFile ( "solid/space_discretization/fiber_fieldname", "fibers");
        std::string sheetFieldName =  dataFile ( "solid/space_discretization/sheet_fieldname", "sheets");
        std::string fiberDir       =  meshPath; //dataFile ( "solid/space_discretization/fiber_dir", "./");
        std::string sheetDir       =  meshPath; //dataFile ( "solid/space_discretization/sheet_dir", "./");
        std::string elementOrder   =  dataFile ( "solid/space_discretization/order", "P1");

        solver.setupFiberVector ( fiberFileName, fiberFieldName, fiberDir, elementOrder );
        solver.setupMechanicalSheetVector ( sheetFileName, sheetFieldName, sheetDir, elementOrder );
    }
    else
    {
        solver.setupFiberVector (1., 0., 0.);
        solver.setupSheetVector (0., 1., 0.);
    }
    
    displayer.leaderPrint ("\ndone!");

    
    //============================================//
    // Initialize electrophysiology
    //============================================//
    displayer.leaderPrint ("\nInitialize electrophysiology ... ");

    solver.initialize();
    
    displayer.leaderPrint ("\ndone!");

    
    //============================================//
    // Building Matrices
    //============================================//
    displayer.leaderPrint ("\nBuilding matrices ... ");

    solver.oneWayCoupling();
    solver.structuralOperatorPtr()->setNewtonParameters(dataFile);
    solver.buildSystem();
    
    displayer.leaderPrint ("\ndone!");

    
    //============================================//
    // Setup exporters for EMSolver
    //============================================//
    displayer.leaderPrint ("\nSetting up exporters ... ");

    solver.setupExporters (problemFolder);
    
    displayer.leaderPrint ("\ndone!");
    
    
    //============================================//
    // Electric stimulus function
    //============================================//
    function_Type stim = &Iapp;
    
    
    //============================================//
    // Kept-normal boundary conditions
    //============================================//
    // Get b.c. flags
    // ID LvFlag =  dataFile ( "solid/boundary_conditions/LvFlag", 0);
    
    // Boundary vector normal in deformed configuration
    // solver.structuralOperatorPtr() -> setBCFlag( LvFlag );
    
    //solver.bcInterfacePtr() -> handler() -> addBC("LvPressure", LVFlag, Natural, Full, *pLvBCVectorPtr, 3); // BC for using function which keeps bc normal
    // Todo: Normal boundary condition!!

    
    //============================================//
    // B.C. endocardia and patches
    //============================================//
    std::vector<vectorPtr_Type> pVecPtrs;
    std::vector<bcVectorPtr_Type> pBCVecPtrs;
    std::vector<vectorPtr_Type> pVecPatchesPtrs;
    std::vector<bcVectorPtr_Type> pBCVecPatchesPtrs;

    std::vector<ID> flagsBC;
    std::vector<ID> flagsBCPatches;

    std::vector<UInt> ventIdx;

    // Endocardia
    UInt nVarBC = dataFile.vector_variable_size ( ( "solid/boundary_conditions/listVariableBC" ) );
    for ( UInt i (0) ; i < nVarBC ; ++i )
    {
        std::string varBCSection = dataFile ( ( "solid/boundary_conditions/listVariableBC" ), " ", i );
        flagsBC.push_back( dataFile ( ("solid/boundary_conditions/" + varBCSection + "/flag").c_str(), 0 ) );
        ventIdx.push_back ( dataFile ( ("solid/boundary_conditions/" + varBCSection + "/index").c_str(), 0 ) );
        
        pVecPtrs.push_back ( vectorPtr_Type ( new vector_Type ( solver.structuralOperatorPtr() -> displacement().map(), Repeated ) ) );
        pBCVecPtrs.push_back ( bcVectorPtr_Type( new bcVector_Type( *pVecPtrs[i], solver.structuralOperatorPtr() -> dispFESpacePtr() -> dof().numTotalDof(), 1 ) ) );
        solver.bcInterfacePtr() -> handler() -> addBC(varBCSection, flagsBC[i], Natural, Full, *pBCVecPtrs[i], 3);
    }
    
    // Force vector patches
    
    auto nGlobalPoints = solver.fullMeshPtr()->numGlobalPoints();
    auto nPoints = solver.fullMeshPtr()->numPoints();
    auto nGlobalFaces = solver.fullMeshPtr()->numGlobalFaces();
    auto nFaces = solver.fullMeshPtr()->numFaces();
    
    if ( 0 == comm->MyPID() )
    {
        std::cout << "\n*****************************************************************";
        std::cout << nGlobalPoints << "\n";
        std::cout << nGlobalFaces << "\n";
        std::cout << nPoints << "\n";
        std::cout << nFaces << "\n";
        std::cout << "\n*****************************************************************\n";
    }
    
    for (int j(0); j < nGlobalFaces; j++)
    {
//        UInt iGID = p1PositionVector.blockMap().GID (j);
//        UInt jGID = p1PositionVector.blockMap().GID (j + p1nCompLocalDof);
//        UInt kGID = p1PositionVector.blockMap().GID (j + 2 * p1nCompLocalDof);
        
        auto face = solver.fullMeshPtr()->face (j); //.x();
        
        for (int k(0); k < 3; ++k)
        {
            auto face_point = solver.fullMeshPtr()->face (j).point(j); //.x();

            auto iGID = j;
            
//            auto X = solver.fullMeshPtr()->point (iGID).x();
//            auto Y = solver.fullMeshPtr()->point (iGID).y();
//            auto Z = solver.fullMeshPtr()->point (iGID).z();
        }
        
    }

    
    UInt nVarPatchesBC = dataFile.vector_variable_size ( ( "solid/boundary_conditions/listForcePatchesBC" ) );
    for ( UInt i (0) ; i < nVarPatchesBC ; ++i )
    {
        std::string varBCPatchesSection = dataFile ( ( "solid/boundary_conditions/listPatchesBC" ), " ", i );
        flagsBCPatches.push_back ( dataFile ( ("solid/boundary_conditions/" + varBCPatchesSection + "/flag").c_str(), 0 ) );
        
        pVecPatchesPtrs.push_back ( vectorPtr_Type ( new vector_Type ( solver.structuralOperatorPtr() -> displacement().map(), Repeated ) ) );
        pBCVecPatchesPtrs.push_back ( bcVectorPtr_Type( new bcVector_Type( *pVecPatchesPtrs[i], solver.structuralOperatorPtr() -> dispFESpacePtr() -> dof().numTotalDof(), 1 ) ) );
        solver.bcInterfacePtr() -> handler() -> addBC(varBCPatchesSection, flagsBCPatches[i], Natural, Full, *pBCVecPatchesPtrs[i], 3);
    }

    // Displacement function patches (not working properly yet)
    BCFunctionBase bcFunction;
    bcFunction.setFunction(patchFunction);
    UInt nVarDispPatchesBC = dataFile.vector_variable_size ( ( "solid/boundary_conditions/listDispPatchesBC" ) );
    for ( UInt i (0) ; i < nVarDispPatchesBC ; ++i )
    {
        std::string varBCPatchesSection = dataFile ( ( "solid/boundary_conditions/listDispPatchesBC" ), " ", i );
        ID flag = dataFile ( ("solid/boundary_conditions/" + varBCPatchesSection + "/flag").c_str(), 0 );
        
        solver.bcInterfacePtr() -> handler() -> addBC(varBCPatchesSection, flag, Essential, Normal, bcFunction);
    }

    
    solver.bcInterfacePtr() -> handler() -> bcUpdate( *solver.structuralOperatorPtr() -> dispFESpacePtr() -> mesh(), solver.structuralOperatorPtr() -> dispFESpacePtr() -> feBd(), solver.structuralOperatorPtr() -> dispFESpacePtr() -> dof() );
    
    // Functions to modify b.c.
    auto modifyFeBC = [&] (const std::vector<Real>& bcValues)
    {
        for ( UInt i (0) ; i < nVarBC ; ++i )
        {
            *pVecPtrs[i] = - bcValues[ ventIdx[i] ] * 1333.224;
            // Check coordinates of pVecPtrs and assign only values to certain cells
            // Implement vector for both natural and essential b.c.
            pBCVecPtrs[i].reset ( ( new bcVector_Type (*pVecPtrs[i], solver.structuralOperatorPtr() -> dispFESpacePtr() -> dof().numTotalDof(), 1) ) );
            solver.bcInterfacePtr() -> handler() -> modifyBC(flagsBC[i], *pBCVecPtrs[i]);
        }
    };

    Real Tmax = dataFile ( "solid/patches/Tmax", 0. );
    Real tmax = dataFile ( "solid/patches/tmax", 0. );
    Real tduration = dataFile ( "solid/patches/tduration", 0. );

    auto modifyFeBCPatches = [&] (const Real& time)
    {
        auto undefPosVec = undeformedPositionVector(solver.fullMeshPtr(), solver.structuralOperatorPtr() -> dispFESpacePtr());
        
        for ( UInt i (0) ; i < nVarPatchesBC ; ++i )
        {
            if ( 0 == comm->MyPID() ) std::cout << "\nPatch force: " << patchForce(time, Tmax, tmax, tduration) << std::endl;
            *pVecPatchesPtrs[i] = - patchForce(time, Tmax, tmax, tduration) * 1333.224;
            
            // Set pVecPatchesPtrs to zero outside of patch area
            
            VectorSmall<3> X;
            Int nLocalDof = undefPosVec.epetraVector().MyLength() / 3;

            for (int ik (0); ik < nLocalDof; ik++)
            {
                UInt iGID = undefPosVec.blockMap().GID (ik);
                UInt jGID = undefPosVec.blockMap().GID (ik + nLocalDof);
                UInt kGID = undefPosVec.blockMap().GID (ik + 2 * nLocalDof);
                
                X[0] = undefPosVec[iGID];
                X[1] = undefPosVec[jGID];
                X[2] = undefPosVec[kGID];
                
                Vector3D center1, center2;
                Real radius1 = 2;
                Real radius2 = 2;
                center1[0] = -0.7;
                center1[1] = -4.7;
                center1[2] = -6;
                center2[0] = 3.8;
                center2[1] = 1.9;
                center2[2] = -6;
                
                bool patch1Area = (X - center1).norm() < radius1;
                bool patch2Area = (X - center2).norm() < radius2;

                if ( (! patch1Area) && (! patch2Area) )
                {
                    (*pVecPatchesPtrs[i]) [iGID] = 0;
                    (*pVecPatchesPtrs[i]) [jGID] = 0;
                    (*pVecPatchesPtrs[i]) [kGID] = 0;
                }
            }
            
            pBCVecPatchesPtrs[i].reset ( ( new bcVector_Type (*pVecPatchesPtrs[i], solver.structuralOperatorPtr() -> dispFESpacePtr() -> dof().numTotalDof(), 1) ) );
            solver.bcInterfacePtr() -> handler() -> modifyBC(flagsBCPatches[i], *pBCVecPatchesPtrs[i]);
        }
    };

    if ( 0 == comm->MyPID() ) solver.bcInterfacePtr() -> handler() -> showMe();
    
    
    //============================================//
    // Volume integrators
    //============================================//
    auto& disp = solver.structuralOperatorPtr() -> displacement();
    auto FESpace = solver.structuralOperatorPtr() -> dispFESpacePtr();
    auto dETFESpace = solver.electroSolverPtr() -> displacementETFESpacePtr();
    auto ETFESpace = solver.electroSolverPtr() -> ETFESpacePtr();
    
    std::vector<int> LVFlags;
    std::vector<int> RVFlags;

    for ( UInt i (0) ; i < nVarBC ; ++i )
    {
        std::string varBCSection = dataFile ( ( "solid/boundary_conditions/listVariableBC" ), " ", i );
        ID flag  =  dataFile ( ("solid/boundary_conditions/" + varBCSection + "/flag").c_str(), 0 );
        ID index =  dataFile ( ("solid/boundary_conditions/" + varBCSection + "/index").c_str(), 0 );
        switch ( index )
        {
            case 0: LVFlags.push_back ( flag );
                break;
            case 1: RVFlags.push_back ( flag );
                break;
            default:
                break;
        }
    }
    
    VolumeIntegrator LV (LVFlags, "Left Ventricle", solver.fullMeshPtr(), solver.localMeshPtr(), ETFESpace, FESpace);
    VolumeIntegrator RV (RVFlags, "Right Ventricle", solver.fullMeshPtr(), solver.localMeshPtr(), ETFESpace, FESpace);

    
    //============================================//
    // Set variables and functions
    //============================================//
    const Real dt_activation = solver.data().electroParameter<Real>("timestep");
    const Real dt_loadstep =  dataFile ( "solid/time_discretization/dt_loadstep", 1.0 );
    const Real activationLimit_loadstep =  dataFile ( "solid/time_discretization/activation_limit_loadstep", 0.0 );
    const Real dt_mechanics = solver.data().solidParameter<Real>("timestep");
    const Real dt_save = dataFile ( "exporter/save", 10. );
    const Real endtime = solver.data().electroParameter<Real>("endtime");
    const UInt mechanicsLoadstepIter = static_cast<UInt>( dt_loadstep / dt_activation );
    const UInt mechanicsCouplingIter = static_cast<UInt>( dt_mechanics / dt_activation );
    const UInt maxiter = static_cast<UInt>( endtime / dt_activation ) ;
    
    const Real pPerturbationFe = dataFile ( "solid/coupling/pPerturbationFe", 1e-2 );
    const Real pPerturbationCirc = dataFile ( "solid/coupling/pPerturbationCirc", 1e-3 );
    const Real couplingError = dataFile ( "solid/coupling/couplingError", 1e-6 );
    const UInt couplingJFeSubIter = dataFile ( "solid/coupling/couplingJFeSubIter", 1 );
    const UInt couplingJFeSubStart = dataFile ( "solid/coupling/couplingJFeSubStart", 1 );
    const UInt couplingJFeIter = dataFile ( "solid/coupling/couplingJFeIter", 1 );
    
    const Real dpMax = dataFile ( "solid/coupling/dpMax", 0.1 );

    std::vector<std::vector<std::string> > bcNames { { "lv" , "p" } , { "rv" , "p" } };
    std::vector<double> bcValues { p ( "lv" ) , p ( "rv") };
    std::vector<double> bcValuesPre ( bcValues );

    VectorSmall<4> ABdplv, ABdprv, ABcoef;
    ABcoef (0) = 55/24; ABcoef (1) = -59/24; ABcoef (2) = 37/24; ABcoef (3) = -3/8;

    VectorSmall<2> VCirc, VCircNew, VCircPert, VFe, VFeNew, VFePert, R, dp;
    MatrixSmall<2,2> JFe, JCirc, JR;

    UInt iter (0);
    Real t (0);
    
    auto printCoupling = [&] ( std::string label ) { if ( 0 == comm->MyPID() )
    {
        std::cout << "\n===============================================================";
        std::cout << "\nCoupling: " << label;
        std::cout << "\nNewton iteration nr. " << iter << " at time " << t;
        std::cout << "\nLV - Pressure: \t\t\t" << bcValues[0];
        std::cout << "\nLV - FE-Volume: \t\t" << VFeNew[0];
        std::cout << "\nLV - Circulation-Volume: \t" << VCircNew[0];
        std::cout << "\nLV - Residual: \t\t\t" << std::abs(VFeNew[0] - VCircNew[0]);
        std::cout << "\nRV - Pressure: \t\t\t" << bcValues[1];
        std::cout << "\nRV - FE-Volume : \t\t" << VFeNew[1];
        std::cout << "\nRV - Circulation-Volume: \t" << VCircNew[1];
        std::cout << "\nRV - Residual: \t\t\t" << std::abs(VFeNew[1] - VCircNew[1]);
        //std::cout << "\nJFe   = " << JFe;
        //std::cout << "\nJCirc = " << JCirc;
        //std::cout << "\nJR    = " << JR;
        std::cout << "\n==============================================================="; }
    };
    
    
    //============================================//
    // Load restart file
    //============================================//
    
    std::string restartInput = command_line.follow ("noRestart", 2, "-r", "--restart");
    const bool restart ( restartInput != "noRestart" );

    if ( restart )
    {
        const std::string restartDir = command_line.follow (problemFolder.c_str(), 2, "-rd", "--restartDir");
        
        Real dtExport = 10.;
        
        // Set time variable
        const unsigned int restartInputStr = std::stoi(restartInput);
        const unsigned int nIter = (restartInputStr - 1) * dtExport / dt_mechanics;
        t = nIter * dt_mechanics;

        // Set time exporter time index
        solver.setTimeIndex(restartInputStr + 1);

        // Load restart solutions from output files
        std::string polynomialDegree = dataFile ( "solid/space_discretization/order", "P2");

        ElectrophysiologyUtility::importVectorField ( solver.structuralOperatorPtr() -> displacementPtr(), "MechanicalSolution" , "displacement", solver.localMeshPtr(), restartDir, polynomialDegree, restartInput );

        for ( unsigned int i = 0; i < solver.electroSolverPtr()->globalSolution().size() ; ++i )
        {
            ElectrophysiologyUtility::importScalarField (solver.electroSolverPtr()->globalSolution().at(i), "ElectroSolution" , ("Variable" + std::to_string(i)), solver.localMeshPtr(), restartDir, polynomialDegree, restartInput );
        }

        ElectrophysiologyUtility::importScalarField (solver.activationModelPtr() -> fiberActivationPtr(), "ActivationSolution" , "Activation", solver.localMeshPtr(), restartDir, polynomialDegree, restartInput );
        //ElectrophysiologyUtility::importScalarField (solver.activationTimePtr(), "ActivationTimeSolution" , "Activation Time", solver.localMeshPtr(), restartDir, polynomialDegree, restartInput );

        circulationSolver.restartFromFile ( restartDir + "solution.dat" , nIter );
        
        // Set boundary mechanics conditions
        bcValues = { p ( "lv" ) , p ( "rv" ) };
        bcValuesPre = { p ( "lv" ) , p ( "rv" ) };
        modifyFeBC(bcValues);
    }

    
    //============================================//
    // Preload
    //============================================//
    
    if ( ! restart )
    {
        solver.structuralOperatorPtr() -> data() -> dataTime() -> setTime(0.0);
        
        const int preloadSteps = dataFile ( "solid/boundary_conditions/numPreloadSteps", 0);
        
        auto preloadPressure = [] (std::vector<double> p, const int& step, const int& steps)
        {
            for (auto& i : p) {i *= double(step) / double(steps);}
            return p;
        };
        
        LifeChrono chronoSave;
        chronoSave.start();

        solver.saveSolution (-1.0);
        
        if ( 0 == comm->MyPID() )
        {
            std::cout << "\n*****************************************************************";
            std::cout << "\nData stored in " << chronoSave.diff() << " s";
            std::cout << "\n*****************************************************************\n";
        }
        
        LifeChrono chronoPreload;
        chronoPreload.start();
        
        for (int i (1); i <= preloadSteps; i++)
        {
            if ( 0 == comm->MyPID() )
            {
                std::cout << "\n*****************************************************************";
                std::cout << "\nPreload step: " << i << " / " << preloadSteps;
                std::cout << "\n*****************************************************************\n";
            }
            
            // Update pressure b.c.
            modifyFeBC(preloadPressure(bcValues, i, preloadSteps));

            // Solve mechanics
            solver.bcInterfacePtr() -> updatePhysicalSolverVariables();
            solver.solveMechanics();
            //solver.saveSolution (i-1);
        }

        auto maxI4fValue ( solver.activationModelPtr()->I4f().maxValue() );
        auto minI4fValue ( solver.activationModelPtr()->I4f().minValue() );
        
        if ( 0 == comm->MyPID() )
        {
            std::cout << "\nI4fmax = " << maxI4fValue;
            std::cout << "\nI4fmin = " << minI4fValue << std::endl;
            
        }

        if ( 0 == comm->MyPID() )
        {
            std::cout << "\n*****************************************************************";
            std::cout << "\nPreload done in: " << chronoPreload.diff();
            std::cout << "\n*****************************************************************\n";
        }

    }

    
    //============================================//
    // Time loop
    //============================================//
    
    VFe[0] = LV.volume(disp, dETFESpace, - 1);
    VFe[1] = RV.volume(disp, dETFESpace, 1);
    VCirc = VFe;
    
    VectorEpetra dispPre ( disp );
    VectorEpetra dispCurrent ( disp );
    ID bdPowerFlag  =  dataFile ( ("solid/boundary_conditions/LVEndo/flag") , 0 );
    
    printCoupling("Initial values");
    
    auto perturbedPressure = [] (std::vector<double> p, const double& dp)
    {
        for (auto& i : p) {i += dp;}
        return p;
    };
    
    auto perturbedPressureComp = [] (std::vector<double> p, const double& dp, int comp)
    {
        p[comp] += dp;
        return p;
    };

    if ( ! restart )
    {
        solver.saveSolution (t);
        circulationSolver.exportSolution( circulationOutputFile );
    }

    for (int k (1); k <= maxiter; k++)
    {
        if ( 0 == comm->MyPID() )
        {
            std::cout << "\n*****************************************************************";
            std::cout << "\nTIME = " << t+dt_activation;
            std::cout << "\n*****************************************************************\n";
        }

        t = t + dt_activation;

        //============================================//
        // Solve electrophysiology and activation
        //============================================//

        auto maxI4fValue ( solver.activationModelPtr()->I4f().maxValue() );
        auto minI4fValue ( solver.activationModelPtr()->I4f().minValue() );
        
        if ( 0 == comm->MyPID() )
        {
            std::cout << "\nI4fmax = " << maxI4fValue;
            std::cout << "\nI4fmin = " << minI4fValue << std::endl;

        }

        
        solver.solveElectrophysiology (stim, t);
        solver.solveActivation (dt_activation);


        //============================================//
        // Load steps mechanics (activation & b.c.)
        //============================================//

        auto minActivationValue ( solver.activationModelPtr() -> fiberActivationPtr() -> minValue() );

        if ( k % mechanicsLoadstepIter == 0 && k % mechanicsCouplingIter != 0 && minActivationValue < activationLimit_loadstep )
        {
            if ( 0 == comm->MyPID() )
            {
                std::cout << "\n*****************************************************************";
                std::cout << "\nLoad step at time = " << t;
                std::cout << "\nMinimal activation value = " << minActivationValue;
                std::cout << "\n*****************************************************************\n";
            }
            
            // Linear b.c. extrapolation
            auto bcValuesLoadstep ( bcValues );
            bcValuesLoadstep[0] = bcValues[0] + ( bcValues[0] - bcValuesPre[0] ) * ( k % mechanicsCouplingIter ) / mechanicsCouplingIter;
            bcValuesLoadstep[1] = bcValues[1] + ( bcValues[1] - bcValuesPre[1] ) * ( k % mechanicsCouplingIter ) / mechanicsCouplingIter;
            
            if ( 0 == comm->MyPID() )
            {
                std::cout << "\n***************************************************************";
                std::cout << "\nLV-Pressure extrapolation from " <<  bcValues[0] << " to " <<  bcValuesLoadstep[0];
                std::cout << "\nRV-Pressure extrapolation from " <<  bcValues[1] << " to " <<  bcValuesLoadstep[1];
                std::cout << "\n***************************************************************\n\n";
            }
            
            // Load step mechanics
            solver.structuralOperatorPtr() -> data() -> dataTime() -> setTime(t);
            modifyFeBC(bcValuesLoadstep);
            solver.bcInterfacePtr() -> updatePhysicalSolverVariables();
            solver.solveMechanics();
        }
        
        
        //============================================//
        // Iterate mechanics / circulation
        //============================================//
        
        if ( k % mechanicsCouplingIter == 0 )
        {
            iter = 0;
            const double dt_circulation ( dt_mechanics / 1000 );
            solver.structuralOperatorPtr() -> data() -> dataTime() -> setTime(t);
            
            modifyFeBCPatches(t);
            
            //============================================//
            // 4th order Adam-Bashforth pressure extrapol.
            //============================================//
            
            for ( unsigned int i = ABcoef.size() - 1; i > 0; --i )
            {
                ABdplv(i) = ABdplv(i-1);
                ABdprv(i) = ABdprv(i-1);
            }
            
            ABdplv(0) = bcValues[0] - bcValuesPre[0];
            ABdprv(0) = bcValues[1] - bcValuesPre[1];

            bcValuesPre = bcValues;
            
            bcValues[0] += std::min( std::max( ABcoef.dot( ABdplv ) , - dpMax ) , dpMax );
            bcValues[1] += std::min( std::max( ABcoef.dot( ABdprv ) , - dpMax ) , dpMax );
            
            if ( 0 == comm->MyPID() )
            {
                std::cout << "\n***************************************************************";
                std::cout << "\nMinimal activation value = " << minActivationValue;
                std::cout << "\nLV-Pressure extrapolation from " <<  bcValuesPre[0] << " to " <<  bcValues[0];
                std::cout << "\nRV-Pressure extrapolation from " <<  bcValuesPre[1] << " to " <<  bcValues[1];
                std::cout << "\n***************************************************************\n\n";
            }

            
            //============================================//
            // Solve mechanics
            //============================================//
            modifyFeBC(bcValues);
            solver.bcInterfacePtr() -> updatePhysicalSolverVariables();
            solver.solveMechanics();
            
            VFeNew[0] = LV.volume(disp, dETFESpace, - 1);
            VFeNew[1] = RV.volume(disp, dETFESpace, 1);

            //============================================//
            // Solve circlation
            //============================================//
            circulationSolver.iterate(dt_circulation, bcNames, bcValues, iter);
            VCircNew[0] = VCirc[0] + dt_circulation * ( Q("la", "lv") - Q("lv", "sa") );
            VCircNew[1] = VCirc[1] + dt_circulation * ( Q("ra", "rv") - Q("rv", "pa") );

            //============================================//
            // Residual computation
            //============================================//
            R = VFeNew - VCircNew;
            printCoupling("Residual Computation");
            
            //============================================//
            // Newton iterations
            //============================================//
            while ( R.norm() > couplingError )
            {
                ++iter;
                
                //============================================//
                // Jacobian circulation
                //============================================//
                
                // Left ventricle
                circulationSolver.iterate(dt_circulation, bcNames, perturbedPressureComp(bcValues, pPerturbationCirc, 0), iter);
                VCircPert[0] = VCirc[0] + dt_circulation * ( Q("la", "lv") - Q("lv", "sa") );
                VCircPert[1] = VCirc[1] + dt_circulation * ( Q("ra", "rv") - Q("rv", "pa") );

                JCirc(0,0) = ( VCircPert[0] - VCircNew[0] ) / pPerturbationCirc;
                JCirc(1,0) = ( VCircPert[1] - VCircNew[1] ) / pPerturbationCirc;
                
                // Right ventricle
                circulationSolver.iterate(dt_circulation, bcNames, perturbedPressureComp(bcValues, pPerturbationCirc, 1), iter);
                VCircPert[0] = VCirc[0] + dt_circulation * ( Q("la", "lv") - Q("lv", "sa") );
                VCircPert[1] = VCirc[1] + dt_circulation * ( Q("ra", "rv") - Q("rv", "pa") );
                
                JCirc(0,1) = ( VCircPert[0] - VCircNew[0] ) / pPerturbationCirc;
                JCirc(1,1) = ( VCircPert[1] - VCircNew[1] ) / pPerturbationCirc;
                

                //============================================//
                // Jacobian fe
                //============================================//

                const bool jFeIter ( ! ( k % (couplingJFeIter * mechanicsCouplingIter) ) );
                const bool jFeSubIter ( ! ( (iter - couplingJFeSubStart) % couplingJFeSubIter) && iter >= couplingJFeSubStart );
                const bool jFeEmpty ( JFe.norm() == 0 );
                
                if ( jFeIter || jFeSubIter || jFeEmpty )
                {
                    JFe *= 0.0;
                    dispCurrent = disp;
                    
                    // Left ventricle
                    modifyFeBC(perturbedPressureComp(bcValues, pPerturbationFe, 0));
                    solver.bcInterfacePtr() -> updatePhysicalSolverVariables();
                    solver.solveMechanicsLin();
                    
                    VFePert[0] = LV.volume(disp, dETFESpace, - 1);
                    VFePert[1] = RV.volume(disp, dETFESpace, 1);

                    JFe(0,0) = ( VFePert[0] - VFeNew[0] ) / pPerturbationFe;
                    JFe(1,0) = ( VFePert[1] - VFeNew[1] ) / pPerturbationFe;
                    
                    disp = dispCurrent;
                    
                    // Right ventricle
                    modifyFeBC(perturbedPressureComp(bcValues, pPerturbationFe, 1));
                    solver.bcInterfacePtr() -> updatePhysicalSolverVariables();
                    solver.solveMechanicsLin();
                    
                    VFePert[0] = LV.volume(disp, dETFESpace, - 1);
                    VFePert[1] = RV.volume(disp, dETFESpace, 1);
                    
                    JFe(0,1) = ( VFePert[0] - VFeNew[0] ) / pPerturbationFe;
                    JFe(1,1) = ( VFePert[1] - VFeNew[1] ) / pPerturbationFe;
                    
                    disp = dispCurrent;
                }
                
                //============================================//
                // Update pressure b.c.
                //============================================//
                JR = JFe - JCirc;

                if ( JR.determinant() != 0 )
                {
                    dp = ( JR | R );
                    if ( iter > 5 ) dp *= 0.7;
                    if ( iter > 20 ) dp *= 0.5;
                    bcValues[0] -= std::min( std::max( dp(0) , - dpMax ) , dpMax );
                    bcValues[1] -= std::min( std::max( dp(1) , - dpMax ) , dpMax );
                }
                
                printCoupling("Pressure Update");

                //============================================//
                // Solve circulation
                //============================================//
                circulationSolver.iterate(dt_circulation, bcNames, bcValues, iter);
                VCircNew[0] = VCirc[0] + dt_circulation * ( Q("la", "lv") - Q("lv", "sa") );
                VCircNew[1] = VCirc[1] + dt_circulation * ( Q("ra", "rv") - Q("rv", "pa") );

                //============================================//
                // Solve mechanics
                //============================================//
                modifyFeBC(bcValues);
                solver.bcInterfacePtr() -> updatePhysicalSolverVariables();
                solver.solveMechanics();
                
                VFeNew[0] = LV.volume(disp, dETFESpace, - 1);
                VFeNew[1] = RV.volume(disp, dETFESpace, 1);

                //============================================//
                // Residual update
                //============================================//
                R = VFeNew - VCircNew;
                printCoupling("Residual Update");
            }
 
            if ( 0 == comm->MyPID() )
            {
                std::cout << "\n*****************************************************************";
                std::cout << "\nCoupling converged after " << iter << " iteration" << ( iter > 1 ? "s" : "" );
                std::cout << "\n*****************************************************************\n\n";
            }
            
            //============================================//
            // Update volume variables
            //============================================//
            VCirc = VCircNew;
            VFe = VFeNew;
            
            //============================================//
            // Export circulation solution
            //============================================//
            if ( 0 == comm->MyPID() ) circulationSolver.exportSolution( circulationOutputFile );
            
        }
        
        //============================================//
        // Export FE-solution
        //============================================//
        bool save ( std::abs(std::remainder(t, dt_save)) < 0.01 );
        if ( save ) solver.saveSolution(t);

    }

    
    //============================================//
    // Close all exporters
    //============================================//
    solver.closeExporters();
    

#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return 0;
}
