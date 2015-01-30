#include <lifev/core/LifeV.hpp>

#include <lifev/bc_interface/3D/bc/BCInterface3D.hpp>
#include <lifev/core/fem/GradientRecovery.hpp>

#include <lifev/core/mesh/MeshLoadingUtility.hpp>
#include <lifev/core/fem/BCManage.hpp>

#include <lifev/eta/fem/ETFESpace.hpp>
#include <lifev/eta/expression/Integrate.hpp>

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

#include <lifev/em/solver/EMData.hpp>


using namespace LifeV;

void createListFromGetPot (Teuchos::ParameterList& solverList, const GetPot& dataFile);


Real pos (const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID& i)
{
	if (i == 0)
		return x;
	else if ( i == 1)
		return y;
	else if (i == 2)
		return z;
	else
		return 0.0;
}

template <typename T> int sgn(T val)
{
    return (  (T(0) < val) - (val < T(0)) );
}


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

    typedef BCHandler                                           bc_Type;
    typedef boost::shared_ptr< bc_Type >                        bcPtr_Type;

    typedef StructuralOperator< RegionMesh<LinearTetra> >      physicalSolver_Type;
    typedef BCInterface3D< bc_Type, physicalSolver_Type >      bcInterface_Type;
    typedef boost::shared_ptr< bcInterface_Type >              bcInterfacePtr_Type;


//    typedef boost::shared_ptr< bcInterface_Type >                  bcInterfacePtr_Type;
//    typedef MeshUtility::MeshTransformer<mesh_Type>                meshTransformer_Type;

    typedef boost::function < Real (const Real&  t,
                                    const Real&  x,
                                    const Real&  y,
                                    const Real&  z,
                                    const ID&    i ) >   function_Type;

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

    if( comm->MyPID() == 0 )
    {
    	std::cout << meshName << "\n";
    }

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

    typedef ETFESpace< RegionMesh<LinearTetra>, MapEpetra, 3, 1 >       scalarETFESpace_Type;
    typedef boost::shared_ptr<scalarETFESpace_Type>                      scalarETFESpacePtr_Type;


    std::string dOrder =  dataFile ( "solid/space_discretization/order", "P1");
    solidFESpacePtr_Type dFESpace ( new solidFESpace_Type (localSolidMesh, dOrder, 3, comm) );
    solidFESpacePtr_Type uFESpace ( new solidFESpace_Type (localSolidMesh, dOrder, 1, comm) );
    solidETFESpacePtr_Type dETFESpace ( new solidETFESpace_Type (localSolidMesh, & (dFESpace->refFE() ), & (dFESpace->fe().geoMap() ), comm) );
    scalarETFESpacePtr_Type uETFESpace ( new scalarETFESpace_Type (localSolidMesh, & (uFESpace->refFE() ), & (uFESpace->fe().geoMap() ), comm) );

    if ( comm->MyPID() == 0 )
    {
        std::cout << " Done!" << std::endl;
    }


    /// STEP I
    /// CREATE FIBERS


    //compute x, y, z
    function_Type coords = &pos;
    vectorPtr_Type position( new vector_Type ( dFESpace -> map() ) );
    dFESpace->interpolate ( static_cast<FESpace<RegionMesh<LinearTetra>, MapEpetra>::function_Type> (coords),*position, 0.0);
    //*************************************************************//
    // We asseble the stiffness matrix using expression template.
    // For more details, look at the ETA tutorial.
    //*************************************************************//

    boost::shared_ptr<matrix_Type> systemMatrix (new matrix_Type ( uFESpace->map() ) );

    *systemMatrix *= 0.0;
    {
        using namespace ExpressionAssembly;


        integrate (  elements (uETFESpace->mesh() ),
                     quadRuleTetra4pt,
                     uETFESpace,
                     uETFESpace,
                     dot ( grad (phi_i) , grad (phi_j) )
                  )
                >> systemMatrix;
    }

    systemMatrix->globalAssemble();

    //set bc
    bcInterfacePtr_Type                     BC ( new bcInterface_Type() );
    BC->createHandler();
    BC->fillHandler ( data_file_name, "problem" );
    BC->handler()->bcUpdate ( *uFESpace->mesh(), uFESpace->feBd(), uFESpace->dof() );
    //preconditioner
    typedef LifeV::Preconditioner             basePrec_Type;
    typedef boost::shared_ptr<basePrec_Type>  basePrecPtr_Type;
    typedef LifeV::PreconditionerIfpack           prec_Type;
    typedef boost::shared_ptr<prec_Type>      precPtr_Type;

    prec_Type* precRawPtr;
    basePrecPtr_Type precPtr;
    precRawPtr = new prec_Type;
    precRawPtr->setDataFromGetPot ( dataFile, "prec" );
    precPtr.reset ( precRawPtr );
    Teuchos::ParameterList solverList;
    createListFromGetPot (solverList, dataFile);

    LinearSolver linearSolver;
    linearSolver.setCommunicator (comm);
    linearSolver.setParameters ( solverList );
    linearSolver.setPreconditioner ( precPtr );

    vectorPtr_Type rhs (new vector_Type ( uFESpace -> map() ) );
    *rhs *= 0.0;
    rhs -> globalAssemble();

    bcManage ( *systemMatrix, *rhs, *uFESpace->mesh(), uFESpace->dof(), *BC -> handler(), uFESpace->feBd(), 1.0, 0.0 );
    vectorPtr_Type solution ( new vector_Type ( uFESpace -> map()) );


    linearSolver.setOperator (systemMatrix);
    linearSolver.setRightHandSide (rhs);
    linearSolver.solve (solution);


    vectorPtr_Type sx (new vector_Type ( uFESpace -> map() ) );
    vectorPtr_Type sy (new vector_Type ( uFESpace -> map() ) );
    vectorPtr_Type sz (new vector_Type ( uFESpace -> map() ) );

    *sx = GradientRecovery::ZZGradient (uFESpace, *solution, 0);
    *sy = GradientRecovery::ZZGradient (uFESpace, *solution, 1);
    *sz = GradientRecovery::ZZGradient (uFESpace, *solution, 2);

    vectorPtr_Type rbFiber ( new vector_Type ( dFESpace -> map() ) );
    vectorPtr_Type rbSheet ( new vector_Type ( dFESpace -> map() ) );
    UInt d =  rbFiber->epetraVector().MyLength();
    UInt nComponentLocalDof = d / 3;

    std::cout << "\nNum dof: " << nComponentLocalDof;
    for ( int l (0); l < nComponentLocalDof; l++)
	{
		UInt iGID = rbFiber->blockMap().GID (l);
        UInt jGID = rbFiber->blockMap().GID (l + nComponentLocalDof);
        UInt kGID = rbFiber->blockMap().GID (l + 2 * nComponentLocalDof);

    	std::cout<< "\nLID: " << l << ", iGID: " <<  iGID << ",jGID: " << jGID << ", kGID: " << kGID ;

//        meshFull->point(iGID).MeshVertex::showMe(true);
		Real x = (*position)[iGID];
		Real y = (*position)[jGID];
		Real z = (*position)[kGID];

		//		std::cout << "x: " << x << ", y: " << y << ", z: " << z << "\n";

		Real t = (*solution)[iGID];

		Real alpha_degrees = 90.0 - 180 * t;
		Real alpha = (alpha_degrees * 3.1415 ) / 180.0;

		Real rl = 1.7 + 0.3 * t;
		Real rs =  0.7 + 0.3 * t;

		Real CosU = z / rl;
		Real SinU = std::sqrt(1.0 - CosU * CosU );

		Real CosV(1.);
		if( std::abs(x) >= 1e-16 )
		{
			Real TgV = y / x;
			Real V = std::atan(TgV);
			CosV = sgn(x) * std::cos(V);
		}
		else if(std::abs(x) < 1e-16 && std::abs(y) > 1e-16)
		{
			CosV = sgn(x);
		}
		else
		{
			CosV = 1.0;
		}


		Real SinV = sgn(y) * std::sqrt(1.0 - CosV * CosV );

		Real dxdu1 = rs * CosU * CosV;
		Real dxdu2 = rs * CosU * SinV;
		Real dxdu3 = - rl * SinU;
		Real norm1 = std::sqrt( dxdu1 * dxdu1 + dxdu2 * dxdu2 + dxdu3 * dxdu3 );

		Real dxdv1 = - SinV;
		Real dxdv2 = CosV;
		Real dxdv3 = 0.0;
		Real norm2 = 1.0;

		Real SinA = std::sin(alpha);
		Real CosA = std::cos(alpha);


		(*rbFiber)[iGID] = dxdu1 * SinA / norm1 + dxdv1 * CosA / norm2;
		(*rbFiber)[jGID] = dxdu2 * SinA / norm1 + dxdv2 * CosA / norm2;
		(*rbFiber)[kGID] = dxdu3 * SinA / norm1 + dxdv3 * CosA / norm2;

		Real fiberNorm = std::sqrt( (*rbFiber)[iGID] * (*rbFiber)[iGID] + (*rbFiber)[jGID] * (*rbFiber)[jGID] + (*rbFiber)[kGID] * (*rbFiber)[kGID]);
		if(fiberNorm >= 1e-13 )
		{
			(*rbFiber)[iGID] = (*rbFiber)[iGID] / fiberNorm;
			(*rbFiber)[jGID] = (*rbFiber)[jGID] / fiberNorm;
			(*rbFiber)[kGID] = (*rbFiber)[kGID] / fiberNorm;
		}
		else
		{
			(*rbFiber)[iGID] = 1.0;
			(*rbFiber)[jGID] = 0.0;
			(*rbFiber)[kGID] = 0.0;
		}


		Real cdot = (*sx) [iGID] * (*rbFiber)[iGID] + (*sy) [iGID] * (*rbFiber)[jGID]  + (*sz) [iGID] * (*rbFiber)[kGID] ;
        (*rbSheet) [iGID] = (*sx) [iGID] - cdot * (*rbFiber)[iGID];
        (*rbSheet) [jGID] = (*sy) [iGID] - cdot * (*rbFiber)[jGID];
        (*rbSheet) [kGID] = (*sz) [iGID] - cdot * (*rbFiber)[kGID];

	}

    rbFiber->spy("vector");
    ElectrophysiologyUtility::normalize(*rbSheet,1);

    //===========================================================
    //===========================================================
    //              BOUNDARY CONDITIONS Part I
    //===========================================================
    //===========================================================
    if ( comm->MyPID() == 0 )
    {
        std::cout << "setup bc ... ";
    }


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

    solid.EMMaterial()->setFiberVectorPtr(rbFiber);
    solid.EMMaterial()->setSheetVectorPtr(rbSheet);

    std::cout << "\nSheets norm = " << solid.EMMaterial()->sheetVectorPtr()->norm2();
    std::cout << "\nFibers norm = " << solid.EMMaterial()->fiberVectorPtr()->norm2();

    if(solid.EMMaterial()->sheetVectorPtr()) std::cout << "I have sheets!\n";
    else
    {
        std::cout << "I don't have sheets!\n";
        return 0.0;
    }



    *( solid.EMMaterial()->fiberActivationPtr() ) = 1.0;

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
    solid.setBCFlag(static_cast<ID>(flag));
    Real FinalPressure =  dataFile ( "solid/boundary_conditions/pressure", 1000.0);

    vectorPtr_Type boundaryVectorPtr(new vector_Type ( solid.displacement().map(), Repeated ) );
//    ComputeBC<solidETFESpace_Type>(localSolidMesh, solid.displacement(), boundaryVectorPtr, dETFESpace, flag);
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
//    exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::ScalarField, "activation", aFESpace, solid.EMMaterial()->activationPtr(), UInt (0) );
	exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, "fibers", dFESpace, solid.EMMaterial()->fiberVectorPtr(), UInt (0) );
	exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, "sheets", dFESpace, solid.EMMaterial()->sheetVectorPtr(), UInt (0) );
        exporter->postProcess ( 0 );

    vectorPtr_Type dispRep( new vector_Type( solid.displacement(), Repeated ) );
    vectorPtr_Type solRep( new vector_Type( *solution, Repeated) );
    vectorPtr_Type fiberRep( new vector_Type( *solid.EMMaterial()->fiberVectorPtr(), Repeated ) );
    vectorPtr_Type sheetRep( new vector_Type( *solid.EMMaterial()->sheetVectorPtr(), Repeated ) );

    ExporterVTK< mesh_Type > exporterVTK;
    exporterVTK.setMeshProcId ( localSolidMesh, comm -> MyPID() );
    exporterVTK.setPostDir (problemFolder);
    exporterVTK.setPrefix ("DeformedConfiguration");
    exporterVTK.addVariable ( ExporterData<mesh_Type>::VectorField,  "Displacement", dFESpace, dispRep, UInt (0) );
    exporterVTK.addVariable ( ExporterData<RegionMesh<LinearTetra> >::ScalarField, "potential", uFESpace, solRep, UInt (0) );
	exporterVTK.addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, "fibers", dFESpace,fiberRep, UInt (0) );
	exporterVTK.addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, "sheets", dFESpace, sheetRep, UInt (0) );
    exporterVTK.postProcess ( 0 );

    //===========================================================
    //===========================================================
    //         SOLVE
    //===========================================================
    //===========================================================
    Real dt =  dataFile ( "solid/time_discretization/timestep", 0.1);
    Real endTime = dataFile ( "solid/time_discretization/endtime", 1.0);
    UInt activeRamp = dataFile ( "importer/activeramp", 1);


    if ( comm->MyPID() == 0 )
    {
		std::cout << "\n=====================================================\n";
		std::cout << "============= SOLVING PART I" ;
		std::cout << "\n=====================================================\n";
    }

    Real time(importTime);
    if(0 == activeRamp) time = 0.0;

    Int iter = 0;
    Int itermax = static_cast<Int>( endTime / dt );
    std::cout << "Starting from time: " << importTime << ", and finishing at: " << endTime << ", with " << itermax << " iterations\n";
    Real Tmax =  dataFile ( "solid/physics/Tmax", 60000);
    for ( ; iter < itermax ; iter++ )
    {

    	time += dt;
    	if( 1 == activeRamp ) emdata.setSolidParameter("MaxActiveTension", Tmax * time );
    	else  emdata.setSolidParameter("MaxActiveTension", Tmax );

        solid.EMMaterial()->setParameters(emdata);
    	if( comm ->MyPID() == 0)
    	{
        solid.EMMaterial()->showMaterialParameters();
    	}

//    	solid.data() -> dataTime() -> updateTime();

    	if( comm ->MyPID() == 0)
    	{
			std::cout << "\n=====================================================\n";
			std::cout << "============= TIME: " << time ;
			std::cout << "\n=====================================================\n";
	        std::cout << "\nPressure: " << FinalPressure * time << ", Active Tension: " << Tmax * time <<"\n";
    	}
//        ComputeBC<solidETFESpace_Type>(localSolidMesh, solid.displacement(), boundaryVectorPtr, dETFESpace, flag);
//        *boundaryVectorPtr *= (FinalPressure * time);
//        Real norm = boundaryVectorPtr -> norm2();
//        bcVectorPtr.reset( new BCVector (*boundaryVectorPtr, dFESpace -> dof().numTotalDof(), 0 ) );
//	    solidBC -> handler() -> modifyBC(flag, *bcVectorPtr);
    	/////
        solid.setLVPressureBC(FinalPressure * time);
//    	solid.iterate ( solidBC -> handler() );
    	solid.iterate ( solidBC -> handler(), true );

        exporter->postProcess ( time );

//        //===========================================================
//        //         SOLVE PART II: Iterate till convergence
//        //===========================================================
//        if ( comm->MyPID() == 0 )
//        {
//    		std::cout << "\n=====================================================\n";
//    		std::cout << "============= SOLVING PART II" ;
//    		std::cout << "\n=====================================================\n";
//        }
//
//        Real tolN = 1e-3;
//	Real resN = 1000.;
//        Real oldres = 4*resN;
//        int iterationN(0);
//        Real timeN(time);
//        Real dtN(dt/200.);
//        while( resN > tolN )
//        {
//	    	oldres = resN;
//
//        	iterationN++;
//        	timeN += dtN;
//
//            if ( comm->MyPID() == 0 )
//            {
//    			std::cout << "\n=====================================================\n";
//    			std::cout << "============= TIME: " << timeN << ", Iteration: " << iterationN;
//    			std::cout << "\n=====================================================\n";
//            }
//
//        	VectorEpetra disp_tn(solid.displacement());
//
//            ComputeBC<solidETFESpace_Type>(localSolidMesh, solid.displacement(), boundaryVectorPtr, dETFESpace, flag);
//            *boundaryVectorPtr *= (FinalPressure * time);
//
//            Real norm = boundaryVectorPtr -> norm2();
//            std::cout << "Norm BC: " << norm << ", Pressure: " << FinalPressure * time << ", Active Tension: " << Tmax * time <<"\n";
//            bcVectorPtr.reset( new BCVector (*boundaryVectorPtr, dFESpace -> dof().numTotalDof(), 0 ) );
//    	    solidBC -> handler() -> modifyBC(flag, *bcVectorPtr);
//
//        	solid.iterate ( solidBC -> handler() );
//
//            exporter->postProcess ( timeN );
//
//            disp_tn -= solid.displacement();
//
//            resN = disp_tn.norm2();
//        Real dispNorm =  solid.displacement().norm2();
//
//        if(dispNorm > 1e-14) resN /= dispNorm;
//
//
//        if(oldres < resN )
//        {
//        	solid.displacement() = solid.displacement() + disp_tn;
//        	break;
//        }
//
//            if ( comm->MyPID() == 0 )
//            {
//    			std::cout << "\n============= Residual: " << resN ;
//    			std::cout << "\n=====================================================\n";
//            }
//        }

    }

//    //===========================================================
//    //         SOLVE PART II: Iterate till convergence
//    //===========================================================
//    if ( comm->MyPID() == 0 )
//    {
//		std::cout << "\n=====================================================\n";
//		std::cout << "============= SOLVING PART II" ;
//		std::cout << "\n=====================================================\n";
//    }
//
//    Real tol = 1e-3;
//    Real res = 1000.;
//    Real oldres = 4*res;
//    int iteration(0);
//    emdata.setParameter("MaxActiveTension", Tmax);
//    solid.EMMaterial()->setParameters(emdata);
//    solid.EMMaterial()->showMaterialParameters();
//
//    while( res > tol )
//    {
//    	oldres = res;
//    	iteration++;
//    	time += dt;
//
//        if ( comm->MyPID() == 0 )
//        {
//			std::cout << "\n=====================================================\n";
//			std::cout << "============= TIME: " << time << ", Iteration: " << iteration;
//			std::cout << "\n=====================================================\n";
//        }
//    	VectorEpetra disp_tn(solid.displacement());
//
//        ComputeBC<solidETFESpace_Type>(localSolidMesh, solid.displacement(), boundaryVectorPtr, dETFESpace, flag);
//        *boundaryVectorPtr *= (FinalPressure);
//
//        Real norm = boundaryVectorPtr -> norm2();
//        std::cout << "Norm BC: " << norm << ", Pressure: " << FinalPressure << ", Active Tension: " << Tmax <<"\n";
//        bcVectorPtr.reset( new BCVector (*boundaryVectorPtr, dFESpace -> dof().numTotalDof(), 0 ) );
//	    solidBC -> handler() -> modifyBC(flag, *bcVectorPtr);
//
//    	solid.iterate ( solidBC -> handler() );
//
//        exporter->postProcess ( time );
//
//        disp_tn -= solid.displacement();
//
//        Real dispNorm =  solid.displacement().norm2();
//
//        res = disp_tn.norm2();
//        if(dispNorm > 1e-14) res /= dispNorm;
//
//
//        if(oldres < res )
//        {
//        	solid.displacement() = solid.displacement() + disp_tn;
//        	break;
//        }
//        if ( comm->MyPID() == 0 )
//        {
//			std::cout << "\n============= Residual: " << res ;
//			std::cout << "\n=====================================================\n";
//        }
//    }

    //===========================================================
    //===========================================================
    //              CLOSE EXPORTER
    //===========================================================
    //===========================================================
    exporter -> closeFile();


    //*************************************************************//
    // We save the potential field just computed on a file
    // using VTK.
    //*************************************************************//
    exporterVTK.closeFile();

#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return 0;
}



//template<typename space> void ComputeBC( const boost::shared_ptr<RegionMesh<LinearTetra> > localMesh,
//						const VectorEpetra& disp,
//						boost::shared_ptr<VectorEpetra> bcVectorPtr,
//						const ETFESpacePtr_Type dETFESpace,
//						int bdFlag)
//{
//
//	*bcVectorPtr *= 0.0;
//
//    MatrixSmall<3,3> Id;
//    Id(0,0) = 1.; Id(0,1) = 0., Id(0,2) = 0.;
//    Id(1,0) = 0.; Id(1,1) = 1., Id(1,2) = 0.;
//    Id(2,0) = 0.; Id(2,1) = 0., Id(2,2) = 1.;
//
//	boost::shared_ptr<VectorEpetra> intergral( new VectorEpetra( disp.map() ) );
//
//	{
//	  	using namespace ExpressionAssembly;
//
//	  	auto I = value(Id);
//	  	auto Grad_u = grad( dETFESpace, disp, 0);
//	  	auto F =  Grad_u + I;
//	  	auto FmT = minusT(F);
//	  	auto J = det(F);
//
//    	QuadratureBoundary myBDQR (buildTetraBDQR (quadRuleTria7pt) );
//
//        integrate ( boundary ( localMesh, bdFlag),
//        		    myBDQR,
//        		    dETFESpace,
//        		    J * dot( FmT * Nface,  phi_i)
//                  ) >> bcVectorPtr;
//
//        bcVectorPtr -> globalAssemble();
//
//
//	}
//}
//Implementation of the predeclared functions:
//
//Read the parameters from a datafile and
// put them in a Teuchos:ParameterList
void createListFromGetPot (Teuchos::ParameterList& solverList, const GetPot& dataFile)
{
    std::string solverName   = dataFile ( "problem/solver/solver_name", "AztecOO");
    std::string solver       = dataFile ( "problem/solver/solver", "gmres");
    std::string conv         = dataFile ( "problem/solver/conv", "rhs");
    std::string scaling      = dataFile ( "problem/solver/scaling", "none");
    std::string output       = dataFile ( "problem/solver/output", "all");
    Int maxIter              = dataFile ( "problem/solver/max_iter", 200);
    Int maxIterForReuse      = dataFile ( "problem/solver/max_iter_reuse", 250);
    Int kspace               = dataFile ( "problem/solver/kspace", 100);
    Int orthog               = dataFile ( "problem/solver/orthog", 0);
    Int auxvec               = dataFile ( "problem/solver/aux_vec", 0);
    double tol               = dataFile ( "problem/solver/tol", 1e-10);
    bool reusePreconditioner = dataFile ( "problem/solver/reuse", true);
    bool quitOnFailure       = dataFile ( "problem/solver/quit", false);
    bool silent              = dataFile ( "problem/solver/silent", false);

    solverList.set ("Solver Type", solverName);
    solverList.set ("Maximum Iterations", maxIter);
    solverList.set ("Max Iterations For Reuse", maxIterForReuse);
    solverList.set ("Reuse Preconditioner", reusePreconditioner);
    solverList.set ("Quit On Failure", quitOnFailure);
    solverList.set ("Silent", silent);
    solverList.sublist ("Solver: Operator List").sublist ("Trilinos: AztecOO List").set ("solver", solver);
    solverList.sublist ("Solver: Operator List").sublist ("Trilinos: AztecOO List").set ("conv", conv);
    solverList.sublist ("Solver: Operator List").sublist ("Trilinos: AztecOO List").set ("scaling", scaling);
    solverList.sublist ("Solver: Operator List").sublist ("Trilinos: AztecOO List").set ("output", output);
    solverList.sublist ("Solver: Operator List").sublist ("Trilinos: AztecOO List").set ("tol", tol);
    solverList.sublist ("Solver: Operator List").sublist ("Trilinos: AztecOO List").set ("max_iter", maxIter);
    solverList.sublist ("Solver: Operator List").sublist ("Trilinos: AztecOO List").set ("kspace", kspace);
    solverList.sublist ("Solver: Operator List").sublist ("Trilinos: AztecOO List").set ("orthog", orthog);
    solverList.sublist ("Solver: Operator List").sublist ("Trilinos: AztecOO List").set ("aux_vec", auxvec);;
}
