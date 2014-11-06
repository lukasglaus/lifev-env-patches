//@HEADER
/*
 *******************************************************************************

    Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
    Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University

    This file is part of LifeV.

    LifeV is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    LifeV is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with LifeV.  If not, see <http://www.gnu.org/licenses/>.

 *******************************************************************************
 */
//@HEADER

/*!
    @file
    @brief Generation muscular fibers and sheets

    @author Simone Rossi <simone.rossi@epfl.ch>
    @maintainer Simone Palamara <palamara.simone@gmail.com>
    @date 31-01-2014


    Generation of the muscular fibers and sheets on a generic
    geometry representing the left or right ventricle, generated
    according to geometrical rules based on anatomical knowledge.
    For more details about the method see [S.Rossi et al,European Journal of
    Mechanics A/Solids (2013), http://dx.doi.org/10.1016/j.euromechsol.2013.10.009]

 */

//#include <Epetra_ConfigDefs.h>

// ------------------------------------------------------------------------------
//  Include MPI for parallel simulations
// ------------------------------------------------------------------------------
#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

// ------------------------------------------------------------------------------
//  To set up the linear solver we need a Teuchos parameter list
// ------------------------------------------------------------------------------
#include <Teuchos_ParameterList.hpp>

// ------------------------------------------------------------------------------
//  Needed to generate the ouput folder
// ------------------------------------------------------------------------------
#include <sys/stat.h>

// ------------------------------------------------------------------------------
// BCInterface is the interface that between the datafile and the
// boundary conditions. We will create a ummy physical solver in order to use it.
// ------------------------------------------------------------------------------
#include <lifev/bc_interface/3D/bc/BCInterface3D.hpp>
#include <lifev/bc_interface/core/solver/EmptyPhysicalSolver.hpp>

// ------------------------------------------------------------------------------
//  Usefule utility to load the mesh in one line.
// Not working with partitioned meshes
// ------------------------------------------------------------------------------
#include <lifev/core/mesh/MeshLoadingUtility.hpp>

// ------------------------------------------------------------------------------
//  The sheets are defined as the gradient of a scalar potential, we therefore
// need to be able to compute the gradient in the nodes.
// ------------------------------------------------------------------------------
#include <lifev/core/fem/GradientRecovery.hpp>
#include <lifev/core/fem/BCManage.hpp>

// ------------------------------------------------------------------------------
// To solve the laplacian, use the linear solver (AztecOO or Belos) with ML
// ------------------------------------------------------------------------------
#include <lifev/core/algorithm/LinearSolver.hpp>
#include <lifev/core/algorithm/PreconditionerML.hpp>

// ------------------------------------------------------------------------------
// The laplacian is assembled using Expression Template Assembly
// ------------------------------------------------------------------------------
#include <lifev/eta/fem/ETFESpace.hpp>
#include <lifev/eta/expression/Integrate.hpp>

// ------------------------------------------------------------------------------
//  Cannot save without HDF5 !!! Please make sure you have it
// ------------------------------------------------------------------------------
#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5.hpp>
#endif

// ------------------------------------------------------------------------------
//  In this file there are a bunch of useful functions.
// Here we use it only to normalize a vector.
// ------------------------------------------------------------------------------
#include <lifev/electrophysiology/util/HeartUtility.hpp>

// ---------------------------------------------------------------
// As usual, we work in the LifeV namespace. Moreover,
// to make the code more readable, we also make typedefs for the mesh type,
// matrix type, vector type, boundary condition
// ---------------------------------------------------------------

using namespace LifeV;

// ---------------------------------------------------------------
// We typedef some common type we will frequently
// ---------------------------------------------------------------

typedef RegionMesh<LinearTetra>                             mesh_Type;
typedef boost::shared_ptr< mesh_Type >                      meshPtr_Type;

typedef MatrixEpetra<Real>                                  matrix_Type;
typedef boost::shared_ptr< matrix_Type >                    matrixPtr_Type;

typedef VectorEpetra                                        vector_Type;
typedef boost::shared_ptr< vector_Type >                    vectorPtr_Type;



typedef FESpace< mesh_Type, MapEpetra >                     fespace_Type;
typedef boost::shared_ptr<fespace_Type >                    fespacePtr_Type;


// ---------------------------------------------------------------
// In order to keep the code more readble I created a couple
// of auxiliary functions.
// In this test we only have the datafile that will be
// a GetPot object. To set up the linear solver we need a Teuchos::ParameterList
// that typically reads  xml files. Therefore,
// the createListFromGetPot function create the requested Teuchos list
// from the datafile we have.
// In the end we will export three vector field in three different
// files. To avoid code repetition, I created this function
// that exports the requested vecotr fields.
// ---------------------------------------------------------------
void createListFromGetPot (Teuchos::ParameterList& solverList, const GetPot& dataFile);
void exportVectorField (boost::shared_ptr<Epetra_Comm> comm,
                        meshPtr_Type mesh,
                        fespacePtr_Type fespace,
                        vectorPtr_Type vector,
                        std::string postDir,
                        std::string outputName,
                        std::string hdf5name);


Real fiberField(const Real& /*t*/, const Real&  X, const Real& Y, const Real& /*Z*/, const ID& i);


// ------------------------------------------------------------------------------
// Dummy class for to be used with the BCInterfac3D.
// We could have used the ElectroETAMonodomainSolver as physical solver,
// but in this way, this test is totally independent.
// ------------------------------------------------------------------------------

Real fzero (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return 0.0;
}


template <typename T> int sgn(T val)
{
    return (T(0) < val) - (val < T(0));
}

// Starting ...
int main ( int argc, char** argv )
{

    // ---------------------------------------------------------------
    //  In parallel? Yes, we dare!!!
    // ---------------------------------------------------------------

#ifdef HAVE_MPI
    MPI_Init (&argc, &argv);
    boost::shared_ptr<Epetra_Comm> Comm (new Epetra_MpiComm (MPI_COMM_WORLD) );
#else
    boost::shared_ptr<Epetra_Comm> Comm (new Epetra_SerialComm);
#endif

    typedef BCHandler                                           bc_Type;
    typedef boost::shared_ptr< bc_Type >                        bcPtr_Type;

    typedef EmptyPhysicalSolver<VectorEpetra>               physicalSolver_Type;
    typedef BCInterface3D< bc_Type, physicalSolver_Type >       bcInterface_Type;

    typedef boost::shared_ptr< bcInterface_Type >                  bcInterfacePtr_Type;
    typedef MeshUtility::MeshTransformer<mesh_Type>                meshTransformer_Type;

    //*************************************************************//
    // We create, as usual, the output folder where
    // we are going to save the solutions of
    // the simulation.
    // It requires to append "-o OutputFolderName" flag
    // when you laucnh the executable, e.g.
    // mpirun -n 2 Electrophysiology_thisTestName -o SolutionFolder
    // By default the name of the folder will be "Output"
    //*************************************************************//
    GetPot commandLine ( argc, argv );
    std::string problemFolder = commandLine.follow ( "Output", 2, "-o", "--output" );
    // Create the problem folder
    if ( problemFolder.compare ("./") )
    {
        problemFolder += "/";

        if ( Comm->MyPID() == 0 )
        {
            mkdir ( problemFolder.c_str(), 0777 );
        }
    }

    //*************************************************************//
    // We create the datafile. The datafile is passed through the flag
    // "-f datafileName" when launching the executable.
    // By default the name of the datafiler is "data".
    // Usually you should not bother about this flag, as the
    // datafile has already the default filename.
    // Remember about this, though, if you are going to change the filename
    // or if you want to add another datafile and run your simulation
    // with that one.
    //*************************************************************//
    GetPot command_line (argc, argv);
    const std::string data_file_name = command_line.follow ("data", 2, "-f", "--file");
    GetPot dataFile (data_file_name);

    //*************************************************************//
    // We specified in the datafile the name and the path of the
    // mesh we are going to create the fibers on.
    //*************************************************************//

    std::string meshName = dataFile ( "problem/space_discretization/mesh_file", "" );
    std::string meshPath = dataFile ( "problem/space_discretization/mesh_dir", "./" );

    //*************************************************************//
    // Here we create a pointer to the mesh and then we load it.
    // Notice how cool is to load the mesh in just one line!!!
    // Note that the pointer will point to the partitioned mesh
    // If you would like to keep informations about the full mesh
    // create another meshPtr_Type and call for example
    // MeshUtility::loadMesh (meshPart, meshFull, meshName, meshPath);
    //*************************************************************//
    meshPtr_Type meshPart (new mesh_Type ( Comm ) );
    meshPtr_Type meshFull (new mesh_Type ( Comm ) );
    MeshUtility::loadMesh (meshPart, meshFull, meshName, meshPath);

    //*************************************************************//
    // Here we define the finite element spaces. In particular
    // - uSpace is the finite element space for Expression Template Assembly
    // - uFESpace is the usual finite element space required for
    //   boundary conditions and for the exporter
    // - vectorESpace is the space for the vectorial fields (sheets and fibers)
    //   we will define. Even if we don't need to do finite element
    //   operation on them, this fespace will be required to export
    //   the solution.
    //*************************************************************//

    std::string order = dataFile ( "problem/space_discretization/order", "P1" );

    boost::shared_ptr<ETFESpace< mesh_Type, MapEpetra, 3, 1 > > uSpace;
    if(order == "P2")
    {
		uSpace.reset( new ETFESpace< mesh_Type, MapEpetra, 3, 1 > (meshPart, &feTetraP2, Comm) );
    }
    else
    {
		uSpace.reset( new ETFESpace< mesh_Type, MapEpetra, 3, 1 > (meshPart, &feTetraP1, Comm) );
    }

    fespacePtr_Type uFESpace ( new FESpace< mesh_Type, MapEpetra > (meshPart, order, 1, Comm) );

    fespacePtr_Type vectorFESpace ( new FESpace< mesh_Type, MapEpetra > (meshPart, order, 3, Comm) );

    //*************************************************************//
    // We asseble the stiffness matrix using expression template.
    // For more details, look at the ETA tutorial.
    //*************************************************************//

    boost::shared_ptr<matrix_Type> systemMatrix (new matrix_Type ( uSpace->map() ) );

    *systemMatrix *= 0.0;
    {
        using namespace ExpressionAssembly;


        integrate (  elements (uSpace->mesh() ),
                     quadRuleTetra4pt,
                     uSpace,
                     uSpace,
                     dot ( grad (phi_i) , grad (phi_j) )
                  )
                >> systemMatrix;
    }

    systemMatrix->globalAssemble();

    //*************************************************************//
    // Setting up the boundary conditions is always an issue.
    // Fortunately Cristiano Malossi implemented a way to read
    // the boundary condition from the datafile. This is achieved
    // by using the BCInterface class.
    // The BC is template on a physicalSolver_Type.
    // We have not created a specific BCInterface for the
    // monodomain model (as usually one imposes homogeneous Neumann
    // conditions). Therefore here I created a default physical solver
    // in the bc interface module which will be used to solve this
    // simple laplacian problem. Check the typedef.
    // Then to set up the boundary conditions, we need to
    // 1 -  create the BCHandler
    // 2 - fill the handler with the boundary conditions specified
    //     in the datafile
    // 3 - update the boundary condtions
    // For more information on how to use the BCInterface refer
    // to that module.
    //*************************************************************//
    bcInterfacePtr_Type                     BC ( new bcInterface_Type() );
    BC->createHandler();
    BC->fillHandler ( data_file_name, "problem" );
    BC->handler()->bcUpdate ( *uFESpace->mesh(), uFESpace->feBd(), uFESpace->dof() );

    //**********************************************************//
    //  We are going to solve the laplace equation with an
    // iterative solver. We precondition the system with Ifpack.
    // The parameters of the preconditioner are in the datafile.
    // If you want to use ML, change the precType to
    // LifeV::PreconditionerML
    //**********************************************************//

    typedef LifeV::Preconditioner             basePrec_Type;
    typedef boost::shared_ptr<basePrec_Type>  basePrecPtr_Type;
    typedef LifeV::PreconditionerIfpack           prec_Type;
    typedef boost::shared_ptr<prec_Type>      precPtr_Type;

    prec_Type* precRawPtr;
    basePrecPtr_Type precPtr;
    precRawPtr = new prec_Type;
    precRawPtr->setDataFromGetPot ( dataFile, "prec" );
    precPtr.reset ( precRawPtr );

    //**********************************************************//
    // As detailed above to setup the linear system we need a
    // Teuchos::ParmeterList list. Usually this is read from
    // an xml file. Since the xml is not available in this test
    // the createListFromGetPot reads the required parameters
    // from the datafile and put the in the Teuchos list.
    // You can see the actual implementation at the end of the
    // test.
    //**********************************************************//
    Teuchos::ParameterList solverList;
    createListFromGetPot (solverList, dataFile);

    LinearSolver linearSolver;
    linearSolver.setCommunicator (Comm);
    linearSolver.setParameters ( solverList );
    linearSolver.setPreconditioner ( precPtr );

    //*************************************************************//
    // We are going to solve the laplace equation. The right hand
    // is zero!
    //*************************************************************//
    vectorPtr_Type rhs (new vector_Type ( uSpace -> map() ) );
    *rhs *= 0.0;
    rhs -> globalAssemble();

    //*************************************************************//
    // We impose the boundary conditions, on our system.
    //*************************************************************//

    bcManage ( *systemMatrix, *rhs, *uSpace->mesh(), uSpace->dof(), *BC -> handler(), uFESpace->feBd(), 1.0, 0.0 );

    //*************************************************************//
    // We declare the vector where we want to put the solution.
    // Here we want to solve the linear system Ax=b;
    // we tell the linear solver which A to use (systemMatrix)
    // and which right hand side to use (rhs).
    // Then we solve telling the solver to put the solution in the
    // vector x (solution).
    //*************************************************************//
    vectorPtr_Type solution ( new vector_Type ( uFESpace -> map() ) );


    linearSolver.setOperator (systemMatrix);
    linearSolver.setRightHandSide (rhs);
    linearSolver.solve (solution);

    //*************************************************************//
    // We save the potential field just computed on a file
    // using HDF5.
    //*************************************************************//
    ExporterHDF5< mesh_Type > exporter;
    exporter.setMeshProcId ( meshPart, Comm -> MyPID() );
    exporter.setPostDir (problemFolder);
    exporter.setPrefix ("Potential");
    exporter.addVariable ( ExporterData<mesh_Type>::ScalarField,  "potential", uFESpace, solution, UInt (0) );
    exporter.postProcess (0.0);
    exporter.closeFile();


    //*************************************************************//
    // The sheets are defined as the gradient of the scalar
    // potential just computed.  We create three vector where we
    // store the components of the sheets direction.
    // We computed the gradient using the superconvergent
    // gradient recovery patch of ZZ (check that file for more infos)
    //*************************************************************//
    vectorPtr_Type sx (new vector_Type ( uSpace -> map() ) );
    vectorPtr_Type sy (new vector_Type ( uSpace -> map() ) );
    vectorPtr_Type sz (new vector_Type ( uSpace -> map() ) );

    *sx = GradientRecovery::ZZGradient (uSpace, *solution, 0);
    *sy = GradientRecovery::ZZGradient (uSpace, *solution, 1);
    *sz = GradientRecovery::ZZGradient (uSpace, *solution, 2);



    //*************************************************************//
    // We save the potential field just computed on a file
    // using HDF5.
    //*************************************************************//
    vectorPtr_Type rbFiber ( new vector_Type ( vectorFESpace -> map() ) );
    vectorPtr_Type rbSheet ( new vector_Type ( vectorFESpace -> map() ) );
    UInt d =  rbFiber->epetraVector().MyLength();
    UInt nComponentLocalDof = d / 3;

    for ( int l (0); l < nComponentLocalDof; l++)
	{
		UInt iGID = rbFiber->blockMap().GID (l);
        UInt jGID = rbFiber->blockMap().GID (l + nComponentLocalDof);
        UInt kGID = rbFiber->blockMap().GID (l + 2 * nComponentLocalDof);

//        meshFull->point(iGID).MeshVertex::showMe(true);
		Real x = meshFull->point(iGID).x();
		Real y = meshFull->point(iGID).y();
		Real z = meshFull->point(iGID).z();

		//		std::cout << "x: " << x << ", y: " << y << ", z: " << z << "\n";

		Real t = (*solution)[iGID];

		Real alpha_degrees = 90.0 - 180 * t;
		Real alpha = (alpha_degrees * 3.1415 ) / 180.0;

		Real rl = 1.7 + 0.3 * t;
		Real rs =  0.7 + 0.3 * t;

		Real CosU = z / rl;
		Real SinU = std::sqrt(1.0 - CosU * CosU );

		Real CosV(1.);
		if( x*x + y*y >= 1e-12 )
		{
			Real TgV = y / x;
			Real V = std::atan(TgV);
			CosV = sgn(x) * std::cos(V);
//			CosV = x / rs / SinU;
		}


		Real SinV = sgn(y) * std::sqrt(1.0 - CosV * CosV );

		Real dxdu1 = rs * CosU * CosV;
		Real dxdu2 = rs * CosU * SinV;
		Real dxdu3 = - rl * SinU;
		Real norm1 = std::sqrt( dxdu1 * dxdu1 + dxdu2 * dxdu2 + dxdu3 * dxdu3 );

//		Real dxdv1 = - rs * SinU * SinV;
//		Real dxdv2 = rs * SinU * CosV;
//		Real dxdv3 = 0.0;
//		Real norm2 = std::sqrt( dxdv1 * dxdv1 + dxdv2 * dxdv2 + dxdv3 * dxdv3 );

		Real dxdv1 = - SinV;
		Real dxdv2 = CosV;
		Real dxdv3 = 0.0;
		Real norm2 = 1.0;
//		std::cout << "dxdu1: " << dxdu1 << ", SinV: " << SinV << ", CosU: " << CosU << "\n";
		Real SinA = std::sin(alpha);
		Real CosA = std::cos(alpha);
		if(norm1 >= 1e-13 && norm2 >= 1e-13)
		{
			(*rbFiber)[iGID] = dxdu1 * SinA / norm1 + dxdv1 * CosA / norm2;
			(*rbFiber)[jGID] = dxdu2 * SinA / norm1 + dxdv2 * CosA / norm2;
			(*rbFiber)[kGID] = dxdu3 * SinA / norm1 + dxdv3 * CosA / norm2;
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

    ElectrophysiologyUtility::normalize(*rbSheet);
//
//    //*************************************************************//
//    // Now that we have computed all the desired vector fields
//    // we export them using HDF5 format.
//    // From datafile we can choose the name of the files of the
//    // fiber and sheet vectors.
//    // Also the field are saved on the h5 file with a identification
//    // name (that is the one you can select on paraview). I called
//    // this fiberHDF5Name and sheetHDF5Name. This name is important
//    // to know if we want to import the computed fields in other
//    // simulations
//    //*************************************************************//
    std::string outputFiberFileName = dataFile ("problem/output_fiber_filename", "FiberDirection");
    std::string fiberHDF5Name = dataFile ("problem/hdf5_fiber_name", "fibers");
    std::string outputSheetFileName = dataFile ("problem/output_sheet_filename", "SheetDirection");
    std::string sheetHDF5Name = dataFile ("problem/hdf5_sheet_name", "sheets");

//    std::string outputSheetsFileName = dataFile ("problem/output_sheets_filename", "SheetsDirection");
//    std::string sheetsHDF5Name = dataFile ("problem/hdf5_sheets_name", "sheets");
//
    exportVectorField (Comm, meshPart, vectorFESpace, rbFiber, problemFolder, outputFiberFileName, fiberHDF5Name );
    exportVectorField (Comm, meshPart, vectorFESpace, rbSheet, problemFolder, outputSheetFileName, sheetHDF5Name );
//    exportVectorField (Comm, meshPart, vectorFESpace, rbSheet, problemFolder, outputSheetsFileName, sheetsHDF5Name );
//    exportVectorField (Comm, meshPart, vectorFESpace, projection, problemFolder, "Projection", "projection" );
//
//
//    //*************************************************************//
//    // We test if the norm of the fuber field and the norm of the
//    // sheet field are the same
//    //*************************************************************//
//
//    Real normS = rbSheet-> norm2();
//    Real normF = rbFiber-> norm2();

#ifdef HAVE_MPI
    MPI_Finalize();
#endif

//    Real err = std::abs (normF - normS) / std::abs (normS);
//    std::cout << std::setprecision (20) << "\nError: " <<  err << "\nFiber Norm: " <<  normF << "\n";
//    std::cout << std::setprecision (20) << "Sheet Norm: " <<  normS << "\n";
//    if ( err > 1e-13 )
//    {
//        return EXIT_FAILURE; // Norm of solution did not match
//    }
//    else
//    {
//        return EXIT_SUCCESS;
//    }
}

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

//Export vector to file using HDF5 exporter
void exportVectorField (boost::shared_ptr<Epetra_Comm> comm,
                        meshPtr_Type mesh,
                        fespacePtr_Type fespace,
                        vectorPtr_Type vector,
                        std::string postDir,
                        std::string outputName,
                        std::string hdf5name)
{
    ExporterHDF5< mesh_Type > exporter;
    exporter.setMeshProcId ( mesh, comm -> MyPID() );
    exporter.setPostDir (postDir);
    exporter.setPrefix (outputName);
    exporter.addVariable ( ExporterData<mesh_Type>::VectorField,  hdf5name, fespace, vector, UInt (0) );
    exporter.postProcess (0);
    exporter.closeFile();
}



