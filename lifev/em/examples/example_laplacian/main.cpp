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
    @brief Tutorial introducing the expression assembly

    @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @date 28-06-2012

    In this first tutorial, we assemble the matrix
    associated to a scalar laplacian problem. The basics
    of the ETA module are explained and are pushed further
    in the next tutorials.

    ETA stands Expression Template Assembly, in reference
    to the metaprogramming technique used.

 */

// ---------------------------------------------------------------
// We include here the MPI headers for the parallel computations.
// The specific "pragma" instructions are used to avoid warning
// coming from the MPI library, that are not useful to us.
// ---------------------------------------------------------------

#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"


// ---------------------------------------------------------------
// We include then the required headers from LifeV. First of all,
// the definition file and mesh related files. We also include
// the MatrixEpetra since this is the kind of object that we want
// to assemble.
// ---------------------------------------------------------------

#include <lifev/core/LifeV.hpp>

#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/core/mesh/RegionMesh3DStructured.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>

#include <lifev/core/array/MatrixEpetra.hpp>

// ---------------------------------------------------------------
// In order to use the ETA framework, a special version of the
// FESpace structure must be used. It is called ETFESpace and
// has basically the same role as the FESpace.
// ---------------------------------------------------------------

#include <lifev/eta/fem/ETFESpace.hpp>


// ---------------------------------------------------------------
// The most important file to include is the Integrate.hpp file
// which contains all the definitions required to perform the
// different integrations.
// ---------------------------------------------------------------

#include <lifev/eta/expression/Integrate.hpp>


// ---------------------------------------------------------------
// Finally, we include shared pointer from boost since we use
// them explicitly in this tutorial.
// ---------------------------------------------------------------

#include <boost/shared_ptr.hpp>

#include <lifev/electrophysiology/util/HeartUtility.hpp>
#include <lifev/core/mesh/MeshUtility.hpp>
#include <lifev/core/mesh/MeshLoadingUtility.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include "Teuchos_XMLParameterListHelpers.hpp"

#include <lifev/bc_interface/3D/bc/BCInterface3D.hpp>
#include <lifev/core/algorithm/LinearSolver.hpp>
#include <lifev/core/algorithm/Preconditioner.hpp>
#include <lifev/core/fem/BCManage.hpp>

#include <lifev/core/filter/ExporterEnsight.hpp>
#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5.hpp>
#endif
#include <lifev/core/filter/ExporterEmpty.hpp>
// ---------------------------------------------------------------
// As usual, we work in the LifeV namespace. For clarity, we also
// make two typedefs for the mesh type and matrix type.
// ---------------------------------------------------------------

using namespace LifeV;

typedef RegionMesh<LinearTetra> mesh_Type;
typedef boost::shared_ptr< mesh_Type > meshPtr_Type;
typedef MatrixEpetra<Real> matrix_Type;
typedef boost::shared_ptr< matrix_Type > matrixPtr_Type;
typedef VectorEpetra vector_Type;
typedef boost::shared_ptr< vector_Type > vectorPtr_Type;
typedef BCHandler                                          bc_Type;
typedef boost::shared_ptr< bc_Type >                       bcPtr_Type;
typedef  EMStructuralOperator< RegionMesh<LinearTetra> >      physicalSolver_Type;
typedef BCInterface3D< bc_Type, physicalSolver_Type >              bcInterface_Type;
typedef boost::shared_ptr< bcInterface_Type >              bcInterfacePtr_Type;


// ---------------------------------------------------------------
// We start the programm by the definition of the communicator
// (as usual) depending on whether MPI is available or not. We
// also define a boolean to allow only one process to display
// messages.
// ---------------------------------------------------------------

Real bcZero (const Real& /*t*/, const Real&  /*X*/, const Real& /*Y*/, const Real& /*Z*/, const ID& /*i*/)
{
    return  0.;
}

Real bcOne (const Real& /*t*/, const Real&  /*X*/, const Real& /*Y*/, const Real& /*Z*/, const ID& /*i*/)
{
    return  1.0;
}

int main ( int argc, char** argv )
{

#ifdef HAVE_MPI
    MPI_Init (&argc, &argv);
    boost::shared_ptr<Epetra_Comm> Comm (new Epetra_MpiComm (MPI_COMM_WORLD) );
#else
    boost::shared_ptr<Epetra_Comm> Comm (new Epetra_SerialComm);
#endif

    const bool verbose (Comm->MyPID() == 0);


    // ---------------------------------------------------------------
    // The next step is to build the mesh. We use here a structured
    // cartesian mesh over the square domain (-1,1)x(-1,1)x(-1,1).
    // The mesh is the partitioned for the parallel computations and
    // the original mesh is deleted.
    // ---------------------------------------------------------------

    if (verbose)
    {
        std::cout << " -- Building and partitioning the mesh ... " << std::flush;
    }
    //********************************************//
    // Import parameters from an xml list. Use    //
    // Teuchos to create a list from a given file //
    // in the execution directory.                //
    //********************************************//

    if ( Comm->MyPID() == 0 )
    {
        std::cout << "Importing parameters list...";
    }
    Teuchos::ParameterList parameterList = * ( Teuchos::getParametersFromXmlFile ( "ParamList.xml" ) );
    if ( Comm->MyPID() == 0 )
    {
        std::cout << " Done!" << endl;
    }


    //********************************************//
    // In the parameter list we need to specify   //
    // the mesh name and the mesh path.           //
    //********************************************//
    if ( Comm->MyPID() == 0 )
    {
        std::cout << "Reading Mesh Name and Path...\n";
    }

    std::string meshName = parameterList.get ("mesh_name", "lid16.mesh");
    std::string meshPath = parameterList.get ("mesh_path", "./");

    meshPtr_Type meshPart (new mesh_Type ( Comm ) );
    MeshUtility::fillWithMesh (meshPart, meshName, meshPath);

    if (verbose)
    {
        std::cout << " done ! " << std::endl;
    }


    // ---------------------------------------------------------------
    // We define now the ETFESpace that we need for the assembly.
    // Remark that we use a shared pointer because other structures
    // will require this ETFESpace to be alive. We can also observe
    // that the ETFESpace has more template parameters than the
    // classical FESpace (this is the main difference). The 3
    // indicates that the problem is in 3D while the 1 indicate that
    // the unknown is scalar.
    //
    // After having constructed the ETFESpace, we display the number
    // of degrees of freedom of the problem.
    // ---------------------------------------------------------------

    if (verbose)
    {
        std::cout << " -- Building ETFESpaces ... " << std::flush;
    }

    boost::shared_ptr<ETFESpace< mesh_Type, MapEpetra, 3, 1 > > uSpace
    ( new ETFESpace< mesh_Type, MapEpetra, 3, 1 > (meshPart, &feTetraP1, Comm) );
    boost::shared_ptr<FESpace< mesh_Type, MapEpetra > > uFESpace
    ( new FESpace< mesh_Type, MapEpetra > (meshPart, "P1", 1, Comm) );

    if (verbose)
    {
        std::cout << " done ! " << std::endl;
    }
    if (verbose)
    {
        std::cout << " ---> Dofs: " << uSpace->dof().numTotalDof() << std::endl;
    }


    // ---------------------------------------------------------------
    // The matrix is then defined using the map of the FE space.
    // ---------------------------------------------------------------

    if (verbose)
    {
        std::cout << " -- Defining the matrix ... " << std::flush;
    }

    boost::shared_ptr<matrix_Type> systemMatrix (new matrix_Type ( uSpace->map() ) );

    *systemMatrix *= 0.0;

    if (verbose)
    {
        std::cout << " done! " << std::endl;
    }

    *systemMatrix *= 0.0;

    if (verbose)
    {
        std::cout << " done! " << std::endl;
    }


    // ---------------------------------------------------------------
    // We start now the assembly of the matrix.
    // ---------------------------------------------------------------

    if (verbose)
    {
        std::cout << " -- Assembling the Laplace matrix ... " << std::flush;
    }


    // ---------------------------------------------------------------
    // To use the ETA framework, it is mandatory to use a special
    // namespace, called ExpressionAssembly. This namespace is useful
    // to avoid collisions with keywords used for the assembly. A
    // special scope is opened to keep only that part of the code
    // in the ExpressionAssembly namespace.
    // ---------------------------------------------------------------

    {
        using namespace ExpressionAssembly;

        // ---------------------------------------------------------------
        // We can now proceed with assembly. The next instruction
        // assembles the laplace operator.
        //
        // The first argument of the integrate function indicates that the
        // integration is done on the elements of the mesh located in the
        // ETFESpace defined earlier.
        //
        // The second argument is simply the quadrature rule to be used.
        //
        // The third argument is the finite element space of the test
        // functions.
        //
        // The fourth argument is the finite element space of the trial
        // functions (those used to represent the solution).
        //
        // The last argument is the expression to be integrated, i.e.
        // that represents the weak formulation of the problem. The
        // keyword phi_i stands for a generic test function and phi_j
        // a generic trial function. The function grad applied to them
        // indicates that the gradient is considered and the dot function
        // indicates a dot product between the two gradients. The
        // expression to be integrated is then the dot product between
        // the gradient of the test function and the gradient of the trial
        // function. This corresponds to the left hand side of the weak
        // formulation of the Laplace problem.
        //
        // Finally, the operator >> indicates that the result of the
        // integration must be added to the systemMatrix.
        // ---------------------------------------------------------------

        integrate (  elements (uSpace->mesh() ),
                     quadRuleTetra4pt,
                     uSpace,
                     uSpace,
                     dot ( grad (phi_i) , grad (phi_j) )
                  )
                >> systemMatrix;
    }

    if (verbose)
    {
        std::cout << " done! " << std::endl;
    }




    // ---------------------------------------------------------------
    // As we are already done with the assembly of the matrix, we
    // finalize it to be able to work on it, e.g. to solve a linear
    // system.
    // ---------------------------------------------------------------

    if (verbose)
    {
        std::cout << " -- Closing the matrix ... " << std::flush;
    }

    systemMatrix->globalAssemble();

    if (verbose)
    {
        std::cout << " done ! " << std::endl;
    }

    if ( Comm->MyPID() == 0 )
    {
        std::cout << "Creating Boundary Conditions...\n";
    }
    //-----------------------
    //  BOundary Conditions
    GetPot command_line (argc, argv);
    const string data_file_name = command_line.follow ("data", 2, "-f", "--file");
    GetPot dataFile (data_file_name);
    bcInterfacePtr_Type                     BC ( new bcInterface_Type() );
    BC->createHandler();
    BC->fillHandler ( data_file_name, "problem" );
    BC -> handler() -> bcUpdate ( *uFESpace->mesh(), uFESpace->feBd(), uFESpace->dof() );

    if ( Comm->MyPID() == 0 )
    {
        std::cout << "Creating Preconditioner...\n";
    }

    //Preconditioner
    typedef LifeV::Preconditioner             basePrec_Type;
    typedef boost::shared_ptr<basePrec_Type>  basePrecPtr_Type;
    typedef LifeV::PreconditionerIfpack       prec_Type;
    typedef boost::shared_ptr<prec_Type>      precPtr_Type;

    prec_Type* precRawPtr;
    basePrecPtr_Type precPtr;
    precRawPtr = new prec_Type;
    precRawPtr->setDataFromGetPot ( dataFile, "prec" );
    precPtr.reset ( precRawPtr );

    if ( Comm->MyPID() == 0 )
    {
        std::cout << "Creating linear solver...\n";
    }

    //linear solver
    Teuchos::RCP< Teuchos::ParameterList > belosList3 = Teuchos::rcp ( new Teuchos::ParameterList );
    belosList3 = Teuchos::getParametersFromXmlFile ( "ParamList.xml" );

    LinearSolver linearSolver;
    linearSolver.setCommunicator (Comm);
    //    .setCommunicator ( Comm );
    linearSolver.setParameters ( *belosList3 );
    linearSolver.setPreconditioner ( precPtr );


    if ( Comm->MyPID() == 0 )
    {
        std::cout << "Creating rhs...\n";
    }
    //Create right hand side
    vectorPtr_Type rhs (new vector_Type ( uSpace -> map() ) );
    *rhs *= 0.0;
    rhs -> globalAssemble();


    if ( Comm->MyPID() == 0 )
    {
        std::cout << "Applying BC...\n";
    }
    bcManage ( *systemMatrix, *rhs, *uSpace->mesh(), uSpace->dof(), *BC -> handler(), uFESpace->feBd(), 1.0, 0.0 );




    linearSolver.setOperator (systemMatrix);
    vectorPtr_Type solution ( new vector_Type ( uFESpace -> map() ) );

    if ( Comm->MyPID() == 0 )
    {
        std::cout << "Exporter...\n";
    }

    ExporterHDF5< RegionMesh <LinearTetra> > exporter;
    exporter.setMeshProcId ( meshPart, Comm->MyPID() );
    exporter.setPrefix ("rescalingGammaf");
    exporter.addVariable ( ExporterData<mesh_Type>::ScalarField,  "rescalingGammaf", uFESpace,
                           solution, UInt (0) );

    //   exporter.postProcess ( 0 );

    if ( Comm->MyPID() == 0 )
    {
        std::cout << "Solve System...\n";
    }
    linearSolver.setRightHandSide (rhs);
    linearSolver.solve (solution);
    //  exporter.postProcess ( 1 );


    if ( Comm->MyPID() == 0 )
    {
        std::cout << "Transforming...\n";
    }
    Real p0 = 3.9981022451448112e-01;
    Real p1 = -1.7249368443147499e+00;
    Real p2 = 6.4695359113991486e+00;
    Real p3 = -1.9192450904013128e+01;
    Real p4 = 4.1297847160139170e+01;
    Real p5 = -5.9665612469298118e+01;
    Real p6 = 5.4040368339225651e+01;
    Real p7 = -2.7524476381745629e+01;
    Real p8 = 5.9999955063690678e+00;

    vector_Type tmp (solution -> map() );
    exporter.postProcess ( 0 );
    tmp *= 0.0;
    tmp += p0;
    tmp += ( p1 * *solution);
    tmp += ( p2 * *solution * *solution);
    tmp += ( p3 * *solution * *solution * *solution);
    tmp += ( p4 * *solution * *solution * *solution * *solution);
    tmp += ( p5 * *solution * *solution * *solution * *solution * *solution);
    tmp += ( p6 * *solution * *solution * *solution * *solution * *solution * *solution);
    tmp += ( p7 * *solution * *solution * *solution * *solution * *solution * *solution * *solution);
    tmp += ( p8 * *solution * *solution * *solution * *solution * *solution * *solution * *solution * *solution);

    *solution = tmp;
    exporter.postProcess ( 1 );
    exporter.closeFile();

    vectorPtr_Type importedSolution ( new vector_Type ( solution -> map() ) );

    std::string filename = parameterList.get ("filename", "rescalingGammaf");
    std::string fieldname = parameterList.get ("fieldname", "rescalingGammaf");
    ElectrophysiologyUtility::importScalarField (importedSolution, filename, fieldname, meshPart);

    ExporterHDF5< RegionMesh <LinearTetra> > exporter2;
    exporter2.setMeshProcId ( meshPart, Comm->MyPID() );
    exporter2.setPrefix ("rescalingGammaf_refined1");
    exporter2.addVariable ( ExporterData<mesh_Type>::ScalarField,  "rescalingGammaf", uFESpace,
                            importedSolution, UInt (0) );
    exporter2.postProcess ( 0 );
    exporter2.closeFile();
    // ---------------------------------------------------------------
    // We finalize the MPI session if MPI was used
    // ---------------------------------------------------------------

#ifdef HAVE_MPI
    MPI_Finalize();
#endif




}


