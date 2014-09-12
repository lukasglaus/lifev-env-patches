#include <lifev/core/LifeV.hpp>


#include <lifev/structure/solver/StructuralOperator.hpp>
#include <lifev/structure/solver/MixedStructuralOperator.hpp>
#include <lifev/structure/solver/StructuralConstitutiveLawData.hpp>
#include <lifev/core/filter/ExporterHDF5.hpp>

#include <lifev/bc_interface/3D/bc/BCInterface3D.hpp>
#include <lifev/core/fem/BCHandler.hpp>
#include <lifev/core/fem/BCManage.hpp>

#include <lifev/eta/fem/ETFESpace.hpp>
#include <lifev/eta/expression/Integrate.hpp>

#include <lifev/core/algorithm/LinearSolver.hpp>
#include <lifev/core/algorithm/Preconditioner.hpp>
#include <lifev/core/algorithm/PreconditionerML.hpp>
#include <lifev/core/algorithm/PreconditionerIfpack.hpp>

#include <lifev/core/array/VectorEpetraStructured.hpp>
#include <lifev/core/array/VectorEpetraStructuredView.hpp>
#include <lifev/core/array/MatrixEpetraStructured.hpp>
#include <lifev/core/array/MatrixEpetraStructuredView.hpp>


#include <lifev/core/array/VectorBlockMonolithicEpetra.hpp>

#include <lifev/core/mesh/MeshLoadingUtility.hpp>




using namespace LifeV;


Real bcZero (const Real& /*t*/, const Real&  /*X*/, const Real& /*Y*/, const Real& /*Z*/, const ID& /*i*/)
{
    return  0.;
}


Real bcOne (const Real& /*t*/, const Real&  /*X*/, const Real& /*Y*/, const Real& /*Z*/, const ID& /*i*/)
{
    return  1.;
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

    typedef VectorEpetra                vector_Type;
    typedef boost::shared_ptr<vector_Type> vectorPtr_Type;

    typedef MatrixEpetra<Real> matrix_Type;
    typedef boost::shared_ptr<matrix_Type> matrixPtr_Type;

    typedef BCHandler                                          bc_Type;
    typedef boost::shared_ptr< bc_Type >                       bcPtr_Type;
    typedef StructuralOperator< RegionMesh<LinearTetra> >      physicalSolver_Type;
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

	meshPtr_Type mesh ( new mesh_Type ( comm ) );
	meshPtr_Type fullMesh ( new mesh_Type ( comm ) );
	MeshUtility::loadMesh( mesh, fullMesh, meshName, meshPath, false, "P1");
//	MeshUtility::loadMesh (mesh, fullMesh, meshName, meshPath);
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

    MixedStructuralOperator<mesh_Type> solver(mesh, comm);
    solver.setup(dataFile);
//
//
//    //! #################################################################################
//    //! BOUNDARY CONDITIONS
//    //! #################################################################################
//    std::vector <ID> compx (1), compy (1), compz (1), compxy (2), compxz (2), compyz (2);
//    compx[0] = 0;
//    compy[0] = 1, compz[0] = 2;
//    compxy[0] = 0;
//    compxy[1] = 1;
//    compxz[0] = 0;
//    compxz[1] = 2;
//    compyz[0] = 1;
//    compyz[1] = 2;
//
//      boost::shared_ptr<BCHandler> BCh ( new BCHandler() );
//
      BCFunctionBase zero (bcZero);
      BCFunctionBase load (bcOne);


      solver.M_bchandlerPtr->addBC ("Top",     600,  Essential, Full,      zero,    3);
      solver.M_bchandlerPtr->addBC ("Bottom",  900,  Natural,   Normal, load);


      solver.iterate();

    std::cout << "\n Test Passed! Congrats!\n\n";

#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return 0;
}
