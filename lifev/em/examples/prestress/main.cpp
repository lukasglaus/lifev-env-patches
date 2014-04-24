#include <lifev/core/LifeV.hpp>
#include <lifev/electrophysiology/solver/ElectroETAMonodomainSolver.hpp>
#include <lifev/electrophysiology/solver/IonicModels/IonicMinimalModel.hpp>

#include <lifev/structure/solver/StructuralConstitutiveLawData.hpp>

#include <lifev/structure/solver/StructuralConstitutiveLaw.hpp>
#include <lifev/structure/solver/StructuralOperator.hpp>
#include <lifev/em/solver/EMStructuralOperator.hpp>

#include <lifev/em/solver/EMNeoHookeanActivatedMaterial.hpp>

#include <lifev/core/filter/ExporterEnsight.hpp>
#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5.hpp>
#endif
#include <lifev/core/filter/ExporterEmpty.hpp>

#include <lifev/eta/fem/ETFESpace.hpp>
#include <lifev/eta/expression/Integrate.hpp>

#include <lifev/core/interpolation/RBFhtp.hpp>
#include <lifev/core/interpolation/RBFhtpVectorial.hpp>
#include <lifev/core/mesh/MeshLoadingUtility.hpp>
#include <lifev/core/mesh/MeshTransformer.hpp>
#include <lifev/core/interpolation/RBFlocallyRescaledVectorial.hpp>
#include <lifev/core/interpolation/RBFlocallyRescaledScalar.hpp>
#include <lifev/core/interpolation/RBFrescaledVectorial.hpp>
#include <lifev/core/interpolation/RBFrescaledScalar.hpp>
//#include <lifev/core/interpolation/RBFscalar.hpp>
#include <lifev/core/interpolation/RBFvectorial.hpp>

#include <lifev/bc_interface/3D/bc/BCInterface3D.hpp>
#include <sys/stat.h>
#include <fstream>

#include <lifev/em/solver/EMSolver.hpp>
#include <lifev/em/solver/EMActiveStrainSolver.hpp>

using namespace LifeV;


void positionVector (const RegionMesh<LinearTetra>& mesh,  VectorEpetra& disp)
{

    //  std:: cout << "\n  mesh list point size : " <<       mesh.storedPoints()  << ",\t local mesh list point size " << localmesh.storedPoints();
    UInt nLocalDof = disp.epetraVector().MyLength();
    UInt nComponentLocalDof = nLocalDof / 3;

    for (UInt k (0); k < nComponentLocalDof; k++)
    {
        //          std::cout << "\n---------";

        UInt iGID = disp.blockMap().GID (k);
        UInt jGID = disp.blockMap().GID (k + nComponentLocalDof);
        UInt kGID = disp.blockMap().GID (k + 2 * nComponentLocalDof);

        //      UInt iLID = disp.blockMap().LID(k);
        //      UInt jLID = disp.blockMap().LID(k + nComponentLocalDof);
        //      UInt kLID = disp.blockMap().LID(k + 2 * nComponentLocalDof);


        //      disp[iGID] = mesh.point(iGID).x();
        //      std:: cout << "\n  i: " <<  k << ",\t j " << k+nComponentLocalDof << ",\t k " << k+nComponentLocalDof*2;
        //      std:: cout << "\n  x: " << mesh.point(k).x() << ",\t y " <<mesh.point(k+nComponentLocalDof).y()  << ",\t z " << mesh.point(k+nComponentLocalDof*2).z();
        //
        //
        //      std:: cout << "\n  i: " <<  iGID << ",\t j " << jGID << ",\t k " << kGID;
        //      std:: cout << "\n  x: " << mesh.point(iGID).x() << ",\t y " <<mesh.point(iGID).y()  << ",\t z " << mesh.point(iGID).z();
        //      std:: cout << "\n  x: " << localmesh.point(iGID).x() << ",\t y " << localmesh.point(iGID).y()  << ",\t z " << localmesh.point(iGID).z();

        disp[iGID] = mesh.point (iGID).x();
        disp[jGID] = mesh.point (iGID).y();
        disp[kGID] = mesh.point (iGID).z();



    }
    //  disp.setMapType(Repeated);
    //    for ( UInt i = 0; i < mesh.pointList.size(); ++i )
    //    {
    //      std:: cout << "\n--------------";
    //      UInt iGID = disp.blockMap().GID(i);
    //      UInt jGID = disp.blockMap().GID(i + nComponentLocalDof);
    //      UInt kGID = disp.blockMap().GID(i + 2 * nComponentLocalDof);
    ////        for ( UInt j = 0; j < nDimensions; ++j )
    ////        {
    //            int globalId = mesh.pointList[i].id();
    //
    //            std:: cout << "\n i " << i << ", GID " << globalId << ", GID from vector: " << disp.blockMap().GID(i);
    ////            ASSERT ( disp.isGlobalIDPresent ( globalId + dim * j ), "global ID missing" );
    //            std:: cout  << " x : " << disp[ iGID ] << ", " << mesh.pointList[ iGID ].coordinate ( 0 );
    //            std:: cout  << " y : " << disp[ jGID ] << ", " << mesh.pointList[ iGID ].coordinate ( 1 );
    //            std:: cout  << " z : " << disp[ kGID ] << ", " << mesh.pointList[ iGID ].coordinate ( 2 );
    //            disp[ iGID ] = mesh.pointList[ iGID ].coordinate ( 0 );
    //            disp[ jGID ] = mesh.pointList[ iGID ].coordinate ( 1 );
    //            disp[ kGID ] = mesh.pointList[ iGID ].coordinate ( 2 );
    //            std:: cout  << " x : " << disp[ iGID ] << ", " << mesh.pointList[ iGID ].coordinate ( 0 );
    //            std:: cout  << " y : " << disp[ jGID ] << ", " << mesh.pointList[ iGID ].coordinate ( 1 );
    //            std:: cout  << " z : " << disp[ kGID ] << ", " << mesh.pointList[ iGID ].coordinate ( 2 );
    ////        }
    //    }
    //  disp.setMapType(Unique);


}

void changeMesh (const RegionMesh<LinearTetra>& /*mesh*/,  VectorEpetra& positionVector)
{
    Int nLocalDof = positionVector.epetraVector().MyLength();
    Int nComponentLocalDof = nLocalDof / 3;

    for (UInt k (0); k < nComponentLocalDof; k++)
    {
        //UInt iGID = positionVector.blockMap().GID (k);
        //UInt jGID = positionVector.blockMap().GID (k + nComponentLocalDof);
        //UInt kGID = positionVector.blockMap().GID (k + 2 * nComponentLocalDof);

        //      mesh.point(iGID).setX(positionVector[iGID]);
        //      mesh.point(iGID).setY(positionVector[jGID]);
        //      mesh.point(iGID).setZ(positionVector[kGID]);

    }
}


int main (int argc, char** argv)
{

    typedef RegionMesh<LinearTetra>                         mesh_Type;
    typedef boost::function < Real (const Real& /*t*/,
                                    const Real &   x,
                                    const Real &   y,
                                    const Real& /*z*/,
                                    const ID&   /*i*/ ) >   function_Type;

    typedef IonicMinimalModel                                        ionicModel_Type;
    typedef boost::shared_ptr<ionicModel_Type>                       ionicModelPtr_Type;
    typedef EMMonodomainSolver< mesh_Type, ionicModel_Type > monodomainSolver_Type;
    typedef boost::shared_ptr< monodomainSolver_Type >               monodomainSolverPtr_Type;

    typedef VectorEpetra                                             vector_Type;
    typedef boost::shared_ptr<vector_Type>                           vectorPtr_Type;

    typedef EMSolver<mesh_Type, ionicModel_Type>    emSolver_Type;
    typedef boost::shared_ptr<emSolver_Type> emSolverPtr_Type;

#ifdef HAVE_MPI
    MPI_Init ( &argc, &argv );
#endif

    boost::shared_ptr<Epetra_Comm> comm (new Epetra_MpiComm (MPI_COMM_WORLD) );



    //*********************************************//
    // creating output folder
    //*********************************************//
    GetPot commandLine ( argc, argv );
    std::string problemFolder = commandLine.follow ( "Output", 2, "-o", "--output" );
    // Create the problem folder
    if ( problemFolder.compare ("./") )
    {
        problemFolder += "/";

        if ( comm->MyPID() == 0 )
        {
            mkdir ( problemFolder.c_str(), 0777 );
        }
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
    Teuchos::ParameterList monodomainList = * ( Teuchos::getParametersFromXmlFile ( "ParamList.xml" ) );
    if ( comm->MyPID() == 0 )
    {
        std::cout << " Done!" << endl;
    }

    //********************************************//
    // We need the GetPot datafile to setup       //
    //                                            //
    //********************************************//
    GetPot command_line (argc, argv);
    const string data_file_name = command_line.follow ("data", 2, "-f",
                                                       "--file");
    GetPot dataFile (data_file_name);

    //********************************************//
    // We create three solvers to solve with:     //
    // 1) Operator Splitting method               //
    //********************************************//
    if ( comm->MyPID() == 0 )
    {
        std::cout << "Building Monodomain Solvers... ";
    }

    emSolverPtr_Type emSolverPtr ( new emSolver_Type (monodomainList, data_file_name, comm) );

    if ( comm->MyPID() == 0 )
    {
        if (emSolverPtr -> monodomainPtr() -> displacementPtr() )
        {
            std::cout << "\nI've set the displacement ptr in monodomain: constructor";
        }
    }


    if ( comm->MyPID() == 0 )
    {
        std::cout << " solver done... ";

    }
    //********************************************//
    // Setting up the initial condition form      //
    // a given function.                          //
    //********************************************//
    if ( comm->MyPID() == 0 )
    {
        cout << "\nInitializing potential and gating variables:  " ;
    }

    //
    //    function_Type stimulus;
    //      stimulus = &PacingProtocolMM;
    //    emSolverPtr -> monodomainPtr() -> setAppliedCurrentFromFunction(stimulus, 0.0);
    //    if ( comm->MyPID() == 0 )
    //    {
    //        if(emSolverPtr -> monodomainPtr() -> displacementPtr()) std::cout << "\nI've set the displacement ptr in monodomain: applied current";
    //    }
    //    if ( comm->MyPID() == 0 )
    //    {
    //        cout << "Done! \n" ;
    //
    //    }

    //********************************************//
    // Create a fiber direction                   //
    //********************************************//
    if ( comm->MyPID() == 0 )
    {
        cout << "\nSetting fibers:  " ;
    }

    VectorSmall<3> fibers;
    fibers[0] = 0.0;
    fibers[1] = 0.0;
    fibers[2] = 1.0;
    emSolverPtr -> monodomainPtr() -> setupFibers ( fibers );
    emSolverPtr -> activationPtr() -> setFiberPtr ( emSolverPtr -> monodomainPtr() -> fiberPtr() );


    VectorSmall<3> sheets;
    sheets[0] = 0.0;
    sheets[1] = 1.0;
    sheets[2] = 0.0;
    emSolverPtr -> solidPtr() -> activeMaterial() -> setupFiberVector (fibers[0], fibers[1], fibers[2]);
    emSolverPtr -> solidPtr() -> activeMaterial() -> setupSheetVector (sheets[0], sheets[1], sheets[2]);
    emSolverPtr -> exportFibersAndSheetsFields (problemFolder);
    //    if ( comm->MyPID() == 0 )
    //    {
    //        if(emSolverPtr -> monodomainPtr() -> displacementPtr()) std::cout << "\nI've set the displacement ptr in monodomain: fibers";
    //    }
    //
    //    if ( comm->MyPID() == 0 )
    //    {
    //        cout << "Done! \n" ;
    //    }
    //
    //    if ( comm->MyPID() == 0 )
    //    {
    //        if(emSolverPtr -> monodomainPtr() -> displacementPtr()) std::cout << "\nI've set the displacement ptr in monodomain";
    //    }
    //********************************************//
    // Create the global matrix: mass + stiffness //
    //********************************************//
    if ( comm->MyPID() == 0 )
    {
        cout << "\nPointer 1:  " << emSolverPtr -> activationPtr() -> gammafPtr() ;
    }
    emSolverPtr -> setup (monodomainList, data_file_name);
    if ( comm->MyPID() == 0 )
    {
        cout << "\nPointer 2:  " << emSolverPtr -> activationPtr() -> gammafPtr() ;
    }
    emSolverPtr -> setupExporters (problemFolder);
    //********************************************//
    // Activation time                            //
    //********************************************//
    emSolverPtr -> registerActivationTime (0.0, 0.8);
    emSolverPtr -> exportSolution (0.0);

    //********************************************//
    // Solving the system                         //
    //********************************************//

    if ( comm->MyPID() == 0 )
    {
        cout << "\nImporting Time data from parametr list:  " ;
    }
    Real dt = monodomainList.get ("timeStep", 0.1);
    Real TF = monodomainList.get ("endTime", 150.0);
    Int iter = monodomainList.get ("saveStep", 1.0) / dt;
    Int reactionSubiter = monodomainList.get ("reaction_subiteration", 100);
    Int solidIter = monodomainList.get ("emdt", 1.0) / dt;
    bool coupling = monodomainList.get ("coupling", false);
    if ( comm->MyPID() == 0 )
    {
        std::cout << "\ndt: " << dt;
        std::cout << "\nTF: " << TF;
        std::cout << "\nIter: " << iter;
        std::cout << "\nsubiter: " << reactionSubiter;
        std::cout << "\nsolid iterations: " << solidIter;
        std::cout << "\nSolving coupled problem: " << coupling;
    }
//    Int k (0);

    ////    Real timeReac = 0.0;
    ////    Real timeDiff = 0.0;
    ////    Real timeReacDiff = 0.0;
    //
    //    std::string solutionMethod = monodomainList.get ("solutionMethod", "SVI");
    //

    //    EMEvaluate util;
    //    util.evaluate(*(emSolverPtr -> monodomainPtr() -> potentialPtr()),
    //                *(emSolverPtr -> activationSolidPtr()),
    //                *(emSolverPtr -> monodomainPtr() -> feSpacePtr()),
    //                emSolverPtr -> activationSolidFESpace(),
    //                *(emSolverPtr -> fullSolidMesh()) );



    //    for ( Real t = 0.0; t < TF - dt; )
    //    {
    //      //register activation time
    //      k++;
    //      t = t + dt;
    //
    //        if ( comm->MyPID() == 0 )
    //        {
    //            cout << "\n==============";
    //            cout << "\nTime: " << t;
    //            cout << "\n==============";
    //        }
    //
    //
    //        if ( comm->MyPID() == 0 )
    //        {
    //            cout << "\nSet applied current";
    //        }
    //
    //        emSolverPtr -> monodomainPtr() -> setAppliedCurrentFromFunction ( stimulus, t );
    //        emSolverPtr -> solveOneMonodomainStep();
    //
    //        if(coupling == true)
    //        {
    //            if ( comm->MyPID() == 0 )
    //            {
    //                cout << "\n---------------------";
    //                cout << "\nActive Strain Solver ";
    //                cout << "\n---------------------";
    //            }
    //          emSolverPtr -> solveOneActivationStep();
    //            if ( comm->MyPID() == 0 )
    //            {
    //                if(emSolverPtr -> monodomainPtr() -> displacementPtr()) std::cout << "\nI've set the displacement ptr in monodomain: activation";
    //            }
    //
    //
    //          if( k % solidIter == 0 )
    //          {
    //                if ( comm->MyPID() == 0 )
    //                {
    //                    cout << "\n---------------------";
    //                    cout << "\nInterpolating gammaf ";
    //                    cout << "\n---------------------";
    //                }
    //              emSolverPtr -> updateSolid();
    //                if ( comm->MyPID() == 0 )
    //                {
    //                    cout << "\n------------------";
    //                    cout << "\nStructural Solver ";
    //                    cout << "\n------------------";
    //                }

    boost::shared_ptr< Exporter<RegionMesh<LinearTetra> > > expc;
    expc.reset ( new ExporterHDF5<RegionMesh<LinearTetra> > ( dataFile, "configurations" ) );
    expc -> setPostDir ( problemFolder );
    expc ->setMeshProcId ( emSolverPtr -> solidPtr() -> mesh(), comm->MyPID() );
    //    exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::ScalarField, "solid_gammaf", solidaFESpace, solidGammaf, UInt (0) );

    vectorPtr_Type Xk ( new vector_Type ( emSolverPtr -> solidPtr() -> displacement().map() ) );
    vectorPtr_Type Xoriginal ( new vector_Type ( Xk-> map() ) );
    vectorPtr_Type Xkp1 ( new vector_Type ( Xk-> map() ) );
    vectorPtr_Type res ( new vector_Type ( Xk-> map() ) );
    expc->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, "Xk", emSolverPtr -> activationPtr() -> FESpacePtrSolid(), Xk, UInt (0) );
    expc->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, "Xor", emSolverPtr -> activationPtr() -> FESpacePtrSolid(), Xoriginal, UInt (0) );
    expc->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, "Xkp1", emSolverPtr -> activationPtr() -> FESpacePtrSolid(), Xkp1, UInt (0) );
    expc->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, "disp", emSolverPtr -> activationPtr() -> FESpacePtrSolid(), emSolverPtr -> solidPtr() -> displacementPtr(), UInt (0) );


//    UInt ntotalDof = emSolverPtr -> solidPtr() -> dispFESpace().dof().numTotalDof();
    UInt nLocalDof = emSolverPtr -> solidPtr() -> displacement().epetraVector().MyLength();
    UInt dim = nLocalDof / 3;

    positionVector ( * (emSolverPtr -> fullSolidMesh() ),  *Xoriginal);
    emSolverPtr -> solidPtr() -> mesh() -> meshTransformer().vertices (*Xkp1, dim);
    *Xk = *Xoriginal;
    //    positionVector( *(emSolverPtr -> fullSolidMesh() ),  *Xoriginal);
    //    positionVector( *(emSolverPtr -> fullSolidMesh() ),  *Xoriginal, *( emSolverPtr -> solidPtr() -> mesh() ));


    //emSolverPtr -> solidPtr() -> mesh() -> meshTransformer().vertices(*Xk, 3); //moveMesh( -emSolverPtr -> solidPtr() -> displacement(), 3);
    // emSolverPtr -> solidPtr() -> mesh() -> meshTransformer().vertices(*Xoriginal, 3); //moveMesh( -emSolverPtr -> solidPtr() -> displacement(), 3);

    expc -> postProcess (0);

    //
    emSolverPtr ->solidBCPtr() -> handler() -> showMe();



    //                emSolverPtr -> solveSolid();
    //
    //                *Xkp1 *= 0.0;
    //                *Xkp1 += *Xoriginal;
    //                *Xkp1 -= emSolverPtr -> solidPtr() -> displacement();
    //
    //                *res = *Xkp1;
    //                *res -= *Xk;
    //
    ////                *res = emSolverPtr -> solidPtr() -> displacement();
    //Real normres = 2.0;
    //Real normdisp = 2.0;
    // res -> norm2();

    Real tol = 1e-6;
    int iterk (0);
    //
    //
    //                  expc -> postProcess(1.0);
    //                  expc -> closeFile();
    //                emSolverPtr -> solidPtr() -> mesh() -> meshTransformer().moveMesh( -emSolverPtr -> solidPtr() -> displacement(), 3);
    //                emSolverPtr -> fullSolidMesh() -> meshTransformer().moveMesh( -( emSolverPtr -> solidPtr() -> displacement()), 3);
    std::ofstream output ("res.txt");




    for ( Real t (0.0) ;  t < TF ; )
    {
        t += dt;
        emSolverPtr ->solidDataPtr() -> dataTime() -> setTime (t);

        if ( comm->MyPID() == 0 )
        {
            std::cout << "\n/////////////////////////////////////////////////////";
            std::cout << "\n Time " << t;
            std::cout << "\n/////////////////////////////////////////////////////";
                  }

                      if ( comm->MyPID() == 0 )
                  {
          output << "\nnormres" << ", " << "normdisp" << ", time: " << t ;
                  }

                      Real normres = 2.0;
                      Real normdisp = 2.0;
                      iterk = 0;
                      while (normres > tol)
                  {
                      iterk++;
                      emSolverPtr -> solidPtr() -> displacement() = *Xoriginal;
                      emSolverPtr -> solidPtr() -> displacement() -= *Xk;

                      emSolverPtr -> solveSolid();
                      emSolverPtr -> exportSolution (t + 0.01 * iterk);

                      *res = *Xk;
                      *res += emSolverPtr -> solidPtr() -> displacement();
                      *res -= *Xoriginal;

                      normres =  res -> normInf();
                      normdisp =  emSolverPtr -> solidPtr() -> displacement().normInf();
                      //
                      if ( comm->MyPID() == 0 )
                  {
          std::cout << "\n RESIDUAL : " << normres << ", Norm disp: " << normdisp ;
                      output << "\n" << normres << ", " << normdisp;
                  }

                      *Xk *= 0.0;
                      *Xk += *Xoriginal;
                      *Xk -= emSolverPtr -> solidPtr() -> displacement();


                      emSolverPtr -> solidPtr() -> mesh() -> meshTransformer().moveMesh ( *Xk, dim, 1);
                      *Xkp1 = emSolverPtr -> solidPtr() -> displacement();
                      //emSolverPtr -> solidPtr() -> mesh() -> meshTransformer().moveMesh( -emSolverPtr -> solidPtr() -> displacement(), dim, 0);

                      expc -> postProcess (t + 0.01 * iterk);

                      if (iterk > 15)
                  {
                      break;
                  }
                  }
                  }

                      output.close();
                      //      //

                      //                if ( comm->MyPID() == 0 )
                      //                {
                      //                    cout << "\n---------------------------";
                      //                    cout << "\nInterpolating displacement ";
                      //                    cout << "\n---------------------------";
                      //                }
                      //                emSolverPtr -> updateMonodomain();
                      //          }
                      //
                      //
                      //        }
                      //      emSolverPtr -> registerActivationTime(t, 0.8);
                      //      matrixPtr_Type broydenMatrix( emSolverPtr -> solidPtr() -> );
                      //        if( k % iter == 0 )
                      //        {
                      //            if ( comm->MyPID() == 0 )
                      //            {
                      //                cout << "\nExporting solutions";
                      //            }


                      //emSolverPtr -> solidPtr() -> mesh() -> meshTransformer().moveMesh( -( emSolverPtr -> solidPtr() -> displacement()), 3);
                      //emSolverPtr -> exportSolution(2.0);

                      //        }
                      //
                      //    }
                      //
                      //    if ( comm->MyPID() == 0 )
                      //    {
                      //        cout << "\nExporting Activation Time";
                      //    }
                      iterk++;
                      emSolverPtr -> exportSolution (iterk);
                      emSolverPtr -> closeExporters();
                      expc -> closeFile();




          std::cout << "\n\nThank you for using EMSolver.\nI hope to meet you again soon!\n All the best for your simulation :P\n  " ;
                          //   }

#ifdef HAVE_MPI
                          MPI_Finalize();
#endif
                          return 0;
                      }
