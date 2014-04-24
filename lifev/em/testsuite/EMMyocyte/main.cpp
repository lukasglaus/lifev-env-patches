#include <lifev/core/LifeV.hpp>
#include <lifev/electrophysiology/solver/ElectroETAMonodomainSolver.hpp>
#include <lifev/electrophysiology/solver/IonicModels/IonicMinimalModel.hpp>

#include <lifev/structure/solver/StructuralConstitutiveLawData.hpp>

#include <lifev/structure/solver/StructuralConstitutiveLaw.hpp>
#include <lifev/structure/solver/StructuralOperator.hpp>
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


Real PacingProtocolMM ( const Real& t, const Real& x, const Real& y, const Real& z, const ID&   /*id*/)
{
    Real pacingSite_X = 0.0045;
    Real pacingSite_Y = 0.0025;
    Real pacingSite_Z = 0.0;
    Real stimulusRadius = 0.0005;
    Real stimulusValue = 10.;
    //Real stimulusPeriod = 2.;

    Real returnValue;

    Real r = std::sqrt ( (x - pacingSite_X ) * (x - pacingSite_X ) + (y - pacingSite_Y ) * (y - pacingSite_Y ) + (z - pacingSite_Z ) * (z - pacingSite_Z ) );
    if ( r <  stimulusRadius  && t < 0.1 )
    {
        returnValue = stimulusValue;
    }
    else
    {
        returnValue = 0.;
    }

    return returnValue;
}

Real initialSphereOnCell (const Real& /*t*/, const Real&  X, const Real& Y, const Real& Z, const ID& /*i*/)
{

    double r = std::sqrt (pow (X - 42, 2) + pow (Y - 46, 2) + pow (Z - 7, 2) );
    double auxexp = 1.0 - 1.0 / (1.0 + exp (-50.0 * (r - 8) ) );

    return 0.1 + 3.5 * auxexp;
}

Real initialSphereOnCylinder (const Real& /*t*/, const Real&  X, const Real& Y, const Real& Z, const ID& /*i*/)
{

    double r = std::sqrt (pow (X, 2) + pow (Y - 60, 2) + pow (Z - 12, 2) );
    double auxexp = 1.0 - 1.0 / (1.0 + exp (-50.0 * (r - 12) ) );

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


    function_Type stimulus;
    stimulus = &PacingProtocolMM;
    emSolverPtr -> monodomainPtr() -> setAppliedCurrentFromFunction (stimulus, 0.0);
    if ( comm->MyPID() == 0 )
    {
        if (emSolverPtr -> monodomainPtr() -> displacementPtr() )
        {
            std::cout << "\nI've set the displacement ptr in monodomain: applied current";
        }
    }
    if ( comm->MyPID() == 0 )
    {
        /*
                std::cout << " Splitting solver done... ";
            }

            monodomain -> setInitialConditions();

            function_Type f = &initialSphereOnCell;
            monodomain -> setPotentialFromFunction(f);

        */
        cout << "Done! \n" ;

    }

    //********************************************//
    // Create a fiber direction                   //
    //********************************************//
    if ( comm->MyPID() == 0 )
    {
        cout << "\nSetting fibers:  " ;
    }

    /*  VectorSmall<3> fibers;
      fibers[0]=0.0;
      fibers[1]=0.0;
      fibers[2]=1.0;
      emSolverPtr -> monodomainPtr() -> setupFibers( fibers );
      emSolverPtr -> activationPtr() -> setFiberPtr( emSolverPtr -> monodomainPtr() -> fiberPtr() );


      VectorSmall<3> sheets;
      sheets[0]=0.0;
      sheets[1]=1.0;
      sheets[2]=0.0;
      emSolverPtr -> solidPtr() -> material() -> setupFiberVector(fibers[0],fibers[1],fibers[2]);
      emSolverPtr -> solidPtr() -> material() -> setupSheetVector(sheets[0],sheets[1],sheets[2]);*/

    emSolverPtr -> setFibersAndSheets (monodomainList);
    emSolverPtr -> exportFibersAndSheetsFields (problemFolder);
    if ( comm->MyPID() == 0 )
    {
        if (emSolverPtr -> monodomainPtr() -> displacementPtr() )
        {
            std::cout << "\nI've set the displacement ptr in monodomain: fibers";
        }
    }

    if ( comm->MyPID() == 0 )
    {
        cout << "Done! \n" ;
    }

    if ( comm->MyPID() == 0 )
    {
        if (emSolverPtr -> monodomainPtr() -> displacementPtr() )
        {
            std::cout << "\nI've set the displacement ptr in monodomain";
        }
    }
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

    emSolverPtr -> activationPtr() -> setActivationType ( "TransverselyIsotropic" );

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
        /*
            {
                std::cout << "\nsetup structural operator" << std::endl;
            }

            //! 1. Constructor of the structuralSolver
             StructuralOperator< RegionMesh<LinearTetra> > solid;
             solid.setup (dataStructure,
                          dFESpace,
                          dETFESpace,
                          solidBC -> handler(),
                          comm);

             //if ( comm->MyPID() == pid )
             if ( comm->MyPID() == 0 )
               {
                 std::cout << "\ninitial guess" << std::endl;
             }

             solid.setDataFromGetPot (dataFile);

                //===========================================================
                //===========================================================
                //              FIBERS
                //===========================================================
                //===========================================================
             // for(int pid(0); pid < 4 ; pid ++){
               // if ( comm->MyPID() == pid )
             if ( comm->MyPID() == 0 )
             {
                 std::cout << "\nreading fibers ... " << std::endl;
             }
             //   }


             vectorPtr_Type solidFibers( new vector_Type( dFESpace -> map() ) );


             std::vector<Real> fvec(3, 0.0);
             fvec.at(0)  = parameterList.get ("fiber_X", 0.0);
             fvec.at(1)  = parameterList.get ("fiber_Y", 1.0);
             fvec.at(2)  = parameterList.get ("fiber_Z", 0.0);
             HeartUtility::setupFibers(*solidFibers, fvec);

             MPI_Barrier(MPI_COMM_WORLD);


             if ( comm->MyPID() == 0 )
             {
                 std::cout << "\nread fibers" << std::endl;
             }

             solid.material() -> setFiberVector( *solidFibers );

        //     monodomain -> setupFibers();

             vectorPtr_Type gammaf( new vector_Type( ( monodomain -> globalSolution().at(0) ) -> map() ) );
             vectorPtr_Type gammas( new vector_Type( ( monodomain -> globalSolution().at(0) ) -> map() ) );
             vectorPtr_Type gamman( new vector_Type( ( monodomain -> globalSolution().at(0) ) -> map() ) );
             vectorPtr_Type solidGammaf;
             vectorPtr_Type emDisp;
             solidFESpacePtr_Type electroFiberFESpace;
             solidETFESpacePtr_Type electrodETFESpace;

            solidGammaf = gammaf;
            monodomain -> setFiberPtr( solidFibers );
            emDisp = solid.displacementPtr();
            electroFiberFESpace = dFESpace;
            electrodETFESpace = dETFESpace;



             monodomain -> exportFiberDirection(problemFolder);
         */    //********************************************//
        // Create the global matrix: mass + stiffness in ELECTROPHYSIOLOGY //
        //********************************************//
        /*    if ( comm->MyPID() == 0 )
            {
                cout << "\nSetup operators:   = " << monodomain -> timeStep() << "\n" ;
            }

            monodomain -> setDisplacementPtr( emDisp );
            monodomain -> setupMassMatrix();
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
            *gammaf *= 0.0;
            *solidGammaf =0;

            solid.material() -> setGammaf( *solidGammaf );

            vectorPtr_Type solidGammas( new vector_Type( solidGammaf -> map() ) );
            vectorPtr_Type solidGamman( new vector_Type( solidGammaf -> map() ) );

         *solidGammas = 1.0;
         *solidGammas /= (1.0 + *solidGammaf);
         EpetraSqrt(*solidGammas);
         *solidGammas -= 1.0;
         solid.material() -> setGamman(*solidGammas);
         *solidGamman = *solidGammas;
         solid.material() -> setGammas(*solidGamman);

         gamman = solidGamman;
         gammas = solidGammas;

            //==================================================================//
            //==================================================================//
            //                 SETUP INTERPOLATION                             //
            //==================================================================//
            //==================================================================//

            if ( comm->MyPID() == 0 )
            {
                std::cout << "\nbuild solid system" << std::endl;
            }

            solid.buildSystem(1.0);
            vectorPtr_Type rhs (new vector_Type (solid.displacement(), Unique) );
            vectorPtr_Type initialDisplacement (new vector_Type (solid.displacement(), Unique) );
            solid.initialize ( initialDisplacement );


            MPI_Barrier (MPI_COMM_WORLD);

            if ( comm->MyPID() == 0 )
            {
                std::cout << "\nsetup solid exporter" << std::endl;
            }

             boost::shared_ptr< Exporter<RegionMesh<LinearTetra> > > exporter;
             exporter.reset ( new ExporterHDF5<RegionMesh<LinearTetra> > ( dataFile, parameterList.get ("StructureOutputFile", "StructureOutput") ) );

             exporter -> setPostDir ( problemFolder );
             exporter->setMeshProcId ( localSolidMesh, comm->MyPID() );

             vectorPtr_Type solidDisp ( new vector_Type (solid.displacement(), exporter->mapType() ) );
             exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, "displacement", dFESpace, solidDisp, UInt (0) );
             exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::ScalarField, "solid_gammaf", solidaFESpace, solidGammaf, UInt (0) );





             //================================================================//
             //================================================================//
             //                    SETUP COUPLING SOLVER                       //
             //                                                                //
             //================================================================//
             //================================================================//
             ExporterHDF5< RegionMesh <LinearTetra> > expGammaf;
        expGammaf.setMeshProcId(monodomain -> localMeshPtr(), comm->MyPID());
        expGammaf.setPrefix(parameterList.get ("ActivationOutputFile", "ActivationOutput"));
        expGammaf.setPostDir ( problemFolder );

           matrixPtr_Type mass(new matrix_Type( monodomain -> massMatrixPtr() -> map() ) ) ;

        boost::shared_ptr<FLRelationship> fl (new FLRelationship);

        boost::shared_ptr<HeavisideFct> H (new HeavisideFct);
        */
    {
        std::cout << "\ndt: " << dt;
        std::cout << "\nTF: " << TF;
        std::cout << "\nIter: " << iter;
        std::cout << "\nsubiter: " << reactionSubiter;
        std::cout << "\nsolid iterations: " << solidIter;
        std::cout << "\nSolving coupled problem: " << coupling;
    }
    Int k (0);


    //    Real timeReac = 0.0;
    //    Real timeDiff = 0.0;
    //    Real timeReacDiff = 0.0;

    std::string solutionMethod = monodomainList.get ("solutionMethod", "SVI");


    //    EMEvaluate util;
    //    util.evaluate(*(emSolverPtr -> monodomainPtr() -> potentialPtr()),
    //                *(emSolverPtr -> activationSolidPtr()),
    //                *(emSolverPtr -> monodomainPtr() -> feSpacePtr()),
    //                emSolverPtr -> activationSolidFESpace(),
    //                *(emSolverPtr -> fullSolidMesh()) );
    emSolverPtr -> preloadRamp (0.1);


    for ( Real t = 0.0; t < TF - dt; )
    {
        //register activation time
        k++;
        t = t + dt;

        if ( comm->MyPID() == 0 )
        {
            cout << "\n==============";
            cout << "\nTime: " << t;
            cout << "\n==============";
        }


        if ( comm->MyPID() == 0 )
        {
            cout << "\nSet applied current";
        }

        emSolverPtr -> monodomainPtr() -> setAppliedCurrentFromFunction ( stimulus, t );
        emSolverPtr -> solveOneMonodomainStep();

        if (coupling == true)
        {
            if ( comm->MyPID() == 0 )
            {
                cout << "\n---------------------";
                cout << "\nActive Strain Solver ";
                cout << "\n---------------------";
            }
            emSolverPtr -> solveOneActivationStep();
            if ( comm->MyPID() == 0 )
            {
                if (emSolverPtr -> monodomainPtr() -> displacementPtr() )
                {
                    std::cout << "\nI've set the displacement ptr in monodomain: activation";
                }
            }


            if ( k % solidIter == 0 )
            {
                if ( comm->MyPID() == 0 )
                {
                    cout << "\n---------------------";
                    cout << "\nInterpolating gammaf ";
                    cout << "\n---------------------";
                }
                emSolverPtr -> updateSolid();
                if ( comm->MyPID() == 0 )
                {
                    cout << "\n------------------";
                    cout << "\nStructural Solver ";
                    cout << "\n------------------";
                }
                emSolverPtr -> solveSolid();
                if ( comm->MyPID() == 0 )
                {
                    cout << "\n---------------------------";
                    cout << "\nInterpolating displacement ";
                    cout << "\n---------------------------";
                }
                emSolverPtr -> updateMonodomain();
            }


        }
        emSolverPtr -> registerActivationTime (t, 0.8);
        //      matrixPtr_Type broydenMatrix( emSolverPtr -> solidPtr() -> );
        if ( k % iter == 0 )
        {
            if ( comm->MyPID() == 0 )
            {
                cout << "\nExporting solutions";
            }

            emSolverPtr -> exportSolution (t);
        }

    }

    if ( comm->MyPID() == 0 )
    {
        cout << "\nExporting Activation Time";
    }

    emSolverPtr -> exportActivationTime (problemFolder);
    emSolverPtr -> closeExporters();
    /*    if ( comm->MyPID() == 0 )

          {
            std::cout << "\nIt does!!!!" << std::endl;
          }




        //Initial condition for gammaf
        *gammaf = -0.015;

        vectorPtr_Type tmpRhsActivation( new vector_Type ( rhsActivation -> map(), Repeated ) );
        solidFESpacePtr_Type emDispFESpace ( new solidFESpace_Type ( monodomain -> localMeshPtr(), "P1", 3, comm) );
        expGammaf.addVariable(ExporterData<mesh_Type>::ScalarField, "gammaf",
                  monodomain -> feSpacePtr(), gammaf, UInt(0));
        expGammaf.addVariable(ExporterData<mesh_Type>::VectorField, "interpolated displacement",
                  emDispFESpace, emDisp, UInt(0));
        expGammaf.addVariable(ExporterData<mesh_Type>::ScalarField, "rhs",
                  monodomain -> feSpacePtr(), rhsActivation, UInt(0));

        expGammaf.postProcess(0.0);

        //===========================================================
        //===========================================================
        //              Initializing solid
        //===========================================================
        //===========================================================

        exporter->postProcess ( 0 );

        vectorPtr_Type emDisp0(new vector_Type( emDisp -> map() ) );
        *emDisp0 = *emDisp;

        //===========================================================
        //===========================================================
        //              TIME LOOP
        //===========================================================
        //===========================================================
        Real emdt = parameterList.get("emdt",1.0);
        int iter((emdt / monodomain -> timeStep()));
        int k(0);
        Real saveStep = parameterList.get("save_step",1.0);
        int saveIter((saveStep / monodomain -> timeStep()));
        int subiter = parameterList.get("subiter",100);

        for( Real t(0.0); t< monodomain -> endTime(); )
          {
        t = t + monodomain -> timeStep();
        k++;

        LifeChrono timer;

        for(int j(0); j<subiter; j++) monodomain -> solveOneReactionStepFE(subiter);

        timer.stop();
    */
    std::cout << "\nExporting fibers: " << std::endl;

    //********************************************//
    // Saving Fiber direction to file             //
    //********************************************//
    emSolverPtr -> monodomainPtr() -> exportFiberDirection (problemFolder);

    /*
        *tmpRhsActivation *= 0;
        if ( comm->MyPID() == 0 )
          {
            std::cout << "\nASSEMBLING ACTIVATION EQUATION!\n" << std::endl;
          }
        {
          using namespace ExpressionAssembly;



          BOOST_AUTO_TPL(Ca ,   value( aETFESpace, *( monodomain -> globalSolution().at(0)  ) ) );
          BOOST_AUTO_TPL(Gammaf,     value( aETFESpace, *gammaf )  );
          //Calcium-dependent activation
          BOOST_AUTO_TPL(activationEquation, value(-0.5)*Ca + value(-2.5)*Gammaf  );

          integrate ( elements ( monodomain -> localMeshPtr() ),
                  monodomain -> feSpacePtr() -> qr() ,
                  monodomain -> ETFESpacePtr(),
                  activationEquation  * phi_i
                  ) >> tmpRhsActivation;
        }

        *rhsActivation *= 0;
        *rhsActivation = ( *(mass) * ( *gammaf ) );
        *rhsActivation += ( ( monodomain -> timeStep() * *tmpRhsActivation ) );

        linearSolver.setRightHandSide(rhsActivation);

        if ( comm->MyPID() == 0 )
          {
            std::cout << "\nSOLVING ACTIVATION EQUATION!\n" << std::endl;
          }

        linearSolver.solve(gammaf);

        if ( k % iter == 0)
          {
            solidGammaf = gammaf;
            solid.material() -> setGammaf( *solidGammaf );
            *solidGammas = 1.0;
            *solidGammas /= (1.0 + *solidGammaf);
            EpetraSqrt(*solidGammas);
            *solidGammas -= 1.0;
            solid.material() -> setGamman(*solidGammas);
            *solidGamman = *solidGammas;
            solid.material() -> setGammas(*solidGamman);

            if ( comm->MyPID() == 0 )
              {
            std::cout << "\nSOLVING STATIC MECHANICS!\n" << std::endl;
              }

            if ( comm->MyPID() == 0 )
              {
            std::cout << "\n*****************************************************";
            std::cout << "\nWE ARE AT TIME: "<< t;
            std::cout << "\n*****************************************************";
              }

            solid.iterate ( solidBC -> handler() );
            *solidDisp = solid.displacement();

            if ( comm->MyPID() == 0 )
              {
            std::cout << "\nREASSEMBLING STIFFNESS MATRIX FOR TOW WAY COUPLING!\n" << std::endl;
              }

            monodomain -> setupStiffnessMatrix();
            monodomain -> setupGlobalMatrix();
          }

        //cout << "\n\n save every " << saveIter << "iteration\n";
        if ( k % saveIter == 0)
          {
            monodomain -> registerActivationTime(*activationTimeVector, t, 1.0);
            monodomain -> exportSolution(expElectro, t);
            expGammaf.postProcess(t);
            exporter->postProcess ( t );
          }
          }*/

    std::cout << "\n\nThank you for using EMSolver.\nI hope to meet you again soon!\n All the best for your simulation :P\n  " ;
    //   }

#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return 0;
}
