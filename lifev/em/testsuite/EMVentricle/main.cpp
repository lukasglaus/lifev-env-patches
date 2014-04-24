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
    Real returnValue = 0;

    if ( t <= 2.0 )
    {
        if ( std::abs ( x ) <= 0.5 &&
                std::abs ( y ) <= 0.5 &&
                std::abs ( z ) <= 0.2 )
        {
            returnValue = 10;
        }
    }

    return returnValue;
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
        cout << "Done! \n" ;
    }

    //********************************************//
    // Create a fiber direction                   //
    //********************************************//
    if ( comm->MyPID() == 0 )
    {
        cout << "\nSetting fibers:  " ;
    }

    emSolverPtr -> setFibersAndSheets (monodomainList);

    emSolverPtr -> exportFibersAndSheetsFields ( problemFolder);
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

    emSolverPtr -> setup (monodomainList, data_file_name);

    {
        cout << "\nPointer 2:  " << emSolverPtr -> activationPtr() -> gammafPtr() ;
    }

    emSolverPtr -> setupExporters (problemFolder);
    //********************************************//
    // Activation time                            //
    //********************************************//
    emSolverPtr -> registerActivationTime (0.0, 0.8);
    emSolverPtr -> exportSolution ( 0.0);

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
    //    if ( comm->MyPID() == 0 )
    //      {
    //        std::cout << "\ndt: " << dt;
    //        std::cout << "\nTF: " << TF;
    //        std::cout << "\nIter: " << iter;
    //        std::cout << "\nsubiter: " << reactionSubiter;
    //        std::cout << "\nsolid iterations: " << solidIter;
    //        std::cout << "\nSolving coupled problem: " << coupling;
    //      }
    Int k (0);

    //    Real timeReac = 0.0;
    //    Real timeDiff = 0.0;
    //    Real timeReacDiff = 0.0;

    std::string solutionMethod = monodomainList.get ("solutionMethod", "splitting");


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
            cout << "\n------------------";
            cout << "\nMonodomain Solver ";
            cout << "\n------------------";
        }
        if ( comm->MyPID() == 0 )
        {
            cout << "\nSet applied current";
        }

        emSolverPtr -> monodomainPtr() -> setAppliedCurrentFromFunction ( stimulus, t );
        if (solutionMethod == "splitting")
        {
            if ( comm->MyPID() == 0 )
            {
                cout << "\nSplitting: Solving reactions";
            }
            for (int j (0); j < reactionSubiter; j++)
            {
                emSolverPtr -> monodomainPtr() -> solveOneReactionStepFE (reactionSubiter);
            }
            //solve diffusion step
            if ( comm->MyPID() == 0 )
            {
                cout << "\nSplitting: : Solving diffusion";
            }
            emSolverPtr -> solveOneDiffusionStep();
        }
        if (solutionMethod == "HLS")
        {
            * (emSolverPtr -> monodomainPtr() -> appliedCurrentPtr() ) *= 10.0;
            if ( comm->MyPID() == 0 )
            {
                cout << "\nHLS: Solving gating variables ";
            }
            emSolverPtr -> monodomainPtr() -> solveOneStepGatingVariablesFE();
            if ( comm->MyPID() == 0 )
            {
                cout << "\nHLS: Solving potential ";
            }
            emSolverPtr -> monodomainPtr() ->solveOneICIStep ( * (emSolverPtr -> activationPtr() -> massMatrixPtr() ) );
        }
        if (solutionMethod == "ICI")
        {
            if ( comm->MyPID() == 0 )
            {
                cout << "\nICI: Solving gating variables ";
            }
            emSolverPtr -> monodomainPtr() -> solveOneStepGatingVariablesFE();
            if ( comm->MyPID() == 0 )
            {
                cout << "\nICI: Solving potential ";
            }
            emSolverPtr -> monodomainPtr() ->solveOneICIStep();
        }
        if (solutionMethod == "SVI")
        {
            if ( comm->MyPID() == 0 )
            {
                if (emSolverPtr -> monodomainPtr() -> displacementPtr() )
                {
                    std::cout << "\nI've set the displacement ptr in monodomain: exporters";
                }
            }
            if ( comm->MyPID() == 0 )
            {
                cout << "\nSVI: Solving gating variables ";
            }
            emSolverPtr -> monodomainPtr() -> solveOneStepGatingVariablesFE();
            if ( comm->MyPID() == 0 )
            {
                cout << "\nSVI: Solving potential ";
            }
            emSolverPtr -> monodomainPtr() -> solveOneSVIStep();
        }

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
        if ( k % iter == 0 )
        {
            if ( comm->MyPID() == 0 )
            {
                cout << "\nExporting solutions";
            }
            emSolverPtr -> registerActivationTime (t, 0.8);
            emSolverPtr -> exportSolution (t);
        }

    }

    if ( comm->MyPID() == 0 )
    {
        cout << "\nExporting Activation Time";
    }

    emSolverPtr -> exportActivationTime (problemFolder);

    emSolverPtr -> closeExporters();
    if ( comm->MyPID() == 0 )
    {
        std::cout << "\nExporting fibers: " << std::endl;
    }

    //********************************************//
    // Saving Fiber direction to file             //
    //********************************************//
    emSolverPtr -> monodomainPtr() -> exportFiberDirection (problemFolder);

    std::cout << "\n\nThank you for using EMSolver.\nI hope to meet you again soon!\n All the best for your simulation :P\n  " ;

#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return 0;
}

