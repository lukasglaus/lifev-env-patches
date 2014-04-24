#include <lifev/core/LifeV.hpp>
#include <lifev/electrophysiology/solver/ElectroETAMonodomainSolver.hpp>
#include <lifev/electrophysiology/solver/IonicModels/IonicMinimalModel.hpp>

#include <lifev/structure/solver/StructuralConstitutiveLawData.hpp>

#include <lifev/structure/solver/StructuralConstitutiveLaw.hpp>
#include <lifev/structure/solver/StructuralOperator.hpp>
#include <lifev/structure/solver/NeoHookeanActivatedMaterial.hpp>

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

void EpetraPow (VectorEpetra& vector, const Real p)
{
    Int size = vector.epetraVector().MyLength();

    for (int j (0); j < size; j++)
    {
        int gid = vector.blockMap().GID (j);
        vector[gid] = std::pow (vector[gid], p);
    }
}
void EpetraSqrt (VectorEpetra& vector)
{
    Int size = vector.epetraVector().MyLength();

    for (int j (0); j < size; j++)
    {
        int gid = vector.blockMap().GID (j);
        vector[gid] = std::sqrt (vector[gid]);
    }
}

Real bcZero (const Real& /*t*/, const Real& /*X*/, const Real& /*Y*/,
             const Real& /*Z*/, const ID& /*i*/)
{
    return 0.;
}
Real d0 (const Real& /*t*/, const Real& /*X*/, const Real& /*Y*/,
         const Real& /*Z*/, const ID& /*i*/)
{
    return 0.;
}

Real initialVlid (const Real& /*t*/, const Real& X, const Real& /*Y*/,
                  const Real& Z, const ID& /*i*/)
{
    if (X > 0.95)
    {
        return 1.0;
    }
    else
    {
        return 0.;
    }
}

Real initialV0left (const Real& /*t*/, const Real& X, const Real& /*Y*/,
                    const Real& /*Z*/, const ID& /*i*/)
{
    if (X == 0)
    {
        return 1.0;
    }
    else
    {
        return 0.;
    }
}

Real fiberRotationRing (const Real& /*t*/, const Real& X, const Real& Y,
                        const Real& /*Z*/, const ID& i)
{
    Real R = std::sqrt (X * X + Y * Y);
    //Real teta = std::atan( Y / X );
    Real fz = 0.0;
    Real fx = Y / R;
    Real fy = -X / R;
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
            return fx * std::cos (theta) + fx * sx * sx * (1.0 - std::cos (theta) )
                   + sx * sy * fy * (1.0 - std::cos (theta) );
            break;
        case 1:
            return sx * sy * fy * (1.0 - std::cos (theta) ) + fy * std::cos (theta)
                   + fy * sy * sy * (1.0 - std::cos (theta) );
            break;
        case 2:
            return sx * fy * std::sin (theta) - sy * fx * std::sin (theta);
            break;
        default:
            ERROR_MSG ("This entry is not allowed: ud_functions.hpp")
            ;
            return 0.;
            break;
    }

}

static Real f0fun (const Real&, const Real& x, const Real&, const Real&,
                   const ID& comp)
{
    Real p = 3.14159265358979;
    Real alpha = -2.0 * p / 3. * x + p / 3.;

    Real compx = 0.0;
    Real compy = std::cos (alpha);
    Real compz = std::sin (alpha);

    //    compx = 0.0;
    //    compy = 0.0;
    //    compz = 1.0;
    if (comp == 0)
    {
        return compx;
    }
    else if (comp == 1)
    {
        return compy;
    }
    else
    {
        return compz;
    }
}

Real rescalingGamma (const Real&, const Real& x, const Real&, const Real&,
                     const ID& /*comp*/)
{
    Real r = (1 - x);
    Real p = 3.14159265358979;
    Real alpha = -2.0 * p / 3. * x + p / 3.;
    Real circumferentialShortening = (0.05 - 0.2) * r + 0.2;
    Real max_gamma = circumferentialShortening / std::cos (alpha);
    return max_gamma;
}

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

void createPositionXVector (const RegionMesh<LinearTetra>& mesh,
                            VectorEpetra& positionVector)
{
    Int nLocalDof = positionVector.epetraVector().MyLength();
    Int nComponentLocalDof = nLocalDof / 3;
    for (int k (0); k < nComponentLocalDof; k++)
    {
        UInt iGID = positionVector.blockMap().GID (k);
        positionVector[iGID] = mesh.point (iGID).x();
    }

}

void computeX (VectorEpetra& positionVector)
{
    Int nLocalDof = positionVector.epetraVector().MyLength();
    Int nComponentLocalDof = nLocalDof / 3;
    for (int k (0); k < nComponentLocalDof; k++)
    {
        UInt iGID = positionVector.blockMap().GID (k);
        UInt jGID = positionVector.blockMap().GID (k + nComponentLocalDof);
        UInt kGID = positionVector.blockMap().GID (k + 2 * nComponentLocalDof);

        positionVector[jGID] = 0.0;
        positionVector[kGID] = 0.0;
    }

}

template<typename space> Real ComputeVolume (
    const boost::shared_ptr<RegionMesh<LinearTetra> > localMesh,
    VectorEpetra positionVector, const VectorEpetra& disp,
    const boost::shared_ptr <
    ETFESpace<RegionMesh<LinearTetra>, MapEpetra, 3, 1> > ETFESpace,
    const boost::shared_ptr<space> dETFESpace, int bdFlag,
    boost::shared_ptr<Epetra_Comm> comm)
{

    Real fluidVolume;

    MatrixSmall<3, 3> Id;
    Id (0, 0) = 1.;
    Id (0, 1) = 0., Id (0, 2) = 0.;
    Id (1, 0) = 0.;
    Id (1, 1) = 1., Id (1, 2) = 0.;
    Id (2, 0) = 0.;
    Id (2, 1) = 0., Id (2, 2) = 1.;
    VectorSmall<3> E1;
    E1 (0) = 1.;
    E1 (1) = 0.;
    E1 (2) = 0.;

    Int nLocalDof = positionVector.epetraVector().MyLength();
    for (int k (0); k < nLocalDof; k++)
    {
        UInt iGID = positionVector.blockMap().GID (k);

        positionVector[iGID] += disp[iGID];
    }

    boost::shared_ptr<VectorEpetra> intergral (
        new VectorEpetra (positionVector.map() ) );

    {
        using namespace ExpressionAssembly;

        BOOST_AUTO_TPL (I, value (Id) );
        BOOST_AUTO_TPL (vE1, value (E1) );
        BOOST_AUTO_TPL (Grad_u, grad (dETFESpace, disp, 0) );
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

Real ComputeCubicVolume (const RegionMesh<LinearTetra> fullMesh,
                         VectorEpetra positionVector, const VectorEpetra& disp, int bdFlag,
                         boost::shared_ptr<Epetra_Comm> comm)
{
    positionVector += disp;
    Real xMin (0.);
    Real xMax (0.);
    Real yMin (0.);
    Real yMax (0.);
    Real zMin (0.);
    Real zMax (0.);
    Int nLocalDof = positionVector.epetraVector().MyLength();
    Int nComponentLocalDof = nLocalDof / 3;

    for (int k (0); k < nComponentLocalDof; k++)
    {
        UInt iGID = positionVector.blockMap().GID (k);
        UInt jGID = positionVector.blockMap().GID (k + nComponentLocalDof);
        UInt kGID = positionVector.blockMap().GID (k + 2 * nComponentLocalDof);
        if (fullMesh.point (iGID).markerID() == bdFlag)
        {
            if (positionVector[iGID] > xMax)
            {
                xMax = positionVector[iGID];
            }
            if (positionVector[iGID] < xMin)
            {
                xMin = positionVector[iGID];
            }
            if (positionVector[jGID] > yMax)
            {
                yMax = positionVector[jGID];
            }
            if (positionVector[jGID] < yMin)
            {
                yMin = positionVector[jGID];
            }
            if (positionVector[kGID] > zMax)
            {
                zMax = positionVector[kGID];
            }
            if (positionVector[kGID] < zMin)
            {
                zMin = positionVector[kGID];
            }

        }
    }

    //  Real p = 3.14159265358979;
    //  Real Volume;
    //  if(xDiameter > yDiameter) Volume = p * xDiameter * xDiameter / 4.0 * (zMax - zMin) / 3.0;
    //  else Volume = p * yDiameter * xDiameter / 4.0 * (zMax - zMin) / 3.0;
    int numProc = comm->NumProc();
    Real xMinGlobal (0.);
    Real xMaxGlobal (0.);
    Real yMinGlobal (0.);
    Real yMaxGlobal (0.);
    Real zMinGlobal (0.);
    Real zMaxGlobal (0.);

    if (numProc == 1)
    {
        Real xDiameter = (xMax - xMin);
        Real yDiameter = (yMax - yMin);
        Real Volume = xDiameter * yDiameter * (zMax - zMin) / 2;
        return Volume;
    }
    else
    {
        comm->MaxAll (&xMax, &xMaxGlobal, 1);
        comm->MinAll (&xMin, &xMinGlobal, 1);
        comm->MaxAll (&yMax, &yMaxGlobal, 1);
        comm->MinAll (&yMin, &yMinGlobal, 1);
        comm->MaxAll (&zMax, &zMaxGlobal, 1);
        comm->MinAll (&zMin, &zMinGlobal, 1);

        Real xDiameter = (xMaxGlobal - xMinGlobal);
        Real yDiameter = (yMaxGlobal - yMinGlobal);
        Real Volume = xDiameter * yDiameter * (zMaxGlobal - zMinGlobal) / 2;
        return Volume;
    }

}

Real evaluatePressure (Real Volume, Real dV, Real pn, Real dp_temporal,
                       Real dV_temporal, Real Cp, Int phase = 0, Real t = 0.0, Real dt = 1.0)
{
    Real pressure;
    Real R = 150 * 1333.22; //mmHg ms / ml
    Real C = 0.9 / 1333.22; //ml / mmHg

    switch (phase)
    {
        case 0:
            pressure = pn + dV_temporal / Cp;
            break;
        case 1:
            pressure = pn - (pn * dt / C / R + dV_temporal / C);
            break;
        case 2:
            pressure = pn + dV_temporal / Cp;
            break;
        case 3:
            R = 100 * 1333;
            C = 0.2 / 1333;
            pressure = pn + (pn * dt / C / R - dV_temporal / C);
            //          if(t > 250)pressure = 12000 + 0.35 * (t-250.0) * (t-250.0);
            //          else pressure = pn;
            break;
        default:
            std::cout
                    << "\nThis case is not yet available\n I'm not gonna change the pressure.\n";
            pressure = pn;
            break;
    }
    return pressure;
}

int main (int argc, char** argv)
{

    typedef RegionMesh<LinearTetra> mesh_Type;
    typedef boost::shared_ptr<mesh_Type> meshPtr_Type;
    typedef boost::function <
    Real (const Real& /*t*/, const Real & x, const Real & y,
          const Real& /*z*/, const ID& /*i*/) > function_Type;
    typedef IonicMinimalModel ionicModel_Type;
    typedef boost::shared_ptr<ionicModel_Type> ionicModelPtr_Type;

    typedef ElectroETAMonodomainSolver<mesh_Type, ionicModel_Type> monodomainSolver_Type;
    typedef boost::shared_ptr<monodomainSolver_Type> monodomainSolverPtr_Type;
    typedef VectorEpetra vector_Type;
    typedef boost::shared_ptr<vector_Type> vectorPtr_Type;

    typedef MatrixEpetra<Real> matrix_Type;
    typedef boost::shared_ptr<matrix_Type> matrixPtr_Type;

    typedef BCHandler bc_Type;
    typedef boost::shared_ptr<bc_Type> bcPtr_Type;
    typedef EMStructuralOperator<RegionMesh<LinearTetra> > physicalSolver_Type;
    typedef BCInterface3D<bc_Type, physicalSolver_Type> bcInterface_Type;
    typedef boost::shared_ptr<bcInterface_Type> bcInterfacePtr_Type;
    typedef EMSolver<mesh_Type, ionicModel_Type> emSolver_Type;
    typedef boost::shared_ptr<emSolver_Type> emSolverPtr_Type;

#ifdef HAVE_MPI
    MPI_Init ( &argc, &argv );
#endif

    boost::shared_ptr<Epetra_Comm> comm (new Epetra_MpiComm (MPI_COMM_WORLD) );
    //*********************************************//
    // creating output folder
    //*********************************************//
    GetPot commandLine (argc, argv);
    std::string problemFolder = commandLine.follow ("Output", 2, "-o",
                                                    "--output");
    // Create the problem folder
    if (problemFolder.compare ("./") )
    {
        problemFolder += "/";

        if (comm->MyPID() == 0)
        {
            mkdir (problemFolder.c_str(), 0777);
        }
    }

    //===========================================================
    //===========================================================
    //              ELECTROPHYSIOLOGY
    //===========================================================
    //===========================================================

    if (comm->MyPID() == 0)
    {
        cout << "% using MPI" << endl;
    }

    //********************************************//
    // Import parameters from an xml list. Use    //
    // Teuchos to create a list from a given file //
    // in the execution directory.                //
    //********************************************//

    if (comm->MyPID() == 0)
    {
        std::cout << "Importing parameters list...";
    }
    Teuchos::ParameterList parameterList = * (Teuchos::getParametersFromXmlFile (
                                                  "ParamList.xml") );
    if (comm->MyPID() == 0)
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


    emSolverPtr_Type emSolverPtr ( new emSolver_Type (parameterList, data_file_name, comm) );

    emSolverPtr -> setPotentialOnBoundary (1.0, 10);

    emSolverPtr -> setFibersAndSheets (parameterList);

    emSolverPtr -> setup (parameterList, data_file_name, comm);

    emSolverPtr -> exportFibersAndSheetsFields ( comm, problemFolder);

    emSolverPtr -> setupExporters (comm, problemFolder);

    emSolverPtr -> registerActivationTime (0.0, 0.2);

    emSolverPtr -> preloadRamp (0.5);

    int iter (0);
    int reactionSubiter = parameterList.get ("subiter", 1);
    for (Real time = 0.0; time < emSolverPtr -> monodomainPtr() -> endTime(); )
    {
        iter++;
        Real electroStep = iter * (emSolverPtr -> monodomainPtr() -> timeStep() );

        emSolverPtr -> updateMonodomain();
        while ( time < electroStep )
        {
            //solve reaction step subiter times
            for (int j (0); j < reactionSubiter; j++)
            {
                emSolverPtr -> M_monodomainPtr -> solveOneReactionStepFE (reactionSubiter);
            }
            //solve diffusion step
            emSolverPtr -> solveOneDiffusionStep();

            //active strain
            emSolverPtr -> solveOneActivationStep();
        }

        emSolverPtr -> updateSolid();
        emSolverPtr -> solveSolid();
        emSolverPtr -> exportSolution (comm, time);
    }



    emSolverPtr -> exportActivationTime (comm, problemFolder);

    emSolverPtr -> closeExporters (comm);
    if (comm->MyPID() == 0)
    {
        for (int i (0); i < emSolverPtr -> lvFlags().size(); i++)
        {
            std::cout << "\nFlags " << i << ": " << emSolverPtr -> lvFlags().at (i) ;
        }
    }

    //    //********************************************//
    //    // In the parameter list we need to specify   //
    //    // the mesh name and the mesh path.           //
    //    //********************************************//

    //    Real Cp = parameterList.get ("Cp", 1e-5);
    //    Real pvtol = parameterList.get ("pvtol", 1e-3);

    //    bool load4restart = parameterList.get("load4restart", false);
    ////    ionicModel -> initialize( monodomain -> globalSolution() );
    //
    //    if(load4restart)
    //    {
    //      std::string V0filename = parameterList.get("V0filename", "V0");
    //      std::string V0fieldname = parameterList.get("V0fieldname", "V0");
    //      HeartUtility::importScalarField(monodomain -> globalSolution().at(0),V0filename,V0fieldname,monodomain -> localMeshPtr() );
    //      std::string V1filename = parameterList.get("V1filename", "V1");
    //      std::string V1fieldname = parameterList.get("V1fieldname", "V1");
    //      HeartUtility::importScalarField(monodomain -> globalSolution().at(1),V1filename,V1fieldname,monodomain -> localMeshPtr() );
    //      std::string V2filename = parameterList.get("V2filename", "V2");
    //      std::string V2fieldname = parameterList.get("V2fieldname", "V2");
    //      HeartUtility::importScalarField(monodomain -> globalSolution().at(2),V2filename,V2fieldname,monodomain -> localMeshPtr() );
    //      std::string V3filename = parameterList.get("V3filename", "V3");
    //      std::string V3fieldname = parameterList.get("V3fieldname", "V3");
    //      HeartUtility::importScalarField(monodomain -> globalSolution().at(3),V3filename,V3fieldname,monodomain -> localMeshPtr() );
    //    }
    //    else{
    //

    //
    //      //function_Type Vlid = &initialVlid;
    //
    //    }
    //
    //     if(load4restart)
    //     {
    //      std::string Dfilename = parameterList.get("Gfilename", "G");
    //      std::string Dfieldname = parameterList.get("Gfieldname", "G");
    //      HeartUtility::importVectorField( solid.displacementPtr(), Dfilename, Dfieldname,  localSolidMesh );
    //     }
    //      //===========================================================
    //      //===========================================================
    //      //              FIBERS
    //      //===========================================================
    //      //===========================================================
    //
    //     monodomain -> exportFiberDirection(problemFolder);
    //     boost::shared_ptr< Exporter<RegionMesh<LinearTetra> > > exporterSheets;
    //     exporterSheets.reset ( new ExporterHDF5<RegionMesh<LinearTetra> > ( dataFile, "SheetsDirection" ) );
    //
    //       //      exporter->setPostDir ( "./" );
    //             exporterSheets -> setPostDir ( problemFolder );
    //       exporterSheets->setMeshProcId ( localSolidMesh, comm->MyPID() );
    //       exporterSheets->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, "sheets", dFESpace, solidSheets, UInt (0) );
    //       exporterSheets->postProcess(0);
    //       exporterSheets->closeFile();
    //       //********************************************//
    //     // Create the global matrix: mass + stiffness in ELECTROPHYSIOLOGY //
    //     //********************************************//
    //     if ( comm->MyPID() == 0 )
    //     {
    //         cout << "\nSetup operators:  dt = " << monodomain -> timeStep() << "\n" ;
    //     }
    //
    //     monodomain -> setDisplacementPtr( emDisp );
    //     monodomain -> setupLumpedMassMatrix();
    //     monodomain -> setupStiffnessMatrix();
    //     monodomain -> setupGlobalMatrix();
    //
    //     if ( comm->MyPID() == 0 )
    //     {
    //         cout << "Done! \n" ;
    //     }


    //
    //
    //
    //
    //      //================================================================//
    //      //================================================================//
    //      //                  SETUP COUPLING SOLVER                       //
    //      //                                                              //
    //      //================================================================//
    //      //================================================================//
    //      ExporterHDF5< RegionMesh <LinearTetra> > expGammaf;
    //  expGammaf.setMeshProcId(monodomain -> localMeshPtr(), comm->MyPID());
    //  expGammaf.setPrefix(parameterList.get ("ActivationOutputFile", "ActivationOutput"));
    //  expGammaf.setPostDir ( problemFolder );
    //
    ////    expGammaf.addVariable(ExporterData<mesh_Type>::ScalarField, "gammaf",
    ////            monodomain -> feSpacePtr(), gammaf, UInt(0));
    ////    expGammaf.postProcess(0.0);
    ////    Real min =  0.2;
    ////    Real max =  0.85;
    ////
    ////    Real beta = -0.3;
    ////

    //
    //
    //
    //      if ( comm->MyPID() == 0 )
    //      {
    //          std::cout << "\nSolve system" << std::endl;
    //      }

    //
    //
    //
    //      vectorPtr_Type tmpRhsActivation( new vector_Type ( rhsActivation -> map(), Repeated ) );
    //    solidFESpacePtr_Type emDispFESpace ( new solidFESpace_Type ( monodomain -> localMeshPtr(), "P1", 3, comm) );
    //      expGammaf.addVariable(ExporterData<mesh_Type>::ScalarField, "gammaf",
    //              monodomain -> feSpacePtr(), gammaf, UInt(0));
    //      expGammaf.addVariable(ExporterData<mesh_Type>::VectorField, "interpolated displacement",
    //              emDispFESpace, emDisp, UInt(0));
    //      expGammaf.addVariable(ExporterData<mesh_Type>::ScalarField, "rhs",
    //              monodomain -> feSpacePtr(), rhsActivation, UInt(0));
    //
    //
    //      expGammaf.postProcess(0.0);
    //
    //    //================================================================//
    //    //================================================================//
    //    //                    SETUP gamma rescaling                         //
    //    //                                                                  //
    //    //================================================================//
    //    //================================================================//
    //      vectorPtr_Type rescaling(new vector_Type( solidGammaf -> map() ) );
    //      function_Type resc = &rescalingGamma;
    //    solidaFESpace -> interpolate ( static_cast< FESpace< RegionMesh<LinearTetra>, MapEpetra >::function_Type > ( resc ), *rescaling , 0);
    //  //  vectorPtr_Type emDisp0(new vector_Type( emDisp -> map() ) );


    //
    //  if(usingDifferentMeshes)
    //  {
    //      if(solid.displacementPtr() -> normInf() > 0)
    //      {
    //          if ( comm->MyPID() == 0 )
    //          {
    //              std::cout << "\nINTERPOLATING FROM COARSE TO FINE!\n" << std::endl;
    //          }
    //
    //          C2F -> updateRhs( solid.displacementPtr() );
    //          C2F -> interpolate();
    //          C2F -> solution( emDisp );
    //      }
    //      else *emDisp *= 0.0;
    //  }

    //  vectorPtr_Type emDisp0(new vector_Type( emDisp -> map() ) );
    //    *emDisp0 = *emDisp;
    //    cout << "\n\naddress emDisp: " << emDisp << ", emDisp0: " << emDisp0 << "\n\n";
    //    //return 0.0;
    //
    //      bool twoWayCoupling = parameterList.get("two_way", false);
    //
    //    if(twoWayCoupling)
    //    {
    //          if ( comm->MyPID() == 0 )
    //          {
    //              std::cout << "\nREASSEMBLING STIFFNESS MATRIX FOR TOW WAY COUPLING!\n" << std::endl;
    //          }
    //
    //           monodomain -> setupStiffnessMatrix();
    //           monodomain -> setupGlobalMatrix();
    //    }

    //===========================================================
    //===========================================================
    //              TIME LOOP
    //===========================================================
    //===========================================================
    //      Real emdt = parameterList.get("emdt",1.0);
    //      int iter((emdt / monodomain -> timeStep()));
    //      int k(0);
    //      Real saveStep = parameterList.get("save_step",1.0);
    //      int saveIter((saveStep / monodomain -> timeStep()));
    //      Real meth = parameterList.get("meth",1.0);
    //      int subiter = parameterList.get("subiter",100);
    //
    //        Real dt_min = 0.01;
    //
    //
    //      Real Ca_diastolic = dataFile( "solid/physics/Ca_diastolic", 0.02155 );
    //
    //      Real p0 = 0;
    //      Real pnm1 = pressure;
    //      Real dp = pressure - p0;
    //
    //      fluidVolume =  ComputeVolume( monodomain -> localMeshPtr(), *referencePositionX, *solidDisp, monodomain->ETFESpacePtr(),dETFESpace,10,comm);
    //      Real V0 = fluidVolume;
    //      Real dV = fluidVolume - V0;
    //      Real dV_temporal(dV);
    //      Real dp_temporal(dp);
    //
    //      bool isoPV = true;
    //      Real res;
    //      std::string outputPV="/"+problemFolder+"/pv.txt";
    //      const char * pvOutputFile = outputPV.c_str();
    //      std::ofstream pv(pvOutputFile);
    //      if(  comm->MyPID() == 0 ){
    //      std::string outputPV="/"+problemFolder+"/pv.txt";
    //      const char * pvOutputFile = outputPV.c_str();
    //      pv << "Pressure, Volume";
    //      pv << "\n" << pressure << ", " << fluidVolume;
    //      }
    //
    //  Int pv_phase(0);
    //     for( Real t(0.0); t< monodomain -> endTime(); )
    //   {
    //        t = t + monodomain -> timeStep();
    //        k++;
    //
    //
    ////            if ( comm->MyPID() == 0 )
    ////            {
    ////                std::cout << "\nSolve REACTION step with ROS3P!\n" << std::endl;
    ////            }
    //          LifeChrono timer;
    //          timer.start();
    //          if(meth == 1.0) monodomain -> solveOneReactionStepROS3P(dt_min);
    //          else
    //              {
    //               for(int j(0); j<subiter; j++) monodomain -> solveOneReactionStepFE(subiter);
    //              }
    //
    //          timer.stop();
    //          for(int pid(0); pid < 4; pid++)
    //          {
    //              if ( comm->MyPID() == pid )
    //              {
    //                      std::cout << "\nDone in " << timer.diff() << std::endl;
    //              }
    //          }
    //
    //      if ( comm->MyPID() == 0 )
    //      {
    //          std::cout << "\nSolve DIFFUSION step with BE!\n" << std::endl;
    //      }
    //      timer.reset();
    //      timer.start();
    //
    //      (*monodomain -> rhsPtrUnique()) *= 0.0;
    //        monodomain -> updateRhs();
    //          monodomain -> solveOneDiffusionStepBE();
    //          timer.stop();
    //        for(int pid(0); pid < 4; pid++)
    //        {
    //              if ( comm->MyPID() == pid )
    //              {
    //                      std::cout << "\nDone in " << timer.diff() << std::endl;
    //              }
    //        }
    //      timer.reset();
    //
    //        *tmpRhsActivation *= 0;
    //          if ( comm->MyPID() == 0 )
    //          {
    //              std::cout << "\nASSEMBLING ACTIVATION EQUATION!\n" << std::endl;
    //          }
    //
    //
    //            if( monodomain -> globalSolution().at(3) -> normInf() >= Ca_diastolic
    //                    ||
    //                 emDisp -> normInf() >=  1e-7)
    //            {
    //
    ////                    if ( comm->MyPID() == 0 )
    ////                    {
    ////                        std::cout << "\nI SHALL SOLVE THE ACTIVATION EQUATION!\n" << std::endl;
    ////                    }
    //
    //                  {
    //                      using namespace ExpressionAssembly;
    //
    //
    //                      BOOST_AUTO_TPL(I,      value(Id) );
    //                      BOOST_AUTO_TPL(Grad_u, grad( electrodETFESpace, *emDisp, 0) );
    //                      BOOST_AUTO_TPL(F,      ( Grad_u + I ) );
    //                      //BOOST_AUTO_TPL(FmT,    minusT(F) );
    //                      BOOST_AUTO_TPL(J,       det(F) );
    //                      BOOST_AUTO_TPL(Jm23,    pow(J, -2./3));
    //                      BOOST_AUTO_TPL(I1,     dot(F, F));
    //
    //                      // Fibres
    //                      BOOST_AUTO_TPL(f0,     value( electrodETFESpace, *( monodomain -> fiberPtr() ) ) );
    //                      BOOST_AUTO_TPL(f,      F * f0 );
    //                      BOOST_AUTO_TPL(I4f,    dot(f, f) );
    //
    //
    //                      BOOST_AUTO_TPL(s0,     value(electrodETFESpace, *(solid.material() -> sheetVectorPtr() ) ) );
    //                      BOOST_AUTO_TPL(I4s,    dot(F * s0, F * s0));
    //
    //                      BOOST_AUTO_TPL(I1iso,   Jm23 * I1);
    //                      BOOST_AUTO_TPL(I4fiso,  Jm23 * I4f);
    //                      BOOST_AUTO_TPL(I4siso,  Jm23 * I4s);
    //
    //
    //                      // Generalised invariants
    //                      BOOST_AUTO_TPL(gf,  value(aETFESpace, *gammaf));
    //                      BOOST_AUTO_TPL(gs,  value(aETFESpace, *gammas));
    //                      BOOST_AUTO_TPL(gn,  value(aETFESpace, *gamman));
    ////                        BOOST_AUTO_TPL(gs,  pow(gf + value(1.0), -0.5) + value(-1.0) );
    ////                        BOOST_AUTO_TPL(gn,  gs);
    //
    //
    //                      BOOST_AUTO_TPL(dI1edI1,   value(1.0) / (   (gn + value(1.0))  *  (gn + value(1.0))  )    );
    //                      BOOST_AUTO_TPL(dI1edI4f,  value(1.0) / (   (gf + value(1.0))  *  (gf + value(1.0))  ) -  value(1.0) / (   (gn + value(1.0))  *  (gn + value(1.0))  ) );
    //                      BOOST_AUTO_TPL(dI1edI4s,  value(1.0) / (   (gs + value(1.0))  *  (gs + value(1.0))  ) -  value(1.0) / (   (gn + value(1.0))  *  (gn + value(1.0))  ) );
    //                      BOOST_AUTO_TPL(dI4fedI4f, value(1.0) / (   (gf + value(1.0)) *   (gf + value(1.0))  )    );
    //
    //                      BOOST_AUTO_TPL(I1eiso,   dI1edI1 * I1iso + dI1edI4f * I4fiso + dI1edI4s * I4siso);
    //                      BOOST_AUTO_TPL(I4feiso,  dI4fedI4f * I4fiso);
    //                      BOOST_AUTO_TPL(I4feisom1, ( I4feiso - value(1.0) ) );
    //
    //                      Real A = dataFile( "solid/physics/a", 4960 );
    //                      Real B = dataFile( "solid/physics/b_activation", 0. );
    //                      Real Af = dataFile( "solid/physics/af", 0. );
    //                      Real Bf = dataFile( "solid/physics/bf", 0. );
    //                      //cout << "\n\nparameters: a: " << A << ", af: " << Af << ", b: " << B << ", bf: " << Bf << "\n\n";
    //                      BOOST_AUTO_TPL(a, value( A ) );
    //                      BOOST_AUTO_TPL(b, value(  B ) );
    //                      BOOST_AUTO_TPL(af, value( Af ) );
    //                      BOOST_AUTO_TPL(bf, value(  Bf ) );
    //                      BOOST_AUTO_TPL(psi_iso_e, a *  Jm23 * pow( eval(EXP, ( I1eiso + value(-3.0) ) ), B ) );
    //                      BOOST_AUTO_TPL(psi_f_e,  value(2.0) * af * eval( H, I4feisom1 ) /* ( I4feiso + value(-1.0) )*/  * pow( eval(EXP2, ( I4feiso + value(-1.0) ) ), Bf ) );
    //
    //                      //initial
    //                      BOOST_AUTO_TPL(Grad_u_i, grad( electrodETFESpace, *emDisp0, 0) );
    //                      BOOST_AUTO_TPL(F_i,      ( Grad_u_i + I ) );
    //                      BOOST_AUTO_TPL(J_i,       det(F_i) );
    //                      BOOST_AUTO_TPL(Jm23_i,    pow(J_i, -2./3));
    //                      BOOST_AUTO_TPL(I1_i,     dot(F_i, F_i));
    //                      BOOST_AUTO_TPL(I1iso_i,   Jm23_i * I1_i);
    //                      // Fibres
    //                      BOOST_AUTO_TPL(f_i,      F_i * f0 );
    //                      BOOST_AUTO_TPL(I4f_i,    dot(f_i, f_i) );
    //                      BOOST_AUTO_TPL(I4fiso_i,  Jm23_i * I4f_i);
    //                      BOOST_AUTO_TPL(I1eiso_i,   I1iso_i );
    //                      BOOST_AUTO_TPL(I4feiso_i,  I4fiso_i);
    //                      BOOST_AUTO_TPL(I4feisom1_i, ( I4feiso_i - value(1.0) ) );
    //                      BOOST_AUTO_TPL(psi_iso_e_i, a *  Jm23 * pow( eval(EXP, ( I1eiso_i + value(-3.0) ) ), B ) );
    //                      BOOST_AUTO_TPL(psi_f_e_i,  value(2.0) * af * eval( H, I4feisom1_i ) /* ( I4feiso + value(-1.0) )*/  * pow( eval(EXP2, ( I4feiso_i + value(-1.0) ) ), Bf ) );
    //
    //                      //BOOST_AUTO_TPL(dW01, value(-1.0) * a + eval( H, I4feisom1_i ) *  value(-2.0) * af * ( I4feiso + value(-1.0) ) ) ;
    ////                        BOOST_AUTO_TPL(dW0, value(-1.0) * ( psi_iso_e_i + psi_f_e_i ) * I4fiso_i ) ;
    //  //                  BOOST_AUTO_TPL(dW, value(-1.0) * ( psi_iso_e + psi_f_e ) * I4fiso * pow(gf + value(1.0), -3) );
    //                      BOOST_AUTO_TPL(dW0, value(-2.0) * I4fiso_i) ;
    //                      BOOST_AUTO_TPL(dW, value(2.0) * I4fiso * ( value(3.0) * gf + value(-6.0) * gf * gf + value(10.0) * gf * gf * gf + value(-15.0) * gf * gf * gf * gf  + value(21.0) * gf * gf * gf * gf * gf) );
    //
    //
    //                      BOOST_AUTO_TPL(Ca,    value( aETFESpace, *( monodomain -> globalSolution().at(3)  ) ) );
    //                      BOOST_AUTO_TPL(Ca2, Ca * Ca );
    //
    //                      BOOST_AUTO_TPL(dCa, ( Ca - value(Ca_diastolic) ) );
    //                  //    Real alpha1 = -2.5;
    //                      Real active_coefficient = dataFile( "solid/physics/active_coefficient", -2.5 );
    //                  //  BOOST_AUTO_TPL(coeff, a * Jm23 * pow( eval(EXP, ( I1iso + value(-3.0) ) ), B )/*+ value(2.0) * af * eval( H, I4feisom1)*/ );
    //                      //BOOST_AUTO_TPL(Pa, value(active_coefficient) * psi_iso_e  * eval(H, dCa) * eval(H, dCa) * eval(fl, I4fiso) + dW0 );
    //                      BOOST_AUTO_TPL(Pa, value(active_coefficient) * eval(H, dCa) * eval(H, dCa) * eval(fl, I4fiso) );
    //                      //BOOST_AUTO_TPL(Pa, value(active_coefficient) * psi_iso_e  * eval(H, dCa) * eval(H, dCa) * eval(fl, I4fiso) + dW0 );
    //
    //                      //BOOST_AUTO_TPL(Pag, value(active_coefficient) * a * eval(H, dCa) * eval(H, dCa) * eval(flg, value( aETFESpace, *gammaf ) ) + dW0 );
    //                    //  Real delta = 0.001;
    //                      Real viscosity = dataFile( "solid/physics/viscosity", 0.0005 );
    //                      BOOST_AUTO_TPL(beta, value(viscosity ) /*/ coeff  /* Jm23 * pow( eval(EXP, ( I1iso + value(-3.0) ) ), -B )*/ );
    //
    ////                        BOOST_AUTO_TPL(dWs, a * pow(g, -1) );
    //                      BOOST_AUTO_TPL(gamma_dot, beta / ( Ca2 ) * ( Pa - dW )  );
    //
    //
    //                      timer.start();
    //                      integrate ( elements ( monodomain -> localMeshPtr() ),
    //                              monodomain -> feSpacePtr() -> qr() ,
    //                              monodomain -> ETFESpacePtr(),
    //                              gamma_dot  * phi_i
    //                      ) >> tmpRhsActivation;
    //
    //                  }
    //                  timer.stop();
    ////                    for(int pid(0); pid < 4; pid++)
    ////                    {
    ////                            if ( comm->MyPID() == pid )
    ////                            {
    ////                                    std::cout << "\nDone in " << timer.diff() << std::endl;
    ////                            }
    ////                    }
    //                  timer.reset();
    //
    //                  *rhsActivation *= 0;
    //                  *rhsActivation = ( *(mass) * ( *gammaf ) );
    //                  *rhsActivation += ( ( monodomain -> timeStep() * *tmpRhsActivation ) );
    //
    //                  linearSolver.setRightHandSide(rhsActivation);
    //
    //
    //
    //                  if ( comm->MyPID() == 0 )
    //                  {
    //                      std::cout << "\nSOLVING ACTIVATION EQUATION!\n" << std::endl;
    //                  }
    //
    //
    //                  linearSolver.solve(gammaf);
    //            }
    //            else *gammaf *= 0.0;
    //
    ////                if( gammaf -> maxValue() > 0.0)
    ////                {
    ////                    std::cout << "\n*****************************************************";
    ////                    std::cout << "\nCannot be positive: " ;
    ////                    std::cout << "\n*****************************************************";
    ////                    return 0;
    ////                    int d = gammaf -> epetraVector().MyLength();
    ////                    int size =  gammaf -> size();
    ////                    for(int l(0); l < d; l++)
    ////                    {
    ////                        if ( comm->MyPID() == 0 )
    ////                        {
    ////                            std::cout << "\n*****************************************************";
    ////                            std::cout << "\nChanging the gamma back to zero: " ;
    ////                            std::cout << "\n*****************************************************";
    ////
    ////                        }
    ////                        int m1 = gammaf -> blockMap().GID(l);
    ////                        //cout << m1 << "\t" << size << "\n";
    ////                        if( (*gammaf)[m1] > 0)
    ////                        {
    ////                            (*gammaf)[m1] = 0.0;
    ////                        }
    ////                        std::cout << "\n";
    ////
    ////                    }
    ////                }
    //
    //
    //            if ( k % iter == 0)
    //            {
    //                  if(usingDifferentMeshes)
    //                  {
    //                      if(gammaf -> normInf() > 0)
    //                      {
    //                          if ( comm->MyPID() == 0 )
    //                          {
    //                              std::cout << "\nINTERPOLATING FROM FINE TO COARSE!\n" << std::endl;
    //                          }
    //
    //                          F2C -> updateRhs( gammaf );
    //                          F2C -> interpolate();
    //                          F2C -> solution( solidGammaf );
    //
    //                      }
    //                      else *solidGammaf *= 0.0;
    //                  }
    //
    //
    //                  solid.material() -> setGammaf( *solidGammaf );
    //

    //
    //                  if ( comm->MyPID() == 0 )
    //                  {
    //                      std::cout << "\nSOLVING STATIC MECHANICS!\n" << std::endl;
    //                  }
    //
    //
    //
    //                   res = 2 * fluidVolume;
    //                   //
    //                   p0 = pressure;
    //                   if( pv_phase == 0 && pressure > 100000)
    //                  {
    //                      isoPV = false;
    //                      pv_phase = 1;
    //                  }
    //                  if( pv_phase == 1 && dV_temporal > 0 )
    //                  {
    //                      isoPV = true;
    //                      pv_phase = 2;
    //                      V0 = fluidVolume;
    //                  }
    //                  if( pv_phase == 2 && pressure < 12000 )
    //                  {
    //                      isoPV = false;
    //                      pv_phase = 3;
    //                  }
    //                  if( pv_phase == 3 && dV_temporal < 0 )
    //                  {
    //                      isoPV = true;
    //                      pv_phase = 0;
    //                      V0 = fluidVolume;
    //                  }
    //                      if ( comm->MyPID() == 0 )
    //                      {
    //                      std::cout  << "\nisoPV: " << isoPV << " \n";
    //                      }
    //
    //
    //
    //                  if ( comm->MyPID() == 0 )
    //                  {
    //                      std::cout << "\n*****************************************************";
    //                      std::cout << "\nWE ARE AT TIME: "<< t << ", PV PHASE: " << pv_phase;
    //                      std::cout << "\n*****************************************************";
    //
    //                  }
    //
    //                   if(isoPV)
    //                   {
    //                      if ( comm->MyPID() == 0 )
    //                      {
    //                      std::cout  << "\nEntering in the loop (dV / V): " << dV / V0 << " \n";
    //                      }
    //                      int iter(0);
    //                      while( res > pvtol )
    //                      {
    //                          iter++;
    //
    //                  if ( comm->MyPID() == 0 )
    //                  {
    //                      std::cout << "\n*****************************************************";
    //                      std::cout << "\nITERATION: "<< iter << ", TIME: " << t ;
    //                      std::cout << "\n*****************************************************\n";
    //
    //                  }
    //
    //                          solid.iterate ( solidBC -> handler() );
    //                          *solidDisp = solid.displacement();
    //
    //                          if ( comm->MyPID() == 0 )
    //                          {
    //                              std::cout << "\nisoPV: Initial Volume of Fluid: ... " <<  V0;
    //                          }
    //                          Real Vn = fluidVolume;
    //                          fluidVolume =  ComputeVolume( monodomain -> localMeshPtr(), *referencePositionX, *solidDisp, monodomain->ETFESpacePtr(),dETFESpace,10,comm);
    //                          //Real newFluidVolume = 1;//ComputeVolume(*fullSolidMesh, *referencePosition, *solidDisp, 10, comm);
    //
    //                          dV = fluidVolume - V0;
    //                          //fluidVolume = newFluidVolume;
    //                          if ( comm->MyPID() == 0 )
    //                          {
    //                              std::cout  << "\nVariation %: " << dV / V0 << " \n";
    //                          }
    //
    //
    //                          if ( comm->MyPID() == 0 )
    //                          {
    //                              std::cout << "\nComputing pressure: ... ";
    //                          }
    //                          Real newPressure = evaluatePressure(fluidVolume, dV, pressure, dp_temporal, dV, Cp, pv_phase);
    //                          dp = newPressure - pressure;
    //                          pressure = newPressure;
    //                          if ( comm->MyPID() == 0 )
    //                          {
    //                              std::cout <<  newPressure <<", Variation %: " << dp / pressure << " \n";
    //                          }
    //                          endoVec = -newPressure;
    //                          pEndo.reset( ( new BCVector (endoVec, dFESpace -> dof().numTotalDof(), 1) ) );
    //                          solidBC -> handler() -> modifyBC(10, *pEndo);
    //
    //                          res = std::abs(dV) / V0;
    //                          if ( comm->MyPID() == 0 )
    //                          {
    //                              std::cout <<  "\nresidual: " << res <<"\n";
    //                          }
    //
    //                          if(pressure > 100000) break;
    //                          if(pressure < 12000) break;
    //                      }
    //                      dV_temporal = fluidVolume  - V0;
    //                      dp_temporal = pressure - p0;
    //                   }
    //                   else
    //                   {
    //                          Real Vn = fluidVolume;
    //                          solid.iterate ( solidBC -> handler() );
    //                          *solidDisp = solid.displacement();
    //
    //                          if ( comm->MyPID() == 0 )
    //                          {
    //                              std::cout << "\nno isoPV: Computing Volume of Fluid: ... ";
    //                          }
    //                          fluidVolume =  ComputeVolume( monodomain -> localMeshPtr(), *referencePositionX, *solidDisp, monodomain->ETFESpacePtr(),dETFESpace,10,comm);//Real newFluidVolume = 1.0;//ComputeVolume(*fullSolidMesh, *referencePosition, *solidDisp, 10, comm);
    //                          dV = fluidVolume - Vn;
    //                          //fluidVolume = newFluidVolume;
    //                          if ( comm->MyPID() == 0 )
    //                          {
    //                              std::cout <<  fluidVolume <<", Variation %: " << dV / fluidVolume << " \n";
    //                          }
    //
    //
    //                          if ( comm->MyPID() == 0 )
    //                          {
    //                              std::cout << "\nComputing pressure: ... ";
    //                          }
    //                          Real newPressure = evaluatePressure(fluidVolume, dV, pressure, dp_temporal, dV_temporal, Cp, pv_phase, t);
    //                          dp = newPressure - pressure;
    //                          pressure = newPressure;
    //                          if ( comm->MyPID() == 0 )
    //                          {
    //                              std::cout <<  newPressure <<", Variation %: " << dp / pressure << " \n";
    //                          }
    //                          endoVec = -newPressure;
    //                          pEndo.reset( ( new BCVector (endoVec, dFESpace -> dof().numTotalDof(), 1) ) );
    //                          solidBC -> handler() -> modifyBC(10, *pEndo);
    //
    //                          dV_temporal = fluidVolume  - Vn;
    //                          dp_temporal = pressure - p0;
    //                   }
    //
    //                      *savePressure = -pressure;
    //                      *saveVolume = fluidVolume;
    //
    //                  //        timeAdvance->shiftRight ( solid.displacement() );
    //
    //
    //                  /*if(parameterList.get("time_prestretch",false))
    //                  {
    //                      if( monodomain -> globalSolution().at(3)-> maxValue() < Ca_diastolic)
    //                      {
    //                          (*emDisp0) = (*emDisp);
    //                      }
    //                      else if( monodomain -> globalSolution().at(3)-> minValue() < Ca_diastolic)
    //                      {
    //                          int d = monodomain -> globalSolution().at(3) -> epetraVector().MyLength();
    ////                            int size =  monodomain -> globalSolution().at(3) -> size();
    //                          for(int l(0); l < d; l++)
    //                          {
    //                              if ( comm->MyPID() == 0 )
    //                              {
    //                                  std::cout << "\n*****************************************************";
    //                                  std::cout << "\nChanging the initial displacement: " ;
    //                                  std::cout << "\n*****************************************************";
    //
    //                              }
    //                              int m1 = monodomain -> globalSolution().at(3) -> blockMap().GID(l);
    //                              //cout << m1 << "\t" << size << "\n";
    //                              if( (*(monodomain -> globalSolution().at(3)))[m1] <= Ca_diastolic)
    //                              {
    //                                  std::cout << "\nchanging at: " << m1 ;
    //                                  int m2 = solid.displacementPtr() -> blockMap().GID(l + d);
    //                                  int m3 = solid.displacementPtr() -> blockMap().GID(l + 2 * d);
    //                                  std::cout << "\n" << (*emDisp0)[m1] << " has become: ";
    //                                  (*emDisp0)[m1] = (*emDisp)[m1];
    //                                  std::cout << (*emDisp0)[m1] << "\n";
    //                                  (*emDisp0)[m2] = (*emDisp)[m2];
    //                                  (*emDisp0)[m3] = (*emDisp)[m3];
    //                              }
    //                              std::cout << "\n";
    //
    //                          }
    //                      }
    //                  }*/
    //                  if(usingDifferentMeshes)
    //                  {
    //                      if(solid.displacementPtr() -> normInf() > 0)
    //                      {
    //                          if ( comm->MyPID() == 0 )
    //                          {
    //                              std::cout << "\nINTERPOLATING FROM COARSE TO FINE!\n" << std::endl;
    //                          }
    //
    //                          C2F -> updateRhs( solid.displacementPtr() );
    //                          C2F -> interpolate();
    //                          C2F -> solution( emDisp );
    //                      }
    //                      else *emDisp *= 0.0;
    //                  }
    //
    //
    //
    //                    if(twoWayCoupling)
    //                    {
    //                          if ( comm->MyPID() == 0 )
    //                          {
    //                              std::cout << "\nREASSEMBLING STIFFNESS MATRIX FOR TWO WAY COUPLING!\n" << std::endl;
    //                          }
    //
    //                           monodomain -> setupStiffnessMatrix();
    //                           monodomain -> setupGlobalMatrix();
    //
    //                    }
    //            }
    //
    //                          if ( comm->MyPID() == 0 )
    //                          {
    //          pv << "\n" << pressure << ", " << fluidVolume;
    //  }
    //        //cout << "\n\n save every " << saveIter << "iteration\n";
    //        if ( k % saveIter == 0)
    //        {
    //            monodomain -> exportSolution(expElectro, t);
    //            expGammaf.postProcess(t);
    //            exporter->postProcess ( t );
    //        }
    //
    //        //emDisp -> spy("interpolatedDisplacement");
    //        //solid.displacementPtr() -> spy("displacement");
    //
    ////          if ( comm->MyPID() == 0 )
    ////          {
    ////              std::cout << "\n===================="
    ////              std::cout << "\nWE ARE AT TIME: " << t ;
    ////              std::cout << "\n===================="
    ////
    ////          }
    //
    //
    //      }
    //  if ( comm->MyPID() == 0 )
    //  {
    //          pv.close();
    //  }
    //
    //      expElectro.closeFile();
    //      expGammaf.closeFile();
    //
    //      exporter -> closeFile();
    //
    //
    //
    //    if ( comm->MyPID() == 0 )
    //    {
    //        std::cout << "\nActive strain example: Passed!" << std::endl;
    //    }
    //
    //
    //    bool save4restart = parameterList.get("save4restart", true);
    //    if(save4restart)
    //    {
    //      ExporterHDF5< RegionMesh <LinearTetra> > restartV0;
    //      ExporterHDF5< RegionMesh <LinearTetra> > restartV1;
    //      ExporterHDF5< RegionMesh <LinearTetra> > restartV2;
    //      ExporterHDF5< RegionMesh <LinearTetra> > restartV3;
    //      ExporterHDF5< RegionMesh <LinearTetra> > restartG;
    //      ExporterHDF5< RegionMesh <LinearTetra> > restartD;
    //
    //      restartV0.setMeshProcId(monodomain -> localMeshPtr(), comm->MyPID());
    //      restartV1.setMeshProcId(monodomain -> localMeshPtr(), comm->MyPID());
    //      restartV2.setMeshProcId(monodomain -> localMeshPtr(), comm->MyPID());
    //      restartV3.setMeshProcId(monodomain -> localMeshPtr(), comm->MyPID());
    //      restartG.setMeshProcId(monodomain -> localMeshPtr(), comm->MyPID());
    //      restartD.setMeshProcId(localSolidMesh, comm->MyPID());
    //
    //      restartV0.setPrefix(parameterList.get ("variable0", "V0"));
    //      restartV1.setPrefix(parameterList.get ("variable1", "V1"));
    //      restartV2.setPrefix(parameterList.get ("variable2", "V2"));
    //      restartV3.setPrefix(parameterList.get ("variable3", "V3"));
    //      restartG.setPrefix(parameterList.get ("gamma_exporter", "A"));
    //      restartD.setPrefix(parameterList.get ("displacement_exporter", "D"));
    //
    //      restartV0.setPostDir ( "./save_restart/" );
    //      restartV1.setPostDir ( "./save_restart/" );
    //      restartV2.setPostDir ( "./save_restart/" );
    //      restartV3.setPostDir ( "./save_restart/" );
    //      restartG.setPostDir ( "./save_restart/" );
    //      restartD.setPostDir ( "./save_restart/" );
    //
    //      restartV0.addVariable ( ExporterData<RegionMesh<LinearTetra> >::ScalarField, "V0", monodomain -> feSpacePtr(), monodomain -> globalSolution().at(0), UInt (0) );
    //      restartV1.addVariable ( ExporterData<RegionMesh<LinearTetra> >::ScalarField, "V1", monodomain -> feSpacePtr(), monodomain -> globalSolution().at(1), UInt (0) );
    //      restartV2.addVariable ( ExporterData<RegionMesh<LinearTetra> >::ScalarField, "V2", monodomain -> feSpacePtr(), monodomain -> globalSolution().at(2), UInt (0) );
    //      restartV3.addVariable ( ExporterData<RegionMesh<LinearTetra> >::ScalarField, "V3", monodomain -> feSpacePtr(), monodomain -> globalSolution().at(3), UInt (0) );
    //      restartG.addVariable ( ExporterData<RegionMesh<LinearTetra> >::ScalarField, "G", monodomain -> feSpacePtr(), gammaf, UInt (0) );
    //      restartD.addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, "D", dFESpace, solidDisp, UInt (0) );
    //
    //      restartV0.closeFile();
    //      restartV1.closeFile();
    //      restartV2.closeFile();
    //      restartV3.closeFile();
    //      restartG.closeFile();
    //      restartD.closeFile();
    //    }

#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return 0;
}

