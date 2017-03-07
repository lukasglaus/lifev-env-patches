/*
 * EMData.cpp
 *
 *  Created on: 26/ott/2014
 *      Author: srossi
 */

#include <lifev/em/solver/EMData.hpp>


namespace LifeV{


void
EMData::setup(GetPot& dataFile)
{
	setupSolver(dataFile);
	setupSolidParameters(dataFile);
	setupActivationParameters(dataFile);
	setupElectrophysiologyParameters(dataFile);
}

void
EMData::setupSolidParameters(GetPot& dataFile, const std::string& section)
{
    std::string EMpassiveMaterialType = dataFile ( ( section + "/physics/EMPassiveMaterialType" ).data(), "NO_DEFAULT_PASSIVE_TYPE" );
    M_solidParametersList.set ("EMPassiveMaterialType", EMpassiveMaterialType);

    double bulkModulus = dataFile ( ( section + "/physics/BulkModulus" ).data(), 35000.0 );
    M_solidParametersList.set ("BulkModulus", bulkModulus);

    double A = dataFile ( ( section + "/physics/a" ).data(), 3330.0 );
    double B = dataFile ( ( section + "/physics/b" ).data(), 9.242 );
    double Af = dataFile ( ( section + "/physics/af" ).data(), 185350.0 );
    double Bf = dataFile ( ( section + "/physics/bf" ).data(), 15.972 );
    double As = dataFile ( ( section + "/physics/as" ).data(), 25640.0 );
    double Bs = dataFile ( ( section + "/physics/bs" ).data(), 10.446 );
    double Afs = dataFile ( ( section + "/physics/afs" ).data(), 4170.0 );
    double Bfs = dataFile ( ( section + "/physics/bfs" ).data(),11.602 );
    M_solidParametersList.set ("a", A);
    M_solidParametersList.set ("b", B);
    M_solidParametersList.set ("af", Af);
    M_solidParametersList.set ("bf", Bf);
    M_solidParametersList.set ("as", As);
    M_solidParametersList.set ("bs", Bs);
    M_solidParametersList.set ("afs", Afs);
    M_solidParametersList.set ("bfs", Bfs);

    double C = dataFile ( ( section + "/physics/C" ).data(), 2000.0 );
    double bff = dataFile ( ( section + "/physics/bff" ).data(), 35.7 );
    double bss = dataFile ( ( section + "/physics/bss" ).data(), 18.9 );
    double bnn = dataFile ( ( section + "/physics/bnn" ).data(), 13.0 );
    double bfs = dataFile ( ( section + "/physics/bfs" ).data(), 10.5 );
    double bfn = dataFile ( ( section + "/physics/bfn" ).data(), 11.4 );
    double bsn = dataFile ( ( section + "/physics/bsn" ).data(), 8.35 );
    M_solidParametersList.set ("C", C);
    M_solidParametersList.set ("bff", bff);
    M_solidParametersList.set ("bss", bss);
    M_solidParametersList.set ("bnn", bnn);
    M_solidParametersList.set ("bfs", bfs);
    M_solidParametersList.set ("bfn", bfn);
    M_solidParametersList.set ("bsn", bsn);
    double C2 = dataFile ( ( section + "/physics/C2" ).data(), 4960.0 );
    M_solidParametersList.set ("C2", C2);


    double mu = dataFile ( ( section + "/physics/mu" ).data(), 4960.0 );
    double lambda = dataFile ( ( section + "/physics/lambda" ).data(), 2000.0 );
    M_solidParametersList.set ("mu", mu);
    M_solidParametersList.set ("lambda", lambda);


    std::string EMactiveStressMaterialType = dataFile ( ( section + "/physics/EMActiveStressMaterialType" ).data(), "NO_DEFAULT_ACTIVESTRESS_TYPE" );
    M_solidParametersList.set ("EMActiveStressMaterialType", EMactiveStressMaterialType);

    std::string EMactiveStrainMaterialType = dataFile ( ( section + "/physics/EMActiveStrainMaterialType" ).data(), "NO_DEFAULT_ACTIVESTRAIN_TYPE" );
    double EMActiveStrainOrthotropicParameter = dataFile ( ( section + "/physics/EMActiveStrainOrthotropicParameter" ).data(), -666. );
    std::string EMactiveStrainType = dataFile ( ( section + "/physics/EMactiveStrainType" ).data(), "TransverselyIsotropic" );
    M_solidParametersList.set ("EMActiveStrainMaterialType", EMactiveStrainMaterialType);
    M_solidParametersList.set ("EMActiveStrainOrthotropicParameter", EMActiveStrainOrthotropicParameter);
    M_solidParametersList.set ("EMactiveStrainType", EMactiveStrainType);

    double Tmax = dataFile ( ( section + "/physics/Tmax" ).data(), 50.0 );
    M_solidParametersList.set ("MaxActiveTension", Tmax);
    double Cs = dataFile ( ( section + "/physics/Cs" ).data(), 0. );
    M_solidParametersList.set ("Cs", Cs);
    double Cn = dataFile ( ( section + "/physics/Cn" ).data(), 0. );
    M_solidParametersList.set ("Cn", Cn);

    double endtime = dataFile ( ( section + "/time_discretization/endtime" ).data(), 1.0 );
    M_solidParametersList.set ("endtime", endtime);
    double timestep = dataFile ( ( section + "/time_discretization/timestep" ).data(), 0.1 );
    M_solidParametersList.set ("timestep", timestep);



}

void
EMData::setupActivationParameters(GetPot& dataFile, const std::string& section)
{
    std::string  ActivationModel = dataFile ( ( section + "/physics/ActivationModel" ).data(), "NO_DEFAULT_ACTIVATION_MODEL" );
    M_activationParametersList.set ("ActivationModel", ActivationModel );

    double  ActiveStressBeta = dataFile ( ( section + "/physics/ActiveStress_Beta" ).data(), 2.279 );
    M_activationParametersList.set ("ActiveStress_Beta", ActiveStressBeta );
    
    double ActiveStressMu = dataFile ( ( section + "/physics/ActiveStressMu" ).data(), 1000.0 );
    M_activationParametersList.set ("ActiveStress_Mu", ActiveStressMu);

    //Devil Number, in case you do not put anything ...
    double EMactiveStrainOrthotropicParameter = dataFile ( ( section + "/physics/EMActiveStrainOrthotropicParameter" ).data(), -666. );
    M_activationParametersList.set ("EMactiveStrainOrthotropicParameter", EMactiveStrainOrthotropicParameter);

    std::string EMactiveStrainType = dataFile ( ( section + "/physics/EMActiveStrainType" ).data(), "Anisotropic" );
    M_activationParametersList.set ("EMActiveStrainType", EMactiveStrainType);
  
    double Tmax = dataFile ( ( section + "/physics/Tmax" ).data(), 50.0 );
    M_activationParametersList.set ("MaxActiveTension", Tmax);

    UInt calciumIndex = dataFile ( ( section + "/physics/CalciumIndex" ).data(), 2 );
    M_activationParametersList.set("CalciumIndex", calciumIndex);

    double inverseViscosity = dataFile ( ( section + "/physics/InverseViscosity" ).data(), 0.0001 );
    M_activationParametersList.set ("InverseViscosity", inverseViscosity);

    double activeForceCoefficient = dataFile ( ( section + "/physics/ActiveForceCoefficient" ).data(), 1.0 );
    M_activationParametersList.set ("ActiveForceCoefficient", activeForceCoefficient);

    double chemicalThreshold = dataFile ( ( section + "/physics/ChemicalThreshold" ).data(), 0.21 );
    M_activationParametersList.set ("ChemicalThreshold", chemicalThreshold);


    double kTa = dataFile ( ( section + "/physics/kTa" ).data(), 49.7 );
    M_activationParametersList.set ("kTa", kTa);

    double epsilon0 = dataFile ( ( section + "/physics/epsilon0" ).data(), 1.0 );
    M_activationParametersList.set ("epsilon0", epsilon0);


    double endtime = dataFile ( ( section + "/time_discretization/endtime" ).data(), 1.0 );
    M_activationParametersList.set ("endtime", endtime);

    double timestep = dataFile ( ( section + "/time_discretization/timestep" ).data(), 0.1 );
    M_activationParametersList.set ("timestep", timestep);
    
    double pathologyX = dataFile ( ( section + "/pathology/infarctX" ).data(), 0. );
    M_activationParametersList.set ("PathologyX", pathologyX);
 
    double pathologyY = dataFile ( ( section + "/pathology/infarctY" ).data(), 0. );
    M_activationParametersList.set ("PathologyY", pathologyY);

    double pathologyZ = dataFile ( ( section + "/pathology/infarctZ" ).data(), 0. );
    M_activationParametersList.set ("PathologyZ", pathologyZ);

    double pathologyR = dataFile ( ( section + "/pathology/radius" ).data(), 0. );
    M_activationParametersList.set ("PathologyRadius", pathologyR);
    
    double pathologyS = dataFile ( ( section + "/pathology/strength" ).data(), 0. );
    M_activationParametersList.set ("PathologyStrength", pathologyS);
    
}


void
EMData::setupElectrophysiologyParameters(GetPot& dataFile, const std::string& section)
{
    double surfaceVolumeRatio = dataFile ( ( section + "/physics/surfaceVolumeRatio" ).data(), 1400.0 );
    M_electroParametersList.set ("surfaceVolumeRatio", surfaceVolumeRatio);

    double fiberDiffusion = dataFile ( ( section + "/physics/fiberDiffusion" ).data(), 3.0 );
    M_electroParametersList.set ("fiberDiffusion", fiberDiffusion);

    double sheetDiffusion = dataFile ( ( section + "/physics/sheetDiffusion" ).data(), 3.0 );
    M_electroParametersList.set ("sheetDiffusion", sheetDiffusion);

    double normalDiffusion = dataFile ( ( section + "/physics/normalDiffusion" ).data(), 3.0 );
    M_electroParametersList.set ("normalDiffusion", normalDiffusion);

    std::string IonicModel = dataFile ( ( section + "/physics/IonicModel" ).data(), "NO_DEFAULT_IONIC_MODEL" );
    M_electroParametersList.set ("IonicModel", IonicModel);

    double initialtime = dataFile ( ( section + "/time_discretization/initialtime" ).data(), 1.0 );
    M_electroParametersList.set ("initialtime", initialtime);

    double endtime = dataFile ( ( section + "/time_discretization/endtime" ).data(), 1.0 );
    M_electroParametersList.set ("endtime", endtime);

    double timestep = dataFile ( ( section + "/time_discretization/timestep" ).data(), 0.1 );
    M_electroParametersList.set ("timestep", timestep);

    std::string elementsOrder = dataFile ( ( std::string("solid/space_discretization/order") ).data(), "P1" );
    M_electroParametersList.set ("elementsOrder", elementsOrder);

    bool LumpedMass = dataFile ( ( section + "/discretization/LumpedMass" ).data(), false);
    M_electroParametersList.set ("LumpedMass", LumpedMass);


}

void
EMData::setupSolver(GetPot& dataFile)
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

    M_solverParametersList.set ("Solver Type", solverName);
    M_solverParametersList.set ("Maximum Iterations", maxIter);
    M_solverParametersList.set ("Max Iterations For Reuse", maxIterForReuse);
    M_solverParametersList.set ("Reuse Preconditioner", reusePreconditioner);
    M_solverParametersList.set ("Quit On Failure", quitOnFailure);
    M_solverParametersList.set ("Silent", silent);
    std::string subListName = "Solver: Operator List";
    std::string subSubListName = "Trilinos: " + solverName + " List";
    M_solverParametersList.sublist (subListName).sublist (subSubListName).set ("solver", solver);
    M_solverParametersList.sublist (subListName).sublist (subSubListName).set ("conv", conv);
    M_solverParametersList.sublist (subListName).sublist (subSubListName).set ("scaling", scaling);
    M_solverParametersList.sublist (subListName).sublist (subSubListName).set ("output", output);
    M_solverParametersList.sublist (subListName).sublist (subSubListName).set ("tol", tol);
    M_solverParametersList.sublist (subListName).sublist (subSubListName).set ("max_iter", maxIter);
    M_solverParametersList.sublist (subListName).sublist (subSubListName).set ("kspace", kspace);
    M_solverParametersList.sublist (subListName).sublist (subSubListName).set ("orthog", orthog);
    M_solverParametersList.sublist (subListName).sublist (subSubListName).set ("aux_vec", auxvec);;

}


} /* namespace LifeV */
