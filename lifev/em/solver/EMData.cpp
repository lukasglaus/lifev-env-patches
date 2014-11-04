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
	setupMaterialParamters(dataFile);
}

void
EMData::setupMaterialParamters(GetPot& dataFile, const std::string& section)
{
    std::string EMpassiveMaterialType = dataFile ( ( section + "/physics/EMPassiveMaterialType" ).data(), "NO_DEFAULT_PASSIVE_TYPE" );
    M_materialParametersList.set ("EMPassiveMaterialType", EMpassiveMaterialType);

    double bulkModulus = dataFile ( ( section + "/physics/BulkModulus" ).data(), 0.0 );
    M_materialParametersList.set ("BulkModulus", bulkModulus);

    double A = dataFile ( ( section + "/physics/a" ).data(), 0.0 );
    double B = dataFile ( ( section + "/physics/b" ).data(), 0.0 );
    double Af = dataFile ( ( section + "/physics/af" ).data(), 0.0 );
    double Bf = dataFile ( ( section + "/physics/bf" ).data(), 0.0 );
    double As = dataFile ( ( section + "/physics/as" ).data(), 0.0 );
    double Bs = dataFile ( ( section + "/physics/bs" ).data(), 0.0 );
    double Afs = dataFile ( ( section + "/physics/afs" ).data(), 0.0 );
    double Bfs = dataFile ( ( section + "/physics/bfs" ).data(), 0.0 );
    M_materialParametersList.set ("a", A);
    M_materialParametersList.set ("b", B);
    M_materialParametersList.set ("af", Af);
    M_materialParametersList.set ("bf", Bf);
    M_materialParametersList.set ("as", As);
    M_materialParametersList.set ("bs", Bs);
    M_materialParametersList.set ("afs", Afs);
    M_materialParametersList.set ("bfs", Bfs);

    double C = dataFile ( ( section + "/physics/C" ).data(), 0.0 );
    double bff = dataFile ( ( section + "/physics/bff" ).data(), 0.0 );
    double bss = dataFile ( ( section + "/physics/bss" ).data(), 0.0 );
    double bnn = dataFile ( ( section + "/physics/bnn" ).data(), 0.0 );
    double bfs = dataFile ( ( section + "/physics/bfs" ).data(), 0.0 );
    double bfn = dataFile ( ( section + "/physics/bfn" ).data(), 0.0 );
    double bsn = dataFile ( ( section + "/physics/bsn" ).data(), 0.0 );
    M_materialParametersList.set ("C", C);
    M_materialParametersList.set ("bff", bff);
    M_materialParametersList.set ("bss", bss);
    M_materialParametersList.set ("bnn", bnn);
    M_materialParametersList.set ("bfs", bfs);
    M_materialParametersList.set ("bfn", bfn);
    M_materialParametersList.set ("bsn", bsn);


    double mu = dataFile ( ( section + "/physics/mu" ).data(), 0.0 );
    double lambda = dataFile ( ( section + "/physics/lambda" ).data(), 0.0 );
    M_materialParametersList.set ("mu", mu);
    M_materialParametersList.set ("lambda", lambda);


    std::string EMactiveStressMaterialType = dataFile ( ( section + "/physics/EMActiveStressMaterialType" ).data(), "NO_DEFAULT_ACTIVESTRESS_TYPE" );
    M_materialParametersList.set ("EMActiveStressMaterialType", EMactiveStressMaterialType);

    double Tmax = dataFile ( ( section + "/physics/Tmax" ).data(), 0.0 );
    M_materialParametersList.set ("MaxActiveTension", Tmax);
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
