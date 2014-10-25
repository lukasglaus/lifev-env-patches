/*
 * RossiModel14.cpp
 *
 *  Created on: 13/mag/2014
 *      Author: srossi
 */

#include <lifev/em/solver/activation/activeStressModels/ActiveStressRossiModel14.hpp>



namespace LifeV
{

ActiveStressRossiModel14::ActiveStressRossiModel14 (const MapEpetra& map, Real beta, Real mu, Real Tmax) :
    super (map),
    M_coefficientBeta (beta),
    M_coefficientMu (mu),
    M_maximumActiveTenstion (Tmax)
{
}

ActiveStressRossiModel14::ActiveStressRossiModel14 (MapEpetra& map, Real beta, Real mu, Real Tmax) :
    super (map),
    M_coefficientBeta (beta),
    M_coefficientMu (mu),
    M_maximumActiveTenstion (Tmax)
{
}

void ActiveStressRossiModel14::solveModel (VectorEpetra& potential, Real timeStep)
{
    VectorEpetra rhs ( potential.map() );
    rhs = potential;
    rhs *= rhs.operator > (0.0);
    rhs *= M_coefficientBeta;
    rhs -= 2.0 * *super::M_activationPtr;
    rhs *= (timeStep * M_maximumActiveTenstion / M_coefficientMu);
    *super::M_activationPtr += rhs;
}

} /* namespace LifeV */
