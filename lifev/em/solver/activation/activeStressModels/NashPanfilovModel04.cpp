/*
 * NashPanfilovModel04.cpp
 *
 *  Created on: 13/mag/2014
 *      Author: srossi
 */

#include <lifev/em/solver/activation/activeStressModels/NashPanfilovModel04.hpp>



namespace LifeV {

NashPanfilovModel04::NashPanfilovModel04(const MapEpetra& map, Real kTa, Real epsilon0) :
		super (map),
	    M_kTa(kTa),
		M_epsilon0(epsilon0)
{
}

NashPanfilovModel04::NashPanfilovModel04(MapEpetra& map, Real kTa, Real epsilon0) :
		super (map),
	    M_kTa(kTa),
		M_epsilon0(epsilon0)
{
}

void NashPanfilovModel04::multiplyByEpsilon(VectorEpetra& potential, VectorEpetra& rhs)
{
	Real threshold = 0.05;
	VectorEpetra eps(potential.map());
	eps = M_epsilon0;
	eps += (9*(potential>threshold));
	rhs *= eps;
}



void NashPanfilovModel04::solveModel(VectorEpetra& potential, Real timeStep)
{
    VectorEpetra rhs( potential.map() );
    rhs = potential;
    rhs *= M_kTa;
    rhs -= *super::M_activationPtr;
    multiplyByEpsilon(potential,rhs);
    rhs *= (timeStep);
    *super::M_activationPtr += rhs;
}

} /* namespace LifeV */
