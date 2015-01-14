/*
 * NashPanfilovModel04.cpp
 *
 *  Created on: 13/mag/2014
 *      Author: srossi
 */

#include <lifev/em/solver/activation/activeStressModels/ActiveStressNashPanfilovModel04.hpp>



namespace LifeV
{


ActiveStressNashPanfilovModel04::ActiveStressNashPanfilovModel04 (Real kTa, Real epsilon0) :
    M_kTa (kTa),
    M_epsilon0 (epsilon0)
{
}

void ActiveStressNashPanfilovModel04::multiplyByEpsilon (VectorEpetra& rhs)
{
    Real threshold = 0.05;
    VectorEpetra eps ( this->M_electroSolution.at(0)->map()  );
    eps = M_epsilon0;
    eps += (9 * (  *( this->M_electroSolution.at(0) ) > threshold) );
    rhs *= eps;
}



void ActiveStressNashPanfilovModel04::solveModel ( Real& timeStep)
{
    VectorEpetra rhs ( *( this->M_electroSolution.at(0) ) );
    rhs *= M_kTa;
    rhs -= *super::M_fiberActivationPtr;
    multiplyByEpsilon (rhs);
    rhs *= (timeStep);
    *super::M_fiberActivationPtr += rhs;
}


void
ActiveStressNashPanfilovModel04::setParameters( EMData& data)
{
	M_kTa = data.activationParameter<Real>("kTa");
	M_epsilon0 = data.activationParameter<Real>("epsilon0");
}


void
ActiveStressNashPanfilovModel04::setupActivationPtrs(	vectorPtr_Type& fiberActivationPtr,
												vectorPtr_Type& sheetActivationPtr,
												vectorPtr_Type& normalActivationPtr )
{
	fiberActivationPtr = this->M_fiberActivationPtr;
	sheetActivationPtr.reset();
	normalActivationPtr.reset();
}

void
ActiveStressNashPanfilovModel04::updateActivation(	vectorPtr_Type& fiberActivationPtr,
											vectorPtr_Type& sheetActivationPtr,
											vectorPtr_Type& normalActivationPtr )
{
	setupActivationPtrs(fiberActivationPtr, sheetActivationPtr, normalActivationPtr);
}


} /* namespace LifeV */
