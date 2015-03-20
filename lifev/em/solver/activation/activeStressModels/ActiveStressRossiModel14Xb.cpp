/*
 * RossiModel14.cpp
 *
 *  Created on: 13/mag/2014
 *      Author: srossi
 */

#include <lifev/em/solver/activation/activeStressModels/ActiveStressRossiModel14.hpp>



namespace LifeV
{

ActiveStressRossiModel14::ActiveStressRossiModel14 (Real beta, Real mu, Real Tmax) :
    M_coefficientBeta (beta),
    M_coefficientMu (mu),
    M_maximumActiveTenstion (Tmax)
{
}



void
ActiveStressRossiModel14::setParameters( EMData& data )
{
	M_maximumActiveTenstion = data.activationParameter<Real>("MaxActiveTension");
	M_coefficientBeta = data.activationParameter<Real>("ActiveStress_Beta");
	M_coefficientMu = data.activationParameter<Real>("ActiveStress_Mu");
}

void
ActiveStressRossiModel14::setup(EMData& data, const MapEpetra& map)
{
	setParameters(data);
    this->M_fiberActivationPtr.reset ( new vector_Type ( map ) );
    super::setup(data, map);

}

void
ActiveStressRossiModel14::setupActivationPtrs(	vectorPtr_Type& fiberActivationPtr,
												vectorPtr_Type& sheetActivationPtr,
												vectorPtr_Type& normalActivationPtr )
{
	fiberActivationPtr = 	this->M_fiberActivationPtr;
	sheetActivationPtr.reset();
	normalActivationPtr.reset();
}

void
ActiveStressRossiModel14::updateActivation(	vectorPtr_Type& fiberActivationPtr,
											vectorPtr_Type& sheetActivationPtr,
											vectorPtr_Type& normalActivationPtr )
{
	setupActivationPtrs(fiberActivationPtr, sheetActivationPtr, normalActivationPtr);
}

void ActiveStressRossiModel14::solveModel (Real& timeStep)
{
    M_XA
    M_M
    M_TnCA

    std::vector< std::shared_ptr<VectorEpetra> > M_states;
    
    XBridge4SM xb;
    M_xb;
    
    M_xb.solveFE(M_states, this->M_electroSolution.at(j), dt);
    
    
    
    VectorEpetra rhs ( *( this->M_electroSolution.at(0) ) );
    rhs *= rhs.operator > (0.0);
    rhs *= M_coefficientBeta;
    rhs -= 2.0 * *super::M_fiberActivationPtr;
    rhs *= (timeStep * M_maximumActiveTenstion / M_coefficientMu);
    *super::M_fiberActivationPtr += rhs;
}

} /* namespace LifeV */
