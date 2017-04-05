/*
 * RossiModel14.cpp
 *
 *  Created on: 13/mag/2014
 *      Author: srossi
 */

#include <lifev/em/solver/activation/activeStrainModels/ActiveStrainRossiModel14.hpp>




namespace LifeV
{



ActiveStrainRossiModel14::ActiveStrainRossiModel14 (  Real inverseViscosity,
		                                              Real coefficient,
                                  		              Real threshold) :
    M_inverseViscosity (inverseViscosity),
    M_activeForceCoefficient(coefficient),
    M_chemicalThreshold(threshold)

{
}


Real
ActiveStrainRossiModel14::computeActiveStress( Real i4f, Real Calcium )
{
    Real dCa = Calcium - M_chemicalThreshold;
    Real fl = EMUtility::FLRelationship( i4f );
    if ( dCa > 0.0 )
    	return M_activeForceCoefficient * dCa * dCa * fl;
    else return 0.0;
}

void
ActiveStrainRossiModel14::solveModel ( Real& timeStep )
{
	if(!M_I4fPtr)
	{
		ASSERT(false, " Cannot solve Active Strain Rossi model without I4f. ")
		exit(-1);
	}
    Int nLocalDof = M_I4fPtr->epetraVector().MyLength();

    for (int ik (0); ik < nLocalDof; ik++)
    {
        int iGID = M_I4fPtr->blockMap().GID (ik);
        Real Pa, dW;
        Real i4f = (*M_I4fPtr)[iGID];
        Real Ca = (*this->M_electroSolution.at(M_calciumIndex) )[iGID];

        Pa = computeActiveStress(i4f, Ca);

        Real g = (*M_fiberActivationPtr) [iGID];
        Real g2 = g * g;
        Real g3 = g * g2;
        Real g4 = g * g3;
        Real g5 = g * g4;

        dW = 2.0 * i4f * ( 3.0 * g - 6.0 * g2 + 10.0 * g3 - 15.0 * g4  + 21.0 * g5 );
        Real grhs = M_inverseViscosity * ( Pa - dW ) / Ca / Ca;
        grhs *= timeStep;
        (*M_fiberActivationPtr) [iGID] += grhs;
    }

}

void
ActiveStrainRossiModel14::solveModelPathology ( Real& timeStep, const VectorEpetra& f0_, const VectorEpetra& disp, const boost::shared_ptr<FESpace<RegionMesh<LinearTetra>, MapEpetra> >& feSpacePtr, boost::shared_ptr<RegionMesh<LinearTetra> > fullMeshPtr )
{
    if(!M_I4fPtr)
    {
        ASSERT(false, " Cannot solve Active Strain Rossi model without I4f. ")
        exit(-1);
    }
    
    computeI4f(I4f(), f0_, disp, feSpacePtr);
    
    Int nLocalDof = M_I4fPtr->epetraVector().MyLength();
    
    VectorSmall<3> X;
    
    for (int ik (0); ik < nLocalDof; ik++)
    {
        int iGID = M_I4fPtr->blockMap().GID (ik);
        
        X[0] = fullMeshPtr -> point(iGID).x();
        X[1] = fullMeshPtr -> point(iGID).y();
        X[2] = fullMeshPtr -> point(iGID).z();
        
        bool infarctZone = (X - M_PathologyCenter).norm() < M_PathologyRadius;
        
        Real Pa, dW;
        Real i4f = (*M_I4fPtr)[iGID];
        Real Ca = (*this->M_electroSolution.at(M_calciumIndex) )[iGID];
        
        Pa = computeActiveStress(i4f, Ca);
        Pa *= (infarctZone ? M_PathologyStrength : 1.0);
        
        Real g = (*M_fiberActivationPtr) [iGID];
        Real g2 = g * g;
        Real g3 = g * g2;
        Real g4 = g * g3;
        Real g5 = g * g4;
        
        dW = 2.0 * i4f * ( 3.0 * g - 6.0 * g2 + 10.0 * g3 - 15.0 * g4  + 21.0 * g5 );
        Real grhs = M_inverseViscosity * ( Pa - dW ) / Ca / Ca;
        grhs *= timeStep;
        (*M_fiberActivationPtr) [iGID] += grhs;

    }
    
}

void
ActiveStrainRossiModel14::computeI4f (VectorEpetra& i4f, const VectorEpetra& f0_, const VectorEpetra& disp, const boost::shared_ptr<FESpace<RegionMesh<LinearTetra>, MapEpetra> >& feSpacePtr)
{
    VectorEpetra dUdx (disp);
    VectorEpetra dUdy (disp);
    VectorEpetra dUdz (disp);
    
    dUdx = GradientRecovery::ZZGradient (feSpacePtr, disp, 0);
    dUdy = GradientRecovery::ZZGradient (feSpacePtr, disp, 1);
    dUdz = GradientRecovery::ZZGradient (feSpacePtr, disp, 2);
    
    int n = i4f.epetraVector().MyLength();
    int i (0); int j (0); int k (0);
    MatrixSmall<3,3> F; VectorSmall<3> f0;
    
    for (int p (0); p < n; p++)
    {
        i = dUdx->blockMap().GID (p);
        j = dUdx->blockMap().GID (p + n);
        k = dUdx->blockMap().GID (p + 2 * n);
        
        F(0,0) = 1.0 + dUdx[i];
        F(0,1) =       dUdy[i];
        F(0,2) =       dUdz[i];
        F(1,0) =       dUdx[j];
        F(1,1) = 1.0 + dUdy[j];
        F(1,2) =       dUdz[j];
        F(2,0) =       dUdx[k];
        F(2,1) =       dUdy[k];
        F(2,2) = 1.0 + dUdz[k];
        
        f0(0) = f0_[i];
        f0(1) = f0_[j];
        f0(2) = f0_[k];
        
        auto f = F * f0;
        i4f[i] = f.dot(f);
    }
}



void
ActiveStrainRossiModel14::setParameters(EMData& data)
{
    
//    M_PathologyCenter[0] = dataFile("activation/pathology/infarctX", 0.0);
//    M_PathologyCenter[1] = dataFile("activation/pathology/infarctY", 0.0);
//    M_PathologyCenter[2] = dataFile("activation/pathology/infarctZ", 0.0);
//    
//    M_PathologyRadius = dataFile("activation/pathology/radius", 0.0);
//    M_PathologyStrength = dataFile("activation/pathology/strength", 1.0);
    

    M_PathologyCenter[0] = data.activationParameter<Real>("PathologyX");
    M_PathologyCenter[1] = data.activationParameter<Real>("PathologyY");
    M_PathologyCenter[2] = data.activationParameter<Real>("PathologyZ");

    M_PathologyRadius = data.activationParameter<Real>("PathologyRadius");
    M_PathologyStrength = data.activationParameter<Real>("PathologyStrength");

    
    M_inverseViscosity = data.activationParameter<Real>("InverseViscosity");
    M_activeForceCoefficient = data.activationParameter<Real>("ActiveForceCoefficient");
    M_chemicalThreshold = data.activationParameter<Real>("ChemicalThreshold");
    M_calciumIndex = data.activationParameter<UInt>("CalciumIndex");
}

} /* namespace LifeV */
