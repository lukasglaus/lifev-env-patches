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
ActiveStrainRossiModel14::solveModelPathology ( Real& timeStep, boost::shared_ptr<RegionMesh<LinearTetra> > fullMeshPtr, const boost::shared_ptr<FESpace<RegionMesh<LinearTetra>, MapEpetra >> dFeSpace)
{
    if(!M_I4fPtr)
    {
        ASSERT(false, " Cannot solve Active Strain Rossi model without I4f. ")
        exit(-1);
    }
        
    Int nLocalDof = M_I4fPtr->epetraVector().MyLength();

    auto positionVector = undeformedPositionVector(fullMeshPtr, dFeSpace);
    
    VectorSmall<3> X;
    
    for (int ik (0); ik < nLocalDof; ik++)
    {
        UInt iGID = positionVector.blockMap().GID (ik);
        UInt jGID = positionVector.blockMap().GID (ik + nLocalDof);
        UInt kGID = positionVector.blockMap().GID (ik + 2 * nLocalDof);
        
        std::cout << iGID << "/" << M_I4fPtr->blockMap().GID(ik) << std::endl;
        
        X[0] = positionVector[iGID];
        X[1] = positionVector[jGID];
        X[2] = positionVector[kGID];
        
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

    
const VectorEpetra
ActiveStrainRossiModel14::undeformedPositionVector (boost::shared_ptr<RegionMesh<LinearTetra> > fullMeshPtr, const boost::shared_ptr<FESpace<RegionMesh<LinearTetra>, MapEpetra >> dFeSpace) const
{
    // New P1 Space
    FESpace<RegionMesh<LinearTetra> , MapEpetra > p1FESpace ( dFeSpace->mesh(), "P1", 3, dFeSpace->map().commPtr() );
    
    // Create P1 VectorEpetra
    VectorEpetra p1PositionVector (p1FESpace.map());
    
    // Fill P1 vector with mesh values
    Int p1nCompLocalDof = p1PositionVector.epetraVector().MyLength() / 3;
    for (int j (0); j < p1nCompLocalDof; j++)
    {
        UInt iGID = p1PositionVector.blockMap().GID (j);
        UInt jGID = p1PositionVector.blockMap().GID (j + p1nCompLocalDof);
        UInt kGID = p1PositionVector.blockMap().GID (j + 2 * p1nCompLocalDof);
        
        p1PositionVector[iGID] = fullMeshPtr->point (iGID).x();
        p1PositionVector[jGID] = fullMeshPtr->point (iGID).y();
        p1PositionVector[kGID] = fullMeshPtr->point (iGID).z();
    }
    
    // Interpolate position vector from P1-space to current space
    VectorEpetra positionVector ( dFeSpace->map() );
    positionVector = dFeSpace -> feToFEInterpolate(p1FESpace, p1PositionVector);
    
    return positionVector;
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
