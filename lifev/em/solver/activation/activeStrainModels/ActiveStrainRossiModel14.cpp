/*
 * RossiModel14.cpp
 *
 *  Created on: 13/mag/2014
 *      Author: srossi
 */

#include <lifev/em/solver/activation/activeStrainModels/ActiveStrainRossiModel14.hpp>
#include <lifev/em/util/EMUtility.hpp>



namespace LifeV
{

ActiveStrainRossiModel14::ActiveStrainRossiModel14 (  const MapEpetra& map,
		                                              Real inverseViscosity,
		                                              Real coefficient,
                                  		              Real threshold ) :
    super (map),
    M_inverseViscosity (inverseViscosity),
    M_activeForceCoefficient(coefficient),
    M_chemicalThreshold(threshold)

{
}

ActiveStrainRossiModel14::ActiveStrainRossiModel14 (  MapEpetra& map,
		                                              Real inverseViscosity,
		                                              Real coefficient,
                                  		              Real threshold ) :
    super (map),
    M_inverseViscosity (inverseViscosity),
    M_activeForceCoefficient(coefficient),
    M_chemicalThreshold(threshold)

{
}


Real
ActiveStrainRossiModel14::computeActiveStress( Real i4f, Real Calcium)
{
    Real dCa = Calcium - M_chemicalThreshold;
    Real fl = EMUtility::FLRelationship( i4f );
    if ( dCa > 0.0 )
    	return M_activeForceCoefficient * dCa * dCa * fl;
    else return 0.0;
}


void
ActiveStrainRossiModel14::solveModel ( VectorEpetra& I4f,
									   VectorEpetra& Calcium,
									   Real timeStep,
									   ActivationAnisotropy type )
{
    Int nLocalDof = I4f.epetraVector().MyLength();

    for (int ik (0); ik < nLocalDof; ik++)
    {
        int iGID = I4f.blockMap().GID (ik);
        Real Pa, dW;
        Real i4f = I4f[iGID];
        Real Ca = Calcium[iGID];

        Pa = computeActiveStress(i4f, Ca);

        Real g = (*M_gammafPtr) [iGID];
        Real g2 = g * g;
        Real g3 = g * g2;
        Real g4 = g * g3;
        Real g5 = g * g4;

        dW = 2.0 * i4f * ( 3.0 * g - 6.0 * g2 + 10.0 * g3 - 15.0 * g4  + 21.0 * g5 );
        Real grhs = M_inverseViscosity * ( Pa - dW ) / Ca / Ca;
        grhs *= timeStep;
        (*M_gammafPtr) [iGID] += grhs;
    }

    computeTransversalActivation(type);


}


void
ActiveStrainRossiModel14::computeTransversalActivation( ActivationAnisotropy  anisotropy )
{

	switch ( anisotropy )
	{
	case TransversilyIsotropic:
		ActiveStrain::transversilyIsotropicActiveStrain( *M_gammafPtr, *M_gammasPtr, *M_gammanPtr);
		break;
	default:
		ActiveStrain::orthotropicActiveStrain(*M_gammafPtr, *M_gammasPtr, *M_gammanPtr, 3.0);
	}
}

} /* namespace LifeV */
