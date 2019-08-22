/*
 * EssentialPatchBCMovingPlane.h
 *
 *  Created on: Apr 5, 2019
 *      Author: pamstad
 */

#ifndef LIFEV_EM_EXAMPLES_EXAMPLE_EMHEART_ESSENTIALPATCHBCMOVINGPLANE_H_
#define LIFEV_EM_EXAMPLES_EXAMPLE_EMHEART_ESSENTIALPATCHBCMOVINGPLANE_H_

#include <stdio.h>
#include <lifev/em/examples/example_EMHeart/EssentialPatchBC.hpp>
#include <lifev/em/examples/example_EMHeart/GenericFactory.hpp>
#include <vector>

#define PI 3.14159265359

namespace LifeV
{

class EssentialPatchBCMovingPlane : public EssentialPatchBC {
public:


	virtual void setup(const GetPot& dataFile, const std::string& name); //In the setup function the basic things come in it, like as name of patch, flag, direction vector, displacement vector and so on

	virtual const bool nodeOnPatch(const Vector3D& coord, const Real& time);

	virtual void modifyPatchArea(EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver,const int& newFlag, const Real& time);

	virtual vectorPtr_Type directionalVectorField (EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver,const boost::shared_ptr<FESpace<RegionMesh<LinearTetra>, MapEpetra >> dFeSpace, Vector3D& direction, const Real& disp, const Real& time);

	 virtual const bool nodeOnPatchCurrent(const Vector3D& coord,const Real& time);

private:
	//vectorPtr_Type p2PositionVector;
	//LifeV::EssentialPatchBC::vector_Type p2PositionVector;
	//vector_Type p2PositionVector;
	vectorPtr_Type m_p2currentPositionVector;
	Real m_maxDisplacement;
	Vector3D normal_vector;
	Vector3D starting_point;

};

REGISTER(EssentialPatchBC, EssentialPatchBCMovingPlane);

}
#endif /* LIFEV_EM_EXAMPLES_EXAMPLE_EMHEART_ESSENTIALPATCHBCMOVINGPLANE_H_ */
