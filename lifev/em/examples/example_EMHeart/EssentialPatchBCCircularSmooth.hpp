//
//  EssentialPatchBCCircularSmooth.hpp
//  lifev-heart
//
//  Created by Thomas Kummer on 30.06.18.
//  Copyright © 2018 Thomas Kummer. All rights reserved.
//

#ifndef EssentialPatchBCCircularSmooth_hpp
#define EssentialPatchBCCircularSmooth_hpp

#include <stdio.h>
#include <lifev/em/examples/example_EMHeart/EssentialPatchBC.hpp>

#define PI 3.14159265359

namespace LifeV
{
    
class EssentialPatchBCCircularSmooth : public EssentialPatchBC
{
public:
    
    EssentialPatchBCCircularSmooth(){}
    ~EssentialPatchBCCircularSmooth(){}
    
    virtual void setup(const GetPot& dataFile, const std::string& name)
    {        
        super::setup(dataFile, name);

        m_Radius= dataFile ( ("solid/boundary_conditions/" + m_Name + "/radius").c_str(), 1.0 );
        
        m_EdgeDispFactor = dataFile ( ("solid/boundary_conditions/" + m_Name + "/edgeDispFactor").c_str(), 0.0 );

        for ( UInt j (0); j < 3; ++j )
        {
            m_Center[j] = dataFile ( ("solid/boundary_conditions/" + m_Name + "/center").c_str(), 0.0, j );
        }
        
        m_tmax = dataFile ( "solid/patches/tmax", 0. );
        m_tduration = dataFile ( "solid/patches/tduration", 0. );
    }
    
    
protected:
    
    virtual vectorPtr_Type directionalVectorField (EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver ,const boost::shared_ptr<FESpace<RegionMesh<LinearTetra>, MapEpetra >> dFeSpace, Vector3D& direction, const Real& disp, const Real& time)
    {
        // auto p2PositionVector = p2PositionVectorDisplaced(dFeSpace);
        auto p2PositionVector = p2PositionVectorInitial(dFeSpace, solver);

        vectorPtr_Type p2PatchDisplacement (new VectorEpetra( dFeSpace->map(), Repeated ));
        auto nCompLocalDof = p2PatchDisplacement->epetraVector().MyLength() / 3;

        direction.normalize();
        
        for (int j (0); j < nCompLocalDof; ++j)
        {
            // Get coordiantes
            UInt iGID = p2PatchDisplacement->blockMap().GID (j);
            UInt jGID = p2PatchDisplacement->blockMap().GID (j + nCompLocalDof);
            UInt kGID = p2PatchDisplacement->blockMap().GID (j + 2 * nCompLocalDof);

/////////////THIS IS MY CHANGE
//
	//   Vector3D coordt = dFeSpace->mesh()->point(j).coordinates();		




            Vector3D coord;
            coord(0) = p2PositionVector[iGID];
            coord(1) = p2PositionVector[jGID];
            coord(2) = p2PositionVector[kGID];
            

	//	std::cout << "These are coordinates in directionali vectorfield with index j: " << coordt[0] << "          " << coordt[1] << "       " << coordt[2] << std::endl;


            // Radial and axial distance to center line
            auto patchAxis = m_Center + 1.0 * direction;
            auto radialDistance = ( (coord - m_Center).cross(coord - patchAxis) ).norm(); // is 1.0 / (m_Center - patchAxis).norm(); //this gives us the radial direction of a point
            auto axialDistanceToCenter = (coord - m_Center).dot(direction); // * direction;
            
            // If coordiantes inside or outside of a certain radius
            auto displacement = (m_EdgeDispFactor * disp - disp) * std::pow(radialDistance / m_Radius, 2.0) + disp;

            // Scale the direction vector // ternary operator; if absolute value of displacement is smaller then 5 then we multiply with displacement otherwise with zero
            auto displacementVec =  direction * ( std::abs(displacement) < 5.0 ? displacement : 0.0 ); //this is ternary operator
            (*p2PatchDisplacement)[iGID] = displacementVec[0];
            (*p2PatchDisplacement)[jGID] = displacementVec[1];
            (*p2PatchDisplacement)[kGID] = displacementVec[2];
        }

        return p2PatchDisplacement;
    }
        
    
    virtual const bool nodeOnPatch(const Vector3D& coord, const Real& time)
    {
        bool pointInCircle = (coord - m_Center).norm() < m_Radius;
        return pointInCircle;
    }
    
    
    Vector3D m_Center;
    Real m_Radius;
    Real m_EdgeDispFactor;

};

REGISTER(EssentialPatchBC, EssentialPatchBCCircularSmooth);

}

#endif /* EssentialPatchBCCircularSmooth_hpp */
