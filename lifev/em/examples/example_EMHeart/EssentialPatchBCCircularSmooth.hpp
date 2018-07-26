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
        
        m_EdgeDispFactor = dataFile ( ("solid/boundary_conditions/" + m_Name + "/edgeDispFactor").c_str(), 0 );

        for ( UInt j (0); j < 3; ++j )
        {
            m_Center[j] = dataFile ( ("solid/boundary_conditions/" + m_Name + "/center").c_str(), 0, j );
        }
        
        m_tmax = dataFile ( "solid/patches/tmax", 0. );
        m_tduration = dataFile ( "solid/patches/tduration", 0. );
    }
    
protected:
    
    virtual vectorPtr_Type directionalVectorField (const boost::shared_ptr<FESpace<RegionMesh<LinearTetra>, MapEpetra >> dFeSpace, Vector3D& direction, const Real& disp, const Real& time) const
    {
        vectorPtr_Type vectorField (new VectorEpetra( dFeSpace->map(), Repeated ));
        auto nCompLocalDof = vectorField->epetraVector().MyLength() / 3;

        direction.normalize();
        // direction *= disp;

        for (int j (0); j < nCompLocalDof; ++j)
        {
            // Get coordiantes
            UInt iGID = vectorField->blockMap().GID (j);
            UInt jGID = vectorField->blockMap().GID (j + nCompLocalDof);
            UInt kGID = vectorField->blockMap().GID (j + 2 * nCompLocalDof);

            Vector3D coord;

            coord(0) = dFeSpace->mesh()->point(iGID).x() + (*m_dispPtr)[iGID];
            coord(1) = dFeSpace->mesh()->point(iGID).y() + (*m_dispPtr)[jGID];
            coord(2) = dFeSpace->mesh()->point(iGID).z() + (*m_dispPtr)[kGID];

            // Radial and axial distance to center line
            auto currentPatchCenter = m_Center + activationFunction(time) * direction;
            auto radialDistance = ( (coord - m_Center).cross(coord - currentPatchCenter) ).norm() / (m_Center - currentPatchCenter).norm();
            auto axialDistance = (coord - currentPatchCenter).dot(direction) * direction;

            // If coordiantes inside or outside of a certain radius
            auto displacement = (m_EdgeDispFactor * disp - disp) * std::pow(radialDistance / m_Radius, 2.0) + disp;

            // If patch inside or outside the structure


            // Scale the direction vector
            auto displacementVec = displacement * direction;
            (*vectorField)[iGID] = direction[0];
            (*vectorField)[jGID] = direction[1];
            (*vectorField)[kGID] = direction[2];
        }

        return vectorField;
    }
    
    virtual const bool nodeOnPatch(const Vector3D& coord) const
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
