//
//  EssentialPatchBCCircularSmooth
//  lifev-heart
//
//  Created by Thomas Kummer on 30.06.18.
//  Copyright © 2018 Thomas Kummer. All rights reserved.
//

#ifndef EssentialPatchBCCircular_hpp
#define EssentialPatchBCCircular_hpp

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
        m_Name = name;
        m_PrevFlag = dataFile ( ("solid/boundary_conditions/" + m_Name + "/flag").c_str(), 0 );
        m_patchDisplacement = dataFile ( ("solid/boundary_conditions/" + m_Name + "/displacement").c_str(), 1.0 );
        m_Radius= dataFile ( ("solid/boundary_conditions/" + m_Name + "/radius").c_str(), 1.0 );
        
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
        
        for (int j (0); j < nCompLocalDof; ++j)
        {
            // Get coordiantes
            UInt iGID = vectorField->blockMap().GID (j);
            UInt jGID = vectorField->blockMap().GID (j + nCompLocalDof);
            UInt kGID = vectorField->blockMap().GID (j + 2 * nCompLocalDof);
            
            Vector3D coord;
            
            coord(0) = dFeSpace->mesh().point(iGID).x() + disp(iGID);
            coord(1) = dFeSpace->mesh().point(iGID).y() + disp(jGID);
            coord(2) = dFeSpace->mesh().point(iGID).z() + disp(kGID);

            // Radial and axial distance to center line
            auto currentPatchCenter = m_Center + activationFunction(time) * direction;
            auto radialDistance = ( (coord - m_Center).cross(coord - currentPatchCenter) ).norm() / (m_Center - currentPatchCenter).norm();
            auto axialDistance = (coord - currentPatchCenter).dot(direction) * direction
            
            // If coordiantes inside or outside of a certain radius
            Real radiusDispFactor (0.5);
            auto displacement = (radiusDispFactor * disp - disp) * std::pow(radialDistance / m_Radius, 2.0) + disp;
            
            
            // If patch inside or outside the structure 
            
            // Scale the direction vector
            
            auto displacementVec = displacement * direction;
            (*vectorField)[iGID] = displacementVec[0];
            (*vectorField)[jGID] = displacementVec[1];
            (*vectorField)[kGID] = displacementVec[2];
        }
        
        
//        // New P1 Space
//        FESpace<RegionMesh<LinearTetra> , MapEpetra > p1FESpace ( M_localMeshPtr, "P1", 3, M_fullMesh.comm() );
//
//        // Create P1 VectorEpetra
//        VectorEpetra p1PositionVector (p1FESpace.map());
//
//        // Fill P1 vector with mesh values
//        Int p1nCompLocalDof = p1PositionVector.epetraVector().MyLength() / 3;
//        for (int j (0); j < p1nCompLocalDof; j++)
//        {
//            UInt iGID = p1PositionVector.blockMap().GID (j);
//            UInt jGID = p1PositionVector.blockMap().GID (j + p1nCompLocalDof);
//            UInt kGID = p1PositionVector.blockMap().GID (j + 2 * p1nCompLocalDof);
//
//            p1PositionVector[iGID] = M_fullMesh.point (iGID).x();
//            p1PositionVector[jGID] = M_fullMesh.point (iGID).y();
//            p1PositionVector[kGID] = M_fullMesh.point (iGID).z();
//        }
//
//        // Interpolate position vector from P1-space to current space
//        VectorEpetra positionVector ( disp.map() );
//        positionVector = M_FESpace -> feToFEInterpolate(p1FESpace, p1PositionVector);
//
//        // Add displacement to position vector
//        positionVector += disp;

        
        
        return vectorField;
    }
    
    virtual const bool nodeOnPatch(Vector3D& coord) const
    {
        bool pointInCircle = (coord - m_Center).norm() < m_Radius;
        return pointInCircle;
    }
    
    virtual Real activationFunction (const Real& time) const
    {
        Real timeInPeriod = fmod(time - m_tmax + 0.5*m_tduration, 800.);
        bool inPeriod ( timeInPeriod < m_tduration && timeInPeriod > 0);
        Real sinusSquared = std::pow( std::sin(timeInPeriod * PI / m_tduration) , 2 ) * m_patchDisplacement;
        return ( inPeriod ? sinusSquared : 0 );
    }
    
    Real m_patchDisplacement;

    Vector3D m_Center;
    Real m_Radius;
    
    Real m_tmax;
    Real m_tduration;
    
};

REGISTER(EssentialPatchBC, EssentialPatchBCCircularSmooth);

}

#endif /* EssentialPatchBCCircularSmooth_hpp */