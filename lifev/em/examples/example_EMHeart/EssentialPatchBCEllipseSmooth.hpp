//
//  EssentialPatchBCEllipseSmooth.hpp
//  lifev-heart
//
//  Created by Thomas Kummer on 30.06.18.
//  Copyright © 2018 Thomas Kummer. All rights reserved.
//

#ifndef EssentialPatchBCEllipseSmooth_hpp
#define EssentialPatchBCEllipseSmooth_hpp

#include <stdio.h>
#include <lifev/em/examples/example_EMHeart/EssentialPatchBC.hpp>

#define PI 3.14159265359

namespace LifeV
{
    
class EssentialPatchBCEllipseSmooth : public EssentialPatchBC
{
public:
    
    EssentialPatchBCEllipseSmooth(){}
    ~EssentialPatchBCEllipseSmooth(){}
    
    virtual void setup(const GetPot& dataFile, const std::string& name)
    {        
        super::setup(dataFile, name);
        
        for ( UInt j (0); j < 3; ++j )
        {
            m_ellipsoidPrincSemiAxesLen[j] = dataFile ( ("solid/boundary_conditions/" + m_Name + "/princSemiAxesLength").c_str(), 0, j );
        }
        
        m_EdgeDispFactor = dataFile ( ("solid/boundary_conditions/" + m_Name + "/edgeDispFactor").c_str(), 0 );

        for ( UInt j (0); j < 3; ++j )
        {
            m_Center[j] = dataFile ( ("solid/boundary_conditions/" + m_Name + "/center").c_str(), 0, j );
        }
    }
    
protected:
    
    virtual vectorPtr_Type directionalVectorField (const boost::shared_ptr<FESpace<RegionMesh<LinearTetra>, MapEpetra >> dFeSpace, Vector3D& direction, const Real& disp, const Real& time) const
    {
        // New P1 Space
        FESpace<RegionMesh<LinearTetra> , MapEpetra > p1FESpace ( dFeSpace->mesh(), "P1", 3, dFeSpace->mesh()->comm() );
        
        // Create P1 VectorEpetra
        VectorEpetra p1PositionVector (p1FESpace.map());
        
        // Fill P1 vector with mesh values
        Int p1nCompLocalDof = p1PositionVector.epetraVector().MyLength() / 3;
        for (int j (0); j < p1nCompLocalDof; j++)
        {
            UInt iGID = p1PositionVector.blockMap().GID (j);
            
            p1PositionVector[iGID] = dFeSpace->mesh()->point (iGID).x();
            p1PositionVector[jGID] = dFeSpace->mesh()->point (iGID).y();
            p1PositionVector[kGID] = dFeSpace->mesh()->point (iGID).z();
        }
        
        // Interpolate position vector from P1-space to current space
        VectorEpetra positionVector ( disp.map() );
        positionVector = dFeSpace->feToFEInterpolate(p1FESpace, p1PositionVector);
        positionVector += disp;
        
        vectorPtr_Type vectorField (new VectorEpetra( dFeSpace->map(), Repeated ));
        auto nCompLocalDof = vectorField->epetraVector().MyLength() / 3;
        
        direction.normalize();

        // Fill P1 vector with mesh values
        Int p1nCompLocalDof = p1PositionVector.epetraVector().MyLength() / 3;
        for (int j (0); j < nCompLocalDof; j++)
        {
            UInt iGID = positionVector->blockMap().GID (j);
            UInt jGID = positionVector->blockMap().GID (j + nCompLocalDof);
            UInt kGID = positionVector->blockMap().GID (j + 2 * nCompLocalDof);
            
            Vector3D coord;
            
            coord(0) = (*positionVector)[iGID];
            coord(1) = (*positionVector)[jGID];
            coord(2) = (*positionVector)[kGID];
            
            // Radial and axial distance to center line
            Vector3D currentPatchCenter = m_Center + activationFunction(time) * direction;
            auto radialDistance = ( (coord - m_Center).cross(coord - currentPatchCenter) ).norm() / (m_Center - currentPatchCenter).norm();
            auto axialDistance = (coord - currentPatchCenter).dot(direction) * direction;
            
            // If coordiantes inside or outside of a certain radius
            auto displacement = disp - disp * (1 - m_EdgeDispFactor) * dispDistributionWeight(coord);
            
            // If patch inside or outside the structure
            
            
            // Scale the direction vector
            auto displacementVec = direction * displacement;
            (*vectorField)[iGID] = displacementVec[0];
            (*vectorField)[jGID] = displacementVec[1];
            (*vectorField)[kGID] = displacementVec[2];
        }
        
        return vectorField;

//        vectorPtr_Type vectorField (new VectorEpetra( dFeSpace->map(), Repeated ));
//        auto nCompLocalDof = vectorField->epetraVector().MyLength() / 3;
//
//        direction.normalize();
//        // direction *= disp;
//
//        for (int j (0); j < nCompLocalDof; ++j)
//        {
//            // Get coordiantes
//            UInt iGID = vectorField->blockMap().GID (j);
//            UInt jGID = vectorField->blockMap().GID (j + nCompLocalDof);
//            UInt kGID = vectorField->blockMap().GID (j + 2 * nCompLocalDof);
//
//            Vector3D coord;
//
//            coord(0) = dFeSpace->mesh()->point(iGID).x() + (*m_dispPtr)[iGID];
//            coord(1) = dFeSpace->mesh()->point(iGID).y() + (*m_dispPtr)[jGID];
//            coord(2) = dFeSpace->mesh()->point(iGID).z() + (*m_dispPtr)[kGID];
//
//            // Radial and axial distance to center line
//            Vector3D currentPatchCenter = m_Center + activationFunction(time) * direction;
//            auto radialDistance = ( (coord - m_Center).cross(coord - currentPatchCenter) ).norm() / (m_Center - currentPatchCenter).norm();
//            auto axialDistance = (coord - currentPatchCenter).dot(direction) * direction;
//
//            // If coordiantes inside or outside of a certain radius
//            auto displacement = disp - disp * (1 - m_EdgeDispFactor) * dispDistributionWeight(coord);
//
//            // If patch inside or outside the structure
//
//
//            // Scale the direction vector
//            auto displacementVec = direction * displacement;
//            (*vectorField)[iGID] = displacementVec[0];
//            (*vectorField)[jGID] = displacementVec[1];
//            (*vectorField)[kGID] = displacementVec[2];
//        }

    }
    
    virtual const bool nodeOnPatch(const Vector3D& coord) const
    {
        auto ellipsoidCoord = transformToLocalEllipsoidCoordinates(coord);
        return nodeInsideEllipsoid(ellipsoidCoord);
    }
    
    virtual const bool nodeInsideEllipsoid(const Vector3D& ellipseCoord) const
    {
        return (ellipsoidFuncEval(ellipseCoord) < 1.0);
    }
    
    virtual const Real dispDistributionWeight(const Vector3D& coord) const
    {
        auto ellipsoidCoord = transformToLocalEllipsoidCoordinates(coord);
        return ( nodeInsideEllipsoid(ellipsoidCoord) ? ellipsoidFuncEval(ellipsoidCoord) : std::pow(1 - m_EdgeDispFactor, -1.0) );
    }
    
    virtual const Real ellipsoidFuncEval(const Vector3D& ellipseCoord) const
    {
        Real ellipsoidFuncValue = std::pow(ellipseCoord(0) / m_ellipsoidPrincSemiAxesLen(0), 2.0) + std::pow(ellipseCoord(1) / m_ellipsoidPrincSemiAxesLen(1), 2.0) + std::pow(ellipseCoord(2) / m_ellipsoidPrincSemiAxesLen(2), 2.0);
        return ellipsoidFuncValue;
    }

    virtual const std::vector<Vector3D> ellipsoidCoordinateSystem(Vector3D patchDirection) const
    {
        auto axis0 = patchDirection.normalized();
        auto axis1 = (Vector3D( 1.0 , 0.0 , - axis0(0) / axis0(2))).normalized();
        auto axis2 = (axis0.cross(axis1)).normalized();
        
        return std::vector<Vector3D> { axis0 , axis1 , axis2 };
    }
    
    virtual const Vector3D transformToLocalEllipsoidCoordinates(const Vector3D& coord) const
    {
        auto ellipsoidCS = ellipsoidCoordinateSystem(m_patchDirection);
        auto localCoord = coord - m_Center;
        Vector3D localEllipsoidCoord( ellipsoidCS[0].dot(localCoord) , ellipsoidCS[1].dot(localCoord) , ellipsoidCS[2].dot(localCoord) );
        return localEllipsoidCoord;
    }
        

    Vector3D m_Center;
    Vector3D m_ellipsoidPrincSemiAxesLen;
    Real m_EdgeDispFactor;
    
};

REGISTER(EssentialPatchBC, EssentialPatchBCEllipseSmooth);

}

#endif /* EssentialPatchBCEllipseSmooth_hpp */
