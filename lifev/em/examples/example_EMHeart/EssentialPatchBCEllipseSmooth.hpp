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
        auto p2PositionVector = p2PositionVectorDisplaced(dFeSpace);
        
        vectorPtr_Type p2PatchDisplacement (new VectorEpetra( dFeSpace->map(), Repeated ));
        auto nCompLocalDof = p2PatchDisplacement->epetraVector().MyLength() / 3;
        
        direction.normalize();

        for (int j (0); j < nCompLocalDof; j++)
        {
            UInt iGID = p2PatchDisplacement->blockMap().GID (j);
            UInt jGID = p2PatchDisplacement->blockMap().GID (j + nCompLocalDof);
            UInt kGID = p2PatchDisplacement->blockMap().GID (j + 2 * nCompLocalDof);
            
            Vector3D coord;
            coord(0) = p2PositionVector[iGID];
            coord(1) = p2PositionVector[jGID];
            coord(2) = p2PositionVector[kGID];
            
            // Radial and axial distance to current patch center and center line
            Vector3D currentPatchCenter = m_Center + activationFunction(time) * direction;
            auto radialDistance = ( (coord - m_Center).cross(coord - currentPatchCenter) ).norm() / (m_Center - currentPatchCenter).norm();
            auto axialDistance = (coord - currentPatchCenter).dot(direction) * direction;
            
            // Determine the patch displacement as a function of patch coordinates
            auto displacement = disp - disp * (1 - m_EdgeDispFactor) * dispDistributionWeight(coord);
            
            // If patch inside or outside the structure
            
            
            // Scale the direction vector
            auto displacementVec = direction * displacement;
            (*p2PatchDisplacement)[iGID] = displacementVec[0];
            (*p2PatchDisplacement)[jGID] = displacementVec[1];
            (*p2PatchDisplacement)[kGID] = displacementVec[2];
        }
        
        return p2PatchDisplacement;

    }
    
    
    virtual vector_Type p2PositionVectorDisplaced(const boost::shared_ptr<FESpace<RegionMesh<LinearTetra>, MapEpetra >> p2dFeSpace) const
    {
        // New P1 Space
        FESpace<RegionMesh<LinearTetra> , MapEpetra > p1dFESpace ( p2dFeSpace->mesh(), "P1", 3, dFeSpace->mesh()->comm() );
        
        // Create P1 VectorEpetra
        VectorEpetra p1PositionVector (p1dFESpace.map());
        
        // Fill P1 vector with mesh values
        Int p1nCompLocalDof = p1PositionVector.epetraVector().MyLength() / 3;
        for (int j (0); j < p1nCompLocalDof; j++)
        {
            UInt iGID = p1PositionVector.blockMap().GID (j);
            UInt jGID = p1PositionVector.blockMap().GID (j + p1nCompLocalDof);
            UInt kGID = p1PositionVector.blockMap().GID (j + 2 * p1nCompLocalDof);
            
            p1PositionVector[iGID] = p2dFeSpace->mesh()->point (iGID).x();
            p1PositionVector[jGID] = p2dFeSpace->mesh()->point (iGID).y();
            p1PositionVector[kGID] = p2dFeSpace->mesh()->point (iGID).z();
        }
        
        // Interpolate position vector from P1-space to current space
        VectorEpetra p2PositionVector ( m_dispPtr->map() );
        p2PositionVector = p2dFeSpace->feToFEInterpolate(p1dFESpace, p1PositionVector);
        p2PositionVector += (*m_dispPtr);
        
        return p2PositionVector;
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
