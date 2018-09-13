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
        // auto p2PositionVector = p2PositionVectorDisplaced(dFeSpace);
        auto p2PositionVector = p2PositionVectorInitial(dFeSpace);

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
            auto currentPatchCenter = m_Center + activationFunction(time) * direction;
            auto radialDistance = ( (coord - m_Center).cross(coord - currentPatchCenter) ).norm() / (m_Center - currentPatchCenter).norm();
            auto axialDistance = (coord - currentPatchCenter).dot(direction) * direction;
            
            // Determine the patch displacement as a function of patch coordinates
            auto displacement = disp - disp * (1 - m_EdgeDispFactor) * dispDistributionWeight(coord);
            
            // Scale the direction vector
            auto displacementVec = direction * displacement;
            (*p2PatchDisplacement)[iGID] = displacementVec[0];
            (*p2PatchDisplacement)[jGID] = displacementVec[1];
            (*p2PatchDisplacement)[kGID] = displacementVec[2];
        }
        
        return p2PatchDisplacement;

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
