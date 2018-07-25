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
        vectorPtr_Type vectorField (new VectorEpetra( dFeSpace->map(), Repeated ));
        auto nCompLocalDof = vectorField->epetraVector().MyLength() / 3;

        direction.normalize();
        direction *= disp;

        for (int j (0); j < nCompLocalDof; ++j)
        {
            // Get coordiantes
            UInt iGID = vectorField->blockMap().GID (j);
            UInt jGID = vectorField->blockMap().GID (j + nCompLocalDof);
            UInt kGID = vectorField->blockMap().GID (j + 2 * nCompLocalDof);

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
//            auto displacement = (m_EdgeDispFactor * disp - disp) * dispDistributionWeight(coord) + disp;

            // If patch inside or outside the structure


            // Scale the direction vector
            auto displacementVec = direction; //* displacement;
            (*vectorField)[iGID] = displacementVec[0];
            (*vectorField)[jGID] = displacementVec[1];
            (*vectorField)[kGID] = displacementVec[2];
        }

        return vectorField;
    }
    
    virtual const bool nodeOnPatch(Vector3D& coord) const
    {
        auto ellipsoidCS = ellipsoidCoordinateSystem(m_patchDirection);
        auto localCoord = coord - m_Center;
        Vector3D ellipsoidCoord( ellipsoidCS[0].dot(localCoord) , ellipsoidCS[1].dot(localCoord) , ellipsoidCS[2].dot(localCoord) );
        
        return nodeInsideEllipsoid(ellipsoidCoord);
    }
    
    virtual const std::vector<Vector3D> ellipsoidCoordinateSystem(Vector3D patchDirection) const
    {
        auto axis0 = patchDirection.normalized();
        auto axis1 = (Vector3D( 1.0 , 0.0 , - axis0(0) / axis0(2))).normalized();
        auto axis2 = (axis0.cross(axis1)).normalized();

        return std::vector<Vector3D> { axis0 , axis1 , axis2 };
    }
    
    virtual const bool nodeInsideEllipsoid(const Vector3D& ellipseCoord) const
    {
        bool pointInsideEllipsoid (ellipsoidFuncEval(ellipseCoord) < 1.0);
        return pointInsideEllipsoid;
    }
    
    virtual const Real ellipsoidFuncEval(const Vector3D& ellipseCoord) const
    {
        Real ellipsoidFuncValue = std::pow(ellipseCoord(0) / m_ellipsoidPrincSemiAxesLen(0), 2.0) + std::pow(ellipseCoord(1) / m_ellipsoidPrincSemiAxesLen(1), 2.0) + std::pow(ellipseCoord(2) / m_ellipsoidPrincSemiAxesLen(2), 2.0);
        return ellipsoidFuncValue;
    }
    
    virtual const Real dispDistributionWeight(Vector3D& coord) const
    {
        auto ellipsoidCS = ellipsoidCoordinateSystem(m_patchDirection);
        auto localCoord = coord - m_Center;
        Vector3D ellipsoidCoord( ellipsoidCS[0].dot(localCoord) , ellipsoidCS[1].dot(localCoord) , ellipsoidCS[2].dot(localCoord) );
        
        Real dispWeight (0);
        if ( nodeInsideEllipsoid(ellipsoidCoord) )
        {
            dispWeight =  = ellipsoidFuncEval(ellipsoidCoord);
        }
        
        return dispWeight;
    }
        

    Vector3D m_Center;
    Vector3D m_ellipsoidPrincSemiAxesLen;
    Real m_EdgeDispFactor;
    
};

REGISTER(EssentialPatchBC, EssentialPatchBCEllipseSmooth);

}

#endif /* EssentialPatchBCEllipseSmooth_hpp */
