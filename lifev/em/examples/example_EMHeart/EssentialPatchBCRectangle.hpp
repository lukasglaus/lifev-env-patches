//
//  EssentialPatchBCRectangle.hpp
//  lifev-heart
//
//  Created by Thomas Kummer on 30.04.18.
//  Copyright © 2018 Thomas Kummer. All rights reserved.
//

#ifndef EssentialPatchBCRectangle_hpp
#define EssentialPatchBCRectangle_hpp

#include <stdio.h>
#include <lifev/em/examples/example_EMHeart/EssentialPatchBC.hpp>

#define PI 3.14159265359

namespace LifeV
{
    
class EssentialPatchBCRectangle : public EssentialPatchBC
{
public:
    
    EssentialPatchBCRectangle(){}
    ~EssentialPatchBCRectangle(){}
    
    virtual void setup(const GetPot& dataFile, const std::string& name)
    {
        super::setup(dataFile, name);

        m_angle= dataFile ( ("solid/boundary_conditions/" + m_Name + "/angle").c_str(), 0.0 );
        m_dAngle = dataFile ( ("solid/boundary_conditions/" + m_Name + "/width").c_str(), 30.0 );
        m_height= dataFile ( ("solid/boundary_conditions/" + m_Name + "/height").c_str(), -7.5 );
        m_width= dataFile ( ("solid/boundary_conditions/" + m_Name + "/width").c_str(), 1.5 );
        
        m_tmax = dataFile ( "solid/patches/tmax", 0. );
        m_tduration = dataFile ( "solid/patches/tduration", 0. );
    }
    
protected:
    
    virtual const bool nodeOnPatch(const Vector3D& coord) const
    {
        Vector3D coordZyl;
        coordZyl(0) = std::sqrt(std::pow(coord(0), 2) + std::pow(coord(2), 2)); // radius
        coordZyl(1) = coord(1); // height
        coordZyl(2) = std::atan2(coord(2), coord(1)); // angle

        angle *= m_angle * PI/180;
        dAngle *= m_dAngle * PI/180;

        const bool inAngleRange ( coordZyl(2) < (angle + dAngle) && coordZyl(2) > (angle - dAngle) );
        const bool inVerticalRange ( coordZyl(1) < (m_height + m_width) && coordZyl(1) > (m_height - m_width) );
        
        return (inAngleRange && inVerticalRange);
    }
    
    
    Real m_angle;
    Real m_dAngle;
    Real m_height;
    Real m_width;

};

REGISTER(EssentialPatchBC, EssentialPatchBCRectangle);

}

#endif /* EssentialPatchBCRectangle_hpp */
