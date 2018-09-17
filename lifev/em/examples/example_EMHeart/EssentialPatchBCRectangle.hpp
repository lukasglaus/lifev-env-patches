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
        m_height= dataFile ( ("solid/boundary_conditions/" + m_Name + "/height").c_str(), 1.0 );
        m_width= dataFile ( ("solid/boundary_conditions/" + m_Name + "/width").c_str(), 1.0 );
        
        m_tmax = dataFile ( "solid/patches/tmax", 0. );
        m_tduration = dataFile ( "solid/patches/tduration", 0. );
    }
    
protected:
    
    virtual const bool nodeOnPatch(const Vector3D& coord) const
    {
        Vector3D center;

        Vector3D direction;
        direction(0) = std::cos(m_angle * PI /180);
        direction(1) = 0.0;
        direction(2) = std::sin(m_angle * PI /180);
        direction.normalize();
        
        
        Real normalDistance = ( (coord - center).cross(coord - direction) ).norm();
        
        return (normalDistance < 1.5);
    }
    
    
    Real m_angle;
    Real m_height;
    Real m_width;

};

REGISTER(EssentialPatchBC, EssentialPatchBCRectangle);

}

#endif /* EssentialPatchBCRectangle_hpp */
