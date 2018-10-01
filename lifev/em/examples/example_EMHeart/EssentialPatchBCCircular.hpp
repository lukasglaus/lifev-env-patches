//
//  EssentialPatchBCCircular.hpp
//  lifev-heart
//
//  Created by Thomas Kummer on 30.04.18.
//  Copyright © 2018 Thomas Kummer. All rights reserved.
//

#ifndef EssentialPatchBCCircular_hpp
#define EssentialPatchBCCircular_hpp

#include <stdio.h>
#include <lifev/em/examples/example_EMHeart/EssentialPatchBC.hpp>

#define PI 3.14159265359

namespace LifeV
{
    
class EssentialPatchBCCircular : public EssentialPatchBC
{
public:
    
    EssentialPatchBCCircular(){}
    ~EssentialPatchBCCircular(){}
    
    virtual void setup(const GetPot& dataFile, const std::string& name)
    {
        super::setup(dataFile, name);

        m_Radius= dataFile ( ("solid/boundary_conditions/" + m_Name + "/radius").c_str(), 1.0 );
        
        for ( UInt j (0); j < 3; ++j )
        {
            m_Center[j] = dataFile ( ("solid/boundary_conditions/" + m_Name + "/center").c_str(), 0, j );
        }
        
        m_tmax = dataFile ( "solid/patches/tmax", 0. );
        m_tduration = dataFile ( "solid/patches/tduration", 0. );
    }
    
protected:
    
    virtual const bool nodeOnPatch(const Vector3D& coord) const
    {
        bool pointInCircle = (coord - m_Center).norm() < m_Radius;
        //std::cout << (coord - m_Center).norm() << "," << coord(0) << "," << coord(1) << "," << coord(2) << "; ";
        return pointInCircle;
    }
    
    
    Vector3D m_Center;
    Real m_Radius;
    
};

REGISTER(EssentialPatchBC, EssentialPatchBCCircular);

}

#endif /* EssentialPatchBCCircular_hpp */
