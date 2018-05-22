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
        m_Name = name;
        m_PrevFlag = dataFile ( ("solid/boundary_conditions/" + m_Name + "/flag").c_str(), 0 );
        m_patchDisplacement = dataFile ( ("solid/boundary_conditions/" + m_Name + "/displacement").c_str(), 1.0 );
        m_Radius= dataFile ( ("solid/boundary_conditions/" + m_Name + "/radius").c_str(), 1.0 );
        
        for ( UInt j (0); j < 3; ++j )
        {
            m_Center[j] = dataFile ( ("solid/boundary_conditions/" + m_Name + "/center").c_str(), 0, j );
        }
        
        Real m_tmax = dataFile ( "solid/patches/tmax", 0. );
        Real m_tduration = dataFile ( "solid/patches/tduration", 0. );
    }
    
protected:
    
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
        std::cout << time << " " << timeInPeriod << " " << inPeriod << " " << sinusSquared << std::endl
        return ( inPeriod ? sinusSquared : 0 );
    }
    
    Real m_patchDisplacement;

    Vector3D m_Center;
    Real m_Radius;
    
    Real m_tmax;
    Real m_tduration;
    
};

REGISTER(EssentialPatchBC, EssentialPatchBCCircular);

}

#endif /* EssentialPatchBCCircular_hpp */
