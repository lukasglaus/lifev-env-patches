//
//  HeartSolver.cpp
//  LifeV
//
//  Created by Thomas Kummer on 03.05.16.
//  Copyright © 2016 Thomas Kummer. All rights reserved.
//

#ifndef _HEARTDATA_H_
#define _HEARTDATA_H_


#include <stdio.h>
#include <lifev/em/solver/EMSolver.hpp>



namespace LifeV
{


class HeartSolver {
   
public:
    
    HeartSolver(EMSolver& emSolver,  CirculationSolver& circulationSolver):
        M_emSolver          (emSolver),
        M_circulationSolver (circulationSolver)
    {};
    
    virtual ~HeartSolver() {};
    

    
    
    
    
    
    
private:
    
    
    EMSolver M_emSolver;
    CirculationSolver M_circulationSolver;
    
    
    
    
    VectorSmall<2> M_pressure;
    VectorSmall<2> M_volume;

    
};

}

#endif
