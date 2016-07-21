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

namespace LifeV
{

class HeartSolver {
   
public:
    
    HeartSolver(FESolver& feSolver,  CirculationSolver& circulationSolver):
        M_feSolver          (feSolver),
        M_circulationSolver (circulationSolver)
    {};
    
    virtual ~HeartSolver() {};
    
    // get data void setup(FESolver& feSolver)
    
    void extrapolatePressure(const double& dt)
    {
        
    }
    
    template<class FESolver, class CirculationSolver>
    void computeJacobian(FESolver& feSolver,  CirculationSolver& circulationSolver)
    {}
    
    void solveCirculation()
    {}
    
    
private:
    
    VectorSmall<2> M_pressure;
    VectorSmall<2> M_volume;

    
};

}

#endif