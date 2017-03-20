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
#include <lifev/em/solver/circulation/Circulation.hpp>



namespace LifeV
{

template<class meshType>
class HeartSolver {
   
public:
    
    
    HeartSolver(){};


    
    virtual ~HeartSolver() {};
    
    
    void setup(EMSolver<meshType,EMMonodomainSolver<meshType> >& emSolver,  Circulation& circulationSolver)
    {
        M_emSolver = emSolver;
        M_circulationSolver = circulationSolver;
    };
    
    
    
    void preloadHeart(const VectorSmall<2>& endocardiaBC);

    
    
    
    
    
    
protected:
    
    
    EMSolver<meshType,EMMonodomainSolver<meshType> > M_emSolver;
    Circulation M_circulationSolver;
    
    
    VectorSmall<2> M_pressure;
    VectorSmall<2> M_volume;

    
    
    
    
    
    
    Real patchForce (const Real& t, const Real& Tmax, const Real& tmax, const Real& tduration) const
    {
        bool time ( fmod(t-tmax+0.5*tduration, 700.) < tduration && fmod(t-tmax+0.5*tduration, 700.) > 0);
        Real force = std::pow( std::sin(fmod(t-tmax+0.5*tduration, 700.)*3.14159265359/tduration) , 2 ) * Tmax;
        return ( time ? force : 0 );
    }
    
    Real patchFunction (const Real& t, const Real&  X, const Real& Y, const Real& Z, const ID& /*i*/) const
    {
        Real disp = std::pow( std::sin(fmod(t, 700.) * 3.14159265359/300) , 2 )*15;
        return disp;
    }
    
    Real Iapp (const Real& t, const Real&  X, const Real& Y, const Real& Z, const ID& /*i*/) const
    {
        bool coords ( Y < -7. );
        //bool coords ( Y > 4. ); //( Y > 1.5 && Y < 3 );
        bool time ( fmod(t, 700.) < 4 && fmod(t, 700.) > 2);
        return ( coords && time ? 30 : 0 );
    }
    
    Real potentialMultiplyerFcn (const Real& t, const Real&  X, const Real& Y, const Real& Z, const ID& /*i*/) const
    {
        bool time ( fmod(t, 700.) < 4 && fmod(t, 700.) > 2);
        return 1.4 * time; // ( Y < 2.5 && Y > 0.5 ? 1.0 : 0.0 );
    }
    
    Real
    externalPower ( const VectorEpetra& dispCurrent,
                    const VectorEpetra& dispPrevious,
                    const boost::shared_ptr<ETFESpace< RegionMesh<LinearTetra>, MapEpetra, 3, 3 > > dispETFESpace,
                    Real pressure,
                    Real dt,
                    const unsigned int bdFlag) const
    {
        VectorEpetra traction ( dispCurrent.map() );
        VectorEpetra velocity ( (dispCurrent - dispPrevious) / dt );
        
        MatrixSmall<3,3> Id;
        Id(0,0) = 1.; Id(0,1) = 0., Id(0,2) = 0.;
        Id(1,0) = 0.; Id(1,1) = 1., Id(1,2) = 0.;
        Id(2,0) = 0.; Id(2,1) = 0., Id(2,2) = 1.;
        
        {
            using namespace ExpressionAssembly;
            
            auto I = value(Id);
            auto Grad_u = grad( dispETFESpace, dispCurrent, 0);
            auto F =  Grad_u + I;
            auto FmT = minusT(F);
            auto J = det(F);
            auto p = value(pressure);
            
            QuadratureBoundary myBDQR (buildTetraBDQR (quadRuleTria7pt) );
            
            integrate ( boundary ( dispETFESpace->mesh(), bdFlag),
                       myBDQR,
                       dispETFESpace,
                       p * J * dot( FmT * Nface,  phi_i)
                       //p * J * dot( FmT * Nface,  phi_i)
                       //value(-1.0) * J * dot (vE1, FmT * Nface) * phi_i) >> intergral
                       ) >> traction;
            
            traction.globalAssemble();
        }
        
        return traction.dot(velocity);
    }
    
    
};

}

#endif
