//
//  main.cpp
//  Circulation
//
//  Created by Thomas Kummer on 10.04.15.
//  Copyright (c) 2015 Thomas Kummer. All rights reserved.
//

#include <iostream>
#include <math.h>

#include "Circulation.hpp"

#define pi 3.14159265

template <typename T> using Vec = std::vector<T>;
typedef std::vector<std::string> VecS;

int main(int argc, const char * argv[]) {
    
    auto Q = [] (Circulation& c, const std::string& N1, const std::string& N2) { return c.solution ( VecS {N1, N2} ); };
    
    auto pCirc = [] (Circulation& c, const std::string& N1) { return c.solution ( N1 ); };
    
    // Heart activation
    auto phi = [] (const double& t, const double& TPhi, const double& TBeat, const double& TStart) { return pow(std::sin(std::fmod(t - TStart, TBeat) * pi / TPhi), 2) * (std::fmod(t - TStart, TBeat) < TPhi ? 1.0 : 0.0) * (t > TStart ? 1.0 : 0.0); };

    //auto pFe = [] (const double& phi, const double& V, const double& P0, const double & kE, const double& Vu, const double& Emax) {return ( phi * Emax * (V - Vu) + (1 - phi) * P0 * (std::exp(kE*V) - 1) ); };
    
    auto V = [] (const double& phi, const double& p, const double& Vn, const double& P0, const double & kE, const double& Vu, const double& Emax) {return ( ( p + phi*Emax*Vu - (1-phi)*P0*(std::exp(kE*Vn)*(1-kE*Vn)-1) ) / ( phi*Emax + (1-phi)*P0*kE*std::exp(kE*Vn) ) ); };
    
    const double P0lv = 1.5;
    const double P0rv = 1.5;
    const double kElv = 0.8*0.014;
    const double kErv = 0.8*0.011;
    const double Vulv = 16.77;
    const double Vurv = 40.8;
    const double Emaxlv = 2.95;
    const double Emaxrv = 1.75;
    
    const double TStart = 0.2;
    const double TPhi = 0.41;
    const double TBeat = 0.833;
    
    const double dt = 0.0025;
    const int nBeats = (argv[1] ? atof(argv[1]) : 1);
    const double tEnd = nBeats * TBeat;

    // Body circulation
    Circulation bc( "inputfile" );
    
    std::vector<std::vector<std::string> > bcNames ;
    bcNames.push_back({ "lv" , "Q" });
    bcNames.push_back({ "rv" , "p" });
    
    std::vector<double> bcValues { 0 , 5 };
    
    std::vector<double> VCirc {111.8, 108.6};
    std::vector<double> VCircNew (VCirc);
    std::vector<double> VFe (VCirc);
    std::vector<double> VFeNew (VCirc);

    bc.exportSolution( "solution.txt" );
    
    
    // Time loop
    for ( unsigned int i (0) ; i < tEnd / dt ; ++i )
    {
        const double t = (i + 1) * dt;
        unsigned int iter (0);
        const double dp (0.00001);

        // Solve circulation
        bc.iterate(dt, bcNames, bcValues, iter);
        const double& plvCirc = pCirc(bc, "lv");
        VCircNew[1] = VCirc[1] + dt * ( Q(bc, "ra", "rv") - Q(bc, "rv", "pa") + Q(bc, "rv", "-") );
        
        // Solve fe-model
        VFeNew[0] = V(phi(t, TPhi, TBeat, TStart), plvCirc, VFe[0], P0lv, kElv, Vulv, Emaxlv);
        VFeNew[1] = V(phi(t, TPhi, TBeat, TStart), bcValues[1], VFe[1], P0rv, kErv, Vurv, Emaxrv);
    
        const std::vector<double> VFeNew0 ( VFeNew );
        
        // Initialize vectors for pressure perturbation
        std::vector<std::vector<double> > pPert { { bcValues[0] , bcValues[1] + dp } };
        std::vector<std::vector<double> > VFEPert (2, std::vector<double> (2));
        
        // Compute VFe (p+dp)
        VFEPert[1][1] = V(phi(t, TPhi, TBeat, TStart), pPert[0][1], VFe[1], P0rv, kErv, Vurv, Emaxrv);
        
        // Iterate between circulation and fe-model
        while ( !bc.coupling().converged(std::vector<double> {VFeNew[1]}, std::vector<double>{VCircNew[1]}, 1e-6) )
        {
            ++iter;
            
            // Initialize vectors for pressure perturbation
            std::vector<std::vector<double> > pPert { { bcValues[0] , bcValues[1] + dp } };
            std::vector<std::vector<double> > VCircPert (2, std::vector<double> (2));
            
            
            // Compute VCirc (p+dp)
            bc.iterate(dt, bcNames, pPert[0], iter);
            VCircPert[1][1] = VCirc[1] + dt * ( Q(bc, "ra", "rv") - Q(bc, "rv", "pa") + Q(bc, "rv", "-") );

            
            // update pressure
            Vec<double> dpPert (1, dp);
            std::vector<double> R (1);
            std::vector<double> prv {bcValues[1]};
            for ( unsigned int k (0) ; k < R.size() ; ++k ) R[0] = VFeNew[1] - VCircNew[1];
            bc.coupling().updatePressure(prv, std::vector<double> {VFeNew0[1]}, std::vector<std::vector<double> > {{VFEPert[1][1]}}, dpPert, std::vector<double> {VCircNew[1]}, std::vector<std::vector<double> > {{VCircPert[1][1]}}, dpPert, 1.0, R);
            bcValues[1] = prv[0];
            
            // Compute VCirc (p)
            bc.iterate(dt, bcNames, bcValues, iter);
            VCircNew[1] = VCirc[1] + dt * ( Q(bc, "ra", "rv") - Q(bc, "rv", "pa") + Q(bc, "rv", "-") );
            
            
            // ( Compute VFe (p) )
            VFeNew[1] = V(phi(t, TPhi, TBeat, TStart), bcValues[1], VFe[1], P0rv, kErv, Vurv, Emaxrv);

        }
        
        bcValues[0] = - ( VFeNew[0] - VFe[0] ) / dt;
        
        VCirc = VCircNew;
        VFe = VFeNew;
 
        if ( i % 1 == 0 ) bc.exportSolution( "solution.txt" );
    }

    return 0;
}
