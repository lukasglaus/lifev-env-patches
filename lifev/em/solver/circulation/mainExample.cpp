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

// template <typename T> using Vec = std::vector<T>;

int main(int argc, const char * argv[]) {
    
    
    //==========================================//
    // Functions                                //
    //==========================================//
    
    // Flow rate between two vertices
    auto Q = [] (Circulation& c, const std::string& N1, const std::string& N2) { return c.solution ( std::vector<std::string> {N1, N2} ); };
    
    // Heart activation
    auto phi = [] (const double& t, const double& TPhi, const double& TBeat, const double& TStart) { return pow(std::sin(std::fmod(t - TStart, TBeat) * pi / TPhi), 2) * (std::fmod(t - TStart, TBeat) < TPhi ? 1.0 : 0.0) * (t > TStart ? 1.0 : 0.0); };

    // Ventricular pressure as a function of volume
    //auto p = [] (const double& phi, const double& V, const double& P0, const double & kE, const double& Vu, const double& Emax) {return ( phi * Emax * (V - Vu) + (1 - phi) * P0 * (std::exp(kE*V) - 1) ); };
    
    // Ventricular volume as a function of pressure
    auto V = [] (const double& phi, const double& p, const double& Vn, const double& P0, const double & kE, const double& Vu, const double& Emax) {return ( ( p + phi*Emax*Vu - (1-phi)*P0*(std::exp(kE*Vn)*(1-kE*Vn)-1) ) / ( phi*Emax + (1-phi)*P0*kE*std::exp(kE*Vn) ) ); };
    
    
    //==========================================//
    // Set heart activation parameters          //
    //==========================================//
    
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
    
    const int nBeats = (argc > 1 ? atof(argv[1]) : 1);
    const double dt = (argc > 2 ? atof(argv[2]) : 0.0025);
    const double tEnd = nBeats * TBeat;

    
    //==========================================//
    // Initialize circulation object and b.c.   //
    //==========================================//
    
    Circulation bc( "inputfile" );
    
    std::vector<std::vector<std::string> > bcNames ;
    bcNames.push_back({ "lv" , "p" });
    bcNames.push_back({ "rv" , "p" });
    
    std::vector<double> bcValues(2, 5);
    
    std::vector<double> VCirc {111.8, 108.6};
    std::vector<double> VCircNew (VCirc);
    std::vector<double> VFe (VCirc);
    std::vector<double> VFeNew (VCirc);

    
    //==========================================//
    // Export initial conditions                //
    //==========================================//

    bc.exportSolution( "solution.txt" );
    
    
    //==========================================//
    // Time loop                                //
    //==========================================//
    
    for ( unsigned int i (0) ; i < tEnd / dt ; ++i )
    {
        const double t = (i + 1) * dt;
        unsigned int iter (0);
        const double dp (0.00001);

        
        //==========================================//
        // Solve circulation                        //
        //==========================================//
        
        bc.iterate(dt, bcNames, bcValues, iter);
        VCircNew[0] = VCirc[0] + dt * ( Q(bc, "la", "lv") - Q(bc, "lv", "sa") );
        VCircNew[1] = VCirc[1] + dt * ( Q(bc, "ra", "rv") - Q(bc, "rv", "pa") );

        
        //==========================================//
        // Solve fe-model                           //
        //==========================================//
        
        VFeNew[0] = V(phi(t, TPhi, TBeat, TStart), bcValues[0], VFe[0], P0lv, kElv, Vulv, Emaxlv);
        VFeNew[1] = V(phi(t, TPhi, TBeat, TStart), bcValues[1], VFe[1], P0rv, kErv, Vurv, Emaxrv);
        
        if ( bc.coupling().converged(VFeNew, VCircNew, 1e-6) ) continue;
        
        
        //==========================================//
        // Perturb fe-model                         //
        //==========================================//
        
        const std::vector<double> VFeNew0 ( VFeNew );
        
        // Initialize vectors for pressure perturbation
        std::vector<std::vector<double> > pPert { { bcValues[0] + dp , bcValues[1] } , { bcValues[0] , bcValues[1] + dp } };
        std::vector<std::vector<double> > VFEPert (2, std::vector<double> (2));
        
        VFEPert[0][0] = V(phi(t, TPhi, TBeat, TStart), pPert[0][0], VFe[0], P0lv, kElv, Vulv, Emaxlv);
        VFEPert[1][0] = V(phi(t, TPhi, TBeat, TStart), pPert[0][1], VFe[1], P0rv, kErv, Vurv, Emaxrv);
        
        VFEPert[0][1] = V(phi(t, TPhi, TBeat, TStart), pPert[1][0], VFe[0], P0lv, kElv, Vulv, Emaxlv);
        VFEPert[1][1] = V(phi(t, TPhi, TBeat, TStart), pPert[1][1], VFe[1], P0rv, kErv, Vurv, Emaxrv);

        
        //==========================================//
        // Iterate between circulation and fe-model //
        //==========================================//
        
        while ( !bc.coupling().converged(VFeNew, VCircNew, 1e-6) )
        {
            ++iter;
            
            // Initialize vectors for pressure perturbation
            std::vector<std::vector<double> > pPert { { bcValues[0] + dp , bcValues[1] } , { bcValues[0] , bcValues[1] + dp } };
            std::vector<std::vector<double> > VCircPert (2, std::vector<double> (2));
            
            
            //==========================================//
            // Perturb circulation                      //
            //==========================================//

            bc.iterate(dt, bcNames, pPert[0], iter);
            VCircPert[0][0] = VCirc[0] + dt * ( Q(bc, "la", "lv") - Q(bc, "lv", "sa") );
            VCircPert[1][0] = VCirc[1] + dt * ( Q(bc, "ra", "rv") - Q(bc, "rv", "pa") );
            
            bc.iterate(dt, bcNames, pPert[1], iter);
            VCircPert[0][1] = VCirc[0] + dt * ( Q(bc, "la", "lv") - Q(bc, "lv", "sa") );
            VCircPert[1][1] = VCirc[1] + dt * ( Q(bc, "ra", "rv") - Q(bc, "rv", "pa") );

            
            //==========================================//
            // Update pressure                          //
            //==========================================//
            
            std::vector<double> dpPert (2, dp);
            std::vector<double> R (2);
            for ( unsigned int k (0) ; k < R.size() ; ++k ) R[k] = VFeNew[k] - VCircNew[k];
            bc.coupling().updatePressure(bcValues, VFeNew0, VFEPert, dpPert, VCircNew, VCircPert, dpPert, 1.0, R);

            
            //==========================================//
            // Solve circulation                        //
            //==========================================//

            bc.iterate(dt, bcNames, bcValues, iter);
            VCircNew[0] = VCirc[0] + dt * ( Q(bc, "la", "lv") - Q(bc, "lv", "sa") );
            VCircNew[1] = VCirc[1] + dt * ( Q(bc, "ra", "rv") - Q(bc, "rv", "pa") );
            
            
            //==========================================//
            // Solve fe-model                           //
            //==========================================//

            VFeNew[0] = V(phi(t, TPhi, TBeat, TStart), bcValues[0], VFe[0], P0lv, kElv, Vulv, Emaxlv);
            VFeNew[1] = V(phi(t, TPhi, TBeat, TStart), bcValues[1], VFe[1], P0rv, kErv, Vurv, Emaxrv);

        }
        
        VCirc = VCircNew;
        VFe = VFeNew;
 
        
        //==========================================//
        // Export solution                          //
        //==========================================//

        bc.exportSolution( "solution.txt" );
    }

    return 0;
}
