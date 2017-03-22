//
//  HeartData.hpp
//  LifeV
//
//  Created by Thomas Kummer on 22.03.17.
//  Copyright © 2017 Thomas Kummer. All rights reserved.
//

#ifndef HeartData_hpp
#define HeartData_hpp

#include <stdio.h>


namespace LifeV
{
    
class HeartData
{
public:
    
    HeartData() {}
    
    virtual ~HeartData() {};
    
    void setup(const GetPot& datafile)
    {
        setGetPotFile(datafile);
        readFile();
    }
    
    const Real& dt_activation () const { return M_dt_activation; }
    const Real& dt_loadstep () const { return M_dt_loadstep; }
    const Real& activationLimit_loadstep () const { return M_activationLimit_loadstep; }
    const Real& dt_mechanics () const { return M_dt_mechanics; }
    const Real& dt_save () const { return M_dt_save; }
    const Real& endtime () const { return M_endtime; }
    const UInt& mechanicsLoadstepIter () const { return M_mechanicsLoadstepIter; }
    const UInt& mechanicsCouplingIter () const { return M_mechanicsCouplingIter; }
    const UInt& maxiter () const { return M_maxiter; }
    const UInt& preloadSteps () const { return M_preloadSteps; }
    const bool& safePreload () const { return M_safePreload; }
    
    const UInt& pPerturbationFe () const { return M_pPerturbationFe; }
    const UInt& pPerturbationCirc () const { return M_pPerturbationCirc; }
    const UInt& couplingError () const { return M_couplingError; }
    const UInt& couplingJFeSubIter () const { return M_couplingJFeSubIter; }
    const UInt& couplingJFeSubStart () const { return M_couplingJFeSubStart; }
    const UInt& couplingJFeIter () const { return M_couplingJFeIter; }
    
    const GetPot& datafile () { return M_datafile; }
    
protected:
    
    
    void setGetPotFile(const GetPot& datafile)
    {
        M_datafile = datafile;
    }
    
    void readFile()
    {
        M_dt_activation = M_datafile ("activation/time_discretization/timestep", 0.05 );
        M_dt_loadstep =  M_datafile ( "solid/time_discretization/dt_loadstep", 1.0 );
        M_activationLimit_loadstep =  M_datafile ( "solid/time_discretization/activation_limit_loadstep", 0.0 );
        M_dt_mechanics = M_datafile ("solid/time_discretization/timestep", 1.0 );
        M_dt_save = M_datafile ( "exporter/save", 10. );
        M_endtime = M_datafile ("solid/time_discretization/endtime", 100000);
        M_mechanicsLoadstepIter = static_cast<UInt>( M_dt_loadstep / M_dt_activation );
        M_mechanicsCouplingIter = static_cast<UInt>( M_dt_mechanics / M_dt_activation );
        M_maxiter = static_cast<UInt>( M_endtime / M_dt_activation ) ;
        M_preloadSteps = M_datafile ( "solid/boundary_conditions/numPreloadSteps", 0);
        M_safePreload = M_datafile ( "exporter/savePreload", false );
        
        M_pPerturbationFe = M_datafile ( "solid/coupling/pPerturbationFe", 1e-2 );
        M_pPerturbationCirc = M_datafile ( "solid/coupling/pPerturbationCirc", 1e-3 );
        M_couplingError = M_datafile ( "solid/coupling/couplingError", 1e-6 );
        M_couplingJFeSubIter = M_datafile ( "solid/coupling/couplingJFeSubIter", 1 );
        M_couplingJFeSubStart = M_datafile ( "solid/coupling/couplingJFeSubStart", 1 );
        M_couplingJFeIter = M_datafile ( "solid/coupling/couplingJFeIter", 1 );
        
        //        const Real dpMax = dataFile ( "solid/coupling/dpMax", 0.1 );
        //
        //        std::vector<std::vector<std::string> > bcNames { { "lv" , "p" } , { "rv" , "p" } };
        //        std::vector<double> bcValues { p ( "lv" ) , p ( "rv") };
        //        std::vector<double> bcValuesPre ( bcValues );
        //
        //        VectorSmall<4> ABdplv, ABdprv, ABcoef;
        //        ABcoef (0) = 55/24; ABcoef (1) = -59/24; ABcoef (2) = 37/24; ABcoef (3) = -3/8;
        //
        //        VectorSmall<2> VCirc, VCircNew, VCircPert, VFe, VFeNew, VFePert, R, dp;
        //        MatrixSmall<2,2> JFe, JCirc, JR;
        //
        //        UInt iter (0);
        //        Real t (0);
        
        
    }
    
    Real M_dt_activation;
    Real M_dt_loadstep;
    Real M_activationLimit_loadstep;
    Real M_dt_mechanics;
    Real M_dt_save;
    Real M_endtime;
    UInt M_mechanicsLoadstepIter;
    UInt M_mechanicsCouplingIter;
    UInt M_maxiter;
    UInt M_preloadSteps;
    bool M_safePreload;
    
    Real M_pPerturbationFe;
    Real M_pPerturbationCirc;
    Real M_couplingError;
    UInt M_couplingJFeSubIter;
    UInt M_couplingJFeSubStart;
    UInt M_couplingJFeIter;
    
    GetPot M_datafile;
    
};
    
}

#endif /* HeartData_hpp */

