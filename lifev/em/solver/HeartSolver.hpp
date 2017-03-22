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

    
class HeartData
{
public:
    
    HeartData() {}
    
    virtual ~HeartData() {};
    
    const Real& dt_activation () const { return M_dt_activation; }
    const Real& dt_loadstep () const { return M_dt_loadstep; }
    const Real& activationLimit_loadstep () const { return M_activationLimit_loadstep; }
    
    
    const GetPot& datafile () { return M_datafile; }
    
protected:
    
    
    
    void setupData(GetPot& datafile)
    {
        M_datafile = datafile;
        
        M_dt_activation = M_datafile ("activation/time_discretization/timestep", 0.05 );
        M_dt_loadstep =  M_datafile ( "solid/time_discretization/dt_loadstep", 1.0 );
        M_activationLimit_loadstep =  M_datafile ( "solid/time_discretization/activation_limit_loadstep", 0.0 );
        M_dt_mechanics = M_datafile ("solid/time_discretization/timestep", 1.0 );
        M_dt_save = M_datafile ( "exporter/save", 10. );
        M_endtime = M_datafile ("solid/time_discretization/endtime", 100000);
        M_mechanicsLoadstepIter = static_cast<UInt>( M_dt_loadstep / M_dt_activation );
        M_mechanicsCouplingIter = static_cast<UInt>( M_dt_mechanics / M_dt_activation );
        M_maxiter = static_cast<UInt>( M_endtime / M_dt_activation ) ;

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
    
    Real M_pPerturbationFe;
    Real M_pPerturbationCirc;
    Real M_couplingError;
    UInt M_couplingJFeSubIter;
    UInt M_couplingJFeSubStart;
    UInt M_couplingJFeIter;
    
    GetPot M_datafile;
    
};
    
    
template <class EmSolver>
class HeartSolver {
   
public:
    
    HeartSolver(EmSolver& emSolver,  Circulation& circulationSolver) :
        M_emSolver          (emSolver),
        M_circulationSolver (circulationSolver),
        M_heartData         (HeartData())
    {}
    
    virtual ~HeartSolver() {}
    
    const HeartData& data()
    {
        return M_heartData;
    }
    
    void preloadHeart(const VectorSmall<2>& endocardiaBC);

    void setupData(const GetPot& datafile)
    {
        M_heartData.setupData(datafile);
    }
    
    
    
    
    
    
protected:
    
    
    EmSolver M_emSolver;
    Circulation M_circulationSolver;
    
    HeartData M_heartData;
    
    
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
