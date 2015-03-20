#include <lifev/em/solver/activation/xBridgeModels/XBridge4SM.hpp>
#include <stdexcept>
#include <cmath>
#include <iostream>

namespace LifeV
{

XBridge4SM::XBridge4SM():
    M_NVar  ( 8 ),
    M_K2    ( 127.0 ),
    M_K3    ( 30.0 ),
    M_K4    ( 3.0 ),
    M_Kd    ( 67.0 ),
    M_Kdd   ( 139 ),
    M_a1    ( 1.0 ),
    M_b1    ( 0.3 ),
    M_aa    ( 1.0 ),
    M_ba    ( 3.0 )
{
}

    
std::vector<double>
XBridge4SM::xbVarVec()
{
    std::vector<Real> xbVarVec ( M_NVar );
    initialize ( xbVarVec );
    return xbVarVec;
}
    
    
void
XBridge4SM::initialize( std::vector<Real>& variables )
{
    if ( variables.size() != M_NVar ) throw std::runtime_error("XBridge4SM::initialize: Number of variables not correct.");
    
    variables.at (0) = 0;     // XA
    variables.at (1) = 70;    // TnCA
    variables.at (2) = 20;    // M
    variables.at (3) = 0;     // Ca TnCA
    variables.at (4) = 0;     // Ca TnCA M
    variables.at (5) = 0;     // TnCA M
    variables.at (6) = M_a1 * std::pow( variables.at (4) + variables.at(5) , 0.5 ) + M_b1; // K1
    variables.at (7) = M_aa * std::pow( variables.at (4) + variables.at(5) , 2.0 ) + M_ba; // Ka
}

    
void
XBridge4SM::computeRHS ( const std::vector<Real>& variables, Real Ca, std::vector<Real>& RHS )
{

    auto XA     = variables.at (0);    // XA
    auto TnCA   = variables.at (1);    // TnCA
    auto M      = variables.at (2);    // M
    auto CaTnCA = variables.at (3);    // Ca TnCA
    auto CaTnCAM= variables.at (4);    // Ca TnCA M
    auto TnCAM  = variables.at (5);    // TnCA M
    auto K1     = variables.at (6);    // K1
    auto Ka     = variables.at (7);    // Ka
    
    RHS.at(0) = CaTnCAM + TnCAM;
    RHS.at(1) = - K1 * Ca * TnCA + M_K3 * CaTnCA + M_Kdd * TnCAM;
    RHS.at(2) = M_Kdd * TnCAM - Ka * CaTnCA * M + M_Kd * CaTnCAM;
    RHS.at(3) = K1 * Ca * TnCA - M_K3 * CaTnCA - Ka * CaTnCA * M + M_Kd * CaTnCAM;
    RHS.at(4) = Ka * CaTnCA * M + M_K2 * Ca * TnCAM - (M_Kd + M_K4) * CaTnCAM;
    RHS.at(5) = M_K4 * CaTnCAM - M_Kdd * TnCAM - M_K2 * Ca *TnCAM;
    RHS.at(6) = M_a1 * std::pow( CaTnCAM + TnCAM , 0.5 ) + M_b1;
    RHS.at(7) = M_aa * std::pow( CaTnCAM + TnCAM , 2.0 ) + M_ba;
}

    
void
XBridge4SM::solveFE( std::vector<Real>& variables, Real Ca, Real dt)
{
    std::vector<Real> RHS( variables.size() );
    computeRHS(variables, Ca, RHS);
    
    for ( unsigned int i (0); i < M_NVar; ++i )
    {
        variables.at(i) += dt * RHS.at(i);
    }
}
    
    
}

