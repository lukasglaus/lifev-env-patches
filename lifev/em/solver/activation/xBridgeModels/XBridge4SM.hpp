#ifndef XBRIDGE4SM_HPP_
#define XBRIDGE4SM_HPP_

#include <vector>

namespace LifeV
{

class XBridge4SM
{
public:
    
    typedef double Real;

    // Empty constructor.
    XBridge4SM ();

    // Return number of x-bridge model variables.
    Real size () { return M_NVar; }
    
    // Return initialized x-bridge model state vector.
    std::vector<Real> xbVarVec ();
    
    // Initialize vector with initial values.
    void initialize( std::vector<Real>& variables );
    
    // Compute the rhs of the x-bridge model
    void computeRHS ( const std::vector<Real>& variables, Real Ca, std::vector<Real>& RHS );

    // Solve system with forward euler.
    void solveFE( std::vector<Real>& variables, Real Ca, Real dt );
    

private:
    
    // Variables
   std::vector<Real> M_variables;
    
    // Constants
    Real M_NVar;
    Real M_K2;
    Real M_K3;
    Real M_K4;
    Real M_Kd;
    Real M_Kdd;
    Real M_a1;
    Real M_b1;
    Real M_aa;
    Real M_ba;

}; // XBridge4SM

}

#endif
