//
//  CirculationCoupling.hpp
//  Circulation
//
//  Created by Thomas Kummer on 03.05.15.
//  Copyright (c) 2015 Thomas Kummer. All rights reserved.
//

#ifndef CIRCULATIONCOUPLING_H_
#define CIRCULATIONCOUPLING_H_

// #include <stdio.h>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <vector>

class CirculationCoupling {
public:
    
    typedef std::vector<double> VectorStdDouble;
    typedef std::vector<VectorStdDouble> MatrixStdDouble;
    typedef Eigen::MatrixXd MatrixEigen;
    typedef Eigen::VectorXd VectorEigen;
    
    CirculationCoupling(){}
    
    virtual ~CirculationCoupling(){}
    
    bool converged(const VectorStdDouble& V1,  const VectorStdDouble& V2, const double& e) const
    {
        VectorEigen R = stdToEigen(V1) - stdToEigen(V2);
        if ( R.norm() < e ) return true;
        else return false;
    }
    
    void updatePressure(VectorStdDouble& p, const VectorStdDouble& V1, const MatrixStdDouble& Vp1, const VectorStdDouble& dp1, const VectorStdDouble& V2, const MatrixStdDouble& Vp2, const VectorStdDouble& dp2, const double& relaxationParam = 1, const VectorStdDouble& Rstd = VectorStdDouble (0)) const
    {
        // Compute Residual
        VectorEigen R = stdToEigen(Rstd);
        //VectorEigen R = stdToEigen(V1) - stdToEigen(V2);

        // Compute Jacobian
        auto J1 = assembleJacobian(V1, Vp1, dp1);
        auto J2 = assembleJacobian(V2, Vp2, dp2);
        auto JR = J1 - J2;
        
        // Update p
        p = solveNewtonStep(p, R, JR, relaxationParam);
    }
    
    const VectorStdDouble solveNewtonStep(const VectorStdDouble& p, const VectorEigen& R, const MatrixEigen& J, const double& relaxationParam = 1) const
    {
        VectorEigen rhs = - relaxationParam * R;
        VectorEigen dp = J.fullPivLu().solve(rhs);
        VectorEigen pNew =  stdToEigen( p ) + dp;
        return eigenToStd( pNew );
    }
    
    const MatrixEigen assembleJacobian(const VectorStdDouble& V, const MatrixStdDouble& Vp, const VectorStdDouble& dp) const
    {
        MatrixEigen J( dp.size() , dp.size() );
        
        for (unsigned int i (0) ; i < dp.size() ; ++i)
        {
            for (unsigned int j (0) ; j < dp.size() ; ++j)
            {
                J(i, j) = ( Vp [i][j] - V[i] ) / dp[j];
            }
        }
        
        return J;
    }
    
    const VectorEigen residual(const VectorEigen& V, const VectorEigen& W) const
    {
        return ( V - W );
    }

    VectorStdDouble eigenToStd(const VectorEigen& u) const
    {
        VectorStdDouble u_( u.size() );
        VectorEigen::Map( &u_[0] , u.size() ) = u;
        return u_;
    }
    
    VectorEigen stdToEigen(const VectorStdDouble& u) const
    {
        VectorEigen u_ = VectorEigen::Map(u.data(), u.size());
        return u_;
    }

    
private:
    
};

#endif /* CIRCULATIONCOUPLING_H_ */