//
//  CirculationTimeIntegrator.cpp
//  Circulation
//
//  Created by Thomas Kummer on 16.05.15.
//  Copyright (c) 2015 Thomas Kummer. All rights reserved.
//

#include <stdio.h>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <vector>
#include <string>

#include "CirculationAssembler.hpp"


class TimeIntegrator {
public:

    TimeIntegrator() {}
    
    virtual ~TimeIntegrator() {}
    
    virtual Eigen::MatrixXd assembleOperator(const double& dt, Eigen::MatrixXd& M, Eigen::MatrixXd& A, Eigen::VectorXd& f) = 0;
    virtual Eigen::VectorXd assembleRhs(const double& dt, Eigen::MatrixXd& M, Eigen::MatrixXd& A, Eigen::VectorXd& f, Eigen::VectorXd& uPrev) = 0;
    
    virtual Eigen::VectorXd solve(Eigen::MatrixXd& A, Eigen::VectorXd& rhs, const CirculationBCHandler& bcHandler)
    {
        bcHandler.addBC(A, rhs);
        return A.fullPivLu().solve(rhs);
    }
    
};


class ImplicitMidpointTimeIntegrator : public TimeIntegrator {
public:
    
    using TimeIntegrator::TimeIntegrator;
    
    virtual Eigen::MatrixXd assembleOperator(const double& dt, Eigen::MatrixXd& M, Eigen::MatrixXd& A, Eigen::VectorXd& f)
    {
        return M/dt + A/2;
    }
    
    virtual Eigen::VectorXd assembleRhs(const double& dt, Eigen::MatrixXd& M, Eigen::MatrixXd& A, Eigen::VectorXd& f, Eigen::VectorXd& uPrev)
    {
        return f + (M/dt - A/2) * uPrev;
    }

};


class ImplicitTimeIntegrator : public TimeIntegrator {
public:
    
    using TimeIntegrator::TimeIntegrator;
    
    virtual Eigen::MatrixXd assembleOperator(const double& dt, Eigen::MatrixXd& M, Eigen::MatrixXd& A, Eigen::VectorXd& f)
    {
        return M/dt + A;
    }
    
    virtual Eigen::VectorXd assembleRhs(const double& dt, Eigen::MatrixXd& M, Eigen::MatrixXd& A, Eigen::VectorXd& f, Eigen::VectorXd& uPrev)
    {
        return f + M/dt * uPrev;
    }
    
};


class ExplicitTimeIntegrator : public TimeIntegrator {
public:
    
    using TimeIntegrator::TimeIntegrator;
    
    virtual Eigen::MatrixXd assembleOperator(const double& dt, Eigen::MatrixXd& M, Eigen::MatrixXd& A, Eigen::VectorXd& f)
    {
        return M/dt;
    }
    
    virtual Eigen::VectorXd assembleRhs(const double& dt, Eigen::MatrixXd& M, Eigen::MatrixXd& A, Eigen::VectorXd& f, Eigen::VectorXd& uPrev)
    {
        return f + (M/dt - A) * uPrev;
    }
    
};
