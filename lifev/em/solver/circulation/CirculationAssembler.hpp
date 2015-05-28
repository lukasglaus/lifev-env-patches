//
//  CirculationAssembler.cpp
//  Circulation
//
//  Created by Thomas Kummer on 16.05.15.
//  Copyright (c) 2015 Thomas Kummer. All rights reserved.
//

#ifndef CIRCULATIONASSEMBLER_HPP_
#define CIRCULATIONASSEMBLER_HPP_

#include <stdio.h>
#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <string>

#include "CirculationDofHandler.hpp"


class CirculationAssembler {
public:
    
    CirculationAssembler() {}
    
    virtual ~CirculationAssembler() {}
    
    enum M_var { Q , dQ , p1 , dp1 , p2 , dp2 };
    
    
    template<class T>
    Eigen::MatrixXd assembleMassMatrix(GridView& gv, const double& time, const T& U) const
    {
        DofHandler dofh(gv);
        
        // Initialize operator matrix
        Eigen::MatrixXd M ( dofh.size() , dofh.size() );
        M.setZero();
        
        for ( auto& element : gv.elements() )
        {
            // Determine global indices of element and vertices
            const unsigned int elementIdx = dofh( element );
            const unsigned int vertex1Idx = dofh( element->node(0) );
            const unsigned int vertex2Idx = dofh( element->node(1) );

            // Current solution of Q and p
            std::vector<double> u { U[elementIdx] , 0.0, 0.0 };
            if ( vertex1Idx < dofh.size() ) u[1] = U[vertex1Idx];
            if ( vertex2Idx < dofh.size() ) u[2] = U[vertex2Idx];
            
            // Add entries for element
            M ( elementIdx , elementIdx ) = element->lhs( dQ , u , time );
            if ( vertex1Idx < dofh.size() ) M ( elementIdx , vertex1Idx ) = element->lhs( dp1 , u , time );
            if ( vertex2Idx < dofh.size() ) M ( elementIdx , vertex2Idx ) = element->lhs( dp2 , u , time );
        }
        
        return M;
    }
    
    
    template<class T>
    Eigen::MatrixXd assembleStiffnessMatrix(GridView& gv, const double& time, const T& U) const
    {
        DofHandler dofh(gv);
        
        // Initialize operator matrix
        Eigen::MatrixXd A ( dofh.size() , dofh.size() );
        A.setZero();
        
        for ( auto& element : gv.elements() )
        {
            // Determine global indices of element and vertices
            const unsigned int elementIdx = dofh( element );
            const unsigned int vertex1Idx = dofh( element->node(0) );
            const unsigned int vertex2Idx = dofh( element->node(1) );
            
            // Current solution of Q and p
            std::vector<double> u { U[elementIdx] , 0.0, 0.0 };
            if ( vertex1Idx < dofh.size() ) u[1] = U[vertex1Idx];
            if ( vertex2Idx < dofh.size() ) u[2] = U[vertex2Idx];
            
            // Add entries for continuity at node
            if ( vertex1Idx < dofh.size() ) A ( vertex1Idx , elementIdx ) = - 1.0;
            if ( vertex2Idx < dofh.size() ) A ( vertex2Idx , elementIdx ) =   1.0;
            
            // Add entries for element
            A ( elementIdx , elementIdx ) = element->lhs( Q , u , time );
            if ( vertex1Idx < dofh.size() ) A ( elementIdx , vertex1Idx ) = element->lhs( p1 , u , time );
            if ( vertex2Idx < dofh.size() ) A ( elementIdx , vertex2Idx ) = element->lhs( p2 , u , time );
        }
        
        return A;
    }

    
    template<class T>
    Eigen::VectorXd assembleSourceVector(GridView& gv, const double& time, T& U) const
    {
        DofHandler dofh(gv);
        
        // Initialize operator matrix
        Eigen::VectorXd f ( dofh.size() );
        f.setZero();
        
        for ( auto& element : gv.elements() )
        {
            // Determine global indices of element and vertices
            const unsigned int elementIdx = dofh( element );
            const unsigned int vertex1Idx = dofh( element->node(0) );
            const unsigned int vertex2Idx = dofh( element->node(1) );
            
            // Current solution of Q and p
            std::vector<double> u { U[elementIdx] , 0.0, 0.0 };
            if ( vertex1Idx < dofh.size() ) u[1] = U[vertex1Idx];
            if ( vertex2Idx < dofh.size() ) u[2] = U[vertex2Idx];
            
            // Add entries for element
            f ( elementIdx ) = element->rhs( u , time );
        }
        
        return f;
    }

    // class TimeIntegrator
        // class BackwardEuler
            // assembleOperator
            // assembleRhs
            // solve method from solver
    
    // class CirculationBCHandler
    
};

#endif