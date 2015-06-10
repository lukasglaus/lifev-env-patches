//
//  CirculationBCHandler.cpp
//  Circulation
//
//  Created by Thomas Kummer on 17.05.15.
//  Copyright (c) 2015 Thomas Kummer. All rights reserved.
//

#include <stdio.h>
#include <Eigen/Dense>

#include "CirculationDofHandler.hpp"


class CirculationBCHandler {
public:
    
    CirculationBCHandler(DofHandler dofh, const std::vector<std::vector<std::string> >& bcNames, const std::vector<double>& bcValues) :
        M_dofh(dofh),
        M_bcNames(bcNames),
        M_bcValues(bcValues)
    {}
    
    virtual ~CirculationBCHandler() {}
    
    void addBC(Eigen::MatrixXd& A, Eigen::VectorXd& rhs) const
    {
        for ( unsigned int bcIdx (0) ; bcIdx < M_bcNames.size() ; ++bcIdx )
        {
            const unsigned int vertexidx = M_dofh( M_bcNames[ bcIdx ][0] );
            if ( M_bcNames[ bcIdx ][1] == "p" )
            {
                setPressureBC(A, rhs, vertexidx, M_bcValues[ bcIdx ]);
            }
            else if ( M_bcNames[ bcIdx ][1] == "Q" )
            {
                setFlowRateBC(A, rhs, vertexidx, M_bcValues[ bcIdx ]);
            }
        }
    }
    
    
private:
    
    void setPressureBC(Eigen::MatrixXd& A, Eigen::VectorXd& rhs, const unsigned int& vertexIdx, const double& bcValue) const
    {
        A.row(vertexIdx).setZero();
        A ( vertexIdx , vertexIdx ) = 1.0;
        rhs ( vertexIdx ) = bcValue;
    }

    
    void setFlowRateBC(Eigen::MatrixXd& A, Eigen::VectorXd& rhs, const unsigned int& vertexIdx, const double& bcValue) const
    {
        rhs ( vertexIdx ) = - bcValue;
    }
    
    
    DofHandler M_dofh;
    
    const std::vector<std::vector<std::string> > M_bcNames;
    const std::vector<double> M_bcValues;
    
};