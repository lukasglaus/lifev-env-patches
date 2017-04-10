//
//  Circulation.hpp
//  Circulation
//
//  Created by Thomas Kummer on 10.04.15.
//  Copyright (c) 2015 Thomas Kummer. All rights reserved.
//

#ifndef CIRCULATION_H_
#define CIRCULATION_H_

#include <iostream>
#include <iomanip>
#include <fstream>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <Eigen/Eigenvalues>
#include <vector>
#include <string>
#include <algorithm>

#include "CirculationIO.hpp"
#include "CirculationGridView.hpp"
#include "CirculationCoupling.hpp"
#include "CirculationDofHandler.hpp"
#include "CirculationAssembler.hpp"
#include "CirculationBCHandler.hpp"
#include "CirculationTimeIntegrator.hpp"
//#include "CirculationElementFactory.hpp"


class Circulation {
public:
    
    typedef std::vector<std::string> VectorStdString;
    typedef std::vector<VectorStdString> MatrixStdString;
    typedef std::vector<double> VectorStdDouble;
    typedef Eigen::MatrixXd MatrixEigen;
    typedef Eigen::VectorXd VectorEigen;
    
    Circulation(){}

    //Circulation(const Circulation&) = delete;

    virtual ~Circulation(){}
    
    Circulation(const std::string& filename) :
        M_time ( 0.0 )
    {
        // Read Grid and initial conditions
        readGrid( filename );
        
        // Initialize solution vector
        initialize();
    }
    
    void readGrid(const std::string& filename)
    {
        CirculationIO::readGrid(filename, M_gv);
    }

    void initialize()
    {
        DofHandler dofh(M_gv);

        // Resize solution vector and set it to zero
        M_u.resize( dofh.size() );
        M_u.setZero();
        
        // Set initial conditions in solution vector
        for ( auto& vertex : M_gv.vertices() ) M_u[ dofh(vertex) ] = vertex->init();
        for ( auto& element : M_gv.elements() ) M_u[ dofh(element) ] = element->init();

        // Initialize previous solution vectors
        M_uPrev0 = M_u;
        M_uPrev1 = M_u;
    }
    
    void solve(const double& dt, const MatrixStdString& bcNames = MatrixStdString(0), const VectorStdDouble& bcValues = VectorStdDouble(0), const unsigned int& iter = 0, const bool plotError = true, const bool plotSystem = false)
    {
        // Update time variable
        if ( iter == 0 ) updateTimeVar(dt);

        // Update solution vectors of previous timesteps
        if ( iter == 0 ) updatePrevSolVectors();

        // assemble mass- & stiffness matrix and source vector
        CirculationAssembler ca;
        auto A = ca.assembleStiffnessMatrix(M_gv, M_time, M_u);
        auto M = ca.assembleMassMatrix(M_gv, M_time, M_u);
        auto f = ca.assembleSourceVector(M_gv, M_time, M_u);
        
//        Eigen::GeneralizedEigenSolver<Eigen::MatrixXd> eigV;
//        eigV.compute(A, M);
//        std::cout << eigV.eigenvalues() << std::endl;

        // b.c. handler
        DofHandler dofh(M_gv);
        CirculationBCHandler bcHandler(dofh, bcNames, bcValues);

        // time integration
        ImplicitTimeIntegrator timeIntegrator;
        auto K = timeIntegrator.assembleOperator(dt, M, A, f);
        auto rhs = timeIntegrator.assembleRhs(dt, M, A, f, M_uPrev0);
        M_u = timeIntegrator.solve(K, rhs, bcHandler);

        // Plot linear sytem
        if ( plotSystem ) plotLinSys(K, M_u, rhs);

        // Plot relative error
        if ( plotError ) plotErrorNorm(K, M_u, rhs);
    }
    
    void iterate(const double& dt, const MatrixStdString& bcNames = MatrixStdString(0), const VectorStdDouble& bcValues = VectorStdDouble(0), const unsigned int& iter = 0, const bool plotError = false, const bool plotSystem = false, const double& error = 1e-6)
    {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        if ( rank == 0 )
        {
            unsigned int subiter (0);
            VectorEigen uPrevIter ( M_u );
            solve(dt, bcNames, bcValues, iter, plotError, plotSystem);
            VectorEigen residuum ( M_u - uPrevIter );

            std::cout << "\n\n=============================================================\n";
            std::cout << "Compute circulation\n";
            std::cout << "t = " << M_time << "\titer = " << subiter << "\tL2-Norm = " << residuum.norm() << std::endl;
            
            while ( residuum.norm() > error && subiter < 10 )
            {
                ++subiter;
                uPrevIter = M_u;
                solve(dt, bcNames, bcValues, 1, plotError, plotSystem);
                residuum = ( M_u - uPrevIter );
                
                std::cout << "\t\titer = " << subiter << "\tL2-Norm = " << residuum.norm() << std::endl;
            }
            
            std::cout << "=============================================================\n\n";
        }
	
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Bcast(M_u.data(), M_u.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
    
    CirculationCoupling coupling()
    {
        return CirculationCoupling();
    }
    
    void updateTimeVar(const double& dt)
    {
        setTime(M_time + dt);
    }
    
    void setTime(const double& time)
    {
        M_time = time;
    }
    
    void updatePrevSolVectors()
    {
        M_uPrev1 = M_uPrev0;
        M_uPrev0 = M_u;
    }
    
    VectorStdDouble eigenToStd(const VectorEigen& u) const
    {
        VectorStdDouble u_( u.size() );
        VectorEigen::Map( &u_[0] , u.size() ) = u;
        return u_;
    }

    Eigen::VectorXd stdToEigen(const VectorStdDouble& u) const
    {
        Eigen::VectorXd u_ = Eigen::VectorXd::Map(u.data(), u.size());
        return u_;
    }
    
    VectorStdDouble solution() const
    {
        return eigenToStd(M_u);
    }
    
    template<class type>
    const double solution(const type& bc)
    {
        DofHandler dofh(M_gv);
        const unsigned int idx = dofh( bc );
        return M_u[ idx ];
    }
    
    void plotLinSys(const MatrixEigen& A, const VectorEigen& u, const VectorEigen& rhs)
    {
        DofHandler dofh(M_gv);
        MatrixEigen ls( dofh.size() , dofh.size() + 2 );
        ls << A, u, rhs;
        std::cout << "t = " << M_time << ":" << std::endl << ls << std::endl;
    }
    
    void plotErrorNorm(const MatrixEigen& A, const VectorEigen& u, const VectorEigen& rhs) const
    {
        double relative_error = (A * u - rhs).norm() / rhs.norm();
        std::cout << "t = " << M_time << ":\t\t" << "rel. error = " << relative_error << std::endl;
    }
    
    void exportSolution(const std::string& filename, const bool& restart = false) const
    {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        
        if ( rank == 0 )
        {
            VectorStdDouble exportRow ({M_time});
            VectorStdDouble u ( eigenToStd(M_u) );
            exportRow.insert(exportRow.end(), u.begin(), u.end());
            
            CirculationIO exporter;
            const bool append = (M_time == 0 && restart == false ? false : true);
            exporter.exportVector(filename, exportRow, append);
        }
    }
    
    void restartFromFile(const std::string& filename, const unsigned int& timestep)
    {
        CirculationIO importer;

        std::vector<double> u ( importer.importVector(filename, timestep) );
        std::vector<double> uPrev0 ( importer.importVector(filename, (timestep < 1 ? 0 : (timestep - 1))) );
        std::vector<double> uPrev1 ( importer.importVector(filename, (timestep < 2 ? 0 : (timestep - 2))) );
        
        M_u = stdToEigen(u).tail( M_u.size() );
        M_uPrev0 = ( stdToEigen(uPrev0) ).tail( M_u.size() );
        M_uPrev1 = ( stdToEigen(uPrev1) ).tail( M_u.size() );
                                    
        M_time = u[0];

        DofHandler dofh (M_gv);
        for ( auto& element : M_gv.elements() )
        {
            const unsigned int elementIdx = dofh( element );
            const unsigned int vertex1Idx = dofh( element->node(0) );
            const unsigned int vertex2Idx = dofh( element->node(1) );
            
            std::vector<double> u {M_u(elementIdx) , 0.0, 0.0 };
            if ( vertex1Idx < dofh.size() ) u[1] = M_u(vertex1Idx);
            if ( vertex2Idx < dofh.size() ) u[2] = M_u(vertex2Idx);
            
            std::vector<double> uPrev0 {M_uPrev0(elementIdx) , 0.0, 0.0 };
            if ( vertex1Idx < dofh.size() ) uPrev0[1] = M_uPrev0(vertex1Idx);
            if ( vertex2Idx < dofh.size() ) uPrev0[2] = M_uPrev0(vertex2Idx);
            
            std::vector<double> uPrev1 {M_uPrev1(elementIdx) , 0.0, 0.0 };
            if ( vertex1Idx < dofh.size() ) uPrev1[1] = M_uPrev1(vertex1Idx);
            if ( vertex2Idx < dofh.size() ) uPrev1[2] = M_uPrev1(vertex2Idx);
            
            element -> initRestart ( u , uPrev0, uPrev1, M_time );
        }
    }
    

private:
    
    // Simulation time
    double M_time;
    
    // Grid
    GridView M_gv;

    // Solution vector
    Eigen::VectorXd M_u;
    
    // Previous solution vectors
    Eigen::VectorXd M_uPrev0;
    Eigen::VectorXd M_uPrev1;

};


#endif /* CIRCULATION_H_ */
