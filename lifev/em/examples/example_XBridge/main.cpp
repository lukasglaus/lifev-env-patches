//@HEADER
/*
*******************************************************************************

    Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
    Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University

    This file is part of LifeV.

    LifeV is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    LifeV is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with LifeV.  If not, see <http://www.gnu.org/licenses/>.

*******************************************************************************
*/
//@HEADER

/*!
    @file
    @brief 0D test with Ten Tusscher et al. 2006 model

    Note that the model is not solved correctly using Rush-Larsen
    with forward Euler. In the original code of Ten Tusscher
    the solution of CaSS and CaSR is computed using a specific
    algorithm (which I do not know).
    Therefore I consider all the variables except for the
    potential as gating variable and use the specific
    method in their original code to make it work.

    @date 08 - 2013
    @author Simone Rossi <simone.rossi@epfl.ch>

    @contributor
    @mantainer Simone Rossi <simone.rossi@epfl.ch>
 */

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

//Tell the compiler to restore the warning previously silented
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"


#include <fstream>


#include <lifev/electrophysiology/solver/IonicModels/IonicTenTusscher06.hpp>
#include <lifev/core/LifeV.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include "Teuchos_XMLParameterListHelpers.hpp"

#include <lifev/em/solver/activation/xBridgeModels/XBridge4SM.hpp>


using namespace LifeV;

#define SolutionTestNorm  -2476.4745158656560307

Int main ( Int argc, char** argv )
{
    //! Initializing Epetra communicator
    MPI_Init (&argc, &argv);
    Epetra_MpiComm Comm (MPI_COMM_WORLD);
    if ( Comm.MyPID() == 0 )
    {
        std::cout << "% using MPI" << std::endl;
    }

    //********************************************//
    // Creates a new model object representing the//
    // model from TenTusscher et al. 2006. The    //
    // model input are the parameters. Pass  the  //
    // parameter list in the constructor          //
    //********************************************//
    std::cout << "Building Constructor for TenTusscher 2006 Model with default parameters ... ";
    IonicTenTusscher06  ionicModel;
    std::cout << " Done!" << std::endl;


    //********************************************//
    // Show the parameters of the model as well as//
    // other informations  about the object.      //
    //********************************************//
    ionicModel.showMe();


    //********************************************//
    // Initialize the solution with the default   //
    // values                                     //
    //********************************************//
    std::cout << "Initializing solution vector...";
    std::vector<Real> states (ionicModel.restingConditions() );
    std::cout << " Done!" << std::endl;


    //********************************************//
    // Initialize the rhs to 0. The rhs is the    //
    // vector containing the numerical values of  //
    // the time derivatives of the state          //
    // variables, that is, the right hand side of //
    // the differential equation.                 //
    //********************************************//
    std::cout << "Initializing rhs..." ;
    std::vector<Real> rhs (ionicModel.Size(), 0);
    std::cout << " Done! "  << std::endl;




    //********************************************//
    // The model needs as external informations   //
    // the contraction velocity and the Calcium   //
    // concentration.                             //
    //********************************************//
    Real Iapp (0.0);

    
    //********************************************//
    // Setup X-Bridge Model.                      //
    //********************************************//
    XBridge4SM xb;
    std::vector<Real> statesXb ( xb.xbVarVec() );


    //********************************************//
    // Simulation starts on t=0 and ends on t=TF. //
    // The timestep is given by dt                //
    //********************************************//
    Real TF (15500.0);
    Real dt (0.001);


    //********************************************//
    // We record the norm of the solution to      //
    // check the failure of the test              //
    //********************************************//
    Real SolutionNorm = states[0];


    //********************************************//
    // Open the file "output.txt" to save the     //
    // solution.                                  //
    //********************************************//
    std::string filename = "output.txt";
    std::ofstream output ("output.txt");
    output << 0 << "\t";

    for ( int j (0); j < ionicModel.Size(); j++)
    {
        output << states[j] << "\t";
    }
    
    for ( int jXb (0); jXb < xb.size(); ++jXb )
    {
        output << statesXb[jXb] << "\t";
    }
    
    output << "\n";
    
    
    //********************************************//
    // Time loop starts.                          //
    //********************************************//
    std::cout << "Time loop starts...\n";


    int savestep ( 10.0 / dt );
    int iter (0);
    int stimulusStep (1000 / dt);
    int stimulusStepBreak (1 / dt);


    std::vector<Real> v (states);
    for ( Real t = 0; t < TF; )
    {

        //********************************************//
        // Compute the applied current. This is a     //
        // simple switch.                             //
        //********************************************//
//        if ( t > 0 && t < 1 && iter % stimulusStep == 0)
        if ( iter % stimulusStep == 0)
        {
            Iapp = 52.;
        }
        if( ( iter - stimulusStepBreak ) % stimulusStep == 0)
        {
            Iapp = 0;
        }
        std::cout << "\r " << t << " ms.       " <<  std::flush;


        iter++;
        //********************************************//
        // Compute the rhs using the model equations  //
        //********************************************//
        ionicModel.setAppliedCurrent (Iapp);
        rhs[0] = ionicModel.computeLocalPotentialRhs (states);
        ionicModel.addAppliedCurrent (rhs);

        //ionicModel.computeGatingVariablesWithRushLarsen ( states, dt);
            Real V = states[0];
            Real m = states[1];
            Real h = states[2];
            Real j = states[3];
            Real d = states[4];
            Real f = states[5];
            Real f2 = states[6];
            Real fcass = states[7];
            Real r = states[8];
            Real s = states[9];
            Real xr1 = states[10];
            Real xr2 = states[11];
            Real xs = states[12];
            Real Nai = states[13];
            Real Ki = states[14];
            Real Cai = states[15];
            Real CaSS = states[16];
            Real CaSR = states[17];
            Real RR = states[18];
        
            states[18] = ionicModel.solveRR(CaSR, CaSS, RR, dt);
            states[17] = ionicModel.solveCaSR(Cai, CaSR, CaSS, RR, dt);
            states[16] = ionicModel.solveCaSS(Cai, CaSR, CaSS, RR, V, d, f, f2, fcass, dt);
            states[15] = ionicModel.solveCai(V, Nai, Cai, CaSR, CaSS, dt);
            states[14] = ionicModel.solveKi(V, r, s, xr1, xr2, xs, Nai, Ki, dt);
            states[13] = ionicModel.solveNai(V, m, h, j, Nai, Cai, dt);
            states[12] = ionicModel.solveXs(V, xs, dt);
            states[11] = ionicModel.solveXr2(V, xr2, dt);
            states[10] = ionicModel.solveXr1(V, xr1, dt);
            states[9] = ionicModel.solveS(V, s, dt);
            states[8] = ionicModel.solveR(V, r, dt);
            states[7] = ionicModel.solveFCaSS(V, fcass, dt);
            states[6] = ionicModel.solveF2(V, f2, dt);
            states[5] = ionicModel.solveF(V, f, dt);
            states[4] = ionicModel.solveD(V, d, dt);
            states[3] = ionicModel.solveJ(V, j, dt);
            states[2] = ionicModel.solveH(V, h, dt);
            states[1] = ionicModel.solveM(V, m, dt);
            states[0] = ionicModel.solveV(V, m, h, j, d, f, f2, fcass, r, s, xr1, xr2, xs, Nai, Ki, Cai, CaSS, dt);
        
        
        //states[0] = states[0]  + dt * rhs[0]; //This should be divided by the membrane capacitance  ionicModel.membraneCapacitance();
        //In the code of Ten Tusscher there is no such a division
        // on the other hand if we don't divide by the membrane capacitance
        // in the three dimensional simulations
        // we get a faster wave (as fast as the LuoRudy ~2 the one in the benchmark!)

        
        //********************************************//
        // Compute X-Bridge Model                     //
        //********************************************//
        xb.solveFE ( statesXb, states[15], dt );
        

        //********************************************//
        // Update the time.                           //
        //********************************************//
        t = t + dt;


        //********************************************//
        // Write solution on file.                    //
        //********************************************//
        if ( iter % savestep == 0  )
        {
            // Write time data
            output << t << "\t";
            
            // Write ionic model data
            for ( int index (0); index < ionicModel.Size(); index++)
            {
                output << states[index] << "\t";
            }
            
            // Write x-bridge model data
            for ( int indexXb (0); indexXb < xb.size(); ++indexXb )
            {
                output << statesXb[indexXb] << "\t";
            }
            
            output << "\n";
            
            
            //********************************************//
            // Update the norm of the solution to check   //
            // test failure                               //
            //********************************************//
            SolutionNorm += states[0];
        }

    }
    std::cout << "\n...Time loop ends.\n";
    std::cout << "Solution written on file: " << filename << "\n";
    //********************************************//
    // Close exported file.                       //
    //********************************************//
    output.close();

    //! Finalizing Epetra communicator
    MPI_Finalize();

    Real returnValue;

    Real err = std::abs (SolutionTestNorm - SolutionNorm) / std::abs (SolutionTestNorm);
    std::cout << std::setprecision (20) << "\nError: " << err << "\nSolution norm: " << SolutionNorm << "\n";
    if ( err > 1e-12 )
    {
        std::cout << "\nTest Failed!\n";
        returnValue = EXIT_FAILURE; // Norm of solution did not match
    }
    else
    {
        returnValue = EXIT_SUCCESS;
    }
    return ( returnValue );
}

#undef SolutionTestNorm
