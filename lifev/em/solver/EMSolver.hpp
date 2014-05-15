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
 @brief Class for solving the Monodomain equations in electrophysiology.

 @date 02-2013
 @author Simone Rossi <simone.rossi@epfl.ch>

 @last update 02-2013

 This class provides interfaces to solve the monodomain equation
 ( reaction diffusion equation ) using the ETA framework.
 The solution can be performed using three different methods:
 -operator splitting method (at this point available only with forward Euler
 for the reaction step and backward Euler for the diffusion step. );
 -Ionic Currents Interpolation (at this point only forward Euler);
 -State Variable interpolation (at this point only forward Euler).
 */

#ifndef _EMSOLVER_H_
#define _EMSOLVER_H_


#include <lifev/em/solver/electrophysiology/EMMonodomainSolver.hpp>
//
//#include <lifev/em/solver/EMStructuralOperator.hpp>
//#include <lifev/em/solver/EMGeneralizedActiveHolzapfelOgdenMaterial.hpp>
//#include <lifev/em/solver/EMActiveStrainSolver.hpp>
//#include <lifev/core/interpolation/RBFlocallyRescaledVectorial.hpp>
//#include <lifev/core/interpolation/RBFlocallyRescaledScalar.hpp>
//#include <lifev/core/interpolation/RBFrescaledVectorial.hpp>
//#include <lifev/core/interpolation/RBFrescaledScalar.hpp>
////#include <lifev/core/interpolation/RBFscalar.hpp>
//#include <lifev/core/interpolation/RBFvectorial.hpp>
//
//#include <lifev/bc_interface/3D/bc/BCInterface3D.hpp>
//
//
//#include <lifev/em/solver/EMEvaluate.hpp>

namespace LifeV
{

//! EMSolver - Class featuring the solution of the electromechanical problem with monodomain equation

template<typename Mesh , typename IonicModel, typename ActivationModel>
class EMSolver
{

	EMSolver();

	EMSolver(const EMSolver& solver);

	EMSolver(EMSolver&& solver);

public:
	EMMonodomainSolver<Mesh, IonicModel> M_monodomain;
	ActivationModel M_activationModel;
};

/////////////////////
// CONSTRUCTORS
template<typename Mesh , typename IonicModel, typename ActivationModel>
EMSolver<Mesh, IonicModel, ActivationModel>::EMSolver()
{
}


/////////////////////
// COPY CONSTRUCTORS
template<typename Mesh , typename IonicModel, typename ActivationModel>
EMSolver<Mesh, IonicModel, ActivationModel>::EMSolver(const EMSolver& solver) :
    M_monodomain(solver.M_monodomain),
    M_activationModel(solver.M_activationModel)
{
}

/////////////////////
// MOVE CONSTRUCTORS
template<typename Mesh , typename IonicModel, typename ActivationModel>
EMSolver<Mesh, IonicModel, ActivationModel>::EMSolver(EMSolver&& solver) :
    M_monodomain(solver.M_monodomain)
{
}



} // namespace LifeV




#endif //_MONODOMAINSOLVER_H_
