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
 *  @file
 *  @brief This file contains solvers for different active materials.
 *  @warning: This is the most important issue related with this class.
 *  At the moment, BC they are applied on the residual directly. This does
 *  not work for nonhomogeneus Dirichlet conditions!!
 *  Check the structural operator class in the structure module
 *
 *  @version 1.0
 *  @date 04-2014
 *  @author Paolo Tricerri
 *
 *  @maintainer  Simone Rossi <simone.rossi@epfl.ch>
*/

#ifndef _EMSTRUCTURALOPERATOR_H_
#define _EMSTRUCTURALOPERATOR_H_ 1

#include <lifev/em/solver/mechanics/EMStructuralConstitutiveLaw.hpp>
#include <lifev/core/fem/BCVector.hpp>
#include <lifev/structure/solver/StructuralOperator.hpp>
#include <lifev/em/util/EMUtility.hpp>
#include <memory>


namespace LifeV
{

using namespace ExpressionAssembly;

/*!
  \class EMStructuralSolver
  \brief
  This class solves the elastodynamics equations for different active materials


*/
template <typename Mesh>
class EMStructuralOperator : public StructuralOperator<Mesh>
{
public:

    //!@name Type definitions
    //@{
    typedef StructuralOperator<Mesh>                    structuralOperator_Type;

    typedef boost::shared_ptr<structuralOperator_Type>  structuralOperatorPtr_Type;

    typedef EMStructuralConstitutiveLaw<Mesh>     material_Type;

    typedef boost::shared_ptr<material_Type>            materialPtr_Type;

    //Old typedefs
    typedef StructuralOperator<Mesh>                    super;

    typedef typename super::data_Type                   data_Type;

    typedef typename super::FESpacePtr_Type             FESpacePtr_Type;

    typedef typename super::ETFESpacePtr_Type           ETFESpacePtr_Type;

    typedef typename super::bcHandler_Type              bcHandler_Type;
    //@}

    //! @name Constructor & Destructor
    //@{

    EMStructuralOperator();

    virtual ~EMStructuralOperator()
    {
        M_EMMaterial.reset();
    };


    //! Setup the created object of the class Venantkirchhof
    /*!
      \param data_file GetPot data file
      \param refFE reference FE for the displacement
      \param BCh boundary conditions for the displacement
      \param comm the Epetra Comunicator
    */
    void setup ( boost::shared_ptr<data_Type>  data,
                 const FESpacePtr_Type&        dFESpace,
                 const ETFESpacePtr_Type&      dETFESpace,
                 bcHandler_Type&       BCh,
                 boost::shared_ptr<Epetra_Comm>&     comm
               );

    /*! Get the offset parameter. It is taken into account when the boundary conditions
      are applied and the matrices are assembled.
    */
    const materialPtr_Type& EMMaterial() const
    {
        return M_EMMaterial;
    }


    //    const materialPtr_Type& EMMaterial() const
    //    {
    //        return M_EMMaterial;
    //    }

    //@}
    FESpacePtr_Type dispFESpacePtr()
    {
        return super::M_dispFESpace;
    }

    void setPressureBC( Real p )
    {
    	M_LVPressure = p;
    }

    Real pressureBC() const
    {
    	return M_LVPressure;
    }

    void setBCFlag(ID flag)
    {
    	M_LVPressureFlag = flag;
    }
    ID BCFlag()
    {
    	return M_LVPressureFlag;
    }

    void solveJac ( vector_Type& step, const vector_Type& res, Real& linear_rel_tol);

    void computePressureBC(const VectorEpetra& disp,
			boost::shared_ptr<VectorEpetra> bcVectorPtr,
			const ETFESpacePtr_Type dETFESpace,
			Real pressure, int bdFlag);

    void computePressureBCJacobian(const VectorEpetra& disp,
    		matrixPtr_Type& jacobian,
			const ETFESpacePtr_Type dETFESpace,
			Real pressure, int bdFlag);


    void evalResidual ( vector_Type& residual, const vector_Type& solution, Int iter);

    void iterate ( const bcHandler_Type& bch, bool pressureBC = false );


    void updateJacobian ( const vector_Type& sol, matrixPtr_Type& jacobian  );

    void setI4fPtr(vectorPtr_Type ptr)
    {
    	M_I4fPtr = ptr;
    }


    void setI4f(vector_Type& i4f)
    {
    	*M_I4fPtr = i4f;
    }

    vectorPtr_Type I4fPtr()
    {
    	return M_I4fPtr;
    }

protected:

    //! Material class
    //    structuralOperatorPtr_Type             M_structuralOperator;
    materialPtr_Type                     M_EMMaterial;

    //left ventricular pressure
    Real                                 M_LVPressure;
    ID								M_LVPressureFlag;
    vectorPtr_Type M_boundaryVectorPtr;
    boost::shared_ptr<BCVector>  M_bcVectorPtr;
    bool M_LVpressureBC;


    vectorPtr_Type  M_I4fPtr;
};

//====================================
// Constructor
//=====================================

template <typename Mesh>
EMStructuralOperator<Mesh>::EMStructuralOperator( ) :
    //    M_structuralOperator(),
	M_LVPressure(),
	M_LVPressureFlag(),
    super(),
    M_EMMaterial(),
    M_boundaryVectorPtr(),
    M_bcVectorPtr()
{
}

template <typename Mesh>
void
EMStructuralOperator<Mesh>::setup (boost::shared_ptr<data_Type>          data,
                                   const FESpacePtr_Type& dFESpace,
                                   const ETFESpacePtr_Type& dETFESpace,
                                   bcHandler_Type&                    BCh,
                                   boost::shared_ptr<Epetra_Comm>&   comm)
{

    this->super::setup (data, dFESpace, dETFESpace, BCh, comm);
    M_EMMaterial = boost::dynamic_pointer_cast<material_Type> (this -> material() );
    M_boundaryVectorPtr.reset(new vector_Type ( this->M_disp->map(), Repeated ) );
    M_I4fPtr.reset(new vector_Type ( M_EMMaterial->scalarETFESpacePtr()->map() ) );
    *M_I4fPtr += 1.0;

}


template <typename Mesh>
void
EMStructuralOperator<Mesh>::computePressureBC(  const VectorEpetra& disp,
										boost::shared_ptr<VectorEpetra> bcVectorPtr,
										const ETFESpacePtr_Type dETFESpace,
										Real pressure, int bdFlag)
{

	*bcVectorPtr *= 0.0;

    MatrixSmall<3,3> Id;
    Id(0,0) = 1.; Id(0,1) = 0., Id(0,2) = 0.;
    Id(1,0) = 0.; Id(1,1) = 1., Id(1,2) = 0.;
    Id(2,0) = 0.; Id(2,1) = 0., Id(2,2) = 1.;

	{
	  	using namespace ExpressionAssembly;

	  	auto I = value(Id);
	  	auto Grad_u = grad( dETFESpace, disp, 0);
	  	auto F =  Grad_u + I;
	  	auto FmT = minusT(F);
	  	auto J = det(F);
	  	auto p = value(pressure);

    	QuadratureBoundary myBDQR (buildTetraBDQR (quadRuleTria7pt) );

        integrate ( boundary ( dETFESpace->mesh(), bdFlag),
        		    myBDQR,
        		    dETFESpace,
        		    p * J * dot( FmT * Nface,  phi_i)
                  ) >> bcVectorPtr;

        bcVectorPtr -> globalAssemble();


	}
}


template <typename Mesh>
void
EMStructuralOperator<Mesh>::evalResidual ( vector_Type& residual, const vector_Type& solution, Int iter)
{
    if(M_LVpressureBC)
    {
		if(M_LVPressure != 0 && M_LVPressureFlag != 0)
		{
			this->M_Displayer->leaderPrint ("\n    S- Updating the pressure boundary conditions: pressure = ", M_LVPressure, "\n");
			computePressureBC(  solution, M_boundaryVectorPtr, this->M_dispETFESpace, M_LVPressure, M_LVPressureFlag);
			M_bcVectorPtr.reset( new BCVector (*M_boundaryVectorPtr, this->M_dispFESpace -> dof().numTotalDof(), 0 ) );
			this-> M_BCh -> modifyBC(M_LVPressureFlag, *M_bcVectorPtr);

		}
		else
		{
			this->M_Displayer->leaderPrint ("\n    S- Not Updating the pressure boundary conditions\n");
		}
    }
    super::evalResidual(residual, solution, iter);
}


//solveJac( const Vector& res, Real& linear_rel_tol, Vector &step)
template <typename Mesh>
void EMStructuralOperator<Mesh>::
solveJac ( vector_Type& step, const vector_Type& res, Real& linear_rel_tol)
{
	*this->M_jacobian *= 0.0;
    if(M_LVpressureBC)
    computePressureBCJacobian(*this->M_disp, this->M_jacobian, this->M_dispETFESpace, M_LVPressure, M_LVPressureFlag);
    updateJacobian ( *this->M_disp, this->M_jacobian );
    this->M_jacobian -> globalAssemble();
    this->solveJacobian (step,  res, linear_rel_tol, this->M_BCh);
}

template <typename Mesh>
void EMStructuralOperator<Mesh>::computePressureBCJacobian(const VectorEpetra& disp,
		matrixPtr_Type& jacobian,
		const ETFESpacePtr_Type dETFESpace,
		Real pressure, int bdFlag)
{

	if(M_LVPressure != 0 && M_LVPressureFlag != 0)
	{
		if(disp.comm().MyPID() == 0)
		{
			std::cout << "\nUpdating the pressure BC Jacobian: pressure = " << pressure << "\n";
		}

		MatrixSmall<3,3> Id;
		Id(0,0) = 1.; Id(0,1) = 0., Id(0,2) = 0.;
		Id(1,0) = 0.; Id(1,1) = 1., Id(1,2) = 0.;
		Id(2,0) = 0.; Id(2,1) = 0., Id(2,2) = 1.;


		{
			using namespace ExpressionAssembly;

			auto dF = grad (phi_j);
			auto I = value(Id);
			auto Grad_u = grad( dETFESpace, disp, 0);
			auto F =  Grad_u + I;
			auto FmT = minusT(F);
			auto J = det(F);
			auto p = value(pressure);
			auto dFmT = value(-1.0) * FmT * transpose( dF ) * FmT;
			auto dJ   = J * dot( FmT, dF ) * I;
			//since it's on the left side we put a minus
			auto dP = value(-1.0) * p * (dJ * FmT + J * dFmT );

			QuadratureBoundary myBDQR (buildTetraBDQR (quadRuleTria7pt) );

			Real normJac1 = jacobian -> norm1();
			jacobian->spy("Jac1");
			jacobian -> openCrsMatrix();
			integrate ( boundary ( dETFESpace->mesh(), bdFlag),
						myBDQR,
						dETFESpace,
						dETFESpace,
						dot( dP * Nface,  phi_i )
					  ) >> jacobian;

			Real normJac = jacobian -> norm1();

			if(disp.comm().MyPID() == 0)
			{
				std::cout << "\nNormJac 1: " << normJac1 << ",\tNormJac 1: " << normJac << "\n";
			}
		}




	}
	else
	{
		if(disp.comm().MyPID() == 0)
		{
			std::cout << "\nI'm not pdating the pressure BC: pressure = " << M_LVPressure << "\n";
		}
	}

}


template <typename Mesh>
void
EMStructuralOperator<Mesh>::iterate ( const bcHandler_Type& bch, bool pressureBC )
{
    LifeChrono chrono;

    M_LVpressureBC = pressureBC;
    // matrix and vector assembling communication
    this->M_Displayer->leaderPrint ("  EMSolver -  Solving the system ... \n");

    this->M_BCh = bch;

    Real abstol  = this->M_nonlinearParameters.M_abstol;
    Real reltol  = this->M_nonlinearParameters.M_reltol;
    UInt maxiter = this->M_nonlinearParameters.M_maxiter;
    Real etamax  = this->M_nonlinearParameters.M_etamax;
    Int NonLinearLineSearch = this->M_nonlinearParameters.M_NonLinearLineSearch;

    Real time = this->M_data->dataTime()->time();

    Int status = 0;

    if ( this->M_data->verbose() )
    {
        status = NonLinearRichardson ( *this->M_disp,
        		                       *this,
        		                       abstol,
        		                       reltol,
        		                       maxiter,
        		                       etamax,
        		                       NonLinearLineSearch,
        		                       0,
        		                       2,
        		                       this->M_out_res,
        		                       this->M_data->dataTime()->time() );
    }
    else
    {
        status = NonLinearRichardson ( *this->M_disp,
        							   *this,
        							   abstol,
									   reltol,
									   maxiter,
									   etamax,
									   NonLinearLineSearch );
    }


    if ( status == 1 )
    {
        std::ostringstream ex;
        ex << "EMStructuralOperator::iterate() Inners nonLinearRichardson iterations failed to converge\n";
        throw std::logic_error ( ex.str() );
    }
    else // if status == 0 NonLinearrRichardson converges
    {
        // std::cout << std::endl;

        // std::cout <<" Number of inner iterations       : " << maxiter <<  std::endl;

        // std::cout <<" We are at the time step          : "  << M_data->dataTime()->time() << std::endl;
        if ( this->M_data->verbose() )
        {
        	this->M_out_iter << time << " " << maxiter << std::endl;
        }
    }

}



//Method UpdateJacobian
template <typename Mesh>
void
EMStructuralOperator<Mesh>::updateJacobian ( const vector_Type& sol, matrixPtr_Type& jacobian  )
{
    this->M_Displayer->leaderPrint ("  S-  Solid: Updating JACOBIAN... ");

    LifeChrono chrono;
    chrono.start();

    this->M_material->updateJacobianMatrix (sol, this->M_data, this->M_mapMarkersVolumes, this->M_mapMarkersIndexes,  this->M_Displayer);

    *jacobian += * (this->M_material->jacobian() );

    //Spying the static part of the Jacobian to check if it is symmetric
    // M_material->jacobian()->spy("staticJacobianMatrix");

    // std::cout << "spyed" << std::endl;
    // int n;
    // std::cin >> n ;


    chrono.stop();
    this->M_Displayer->leaderPrintMax ("   ... done in ", chrono.diff() );

}


}
#endif
