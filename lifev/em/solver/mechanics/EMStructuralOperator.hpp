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
#include <lifev/structure/solver/StructuralOperator.hpp>


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
	typedef StructuralOperator<Mesh>					structuralOperator_Type;

	typedef boost::shared_ptr<structuralOperator_Type>  structuralOperatorPtr_Type;

    typedef EMStructuralConstitutiveLaw<Mesh>     material_Type;

    typedef boost::shared_ptr<material_Type>			materialPtr_Type;

 //Old typedefs
    typedef StructuralOperator<Mesh>					super;

    typedef typename super::data_Type					data_Type;

    typedef typename super::FESpacePtr_Type				FESpacePtr_Type;

    typedef typename super::ETFESpacePtr_Type			ETFESpacePtr_Type;

    typedef typename super::bcHandler_Type				bcHandler_Type;
    //@}

    //! @name Constructor & Destructor
    //@{

    EMStructuralOperator();

    virtual ~EMStructuralOperator()
    {
    	M_EMMaterial.reset();
    };


    //!@name Type definitions
    //@{
//    typedef Real ( *function ) ( const Real&, const Real&, const Real&, const Real&, const ID& );
//
//    typedef boost::function<Real ( Real const&, Real const&, Real const&, Real const&, ID const& ) > source_Type;
//
//    typedef StructuralConstitutiveLaw<Mesh>               material_Type;
//    typedef boost::shared_ptr<material_Type>              materialPtr_Type;
//
//    typedef BCHandler                                     bcHandlerRaw_Type;
//    typedef boost::shared_ptr<bcHandlerRaw_Type>          bcHandler_Type;
//
//    typedef LinearSolver                                  solver_Type;
//
//    typedef typename solver_Type::matrix_Type             matrix_Type;
//    typedef boost::shared_ptr<matrix_Type>                matrixPtr_Type;
//    typedef typename solver_Type::vector_Type             vector_Type;
//    typedef boost::shared_ptr<vector_Type>                vectorPtr_Type;
//    typedef vector_Type                                   solution_Type;
//    typedef boost::shared_ptr<solution_Type>              solutionPtr_Type;
//
//    typedef StructuralConstitutiveLawData                 data_Type;
//
//    typedef RegionMesh<LinearTetra >                      mesh_Type;
//    typedef  boost::shared_ptr<mesh_Type>                 meshPtr_Type;
//    typedef std::vector< mesh_Type::element_Type* >       vectorVolumes_Type;
//    typedef std::vector< UInt >                           vectorIndexes_Type;
//
//    typedef std::map< UInt, vectorVolumes_Type>           mapMarkerVolumes_Type;
//    typedef std::map< UInt, vectorIndexes_Type>           mapMarkerIndexes_Type;
//    typedef boost::shared_ptr<mapMarkerVolumes_Type>      mapMarkerVolumesPtr_Type;
//    typedef boost::shared_ptr<mapMarkerIndexes_Type>      mapMarkerIndexesPtr_Type;
//    typedef mapMarkerVolumes_Type::const_iterator         mapIterator_Type;
//
//    typedef typename mesh_Type::element_Type              meshEntity_Type;
//
//    typedef typename boost::function2<bool, const UInt, const UInt> comparisonPolicy_Type;
//
//    typedef MarkerSelector<meshEntity_Type, comparisonPolicy_Type> markerSelector_Type;
//    typedef boost::scoped_ptr<markerSelector_Type>          markerSelectorPtr_Type;
//
//    typedef ETFESpace< RegionMesh<LinearTetra>, MapEpetra, 3, 3 >  ETFESpace_Type;
//    typedef boost::shared_ptr<ETFESpace_Type>                      ETFESpacePtr_Type;
//
//    typedef FESpace< RegionMesh<LinearTetra>, MapEpetra >          FESpace_Type;
//    typedef boost::shared_ptr<FESpace_Type>                        FESpacePtr_Type;
//
//    //Preconditioners typedef
//    typedef LifeV::Preconditioner                   basePrec_Type;
//    typedef boost::shared_ptr<basePrec_Type>        basePrecPtr_Type;
//    typedef LifeV::PreconditionerIfpack             precIfpack_Type;
//    typedef boost::shared_ptr<precIfpack_Type>      precIfpackPtr_Type;
//    typedef LifeV::PreconditionerML                 precML_Type;
//    typedef boost::shared_ptr<precML_Type>          precMLPtr_Type;
//
//    // Time advance
//    typedef TimeAdvance< vector_Type >                                  timeAdvance_Type;
//    typedef boost::shared_ptr< timeAdvance_Type >                       timeAdvancePtr_Type;

    //@}

    //@}

    //!@name Methods
    //@{

//    //! Updates the system at the end of each time step
//    inline void updateSystem ( void )
//    {
//    	M_structuralOperator -> updateSystem();
//    }
//
//    //! Updates the system at the end of each time step when the matrix is passed from outside
//    /*!
//      \param stiff stiffness matrix provided from outside
//    */
//    inline void updateSystem ( matrixPtr_Type& mat_stiff)
//    {
//    	M_structuralOperator -> updateSystem(mat_stiff);
//    }
//
//    //! Updates the system at the end of each time step given a source term
//    /*!
//      \param source volumic source
//      \param time present time
//    */
//    inline void updateSystem ( source_Type const& source )
//    {
//    	M_structuralOperator -> updateSystem(source);
//    }
//
//
//
//    //! Updates the system at the end of each time step given a source term
//    /*!
//      \param source volumic source
//      \param time present time
//    */
//    inline void updateSourceTerm ( source_Type const& source )
//    {
//    	M_structuralOperator -> updateSourceTerm(source);
//    }
//
//
//    //! Updates the rhs at the start of each time step
//    /*!
//      \param rhs: solid  right hand side
//      !*/
//    inline void setRightHandSide (const vector_Type& rightHandSide)
//    {
//    	M_structuralOperator -> setRightHandSide(rightHandSide);
//    }
//
//    //! Comuptes the right hand side in the updateSystem methods
//    void computeRHSNoBC ( void );
//
//    //! Compute the mass matrix and it calls the method to build the linear part of the stiffness matrix of the material class
//    void buildSystem ( const Real coefficient );
//
//    //void buildSystem(matrix_Type & bigMatrixStokes, const Real& timeAdvanceCoefficient, const Real& factor); // used for monolithic
//
//    //! Compute the mass matrix and the linear part of the stiffness matrix
//    /*!
//      \param matrix the matrix containing the mass matrix and the linear part of he stiffness matrix
//      \param rescale factor for FSI problems
//    */
//    void computeMassMatrix ( const Real factor = 1.);
//
//    //! Solve the non-linear system
//    /*!
//      \param bch BCHander object containing the boundary conditions
//    */
//    void iterate ( const bcHandler_Type& bch );
//
//    //! Solve the linearized problem. Used in FSI segregated in ExactJacobian
//    /*!
//      \param bch BCHander object containing the boundary conditions
//    */
//    void iterateLin ( bcHandler_Type& bch );
//
//
//    //! Output
//    /*!
//      \param c output file
//    */
//    void showMe ( std::ostream& c = std::cout ) const;
//
//    //! Update the Jacobian Matrix at each iteration of the nonLinearRichardson method
//    /*!
//      \param solution the current solution at each iteration of the nonLinearRichardson method
//      \param jacobian the Jacobian matrix that must be updated
//    */
//    void updateJacobian ( const vector_Type& solution, matrixPtr_Type& jacobian );
//
//    //! Update the Jacobian Matrix at each iteration of the nonLinearRichardson method
//    /*! Note: this method is used in FSIExactJacobian
//      \param solution the current solution at each iteration of the nonLinearRichardson method
//    */
//    void updateJacobian (const vector_Type& solution );
//
//    //! Solves the tangent problem for newton iterations
//    /*!
//      \param step the vector containing the solution of the sistem J*step=-Res
//      \param res the vector conteining the residual
//      \param lin_res_tol linear_rel_tol send for the relative tolerance to the linear solver is therefore eta.
//      eta is determined by the modified Eisenstat-Walker formula
//    */
//    void solveJac ( vector_Type&       step,
//                    const vector_Type& residual,
//                    Real&            linear_rel_tol ) ;
//    //    void solveJac( const Vector& res, Real& linear_rel_tol, Vector &step);
//
//    //! Solves the tangent problem with custom BC
//    /*!
//      \param step the vector containing the solution of the sistem J*step=-Res
//      \param res the vector conteining the residual
//      \param lin_res_tol linear_rel_tol send for the relative tolerance to the linear solver is therefore eta.
//      eta is determined by the modified Eisenstat-Walker formula
//      \param BCd BCHandler object containing the boundary condition
//    */
//    void solveJacobian ( vector_Type&       step,
//                         const vector_Type& residual,
//                         Real&            linear_rel_tol,
//                         bcHandler_Type&    BCd ) ;
//
//
//    //! Evaluates residual for newton interations
//    /*!
//      \param res residal vector that is update every time the method is called
//      \param sol solution vector from which the residual is computed
//      \param iter iteration of the nonLinearRichardson method
//    */
//    void evalResidual ( vector_Type& residual, const vector_Type& solution, Int iter);
//
//    //! Evaluates residual of the displacement for FSI problems
//    /*!
//      \param sol, the current displacement of he sturcture
//    */
//    void evalResidualDisplacement ( const vector_Type& solution );
//
//    //! Evaluates residual of the displacement in the Linearized problem of ExactJcobian. FSI problems
//    /*!
//      \param sol, the current displacement of he sturcture
//    */
//    void evalResidualDisplacementLin ( const vector_Type& solution );
//
//    //! Sets the initial displacement, velocity, acceleration
//    /*!
//      \param d0 space function describing the initial displacement
//      \param w0 space function describing the initial velocity
//      \param a0 space function describing the initial acceleration
//    */
//    void initialize ( const function& d0 );
//
//    //! Sets the initial displacement, velocity, acceleration
//    /*!
//      \param d0 space function describing the initial displacement
//      \param w0 empty vector
//      \param a0 empty vector
//    */
//    void initialize ( vectorPtr_Type d0 );
//
//    //! Computes the velocity and acceleration vector at the n-th time step
//    //void updateVelAndAcceleration();
//
//    //! Reduce the complete solution to the solution on the pocessor with rank 0
//    /*!
//      \param disp displacement solution
//      \param vel velocity solution
//    */
//    void reduceSolution ( Vector& displacement, Vector& velocity );
//
//    //! Multiply the mass matrix and the linear stiffness matrix by the rescaleFactor
//    //  void rescaleMatrices(); // used for monolithic
//
//    /**
//       in the linear case the solid matrix is constant, thus it does not need to be recomputed.
//    */
//
//    //! Update (in the case of nonlinear material) the solid matrix
//    /*!
//      \param stiff stiffness matrix
//      \param sol the current solution
//      \param factor the rescaleFactor
//    */
//    void computeMatrix ( matrixPtr_Type& stiff, const vector_Type& sol, Real const& factor );
//
//
//#ifdef COMPUTATION_JACOBIAN
//    //! compute the value of the determinant of F in all the volumes of the mesh
//    /*!
//      \param displacement the solution at a certain time
//      \return the vector with the values for J
//    */
//    void jacobianDistribution ( vectorPtr_Type displacement, vector_Type& jacobianDistribution );
//#endif
//
//
//#ifdef COLORING_MESH
//    //! compute the value of the determinant of F in all the volumes of the mesh
//    /*!
//      \param displacement the solution at a certain time
//      \return the vector with the values for J
//    */
//    void colorMesh ( vector_Type& meshColors );
//#endif
//
//    //void updateMatrix(matrix_Type & bigMatrixStokes);// used for monolithic
//    //void updateCoupling(matrix_Type couplingMatrix);// used for monolithic
//
//    //@}
//
//    //! @name Set Methods
//    //@{
//
//    //!Setters
//    //! Set the BCHandler object
//    void setBC (const bcHandler_Type& BCd)
//    {
//        M_BCh = BCd;
//    }
//
//    //! Set the source object
//    void setSourceTerm ( source_Type const& s )
//    {
//        M_source = s;
//    }
//
//    // //! Set the preconditioner
//    // void resetPrec(bool reset = true) { if (reset) M_linearSolver.precReset(); }
//
//    // //! Set the displacement
//    // virtual void setDisp(const vector_Type& disp) {*M_disp = disp;} // used for monolithic
//
//    //! Set the recur parameter
//    void setRecur (UInt recur)
//    {
//        M_recur = recur;
//    }
//
//    //! Set the data fields with the Getpot data file for preconditioners and solver
//    void setDataFromGetPot ( const GetPot& dataFile );
//
//    void setTimeAdvance ( const timeAdvancePtr_Type& timeAdvancePtr )
//    {
//        M_timeAdvance = timeAdvancePtr;
//    }
//
//    //@}
//
//
//    //! @name Get Methods
//    //@{
//
//    //! Getters
//    //! Get the Epetramap
//    MapEpetra   const& map()       const
//    {
//        return *M_localMap;
//    }
//
//    //! Get the Displayer object
//    Displayer   const& displayer() const
//    {
//        return *M_Displayer;
//    }
//
//    boost::shared_ptr<const Displayer>   const& displayerPtr() const
//    {
//        return M_Displayer;
//    }
//
//    //! Get the matrix containing the mass mtrix and the linear part of the stiffness matrix
//    //matrixPtr_Type const MassStiff() const {return M_massStiff; }
//
//    //! Get the mass matrix
//    matrixPtr_Type const massMatrix() const
//    {
//        return M_massMatrix;
//    }
//
//    //! Get the FESpace object
//    FESpace_Type& dispFESpace()
//    {
//        return *M_dispFESpace;
//    }
//    FESpacePtr_Type dispFESpacePtr()
//    {
//        return M_dispFESpace;
//    }
//
//    //! Get the ETFESpace object
//    ETFESpace_Type& dispETFESpace()
//    {
//        return *M_dispETFESpace;
//    }
//
//    ETFESpacePtr_Type dispETFESpacePtr()
//    {
//        return M_dispETFESpace;
//    }
//
//    //! Get the bCHandler object
//    bcHandler_Type const& bcHandler() const
//    {
//        return M_BCh;
//    }
//
//    //! Get the residual
//    vector_Type& residual()
//    {
//        return *M_residual_d;
//    }
//
//    //! Get the source term
//    source_Type const& sourceTerm() const
//    {
//        return M_source;
//    }
//
//    //! Get the displacement
//    vector_Type& displacement()
//    {
//        return *M_disp;
//    }
//
//    vectorPtr_Type displacementPtr()
//    {
//        return M_disp;
//    }
//
//    //! Get the right hand sde without BC
//    vectorPtr_Type& rhsWithoutBC()
//    {
//        return M_rhsNoBC;
//    }
//
//    //! Get the right hand. The member rhsCopy is used for Debug purposes!
//    vector_Type& rhsCopy()
//    {
//        return *M_rhsCopy;
//    }
//    vector_Type& residualCopy()
//    {
//        return *M_residualCopy;
//    }
//
//    //! Get the comunicator object
//    boost::shared_ptr<Epetra_Comm> const& comunicator() const
//    {
//        return M_Displayer->comm();
//    }
//
//    //! Get the rescaleFactor
//    Real rescaleFactor()
//    {
//        return M_rescaleFactor;
//    }
//
//    /*! Get the offset parameter. It is taken into account when the boundary conditions
//      are applied and the matrices are assembled.
//    */
//    const UInt& offset() const
//    {
//        return M_offset;
//    }
//
//    /*! Get the offset parameter. It is taken into account when the boundary conditions
//      are applied and the matrices are assembled.
//    */
//    virtual const materialPtr_Type& material() const
//    {
//        return M_material;
//    }
//
//    /**
//       Do nothing in the linear case: the matrix remains constant. Otherwise substitute the matrix with an updated one
//    */
//    //! Get the Solid Matrix
//    void solidMatrix ( matrixPtr_Type& /*matrix*/ )
//    {
//    }
//
//    // Physic constant
//    //! Get the thickness
//    Real thickness() const
//    {
//        return M_data->thickness();
//    }
//
//    //! Get the Young modulus
//    Real young ( UInt material = 1)            const
//    {
//        return M_data->young ( material );
//    }
//
//    //! Get the Poisson coefficient
//    Real poisson ( UInt material = 1 )          const
//    {
//        return M_data->poisson ( material );
//    }
//
//    //! Get the density
//    Real rho()       const
//    {
//        return M_data->rho();
//    }
//
//    //! Get the data container
//    const boost::shared_ptr<data_Type>& data() const
//    {
//        return M_data;
//    }
//
//    void apply ( const vector_Type& sol, vector_Type& res) const;
//
//    //! Get the density
//    mapMarkerVolumesPtr_Type mapMarkersVolumes() const
//    {
//        return M_mapMarkersVolumes;
//    }
//
//    //! Get the density
//    mapMarkerIndexesPtr_Type mapMarkersIndexes() const
//    {
//        return M_mapMarkersIndexes;
//    }
//
//    const timeAdvancePtr_Type& timeAdvancePtr() const
//    {
//        return M_timeAdvance;
//    }
//
//    inline meshPtr_Type mesh() const
//    {
//        return M_dispFESpace -> mesh();
//    }
//
//    inline Real res1()
//    {
//        return M_res1;
//    }
//    //@}
//
//protected:
//
//    //! Apply boundary condition
//    /*!
//      \param matrix the matrix of the system
//      \param rhs the right hand side of the system
//      \param BCh BCHandler object
//      \param offset the offset parameter
//    */
//    void applyBoundaryConditions (matrix_Type& matrix,
//                                  vector_Type& rhs,
//                                  bcHandler_Type& BCh,
//                                  UInt         offset = 0);
//
//
//    UInt dim() const
//    {
//        return M_dispFESpace->dim();
//    }
//
//
//    //! construct the map between the markers and the volumes
//    /*!
//      \param VOID
//      \return VOID
//    */
//    void setupMapMarkersVolumes ( void );
//
//
//    //@}
//
    //!@name Methods
    //@{

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
protected:

    //! Material class
//    structuralOperatorPtr_Type 			 M_structuralOperator;
    materialPtr_Type                     M_EMMaterial;

};

//====================================
// Constructor
//=====================================

template <typename Mesh>
EMStructuralOperator<Mesh>::EMStructuralOperator( ) :
//    M_structuralOperator(),
super(),
M_EMMaterial()
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
   M_EMMaterial = dynamic_pointer_cast<material_Type> (this -> material());
}




}
#endif
