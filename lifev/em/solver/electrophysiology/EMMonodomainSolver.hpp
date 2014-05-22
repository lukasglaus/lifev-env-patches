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
 @brief Class for solving the Monodomain model in electrophysiology.

 @date 02 - 2013
 @author Simone Rossi <simone.rossi@epfl.ch>

 @last update 04 - 2014

 This class provides interfaces to solve the monodomain equation
 ( reaction diffusion equation ) using the ETA framework.
 The solution can be performed using three different methods:
 -operator splitting method (at this point available only with forward Euler
 for the reaction step and backward Euler for the diffusion step.
 Second order splitting is available but still experimental. );
 -Ionic Currents Interpolation (at this point only forward Euler);
 -State Variable interpolation (at this point only forward Euler).
 */

#ifndef _EMMONODOMAINSOLVER_H_
#define _EMMONODOMAINSOLVER_H_

#include <lifev/electrophysiology/solver/ElectroETAMonodomainSolver.hpp>
#include <lifev/core/fem/GradientRecovery.hpp>

namespace LifeV
{

//! monodomainSolver - Class featuring the solver for monodomain equations

template<typename Mesh, typename EMIonicModel>
class EMMonodomainSolver : public virtual ElectroETAMonodomainSolver<Mesh, EMIonicModel>
{

    //!Monodomain Solver
    /*!
     The monodomain equation reads
     \f \Chi

     */

public:

    //! @name Type definitions
    //@{

	//! Mesh
    typedef Mesh 														mesh_Type;

    typedef boost::shared_ptr<mesh_Type>                    			meshPtr_Type;

    //! Distributed vector // For parallel usage
    typedef VectorEpetra 												vector_Type;

    typedef boost::shared_ptr<VectorEpetra> 							vectorPtr_Type;

    typedef std::vector<vectorPtr_Type> 								vectorOfPtr_Type;

    //! Distributed Matrix // For parallel usage
    typedef MatrixEpetra<Real> 											matrix_Type;

    typedef boost::shared_ptr<matrix_Type> 								matrixPtr_Type;

    //! Communicator to exchange informations among processes
    typedef Epetra_Comm 												comm_Type;

    typedef boost::shared_ptr<comm_Type> 								commPtr_Type;

    //! Expression template  scalar finite element space
    //! To be used in the expression assembly namespace
    typedef ETFESpace<mesh_Type, MapEpetra, 3, 1> 						ETFESpace_Type;

    typedef boost::shared_ptr<ETFESpace<mesh_Type, MapEpetra, 3, 1> > 	ETFESpacePtr_Type;

    //! Expression template vectorial finite element space
    //! To be used in the expression assembly namespace
    typedef ETFESpace<mesh_Type, MapEpetra, 3, 3> 						ETFESpaceVectorial_Type;

    typedef boost::shared_ptr<ETFESpaceVectorial_Type> 					ETFESpaceVectorialPtr_Type;

    //! Finite element space
    typedef FESpace<mesh_Type, MapEpetra> 								feSpace_Type;

    typedef boost::shared_ptr<feSpace_Type> 							feSpacePtr_Type;

    //! Linear Solver
    typedef LinearSolver 												linearSolver_Type;

    typedef boost::shared_ptr<LinearSolver> 							linearSolverPtr_Type;

    //! Exporter to save the solution
    typedef Exporter<mesh_Type> 										IOFile_Type;

    typedef boost::shared_ptr<IOFile_Type> 								IOFilePtr_Type;

    //! Exporter data
    typedef ExporterData<mesh_Type> 									IOData_Type;

    typedef ExporterEnsight<mesh_Type> 									ensightIOFile_Type;

#ifdef HAVE_HDF5
    typedef ExporterHDF5< mesh_Type > 									hdf5IOFile_Type;
#endif

    //! Preconditioner
    typedef LifeV::Preconditioner 										basePrec_Type;

    typedef boost::shared_ptr<basePrec_Type> 							basePrecPtr_Type;

    //! MultiLevel Preconditioner
    typedef LifeV::PreconditionerML 									prec_Type;

    typedef boost::shared_ptr<prec_Type> 								precPtr_Type;

    //! Ionic model
    typedef EMIonicModel 													ionicModel_Type;

    //! Base class of the ionic model
    typedef ElectroIonicModel 											superIonicModel;

    typedef boost::shared_ptr<ionicModel_Type> 							ionicModelPtr_Type;

    //! xml list to read parameters
    typedef Teuchos::ParameterList 										list_Type;

    //! boost function
    typedef boost::function <Real (const Real& t,
    							   const Real& x,
    							   const Real& y,
    							   const Real& z,
    							   const ID&   i) > 					function_Type;

    //! 3x3 matrix
    typedef MatrixSmall<3, 3> 											matrixSmall_Type;

    typedef ElectroETAMonodomainSolver<Mesh, EMIonicModel>              super;
    //@}

    //! @name Constructors & Destructor
    //@{

    //!Empty Constructor
    /*!
     */
    EMMonodomainSolver();
    
        //! Constructor
    /*!
     * @param GetPot datafile (for preconditioner)
     * @param boost::shared_ptr<IonicModel>  chosen ionic model pointer
     * @param boost::shared_ptr<Mesh> Pointer to the partitioned mesh
     */
    EMMonodomainSolver (GetPot&            dataFile,
    		            ionicModelPtr_Type model,
                        meshPtr_Type       meshPtr);
    //! Constructor
    /*!
     * @param meshName file name of the mesh
     * @param meshPath path to the mesh
     * @param datafile  GetPot file to setup the preconditioner
     * @param model shared pointer to the chosen ionic model
     */
    EMMonodomainSolver (std::string        meshName,
    		            std::string        meshPath,
                        GetPot&            dataFile,
                        ionicModelPtr_Type model);

    //! Constructor
    /*!
     * @param string file name of the mesh
     * @param string path to the mesh
     * @param GetPot datafile (for preconditioner)
     * @param boost::shared_ptr<IonicModel> chosen ionic model pointer
     * @param boost::shared_ptr<Epetra_Comm> Epetra communicator
     */
    EMMonodomainSolver (std::string        meshName,
    		            std::string        meshPath,
                        GetPot&            dataFile,
                        ionicModelPtr_Type model,
                        commPtr_Type       comm);
    //! Copy Constructor
    /*!
     * @param ElectroETAmonodomainSolver object
     */
    EMMonodomainSolver (const EMMonodomainSolver& solver);

    //!Operator=()
    /*!
     * @param ElectroETAmonodomainSolver object
     */
    EMMonodomainSolver<Mesh, EMIonicModel>& operator= ( const EMMonodomainSolver& solver);

    //! Destructor
    virtual ~EMMonodomainSolver() {}
    //@}

    //! @name Get Methods
    //@{

    //! get the pointer to the transmembrane potential
    inline vectorPtr_Type displacementPtr() const
    {
        return M_displacementPtr;
    }

    //! getter of the ET FESpace associated with the displacement
    inline ETFESpaceVectorialPtr_Type displacementETFESpacePtr() const
    {
        return M_displacementETFESpacePtr;
    }


    //! set the pointer to the displacement for mechanical feedback
    /*!
     @param displacementPtr pointer to the displacement vector
     */
    inline void setDisplacementPtr (const vectorPtr_Type displacementPtr)
    {
        this->M_displacementPtr = displacementPtr;
    }

    //! set the displacement for mechanical feedback
    /*!
     @param displacement the displacement vector
     */
    inline void setDisplacement (const vector_Type& displacement)
    {
        (* (M_displacementPtr) ) = displacement;
    }


    //@}

    //! @name Methods
    //@{

    //! setup method used in the constructor
    /*!
     @param dataFile needed to set up the preconditioner
     @param ionicSize number of equation in the ionic model
     */
//    void setup (GetPot& dataFile, short int ionicSize);
//
//
//    void setup (std::string meshName, std::string meshPath, GetPot& dataFile,
//                short int ionicSize);
    //! create mass matrix
    /*!
     * Computes the mass matrix calling different methods if the mass should be lumped
     * and if the mechanical feedback (the displacement vector) is set.
     */
    void setupMassMatrix(); //! create mass matrix

    //! create mass matrix
    /*!
     * Computes the mass matrix calling with the mechanical feedback (the displacement vector).
     */
    void setupMassMatrixWithMehcanicalFeedback ();

    //! create mass matrix
    /*!
     * Computes the lumped mass matrix by nodal integration and
     * with the mechanical feedback (the displacement vector).
     @param disp vector of the displacement of the domain (for electromechanical simulations)
     */
    void setupLumpedMassMatrixWithMehcanicalFeedback ();

    //! create stiffness matrix
    /*!
      * Computes the stiffness matrix calling different methods
      * if the mechanical feedback (the displacement vector) is set.
      */
    void setupStiffnessMatrix();

    //! create stiffness matrix in a moving domain
    /*!
     * Computes the stiffenes matrix
     * with the mechanical feedback (the displacement vector).
     @param disp vector of the displacement of the domain (for electromechanical simulations)
     */
    void setupStiffnessMatrixWithMehcanicalFeedback ();

    //! setup all the matrices of the system
    /*!
     * Computes all the matrices
     * with the mechanical feedback (the displacement vector) if needed.
     */
    void setupMatrices();

    //! update all the matrices of the system
    /*!
     * Computes all the matrices
     * with the mechanical feedback (the displacement vector).
     @param disp vector of the displacement of the domain (for electromechanical simulations)
     */
    void updateMatrices();

//
//
//
//
//
//    //! Solves one diffusion step using the backward Euler scheme
//    /*!
//     * \f[
//     * A\mathbf{V}^{n+1} = \left( C_m \frac{M}{\Delta t} + K(\mathbf{f}) \right)\mathbf{V}^{n+1} =\frac{M}{\Delta t} \mathbf{V}^*.
//     * \f]
//     */
//    void solveOneDiffusionStepBE();
//
//
//    //! add to a given exporter the pointer to the potential and to the gating variables saved with name fileName
//    /*!
//     @param exporter where you want to save the solution
//     @param fileName name of the file we wish to export
//     @param folder directory where to save the solution
//     */
//    void setupExporter (IOFile_Type& exporter, std::string fileName = "output",
//                        std::string folder = "./");
//
//
//    //! Solves the gating variables with forward Euler
//    void solveOneStepGatingVariablesFE();
//
//    //! Solves the gating variables with Rush-Larsen scheme
//    void solveOneStepGatingVariablesRL();
//
//    //! Compute the rhs using state variable interpolation
//    void computeRhsSVI();
//
    //! Compute the rhs using ionic current interpolation
    //void computeRhsICI();

//    //! Compute the rhs using ionic current interpolation
//    /*!
//     * This method is useful to solve ICI without lumping the mass matrix
//     * in fron of the reaction term.
//     * Lump the mass matrix, and pass as argument a full mass matrix
//     */
//    void computeRhsICIWithFullMass ();
//
    //!Solve one full step with ionic current interpolation
    /*!
     * \f[
     * A\mathbf{V}^{n+1} = \left( \frac{M}{\Delta t} + K(\mathbf{f}) \right)\mathbf{V}^{n+1} =\frac{M}{\Delta t} \mathbf{V}^n+M\mathbf{I},
     * \f]
     * where $\mathbf{I}$ is the vector of the ionic currents $I_j = I_{ion}(V_j^n)$
     */
    void solveOneICIStep();

    virtual void setup( GetPot& dataFile, short int ionicModelSize);
//
//    //! solves using ionic current interpolation
//    /*!
//     * This method is useful to solve ICI without lumping the mass matrix
//     * in fron of the reaction term.
//     * Lump the mass matrix, and pass as argument a full mass matrix
//	@param mass mass matrix
//     */
//    void solveOneICIStep (matrix_Type& mass);
//
//    //!Solve one full step with state variable interpolation
//    /*!
//     * \f[
//     * A\mathbf{V}^{n+1} = \left( \frac{M}{\Delta t} + K(\mathbf{f}) \right)\mathbf{V}^{n+1} =\frac{M}{\Delta t} \mathbf{V}^n+\mathbf{I}_{ion}(\mathbf{V}^n).
//     * \f]
//     */
//    void solveOneSVIStep();
    inline bool oneWayCoupling()
    {
    	std::cout << "\nI'm " << M_oneWayCoupling << "\n";
    	return M_oneWayCoupling;
    }

    inline bool mechanicsModifiesConductivity()
    {
    	return M_mechanicsModifiesConductivity;
    }

    inline void setOneWayCoupling(bool oneWay)
    {
    	M_oneWayCoupling = oneWay;
    }

    inline void setMechanicsModifiesConductivity(bool modifiesConductivity)
    {
    	M_mechanicsModifiesConductivity = modifiesConductivity;
    }
    //@}

private:
    //pointer to the displacement
    vectorPtr_Type M_displacementPtr;

    ETFESpaceVectorialPtr_Type M_displacementETFESpacePtr;

    //Defines if you want to consider the mechanical feedback
    bool M_oneWayCoupling;
    //true if the mechanical feedback changes the conductivity tensor
    bool M_mechanicsModifiesConductivity;

};
// class MonodomainSolver

//
// IMPLEMENTATION
//
// ===================================================
//! Constructors
// ===================================================
//!Empty constructor
template<typename Mesh, typename IonicModel>
EMMonodomainSolver<Mesh, IonicModel>::EMMonodomainSolver() :
    super(),
    M_displacementETFESpacePtr(),
    M_oneWayCoupling(false),
    M_mechanicsModifiesConductivity(true)
{
//    M_oneWayCoupling = false;
//    M_mechanicsModifiesConductivity = true;
//    M_displacementETFESpacePtr.reset ( new ETFESpaceVectorial_Type (this->M_localMeshPtr,
//    		                                                       & (this->M_feSpacePtr -> refFE() ),
//    		                                                       this -> M_commPtr) );
}


template<typename Mesh, typename IonicModel>
EMMonodomainSolver<Mesh, IonicModel>::EMMonodomainSolver (std::string        meshName,
														  std::string        meshPath,
													      GetPot&            dataFile,
														  ionicModelPtr_Type model) :
	super                           (meshName, meshPath, dataFile, model),
	M_oneWayCoupling                (false),
	M_mechanicsModifiesConductivity (true)
{
    M_displacementETFESpacePtr.reset ( new ETFESpaceVectorial_Type (this->M_localMeshPtr,
    		                                                       & (this->M_feSpacePtr -> refFE() ),
    		                                                       this -> M_commPtr) );
}

template<typename Mesh, typename IonicModel>
EMMonodomainSolver<Mesh, IonicModel>::EMMonodomainSolver (GetPot&            dataFile,
														  ionicModelPtr_Type model,
														  meshPtr_Type       meshPtr) :
	super                           (dataFile, model, meshPtr),
	M_oneWayCoupling                (false),
	M_mechanicsModifiesConductivity (true)
{
    M_displacementETFESpacePtr.reset ( new ETFESpaceVectorial_Type (this->M_localMeshPtr,
    		                                                       & (this->M_feSpacePtr -> refFE() ),
    		                                                       this -> M_commPtr) );
}

template<typename Mesh, typename IonicModel>
EMMonodomainSolver<Mesh, IonicModel>::EMMonodomainSolver (std::string        meshName,
														  std::string        meshPath,
														  GetPot&            dataFile,
														  ionicModelPtr_Type model,
														  commPtr_Type       comm):
	super                           (meshName, meshPath, dataFile, model, comm),
	M_oneWayCoupling                (false),
	M_mechanicsModifiesConductivity (true)
{
    M_displacementETFESpacePtr.reset ( new ETFESpaceVectorial_Type (this->M_localMeshPtr,
    		                                                       & (this->M_feSpacePtr -> refFE() ),
    		                                                       this -> M_commPtr) );
}


//! Copy constructor
template<typename Mesh, typename IonicModel>
EMMonodomainSolver<Mesh, IonicModel>::EMMonodomainSolver (
    const EMMonodomainSolver& solver) :
    super::ElectroETAMonodomainSolver(solver)
{
    if(solver.M_displacementPtr)
    {
    	M_displacementPtr.reset(new vector_Type(*solver.M_displacementPtr));
        M_displacementETFESpacePtr.reset(new ETFESpaceVectorial_Type(*solver.ETFESpaceVectorialPtr_Type));
    }

    M_oneWayCoupling = solver.M_oneWayCoupling;
    M_mechanicsModifiesConductivity = solver.M_mechanicsModifiesConductivity;
}

//! Assignment operator
template<typename Mesh, typename IonicModel>
EMMonodomainSolver<Mesh, IonicModel>& EMMonodomainSolver < Mesh,
                           IonicModel >::operator= (const EMMonodomainSolver& solver)
{
	super::operator =(solver);
	if(solver.M_displacementPtr)
    {
    	M_displacementPtr.reset(new vector_Type(*solver.M_displacementPtr));
        M_displacementETFESpacePtr.reset(new ETFESpaceVectorial_Type(*solver.ETFESpaceVectorialPtr_Type));
    }

    M_oneWayCoupling = solver.M_oneWayCoupling;
    M_mechanicsModifiesConductivity = solver.M_mechanicsModifiesConductivity;

    return *this;
}

/********* SETUP METHODS */ //////

//template<typename Mesh, typename IonicModel>
//void EMMonodomainSolver<Mesh, IonicModel>::setup (GetPot& dataFile,
//                                                          short int ionicSize)
//{
//	super::setup(dataFile, ionicSize);
//
//    M_displacementETFESpacePtr.reset ( new ETFESpaceVectorial_Type (this->M_localMeshPtr,
//    		                                                       & (this->M_feSpacePtr -> refFE() ),
//    		                                                       this -> M_commPtr) );
//}
//
//template<typename Mesh, typename IonicModel>
//void EMMonodomainSolver<Mesh, IonicModel>::setup (std::string meshName, std::string meshPath, GetPot& dataFile,
//            short int ionicSize)
//{
//	super::setup(meshName, meshPath, dataFile, ionicSize);
//
//    M_displacementETFESpacePtr.reset ( new ETFESpaceVectorial_Type (this->M_localMeshPtr,
//    		                                                       & (this->M_feSpacePtr -> refFE() ),
//    		                                                       this -> M_commPtr) );
//}
///////////////////////////////////////////////////////////
//                MATRICES SETUP
///////////////////////////////////////////////////////////
template<typename Mesh, typename IonicModel>
void EMMonodomainSolver<Mesh, IonicModel>::setup( GetPot& dataFile, short int ionicModelSize)
{
    super::setup(dataFile, ionicModelSize);
    M_displacementETFESpacePtr.reset ( new ETFESpaceVectorial_Type (this->M_localMeshPtr,
    		                                                       & (this->M_feSpacePtr -> refFE() ),
    		                                                       this -> M_commPtr) );
}

///////////////////////////////////////////////////////////
//                MASS MATRIX
///////////////////////////////////////////////////////////
//!setup mass matrix
// choose between the one with mechanical feedback
// or the one from the monodomain model
template<typename Mesh, typename IonicModel>
void EMMonodomainSolver<Mesh, IonicModel>::setupMassMatrix()
{
	if(M_displacementPtr && !M_oneWayCoupling)
	{
	    if (this->M_lumpedMassMatrix)
	    {
	    	setupLumpedMassMatrixWithMehcanicalFeedback ();
	    }
	    else
	    {
            setupMassMatrixWithMehcanicalFeedback ();
	    }
	}
	else super::setupMassMatrix();
}

//setup mass matrix with mechanical fedback
template<typename Mesh, typename IonicModel>
void EMMonodomainSolver<Mesh, IonicModel>::setupMassMatrixWithMehcanicalFeedback ()
{
    if (this->M_verbose && this->M_commPtr->MyPID() == 0)
    {
        std::cout << "\nEM Monodomain Solver: Setting up mass matrix with coupling with mechanics ";
    }

    *this->M_massMatrixPtr *= 0.0;
    {
        using namespace ExpressionAssembly;

        if(M_displacementETFESpacePtr)std::cout <<"Disp ETFE space available\n";
        if(M_displacementPtr)std::cout <<"Disp available\n";
        if(super::M_localMeshPtr)std::cout <<"Local Mesh available\n";
        if(super::M_massMatrixPtr)std::cout <<"Mass Matrix available\n";
        BOOST_AUTO_TPL (I, value (super::M_identity) );
        BOOST_AUTO_TPL (Grad_u, grad (M_displacementETFESpacePtr, *M_displacementPtr) );
        BOOST_AUTO_TPL (F, (Grad_u + I) );
        BOOST_AUTO_TPL (J, det (F) );

        integrate (elements (super::M_localMeshPtr),
        		   this->M_feSpacePtr->qr(),
        		   this->M_ETFESpacePtr,
        		   this->M_ETFESpacePtr,
                   J * phi_i * phi_j ) >> this->M_massMatrixPtr;

    }
    this->M_massMatrixPtr->globalAssemble();
}

//setup lumped mass matrix with mechanical fedback
template<typename Mesh, typename IonicModel>
void EMMonodomainSolver<Mesh, IonicModel>::setupLumpedMassMatrixWithMehcanicalFeedback ()
{

    *this->M_massMatrixPtr *= 0.0;

    vectorPtr_Type dUdx (new vector_Type (M_displacementPtr->map() ) );
    vectorPtr_Type dUdy (new vector_Type (M_displacementPtr->map() ) );
    vectorPtr_Type dUdz (new vector_Type (M_displacementPtr->map() ) );

    *dUdx = GradientRecovery::ZZGradient (M_displacementETFESpacePtr, *M_displacementPtr, 0);
    *dUdy = GradientRecovery::ZZGradient (M_displacementETFESpacePtr, *M_displacementPtr, 1);
    *dUdz = GradientRecovery::ZZGradient (M_displacementETFESpacePtr, *M_displacementPtr, 2);

    vectorPtr_Type J (new vector_Type (this->M_potentialPtr->map() ) );
    int n = J->epetraVector().MyLength();
    int i (0);
    int j (0);
    int k (0);
    for (int p (0); p < n; p++)
    {
        i = dUdx->blockMap().GID (p);
        j = dUdx->blockMap().GID (p + n);
        k = dUdx->blockMap().GID (p + 2 * n);

        Real F11 = 1.0 + (*dUdx) [i];
        Real F12 =       (*dUdy) [i];
        Real F13 =       (*dUdz) [i];
        Real F21 =       (*dUdx) [j];
        Real F22 = 1.0 + (*dUdy) [j];
        Real F23 =       (*dUdz) [j];
        Real F31 =       (*dUdx) [k];
        Real F32 =       (*dUdy) [k];
        Real F33 = 1.0 + (*dUdz) [k];

        (*J) [i] = F11 * ( F22 * F33 - F32 * F23 )
                   - F12 * ( F21 * F33 - F31 * F23 )
                   + F13 * ( F21 * F32 - F31 * F22 );
    }

    if (this->M_verbose && this->M_localMeshPtr->comm()->MyPID() == 0)
    {
        std::cout << "\nEM Monodomain Solver: Setting up lumped mass matrix coupling with mechanics";
    }

    {
        using namespace ExpressionAssembly;

        integrate ( elements (this->M_localMeshPtr),
        		    quadRuleTetra4ptNodal,
        		    this->M_ETFESpacePtr,
        		    this->M_ETFESpacePtr,
                    value (this->M_ETFESpacePtr, *J) * phi_i * phi_j
                    ) >> this->M_massMatrixPtr;

    }

    this->M_massMatrixPtr->globalAssemble();

}

///////////////////////////////////////////////////////////
//                STIFFNESS MATRIX
///////////////////////////////////////////////////////////

//setup the stiffness matrix checking if we have mechano electrical feedback
template<typename Mesh, typename IonicModel>
void EMMonodomainSolver<Mesh, IonicModel>::setupStiffnessMatrix()
{

    if (M_displacementPtr && !M_oneWayCoupling && M_mechanicsModifiesConductivity)
    {
    	std::cout << "=========================\n";
    	std::cout << "Using mechanical feedback\n";
    	std::cout << "=========================\n";
//        setupStiffnessMatrix (M_displacementPtr);
        setupStiffnessMatrixWithMehcanicalFeedback ();
    }
    else
    {
        super::setupStiffnessMatrix();
    }
}

//setup the stiffness matrix modifying the conductivity
template<typename Mesh, typename IonicModel>
void EMMonodomainSolver<Mesh, IonicModel>::setupStiffnessMatrixWithMehcanicalFeedback ()
{
    if (this->M_verbose && this->M_localMeshPtr->comm()->MyPID() == 0)
    {
        std::cout
                << "\nETA Monodomain Solver: Setting up stiffness matrix  coupling with mechanics";
    }

    *this->M_stiffnessMatrixPtr *= 0.0;
    Real sigmal = this->M_diffusionTensor[0];
    Real sigmat = this->M_diffusionTensor[1];

    {
        using namespace ExpressionAssembly;

        BOOST_AUTO_TPL (I, value (this->M_identity) );
        BOOST_AUTO_TPL (Grad_u, grad (M_displacementETFESpacePtr, *M_displacementPtr) );
        BOOST_AUTO_TPL (F, (Grad_u + I) );
        BOOST_AUTO_TPL (FmT, minusT (F) );
        BOOST_AUTO_TPL (Fm1, transpose (FmT) );
        BOOST_AUTO_TPL (J, det (F) );
        BOOST_AUTO_TPL (Jm23, pow (J, -2. / 3) );
        BOOST_AUTO_TPL (f0, value (M_displacementETFESpacePtr, *this->M_fiberPtr) );
        BOOST_AUTO_TPL (D,
                        value (sigmat) * I
                        + (value (sigmal) - value (sigmat) )
                        * outerProduct (f0, f0) );

        integrate ( elements (this->M_localMeshPtr),
        		    this->M_feSpacePtr->qr(),
        		    this->M_ETFESpacePtr,
        		    this->M_ETFESpacePtr,
                    dot (J * Fm1 * D * FmT * grad (phi_i), grad (phi_j) )
                    )
                    >> this->M_stiffnessMatrixPtr;

    }

    this->M_stiffnessMatrixPtr->globalAssemble();

}

//update matrices
template<typename Mesh, typename IonicModel>
void EMMonodomainSolver<Mesh, IonicModel>::setupMatrices()
{
    	setupMassMatrix();
    	setupStiffnessMatrix();
    	super::setupGlobalMatrix();
}


//update matrices
template<typename Mesh, typename IonicModel>
void EMMonodomainSolver<Mesh, IonicModel>::updateMatrices()
{
    if (M_displacementPtr && !M_oneWayCoupling && M_mechanicsModifiesConductivity)
    {
    	setupMassMatrixWithMehcanicalFeedback();
    	setupStiffnessMatrixWithMehcanicalFeedback();
    	super::setupGlobalMatrix();
    }
}
//
///********* SOLVING METHODS */    ////////////////////////
//
//
//
//template<typename Mesh, typename IonicModel>
//void EMMonodomainSolver<Mesh, IonicModel>::solveOneDiffusionStepBE()
//{
//    if (M_displacementPtr)
//    {
//        M_linearSolverPtr->setOperator (M_globalMatrixPtr);
//    }
//    M_linearSolverPtr->setRightHandSide (M_rhsPtrUnique);
//    M_linearSolverPtr->solve (M_potentialPtr);
//}
//
//template<typename Mesh, typename IonicModel>
//void EMMonodomainSolver<Mesh, IonicModel>::solveOneSplittingStep()
//{
//    solveOneReactionStepFE();
//    (*M_rhsPtrUnique) *= 0;
//    updateRhs();
//    solveOneDiffusionStepBE();
//}
//
//template<typename Mesh, typename IonicModel>
//void EMMonodomainSolver<Mesh, IonicModel>::solveOneSplittingStep (
//    IOFile_Type& exporter, Real t)
//{
//    solveOneSplittingStep();
//    exportSolution (exporter, t);
//}
//
//template<typename Mesh, typename IonicModel>
//void EMMonodomainSolver<Mesh, IonicModel>::solveSplitting()
//{
//    for (Real t = M_initialTime; t < M_endTime;)
//    {
//        t = t + M_timeStep;
//        solveOneSplittingStep();
//    }
//}
//
//template<typename Mesh, typename IonicModel>
//void EMMonodomainSolver<Mesh, IonicModel>::solveSplitting (
//    IOFile_Type& exporter)
//{
//    if (M_endTime > M_timeStep)
//    {
//        for (Real t = M_initialTime; t < M_endTime;)
//        {
//            t = t + M_timeStep;
//            solveOneSplittingStep (exporter, t);
//        }
//    }
//}
//
//
//template<typename Mesh, typename IonicModel>
//void EMMonodomainSolver<Mesh, IonicModel>::solveOneStepGatingVariablesFE()
//{
//    M_ionicModelPtr->superIonicModel::computeGatingRhs (M_globalSolution,
//                                                        M_globalRhs);
//
//    for (int i = 1; i < M_ionicModelPtr->Size(); i++)
//    {
//        * (M_globalSolution.at (i) ) = * (M_globalSolution.at (i) )
//                                       + M_timeStep * (* (M_globalRhs.at (i) ) );
//    }
//}
//template<typename Mesh, typename IonicModel>
//void EMMonodomainSolver<Mesh, IonicModel>::solveOneStepGatingVariablesRL()
//{
//
//    M_ionicModelPtr->superIonicModel::computeGatingVariablesWithRushLarsen (
//        M_globalSolution, M_timeStep);
//    M_ionicModelPtr->superIonicModel::computeNonGatingRhs (M_globalSolution,
//                                                           M_globalRhs);
//    int offset = M_ionicModelPtr->numberOfGatingVariables() + 1;
//    for (int i = offset; i < M_ionicModelPtr->Size(); i++)
//    {
//    	*(M_globalRhs[i]) *= M_timeStep;
//    	* (M_globalSolution[i]) += *(M_globalRhs[i]);
//    }
//}
//
//template<typename Mesh, typename IonicModel>
//void EMMonodomainSolver<Mesh, IonicModel>::computeRhsICI()
//{
//    M_ionicModelPtr->superIonicModel::computePotentialRhsICI (M_globalSolution,
//                                                              M_globalRhs, (*M_massMatrixPtr) );
//    updateRhs();
//}

//template<typename Mesh, typename IonicModel>
//void EMMonodomainSolver<Mesh, IonicModel>::computeRhsICIWithFullMass ()
//{
//    M_ionicModelPtr->superIonicModel::computePotentialRhsICI (M_globalSolution,
//                                                              M_globalRhs, *M_fullMassMatrixPtr);
//    updateRhs();
//}
//
//template<typename Mesh, typename IonicModel>
//void EMMonodomainSolver<Mesh, IonicModel>::computeRhsSVI()
//{
//    if (M_displacementPtr)
//    {
//        if (M_verbose && M_commPtr -> MyPID() == 0)
//        {
//            std::cout << "\nETA Monodomain Solver: updating rhs with SVI with mechanical coupling\n";
//        }
//        boost::shared_ptr<FESpace<mesh_Type, MapEpetra> > vectorialSpace (
//            new FESpace<mesh_Type, MapEpetra> (M_localMeshPtr, M_elementsOrder,
//                                               3, M_commPtr) );
//        M_ionicModelPtr -> computePotentialRhsSVI (M_globalSolution,
//                                                   M_globalRhs, (*M_feSpacePtr), *M_displacementPtr, vectorialSpace);
//
//    }
//    else
//    {
//        if (M_verbose && M_commPtr -> MyPID() == 0)
//        {
//            std::cout << "\nETA Monodomain Solver: updating rhs with SVI";
//        }
//        M_ionicModelPtr->superIonicModel::computePotentialRhsSVI (M_globalSolution,
//                                                                  M_globalRhs, (*M_feSpacePtr) );
//    }
//    updateRhs();
//}
//
template<typename Mesh, typename IonicModel>
void EMMonodomainSolver<Mesh, IonicModel>::solveOneICIStep()
{
	updateMatrices();
	if (M_displacementPtr  && !M_oneWayCoupling && M_mechanicsModifiesConductivity)
	{
        this->M_linearSolverPtr->setOperator (this->M_globalMatrixPtr);
    }
	super::computeRhsICI();
    this->M_linearSolverPtr->setRightHandSide (this->M_rhsPtrUnique);
    this->M_linearSolverPtr->solve (this->M_potentialPtr);
}
//
//template<typename Mesh, typename IonicModel>
//void EMMonodomainSolver<Mesh, IonicModel>::solveOneICIStep (matrix_Type& mass)
//{
//    computeRhsICI (mass);
//    if (M_displacementPtr)
//    {
//        M_linearSolverPtr->setOperator (M_globalMatrixPtr);
//    }
//    M_linearSolverPtr->setRightHandSide (M_rhsPtrUnique);
//    M_linearSolverPtr->solve (M_potentialPtr);
//}
//
//template<typename Mesh, typename IonicModel>
//void EMMonodomainSolver<Mesh, IonicModel>::solveOneSVIStep()
//{
//    computeRhsSVI();
//    if (M_displacementPtr)
//    {
//        M_linearSolverPtr->setOperator (M_globalMatrixPtr);
//    }
//    M_linearSolverPtr->setRightHandSide (M_rhsPtrUnique);
//    M_linearSolverPtr->solve (M_potentialPtr);
//}
//













} // namespace LifeV

#endif //_MONODOMAINSOLVER_H_
