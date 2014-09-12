/*
 * MixedStructuralOperator.h
 *
 *  Created on: 15/ago/2014
 *      Author: srossi
 */

#ifndef MIXEDSTRUCTURALOPERATOR_H_
#define MIXEDSTRUCTURALOPERATOR_H_

#include <lifev/core/filter/ExporterHDF5.hpp>

#include <lifev/bc_interface/3D/bc/BCInterface3D.hpp>
#include <lifev/core/fem/BCHandler.hpp>
#include <lifev/core/fem/BCManage.hpp>

#include <lifev/eta/fem/ETFESpace.hpp>
#include <lifev/eta/expression/Integrate.hpp>

#include <lifev/core/algorithm/LinearSolver.hpp>
#include <lifev/core/algorithm/Preconditioner.hpp>
#include <lifev/core/algorithm/PreconditionerML.hpp>
#include <lifev/core/algorithm/PreconditionerIfpack.hpp>

#include <lifev/core/array/VectorEpetraStructured.hpp>
#include <lifev/core/array/VectorEpetraStructuredView.hpp>
#include <lifev/core/array/MatrixEpetraStructured.hpp>
#include <lifev/core/array/MatrixEpetraStructuredView.hpp>
#include <lifev/core/array/MatrixEpetraStructuredUtility.hpp>

#include <lifev/core/algorithm/NonLinearRichardson.hpp>

#include <lifev/structure/solver/NeoHookeanMaterialNonLinear.hpp>
#include <lifev/structure/solver/StructuralConstitutiveLaw.hpp>
#include <lifev/structure/solver/StructuralConstitutiveLawData.hpp>


namespace LifeV {

template <typename Mesh>
class MixedStructuralOperator {

	typedef VectorEpetra                                vector_Type;
	typedef boost::shared_ptr<vector_Type>           vectorPtr_Type;
    typedef VectorEpetraStructured                 vectorBlock_Type;
    typedef boost::shared_ptr<vectorBlock_Type> vectorBlockPtr_Type;

	typedef MatrixEpetra<Real>                          matrix_Type;
	typedef boost::shared_ptr<matrix_Type>           matrixPtr_Type;
    typedef MatrixEpetraStructured<Real>           matrixBlock_Type;
    typedef boost::shared_ptr<matrixBlock_Type> matrixBlockPtr_Type;

   typedef Real ( *function ) ( const Real&, const Real&, const Real&, const Real&, const ID& );

	typedef boost::function<Real ( Real const&, Real const&, Real const&, Real const&, ID const& ) > source_Type;

	typedef StructuralConstitutiveLaw<Mesh>               material_Type;
	typedef boost::shared_ptr<material_Type>              materialPtr_Type;

	typedef BCHandler                                     bcHandler_Type;
	typedef boost::shared_ptr<bcHandler_Type>          bcHandlerPtr_Type;

	typedef StructuralConstitutiveLawData                 data_Type;

	typedef Mesh                                          mesh_Type;
	typedef boost::shared_ptr<Mesh>					   meshPtr_Type;


    typedef LinearSolver                            linearSolver_Type;
    typedef boost::shared_ptr<LinearSolver>         linearSolverPtr_Type;
    typedef LifeV::Preconditioner                   basePrec_Type;
    typedef boost::shared_ptr<basePrec_Type>        basePrecPtr_Type;
	typedef LifeV::PreconditionerIfpack             precIfpack_Type;
	typedef boost::shared_ptr<precIfpack_Type>      precIfpackPtr_Type;
	typedef LifeV::PreconditionerML                 precML_Type;
	typedef boost::shared_ptr<precML_Type>          precMLPtr_Type;

	    // Time advance
//	typedef TimeAdvance< vector_Type >                                  timeAdvance_Type;
//	typedef boost::shared_ptr< timeAdvance_Type >                       timeAdvancePtr_Type;

	typedef FESpace< RegionMesh<LinearTetra>, MapEpetra >               FESpace_Type;
	typedef boost::shared_ptr<FESpace_Type>                             FESpacePtr_Type;

	typedef ETFESpace< RegionMesh<LinearTetra>, MapEpetra, 3, 1 >       scalarETFESpace_Type;
	typedef boost::shared_ptr<scalarETFESpace_Type>                     scalarETFESpacePtr_Type;
	typedef ETFESpace< RegionMesh<LinearTetra>, MapEpetra, 3, 3 >       vectorETFESpace_Type;
	typedef boost::shared_ptr<vectorETFESpace_Type>                     vectorETFESpacePtr_Type;

	typedef Epetra_Comm                        comm_Type;
	typedef boost::shared_ptr<Epetra_Comm>  commPtr_Type;

	typedef ExporterHDF5< Mesh >   exporter_Type;
	typedef boost::shared_ptr<exporter_Type> exportePtr_Type;

	typedef TimeData        dataTime_Type;
	typedef boost::shared_ptr<dataTime_Type> dataTimePtr_Type;
public:


	MixedStructuralOperator(meshPtr_Type meshPtr, commPtr_Type comm);

	virtual ~MixedStructuralOperator() {}

	void setup( GetPot& datafile, std::string problemFolder = "./");

	void evalResidual( vectorBlock_Type& residual, vectorBlock_Type& solution, UInt iter);

	void solveJac( vectorBlock_Type& step, vectorBlock_Type& residual, Real linearRelTol);

	void assembleVelocityMassMatrix ();
	void assemblePressureMassMatrix ();
	void assemblePressureBlock12();
	void assemblePressureBlock21();
	void addMaterialResidual();

	void iterate();

	void updateOldSolution();

	commPtr_Type            M_commPtr;
	meshPtr_Type            M_meshPtr;
	FESpacePtr_Type         M_pressureFESpacePtr;
	FESpacePtr_Type         M_velocityFESpacePtr;
	scalarETFESpacePtr_Type M_pressureETFESpacePtr;
	vectorETFESpacePtr_Type M_velocityETFESpacePtr;

	vectorPtr_Type           M_displacementPtr;
	vectorPtr_Type           M_oldDisplacementPtr;
	vectorPtr_Type           M_velocityPtr;
	vectorPtr_Type           M_pressurePtr;

	vectorBlockPtr_Type        M_solutionPtr;
	vectorBlockPtr_Type        M_oldSolutionPtr;
	vectorBlockPtr_Type        M_blockRhsRepPtr;
	vectorBlockPtr_Type        M_blockRhsUniquePtr;

	matrixBlockPtr_Type            M_velocityMassMatrixPtr;
	matrixBlockPtr_Type            M_pressureMassMatrixPtr;
	matrixBlockPtr_Type     M_systemMatrixPtr;

	bcHandlerPtr_Type          M_bchandlerPtr;

	linearSolverPtr_Type   M_linearSolverPtr;
	basePrecPtr_Type       M_preconditionerPtr;

    boost::shared_ptr<data_Type>         M_dataPtr;
    materialPtr_Type                     M_materialPtr;

    exportePtr_Type        M_exporterPtr;

    dataTimePtr_Type       M_dataTimePtr;

    UInt 					M_counter;

};


template <typename Mesh>
MixedStructuralOperator<Mesh>::MixedStructuralOperator( meshPtr_Type meshPtr, commPtr_Type comm ):
M_meshPtr ( meshPtr),
M_commPtr ( comm ),
M_pressureFESpacePtr( new FESpace_Type( M_meshPtr, "P1", 1, M_commPtr) ),
M_velocityFESpacePtr( new FESpace_Type( M_meshPtr, "P2", 3, M_commPtr) ),
M_pressureETFESpacePtr( new scalarETFESpace_Type( M_meshPtr, & (M_pressureFESpacePtr->refFE() ), M_commPtr) ),
M_velocityETFESpacePtr( new vectorETFESpace_Type( M_meshPtr, & (M_velocityFESpacePtr->refFE() ), M_commPtr) ),
M_displacementPtr ( new vector_Type( M_velocityFESpacePtr -> map() ) ),
M_oldDisplacementPtr ( new vector_Type( M_velocityFESpacePtr -> map() ) ),
M_velocityPtr     ( new vector_Type( M_velocityFESpacePtr -> map() ) ),
M_pressurePtr     ( new vector_Type( M_pressureFESpacePtr -> map() ) ),
M_solutionPtr     ( new vectorBlock_Type( M_velocityFESpacePtr -> map() | M_pressureFESpacePtr -> map() ) ),
M_oldSolutionPtr  ( new vectorBlock_Type( M_velocityFESpacePtr -> map() | M_pressureFESpacePtr -> map() ) ),
M_blockRhsUniquePtr( new vectorBlock_Type( M_velocityFESpacePtr -> map() | M_pressureFESpacePtr -> map() ) ),
M_blockRhsRepPtr  ( new vectorBlock_Type( M_velocityFESpacePtr -> map() | M_pressureFESpacePtr -> map(), Repeated ) ),
M_velocityMassMatrixPtr( new matrixBlock_Type( M_velocityFESpacePtr -> map() ) ),
M_pressureMassMatrixPtr( new matrixBlock_Type( M_pressureFESpacePtr -> map() ) ),
M_systemMatrixPtr ( new matrixBlock_Type( M_velocityFESpacePtr -> map() | M_pressureFESpacePtr -> map() ) ),
M_bchandlerPtr    ( new BCHandler() ),
M_linearSolverPtr ( new linearSolver_Type( M_commPtr ) ),
M_preconditionerPtr(),
M_dataPtr( new StructuralConstitutiveLawData( ) ),
M_materialPtr(),
M_exporterPtr ( new exporter_Type() ),
M_dataTimePtr ( new dataTime_Type() ),
M_counter (0)
{
}


template<typename Mesh>
void
MixedStructuralOperator<Mesh>::setup( GetPot& datafile, std::string problemFolder )
{
	//setting up the container of the time data, i.e. initial time, end time, timestep etc. etc.
	M_dataTimePtr->setup(datafile, "solid/time_discretization" );
	//set up the data of the material law, that is the parameters of the constitutive law
    M_dataPtr->setup (datafile);
    //Create a new material of type as specified in the M_dataPtr.
    // Currently working only with the NeoHookean
    M_materialPtr.reset ( material_Type::StructureMaterialFactory::instance().createObject ( M_dataPtr->solidType() ) );
    boost::shared_ptr<Displayer> displayerPtr(new Displayer (M_commPtr) );
    M_materialPtr->setup ( M_velocityFESpacePtr,
    		               M_velocityETFESpacePtr,
    		               M_velocityFESpacePtr->mapPtr(), 0, M_dataPtr, displayerPtr );
	//To use a mixed fomrulation we had to add some additional information in the material,
    // in particular the pointers to the pressure spaces. We set those separately.
    M_materialPtr -> setPressureFESpace(M_pressureFESpacePtr);
	M_materialPtr -> setPressureETFESpace(M_pressureETFESpacePtr);

	//Setting up the preconditioner
	precIfpack_Type* precIfpack = new precIfpack_Type;
	precIfpack->setDataFromGetPot (datafile, "prec");
	M_preconditionerPtr.reset(precIfpack);

	Teuchos::RCP < Teuchos::ParameterList > solverParamList = Teuchos::rcp (
	                                                                  new Teuchos::ParameterList);
	//It looks in the electrophysiology path. This should be changed
	//It is like this because it was copied from a file in the electrophysiology module.
	std::string xmlpath = datafile ("electrophysiology/monodomain_xml_path",
	                                    "./");
	std::string xmlfile = datafile ("electrophysiology/monodomain_xml_file",
	                                    "MonodomainSolverParamList.xml");

	solverParamList = Teuchos::getParametersFromXmlFile (xmlpath + xmlfile);

	//setup the linear solver
	M_linearSolverPtr->setParameters ( *solverParamList );
    M_linearSolverPtr->setPreconditioner ( M_preconditionerPtr );
    M_linearSolverPtr->setOperator ( M_systemMatrixPtr );

    //Initialize all the vectors and matrices to zero
    *M_displacementPtr *= 0.0;
    *M_velocityPtr *= 0.0;
    *M_pressurePtr *= 0.0;
    UInt offset = M_velocityPtr -> size();
    M_solutionPtr->replace(*M_pressurePtr, offset);
    *M_oldDisplacementPtr *= *M_displacementPtr;
    updateOldSolution();
    *M_systemMatrixPtr *= 0.0;
    *M_pressureMassMatrixPtr *= 0.0;

    //setup exporter and export initial solution
    M_exporterPtr->setMeshProcId (M_meshPtr, M_commPtr->MyPID() );
    M_exporterPtr->setPrefix ("solution");
    M_exporterPtr->setPostDir(problemFolder);

    M_exporterPtr->addVariable(ExporterData<mesh_Type>::VectorField, "displacement",
             M_velocityFESpacePtr, M_displacementPtr, UInt(0));
    M_exporterPtr->addVariable(ExporterData<mesh_Type>::VectorField, "velocity",
    		 M_velocityFESpacePtr, M_velocityPtr, UInt(0));
    M_exporterPtr->addVariable(ExporterData<mesh_Type>::ScalarField, "pressure",
    		 M_pressureFESpacePtr, M_pressurePtr, UInt(0));

    M_exporterPtr->postProcess(M_dataTimePtr->time());
}

template< typename Mesh >
void
MixedStructuralOperator<Mesh>::updateOldSolution()
{
	//copy the vecotr velocity and pressure (v, p) in the monolithic old solution vector
    UInt offset = M_velocityPtr -> size();
	M_oldSolutionPtr->replace(*M_velocityPtr, 0);
	M_oldSolutionPtr->replace(*M_pressurePtr, offset);
}

template< typename Mesh >
void
MixedStructuralOperator<Mesh>::assembleVelocityMassMatrix()
{
	*M_velocityMassMatrixPtr *= 0.0;
	M_velocityMassMatrixPtr->openCrsMatrix();
	{
    	using namespace ExpressionAssembly;
    	integrate (
            elements ( M_meshPtr ), // Mesh
            M_velocityFESpacePtr->qr(),
            M_velocityETFESpacePtr,
            M_velocityETFESpacePtr,
            dot ( phi_i , phi_j )
            )
            >> M_velocityMassMatrixPtr;
	}

	M_velocityMassMatrixPtr->globalAssemble();
}

template< typename Mesh >
void
MixedStructuralOperator<Mesh>::assemblePressureMassMatrix()
{
	*M_pressureMassMatrixPtr *= 0.0;
	M_pressureMassMatrixPtr->openCrsMatrix();
	{
    	using namespace ExpressionAssembly;
    	integrate (
            elements ( M_meshPtr ), // Mesh
            M_pressureFESpacePtr->qr(),
            M_pressureETFESpacePtr,
            M_pressureETFESpacePtr,
            phi_i * phi_j
            )
            >> M_pressureMassMatrixPtr;

    	integrate (
            elements ( M_meshPtr ), // Mesh
            M_pressureFESpacePtr->qr(),
            M_pressureETFESpacePtr,
            M_pressureETFESpacePtr,
            phi_i * phi_j
            )
            >> M_systemMatrixPtr->block(1,1);
	}

	M_pressureMassMatrixPtr->globalAssemble();
}


template< typename Mesh >
void
MixedStructuralOperator<Mesh>::assemblePressureBlock12()
{
	using namespace ExpressionAssembly;

	MatrixSmall<3,3> Id;
    Id (0, 0) = 1.0;
    Id (0, 1) = 0.0;
    Id (0, 2) = 0.0;
    Id (1, 0) = 0.0;
    Id (1, 1) = 1.0;
    Id (1, 2) = 0.0;
    Id (2, 0) = 0.0;
    Id (2, 1) = 0.0;
    Id (2, 2) = 1.0;

    auto dt =  value( M_dataTimePtr->timeStep() );
    auto F = grad( M_velocityETFESpacePtr,  *M_displacementPtr, 0 ) + value(Id);
    auto FmT = minusT(F);
    auto J = det(F);
    auto dp = phi_j;

    integrate ( elements ( M_meshPtr ) ,
                quadRuleTetra4pt,
                M_velocityETFESpacePtr,
                M_pressureETFESpacePtr,
               dt *  dot ( dp * J * FmT, grad (phi_i) )
              ) >> M_systemMatrixPtr->block (0, 1);
}


template< typename Mesh >
void
MixedStructuralOperator<Mesh>::assemblePressureBlock21()
{
	using namespace ExpressionAssembly;

	MatrixSmall<3,3> Id;
    Id (0, 0) = 1.0;
    Id (0, 1) = 0.0;
    Id (0, 2) = 0.0;
    Id (1, 0) = 0.0;
    Id (1, 1) = 1.0;
    Id (1, 2) = 0.0;
    Id (2, 0) = 0.0;
    Id (2, 1) = 0.0;
    Id (2, 2) = 1.0;

    auto dt = value( M_dataTimePtr->timeStep() );
    auto dv = grad( M_velocityETFESpacePtr,  *M_velocityPtr, 0 );
    auto div_v = trace(dv);
    auto F = grad( M_velocityETFESpacePtr,  *M_displacementPtr, 0 ) + value(Id);
    auto FmT = minusT(F);
    auto J = det(F);
    auto K = value(100000);

    auto dF = grad(phi_j);

    auto Jd2U = value(2.) * K * (1.0 - log(J) ) / J;
    auto dJd2UdF = value(2.0) * K * ( log(J) - value(2.0) ) / J * dot(FmT, dF);

    auto dpdF = dt * dt * dJd2UdF * div_v + dt * Jd2U * div(phi_j);
//    auto dpdF = dt * div(phi_j);

    integrate ( elements ( M_meshPtr ) ,
                quadRuleTetra4pt,
                M_pressureETFESpacePtr,
                M_velocityETFESpacePtr,
                value(-1.0) * dpdF * phi_i
              ) >> M_systemMatrixPtr->block (1, 0);

}


template< typename Mesh >
void
MixedStructuralOperator<Mesh>::solveJac( vectorBlock_Type& step, vectorBlock_Type& residual, Real linearRelTol)
{

	*M_systemMatrixPtr *= 0.0;
	M_systemMatrixPtr->openCrsMatrix();
    vectorBlockPtr_Type pointerToRes ( new vectorBlock_Type (residual) );
    vectorBlockPtr_Type pointerToStep ( new vectorBlock_Type (step) );

	assembleVelocityMassMatrix ();
	assemblePressureMassMatrix ();
	assemblePressureBlock12();
	assemblePressureBlock21();

    M_materialPtr->updateJacobianMatrix (*M_displacementPtr, *M_pressurePtr, M_dataPtr);
	Real dt2 = M_dataTimePtr->timeStep();
	dt2 *= dt2;
	matrixBlock_Type tmpMatrix( *M_velocityMassMatrixPtr );
	*(M_materialPtr->stiffMatrix() ) *= dt2;
	tmpMatrix +=  *(M_materialPtr->stiffMatrix());

    MatrixEpetraStructuredView<double> A00, A;
    M_systemMatrixPtr->blockView ( 0, 0, A00 );
    tmpMatrix.blockView ( 0, 0, A );

	MatrixEpetraStructuredUtility::copyBlock(A, A00);

	M_systemMatrixPtr->globalAssemble();

	if( !M_bchandlerPtr->bcUpdateDone() )
	{
		M_bchandlerPtr->bcUpdate ( *M_meshPtr, M_velocityFESpacePtr->feBd(), M_velocityFESpacePtr->dof() );
	}

    bcManageMatrix ( *M_systemMatrixPtr,
    		         *M_meshPtr,
    		         M_velocityFESpacePtr->dof(),
    		         *M_bchandlerPtr,
    		         M_velocityFESpacePtr->feBd(),
    		         1.0 );



        ++M_counter;

    M_linearSolverPtr->setOperator (M_systemMatrixPtr);
	M_linearSolverPtr->setTolerance(linearRelTol);
	M_linearSolverPtr->setRightHandSide (pointerToRes);
	M_linearSolverPtr->solve (pointerToStep);

	step = *pointerToStep;
}

template< typename Mesh >
void
MixedStructuralOperator<Mesh>::evalResidual( vectorBlock_Type& residual,
											 vectorBlock_Type& solution,
											 UInt iter)
{
	Real dt = M_dataTimePtr->timeStep();

	M_velocityPtr->subset(solution);
	UInt offset = M_velocityPtr->size();
	M_pressurePtr->subset(solution, offset);
	vectorPtr_Type oldVelocity(new vector_Type( M_velocityPtr->map() ) );
	vectorPtr_Type residualVelocity( new vector_Type( M_velocityPtr->map(), Repeated ) );
	vectorPtr_Type oldPressure( new vector_Type( M_pressurePtr->map() ) ) ;
	vectorPtr_Type residualPressure( new vector_Type( M_pressurePtr->map(), Repeated ) );
	*residualVelocity *= 0.;
	*residualPressure *= 0.;

	*M_displacementPtr = *M_oldDisplacementPtr;
	*M_displacementPtr += ( dt * *M_velocityPtr);

	oldVelocity->subset( *M_oldSolutionPtr );
	oldPressure->subset( *M_oldSolutionPtr, offset );


	M_materialPtr->computeStiffness(*M_displacementPtr, *M_pressurePtr, M_dataPtr);

	*residualVelocity += *(M_materialPtr->stiffVector() );
	*residualVelocity *= dt;

	assembleVelocityMassMatrix();

	vector_Type tmpVector( *M_velocityPtr );
	vector_Type tmpVector2( M_velocityPtr->map() );
	tmpVector -= *oldVelocity;
	tmpVector2 *= 0.0;

	M_velocityMassMatrixPtr->multiply(false, tmpVector, tmpVector2);

	*residualVelocity += tmpVector2;

	{
		using namespace ExpressionAssembly;

		MatrixSmall<3,3> Id;
	    Id (0, 0) = 1.0;
	    Id (0, 1) = 0.0;
	    Id (0, 2) = 0.0;
	    Id (1, 0) = 0.0;
	    Id (1, 1) = 1.0;
	    Id (1, 2) = 0.0;
	    Id (2, 0) = 0.0;
	    Id (2, 1) = 0.0;
	    Id (2, 2) = 1.0;

	    auto dt = value( M_dataTimePtr->timeStep() );
	    auto dv = grad( M_velocityETFESpacePtr,  *M_velocityPtr, 0 );
	    auto div_v = trace(dv);
	    auto F = grad( M_velocityETFESpacePtr,  *M_displacementPtr, 0 ) + value(Id);
	    auto FmT = minusT(F);
	    auto J = det(F);
	    auto K = value(100000);

	    auto Jd2U = value(2.) * K * (1.0 - log(J) ) / J;

	    integrate ( elements ( M_meshPtr ) ,
	                quadRuleTetra4pt,
	                M_pressureETFESpacePtr,
	                value(-1.0) * Jd2U * div_v * phi_i
	              ) >> residualPressure; //->block (1, 0);

		residualPressure->globalAssemble();
	}

	*residualPressure *= dt;

	assemblePressureMassMatrix();

	vector_Type tmpVector3( *M_pressurePtr );
	vector_Type tmpVector4( M_pressurePtr->map() );
	tmpVector3 -= *oldPressure;
	M_pressureMassMatrixPtr->multiply(false, tmpVector3, tmpVector4);
	*residualPressure += tmpVector4;

	vector_Type resV(*residualVelocity, Unique);
	vector_Type resP(*residualPressure, Unique);
	residual.replace(resV, 0);
	residual.replace(resP, offset);



    vectorBlock_Type solRep (solution, Repeated);
    vectorBlock_Type rhs(solution);
    rhs *= 0.0;

    if ( !M_bchandlerPtr->bcUpdateDone() )
    {
     M_bchandlerPtr->bcUpdate ( *M_meshPtr, M_velocityFESpacePtr->feBd(), M_velocityFESpacePtr->dof() );
    }
//    M_bchandlerPtr -> showMe();
    bcManageResidual ( residual,
    		           rhs,
    		           solRep,
    		           *M_meshPtr,
    		           M_velocityFESpacePtr->dof(),
    		           *M_bchandlerPtr,
    		           M_velocityFESpacePtr->feBd(),
    		           M_dataTimePtr->time(),
    		           1.0 );

	residual -= rhs;
}


template< typename Mesh >
void
MixedStructuralOperator<Mesh>::iterate()
{
    Real time = M_dataTimePtr->time();
    std::cout << "\nTime : " << time << "\n";
    std::cout << "\nEnd Time : " << M_dataTimePtr->endTime() << "\n";
    for( ; time < M_dataTimePtr->endTime() ;)
    {

        Real abstol  = 1.e-9;
        Real reltol  = 1.e-9;
        UInt maxiter = 10;
        Real etamax  = 1e-7;
        Int NonLinearLineSearch = 0;
		Int status = 0;

		status = NonLinearRichardson ( *M_solutionPtr, *this, abstol, reltol, maxiter, etamax, NonLinearLineSearch, 0, 2 );

		updateOldSolution();
		*M_oldDisplacementPtr += ( M_dataTimePtr->timeStep() * *M_velocityPtr);
		M_dataTimePtr->updateTime();
		time = M_dataTimePtr->time();
		M_exporterPtr->postProcess(time);

    }
	M_exporterPtr->closeFile();

}



} /* namespace LifeV */

#endif /* MIXEDSTRUCTURALOPERATOR_H_ */
