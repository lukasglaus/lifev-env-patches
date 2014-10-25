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

#ifndef _EMACTIVESTRAINSOLVER_H_
#define _EMACTIVESTRAINSOLVER_H_

#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <lifev/core/filter/ExporterEnsight.hpp>
#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5.hpp>
#endif
#include <lifev/core/filter/ExporterEmpty.hpp>

#include <lifev/eta/fem/ETFESpace.hpp>
#include <lifev/eta/expression/Integrate.hpp>
#include <lifev/structure/solver/StructuralOperator.hpp>
#include <boost/typeof/typeof.hpp>
#include <lifev/em/solver/EMETAFunctors.hpp>
#include <lifev/em/util/EMUtility.hpp>

namespace LifeV
{

//! EMSolver - Class featuring the solution of the electromechanical problem with monodomain equation

template<typename Mesh>
class EMActiveStrainSolver
{

    //!Monodomain Solver
    /*!
     The monodomain equation reads
     \f \Chi

     */

public:

    //! @name Type definitions
    //@{

    typedef Mesh mesh_Type;
    typedef boost::shared_ptr<mesh_Type> meshPtr_Type;

    typedef VectorEpetra vector_Type;
    typedef boost::shared_ptr<VectorEpetra> vectorPtr_Type;

    typedef std::vector<vectorPtr_Type> vectorOfPtr_Type;

    typedef MatrixEpetra<Real> matrix_Type;
    typedef boost::shared_ptr<matrix_Type> matrixPtr_Type;

    typedef Epetra_Comm comm_Type;
    typedef boost::shared_ptr<comm_Type> commPtr_Type;

    typedef ETFESpace<mesh_Type, MapEpetra, 3, 1> ETFESpace_Type;
    typedef boost::shared_ptr<ETFESpace<mesh_Type, MapEpetra, 3, 1> > ETFESpacePtr_Type;

    typedef ETFESpace<mesh_Type, MapEpetra, 3, 3> ETFESpaceVectorial_Type;
    typedef boost::shared_ptr< ETFESpaceVectorial_Type > ETFESpaceVectorialPtr_Type;

    typedef FESpace<mesh_Type, MapEpetra> feSpace_Type;
    typedef boost::shared_ptr<feSpace_Type> feSpacePtr_Type;

    typedef LinearSolver linearSolver_Type;
    typedef boost::shared_ptr<LinearSolver> linearSolverPtr_Type;

    //    typedef ExporterHDF5< mesh_Type >          exporter_Type;
    //    typedef boost::shared_ptr<exporter_Type>                       exporterPtr_Type;
    //  typedef Exporter<mesh_Type> exporter_Type;    //                IOFile_Type;
    //  typedef boost::shared_ptr<exporter_Type> exporterPtr_Type; //                IOFilePtr_Type;

    typedef LifeV::Preconditioner basePrec_Type;
    typedef boost::shared_ptr<basePrec_Type> basePrecPtr_Type;
    typedef LifeV::PreconditionerIfpack prec_Type;
    typedef boost::shared_ptr<prec_Type> precPtr_Type;

    typedef Teuchos::ParameterList list_Type;

    typedef boost::function <
    Real (const Real& t, const Real& x, const Real& y, const Real& z,
          const ID& i) > function_Type;

    typedef MatrixSmall<3, 3>                          matrixSmall_Type;

    //  typedef ElectroETAMonodomainSolver<mesh_Type, ionicModel_Type>                  monodomainSolver_Type;
    //  typedef boost::shared_ptr<monodomainSolver_Type>                    monodomainSolverPtr_Type;

    typedef FESpace< RegionMesh<LinearTetra>, MapEpetra >               FESpace_Type;
    typedef boost::shared_ptr<FESpace_Type>                        FESpacePtr_Type;

    typedef ETFESpace< RegionMesh<LinearTetra>, MapEpetra, 3, 1 >       scalarETFESpace_Type;
    typedef boost::shared_ptr<scalarETFESpace_Type>                      scalarETFESpacePtr_Type;
    typedef ETFESpace< RegionMesh<LinearTetra>, MapEpetra, 3, 3 >       solidETFESpace_Type;
    typedef boost::shared_ptr<solidETFESpace_Type>                      solidETFESpacePtr_Type;

    typedef StructuralConstitutiveLawData           structureData_Type;
    typedef boost::shared_ptr<structureData_Type>           structureDataPtr_Type;

    typedef EMStructuralOperator< RegionMesh<LinearTetra> > structuralOperator_Type;
    typedef boost::shared_ptr< structuralOperator_Type > structuralOperatorPtr_Type;

    //
    //    typedef BCHandler                                          bc_Type;
    //    typedef boost::shared_ptr< bc_Type >                       bcPtr_Type;
    //    typedef StructuralOperator< RegionMesh<LinearTetra> >     physicalSolver_Type;
    //    typedef BCInterface3D< bc_Type, physicalSolver_Type >      bcInterface_Type;
    //    typedef boost::shared_ptr< bcInterface_Type >              bcInterfacePtr_Type;

    typedef ExporterHDF5< RegionMesh <LinearTetra> >        exporter_Type;
    typedef boost::shared_ptr< exporter_Type >             exporterPtr_Type;
    //@}


    enum activation_Type { Orthotropic, TransverselyIsotropic };
    activation_Type                 M_activationType;



    ////getters
    inline vectorPtr_Type   gammafPtr()
    {
        return M_gammafPtr;
    }
    inline vector_Type&     gammaf()
    {
        return *M_gammafPtr;
    }
    inline vectorPtr_Type   gammasPtr()
    {
        return M_gammasPtr;
    }
    inline  vectorPtr_Type  gammanPtr()
    {
        return M_gammanPtr;
    }

    inline  meshPtr_Type    meshPtr()
    {
        return M_meshPtr;
    }

    inline  scalarETFESpacePtr_Type ETFESpacePtr()
    {
        return M_ETFESpacePtr;
    }
    inline  solidETFESpacePtr_Type  solidETFESpacePtr()
    {
        return M_solidETFESpacePtr;
    }
    inline  FESpacePtr_Type     FESpacePtr()
    {
        return M_FESpacePtr;
    }
    inline  FESpacePtr_Type     FESpacePtrSolid()
    {
        return M_FESpacePtrSolid;
    }

    inline  activation_Type         activationType()
    {
        return M_activationType;
    }
    inline  Real    orthotropicActivationAnisotropyRatio()
    {
        return M_orthotropicActivationAnisotropyRatio;
    }
    inline  Real    activeViscosity()
    {
        return M_activeViscosity;
    }
    inline  Real    CaDiastolic()
    {
        return M_CaDiastolic;
    }
    inline  Real    activeCoefficient()
    {
        return M_activeCoefficient;
    }

    inline      matrixPtr_Type  massMatrixPtr()
    {
        return M_massMatrixPtr;
    }
    inline      vectorPtr_Type  rhsPtr()
    {
        return M_rhsPtr;
    }
    inline      vectorPtr_Type  fiberPtr()
    {
        return M_fiberPtr;
    }

    inline      basePrecPtr_Type    precPtr()
    {
        return M_precPtr;
    }
    inline      linearSolver_Type&  linearSolver()
    {
        return M_linearSolver;
    }
    inline      exporterPtr_Type    exporterPtr()
    {
        return M_exporterPtr;
    }


    ////setters
    inline void     setGammafPtr (vectorPtr_Type p)
    {
        M_gammafPtr = p;
    }
    inline void     setGammafPtr (vector_Type& p)
    {
        *M_gammafPtr = p;
    }

    inline  void    setMeshPtr ( meshPtr_Type p)
    {
        M_meshPtr = p;
    }
    inline  void    setMeshPtr ( mesh_Type& p)
    {
        *M_meshPtr = p;
    }

    inline  void setETFESpacePtr ( ETFESpacePtr_Type p)
    {
        M_ETFESpacePtr = p;
    }
    inline  void setSolidETFESpacePtr ( solidETFESpacePtr_Type p)
    {
        M_solidETFESpacePtr = p;
    }
    inline  void setFESpacePtr ( FESpacePtr_Type p)
    {
        M_FESpacePtr = p;
    }

    inline  void setActivationType ( std::string p)
    {
        std::map< std::string, activation_Type > activationTypeMap;
        activationTypeMap["TransverselyIsotropic"] = TransverselyIsotropic;
        activationTypeMap["Orthotropic"]           = Orthotropic;
        M_activationType   = activationTypeMap[p];
    }

    inline  void    setOrthotropicActivationAnisotropyRatio (Real p)
    {
        M_orthotropicActivationAnisotropyRatio = p;
    }
    inline  void    setActiveViscosity (Real p)
    {
        M_activeViscosity = p;
    }
    inline  void    setCaDiastolic (Real p)
    {
        M_CaDiastolic = p;
    }
    inline  void    setActiveCoefficient (Real p)
    {
        M_activeCoefficient = p;
    }

    inline  void    setMassMatrixPtr (matrixPtr_Type p )
    {
        M_massMatrixPtr = p;
        M_linearSolver.setOperator (M_massMatrixPtr);
    }
    inline  void    setRhsPtr (vectorPtr_Type p)
    {
        M_rhsPtr = p;
    }
    inline  void    setFiberPtr (vectorPtr_Type p )
    {
        M_fiberPtr = p;
    }
    inline  void    setFiberPtr (vector_Type& p )
    {
        *M_fiberPtr = p;
    }



    //! @name Constructors & Destructor
    //@{
    EMActiveStrainSolver();

    EMActiveStrainSolver ( Teuchos::ParameterList parameterList, GetPot& dataFile, commPtr_Type comm ); //!Empty Constructor

    EMActiveStrainSolver ( Teuchos::ParameterList parameterList, GetPot& dataFile, meshPtr_Type meshPtr, commPtr_Type comm );   //!Empty Constructor

    virtual ~EMActiveStrainSolver() {};

    void setupMassMatrix();

    void setupLumpedMassMatrix();

    void setupLinearSolver (GetPot& dataFile, commPtr_Type comm);

    void setupExporter (exporter_Type& exporter, commPtr_Type comm, std::string outputFolder = "./");

    void setupExporter (commPtr_Type comm, std::string outputFolder = "./");

    void computeGammasAndGamman();

    void computeGammasAndGamman (const vector_Type& gammaf, vector_Type& gammas, vector_Type& gamman);

    inline void computeGammasAndGamman (const vectorPtr_Type gammaf, vectorPtr_Type gammas, vectorPtr_Type gamman)
    {
        computeGammasAndGamman (*gammaf, *gammas, *gamman);
    }

    void setupRhs ( const vector_Type& disp, solidETFESpacePtr_Type solidETFESpacePtr,
                    const vector_Type& calcium, scalarETFESpacePtr_Type ETFESpacePtr, Real dt);

    void solve();

    void readMeshFromList (Teuchos::ParameterList parameterList);

    void exportSolution (Real t = 0.0);

    void showMe ( commPtr_Type comm );
    /*!
     */

    //  ionicModelPtr_Type          M_ionicPtr;
    vectorPtr_Type              M_gammafPtr;
    vectorPtr_Type              M_gammasPtr;
    vectorPtr_Type              M_gammanPtr;

    meshPtr_Type                M_meshPtr;

    scalarETFESpacePtr_Type     M_ETFESpacePtr;
    solidETFESpacePtr_Type      M_solidETFESpacePtr;
    FESpacePtr_Type             M_FESpacePtr;
    FESpacePtr_Type             M_FESpacePtrSolid;


    Real                        M_orthotropicActivationAnisotropyRatio;
    Real                        M_activeViscosity;
    Real                        M_CaDiastolic;
    Real                        M_activeCoefficient;

    matrixPtr_Type              M_massMatrixPtr;
    vectorPtr_Type              M_rhsRepeatedPtr;
    vectorPtr_Type              M_rhsPtr;
    vectorPtr_Type              M_fiberPtr;

    basePrecPtr_Type            M_precPtr;
    linearSolver_Type           M_linearSolver;
    exporterPtr_Type            M_exporterPtr;

    commPtr_Type                M_commPtr;

private:
    void init ( Teuchos::ParameterList parameterList, commPtr_Type comm );
    void init ( Teuchos::ParameterList parameterList, meshPtr_Type mesh, commPtr_Type comm );

    //     solidFESpacePtr_Type dFESpace ( new solidFESpace_Type (localSolidMesh, dOrder, 3, comm) );
    //      solidFESpacePtr_Type aFESpace ( new solidFESpace_Type (monodomain -> localMeshPtr(), dOrder, 1, comm) );
    //      solidETFESpacePtr_Type dETFESpace ( new solidETFESpace_Type (localSolidMesh, & (dFESpace->refFE() ), & (dFESpace->fe().geoMap() ), comm) );
    //      scalarETFESpacePtr_Type aETFESpace ( new scalarETFESpace_Type (monodomain -> localMeshPtr(), & (aFESpace->refFE() ), & (aFESpace->fe().geoMap() ), comm) );
    //      solidFESpacePtr_Type solidaFESpace ( new solidFESpace_Type (localSolidMesh, "P1", 1, comm) );



};

// ===================================================
//! Constructors
// ===================================================
template<typename Mesh>
EMActiveStrainSolver<Mesh>::EMActiveStrainSolver() :
    M_gammafPtr(),
    M_gammasPtr(),
    M_gammanPtr(),
    M_meshPtr(),
    M_ETFESpacePtr(),
    M_solidETFESpacePtr(),
    M_FESpacePtr(),
    M_FESpacePtrSolid(),
    M_activationType (Orthotropic),
    M_orthotropicActivationAnisotropyRatio (4.0),
    M_CaDiastolic (0.02155),
    M_activeViscosity (4000.0),
    M_activeCoefficient (-4.0),
    M_massMatrixPtr(),
    M_rhsRepeatedPtr(),
    M_rhsPtr(),
    M_fiberPtr(),
    M_precPtr(),
    M_linearSolver(),
    M_exporterPtr(),
    M_commPtr()
{}

template<typename Mesh>
EMActiveStrainSolver<Mesh>::EMActiveStrainSolver (  Teuchos::ParameterList parameterList, GetPot& dataFile, commPtr_Type comm )
{
    if (comm->MyPID() == 0)
    {
        std::cout << "\n==========================================";
        std::cout << "\n\t Active Strain Solver Constructor";
        std::cout << "\n==========================================";
        std::cout << "\n\t Initializing ... ";
    }
    init (parameterList, comm);

    if (comm->MyPID() == 0)
    {
        std::cout << "\t Done. \n";
        std::cout << "\t Setting up matrix, preconditioner and linear solver ... ";

    }
    if (parameterList.get ("activationMassLumped", false) )
    {
        setupLumpedMassMatrix();
    }
    else
    {
        setupMassMatrix();
    }
    setupLinearSolver (dataFile, comm);
    showMe ( comm );
    if (comm->MyPID() == 0)
    {
        std::cout << "\t Done.";

    }

}

template<typename Mesh>
EMActiveStrainSolver<Mesh>::EMActiveStrainSolver (  Teuchos::ParameterList parameterList, GetPot& dataFile, meshPtr_Type meshPtr, commPtr_Type comm )
{
    if (comm->MyPID() == 0)
    {
        std::cout << "\n==========================================";
        std::cout << "\n\t Active Strain Solver Constructor";
        std::cout << "\n==========================================";
        std::cout << "\n\t Initializing ... ";
    }
    M_meshPtr = meshPtr;
    init (parameterList, meshPtr, comm);

    if (comm->MyPID() == 0)
    {
        std::cout << "\t Done. \n";
        std::cout << "\t Setting up matrix, preconditioner and linear solver ... ";

    }
    setupMassMatrix();
    setupLinearSolver (dataFile, comm);
    showMe ( comm );

    if (comm->MyPID() == 0)
    {
        std::cout << "\t Done.";

    }

}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
template<typename Mesh>
void EMActiveStrainSolver<Mesh>::setupExporter (exporter_Type& exporter, commPtr_Type comm, std::string outputFolder)
{
    exporter.setMeshProcId ( M_meshPtr, comm->MyPID() );
    exporter.setPrefix ( "ActivationOutput" );
    exporter.setPostDir ( outputFolder );
    exporter.addVariable (ExporterData<mesh_Type>::ScalarField, "gammaf",
                          M_FESpacePtr, M_gammafPtr, UInt (0) );
}
template<typename Mesh>
void EMActiveStrainSolver<Mesh>::setupExporter (commPtr_Type comm, std::string outputFolder)
{
    M_exporterPtr.reset (new exporter_Type() );
    setupExporter ( *M_exporterPtr, comm, outputFolder);
}


template<typename Mesh>
void EMActiveStrainSolver<Mesh>::solve()
{
    if (M_commPtr->MyPID() == 0)
    {
        std::cout << "\n Active Strain Solver: Solving";
    }

    M_linearSolver.solve (M_gammafPtr);
}

template<typename Mesh>
void EMActiveStrainSolver<Mesh>::exportSolution (Real t)
{
    M_exporterPtr -> postProcess (t);
}

template<typename Mesh>
void EMActiveStrainSolver<Mesh>::setupLumpedMassMatrix()
{
    if (M_commPtr -> MyPID() == 0)
    {
        std::cout << "\n Active Strain Solver: Assembling Lumped Mass Matrix";
    }

    M_massMatrixPtr.reset (new matrix_Type ( M_FESpacePtr -> map() ) ) ;
    {
        using namespace ExpressionAssembly;

        integrate ( elements ( M_meshPtr ),
                    quadRuleTetra4ptNodal,
                    M_ETFESpacePtr,
                    M_ETFESpacePtr,
                    phi_i * phi_j ) >> M_massMatrixPtr;

    }
    std::cout << "\nMatrix Assembled!!!\n";
    M_massMatrixPtr -> globalAssemble();
}


template<typename Mesh>
void EMActiveStrainSolver<Mesh>::setupMassMatrix()
{
    if (M_commPtr -> MyPID() == 0)
    {
        std::cout << "\n Active Strain Solver: Assembling Mass Matrix";
    }

    M_massMatrixPtr.reset (new matrix_Type ( M_FESpacePtr -> map() ) ) ;
    {
        using namespace ExpressionAssembly;

        integrate ( elements ( M_meshPtr ),
                    M_FESpacePtr -> qr(),
                    M_ETFESpacePtr,
                    M_ETFESpacePtr,
                    phi_i * phi_j ) >> M_massMatrixPtr;

    }
    std::cout << "\nMatrix Assembled!!!\n";
    M_massMatrixPtr -> globalAssemble();
}

template<typename Mesh>
void EMActiveStrainSolver<Mesh>::setupLinearSolver (GetPot& dataFile, commPtr_Type comm)
{
    if (comm -> MyPID() == 0)
    {
        std::cout << "\nActive Strain Solver: Setting Activation Linear Solver";
    }
    prec_Type* precRawPtr;
    precRawPtr = new prec_Type;
    precRawPtr->setDataFromGetPot (dataFile, "prec");
    M_precPtr.reset (precRawPtr);

    Teuchos::RCP < Teuchos::ParameterList > solverParamList = Teuchos::rcp (
                                                                  new Teuchos::ParameterList);

    std::string xmlpath = dataFile ("activation/activation_xml_path",
                                    "./");
    std::string xmlfile = dataFile ("activation/activation_xml_file",
                                    "ParamList.xml");

    solverParamList = Teuchos::getParametersFromXmlFile (xmlpath + xmlfile);
    //linearSolver_Type linearSolver;
    M_linearSolver.setCommunicator ( comm );
    M_linearSolver.setParameters ( *solverParamList );
    M_linearSolver.setPreconditioner ( M_precPtr );
    M_linearSolver.setOperator ( M_massMatrixPtr );
}

template<typename Mesh>
void EMActiveStrainSolver<Mesh>::computeGammasAndGamman()
{
    computeGammasAndGamman (*M_gammafPtr, *M_gammasPtr, *M_gammanPtr);
}

template<typename Mesh>
void EMActiveStrainSolver<Mesh>::computeGammasAndGamman (const vector_Type& gammaf, vector_Type& gammas, vector_Type& gamman)
{
    switch (M_activationType)
    {
        case Orthotropic:
            gamman = M_orthotropicActivationAnisotropyRatio * gammaf;
            gammas = 1.0;
            gammas /= (1.0 + gammaf);
            gammas /= (1.0 + gamman);
            gammas -= 1.0;
            break;
        case TransverselyIsotropic:
            gammas = 1.0;
            gammas /= (1.0 + gammaf);
            //Computing the square root
            EMUtility::EpetraSqrt (gammas);
            gammas -= 1.0;
            gamman = gammas;
            break;
    }
}


template<typename Mesh>
void EMActiveStrainSolver<Mesh>::setupRhs ( const vector_Type& disp, solidETFESpacePtr_Type solidETFESpacePtr,
                                            const vector_Type& calcium, scalarETFESpacePtr_Type ETFESpacePtr, Real dt)
{
    Real maxCa = calcium.maxValue();
    Real minCa = calcium.minValue();
    if (M_commPtr->MyPID() == 0)
    {
        std::cout << "\n Active Strain Solver: Assembling rhs";
        if (!M_fiberPtr)
        {
            std::cout << "\n Active Strain Solver: Fibers not set cannot contract";
        }
        std::cout << "\n Activae Strain Solver: Max Calcium = " << maxCa << ", min Ca2+ : " << minCa;

    }
    *M_rhsRepeatedPtr *= 0.0;

    {
        using namespace ExpressionAssembly;

        MatrixSmall<3, 3> Id;
        Id (0, 0) = 1.;
        Id (0, 1) = 0., Id (0, 2) = 0.;
        Id (1, 0) = 0.;
        Id (1, 1) = 1., Id (1, 2) = 0.;
        Id (2, 0) = 0.;
        Id (2, 1) = 0., Id (2, 2) = 1.;

        boost::shared_ptr<FLRelationship> fl (new FLRelationship);
        boost::shared_ptr<HeavisideFct> H (new HeavisideFct);


        BOOST_AUTO_TPL (I,      value (Id) );
        BOOST_AUTO_TPL (Grad_u, grad ( solidETFESpacePtr, disp, 0) );
        BOOST_AUTO_TPL (F,      ( Grad_u + I ) );
        BOOST_AUTO_TPL (J,       det (F) );
        BOOST_AUTO_TPL (Jm23,    pow (J, -2. / 3) );
        // Fibres
        BOOST_AUTO_TPL (f0,     value ( M_solidETFESpacePtr, *M_fiberPtr ) );
        BOOST_AUTO_TPL (f,      F * f0 );
        BOOST_AUTO_TPL (I4f,    dot (f, f) );
        BOOST_AUTO_TPL (I4fiso,  Jm23 * I4f);
        // Generalised invariants
        BOOST_AUTO_TPL (gf,  value (M_ETFESpacePtr, *M_gammafPtr) );

        BOOST_AUTO_TPL (dW, value (2.0) * I4fiso * ( value (3.0) * gf + value (-6.0) * gf * gf + value (10.0) * gf * gf * gf + value (-15.0) * gf * gf * gf * gf  + value (21.0) * gf * gf * gf * gf * gf) );
        BOOST_AUTO_TPL (Ca,    value ( ETFESpacePtr, calcium) );
        BOOST_AUTO_TPL (Ca2, Ca * Ca );
        BOOST_AUTO_TPL (dCa, ( Ca - value (M_CaDiastolic) ) );
        BOOST_AUTO_TPL (Pa, value (M_activeCoefficient) * eval (H, dCa) * eval (H, dCa) * eval (fl, I4fiso) );
        BOOST_AUTO_TPL (beta, value (1.0 / M_activeViscosity ) /*/ coeff  /* Jm23 * pow( eval(EXP, ( I1iso + value(-3.0) ) ), -B )*/ );
        BOOST_AUTO_TPL (gamma_dot, beta / ( Ca2 ) * ( Pa - dW )  );

        integrate ( elements ( M_meshPtr ),
                    M_FESpacePtr -> qr() ,
                    M_ETFESpacePtr,
                    gamma_dot  * phi_i
                  ) >> M_rhsRepeatedPtr;

    }

    *M_rhsPtr *= 0;
    *M_rhsPtr = ( * (M_massMatrixPtr) * ( *M_gammafPtr ) );
    //M_rhsRepeatedPtr -> globalAssemble();
    *M_rhsPtr += ( dt  * *M_rhsRepeatedPtr );

    M_rhsPtr -> globalAssemble();
    //  M_rhsRepeatedPtr.reset( new vector_Type( M_) )
    M_linearSolver.setRightHandSide (M_rhsPtr);

}
///////////////////////////////////////
//////////////////////////////////////////

template<typename Mesh>
void EMActiveStrainSolver<Mesh>::init ( Teuchos::ParameterList parameterList, commPtr_Type comm)
{
    if (comm->MyPID() == 0)
    {
        std::cout << "\nActivation Solver: Importing mesh";
    }
    readMeshFromList (parameterList);
    init (parameterList, M_meshPtr, comm);
}

template<typename Mesh>
void EMActiveStrainSolver<Mesh>::init ( Teuchos::ParameterList parameterList, meshPtr_Type meshPtr, commPtr_Type comm)
{
    M_commPtr = comm;
    std::string activationModelType = parameterList.get ( "activation_type", "Orthotropic" );
    setActivationType ( activationModelType );
    M_orthotropicActivationAnisotropyRatio = parameterList.get ("active_orthtropic_ratio", 4.0);
    M_activeViscosity  = parameterList.get ("active_viscosity", 4000.0);//(5000.0),
    M_CaDiastolic  = parameterList.get ("Ca_diastolic", 0.02155);//(0.02155),
    M_activeCoefficient  = parameterList.get ("active_coefficient", -4.0); //(-4.0),
    M_ETFESpacePtr.reset ( new scalarETFESpace_Type (meshPtr,  &feTetraP1, comm) );
    M_solidETFESpacePtr.reset ( new solidETFESpace_Type (meshPtr, &feTetraP1, comm) );
    M_FESpacePtr.reset ( new FESpace_Type (meshPtr, "P1", 1, comm) );
    M_FESpacePtrSolid.reset ( new FESpace_Type (meshPtr, "P1", 3, comm) );

    M_gammafPtr.reset (new vector_Type ( M_FESpacePtr -> map() ) );
    M_gammasPtr.reset (new vector_Type ( M_FESpacePtr -> map() ) );
    M_gammanPtr.reset (new vector_Type ( M_FESpacePtr -> map() ) );
    M_rhsRepeatedPtr.reset (new vector_Type ( M_FESpacePtr -> map(), Repeated ) );
    M_rhsPtr.reset (new vector_Type ( M_FESpacePtr -> map() ) );

    M_massMatrixPtr.reset (new matrix_Type ( M_FESpacePtr -> map() ) ) ;
}

template<typename Mesh>
void EMActiveStrainSolver<Mesh>::readMeshFromList (Teuchos::ParameterList parameterList)
{
    std::string meshName = parameterList.get ("mesh_name", "lid16.mesh");
    std::string meshPath = parameterList.get ("mesh_path", "./");
    MeshUtility::loadMesh (M_meshPtr, meshName,  meshPath );
}

template<typename Mesh>
void EMActiveStrainSolver<Mesh>::showMe ( commPtr_Type comm )
{
    if (comm->MyPID() == 0)
    {
        std::cout << "\n==========================================";
        std::cout << "\n\t ACTIVE STRAIN SOLVER: infos";
        std::cout << "\n==========================================";

        std::cout << "\nOrthotropic activation anisotropy ratio: " <<   M_orthotropicActivationAnisotropyRatio;
        std::cout << "\nActive viscosity: " <<  M_activeViscosity;
        std::cout << "\nDiastolic calcium: " << M_CaDiastolic;
        std::cout << "\nActive coefficient: " <<    M_activeCoefficient;
        std::cout << "\nMesh size: " << M_meshPtr -> numPoints();
        std::cout << "\nRHS: " <<   M_rhsPtr -> size();
        std::cout << "\n==========================================";

    }
}

} // namespace LifeV

#endif //_MONODOMAINSOLVER_H_
