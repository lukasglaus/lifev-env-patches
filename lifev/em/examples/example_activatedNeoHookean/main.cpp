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
/**
   \file main.cpp

   This test is the case of traction of a cube. It does not use the symmetry BCs
   This test uses the FESpace which is the standard in LifeV and the ETFESpace
   The FESpace is used for the BCs of Neumann type since in the ET branch there
   is not the integration over the boundary faces.

   \author Paolo Tricerri <paolo.tricerri@epfl.ch>
   \date 1861-03-17
 */

#ifdef TWODIM
#error test_structure cannot be compiled in 2D
#endif

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

#include <lifev/core/LifeV.hpp>
#include <lifev/core/algorithm/PreconditionerIfpack.hpp>
#include <lifev/core/algorithm/PreconditionerML.hpp>

#include <lifev/core/array/MapEpetra.hpp>

#include <lifev/core/fem/TimeAdvance.hpp>
#include <lifev/core/fem/TimeAdvanceNewmark.hpp>
#include <lifev/core/fem/TimeAdvanceBDF.hpp>

#include <lifev/core/mesh/MeshData.hpp>
#include <lifev/core/mesh/MeshPartitioner.hpp>

#include <lifev/structure/solver/StructuralConstitutiveLawData.hpp>

#include <lifev/structure/solver/StructuralConstitutiveLaw.hpp>
#include <lifev/structure/solver/StructuralOperator.hpp>
#include <lifev/structure/solver/NeoHookeanMaterialNonLinear.hpp>
#include <lifev/structure/solver/NeoHookeanActivatedMaterial.hpp>

#include <lifev/core/filter/ExporterEnsight.hpp>
#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5.hpp>
#endif
#include <lifev/core/filter/ExporterEmpty.hpp>

//Includes for the Expression Template
#include <lifev/eta/fem/ETFESpace.hpp>

#include <iostream>


using namespace LifeV;

int returnValue = EXIT_FAILURE;
enum TimeScheme { BDF_ORDER_ONE = 1, BDF_ORDER_TWO, BDF_ORDER_THREE };

namespace
{
static bool regIF = (PRECFactory::instance().registerProduct ( "Ifpack", &createIfpack ) );
static bool regML = (PRECFactory::instance().registerProduct ( "ML", &createML ) );
}


std::set<UInt> parseList ( const std::string& list )
{
    std::string stringList = list;
    std::set<UInt> setList;
    if ( list == "" )
    {
        return setList;
    }
    size_t commaPos = 0;
    while ( commaPos != std::string::npos )
    {
        commaPos = stringList.find ( "," );
        setList.insert ( atoi ( stringList.substr ( 0, commaPos ).c_str() ) );
        stringList = stringList.substr ( commaPos + 1 );
    }
    setList.insert ( atoi ( stringList.c_str() ) );
    return setList;
}


class Structure
{
public:

    /** @name Constructors, destructor
     */
    //@{
    Structure ( int                                   argc,
                char**                                argv,
                boost::shared_ptr<Epetra_Comm>        structComm );

    ~Structure()
    {}
    //@}

    //@{
    void run();

    //@}

protected:

private:




private:
    struct Private;
    boost::shared_ptr<Private> parameters;

};



struct Structure::Private
{
    Private() :
        rho (1), poisson (1), young (1), bulk (1), alpha (1), gamma (1)
    {}
    typedef boost::function<Real ( Real const&, Real const&, Real const&, Real const&, ID const& ) > fct_type;
    double rho, poisson, young, bulk, alpha, gamma, gammaf;

    std::string data_file_name;

    boost::shared_ptr<Epetra_Comm>     comm;

    static Real bcZero (const Real& /*t*/, const Real&  /*X*/, const Real& /*Y*/, const Real& /*Z*/, const ID& /*i*/)
    {
        return  0.;
    }

    static Real d0 (const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID& i)
    {
        Teuchos::ParameterList list = * ( Teuchos::getParametersFromXmlFile ( "xmlParameters.xml" ) );
        Real M_gammaf = list.get ( "gammaf", 0.0);
        switch (i)
        {
            case 0:
                return  M_gammaf * x;
                break;
            case 1:
                return ( std::sqrt ( 1.0 / ( 1.0 + M_gammaf ) ) - 1.0 ) * y;
                break;
            case 2:
                return ( std::sqrt ( 1.0 / ( 1.0 + M_gammaf ) ) - 1.0 ) * z;
                break;
            default:
                ERROR_MSG ("This entry is not allowed: ud_functions.hpp");
                return 0.;
                break;
        }

    }

    static Real boundaryLoad (const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID& i)
    {
        Teuchos::ParameterList list = * ( Teuchos::getParametersFromXmlFile ( "xmlParameters.xml" ) );
        Real traction = list.get ( "traction", 0.0);
        switch (i)
        {
            case 0:
                return  traction;
                break;
            case 1:
                return 0.0;
                break;
            case 2:
                return 0.0;
                break;
            default:
                ERROR_MSG ("This entry is not allowed: ud_functions.hpp");
                return 0.;
                break;
        }

    }

    static Real gf (const Real& /*t*/, const Real& x, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
    {
        Teuchos::ParameterList list = * ( Teuchos::getParametersFromXmlFile ( "xmlParameters.xml" ) );
        Real M_gammaf = list.get ( "gammaf", 0.0 );
        return  ( M_gammaf );
    }

};



Structure::Structure ( int                                   argc,
                       char**                                argv,
                       boost::shared_ptr<Epetra_Comm>        structComm) :
    parameters ( new Private() )
{
    GetPot command_line (argc, argv);
    string data_file_name = command_line.follow ("data", 2, "-f", "--file");
    GetPot dataFile ( data_file_name );

    parameters->data_file_name = data_file_name;

    parameters->rho     = dataFile ( "solid/physics/density", 1. );
    parameters->young   = dataFile ( "solid/physics/young",   1. );
    parameters->poisson = dataFile ( "solid/physics/poisson", 1. );
    parameters->bulk    = dataFile ( "solid/physics/bulk",    1. );
    parameters->alpha   = dataFile ( "solid/physics/alpha",   1. );
    parameters->gamma   = dataFile ( "solid/physics/gamma",   1. );
    //  M_gammaf  = dataFile ( "solid/physics/gammaf",  0. );

    std::cout << "density = " << parameters->rho     << std::endl
              << "young   = " << parameters->young   << std::endl
              << "poisson = " << parameters->poisson << std::endl
              << "bulk    = " << parameters->bulk    << std::endl
              << "alpha   = " << parameters->alpha   << std::endl
              << "gamma   = " << parameters->gamma   << std::endl;

    parameters->comm = structComm;
    int ntasks = parameters->comm->NumProc();

    if (!parameters->comm->MyPID() )
    {
        std::cout << "My PID = " << parameters->comm->MyPID() << " out of " << ntasks << " running." << std::endl;
    }
}




void
Structure::run()
{
    typedef EMStructuralOperator< RegionMesh<LinearTetra> >::vector_Type  vector_Type;
    typedef boost::shared_ptr<vector_Type>                              vectorPtr_Type;
    typedef boost::shared_ptr< TimeAdvance< vector_Type > >             timeAdvance_Type;
    typedef FESpace< RegionMesh<LinearTetra>, MapEpetra >               solidFESpace_Type;
    typedef boost::shared_ptr<solidFESpace_Type>                        solidFESpacePtr_Type;

    typedef ETFESpace< RegionMesh<LinearTetra>, MapEpetra, 3, 3 >       solidETFESpace_Type;
    typedef boost::shared_ptr<solidETFESpace_Type>                      solidETFESpacePtr_Type;


    bool verbose = (parameters->comm->MyPID() == 0);

    //! Number of boundary conditions for the velocity and mesh motion
    boost::shared_ptr<BCHandler> BCh ( new BCHandler() );

    //! dataElasticStructure
    GetPot dataFile ( parameters->data_file_name.c_str() );

    boost::shared_ptr<StructuralConstitutiveLawData> dataStructure (new StructuralConstitutiveLawData( ) );
    dataStructure->setup (dataFile);

    MeshData             meshData;
    meshData.setup (dataFile, "solid/space_discretization");

    boost::shared_ptr<RegionMesh<LinearTetra> > fullMeshPtr (new RegionMesh<LinearTetra> (  parameters->comm  ) );
    readMesh (*fullMeshPtr, meshData);

    MeshPartitioner< RegionMesh<LinearTetra> > meshPart ( fullMeshPtr, parameters->comm );

    std::string dOrder =  dataFile ( "solid/space_discretization/order", "P1");

    //Mainly used for BCs assembling (Neumann type)
    solidFESpacePtr_Type dFESpace ( new solidFESpace_Type (meshPart, dOrder, 3, parameters->comm) );
    solidETFESpacePtr_Type dETFESpace ( new solidETFESpace_Type (meshPart, & (dFESpace->refFE() ), & (dFESpace->fe().geoMap() ), parameters->comm) );
	

	///////////////////////////////////////////////////////////////////////7
	//HERE I MAKE A CHANGE


	solidFESpace_Type p1FESpace (meshPart, dOrder, 1, parameters->comm);
	

	VectorEpetra p1Scalarfield (p1FESpace.map());
	p1Scalarfield *=0.0;


	int p1DOF = p1Scalarfield.epetraVector().MyLength()/3;

	for(int j(0); j < p1DOF; j++)
{
	int GID = p1Scalarfield.blockMap().GID(j);
	Vector3D coord meshPart.point(GID).coordinates();

	std::cout << "THESE ARE THE NEWEST COORDINATES: " << coord[0] << "\t" << coord[1] << "\t" << coord[2] << std::endl;

}

/////////
///////////////////////////////////////////////HERE THE CHANGE ENDS

    if (verbose)
    {
        std::cout << std::endl;
    }

    std::string timeAdvanceMethod =  dataFile ( "solid/time_discretization/method", "Newmark");

    timeAdvance_Type  timeAdvance ( TimeAdvanceFactory::instance().createObject ( timeAdvanceMethod ) );

    UInt OrderDev = 2;

    //! initialization of parameters of time Advance method:
    if (timeAdvanceMethod == "Newmark")
    {
        timeAdvance->setup ( dataStructure->dataTimeAdvance()->coefficientsNewmark() , OrderDev);
    }

    if (timeAdvanceMethod == "BDF")
    {
        timeAdvance->setup (dataStructure->dataTimeAdvance()->orderBDF() , OrderDev);
    }

    timeAdvance->setTimeStep (dataStructure->dataTime()->timeStep() );
    //timeAdvance->showMe();

    //! #################################################################################
    //! BOUNDARY CONDITIONS
    //! #################################################################################
    vector <ID> compx (1), compy (1), compz (1), compxy (2), compxz (2), compyz (2);
    compx[0] = 0;
    compy[0] = 1, compz[0] = 2;
    compxy[0] = 0;
    compxy[1] = 1;
    compxz[0] = 0;
    compxz[1] = 2;
    compyz[0] = 1;
    compyz[1] = 2;

    BCFunctionBase zero (Private::bcZero);
    BCFunctionBase load (Private::boundaryLoad);


    //! =================================================================================
    //! BC for SymmCube.mesh
    //! =================================================================================
    BCh->addBC ("EdgesIn",      100,  Essential, Component, zero,    compx);
    BCh->addBC ("SymmetryX",    200,  Essential, Component, zero,    compz);
    BCh->addBC ("SymmetryY",    300,  Essential, Component, zero,    compy);

    Teuchos::ParameterList list = * ( Teuchos::getParametersFromXmlFile ( "xmlParameters.xml" ) );
    Real tract = list.get ( "traction", 0.0 );
    Real Rgammaf = list.get ("gammaf", 0.0);
    if (tract != 0 )
    {
        BCh->addBC ("EdgesIn",      600,  Natural,   Component, load, compx);
    }
    BCh->addBC ("edgetwo",      10,  EssentialEdges, Component, zero,    compxz);
    BCh->addBC ("edgetwo",      20,  EssentialEdges, Component, zero,    compxy);
    BCh->addBC ("edgetwo",      30,  EssentialEdges, Component, zero,    compyz);
    //    BCh->addBC ("edgetwo",      10,  EssentialVertices, Component, zero,    compxz);
    //    BCh->addBC ("edgetwo",      20,  EssentialVertices, Component, zero,    compxy);
    //    BCh->addBC ("edgetwo",      30,  EssentialVertices, Component, zero,    compyz);
    //! =================================================================================


    //! 1. Constructor of the structuralSolver
    EMStructuralOperator< RegionMesh<LinearTetra> > solid;

    //    //! 2. Setup of the structuralSolver
    //    solid.setup (dataStructure,
    //                 dFESpace,
    //                 dETFESpace,
    //                 BCh,
    //                 parameters->comm);

    //! 2. Setup of the structuralSolver
    solid.setup (dataStructure,
                 dFESpace,
                 dETFESpace,
                 BCh,
                 parameters->comm);


    solidFESpacePtr_Type gammaFESpace ( new solidFESpace_Type (meshPart, dOrder, 1, parameters->comm) );
    VectorEpetra gammaf ( gammaFESpace -> map() );
    typedef boost::function<Real ( Real const&, Real const&, Real const&, Real const&, ID const& ) > fct_type;
    fct_type fg = & (Private::gf);

    gammaFESpace -> interpolate ( static_cast< FESpace< RegionMesh<LinearTetra>, MapEpetra >::function_Type > ( fg ), gammaf , 0);
    //solidFESpacePtr_Type initd( new solidFESpace_Type (meshPart, dOrder, 3, parameters->comm) );

    fct_type exsol = & (Private::d0);
    vectorPtr_Type initd ( new vector_Type ( dFESpace -> map() ) );
    dFESpace -> interpolate ( static_cast< FESpace< RegionMesh<LinearTetra>, MapEpetra >::function_Type > ( exsol ), *initd , 0);


    if ( dataStructure->solidType() == "neoHookeanActivated" )
    {
        solid.material() -> setGammaf (gammaf);
    }
    //! 3. Setting data from getPot
    solid.setDataFromGetPot (dataFile);


    Real fx = 1.0;
    Real fy = 0.0;
    Real fz = 0.0;
    if ( dataStructure->solidType() == "neoHookeanActivated" )
    {
        solid.material() -> setupFiberVector (fx, fy, fz);
    }



    //! 4. Building system using TimeAdvance class
    double timeAdvanceCoefficient = timeAdvance->coefficientSecondDerivative ( 0 ) / (dataStructure->dataTime()->timeStep() * dataStructure->dataTime()->timeStep() );
    solid.buildSystem (timeAdvanceCoefficient);


    //dataStructure->showMe();
    //! =================================================================================
    //! Temporal data and initial conditions
    //! =================================================================================

    //! 5. Initial data
    Real dt = dataStructure->dataTime()->timeStep();
    // Real T  = dataStructure->dataTime()->endTime();

    vectorPtr_Type rhs (new vector_Type (solid.displacement(), Unique) );
    vectorPtr_Type disp (new vector_Type (solid.displacement(), Unique) );
    vectorPtr_Type vel (new vector_Type (solid.displacement(), Unique) );
    vectorPtr_Type acc (new vector_Type (solid.displacement(), Unique) );

    if (verbose)
    {
        std::cout << "S- initialization ... ";
    }

    std::vector<vectorPtr_Type> uv0;

    if (timeAdvanceMethod == "Newmark")
    {
        uv0.push_back (disp);
        uv0.push_back (vel);
        uv0.push_back (acc);
    }

    vectorPtr_Type initialDisplacement (new vector_Type (solid.displacement(), Unique) );

    if ( !dataStructure->solidType().compare ("secondOrderExponential") )
    {
        dFESpace->interpolate ( static_cast<solidFESpace_Type::function_Type> ( Private::d0 ), *initialDisplacement, 0.0 );
    }

    if (timeAdvanceMethod == "BDF")
    {
        Real tZero = dataStructure->dataTime()->initialTime();

        for ( UInt previousPass = 0; previousPass < timeAdvance->size() ; previousPass++)
        {
            Real previousTimeStep = tZero - previousPass * dt;
            std::cout << "BDF " << previousTimeStep << "\n";
            if ( !dataStructure->solidType().compare ("secondOrderExponential") )
            {
                uv0.push_back (initialDisplacement);
            }
            else
            {
                uv0.push_back (disp);
            }
        }
    }

    //    timeAdvance->setInitialCondition (uv0);

    //    timeAdvance->setTimeStep ( dt );

    //    timeAdvance->updateRHSContribution ( dt );

    //    if ( !dataStructure->solidType().compare ("secondOrderExponential") )
    //    {
    solid.initialize ( initialDisplacement );
    //    }
    //    else
    //    {
    //        solid.initialize ( disp );
    //    }

    MPI_Barrier (MPI_COMM_WORLD);

    if (verbose )
    {
        std::cout << "ok." << std::endl;
    }

    boost::shared_ptr< Exporter<RegionMesh<LinearTetra> > > exporter;

    std::string const exporterType =  dataFile ( "exporter/type", "ensight");
#ifdef HAVE_HDF5
    if (exporterType.compare ("hdf5") == 0)
    {
        exporter.reset ( new ExporterHDF5<RegionMesh<LinearTetra> > ( dataFile, "structure" ) );
    }
    else
#endif
    {
        if (exporterType.compare ("none") == 0)
        {
            exporter.reset ( new ExporterEmpty<RegionMesh<LinearTetra> > ( dataFile, meshPart.meshPartition(), "structure", parameters->comm->MyPID() ) );
        }

        else
        {
            exporter.reset ( new ExporterEnsight<RegionMesh<LinearTetra> > ( dataFile, meshPart.meshPartition(), "structure", parameters->comm->MyPID() ) );
        }
    }

    exporter->setPostDir ( "./" );
    exporter->setMeshProcId ( meshPart.meshPartition(), parameters->comm->MyPID() );

    vectorPtr_Type solidDisp ( new vector_Type (solid.displacement(), exporter->mapType() ) );
    //    vectorPtr_Type solidVel  ( new vector_Type (solid.displacement(), exporter->mapType() ) );
    //    vectorPtr_Type solidAcc  ( new vector_Type (solid.displacement(), exporter->mapType() ) );


    exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, "displacement", dFESpace, solidDisp, UInt (0) );
    exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, "initd", dFESpace, initd, UInt (0) );
    //    exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, "velocity",     dFESpace, solidVel,  UInt (0) );
    //    exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, "acceleration", dFESpace, solidAcc,  UInt (0) );
    if ( dataStructure->solidType() == "neoHookeanActivated" )
    {
        vectorPtr_Type solidgamma  ( new vector_Type ( * ( solid.material() -> gammaf() ), exporter->mapType() ) );
        exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::ScalarField, "gamma_function", gammaFESpace, solidgamma,  UInt (0) );
    }


    exporter->postProcess ( 0 );


    Real normVect;
    normVect =  solid.displacement().norm2();
    std::cout << "The norm 2 of the displacement field is: " << normVect << std::endl;

    //! =============================================================================
    //! Solving loop
    //! =============================================================================

    //! 7. Iterate --> Calling Newton
    solid.iterate ( BCh );

    //        timeAdvance->shiftRight ( solid.displacement() );

    *solidDisp = solid.displacement();
    //*solidVel  = timeAdvance->firstDerivative();
    //*solidAcc  = timeAdvance->secondDerivative();

    exporter->postProcess ( 1.0 );

    normVect =  solid.displacement().norm2();
    std::cout << "The norm 2 of the displacement field is: " << normVect << std::endl;



    //!--------------------------------------------------------------------------------------------------
    MPI_Barrier (MPI_COMM_WORLD);
    //}
}

int
main ( int argc, char** argv )
{

#ifdef HAVE_MPI
    MPI_Init (&argc, &argv);
    boost::shared_ptr<Epetra_MpiComm> Comm (new Epetra_MpiComm ( MPI_COMM_WORLD ) );
    if ( Comm->MyPID() == 0 )
    {
        cout << "% using MPI" << endl;
    }
#else
    boost::shared_ptr<Epetra_SerialComm> Comm ( new Epetra_SerialComm() );
    cout << "% using serial Version" << endl;
#endif


    Structure structure ( argc, argv, Comm );
    structure.run();

#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return returnValue ;
}
