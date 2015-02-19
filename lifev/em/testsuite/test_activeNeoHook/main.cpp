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
    @brief Simple test using the electromechanical passive constitutive laws


    @date 12-2014
    @author Simone Rossi <simone.rossi@epfl.ch>

    @contributor
    @mantainer Simone Rossi <simone.rossi@epfl.ch>
 */


#include <lifev/core/LifeV.hpp>
#include <lifev/electrophysiology/solver/ElectroETAMonodomainSolver.hpp>
#include <lifev/electrophysiology/solver/IonicModels/IonicMinimalModel.hpp>

#include <lifev/structure/solver/StructuralConstitutiveLawData.hpp>

#include <lifev/structure/solver/StructuralConstitutiveLaw.hpp>
#include <lifev/structure/solver/StructuralOperator.hpp>

#include <lifev/em/solver/mechanics/EMStructuralOperator.hpp>
#include <lifev/em/solver/mechanics/EMStructuralConstitutiveLaw.hpp>
#include <lifev/em/solver/EMETAFunctors.hpp>

#include <lifev/core/filter/ExporterEnsight.hpp>
#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5.hpp>
#endif
#include <lifev/core/filter/ExporterEmpty.hpp>

#include <lifev/eta/fem/ETFESpace.hpp>
#include <lifev/eta/expression/Integrate.hpp>

//#include <lifev/em/solver/mechanics/materials/EMMaterial.hpp>


#include <lifev/bc_interface/3D/bc/BCInterface3D.hpp>
#include <lifev/structure/solver/NeoHookeanMaterialNonLinear.hpp>

using namespace LifeV;



int main (int argc, char** argv)
{

#ifdef HAVE_MPI
    MPI_Init ( &argc, &argv );
#endif

    //start with the communicator
    boost::shared_ptr<Epetra_Comm>  comm ( new Epetra_MpiComm (MPI_COMM_WORLD) );
    if ( comm->MyPID() == 0 )
    {
        cout << "% using MPI" << std::endl;
    }

    //===========================================================
    //===========================================================
    //              READ DATAFILE and CREATE OUTPUT FOLDER
    //===========================================================
    //===========================================================
    GetPot command_line (argc, argv);
    const string data_file_name = command_line.follow ("data", 2, "-f", "--file");
    GetPot dataFile (data_file_name);

    EMData emdata;
    emdata.setup (dataFile);

    //When launching the executable use the flag:  -o OutputFolderName
    std::string problemFolder = EMUtility::createOutputFolder (command_line, *comm);


    //===========================================================
    //===========================================================
    //              LOAD MESH
    //===========================================================
    //===========================================================
    if ( comm->MyPID() == 0 )
    {
        std::cout << "Reading Mesh Name and Path...\n";
    }

    std::string meshName = dataFile ( "solid/space_discretization/mesh_file", "" );
    std::string meshPath = dataFile ( "solid/space_discretization/mesh_dir", "./" );

    typedef RegionMesh<LinearTetra>                         mesh_Type;
    typedef boost::shared_ptr<mesh_Type>                    meshPtr_Type;

    meshPtr_Type localSolidMesh ( new mesh_Type ( comm ) );
    meshPtr_Type fullSolidMesh ( new mesh_Type ( comm ) );

    MeshUtility::loadMesh (localSolidMesh, fullSolidMesh, meshName, meshPath);

    if ( comm->MyPID() == 0 )
    {
        std::cout << " Done!" << std::endl;
    }

    //===========================================================
    //===========================================================
    //              FINITE ELEMENT SPACES
    //===========================================================
    //===========================================================

    //Define the finite element space for bc and exporter
    // and the ET finite element space for assembly
    if ( comm->MyPID() == 0 )
    {
        std::cout << "setup spaces ... ";
    }

    typedef FESpace< RegionMesh<LinearTetra>, MapEpetra >               solidFESpace_Type;
    typedef boost::shared_ptr<solidFESpace_Type>                        solidFESpacePtr_Type;

    typedef ETFESpace< RegionMesh<LinearTetra>, MapEpetra, 3, 3 >       solidETFESpace_Type;
    typedef boost::shared_ptr<solidETFESpace_Type>                      solidETFESpacePtr_Type;



    std::string dOrder =  dataFile ( "solid/space_discretization/order", "P1");
    solidFESpacePtr_Type dFESpace ( new solidFESpace_Type (localSolidMesh, dOrder, 3, comm) );
    solidETFESpacePtr_Type dETFESpace ( new solidETFESpace_Type (localSolidMesh, & (dFESpace->refFE() ), & (dFESpace->fe().geoMap() ), comm) );
    solidFESpacePtr_Type activationFESpace ( new solidFESpace_Type (localSolidMesh, dOrder, 1, comm) );
    solidETFESpacePtr_Type activationETFESpace ( new solidETFESpace_Type (localSolidMesh, & (activationFESpace->refFE() ), & (activationFESpace->fe().geoMap() ), comm) );

    if ( comm->MyPID() == 0 )
    {
        std::cout << " Done!" << std::endl;
    }

    //===========================================================
    //===========================================================
    //              BOUNDARY CONDITIONS
    //===========================================================
    //===========================================================
    if ( comm->MyPID() == 0 )
    {
        std::cout << "setup bc ... ";
    }

    typedef BCHandler                                          bc_Type;
    typedef StructuralOperator< RegionMesh<LinearTetra> >      physicalSolver_Type;
    typedef BCInterface3D< bc_Type, physicalSolver_Type >      bcInterface_Type;
    typedef boost::shared_ptr< bcInterface_Type >              bcInterfacePtr_Type;

    bcInterfacePtr_Type                     solidBC ( new bcInterface_Type() );
    solidBC->createHandler();
    solidBC->fillHandler ( data_file_name, "solid" );
    solidBC->handler()->bcUpdate ( *dFESpace->mesh(), dFESpace->feBd(), dFESpace->dof() );

    if ( comm->MyPID() == 0 )
    {
        std::cout << " Done!" << std::endl;
    }

    //===========================================================
    //===========================================================
    //              SOLID MECHANICS
    //===========================================================
    //===========================================================

    if ( comm->MyPID() == 0 )
    {
        std::cout << "setup constitutive law data ... ";
    }

    //Setup the data of the constitutive law
    boost::shared_ptr<StructuralConstitutiveLawData> dataStructure (new StructuralConstitutiveLawData( ) );
    dataStructure->setup (dataFile);

    if ( comm->MyPID() == 0 )
    {
        std::cout << " Done!" << std::endl;
    }

    //Setup the structural operator
    if ( comm->MyPID() == 0 )
    {
        std::cout << "setup structural operator ... " << std::endl;
    }
    //! 1. Constructor of the structuralSolver
    EMStructuralOperator< RegionMesh<LinearTetra> > solid;

    // solidBC is copied inside the StructuralOperator
    // any changes to it don't affect the M_BCh inside the solver
    // after a change in solidBC add it in the solver!
    solid.setup ( dataStructure, dFESpace, dETFESpace, solidBC -> handler(), comm);
    solid.setDataFromGetPot (dataFile);
    solid.setOutputResStream(problemFolder);
    solid.setOutputIterStream(problemFolder);
    solid.EMMaterial()->setParameters(emdata);
    solid.EMMaterial()->showMaterialParameters();

    solid.EMMaterial() -> setupFiberVector( 0.0, 0.0, 1.0);
    // solid.EMMaterial() -> setupSheetVector( 0.0, 1.0, 0.0);

    solid.setNewtonParameters(dataFile);
    solid.buildSystem (1.0);

    if ( comm->MyPID() == 0 )
    {
        std::cout << " Done!" << std::endl;
    }

    //===========================================================
    //===========================================================
    //              CREATE EXPORTER and SAVE INITIAL SOLUTION
    //===========================================================
    //===========================================================

    boost::shared_ptr< Exporter<RegionMesh<LinearTetra> > > exporter;
    exporter.reset ( new ExporterHDF5<RegionMesh<LinearTetra> > ( dataFile, "structure" ) );
    exporter->setPostDir ( problemFolder );
    exporter->setMeshProcId ( localSolidMesh, comm->MyPID() );
    exporter->addVariable ( ExporterData<RegionMesh<LinearTetra> >::VectorField, "displacement", dFESpace, solid.displacementPtr(), UInt (0) );
    exporter->postProcess ( 0 );

    //===========================================================
    //===========================================================
    //         SOLVE
    //===========================================================
    //===========================================================
    Real dt =  dataFile ( "solid/time_discretization/timestep", 0.1);
    Real endTime =  dataFile ( "solid/time_discretization/endtime", 1.0);
    ID LVFlag =  dataFile ( "solid/boundary_conditions/LV_flag", 0);
    Real LVPreloadPressure =  dataFile ( "solid/boundary_conditions/LV_preload_pressure", 0.);
    bool deformedPressure =  dataFile ( "solid/boundary_conditions/deformed_pressure", 0);
    
    solid.setBCFlag( LVFlag );

    LifeV::VectorEpetra gammaf( activationFESpace->map() );
    // LifeV::VectorEpetra gammas( activationFESpace->map() );
    // LifeV::VectorEpetra gamman( activationFESpace->map() );

    std::cout << "Starting Preload Ramp\n";
    // pressure ramp
    for (Real time (0.0); time < endTime;)
    {
        time += dt;
        if ( time > endTime )
        {
            break;
        }
        std::cout << "----- Preload Time: " << time << std::endl;
        solid.data() -> dataTime() -> updateTime();

        solidBC -> updatePhysicalSolverVariables();

        // solid.bcH() = solidBC -> handler();
        if ( deformedPressure )
        {
            std::cout << "----- Preload Pressure: " << time/endTime*LVPreloadPressure << std::endl;
            solid.setLVPressureBC( time/endTime*LVPreloadPressure );
        }
	solid.iterate ( solidBC -> handler() , deformedPressure );

        exporter->postProcess ( time );
    }


    Real maxShortening   = dataFile ( "solid/activation/max_shortening", -0.2);
    UInt activationSteps = dataFile ( "solid/activation/activation_steps", 20);

    if ( maxShortening != 0. )
    {
        std::cout << "Starting Activation Ramp\n";
        // activation ramp
        for (UInt i (1); i <= activationSteps; ++i )
        {
            std::cout << "----- Activation Step: " << i << std::endl;
            // activation
            gammaf = i*maxShortening/static_cast<Real>(activationSteps);
            // gammas = 1/std::sqrt(time/endTime * maxShortening);
            // gamman = 1/std::sqrt(time/endTime * maxShortening);
            std::cout << "----- gamma_f: " << i*maxShortening/activationSteps << std::endl;
            *(solid.EMMaterial()->fiberActivationPtr())  = gammaf;
            // *(solid.EMMaterial()->sheetActivationPtr())  = gammas;
            // *(solid.EMMaterial()->normalActivationPtr()) = gamman;

            solid.iterate ( solidBC -> handler() , deformedPressure );
        
            exporter->postProcess ( solid.data() -> dataTime() -> time() + i*dt );
        }
    }

    //===========================================================
    //===========================================================
    //              CLOSE EXPORTER
    //===========================================================
    //===========================================================
    exporter -> closeFile();


#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return 0;
}
