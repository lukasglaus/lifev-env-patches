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

#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include "mpi.h"
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_ParameterXMLFileReader.hpp>

#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/core/filter/ParserGmsh.hpp>
#include <lifev/core/filter/ParserINRIAMesh.hpp>
#include <lifev/core/mesh/ConvertBareMesh.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>

#include <lifev/core/mesh/MeshEntity.hpp>

#include <lifev/core/fem/FESpace.hpp>
#include <lifev/eta/fem/ETFESpace.hpp>

#include <lifev/core/filter/GetPot.hpp>
#include <lifev/structure/solver/StructuralConstitutiveLawData.hpp>
#include <lifev/structure/solver/StructuralOperator.hpp>

#include <lifev/structure/solver/HolzapfelOgdenMaterial.hpp>

#include <lifev/core/algorithm/PreconditionerIfpack.hpp>
#include <lifev/core/algorithm/PreconditionerML.hpp>

#include <lifev/core/LifeV.hpp>

#include <boost/shared_ptr.hpp>

#include <iostream>
#include <cstdlib>

namespace
{
typedef LifeV::Real real_t;

static real_t bcZero (const real_t&, const real_t&, const real_t&, const real_t&, const LifeV::ID&)
{
    return 0.;
}

static real_t u0fun (const real_t&, const real_t&, const real_t&, const real_t&, const LifeV::ID&)
{
    return 0.;
}

static real_t f0fun (const real_t&, const real_t&, const real_t&, const real_t& /*z*/, const LifeV::ID& comp)
{
    /*
    real_t alpha = - M_PI/3. * z + M_PI/3. * (1-z);

    if (comp == 0)
        return cos(alpha);
    else if (comp == 1)
        return sin(alpha);
    else
        return 0.;
    */
    return (comp == 0) ? 1.0 : 0.0;
}

static real_t gamma0fun (const real_t&, const real_t&, const real_t&, const real_t&, const LifeV::ID&)
{
    return -0.1;
}

using namespace LifeV;

static bool regIF = (PRECFactory::instance().registerProduct ( "Ifpack", &createIfpack ) );
static bool regML = (PRECFactory::instance().registerProduct ( "ML", &createML ) );

}

class ActiveStructure
{
private:
    typedef boost::shared_ptr<Epetra_Comm>   comm_t;
    typedef Teuchos::ParameterList           param_t;

    typedef LifeV::LinearTetra              shape_t;
    typedef LifeV::BareMesh<shape_t>        baremesh_t;
    typedef LifeV::RegionMesh<shape_t>      mesh_t;
    typedef LifeV::MeshPartitioner<mesh_t>  pmesh_t;

    typedef LifeV::FESpace<mesh_t, LifeV::MapEpetra>          fespace_t;
    typedef LifeV::ETFESpace<mesh_t, LifeV::MapEpetra, 3, 3>  etfespace_t;


    param_t                    params;
    boost::shared_ptr<mesh_t>  mesh;

public:
    ActiveStructure (char* dataxml, comm_t comm)
    {

        LifeV::LifeChrono chrono;
        Real timeToAssemble = 0;
        Real timeToLinearSolve = 0;

        bool ilead = (comm->MyPID() == 0);
        if (ilead)
        {
            std::cout << "== Initialise ActiveStructure ==" << std::endl;
        }
        if (ilead)
        {
            std::cout << "   --- Using datafile: " << dataxml << std::endl;
        }

        // ==========
        // Parameters
        // ----------
        if (ilead)
        {
            std::cout << "   --- START reading parameters ... " << std::endl;
        }
        {
            // XML datafile
            Teuchos::ParameterXMLFileReader reader (dataxml);
            this->params = reader.getParameters();
        }
        if (ilead)
        {
            std::cout << "   ---   END reading parameters\n" << std::endl;
        }

        // ====
        // Mesh
        // ----
        if (ilead)
        {
            std::cout << "   --- START reading mesh ... " << std::endl;
        }
        {
            // Start reading a bare mesh
            baremesh_t baremesh;
            std::string meshfile = this->params.get ("meshfile", "");
            //LifeV::MeshIO::ReadGmshFile(meshfile, baremesh, 1, true);
            LifeV::MeshIO::ReadINRIAMeshFile (baremesh, meshfile, 1, true);
            // Conversion to a region mesh
            boost::shared_ptr<mesh_t> fullmesh (new mesh_t (comm) );
            LifeV::convertBareMesh (baremesh, *fullmesh, true);
            //std::cout << fullmesh->faceList.countElementsWithMarkerID(10, std::equal_to<LifeV::markerID_Type>()) << std::endl;
            // Partition
            fullmesh->updateElementFacets (true, true);
            fullmesh->updateElementRidges (true, true);
            // Retrieve perprocess mesh
            pmesh_t pmesh (fullmesh, comm);
            fullmesh.reset();
            mesh = pmesh.meshPartition();
        }
        if (ilead)
        {
            std::cout << "   ---   END reading mesh\n" << std::endl;
        }


        // =======
        // FESpace
        // -------
        // GetPot datafile
        GetPot datatxt ("data.txt");
        std::string order (datatxt ("solid/space_discretization/order", "P1") );
        boost::shared_ptr<fespace_t> u_space (new fespace_t (mesh, order, 3, comm) );
        boost::shared_ptr<etfespace_t> u_etspace (new etfespace_t (mesh,
                                                                   & (u_space->refFE() ),
                                                                   & (u_space->fe().geoMap() ),
                                                                   comm) );

        // ===================
        // Boundary conditions
        // -------------------
        boost::shared_ptr<LifeV::BCHandler> bcs (new LifeV::BCHandler() );

        LifeV::BCFunctionBase zero (bcZero);

        std::vector<LifeV::ID> compx (1, 0), compy (1, 1), compz (1, 2);

        // cube.msh
        //   10: top,  20: bottom
        //   30: back, 40: front
        //   50: left, 60: right
        bcs->addBC ("SymmetryX", 30, LifeV::Essential, LifeV::Component, zero, compx);
        bcs->addBC ("SymmetryY", 50, LifeV::Essential, LifeV::Component, zero, compy);
        bcs->addBC ("SymmetryZ", 10, LifeV::Essential, LifeV::Component, zero, compz);
        //bcs->addBC("SymmetryZ2", 20, LifeV::Essential, LifeV::Component, zero, compz);
        //bcs->addBC("LockedZ", 10, LifeV::Essential, LifeV::Full, zero, 3);

        typedef LifeV::StructuralConstitutiveLawData datasolid_t;
        typedef LifeV::EMStructuralOperator<mesh_t>    solid_t;

        boost::shared_ptr<datasolid_t> datasolid (new datasolid_t() );
        datasolid->setup (datatxt);

        solid_t solid;
        solid.setup (datasolid, u_space, u_etspace, bcs, comm);
        solid.setDataFromGetPot (datatxt);

        double ta_coeff = 0.0;
        solid.buildSystem (ta_coeff);

        typedef LifeV::VectorEpetra vector_t;
        boost::shared_ptr<vector_t> u0 (new vector_t (solid.displacement(), LifeV::Unique) );
        u_space->interpolate (static_cast<fespace_t::function_Type> (u0fun), *u0, 0.0);

        solid.initialize (u0);

        // Fibers
        {
            vector_t f0 (u_space->map() );
            u_space->interpolate (static_cast<fespace_t::function_Type> (f0fun), f0, 0.0);
            solid.material()->setFiberVector (f0);
        }
        // Sheets
        {
            Real sx = 0.0, sy = 0.0, sz = 1.0;
            solid.material()->setupSheetVector (sx, sy, sz);
        }
        // Activation
        {
            fespace_t gamma_space (mesh, order, 1, comm);
            vector_t  gammaf (gamma_space.map() );
            gamma_space.interpolate (static_cast<fespace_t::function_Type> (gamma0fun), gammaf, 0.0);
            solid.material()->setGammaf (gammaf);
        }

        MPI_Barrier (MPI_COMM_WORLD);

        // Exporter
        typedef LifeV::Exporter<mesh_t> exporter_t;
        boost::shared_ptr<exporter_t> exporter (new LifeV::ExporterHDF5<mesh_t> (datatxt, "structure") );
        exporter->setPostDir ("./");
        exporter->setMeshProcId (mesh, comm->MyPID() );
        boost::shared_ptr<vector_t> u0exp (new vector_t (solid.displacement(), exporter->mapType() ) );
        exporter->addVariable (LifeV::ExporterData<mesh_t>::VectorField, "displacement", u_space, u0exp, UInt (0) );

        //if (ilead) std::cout << "Solid displacement = " << solid.displacement().norm2() << std::endl;
        //exporter->postProcess(0.0);

        // ================
        // Nonlinear solver
        // ----------------

        chrono.reset();
        chrono.start();
        solid.iterate (bcs);
        chrono.stop();
        timeToLinearSolve += chrono.globalDiff ( *comm );

        MPI_Barrier (MPI_COMM_WORLD);

        //if (ilead) std::cout << solid.displacement().norm2() << std::endl;

        *u0exp = solid.displacement();
        //exporter->postProcess(1.0);

        if (ilead)
        {
            std::cout << "Total time to solution: " << timeToLinearSolve << std::endl;
        }
    }

};

int main (int argc, char** argv)
{

#ifdef HAVE_MPI
    MPI_Init (&argc, &argv);
    boost::shared_ptr<Epetra_Comm> comm (new Epetra_MpiComm (MPI_COMM_WORLD) );
#else
    boost::shared_ptr<Epetra_Comm> comm (new Epetra_SerialComm() );
#endif

    // Parse command-line
    if (argc < 2)
    {
        std::cout << "Usage: " << argv[0] << " [datafile.xml]" << std::endl;
        return EXIT_SUCCESS;
    }

    ActiveStructure structure (argv[1], comm);

#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    return EXIT_SUCCESS;
}
