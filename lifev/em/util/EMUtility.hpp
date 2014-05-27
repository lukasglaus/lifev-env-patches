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
    @brief Utilities

    @contributor Simone Rossi <simone.rossi@epfl.ch>
    @maintainer Simone Rossi <simone.rossi@epfl.ch>

    This file contains a set of base utilities used to read vectorial data (mainly fiber
    and sheet directions) from different formats to VectorEpetra objects.
 */

#ifndef EMUTILITY_H
#define EMUTILITY_H 1

#include <lifev/core/LifeV.hpp>
//#include <lifev/core/array/VectorEpetra.hpp>
//#include <lifev/core/array/MatrixEpetra.hpp>
//#include <lifev/core/array/MapEpetra.hpp>
//#include <lifev/core/filter/ExporterEnsight.hpp>
#include <lifev/core/filter/ExporterHDF5.hpp>
//#include <lifev/core/filter/ExporterEmpty.hpp>
//#include <lifev/core/filter/Exporter.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>
#include <lifev/core/fem/FESpace.hpp>
#include <Epetra_MpiComm.h>
#include <lifev/core/filter/ExporterHDF5.hpp>
#include <lifev/core/filter/GetPot.hpp>
#include <boost/shared_ptr.hpp>



#include <sys/stat.h>

namespace LifeV
{

// Predeclaration

namespace EMUtility
{
//! EMUtility -
/*!
 *  @author(s) Simone Rossi
 *
 *
 */

std::string createOutputFolder (GetPot& command_line, Epetra_Comm& comm)
{
    //  std::string problemFolder = "Output";//command_line.follow("Output", 2, "-o","--output");//) command_line.follow ( "Output", 2, "-o", "--output" );
    std::string problemFolder = command_line.follow ("Output", 2, "-o", "--output"); //) command_line.follow ( "Output", 2, "-o", "--output" );
    // Create the problem folder
    std::cout << "EMU - creating output folder: " << problemFolder << "\n";
    if ( problemFolder.compare ("./") )
    {
        problemFolder += "/";

        if ( comm.MyPID() == 0 )
        {
            mkdir ( problemFolder.c_str(), 0777 );
        }
    }
    return problemFolder;
}

template<class Mesh>
void setupExporter ( ExporterHDF5<Mesh>& exporter,
                     boost::shared_ptr<Mesh> localMeshPtr,
                     boost::shared_ptr<Epetra_Comm> commPtr,
                     std::string fileName,
                     std::string folder)
{
    exporter.setMeshProcId (localMeshPtr, commPtr->MyPID() );
    exporter.setPrefix (fileName);
    exporter.exportPID (localMeshPtr, commPtr);
    exporter.setPostDir (folder);
}









void EpetraSqrt ( VectorEpetra& vec)
{
    Int size = vec.epetraVector().MyLength();
    for (int j (0); j < size; j++ )
    {
        int gid = vec.blockMap().GID (j);
        vec[gid] = std::sqrt (vec[gid]);
    }
}


MatrixSmall<3, 3> identity()
{
    MatrixSmall<3, 3> I;
    I (0, 0) = 1.0;
    I (0, 1) = 0.0;
    I (0, 2) = 0.0;
    I (1, 0) = 0.0;
    I (1, 1) = 1.0;
    I (1, 2) = 0.0;
    I (2, 0) = 0.0;
    I (2, 1) = 0.0;
    I (2, 2) = 1.0;
    return I;
}



//void normalize(VectorSmall<3>& v, Real component = 0)
//{
//  normalize(v, static_cast<int>(component) );
//}

void normalize (VectorSmall<3>& v, int component = 0)
{
    LifeV::Real norm = std::sqrt (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    if ( norm >= 1e-13 )
    {
        v[0] = v[0] / norm;
        v[1] = v[1] / norm;
        v[2] = v[2] / norm;
    }
    else
    {
        v *= 0.0;
        v[component] = 1.0;
    }
}

void orthonormalize (VectorSmall<3>& s, VectorSmall<3>& f, int component = 1)
{
    EMUtility::normalize (f);
    s = s - s.dot (f) *  f;
    EMUtility::normalize (s, component);
}




} // namespace EMUtility

} // namespace LifeV

#endif /* HEARTUTILITY_H */

