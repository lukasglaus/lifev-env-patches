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
#include <lifev/core/fem/GradientRecovery.hpp>
#include <lifev/core/array/MatrixSmall.hpp>



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



Real FLRelationship(Real I4f)
{
    if (I4f > 0.87277 && I4f < 1.334)
    {
        Real d0 = -4.333618335582119e3;
        Real d1 = 2.570395355352195e3;
        Real e1 = -2.051827278991976e3;
        Real d2 = 1.329536116891330e3;
        Real e2 = 0.302216784558222e3;
        Real d3 = 0.104943770305116e3;
        Real e3 = 0.218375174229422e3;
        Real l0 = 1.95;

        Real Force = d0 / 2
        		   + d1 * std::sin (I4f * l0)
                   + e1 * std::cos (I4f * l0)
                   + d2 * std::sin (2 * I4f * l0)
                   + e2 * std::cos (2 * I4f * l0)
                   + d3 * std::sin (3 * I4f * l0)
                   + e3 * std::cos (3 * I4f * l0);
        return Force;
    }
    else
    {
        return 0.0;
    }
}

template<typename DispVectorPtr, typename FESpaceType>
void computeZZGradient(VectorEpetra& displacement, std::vector<DispVectorPtr> gradientPtr, boost::shared_ptr<FESpaceType>  dFESpace)
{
    gradientPtr[0].reset( &GradientRecovery::ZZGradient(dFESpace, displacement, 0) );
    gradientPtr[1].reset( &GradientRecovery::ZZGradient(dFESpace, displacement, 1) );
    gradientPtr[2].reset( &GradientRecovery::ZZGradient(dFESpace, displacement, 2) );
}


template< typename FESpaceType >
void computeI4 ( VectorEpetra& I4, VectorEpetra& displacement, VectorEpetra& fibers, boost::shared_ptr<FESpaceType> dFESpace )
{

//    VectorEpetra sx = GradientRecovery::ZZGradient (dFESpace, displacement, 0);
//    VectorEpetra sy = GradientRecovery::ZZGradient (dFESpace, displacement, 1);
//    VectorEpetra sz = GradientRecovery::ZZGradient (dFESpace, displacement, 2);
    std::vector<boost::shared_ptr<VectorEpetra> >gradientPtr(3);
    EMUtility::computeZZGradient(displacement, gradientPtr, dFESpace);

    I4 *= 0.0;
    Int nLocalDof = I4.epetraVector().MyLength();
    //Int nComponentLocalDof = nLocalDof / 3;
    for (int k (0); k < nLocalDof; k++)
    {
        UInt iGID = ( *gradientPtr[0] ).blockMap().GID (k);
        UInt jGID = ( *gradientPtr[0] ).blockMap().GID (k + nLocalDof);
        UInt kGID = ( *gradientPtr[0] ).blockMap().GID (k + 2 * nLocalDof);

        Real fx, fy, fz;
        Real F11 = ( *gradientPtr[0] ) [iGID] + 1.0;
        Real F12 = ( *gradientPtr[1] ) [iGID];
        Real F13 = ( *gradientPtr[2] ) [iGID];
        Real F21 = ( *gradientPtr[0] ) [jGID];
        Real F22 = ( *gradientPtr[1] ) [jGID] + 1.0;
        Real F23 = ( *gradientPtr[2] ) [jGID];
        Real F31 = ( *gradientPtr[0] ) [kGID];
        Real F32 = ( *gradientPtr[1] ) [kGID];
        Real F33 = ( *gradientPtr[2] ) [kGID] + 1.0;
        fx = F11 * fibers[iGID];
        fx += ( F12 * fibers[jGID] );
        fx += ( F13 * fibers[kGID] );
        fy = F21 * fibers[iGID];
        fy += ( F22 * fibers[jGID] );
        fy += ( F23 * fibers[kGID] );
        fz = F31 * fibers[iGID];
        fz += ( F32 * fibers[jGID] );
        fz += ( F33 * fibers[kGID] );
        Real J = F11 * (F22 * F33 - F32 * F23) - F22 * (F21 * F33 - F31 * F23) + F33 * (F21 * F32 - F31 * F22);

        I4[iGID] = fx * fx + fy * fy + fz * fz;
//        I4[iGID] *= std::pow (J, -2.0 / 3.0);

    }

}


} // namespace EMUtility

} // namespace LifeV

#endif /* HEARTUTILITY_H */

