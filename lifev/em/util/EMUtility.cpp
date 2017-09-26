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


#include <lifev/em/util/EMUtility.hpp>


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
    
    if (comm.MyPID() == 0)
    {
        std::cout << "\nEMUtility: createOutputFolder (with name: " << problemFolder << ") ... " << '\r' << std::flush;
    }
    
    if ( problemFolder.compare ("./") )
    {
        problemFolder += "/";

        if ( comm.MyPID() == 0 )
        {
            mkdir ( problemFolder.c_str(), 0777 );
        }
    }
    
    if (comm.MyPID() == 0)
    {
        std::cout << "EMUtility: createOutputFolder (with name: " << problemFolder << ") - done " << '\r' << std::flush;
    }
    
    return problemFolder;
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


void normalize (VectorSmall<3>& v, int component )
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

void orthonormalize (VectorSmall<3>& s, VectorSmall<3>& f, int component )
{
    EMUtility::normalize (s, component);
    s = s - s.dot (f) *  f;
    EMUtility::normalize (s, ++component);
}



Real FLRelationship(Real I4f)
{
	Real lambda = std::sqrt(I4f);
    if ( lambda > 0.87277 && lambda < 1.334)
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
        		   + d1 * std::sin (lambda * l0)
                   + e1 * std::cos (lambda * l0)
                   + d2 * std::sin (2 * lambda * l0)
                   + e2 * std::cos (2 * lambda * l0)
                   + d3 * std::sin (3 * lambda * l0)
                   + e3 * std::cos (3 * lambda * l0);
        return Force;
    }
    else
    {
        return 0.0;
    }
}




} // namespace EMUtility

} // namespace LifeV

