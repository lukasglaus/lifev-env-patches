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
#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/array/MapEpetra.hpp>
#include <lifev/core/filter/ExporterEnsight.hpp>
#include <lifev/core/filter/ExporterHDF5.hpp>
#include <lifev/core/filter/ExporterEmpty.hpp>
#include <lifev/core/filter/Exporter.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>
#include <lifev/core/fem/FESpace.hpp>

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

} // namespace EMUtility

} // namespace LifeV

#endif /* HEARTUTILITY_H */

