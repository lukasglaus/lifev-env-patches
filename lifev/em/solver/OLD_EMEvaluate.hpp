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

 */

#ifndef EMEVALUATE_H
#define EMEVALUATE_H 1

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

class EMEvaluate
{

public:

    EMEvaluate()
    {
    }

    virtual ~EMEvaluate()
    {
    }

    void evaluate (VectorEpetra& v1, VectorEpetra& v2,
                   FESpace<RegionMesh<LinearTetra>, MapEpetra>& fespace1,
                   FESpace<RegionMesh<LinearTetra>, MapEpetra>& fespace2,
                   RegionMesh<LinearTetra>& finalMesh);

    MatrixSmall<3, 3>& BK (VectorSmall<3>& A, VectorSmall<3>& B,
                           VectorSmall<3>& C, VectorSmall<3>& D);

    Real DetBK (MatrixSmall<3, 3>& Bk);

    MatrixSmall<3, 3>& invBK (MatrixSmall<3, 3>& Bk);

    VectorSmall<3>& TransformPoint (Real x, Real y, Real z,
                                    MatrixSmall<3, 3>& invBk, VectorSmall<3>& A);

    bool checkPoint (VectorSmall<3>& X);

    Real N1 (VectorSmall<3>& X);
    Real N2 (VectorSmall<3>& X);
    Real N3 (VectorSmall<3>& X);
    Real N4 (VectorSmall<3>& X);

};
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////

// v1 known vector on fespace fespace1 to be evaluated on fespace2;
// v2 unknown
void EMEvaluate::evaluate (VectorEpetra& v1, VectorEpetra& v2,
                           FESpace<RegionMesh<LinearTetra>, MapEpetra>& fespace1,
                           FESpace<RegionMesh<LinearTetra>, MapEpetra>& fespace2,
                           RegionMesh<LinearTetra>& finalMesh)
{
    LifeChrono timer;
    timer.start();
    int counter (0);
    if ( v2.map().comm().MyPID() == 0)
    {
        std::cout << "\n++++++++++++++++++++++++++++++++";
        std::cout << "\nEVALUATION ON DIFFERENT FESPACE";
    }
    UInt n2 = v2.epetraVector().MyLength();
    std::set<UInt> nodeLID;
    for (UInt i (0); i < n2; i++)
    {
        nodeLID.insert (i);
    }

    std::cout << "\n num ID: " << nodeLID.size();
    std::set<UInt>::iterator itLID;
    UInt ne1 = fespace1.mesh()->numVolumes();

    for (UInt i (0); i < 8; i++)
    {
        Real x = fespace2.mesh() -> point (i).x();
        Real y = fespace2.mesh() -> point (i).y();
        Real z = fespace2.mesh() -> point (i).z();

        std::cout << "\n i: " << i << ", x: " << x << ", y: " << y << ", z: " << z ;
    }
    //UInt ne2 = fespace1.mesh() ->numVolumes();

    if ( v2.map().comm().MyPID() == 0)
    {
        std::cout << "\nloop on volumes of known mesh: " << ne1 << ", unknown mesh nodes" << n2;
    }

    fespace1.fe().update (fespace1.mesh()->volumeList (2) );


    for (UInt i1LocVol (0); i1LocVol < ne1; i1LocVol++)
    {
        //std::cout << "\niLocalVolume: " << i1LocVol;
        fespace1.fe().update (fespace1.mesh()->volumeList (i1LocVol) );

        UInt eleIDu1 = fespace1.fe().currentLocalId();
        UInt eleGIDu1 = fespace1.fe().currentId();
        UInt nbNode1 = (UInt) fespace1.fe().nbFEDof();
        //      if( v2.map().comm().MyPID() == 0)
        //      {
        //          std::cout << "\niLocVol: " << i1LocVol << ", eleIDu1: " <<eleIDu1 << ", eleGIDu1: " <<eleGIDu1 << ", nbNode1: " << nbNode1;
        //      }

        std::vector<ID> ids1;
        for (UInt i1Node = 0; i1Node < nbNode1; i1Node++)
        {
            ID GID1 = fespace1.mesh() -> element ( eleIDu1 ).point ( i1Node ).id();
            //  std::cout << "\n GID1: " << GID1;
            ids1.push_back (GID1);
        }
        // ids2.size() Assume is 4 for P1
        VectorSmall<3> A, B, C, D;
        A[0] = fespace1.mesh()-> element (eleIDu1).point (0).x();
        A[1] = fespace1.mesh()-> element (eleIDu1).point (0).y();
        A[2] = fespace1.mesh()-> element (eleIDu1).point (0).z();
        B[0] = fespace1.mesh()-> element (eleIDu1).point (1).x();
        B[1] = fespace1.mesh()-> element (eleIDu1).point (1).y();
        B[2] = fespace1.mesh()-> element (eleIDu1).point (1).z();
        C[0] = fespace1.mesh()-> element (eleIDu1).point (2).x();
        C[1] = fespace1.mesh()-> element (eleIDu1).point (2).y();
        C[2] = fespace1.mesh()-> element (eleIDu1).point (2).z();
        D[0] = fespace1.mesh()-> element (eleIDu1).point (3).x();
        D[1] = fespace1.mesh()-> element (eleIDu1).point (3).y();
        D[2] = fespace1.mesh()-> element (eleIDu1).point (3).z();

        MatrixSmall<3, 3> bk (BK (A, B, C, D) );
        MatrixSmall<3, 3> invbk = invBK (bk);

        //      if( v2.map().comm().MyPID() == 0)
        //      {
        //          std::cout << "\n!!!!!!!!!!!!!!!!!    BK !!!!!!!!!!!!!!!!!!!!!!!";
        //          std::cout << "\n" << bk(0,0) <<  ", " << bk(0,1) <<  ", "<< bk(0,2);
        //          std::cout << "\n" << bk(1,0) <<  ", " << bk(1,1) <<  ", "<< bk(1,2);
        //          std::cout << "\n" << bk(2,0) <<  ", " << bk(2,1) <<  ", "<< bk(2,2);
        //          std::cout << "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!";
        //          std::cout << "\n!!!!!!!!!!!!!!!!!    invBK !!!!!!!!!!!!!!!!!!!!!!!";
        //          std::cout << "\n" << invbk(0,0) <<  ", " << invbk(0,1) <<  ", "<< invbk(0,2);
        //          std::cout << "\n" << invbk(1,0) <<  ", " << invbk(1,1) <<  ", "<< invbk(1,2);
        //          std::cout << "\n" << invbk(2,0) <<  ", " << invbk(2,1) <<  ", "<< invbk(2,2);
        //          std::cout << "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!";
        //      }

        for (itLID = nodeLID.begin(); itLID != nodeLID.end(); ++itLID)
        {


            //  std::cout << "\n\tiLocalNodes: " << iLocNodes;
            UInt i2GID = v2.blockMap().GID (*itLID);

            Real xp2 = finalMesh.point (i2GID).x();
            Real yp2 = finalMesh.point (i2GID).y();
            Real zp2 = finalMesh.point (i2GID).z();
            //  std::cout << "\nFrom: x: " << xp2 <<  ", y: "<< yp2 <<  ", z: "<< zp2;

            VectorSmall<3> X = TransformPoint (xp2, yp2, zp2, invbk, A);
            //VectorSmall<3> tmp = TransformPoint(xp2, yp2, zp2, invbk, A);
            //std::cout << "\nFrom: x: " << xp2 <<  ", y: "<< yp2 <<  ", z: "<< zp2;
            //std::cout << "\nTransformed in! X: " << X[0] <<  ", Y: "<< X[1] <<  ", Z: "<< X[2];
            bool inCurrentVol = checkPoint (X);


            //          if(*itLID == 11 && i1LocVol == 16 )
            //          {
            //
            //
            //              if( v2.map().comm().MyPID() == 0)
            //              {
            //                  std::cout << "\n===========================================================";
            //                  std::cout << "\nFound a point! x: " << xp2 <<  ", y: "<< yp2 <<  ", z: "<< zp2;
            //                  std::cout << "\nTransformed in! X: " << X[0] <<  ", Y: "<< X[1] <<  ", Z: "<< X[2];
            //                  std::cout << "\nin element: " << i1LocVol << ". GID: " << i2GID;
            //                  std::cout << "\nTetra was: x1: " << A[0] <<  ", y1: "<< A[1] <<  ", z1: "<< A[2];
            //                  std::cout << "\nTetra was: x2: " << B[0] <<  ", y2: "<< B[1] <<  ", z2: "<< B[2];
            //                  std::cout << "\nTetra was: x3: " << C[0] <<  ", y3: "<< C[1] <<  ", z3: "<< C[2];
            //                  std::cout << "\nTetra was: x4: " << D[0] <<  ", y4: "<< D[1] <<  ", z4: "<< D[2];
            //                  std::cout << "\nOriginal Values: v1: " << v1[ids1[0]] <<  ", v2: "<< v1[ids1[1]]  << ", v3: " << v1[ids1[2]] <<  ", v4: "<< v1[ids1[3]];
            //                  std::cout << "\nBasis: N1: " << N1(X) <<  ", N2: "<< N2(X) <<  ", N3: "<< N3(X) <<  ", N4: "<< N4(X);
            //
            //                  std::cout << "\n==========================================================";
            //              }
            //          }

            if (inCurrentVol)
            {
                counter++;
                //std::vector<ID> ids = sortIDs(A, B, C, D, invbk, ids1);
                //std::cout << "\nLID vol: " << eleIDu1  << ", GID vol: " << eleGIDu1 <<  ", LID vertex: "<<  *itLID <<  ", GID vertex: "<< i2GID << ", ids1 size:" << ids1.size();



                Real aux = v1[ids1[0]] * N1 (X) + v1[ids1[1]] * N2 (X)
                           + v1[ids1[2]] * N3 (X) + v1[ids1[3]] * N4 (X);
                //              if( v2.map().comm().MyPID() == 0)
                //              {
                //                  std::cout << "\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++";
                //                  std::cout << "\nFound a point! x: " << xp2 <<  ", y: "<< yp2 <<  ", z: "<< zp2;
                //                  std::cout << "\nTransformed in! X: " << X[0] <<  ", Y: "<< X[1] <<  ", Z: "<< X[2];
                //                  std::cout << "\nin element: " << i1LocVol << ". GID: " << i2GID;
                //                  std::cout << "\nTetra was: x1: " << A[0] <<  ", y1: "<< A[1] <<  ", z1: "<< A[2];
                //                  std::cout << "\nTetra was: x2: " << B[0] <<  ", y2: "<< B[1] <<  ", z2: "<< B[2];
                //                  std::cout << "\nTetra was: x3: " << C[0] <<  ", y3: "<< C[1] <<  ", z3: "<< C[2];
                //                  std::cout << "\nTetra was: x4: " << D[0] <<  ", y4: "<< D[1] <<  ", z4: "<< D[2];
                //                  std::cout << "\nOriginal Values: v1: " << v1[ids1[0]] <<  ", v2: "<< v1[ids1[1]]  << ", v3: " << v1[ids1[2]] <<  ", v4: "<< v1[ids1[3]];
                //                  //std::cout << "\nSorted Values: v1: " << v1[ids[0]] <<  ", v2: "<< v1[ids[1]]  << ", v3: " << v1[ids[2]] <<  ", v4: "<< v1[ids[3]];
                //                  std::cout << "\nBasis: N1: " << N1(X) <<  ", N2: "<< N2(X) <<  ", N3: "<< N3(X) <<  ", N4: "<< N4(X);
                //                  std::cout << "\nValue: " << aux;
                //                  std::cout << "\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++";
                //              }

                v2[*itLID] = aux;
                nodeLID.erase (itLID);
                inCurrentVol = false;
            }
            //          else
            //          {
            //              std::cout << "\n---------------------------------------------------------------";
            //              std::cout << "\nFound a point! x: " << xp2 <<  ", y: "<< yp2 <<  ", z: "<< zp2;
            //              std::cout << "\nTransformed in! X: " << X[0] <<  ", Y: "<< X[1] <<  ", Z: "<< X[2];
            //              std::cout << "\nin element: " << i1LocVol << ". GID: " << i2GID;
            //              std::cout << "\nTetra was: x1: " << A[0] <<  ", y1: "<< A[1] <<  ", z1: "<< A[2];
            //              std::cout << "\nTetra was: x2: " << B[0] <<  ", y2: "<< B[1] <<  ", z2: "<< B[2];
            //              std::cout << "\nTetra was: x3: " << C[0] <<  ", y3: "<< C[1] <<  ", z3: "<< C[2];
            //              std::cout << "\nTetra was: x4: " << D[0] <<  ", y4: "<< D[1] <<  ", z4: "<< D[2];
            //              std::cout << "\nOriginal Values: v1: " << v1[ids1[0]] <<  ", v2: "<< v1[ids1[1]]  << ", v3: " << v1[ids1[2]] <<  ", v4: "<< v1[ids1[3]];
            //              std::cout << "\nBasis: N1: " << N1(X) <<  ", N2: "<< N2(X) <<  ", N3: "<< N3(X) <<  ", N4: "<< N4(X);
            //              std::cout << "\n---------------------------------------------------------------";
            //
            //          }
            //if()
        }
    }

    timer.stop();
    if ( v2.map().comm().MyPID() == 0)
    {
        std::cout << "\nEVALUATION DONE IN " << timer.diff() << " s, found " << counter << " nodes.";
        std::cout << "\n++++++++++++++++++++++++++++++++";
    }

}

MatrixSmall<3, 3>& EMEvaluate::BK (VectorSmall<3>& A, VectorSmall<3>& B,
                                   VectorSmall<3>& C, VectorSmall<3>& D)
{
    Real x1 = A[0];
    Real y1 = A[1];
    Real z1 = A[2];
    Real x2 = B[0];
    Real y2 = B[1];
    Real z2 = B[2];
    Real x3 = C[0];
    Real y3 = C[1];
    Real z3 = C[2];
    Real x4 = D[0];
    Real y4 = D[1];
    Real z4 = D[2];

    MatrixSmall<3, 3> Bk;
    Bk (0, 0) = x2 - x1;
    Bk (0, 1) = x3 - x1;
    Bk (0, 2) = x4 - x1;
    Bk (1, 0) = y2 - y1;
    Bk (1, 1) = y3 - y1;
    Bk (1, 2) = y4 - y1;
    Bk (2, 0) = z2 - z1;
    Bk (2, 1) = z3 - z1;
    Bk (2, 2) = z4 - z1;

    return Bk;
}

Real EMEvaluate::DetBK (MatrixSmall<3, 3>& Bk)
{
    return Bk (0, 0) * (Bk (1, 1) * Bk (2, 2) - Bk (2, 1) * Bk (1, 2) )
           - Bk (0, 1) * (Bk (1, 0) * Bk (2, 2) - Bk (2, 0) * Bk (1, 2) )
           + Bk (0, 2) * (Bk (1, 0) * Bk (2, 1) - Bk (2, 0) * Bk (1, 1) );
}

MatrixSmall<3, 3>& EMEvaluate::invBK (MatrixSmall<3, 3>& Bk)
{
    Real det = DetBK (Bk);

    MatrixSmall<3, 3> invBk;
    invBk (0, 0) = (  Bk (1, 1) * Bk (2, 2) - Bk (2, 1) * Bk (1, 2) )  / det;
    invBk (1, 0) = - (  Bk (1, 0) * Bk (2, 2) - Bk (2, 0) * Bk (1, 2) )  / det;
    invBk (2, 0) = (  Bk (1, 0) * Bk (2, 1) - Bk (2, 0) * Bk (1, 1) )  / det;
    invBk (0, 1) = - (  Bk (0, 1) * Bk (2, 2) - Bk (2, 1) * Bk (0, 2) )  / det;
    invBk (1, 1) = (  Bk (0, 0) * Bk (2, 2) - Bk (2, 0) * Bk (0, 2) )  / det;
    invBk (2, 1) = - (  Bk (0, 0) * Bk (2, 1) - Bk (2, 0) * Bk (0, 1) )  / det;
    invBk (0, 2) = (  Bk (0, 1) * Bk (1, 2) - Bk (1, 1) * Bk (0, 2) )  / det;
    invBk (1, 2) = - (  Bk (0, 0) * Bk (1, 2) - Bk (1, 0) * Bk (0, 2) )  / det;
    invBk (2, 2) = (  Bk (0, 0) * Bk (1, 1) - Bk (1, 0) * Bk (0, 1) )  / det;

    return invBk;
}

//transform the point to the reference space
VectorSmall<3>& EMEvaluate::TransformPoint (Real x, Real y, Real z,
                                            MatrixSmall<3, 3>& invBk, VectorSmall<3>& A)
{
    VectorSmall<3> X;
    X[0] = invBk (0, 0) * x + invBk (0, 1) * y + invBk (0, 2) * z;
    X[1] = invBk (1, 0) * x + invBk (1, 1) * y + invBk (1, 2) * z;
    X[2] = invBk (2, 0) * x + invBk (2, 1) * y + invBk (2, 2) * z;

    VectorSmall<3> aux;
    aux[0] = invBk (0, 0) * A[0] + invBk (0, 1) * A[1] + invBk (0, 2) * A[2];
    aux[1] = invBk (1, 0) * A[0] + invBk (1, 1) * A[1] + invBk (1, 2) * A[2];
    aux[2] = invBk (2, 0) * A[0] + invBk (2, 1) * A[1] + invBk (2, 2) * A[2];

    X[0] -= aux[0];
    X[1] -= aux[1];
    X[2] -= aux[2];

    return X;
}

//check if the transformed point is inside the reference element
bool EMEvaluate::checkPoint (VectorSmall<3>& X)
{

    Real eps = 1e-12;
    Real x = X[0];
    Real y = X[1];
    Real z = X[2];

    if (x >= -eps && y >= -eps && z >= -eps && (1 - x - y - z) >= -eps)
    {
        return true;
    }
    else
    {
        return false;
    }
}

//basis function for P1 elements
Real EMEvaluate::N1 (VectorSmall<3>& X)
{
    return 1 - X[0] - X[1] - X[2];
}
Real EMEvaluate::N2 (VectorSmall<3>& X)
{
    return X[0];
}
Real EMEvaluate::N3 (VectorSmall<3>& X)
{
    return X[1];
}
Real EMEvaluate::N4 (VectorSmall<3>& X)
{
    return X[2];
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
//@}

}// namespace LifeV

#endif /* HEARTUTILITY_H */

