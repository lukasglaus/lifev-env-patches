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
    @brief  Heart Functors for the Luo-Rudy Kinetics

    @date 04−2010
    @author

    @contributor Simone Rossi <simone.rossi@epfl.ch>, Ricardo Ruiz-Baier <ricardo.ruiz@epfl.ch>
    @mantainer Simone Rossi <simone.rossi@epfl.ch>
 */


#ifndef _EMFUNCTORS_H_
#define _EMFUNCTORS_H_

#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <boost/shared_ptr.hpp>
#include <boost/bind.hpp>

#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"


#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/eta/fem/ETFESpace.hpp>
#include <lifev/eta/expression/Integrate.hpp>
#include <lifev/core/fem/FESpace.hpp>

#include <lifev/em/util/EMUtility.hpp>
#include <lifev/em/solver/mechanics/EMElasticityFunctions.hpp>



namespace LifeV
{

//==================================================================
// FORCE LENGTH RELATIONSHIP WITH I4f
//==================================================================
class FLRelationship
{
public:
    typedef Real return_Type;

    return_Type operator() (const VectorSmall<1>& I4f)
    {
        Real i4f = I4f[0];
        Real Force = this->operator() (i4f);
        return Force;
    }

    return_Type operator() (const Real& I4f)
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

            Real Force = d0 / 2 + d1 * std::sin (I4f * l0)
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

    FLRelationship() {}
    FLRelationship (const FLRelationship&) {}
    ~FLRelationship() {}
};




//==================================================================
// FORCE LENGTH RELATIONSHIP WITH GAMMA
//==================================================================
class FLRelationshipGamma
{
public:
    typedef Real return_Type;

    return_Type operator() (const VectorSmall<1>& gamma)
    {
        Real g = gamma[0];
        Real Force = this->operator() (g);
        return Force;
    }

    return_Type operator() (const Real& g)
    {

        if (g > -0.0657788 && g < 0.154989)
        {
            Real d0 = -4.333618335582119e3;
            Real d1 = 2.570395355352195e3;
            Real e1 = -2.051827278991976e3;
            Real d2 = 1.329536116891330e3;
            Real e2 = 0.302216784558222e3;
            Real d3 = 0.104943770305116e3;
            Real e3 = 0.218375174229422e3;
            Real l0 = 1.95;

            Real Force = d0 / 2 + d1 * std::sin ( ( 1.0 + g ) * ( 1.0 + g ) * l0)
                         + e1 * std::cos ( ( 1.0 + g ) * ( 1.0 + g ) * l0)
                         + d2 * std::sin ( 2 * ( 1.0 + g ) * ( 1.0 + g ) * l0)
                         + e2 * std::cos ( 2 * ( 1.0 + g ) * ( 1.0 + g ) * l0)
                         + d3 * std::sin ( 3 * ( 1.0 + g ) * ( 1.0 + g ) * l0)
                         + e3 * std::cos ( 3 * ( 1.0 + g ) * ( 1.0 + g ) * l0);
            return Force;
        }
        else
        {
            return 0.0;
        }
    }

    FLRelationshipGamma() {}
    FLRelationshipGamma (const FLRelationshipGamma&) {}
    ~FLRelationshipGamma() {}
};

//==================================================================
// HEAVISIDE
//==================================================================
class HeavisideFct
{
public:
    typedef Real return_Type;

    return_Type operator() (const VectorSmall<1>& g)
    {
        Real p = g[0];
        Real F = this->operator() (p);
        return F;
    }

    return_Type operator() (const Real& p)
    {
        if (p > 0.0 )
        {
            //cout << "\n\nSono positivo\n\n";
            return p;
        }
        else
        {
            return 0.0;
        }
    }

    HeavisideFct() {}
    HeavisideFct (const HeavisideFct&) {}
    ~HeavisideFct() {}
};

//==================================================================
// Holzapfel Ogden
//==================================================================
class Exp
{
public:
    typedef Real return_Type;

    return_Type operator() (const Real& a)
    {
        if (a != 0)
        {
            return std::exp ( a );
        }
        else
        {
            return 1.0;
        }
    }



    Exp() {}
    Exp (const Exp&) {}
    ~Exp() {}
};

class Exp2
{
public:
    typedef Real return_Type;

    return_Type operator() (const Real& a)
    {
        if (a != 0)
        {
            return std::exp ( a * a );
        }
        else
        {
            return 1.0;
        }
    }



    Exp2() {}
    Exp2 (const Exp2&) {}
    ~Exp2() {}
};

class ShowValue
{
public:
    typedef Real return_Type;

    return_Type operator() (const Real& a)
    {
        std::cout.precision (15);
        std::cout << "value is: " << a << " \n";
        return 1.0;
    }

    return_Type operator() (const VectorSmall<3>& a)
    {
        std::cout << "value is: " << a[0]  << " \n" << a[1]  << " \n" << a[2]  << " \n";
        return 1.0;
    }

    return_Type operator() (const MatrixSmall<3, 3>& a)
    {
        std::cout << "value is\n";
        std::cout << a[0][0]  << " \t" << a[0][1]  << " \t" << a[0][2]  << " \n";
        std::cout << a[1][0]  << " \t" << a[1][1]  << " \t" << a[1][2]  << " \n";
        std::cout << a[2][0]  << " \t" << a[2][1]  << " \t" << a[2][2]  << " \n";

        return 1.0;
    }



    ShowValue() {}
    ShowValue (const ShowValue&) {}
    ~ShowValue() {}
};


////////////////////////////////////////////////////////////////////////////////////////
class orthonormalizeFibers
{
public:
    typedef LifeV::VectorSmall<3> return_Type;

    //    return_Type operator() (const LifeV::VectorSmall<3>& v)
    //    {
    //      auto a(v);
    //      EMUtility::normalize(a);
    //      return a;
    //    }

    return_Type operator() (const LifeV::VectorSmall<3>& v)
    {
        auto a (v);
        EMUtility::normalize (a, M_component);
        return a;
    }

    return_Type operator() (const LifeV::VectorSmall<3>& v, const LifeV::VectorSmall<3>& w )
    {
        auto a (v);
        auto b (w);
        EMUtility::normalize (a);
        EMUtility::orthonormalize (b, a);
        return b;
    }

    orthonormalizeFibers (LifeV::Int component = 0) : M_component (component) {}
    ~orthonormalizeFibers() {}

    LifeV::Int M_component;
};


class CrossProduct
{
public:
    typedef LifeV::VectorSmall<3> return_Type;

    return_Type operator() (const LifeV::VectorSmall<3>& v1, const LifeV::VectorSmall<3>& v2)
    {
        return Elasticity::crossProduct (v1, v2);
    }

    CrossProduct() {}
    ~CrossProduct() {}
};


}



#endif
