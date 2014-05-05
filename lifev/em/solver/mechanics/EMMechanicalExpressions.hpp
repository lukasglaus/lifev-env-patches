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
    @brief File to define the expressions for assembly

    @author Simone Rossi <simone.rossi@epfl.ch>

    @date 04-2014
 */


#ifndef EMMECHANICALEXPRESSIONS_HPP_
#define EMMECHANICALEXPRESSIONS_HPP_


//w = w
//x = identity matrix
//y = y
//z = z



///////////////////////////////////////////////////////////////////////////
// BASIC DEFINITIONS
///////////////////////////////////////////////////////////////////////////
// I = identity matrix
//matrix
#define _I(x)                 ( value (x) )
// \nabla u = gradient of the zlacement
//matrix
#define _Grad_u(y, z, w)      ( grad (y, z, w) )
// F = deformation gradient tensor
//matrix
#define _F(x, y, z, w)        ( _Grad_u(y, z, w) + _I(x) )
// F^-T = inverse transpose of F
//matrix
#define _FmT(x, y, z, w)      ( minusT (_F(x, y, z, w)) )
// C = F^TF = Right Cauchy Green tensor
//matrix
#define _C(x, y, z, w)        ( transpose(_F(x, y, z, w)) * _F(x, y, z, w) )

// dF = \delta \nabla u = variation of F
//matrix
#define _dF   ( grad (phi_j) )
// dFmTdF = dFmT : dF = Derivative of FmT with respect to F in direction dF
//matrix
#define _dFmTdF(x, y, z, w)  ( value (-1.0) * _FmT(x, y, z, w) * transpose (_dF) * _FmT(x, y, z, w) )


///////////////////////////////////////////////////////////////////////////
// J terms
///////////////////////////////////////////////////////////////////////////
// J = determinant of F
//scalar
#define _J(x, y, z, w)      ( det (_F(x, y, z, w)) )// J^{-2/3}
// Jm23 = J^{-2/3}
//scalar
#define _Jm23(x, y, z, w)   ( pow (_J(x, y, z, w), -2. / 3) )
// dJ = Derivative of J with respect to F = J F^{-T}
//matrix
#define _dJ(x, y, z, w)     ( _J(x, y, z, w) * _FmT(x, y, z, w) )
// dJdF = dJ : dF = Derivative of J in direction dF
//scalar
#define _dJdF(x, y, z, w)   ( dot( _dJ(x, y, z, w), _dF) )
// d2JdF = d2J : dF = Second derivative of J in direction dF
//matrix
#define _d2JdF(x, y, z, w)  ( _dJdF(x, y, z, w) * _FmT(x, y, z, w) \
		                                   + _J(x, y, z, w) * _dFmTdF(x, y, z, w) )
// dJm23 = Derivative of J^{-2/3}
//matrix
#define _dJm23(x, y, z, w)  ( value(-2.0/3.0) * _Jm23(x, y, z, w) * _FmT(x, y, z, w) )
// dJm23dF = dJm23 : dF = Derivative of Jm23 with respect to F in direction dF
//scalar
#define _dJm23dF(x, y, z, w)  ( dot( _dJm23(x, y, z, w), _dF)  )
// d2Jm23dF = d2Jm23 : dF = Second derivative of Jm23 with respect to F in direction dF
//matrix
#define _d2Jm23dF(x, y, z, w)  ( value(-2.0/3.0) * _Jm23(x, y, z, w) * ( _dFmTdF(x, y, z, w) \
		                                      + _dJm23dF(x, y, z, w) * _FmT(x, y, z, w) ) )

///////////////////////////////////////////////////////////////////////////
// VOLUMETRIC PART
///////////////////////////////////////////////////////////////////////////
// Wvol = Derivative of the energy with respect to J
//  scalar
#define _Wvol ( eval (M_materialPtr ->Wvol(), _F) )
// dWvol = Second derivative of the energy with respect to J
// scalar
#define _dWvol ( eval (M_materialPtr -> dWvol(), _F) )
// K = bulk modulus
// scalar
#define _K   ( value(M_bulkModulus) )
// Pvol = volumetric part of the stress tensor
// matrix
#define _Pvol ( _K * _Wvol * _dJ )
// dPvol = Derivative of the volumetric part of the stress tensor in direction dF
// matrix
#define  _dPvol ( _K * ( _Wvol * _d2JdF + _dWvol * _dJdF * _dJ ) )
///////////////////////////////////////////////////////////////////////////
// ISOTROPIC TERMS DEPENDING ON I1
///////////////////////////////////////////////////////////////////////////
// I_1 = first invariant of C
#define _I1(x, y, z, w)           ( dot( _F(x, y, z, w), _F(x, y, z, w)) )
// dI1 = Derivative of I_1 with respect of F = 2 F
//matrix
#define _dI1(x, y, z, w)          ( value(2.0) * _F(x, y, z, w) )

// \bar{I}_1 = modified invariant of C
#define _I1bar(x, y, z, w)        ( _Jm23(x, y, z, w) * _I1(x, y, z, w) )
// dI1bar = Derivative of \bar{I}_1 with respect to F
//matrix
#define _dI1bar(x, y, z, w)       ( _Jm23(x, y, z, w) * _dI1(x, y, z, w) + _I1(x, y, z, w) * _dJm23(x, y, z, w) )


// dI1dF = dI1 : dF = Derivative of I_1 with respect to F in direction dF
//scalar
#define _dI1dF(x, y, z, w)        ( dot( _dI1(x, y, z, w), _dF ) )
// dI1bardF = dI1bar : dF = Derivative of \bar{I}_1 with respect to F in direction dF
//scalar
#define _dI1bardF(x, y, z, w)     ( dot( _dI1bar(x, y, z, w), _dF ) )

// d2I1dF = d2I1 : dF = Second derivative of I_1 with respect to F in direction dF
//matrix
#define _d2I1dF       ( value(2.0) * _dF )

// d2I1bardF = d2I1bar : dF = Second derivative of I1bar with respect to F in direction dF
//matrix
#define _d2I1bardF(x, y, z, w)    ( _dJm23dF(x, y, z, w) * _dI1(x, y, z, w) \
		                                         + _Jm23(x, y, z, w) * _d2I1dF                            \
		                                         + _I1(x, y, z, w) * _d2Jm23dF(x, y, z, w) \
		                                         + _dI1dF(x, y, z, w) * _dJm23(x, y, z, w) )

// WI = Derivative of the energy W with respect to I1bar
//scalar
#define _W1           ( eval (M_materialPtr ->W1(), _F ) )
// dWI = Second derivative of the energy W with respect to I1bar
//scalar
#define _dW1          ( eval (M_materialPtr ->dW1(), _F ) )

// PI = First Piola Stress tensor (I1 component)
//matrix
#define _P1           ( _W1 * _dI1bar)

// dPdF =  dP : dF = Derivative of the stress tensor P with respect to F in direction dF
//matrix
#define _dP1          ( _dW1 * _dI1bardF * _dI1bar + _W1 * _d2I1bardF )





#endif /* EMMECHANICALEXPRESSIONS_HPP_ */
