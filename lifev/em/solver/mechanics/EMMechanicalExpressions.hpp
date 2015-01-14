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


#include <lifev/core/LifeV.hpp>
#include <lifev/em/util/EMUtility.hpp>

//w = w
//x = identity matrix
//y = y
//z = z



///////////////////////////////////////////////////////////////////////////
// BASIC DEFINITIONS
///////////////////////////////////////////////////////////////////////////
// I = identity matrix
//matrix
#define _I                 ( value (EMUtility::identity()) )
// \nabla u = gradient of the zlacement
//matrix
#define _Grad_u(y, z, w)      ( grad (y, z, w) )
// F = deformation gradient tensor
//matrix
#define _F(y, z, w)        ( ( _Grad_u(y, z, w) + _I ) )
// FT = transpose of the deformation gradient tensor
//matrix
#define _FT(y, z, w)        ( transpose( _F(y, z, w) ) )
// F^-T = inverse transpose of F
//matrix
#define _FmT(y, z, w)      ( minusT ( _F(y, z, w) ) )
// C = F^TF = Right Cauchy Green tensor
//matrix
#define _C(y, z, w)        ( _FT(y, z, w) * _F(y, z, w) )
// B = FF^T = Left Cauchy Green tensor
//matrix
#define _B(y, z, w)        ( _F(y, z, w) * _FT(y, z, w) )

// dF = \delta \nabla u = variation of F
//matrix
#define _dF   ( grad (phi_j) )
// dF = \delta \nabla u = variation of F
//matrix
#define _dFT   ( transpose( _dF ) )

// dFmTdF = dFmT : dF = Derivative of FmT with respect to F in direction dF
//matrix
#define _dFmTdF(y, z, w)  ( value (-1.0) * ( _FmT(y, z, w) ) * ( _dFT ) * ( _FmT(y, z, w) ) )

// dBdF = dB : dF = Derivative of B with respect to F in direction dF
//matrix
#define _dBdF(y, z, w)  ( ( _F(y, z, w) * _dFT + _dF * _FT(y, z, w) )  )
///////////////////////////////////////////////////////////////////////////
// J terms
///////////////////////////////////////////////////////////////////////////
// J = determinant of F
//scalar
#define _J(y, z, w)      ( det ( _F(y, z, w) ) )// J^{-2/3}
// Jm23 = J^{-2/3}
//scalar
#define _Jm23(y, z, w)   ( pow ( _J(y, z, w), 2 / -3.0 ) )
// dJ = Derivative of J with respect to F = J F^{-T}
//matrix
#define _dJ(y, z, w)     ( _J(y, z, w) * _FmT(y, z, w) )
// dJdF = dJ : dF = Derivative of J in direction dF
//scalar
#define _dJdF(y, z, w)   ( dot( _dJ(y, z, w), _dF ) )
// d2JdF = d2J : dF = Second derivative of J in direction dF
//matrix
#define _d2JdF(y, z, w)  ( _dJdF(y, z, w) * _FmT(y, z, w)    \
                         + _J(y, z, w)    * _dFmTdF(y, z, w) )
// dJm23 = Derivative of J^{-2/3}
//matrix
#define _dJm23(y, z, w)  ( value(-2.0/3.0) * _Jm23(y, z, w) * _FmT(y, z, w) )
// dJm23dF = dJm23 : dF = Derivative of Jm23 with respect to F in direction dF
//scalar
#define _dJm23dF(y, z, w)  ( dot( _dJm23(y, z, w), _dF)  )
// d2Jm23dF = d2Jm23 : dF = Second derivative of Jm23 with respect to F in direction dF
//matrix
#define _d2Jm23dF(y, z, w)  ( value(-2.0/3.0) * ( _Jm23(y, z, w)    *  _dFmTdF(y, z, w) \
                                                + _dJm23dF(y, z, w) *  _FmT(y, z, w)    ) )


// Jm43 = J^{-4/3}
//scalar
#define _Jm43(y, z, w)   ( pow (_J(y, z, w), 4 / -3.0 ) )
// dJm43 = Derivative of J^{-4/3}
//matrix
#define _dJm43(y, z, w)  ( value(-4.0/3.0) * _Jm43(y, z, w) * _FmT(y, z, w) )
// dJm43dF = dJm43 : dF = Derivative of Jm43 with respect to F in direction dF
//scalar
#define _dJm43dF(y, z, w)  ( dot( _dJm43(y, z, w), _dF)  )
// d2Jm43dF = d2Jm43 : dF = Second derivative of Jm43 with respect to F in direction dF
//matrix
#define _d2Jm43dF(y, z, w)  ( value(-4.0/3.0) * _Jm43(y, z, w) * ( _dFmTdF(y, z, w) \
                                              + _dJm43dF(y, z, w) * _FmT(y, z, w) ) )

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
#define _I1(y, z, w)           ( dot( _F(y, z, w), _F(y, z, w)) )
// dI1 = Derivative of I_1 with respect of F = 2 F
//matrix
#define _dI1(y, z, w)          ( value(2.0) * _F(y, z, w) )

// \bar{I}_1 = modified invariant of C
#define _I1bar(y, z, w)        ( _Jm23(y, z, w) * _I1(y, z, w) )
// dI1bar = Derivative of \bar{I}_1 with respect to F
//matrix
#define _dI1bar(y, z, w)       ( _Jm23(y, z, w) * _dI1(y, z, w) + _I1(y, z, w) * _dJm23(y, z, w) )


// dI1dF = dI1 : dF = Derivative of I_1 with respect to F in direction dF
//scalar
#define _dI1dF(y, z, w)        ( dot( _dI1(y, z, w), _dF ) )
// dI1bardF = dI1bar : dF = Derivative of \bar{I}_1 with respect to F in direction dF
//scalar
#define _dI1bardF(y, z, w)     ( dot( _dI1bar(y, z, w), _dF ) )

// d2I1dF = d2I1 : dF = Second derivative of I_1 with respect to F in direction dF
//matrix
#define _d2I1dF       ( value(2.0) * _dF )

// d2I1bardF = d2I1bar : dF = Second derivative of I1bar with respect to F in direction dF
//matrix
#define _d2I1bardF(y, z, w)    ( ( _dJm23dF(y, z, w) * _dI1(y, z, w) ) \
                               + ( _Jm23(y, z, w) * _d2I1dF          ) \
                               + ( _I1(y, z, w) * _d2Jm23dF(y, z, w) ) \
                               + ( _dI1dF(y, z, w) * _dJm23(y, z, w) ) )

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


///////////////////////////////////////////////////////////////////////////
// ISOTROPIC TERMS DEPENDING ON I2
///////////////////////////////////////////////////////////////////////////
// I_2 = first invariant of C
#define _I2(y, z, w)           ( value(0.5) * ( _I1(y, z, w) * _I1(y, z, w) - dot( _C(y, z, w), _C(y, z, w) ) ) )
// dI2 = Derivative of I_2 with respect of F = 2 F
//matrix
#define _dI2(y, z, w)          (  value(2.0) * ( _I1(y, z, w) * _I - _B(y, z, w) ) * _F(y, z, w) )

// \bar{I}_2 = modified invariant of C
#define _I2bar(y, z, w)        ( _Jm43(y, z, w) * _I2(y, z, w) )
// dI2bar = Derivative of \bar{I}_2 with respect to F
//matrix
#define _dI2bar(y, z, w)       ( _Jm43(y, z, w) * _dI2(y, z, w) + _I2(y, z, w) * _dJm43(y, z, w) )


// dI2dF = dI2 : dF = Derivative of I_2 with respect to F in direction dF
//scalar
#define _dI2dF(y, z, w)        ( dot( _dI2(y, z, w), _dF ) )
// dI2bardF = dI2bar : dF = Derivative of \bar{I}_2 with respect to F in direction dF
//scalar
#define _dI2bardF(y, z, w)     ( dot( _dI2bar(y, z, w), _dF ) )

// d2I2dF = d2I2 : dF = Second derivative of I_2 with respect to F in direction dF
//matrix
//#define _d2I2dF(y, z, w)       ( value(2.0) * ( _dI1dF(y, z, w)  *  _F(y, z, w) + _I1(y, z, w) * _dF )     \
//		                       + value(-1.0) * ( _dBdF(y, z, w) * _F(y, z, w)  + _B(y, z, w) * _dF  ) )

#define _d2I2dF(y, z, w)       ( value(2.0) * ( (  _dI1dF(y, z, w) * _F(y, z, w) - _dBdF(y, z, w) * _F(y, z, w) )   \
		                                      + (  _I1(y, z, w) * _I - _B(y, z, w) ) *  _dF )  )

// d2I2bardF = d2I2bar : dF = Second derivative of I2bar with respect to F in direction dF
//matrix
#define _d2I2bardF(y, z, w)    ( ( _dJm43dF(y, z, w) * _dI2(y, z, w) )  \
                               + ( _Jm43(y, z, w) * _d2I2dF(y, z, w) )  \
                               + ( _I2(y, z, w) * _d2Jm43dF(y, z, w) )  \
                               + ( _dI2dF(y, z, w) * _dJm43(y, z, w) ) )

// W2 = Derivative of the energy W with respect to I2bar
//scalar
#define _W2           ( eval (M_materialPtr ->W2(), _F ) )
// dW2 = Second derivative of the energy W with respect to I2bar
//scalar
#define _dW2          ( eval (M_materialPtr ->dW2(), _F ) )

// PI = First Piola Stress tensor (I1 component)
//matrix
//#define _P2           ( _W2 * _dI2bar)

// dPdF =  dP : dF = Derivative of the stress tensor P with respect to F in direction dF
//matrix
//#define _dP1          ( _dW1 * _dI1bardF * _dI1bar + _W1 * _d2I1bardF )


///////////////////////////////////////////////////////////////////////////
// ANISOTROPIC TERMS DEPENDING ON I4
///////////////////////////////////////////////////////////////////////////
//f fiber vector
//_
#define _v0(y, a)            ( value(y, a) )

#define _v(y, z, w, v0)   ( _F(y, z, w) * v0 )

#define _I4(y, z, w, v0)   ( dot( _v(y, z, w, v0), _v(y, z, w, v0) ) )

#define _dI4(y, z, w, v0) ( ( value(2.0) ) * ( outerProduct( _v(y, z, w, v0), v0 ) ) )

#define _I4bar(y, z, w, v0)  ( _Jm23(y, z, w) * _I4(y, z, w, v0) )

#define _dI4bar(y, z, w, v0) ( _Jm23(y, z, w) * _dI4(y, z, w, v0) \
							 + _dJm23(y, z, w) * _I4(y, z, w, v0) )

#define _dv(v0)   ( _dF * v0 )

#define _dI4dF(y, z, w, v0)   ( dot( ( _dI4(y, z, w, v0) ), ( _dF )  ) )

#define _d2I4dF(v0)   ( value(2.0) * outerProduct( _dv(v0), v0 )  )

#define _dI4bardF(y, z, w, v0)   ( dot( _dI4bar(y, z, w, v0), _dF  ) )


#define _d2I4bardF(y, z, w, v0) ( ( _Jm23(y, z, w)     * _d2I4dF(v0)         ) \
                                + ( _dJm23dF(y, z, w)  * _dI4(y, z, w, v0)   ) \
                                + ( _dJm23(y, z, w)    * _dI4dF(y, z, w, v0) ) \
                                + ( _d2Jm23dF(y, z, w) * _I4(y, z, w, v0)    ) )



///////////////////////////////////////////////////////////////////////////
// SHEAR TERMS DEPENDING ON I8fs
///////////////////////////////////////////////////////////////////////////

#define _I8(y, z, w, v0, w0)   ( dot ( _F(y, z, w) * v0, _F(y, z, w) * w0 ) )

#define _dI8(y, z, w, v0, w0)   ( _F(y, z, w) * ( outerProduct( v0, w0 ) + outerProduct( w0, v0 ) ) )
//#define _dI8fs(y, z, w, a, b)   ( _F(y, z, w) )* ( outerProduct( _v0(y, a), _w0(y, b) ) + outerProduct( _w0(y, b), _v0(y, a) ) )

#define _I8bar(y, z, w, v0, w0)            ( _Jm23(y, z, w) * _I8(y, z, w, v0, w0) )

// dI1bar = Derivative of \bar{I}_1 with respect to F
//matrix
#define _dI8bar(y, z, w, v0, w0)       ( _Jm23(y, z, w) * _dI8(y, z, w, v0, w0) + _I8(y, z, w, v0, w0) * _dJm23(y, z, w) )

#define _dI8dF(y, z, w, v0, w0)   ( dot( _dI8(y, z, w, v0, w0), _dF  ) )
#define _dI8bardF(y, z, w, v0, w0)   ( dot( _dI8bar(y, z, w, v0, w0), _dF  ) )

#define _d2I8dF(y, v0, w0)   ( _dF * ( outerProduct( v0, w0 ) + outerProduct( w0, v0 ) ) )

#define _d2I8bardF(y, z, w, v0, w0) ( ( _Jm23(y, z, w) * _d2I8dF(y, v0, w0) )        \
                                      + ( _dJm23dF(y, z, w) * _dI8(y, z, w, v0, w0)  ) \
                                      + ( _dI8dF(y, z, w, v0, w0) * _dJm23(y, z, w)  ) \
                                      + ( _I8(y, z, w, v0, w0)    * _d2Jm23dF(y, z, w) ) )


///////////////////////////////////////////////////////////////////////////
// ACTIVE STRAIN TERMS
///////////////////////////////////////////////////////////////////////////

// F_A = active deformation gradient tensor
//matrix
#define _FA( gf, gs, gn, f0, s0, n0 )        ( _I + gf * outerProduct(f0, f0) + gs * outerProduct(s0, s0) + gn * outerProduct(n0, n0) )
//active function
//scalar
#define _gm( g )  ( value(-1.0) * ( g ) / ( ( g ) + 1 ) )

//scalar
#define _gti( g )  ( pow( ( 1 + ( g ) ), 0.5 ) - value(1.0) )


//scalar
#define _go( g, k )  ( g * ( k + g * k + value(1.0) ) )

// C_A^-1_{anisotropic} =  inverse of the active right Cauchy-Green tensor
//matrix
#define _FAinvA( gf, gs, gn, f0, s0, n0 )        ( _I + _gm( gf ) * outerProduct(f0, f0) \
		 	 	 	 	 	 	 	 	 	 	 	  + _gm( gs ) * outerProduct(s0, s0) \
		 	 	 	 	 	 	 	 	 	 	 	  + _gm( gn ) * outerProduct(n0, n0) )
// C_A^-1_{anisotropic} =  inverse of the active right Cauchy-Green tensor
//matrix
#define _FAinvO( gf, k, f0, s0, n0 )        ( _I + _gm( gf ) * outerProduct(f0, f0) \
		 	 	 	 	 	 	 	 	 	 	 + _go( gf, k ) * outerProduct(s0, s0) \
		 	 	 	 	 	 	 	 	 	 	 + _gm( (k * gf ) ) * outerProduct(n0, n0) )
// C_A^-1_{anisotropic} =  inverse of the active right Cauchy-Green tensor
//matrix
#define _FAinvTI( gf, f0, s0, n0 )       ( _I + _gm( gf ) * outerProduct(f0, f0) \
		 	 	 	 	 	 	 	 	 	  + _gti( gf ) * outerProduct(s0, s0) \
		 	 	 	 	 	 	 	 	 	  + _gti( gf ) * outerProduct(n0, n0) )

//helper macros
#define GET_FAinv(_1,_2,_3,_4,_5,_6,NAME,...) NAME
#define _FAinv(...) GET_FAinv(__VA_ARGS__, _FAinvA, _FAinvO, _FAinvTI ) (__VA_ARGS__)

//active function
//scalar
#define _gfun( g )  ( value(-1.0) * g * ( g + 2 ) / ( g + 1 ) / ( g + 1 ) )


// C_A^-1_{anisotropic} =  inverse of the active right Cauchy-Green tensor
//matrix
#define _CAinvA( gf, gs, gn, f0, s0, n0 )        ( _I + _gfun( gf ) * outerProduct(f0, f0) \
		 	 	 	 	 	 	 	 	 	 	 	 + _gfun( gs ) * outerProduct(s0, s0) \
		 	 	 	 	 	 	 	 	 	 	 	 + _gfun( gn ) * outerProduct(n0, n0) )
// C_A^-1_{transversely isotropic} =  inverse of the active right Cauchy-Green tensor in case of trasversely isotropi activation
//matrix
#define _CAinvTI( gf, f0, s0, n0 )        ( _I + _gfun( gf ) * outerProduct(f0, f0) \
		 	 	 	 	 	 	 	 	 	 	 	 + gf   * outerProduct(s0, s0) \
		 	 	 	 	 	 	 	 	 	 	 	 + gf   * outerProduct(n0, n0) )

//active function
//scalar
#define _gfun2( g, k )  ( g * ( k + g * k + value(1.) ) * ( g + g * k + g * g * k + value(2.0) ) )


// C_A^-1_{orthotropic} =  inverse of the active right Cauchy-Green tensor in case of orthotropic activation
//matrix
#define _CAinvO( gf, k, f0, s0, n0 )        ( _I + _gfun( gf ) * outerProduct(f0, f0) \
		 	 	 	 	 	 	 	 	 	 	 + _gfun2( gf, k ) * outerProduct(s0, s0) \
		 	 	 	 	 	 	 	 	 	 	 + _gfun( ( k * gf ) ) * outerProduct(n0, n0) )

//helper macros
#define GET_CAinv(_1,_2,_3,_4,_5,_6,NAME,...) NAME
#define _CAinv(...) GET_CAinv(__VA_ARGS__, _CAinvA, _CAinvO, _CAinvTI ) (__VA_ARGS__)







#endif /* EMMECHANICALEXPRESSIONS_HPP_ */
