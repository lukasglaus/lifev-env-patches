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
#define _FT(F)        ( transpose( F ) )
// F^-T = inverse transpose of F
//matrix
#define _FmT(F)      ( minusT ( F ) )
// C = F^TF = Right Cauchy Green tensor
//matrix
#define _C(F)        ( _FT(F) * F )
// B = FF^T = Left Cauchy Green tensor
//matrix
#define _B(F)        ( F * _FT(F) )

//// F = deformation gradient tensor
////matrix
//#define _F(y, z, w)        ( ( _Grad_u(y, z, w) + _I ) )
//// FT = transpose of the deformation gradient tensor
////matrix
//#define _FT(y, z, w)        ( transpose( _F(y, z, w) ) )
//// F^-T = inverse transpose of F
////matrix
//#define _FmT(y, z, w)      ( minusT ( _F(y, z, w) ) )
//// C = F^TF = Right Cauchy Green tensor
////matrix
//#define _C(y, z, w)        ( _FT(y, z, w) * _F(y, z, w) )
//// B = FF^T = Left Cauchy Green tensor
////matrix
//#define _B(y, z, w)        ( _F(y, z, w) * _FT(y, z, w) )

// dF = \delta \nabla u = variation of F
//matrix
#define _dF   ( grad (phi_j) )
// dF = \delta \nabla u = variation of F
//matrix
#define _dFT(dF)   ( transpose( dF ) )

// dFmTdF = dFmT : dF = Derivative of FmT with respect to F in direction dF
//matrix
#define _dFmTdF(F, dF)  ( value (-1.0) * ( _FmT(F) ) * ( _dFT(dF) ) * ( _FmT(F) ) )

// dBdF = dB : dF = Derivative of B with respect to F in direction dF
//matrix
#define _dBdF(F, dF)  ( ( F * _dFT(dF) + dF * _FT(F) )  )
///////////////////////////////////////////////////////////////////////////
// J terms
///////////////////////////////////////////////////////////////////////////
// J = determinant of F
//scalar
#define _J(F)      ( det ( F ) )// J^{-2/3}
// Jm23 = J^{-2/3}
//scalar
#define _Jm23(F)   ( pow ( _J(F), 2 / -3.0 ) )
// dJ = Derivative of J with respect to F = J F^{-T}
//matrix
#define _dJ(F)     ( _J(F) * _FmT(F) )
// dJdF = dJ : dF = Derivative of J in direction dF
//scalar
#define _dJdF(F, dF)   ( dot( _dJ(F), _dF ) )
// d2JdF = d2J : dF = Second derivative of J in direction dF
//matrix
#define _d2JdF(F, dF)  ( _dJdF(F, dF) * _FmT(F)    \
    			   + _J(F)    * _dFmTdF(F, dF) )
// dJm23 = Derivative of J^{-2/3}
//matrix
#define _dJm23(F)  ( value(-2.0/3.0) * _Jm23(F) * _FmT(F) )
// dJm23dF = dJm23 : dF = Derivative of Jm23 with respect to F in direction dF
//scalar
#define _dJm23dF(F, dF)  ( dot( _dJm23(F), dF)  )
// d2Jm23dF = d2Jm23 : dF = Second derivative of Jm23 with respect to F in direction dF
//matrix
#define _d2Jm23dF(F, dF)  ( value(-2.0/3.0) * ( _Jm23(F)    *  _dFmTdF(F, dF) \
                                              + _dJm23dF(F, dF) *  _FmT(F)    ) )


// Jm43 = J^{-4/3}
//scalar
#define _Jm43(F)   ( pow (_J(F), 4 / -3.0 ) )
// dJm43 = Derivative of J^{-4/3}
//matrix
#define _dJm43(F)  ( value(-4.0/3.0) * _Jm43(F) * _FmT(F) )
// dJm43dF = dJm43 : dF = Derivative of Jm43 with respect to F in direction dF
//scalar
#define _dJm43dF(F, dF)  ( dot( _dJm43(F), dF)  )
// d2Jm43dF = d2Jm43 : dF = Second derivative of Jm43 with respect to F in direction dF
//matrix
#define _d2Jm43dF(F, dF)  ( value(-4.0/3.0) * _Jm43(F) * ( _dFmTdF(F, dF) \
                                            + _dJm43dF(F, dF) * _FmT(F) ) )

///////////////////////////////////////////////////////////////////////////
// VOLUMETRIC PART
///////////////////////////////////////////////////////////////////////////
//// Wvol = Derivative of the energy with respect to J
////  scalar
//#define _Wvol ( eval (M_materialPtr ->Wvol(), _F) )
//// dWvol = Second derivative of the energy with respect to J
//// scalar
//#define _dWvol ( eval (M_materialPtr -> dWvol(), _F) )
//// K = bulk modulus
//// scalar
//#define _K   ( value(M_bulkModulus) )
//// Pvol = volumetric part of the stress tensor
//// matrix
//#define _Pvol ( _K * _Wvol * _dJ )
//// dPvol = Derivative of the volumetric part of the stress tensor in direction dF
//// matrix
//#define  _dPvol ( _K * ( _Wvol * _d2JdF + _dWvol * _dJdF * _dJ ) )
///////////////////////////////////////////////////////////////////////////
// ISOTROPIC TERMS DEPENDING ON I1
///////////////////////////////////////////////////////////////////////////
// I_1 = first invariant of C
#define _I1(F)           ( dot( F, F ) )
// dI1 = Derivative of I_1 with respect of F = 2 F
//matrix
#define _dI1(F)          ( value(2.0) * F )

// \bar{I}_1 = modified invariant of C
#define _I1bar(F)        ( _Jm23(F) * _I1(F) )
// dI1bar = Derivative of \bar{I}_1 with respect to F
//matrix
#define _dI1bar(F)       ( _Jm23(F) * _dI1(F) + _I1(F) * _dJm23(F) )


// dI1dF = dI1 : dF = Derivative of I_1 with respect to F in direction dF
//scalar
#define _dI1dF(F, dF)        ( dot( _dI1(F), dF ) )
// dI1bardF = dI1bar : dF = Derivative of \bar{I}_1 with respect to F in direction dF
//scalar
#define _dI1bardF(F, dF)     ( dot( _dI1bar(F), dF ) )

// d2I1dF = d2I1 : dF = Second derivative of I_1 with respect to F in direction dF
//matrix
#define _d2I1dF(dF)       ( value(2.0) * dF )

// d2I1bardF = d2I1bar : dF = Second derivative of I1bar with respect to F in direction dF
//matrix
#define _d2I1bardF(F, dF)    ( ( _dJm23dF(F, dF) * _dI1(F)          ) \
						     + ( _Jm23(F)        * _d2I1dF(dF)      ) \
						     + ( _I1(F)          * _d2Jm23dF(F, dF) ) \
						     + ( _dI1dF(F, dF)   * _dJm23(F) )      )

//// WI = Derivative of the energy W with respect to I1bar
////scalar
//#define _W1           ( eval (M_materialPtr ->W1(), _F ) )
//// dWI = Second derivative of the energy W with respect to I1bar
////scalar
//#define _dW1          ( eval (M_materialPtr ->dW1(), _F ) )
//
//// PI = First Piola Stress tensor (I1 component)
////matrix
//#define _P1           ( _W1 * _dI1bar)
//
//// dPdF =  dP : dF = Derivative of the stress tensor P with respect to F in direction dF
////matrix
//#define _dP1          ( _dW1 * _dI1bardF * _dI1bar + _W1 * _d2I1bardF )


///////////////////////////////////////////////////////////////////////////
// ISOTROPIC TERMS DEPENDING ON I2
///////////////////////////////////////////////////////////////////////////
// I_2 = first invariant of C
#define _I2(F)           ( value(0.5) * ( _I1(F) * _I1(F) - dot( _C(F), _C(F) ) ) )
// dI2 = Derivative of I_2 with respect of F = 2 F
//matrix
#define _dI2(F)          (  value(2.0) * ( _I1(F) * _I - _B(F) ) * F )

// \bar{I}_2 = modified invariant of C
#define _I2bar(F)        ( _Jm43(F) * _I2(F) )
// dI2bar = Derivative of \bar{I}_2 with respect to F
//matrix
#define _dI2bar(F)       ( _Jm43(F) * _dI2(F) + _I2(F) * _dJm43(F) )


// dI2dF = dI2 : dF = Derivative of I_2 with respect to F in direction dF
//scalar
#define _dI2dF(F, dF)        ( dot( _dI2(F), dF ) )
// dI2bardF = dI2bar : dF = Derivative of \bar{I}_2 with respect to F in direction dF
//scalar
#define _dI2bardF(F, dF)     ( dot( _dI2bar(F), dF ) )

// d2I2dF = d2I2 : dF = Second derivative of I_2 with respect to F in direction dF
//matrix
//#define _d2I2dF(y, z, w)       ( value(2.0) * ( _dI1dF(y, z, w)  *  _F(y, z, w) + _I1(y, z, w) * _dF )     \
//		                       + value(-1.0) * ( _dBdF(y, z, w) * _F(y, z, w)  + _B(y, z, w) * _dF  ) )

#define _d2I2dF(F, dF)       ( value(2.0) * ( (  _dI1dF(F, dF) * F - _dBdF(F, dF) * F )   \
										    + (  _I1(F) * _I - _B(F) ) *  dF )  )

// d2I2bardF = d2I2bar : dF = Second derivative of I2bar with respect to F in direction dF
//matrix
#define _d2I2bardF(F, dF)  ( ( _dJm43dF(F, dF) * _dI2(F)          )  \
						   + ( _Jm43(F)        * _d2I2dF(F, dF)   )  \
						   + ( _I2(F)          * _d2Jm43dF(F, dF) )  \
						   + ( _dI2dF(F, dF)   * _dJm43(F)        )  )

//// W2 = Derivative of the energy W with respect to I2bar
////scalar
//#define _W2           ( eval (M_materialPtr ->W2(), _F ) )
//// dW2 = Second derivative of the energy W with respect to I2bar
////scalar
//#define _dW2          ( eval (M_materialPtr ->dW2(), _F ) )

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

#define _v(F, v0)   ( F * v0 )

#define _I4(F, v0)   ( dot( _v(F, v0), _v(F, v0) ) )

#define _dI4(F, v0) ( ( value(2.0) ) * ( outerProduct( _v(F, v0), v0 ) ) )

#define _I4bar(F, v0)  ( _Jm23(F) * _I4(F, v0) )

#define _dI4bar(F, v0) ( _Jm23(F) * _dI4(F, v0) \
					   + _dJm23(F) * _I4(F, v0) )

#define _dv(v0, dF)   ( dF * v0 )

#define _dI4dF(F, v0, dF)   ( dot( ( _dI4(F, v0) ), ( dF )  ) )

#define _d2I4dF(v0, dF)   ( value(2.0) * outerProduct( _dv(v0, dF), v0 )  )

#define _dI4bardF(F, v0, dF)   ( dot( _dI4bar(F, v0), dF  ) )


#define _d2I4bardF(F, v0, dF) ( ( _Jm23(F)         * _d2I4dF(v0, dF)   ) \
						      + ( _dJm23dF(F, dF)  * _dI4(F, v0)       ) \
						      + ( _dJm23(F)        * _dI4dF(F, v0, dF) ) \
						      + ( _d2Jm23dF(F, dF) * _I4(F, v0)        ) )



///////////////////////////////////////////////////////////////////////////
// SHEAR TERMS DEPENDING ON I8fs
///////////////////////////////////////////////////////////////////////////

#define _I8(F, v0, w0)   ( dot ( F * v0, F * w0 ) )

#define _dI8(F, v0, w0)   ( F * ( outerProduct( v0, w0 ) + outerProduct( w0, v0 ) ) )
//#define _dI8fs(y, z, w, a, b)   ( _F(y, z, w) )* ( outerProduct( _v0(y, a), _w0(y, b) ) + outerProduct( _w0(y, b), _v0(y, a) ) )

#define _I8bar(F, v0, w0)            ( _Jm23(F) * _I8(F, v0, w0) )

// dI1bar = Derivative of \bar{I}_1 with respect to F
//matrix
#define _dI8bar(F, v0, w0)       ( _Jm23(F) * _dI8(F, v0, w0) + _I8(F, v0, w0) * _dJm23(F) )

#define _dI8dF(F, v0, w0, dF)   ( dot( _dI8(F, v0, w0), dF  ) )
#define _dI8bardF(F, v0, w0, dF)   ( dot( _dI8bar(F, v0, w0), dF  ) )

#define _d2I8dF(v0, w0, dF)   ( dF * ( outerProduct( v0, w0 ) + outerProduct( w0, v0 ) ) )

#define _d2I8bardF(F, v0, w0, dF) ( ( _Jm23(F)              * _d2I8dF(v0, w0, dF) ) \
							      + ( _dJm23dF(F, dF)       * _dI8(F, v0, w0)     ) \
							      + ( _dI8dF(F, v0, w0, dF) * _dJm23(F)           ) \
							      + ( _I8(F, v0, w0)        * _d2Jm23dF(F, dF)    ) )


///////////////////////////////////////////////////////////////////////////
// ACTIVE STRAIN TERMS
///////////////////////////////////////////////////////////////////////////

// F_A = active deformation gradient tensor
//matrix
#define _FA( gf, gs, gn, f0, s0, n0 )        ( _I + gf * outerProduct(f0, f0) + gs * outerProduct(s0, s0) + gn * outerProduct(n0, n0) )

//active function
//scalar
#define _gm( g )  ( value(-1.0) * ( g ) / ( ( g ) + 1.0 ) )

//scalar
#define _gti( g )  ( pow( ( 1. + ( g ) ), 0.5 ) - value(1.0) )
//scalar
#define _gtip1( g )  ( pow( ( 1. + ( g ) ), 0.5 ) )


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
//// C_A^-1_{anisotropic} =  inverse of the active right Cauchy-Green tensor
////matrix
#define _FAinvTI( gf, f0, s0, n0 )       ( _I + _gm( gf ) * outerProduct(f0, f0) \
		 	 	 	 	 	 	 	 	 	  + _gti( gf ) * outerProduct(s0, s0) \
		 	 	 	 	 	 	 	 	 	  + _gti( gf ) * outerProduct(n0, n0) )
// C_A^-1_{anisotropic} =  inverse of the active right Cauchy-Green tensor
//matrix
#define _FAinvTI2( gf, f0 )       ( _gtip1( gf ) * _I + ( _gm( gf ) - _gti( gf ) ) * outerProduct(f0, f0) )

//helper macros
#define GET_FAinv(_1,_2,_3,_4,_5,_6,NAME,...) NAME
#define _FAinv(...) GET_FAinv(__VA_ARGS__, _FAinvA, _FAinvO, _FAinvTI, _FAinvTI, _FAinvTI2 ) (__VA_ARGS__)


// dF_E = dF * FAinv
//matrix
#define _FE( F, FAinv )  ( F * ( FAinv ) )
// dF_E = dF * FAinv
//matrix
#define _dFE( FAinv )  ( _dF * ( FAinv ) )


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


/////////////////////////////////////////////////////////////////////////////
//// ISOTROPIC TERMS DEPENDING ON I1
/////////////////////////////////////////////////////////////////////////////
//// I_1 = first invariant of C
//#define _I1E(F, FAinv)           ( dot( _FE( F, FAinv ) , _FE( F, FAinv ) ) )
//// dI1 = Derivative of I_1 with respect of F = 2 F
////matrix
//#define _dI1E(F, FAinv)          ( value(2.0) * _FE( F, FAinv ) )
//
//// \bar{I}_1 = modified invariant of C
//#define _I1Ebar(F, FAinv)        ( _Jm23(F) * _IE1(F, FAinv) )
//// dI1bar = Derivative of \bar{I}_1 with respect to F
////matrix
//#define _dI1Ebar(F, FAinv)       ( _Jm23(F) * _dI1E(F, FAinv) + _I1E(F, FAinv) * _dJm23(F) )
//
//
//// dI1dF = dI1 : dF = Derivative of I_1 with respect to F in direction dF
////scalar
//#define _dI1EdF(F, FAinv)        ( dot( _dI1E(F, FAinv), _dFE(FAinv) ) )
//// dI1bardF = dI1bar : dF = Derivative of \bar{I}_1 with respect to F in direction dF
////scalar
//#define _dI1EbardF(F, FAinv)     ( dot( _dI1Ebar(F, FAinv), _dFE(FAinv) ) )
//
//// d2I1dF = d2I1 : dF = Second derivative of I_1 with respect to F in direction dF
////matrix
//#define _d2I1EdF (FAinv)      ( value(2.0) * ( _dFE( FAinv ) ) )
//
//// d2I1bardF = d2I1bar : dF = Second derivative of I1bar with respect to F in direction dF
////matrix
//#define _d2I1EbardF(F, FAinv)    ( ( _dJm23dF(F)       * _dI1E(F, FAinv) ) \
//						         + ( _Jm23(F)          * _d2I1EdF(FAinv) ) \
//						         + ( _I1E(F, FAinv)    * _d2Jm23dF(F) ) \
//						         + ( _dI1EdF(F, FAinv) * _dJm23(F) ) )




#endif /* EMMECHANICALEXPRESSIONS_HPP_ */
