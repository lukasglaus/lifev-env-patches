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
  @file IonicOharaRudy
  @brief Ionic model O'hara Rudy

  @date 03-2015
  @author Thomas Kummer <kummerth@ethz.ch>

  @contributors
  @mantainer Thomas Kummer <kummerth@ethz.ch>
  @last update 03-2015
 */


#ifndef _IONICOHARARUDY_H_
#define _IONICOHARARUDY_H_

#include <lifev/electrophysiology/solver/IonicModels/ElectroIonicModel.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include "Teuchos_XMLParameterListHelpers.hpp"

namespace LifeV
{
//! IonicModel - This class implements an ionic model.

class IonicOharaRudy : public virtual ElectroIonicModel
{

public:
    //! @name Type definitions
    //@{
    typedef ElectroIonicModel                         super;
    typedef boost::shared_ptr<VectorEpetra>           vectorPtr_Type;
    typedef boost::shared_ptr<VectorElemental>        elvecPtr_Type;
    typedef RegionMesh<LinearTetra>                   mesh_Type;
    //@}



    //! @name Constructors & Destructor
    //@{

    //! Constructor
    IonicOharaRudy();

    /*!
     * @param Epetra communicator
     * @param list of parameters in an xml file
     */
    IonicOharaRudy ( Teuchos::ParameterList& parameterList );

    /*!
     * @param IonicLuoRudyI object
     */
    IonicOharaRudy ( const IonicOharaRudy& model );
    //! Destructor
    virtual ~IonicOharaRudy() {}

    //@}

    //! @name Overloads
    //@{

    IonicOharaRudy& operator= ( const IonicOharaRudy& model );

    //@}

    
    
    //! @name Setters and getters
    //@{

    //parameters getters and setters
    
    //inline const short int& Size() const { return M_numberOfEquations; }
    //@}

    
    
    //! @name Methods
    //@{

    // Compute reversal potentials
    void revpots ( Real nai,
                  Real ki );
    
    // CaMKt
    Real dCaMKt ( Real CaMKt,
                 Real cass );
    
    }
   
    
       //Total current
    //=============
    inline Real Itot (Real V, Real m, Real h, Real j, Real d, Real f, Real X, Real Ca)
    {
        return ( INa (V, m, h, j) + Isi (V, d, f, Ca) + IK (V, X) + IK1 (V) + IKp (V) + Ib (V) );
    }
    
    //Compute the rhs on a single node or for the 0D case
    void computeGatingRhs ( const std::vector<Real>& v, std::vector<Real>& rhs);

    void computeNonGatingRhs ( const std::vector<Real>& v, std::vector<Real>& rhs );

    void computeRhs ( const std::vector<Real>& v, std::vector<Real>& rhs);

    void computeGatingVariablesWithRushLarsen ( std::vector<Real>& v, const Real dt );

    // compute the rhs with state variable interpolation
    Real computeLocalPotentialRhs ( const std::vector<Real>& v );

    //! Display information about the model
    void showMe();

    //! Solves the ionic model
    //virtual void solveXbModel( const vector_Type& Calcium,
    //                           const Real timeStep )=0;
    //@}

private:
    //! Model Parameters
    
    //constants
    Real M_nao;
    Real M_cao;
    Real M_ko;
    
    //buffer paramaters
    Real M_BSRmax;
    Real M_KmBSR;
    Real M_BSLmax;
    Real M_KmBSL;
    Real M_cmdnmax;
    Real M_kmcmdn;
    Real M_trpnmax;
    Real M_kmtrpn;
    Real M_csqnmax;
    Real M_kmcsqn;
    
    //CaMK paramaters
    Real M_aCaMK;
    Real M_bCaMK;
    Real M_CaMKo;
    Real M_KmCaM;
    Real M_KmCaMK;
    
    //physical constants
    Real M_R;
    Real M_T;
    Real M_F;
    
    //cell geometry
    Real M_L;
    Real M_rad;
    Real M_vcell;
    Real M_Ageo;
    Real M_Acap;
    Real M_vmyo;
    Real M_vmito;
    Real M_vsr;
    Real M_vnsr;
    Real M_vjsr;
    Real M_vss;
    
    //introduce varaibles for reversal potentials, currents, fluxes, and CaMK
    Real M_ENa,M_EK,M_EKs;
    Real M_INa,M_INaL,M_Ito,M_ICaL,M_ICaNa,M_ICaK,M_IKr,M_IKs,M_IK1,M_INaCa_i,M_INaCa_ss,M_INaCa,M_INaK,M_IKb,M_INab,M_IpCa,M_ICab,M_Ist;
    Real M_Jrel,M_Jup,M_Jtr,M_Jdiff,M_JdiffNa,M_JdiffK,M_Jleak;
    Real M_CaMKa,M_CaMKb;
    
    //introduce APD, timing, and counting parameters
    int M_APD_flag;
    Real M_APD;
    Real M_t_vdot_max;
    Real M_vrest;
    Real M_vo;
    Real M_dt;
    Real M_t0;
    Real M_t;
    Real M_dto;
    Real M_vdot_old;
    Real M_vdot;
    Real M_vdot_max;
    int M_p;
    int M_n;
    int M_count;
    
    //! Xb states == equivalent to the number of equations
    //short int M_numberOfEquations;

    //@}

}; // class IonicLuoRudyI

inline ElectroIonicModel* createIonicOharaRudy()
{
    return new IonicOharaRudy();
}

namespace
{
    static bool register_IonicOharaRudy = ElectroIonicModel::IonicModelFactory::instance().registerProduct ("OharaRudy", &createIonicOharaRudy );
}


}

#endif
