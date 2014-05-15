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
  @file IonicAlievPanfilov
  @brief Ionic model of Aliev-Panfilov

  Model from:
  Nash, Martyn P., and Alexander V. Panfilov.
  "Electromechanical model of excitable tissue to study reentrant cardiac arrhythmias."
  Progress in biophysics and molecular biology 85.2 (2004): 501-522.

  @date 01-2013
  @author Simone Rossi <simone.rossi@epfl.ch>

  @contributors
  @mantainer Simone Rossi <simone.rossi@epfl.ch>
  @last update 01-2013
 */


#ifndef _EMIONICALIEVPANFILOV_H_
#define _EMIONICALIEVPANFILOV_H_

#include <lifev/electrophysiology/solver/IonicModels/IonicAlievPanfilov.hpp>

namespace LifeV
{
//! IonicModel - This class implements an ionic model.


class EMIonicAlievPanfilov : public virtual IonicAlievPanfilov
{

public:
    //! @name Type definitions
    //@{
    typedef IonicAlievPanfilov  super;
    typedef boost::shared_ptr<VectorEpetra> vectorPtr_Type;
    typedef RegionMesh<LinearTetra> mesh_Type;
    //@}



    //! @name Constructors & Destructor
    //@{

    //! Constructor
    EMIonicAlievPanfilov();

    /*!
     * @param Epetra communicator
     * @param list of parameters in an xml file
     */
    EMIonicAlievPanfilov ( Teuchos::ParameterList& parameterList );

    /*!
     * @param IonicAlievPanfilov object
     */
    EMIonicAlievPanfilov ( const EMIonicAlievPanfilov& model );
    //! Destructor
    virtual ~EMIonicAlievPanfilov() {}

    //@}

    //! @name Overloads
    //@{

    EMIonicAlievPanfilov& operator= ( const EMIonicAlievPanfilov& model );

    //@}

    //! @name Setters and getters
    //@{


    // inline const short int& Size() const { return M_numberOfEquations; }
    //@}

    //! @name Methods
    //@{

    //Compute the rhs on a single node or for the 0D case
    void computeGatingRhs ( const std::vector<Real>& v, std::vector<Real>& rhs);

    void computeRhs ( const std::vector<Real>& v, std::vector<Real>& rhs);

    //Compute the rhs on a mesh/ 3D case
    //    void computeRhs( const std::vector<vectorPtr_Type>& v, std::vector<vectorPtr_Type>& rhs );
    //
    //    void computeRhs( const std::vector<vectorPtr_Type>& v, const VectorEpetra& Iapp, std::vector<vectorPtr_Type>& rhs );
    //

    // compute the rhs with state variable interpolation
    Real computeLocalPotentialRhs ( const std::vector<Real>& v);

    //    void computePotentialRhs(     const std::vector<vectorPtr_Type>& v,
    //                      const VectorEpetra& Iapp,
    //                      std::vector<vectorPtr_Type>& rhs,
    //                      FESpace<mesh_Type, MapEpetra>& uFESpace );

    //! Display information about the model
    void showMe();

    //! Solves the ionic model
    //virtual void solveXbModel( const vector_Type& Calcium,
    //                           const Real timeStep )=0;
    //@}

private:
    //! Model Parameters

    //! Chemical kinetics parameters
    Real M_mu1;
    Real M_mu2;
    Real M_epsilon;
    Real M_k;
    Real M_a;


    //@}

}; // class EMIonicAlievPanfilov


}

#endif
