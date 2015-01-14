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
  @brief Base class for ionic models

  This class is a  abstract base class for all ionic models.

  Assume that the ionic model is written in the form

  \f[ \dfrac{\partial \mathbf{v} }{\partial t} = f(\mathbf{v), t\f]

  If you wish to implement a new ionic model you should create a new
  class which inherits from this class and implement the abstract methods:

  Given \f[ \mathbf{v} \f] compute \f[ f(\mathbf{v}, t) \f]
  void computeRhs ( const std::vector<Real>& v, std::vector<Real>& rhs);
  In principle this method should call the following two methods:

  Given \f[ \mathbf{v} \f] compute \f[ f(\mathbf{v}, t) \f]  only for the first
  equation, that is only for the  transmembrane potential
  Real computeLocalPotentialRhs ( const std::vector<Real>& v );

  Given \f[ \mathbf{v} \f] compute \f[ f(\mathbf{v}, t) \f] for all the variables
  excpet the first one, that is all variables except the transmembrane potential
  void computeGatingRhs ( const std::vector<Real>& v, std::vector<Real>& rhs );

  Print out informations about the ionic model
  void showMe();

  Some ionic models may be solved with the Rush Larsen scheme.
  The implementation of this method is not mandatory.
  For biophysically detailed models, it may be convenient,
  to have it.
  void computeGatingVariablesWithRushLarsen ( std::vector<Real>& v, const Real dt ) {}

  @date 01-2013
  @author Simone Rossi <simone.rossi@epfl.ch>

  @contributors
  @mantainer Simone Rossi <simone.rossi@epfl.ch>
  @last update 02-2012
 */


#ifndef _EMIonicModel_H_
#define _EMIonicModel_H_


#include <lifev/core/array/MapEpetra.hpp>
#include <lifev/core/array/VectorEpetra.hpp>

//ET include for assemblings
#include <lifev/eta/fem/ETFESpace.hpp>
#include <lifev/eta/expression/Integrate.hpp>

#include <lifev/electrophysiology/solver/IonicModels/ElectroIonicModel.hpp>


namespace LifeV
{

class EMIonicModel : public virtual ElectroIonicModel
{

public:
    //! @name Type definitions
    //@{

    typedef VectorEpetra                            vector_Type;
    typedef boost::shared_ptr<VectorEpetra>         vectorPtr_Type;
    typedef boost::shared_ptr<VectorElemental>  elvecPtr_Type;
    typedef RegionMesh<LinearTetra>                 mesh_Type;
    typedef MatrixEpetra<Real> matrix_Type;
    typedef boost::shared_ptr<matrix_Type> matrixPtr_Type;
    typedef FESpace<mesh_Type, MapEpetra> feSpace_Type;
    typedef boost::shared_ptr<feSpace_Type> feSpacePtr_Type;
    typedef boost::function <
    Real (const Real& t, const Real& x, const Real& y, const Real& z,
          const ID& i) > function_Type;
    //@}

    //! @name Constructors & Destructor
    //@{

    //! Empty Constructor
    /*!
     */

    //@}
    template < class Mesh >
    void addSAC( vector_Type& V,
                 boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 1 > > potentialETFESpace,
                 vector_Type& disp,
                 boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > > dispETFESpace ){}


    void setupParameters(EMData& data){}

protected:


};


}

#endif
