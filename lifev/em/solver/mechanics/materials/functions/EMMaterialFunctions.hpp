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
 *  @file EMMaterialFunctions
 *  @brief File containing the base class to introduce new material functions
 *
 *  @date 04-2014
 *  @author Simone Rossi <simone.rossi@epfl.ch>
 *
 */

#ifndef EMMATERIALFUNCTIONS_HPP_
#define EMMATERIALFUNCTIONS_HPP_

#include <boost/smart_ptr/enable_shared_from_this.hpp>
#include <lifev/em/solver/EMData.hpp>

#include <lifev/em/solver/mechanics/EMElasticityFunctions.hpp>

#include <lifev/em/solver/mechanics/EMETAJacobianAssembler.hpp>
#include <lifev/em/solver/mechanics/EMETAResidualAssembler.hpp>
#include <lifev/em/solver/mechanics/EMETAActiveStrainResidualAssembler.hpp>
#include <lifev/em/solver/mechanics/EMETAJacobianIsotropicAssembler.hpp>
#include <lifev/em/solver/mechanics/EMETAResidualIsotropicAssembler.hpp>
#include <lifev/em/solver/mechanics/EMETAActiveStrainJacobianAssembler.hpp>

//namespace MaterialFunctions
namespace LifeV
{

class EMData;
//! EMMaterialFunctions
/*!
 *  @author Simone Rossi
 *
 *  This class allows to define the functions corresponding to a given constitutive law.
 *  To introduce a new constitutive law  you should defined by the energy function
 *  \f[ \mathcal{W} = \mathcal{W}(\mathbf{F})\f]
 *  you need to specify the derivatives of the energy.
 *
 *  Without loss of generality we can assume that the energy may be written as
 *  \f[ \mathcal{W} = \sum_j\mathcal{W}_j(\mathbf{F})\f]
 *
 *  For example consider the case in which
 *  \f[ \mathcal{W} = \mathcal{W}_I(\mathcal{I}_1) + \mathcal{W}_{vol}(J)\f]
 *
 *  Then the stress tensor is computed as
 *  \f[ \mathbf{P} = W_I \dfrac{\partial \mathcal{I}_1}{\partial \mathbf{F}}} + W_{vol} \dfrac{\partial J}{\partial \mathbf{F}} \f]
 *
 *  Therefore, if we would like to assemble, for example, the residual of different constitutive
 *  laws of this form we only need to change the functions \f[W_I\f] and \f[W_{vol}\f].
 *
 *  You can derive this class in order to specify \f[W_I\f] and \f[W_{vol}\f].
 *  In order to deal with particular cases, this class needs the implementation
 *  of the computation of the residual and of the jacobian.
 *  Since the implementation is not function specific the general implementation is
 *  inserted in the EMETAAssembler. So material specific-function (you will implement)
 *  will call to the functions in the EMETAAssembler.
 *  If your function needs to assemble some terms which are not  included
 *  in the EMETAAssembler, you should add your own assembly  in that file.
 *
 *
 */


namespace MaterialFunctions
{

typedef VectorEpetra           vector_Type;
typedef boost::shared_ptr<vector_Type>         vectorPtr_Type;

typedef MatrixEpetra<Real>           matrix_Type;
typedef boost::shared_ptr<matrix_Type>         matrixPtr_Type;


template <class Mesh>
class EMMaterialFunctions : public boost::enable_shared_from_this<EMMaterialFunctions<Mesh> >
{
public:


    typedef  boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > > ETFESpacePtr_Type;
    typedef boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 1 > > scalarETFESpacePtr_Type;
    typedef Real return_Type;
    typedef MatrixSmall<3, 3> matrix_return_Type;
    typedef boost::enable_shared_from_this<EMMaterialFunctions< Mesh > > boostShared_Type;

    typedef EMData          data_Type;
    typedef typename boost::shared_ptr<data_Type>  dataPtr_Type;

    virtual return_Type operator() (const MatrixSmall<3, 3>& F)
    {
        return 0.0;
    }

    virtual return_Type operator() (const VectorSmall<3>& f)
    {
        return 0.0;
    }

    virtual return_Type operator() (const MatrixSmall<3, 3>& F, const VectorSmall<3>& f, const VectorSmall<3>& s)
    {
        return 0.0;
    }

    virtual return_Type operator() (const VectorSmall<3>& f, const VectorSmall<3>& s)
    {
        return 0.0;
    }

    virtual return_Type operator() (const Real& H)
    {
        return 0.0;
    }

    virtual matrix_return_Type operator() (const MatrixSmall<3, 3>& F, bool isReturnTypeMatrix)
    {
        return MatrixSmall<3,3>();
    }

    EMMaterialFunctions() {}
    EMMaterialFunctions (const EMMaterialFunctions&) {}
    virtual ~EMMaterialFunctions() {}

    inline virtual boost::shared_ptr<EMMaterialFunctions<Mesh> > getMe()
    {
        //      shared_from_this();
        return boostShared_Type::shared_from_this();
    }

    inline virtual void computeJacobian ( const vector_Type& disp,
                                          ETFESpacePtr_Type  dispETFESpace,
                                          const vector_Type& fibers,
                                          const vector_Type& sheets,
                                          matrixPtr_Type           jacobianPtr) {}

    inline virtual void computeJacobian ( const vector_Type& disp,
                                          ETFESpacePtr_Type  dispETFESpace,
                                          const vector_Type& fibers,
                                          const vector_Type& sheets,
                                          const vectorPtr_Type& fiberActivation,
                                          const vectorPtr_Type& sheetActivation,
                                          const vectorPtr_Type& normalActivation,
                                          scalarETFESpacePtr_Type activationETFESpace,
                                          matrixPtr_Type     jacobianPtr) {}



    inline virtual void computeResidual ( const vector_Type& disp,
                                          ETFESpacePtr_Type dispETFESpace,
                                          const vector_Type& fibers,
                                          const vector_Type& sheets,
                                          vectorPtr_Type           residualVectorPtr) {}

    inline virtual void computeResidual ( const vector_Type& disp,
                                          ETFESpacePtr_Type dispETFESpace,
                                          const vector_Type& fibers,
                                          const vector_Type& sheets,
                                          const vectorPtr_Type& fiberActivation,
                                          const vectorPtr_Type& sheetActivation,
                                          const vectorPtr_Type& normalActivation,
                                          scalarETFESpacePtr_Type  activationETFESpace,
                                          vectorPtr_Type           residualVectorPtr) {}


    virtual void showMe() { std::cout << "\nYuo should implement the showMe() method for your constitutive law!\n"; }

    virtual void setParameters (data_Type& data) = 0;

};

} //EMMaterialFunctions

} //LifeV
#endif /* EMMATERIALFUNCTIONS_HPP_ */
