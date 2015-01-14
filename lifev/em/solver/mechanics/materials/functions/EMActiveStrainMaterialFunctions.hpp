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
 *  @file EMActiveStrainMaterialFunctions
 *  @brief File containing the base class to introduce new material functions
 *
 *  @date 04-2014
 *  @author Simone Rossi <simone.rossi@epfl.ch>
 *
 */

#ifndef EMACTIVESTRAINMATERIALFUNCTIONS_HPP_
#define EMACTIVESTRAINMATERIALFUNCTIONS_HPP_


//namespace MaterialFunctions
namespace LifeV
{
//! EMActiveStrainMaterialFunctions
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

template <class Mesh>
class EMActiveStrainMaterialFunctions : public virtual EMMaterialFunctions<Mesh>
{
public:


    typedef EMData          data_Type;
    typedef typename boost::shared_ptr<data_Type>  dataPtr_Type;



    EMActiveStrainMaterialFunctions() {}
    EMActiveStrainMaterialFunctions (const EMActiveStrainMaterialFunctions& function)
    {
    	M_activeStrainOrthotropicParameter = function.M_activeStrainOrthotropicParameter;
    	M_activeStrainType = function.M_activeStrainType;
    }
    virtual ~EMActiveStrainMaterialFunctions() {}

    virtual void setParameters (data_Type& data)
    {
    	M_activeStrainOrthotropicParameter = data.solidParameter<Real>("EMActiveStrainOrthotropicParameter");
    	selectActiveStrainType( data.solidParameter<std::string>("EMactiveStrainType") );

    }

    void selectActiveStrainType(std::string type)
    {
    	if( "TransverselyIsotropic" == type)
    	{
			M_activeStrainType = TransverselyIsotropic;
    	}
    	else if ("Orthotropic" == type )
    	{
    		M_activeStrainType = Orthotropic;
    	}
    	else
    	{
    		M_activeStrainType = Anisotropic;
    	}
    }

    enum ActivationType { Anisotropic, Orthotropic, TransverselyIsotropic };

    ActivationType M_activeStrainType;
    Real M_activeStrainOrthotropicParameter;
};

} //EMActiveStrainMaterialFunctions

} //LifeV
#endif /* EMActiveStrainMaterialFunctions_HPP_ */
