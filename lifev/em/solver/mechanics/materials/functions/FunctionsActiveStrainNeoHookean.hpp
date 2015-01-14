/*
 * EMMaterialFunctions.hpp
 *
 *  Created on: 28/apr/2014
 *      Author: srossi
 */

#ifndef FUNCTIONActiveStrainNeoHookean_HPP_
#define FUNCTIONActiveStrainNeoHookean_HPP_

#include <lifev/em/solver/mechanics/EMElasticityFunctions.hpp>

//#include <lifev/em/solver/mechanics/EMETAAssembler.hpp>
#include <lifev/em/solver/mechanics/materials/functions/FunctionsNeoHookean.hpp>
#include <lifev/em/solver/mechanics/materials/functions/EMActiveStrainMaterialFunctions.hpp>
//using namespace LifeV;

//namespace MaterialFunctions
namespace LifeV
{

namespace MaterialFunctions
{



typedef VectorEpetra           vector_Type;
typedef boost::shared_ptr<vector_Type>         vectorPtr_Type;

typedef MatrixEpetra<Real>           matrix_Type;
typedef boost::shared_ptr<matrix_Type>         matrixPtr_Type;


////////////////////////////////////////////////////////////////////////
//  ACTIVE STRAIN NEO HOOKEAN FUNCTIONS
////////////////////////////////////////////////////////////////////////
template <class Mesh>
class ActiveStrainNeoHookean : public virtual NeoHookean<Mesh>,
                               public virtual EMActiveStrainMaterialFunctions<Mesh>
{
public:

	typedef NeoHookean<Mesh> super;
	typedef typename super::data_Type data_Type;
	typedef EMActiveStrainMaterialFunctions<Mesh> activeStrain;
	ActiveStrainNeoHookean(Real mu = 4960) : super(mu) {}


    void computeJacobian ( const vector_Type& disp,
                                  boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > >  dispETFESpace,
                                  const vector_Type& fibers,
                                  const vector_Type& sheets,
                                  const vectorPtr_Type& fiberActivation,
                                  const vectorPtr_Type& sheetActivation,
                                  const vectorPtr_Type& normalActivation,
                                  boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 1 > >  activationETFESpace,
                                  matrixPtr_Type           jacobianPtr)
    {
        EMAssembler::computeActiveStrainI1JacobianTerms (disp,
														dispETFESpace,
														fibers,
														sheets,
														fiberActivation,
														sheetActivation,
														normalActivation,
														activationETFESpace,
														jacobianPtr,
														this->M_W1 );
    }

    virtual void computeResidual ( const vector_Type& disp,
                                          boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > >  dispETFESpace,
                                          const vector_Type& fibers,
                                          const vector_Type& sheets,
                                          const vectorPtr_Type& fiberActivation,
                                          const vectorPtr_Type& sheetActivation,
                                          const vectorPtr_Type& normalActivation,
                                          boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 1 > >  activationETFESpace,
                                          vectorPtr_Type           residualVectorPtr)
    {
    	EMAssembler::computeActiveStrainI1ResidualTerms(  disp,
    			 	 	 	 	 	 	 	 	 	 	 	 	 dispETFESpace,
    			 	 	 	 	 	 	 	 	 	 	 	 	 fibers,
    			 	 	 	 	 	 	 	 	 	 	 	 	 sheets,
    			 	 	 	 	 	 	 	 	 	 	 	 	 fiberActivation,
    			 	 	 	 	 	 	 	 	 	 	 	 	 sheetActivation,
    			 	 	 	 	 	 	 	 	 	 	 	 	 normalActivation,
    			 	 	 	 	 	 	 	 	 	 	 	 	 activationETFESpace,
    			 	 	 	 	 	 	 	 	 	 	 	 	 residualVectorPtr,
    			 	 	 	 	 	 	 	 	 	 	 	 	 this->M_W1 );
    }

    virtual void setParameters (data_Type& data)
    {
    	activeStrain::setParameters(data);
    	super::setParameters(data);
    }

private:

};


} //EMMaterialFunctions

} //LifeV
#endif /* EMMATERIALFUNCTIONS_HPP_ */
