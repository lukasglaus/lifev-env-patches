/*
 * EMETAAssmebler.hpp
 *
 *  Created on: 29/apr/2014
 *      Author: srossi
 */

#ifndef EMETAACTIVESTRAINJACOBIANASSMEBLER_HPP_
#define EMETAACTIVESTRAINJACOBIANASSMEBLER_HPP_


#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/array/MatrixEpetra.hpp>

//ET include for assemblings
#include <lifev/eta/fem/ETFESpace.hpp>
#include <lifev/eta/expression/Integrate.hpp>

#include <lifev/em/solver/mechanics/EMMechanicalExpressions.hpp>
//#include <lifev/em/solver/mechanics/materials/EMMaterialFunctions.hpp>

#include <lifev/em/util/EMUtility.hpp>
#include <lifev/em/solver/EMETAFunctors.hpp>

//#include <boost/typeof/typeof.hpp>

namespace LifeV
{

typedef VectorEpetra           vector_Type;
typedef boost::shared_ptr<vector_Type>         vectorPtr_Type;

typedef MatrixEpetra<Real>           matrix_Type;
typedef boost::shared_ptr<matrix_Type>         matrixPtr_Type;


namespace EMAssembler
{

template <typename Mesh, typename FunctorPtr >
void
computeActiveStrainI1JacobianTerms (    const vector_Type& disp,
                                        boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > > dispETFESpace,
                                        const vector_Type& fibers,
                                        const vector_Type& sheets,
										const vectorPtr_Type& gammaf,
										const vectorPtr_Type& gammas,
										const vectorPtr_Type& gamman,
										boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 1 > >  activationETFESpace,
										matrixPtr_Type           jacobianPtr,
                                        FunctorPtr               W1,
                                        Real orthotropicParameter = -666.)
{
    //
	if(disp.comm().MyPID() == 0)
    std::cout << "EMETA - Computing Isotropic Active Strain jacobian terms: ";

    using namespace ExpressionAssembly;

    auto F = _F (dispETFESpace, disp, 0);
    auto I = value (EMUtility::identity() );

    auto f_0 = _v0 (dispETFESpace, fibers);
    auto s_0 = _v0 (dispETFESpace, sheets);

    boost::shared_ptr<orthonormalizeFibers> normalize0 (new orthonormalizeFibers);
    auto f0 = eval (normalize0, f_0);

    boost::shared_ptr<orthonormalizeFibers> normalize1 (new orthonormalizeFibers (1) );
    auto s0 = eval (normalize1, f0, s_0);

    boost::shared_ptr<CrossProduct> wedge (new CrossProduct);
    auto n0 = eval ( wedge, f0, s0);


    if(gammas && gamman)
    {
    	if(disp.comm().MyPID() == 0)
        std::cout << " Anisotropic case ... \n";

		auto gf = value (activationETFESpace, *gammaf);
		auto gs = value (activationETFESpace, *gammas);
		auto gn = value (activationETFESpace, *gamman);


		auto FAinv = _FAinv(gf, gs, gn, f0, s0, n0);
		auto CAinv = _CAinv(gf, gs, gn, f0, s0, n0);

		auto FE =  F * FAinv;

		auto dP = eval (W1, FE ) * (_d2I1bardF (dispETFESpace, disp, 0) ) * CAinv;

	    integrate ( elements ( dispETFESpace->mesh() ) ,
	                quadRuleTetra4pt,
	                dispETFESpace,
	                dispETFESpace,
	                dot ( dP , grad (phi_i) )
	              ) >> jacobianPtr;
    }
    else
    {
		auto gf = value (activationETFESpace, *gammaf);

    	if(orthotropicParameter < 0 )
    	{
        	if(disp.comm().MyPID() == 0)
            std::cout << " Orthotropic case ... \n";

    		auto k = value(orthotropicParameter);
    		auto FAinv = _FAinv(gf, k, f0, s0, n0);
    		auto CAinv = _CAinv(gf, k, f0, s0, n0);

    		auto FE =  F * FAinv;
    		auto dP = eval (W1, FE ) * (_d2I1bardF (dispETFESpace, disp, 0) ) * CAinv;

    	    integrate ( elements ( dispETFESpace->mesh() ) ,
    	                quadRuleTetra4pt,
    	                dispETFESpace,
    	                dispETFESpace,
    	                dot ( dP , grad (phi_i) )
    	              ) >> jacobianPtr;
    	}
    	else
    	{
        	if(disp.comm().MyPID() == 0)
            std::cout << " Transversely isotropic case ... \n";

    		auto FAinv = _FAinv(gf, f0, s0, n0);
    		auto CAinv = _CAinv(gf, f0, s0, n0);

    		auto FE =  F * FAinv;
    		auto dP = eval (W1, FE ) * (_d2I1bardF (dispETFESpace, disp, 0) ) * CAinv;

    	    integrate ( elements ( dispETFESpace->mesh() ) ,
    	                quadRuleTetra4pt,
    	                dispETFESpace,
    	                dispETFESpace,
    	                dot ( dP , grad (phi_i) )
    	              ) >> jacobianPtr;
    	}
    }

}




}//EMAssembler

}//LifeV

#endif /* EMETAASSMEBLER_HPP_ */
