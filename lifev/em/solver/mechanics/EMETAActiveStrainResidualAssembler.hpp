/*
 * EMETAAssmebler.hpp
 *
 *  Created on: 29/apr/2014
 *      Author: srossi
 */

#ifndef EMETAACTIVESTRAINRESIDUALASSMEBLER_HPP_
#define EMETAACTIVESTRAINRESIDUALASSMEBLER_HPP_

#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/array/MatrixEpetra.hpp>

//ET include for assemblings
#include <lifev/eta/fem/ETFESpace.hpp>
#include <lifev/eta/expression/Integrate.hpp>

#include <lifev/em/solver/mechanics/EMMechanicalExpressions.hpp>

#include <lifev/em/solver/EMETAFunctors.hpp>

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
computeActiveStrainI1ResidualTerms (  const vector_Type& disp,
													boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > >  dispETFESpace,
													const vector_Type& fibers,
													const vector_Type& sheets,
													const vectorPtr_Type& gammaf,
													const vectorPtr_Type& gammas,
													const vectorPtr_Type& gamman,
													boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 1 > >  activationETFESpace,
													vectorPtr_Type           residualVectorPtr,
													FunctorPtr               W1,
													Real orthotropicParameter = -666.)
{
    //
	if(disp.comm().MyPID() == 0)
    std::cout << "EMETA - Computing Isotropic Active Strain residual terms: ";

    using namespace ExpressionAssembly;

    auto F = _F (dispETFESpace, disp, 0);
    auto I = _I;

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
		auto P = eval (W1, FE ) * _dI1bar (dispETFESpace, disp, 0) * CAinv;

		integrate ( elements ( dispETFESpace->mesh() ) ,
					quadRuleTetra4pt,
					dispETFESpace,
					dot ( P, grad (phi_i) )
				  ) >> residualVectorPtr;
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
    		auto P = eval (W1, FE ) * _dI1bar (dispETFESpace, disp, 0) * CAinv;

    		integrate ( elements ( dispETFESpace->mesh() ) ,
    					quadRuleTetra4pt,
    					dispETFESpace,
    					dot ( P, grad (phi_i) )
    				  ) >> residualVectorPtr;
    	}
    	else
    	{
        	if(disp.comm().MyPID() == 0)
            std::cout << " Transversely isotropic case ... \n";

    		auto FAinv = _FAinv(gf, f0, s0, n0);
    		auto CAinv = _CAinv(gf, f0, s0, n0);

    		auto FE =  F * FAinv;
    		auto P = eval (W1, FE ) * _dI1bar (dispETFESpace, disp, 0) * CAinv;

    		integrate ( elements ( dispETFESpace->mesh() ) ,
    					quadRuleTetra4pt,
    					dispETFESpace,
    					dot ( P, grad (phi_i) )
    				  ) >> residualVectorPtr;
    	}
    }

}





}//EMAssembler

}//LifeV

#endif /* EMETAASSMEBLER_HPP_ */
