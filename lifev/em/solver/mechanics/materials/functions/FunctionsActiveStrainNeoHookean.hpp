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
    	std::cout << "Norm of the activation = " <<  fiberActivation->norm2() << std::endl;

    	using namespace ExpressionAssembly;

        MatrixSmall<3,3> Id;
        Id (0, 0) = 1.; Id (0, 1) = 0., Id (0, 2) = 0.;
        Id (1, 0) = 0.; Id (1, 1) = 1., Id (1, 2) = 0.;
        Id (2, 0) = 0.; Id (2, 1) = 0., Id (2, 2) = 1.;

        auto I = value(Id);
        auto F = I + grad(dispETFESpace, disp, 0);
        auto dF = grad( phi_j );
        auto gf = value(activationETFESpace, *fiberActivation);
        // auto gf = value(0.);
        // auto gs = value(activationETFESpace, *sheetActivation);
        // auto gn = value(activationETFESpace, *normalActivation);
    	auto f0 = value(dispETFESpace, fibers);
        auto gtip1 = pow( ( 1. + ( gf ) ), 0.5 );
        auto ggf = value(1.) - gtip1 - gf / (gf +value(1.));
        // auto _gm = value(-1.0) * ( gf ) / ( ( gf ) + 1.0 );
        // auto _gti = pow( ( 1. + ( gf ) ), 0.5 ) - value(1.0);
    	// auto FAinv = ( _gtip1 * I + ( _gm - _gti ) * outerProduct(f0, f0) );
        // auto FA = I + (gf)*outerProduct(f0, f0);

        auto gn = pow( 1.+gf, -0.5 )+value(-1.);

//		auto FAinv = gtip1* I + ggf * outerProduct(f0, f0);

//        auto FAinv = ( 1 - gn / (gn + value(1.) ) ) * I + ( gn / (gn + value(1.) ) - gf/(value(1.)+gf) )* outerProduct(f0, f0);
		auto FAinv = _FAinv(gf, f0);

//        auto FAinv = I - gf/(value(1.)+gf)*outerProduct(f0, f0);
        // auto k  = value(1.e+5);
        // auto dW1 = eval(this->M_W1, F);

        auto J = det(F);
        auto Jm23 = pow( J, -2./3.);
        auto Fe = F * FAinv;
        auto Je = det(Fe);
        auto Jem23 = pow( Je, -2./3.);
        auto FmT = minusT(F);
        auto FemT = minusT(Fe);
        auto dJdF =  J * FmT;
        auto dJedF = Je * FemT;
        // // VOLUMETRIC PART
        // auto dWvoldJ = value(0.5) * k * ( J - value(1.)/J );
        // auto dWvoldF = dWvoldJ * dJdF;
        // auto d2Wvold2J = value(0.5) * k * ( value(1.) + value(1.)/(J*J) );
        // auto d2Wvold2F = d2Wvold2J * dot( dJdF, dF ) * dJdF + dWvoldJ * ( dJdF * dot(FmT, dF)  - J * FmT * transpose(dF) * FmT );
        auto dW1dI1e = eval(this->M_W1, F);
        // auto mu = value(10.);
        // auto dW1dI1e = mu*value(0.5);
        // ISOTROPIC PART
        // auto dI1edI1 =   value (1.0) - gn * ( gn + value (2.0) ) * pow (gn + value (1.0), -2.0);
        // auto dI1edI4f =   gn * ( gn + value (2.0) ) * pow (gn + value (1.0), -2.0) - gf * ( gf + value (2.0) ) * pow (gf + value (1.0), -2.0);
        // auto dI1edI4s =   gn * ( gn + value (2.0) ) * pow (gn + value (1.0), -2.0) - gs * ( gs + value (2.0) ) * pow (gs + value (1.0), -2.0);
        // FIRST ELASTIC INVARIANT DEPENDS ON FIRST AND FOURTH INVARIANTS
        // FIRST BAR-INVARIANT DERIVATIVE
        auto I1Ce = dot( Fe, Fe );
        auto dFe = dF * FAinv;
        auto dI1Ce = value(2.0) * Fe;
        auto dI1CedF = value(2.0) * dot( Fe, dF );
        auto dI1CedFe = value(2.0) * dot( Fe, dFe );
        auto d2I1CedF = value(2.0) * dF;
        auto d2I1CedFe = value(2.0) * dFe;
        auto dJem23 = value(-2.0/3.0) * Jem23 * FemT;
        auto dJem23dF = value(-2.0/3.0) * Jem23 * dot( FemT, dF );
        auto dJem23dFe = value(-2.0/3.0) * Jem23 * dot( FemT, dFe );
        auto dFemTdF = value(-1.) * FemT * transpose( dF ) * FemT;
        auto dFemTdFe = value(-1.) * FemT * transpose( dFe ) * FemT;
        auto d2Jem23dF = value(-2./3.) * dJem23dF * FemT + value(-2./3.) * Jem23 * dFemTdF;
        auto d2Jem23dFe = value(-2./3.) * dJem23dF * FemT + value(-2./3.) * Jem23 * dFemTdFe;
        auto d2I1CebardFe =
            dJem23dFe * dI1Ce // value(-4./3.) * Jem23 * dot( FemT, dF ) * Fe
            + Jem23 * d2I1CedFe // + value(2.) * Jem23 * dF
            + I1Ce * d2Jem23dFe // + value(-2./3.) * I1Ce * dJem23dF * FemT // + value(4./9.) * Jem23 * dot( FemT, dF ) * I1Ce * FemT + value(-2./3.) * I1Ce * Jem23 * dFemTdF // + value(2./3.) * I1Ce * Jem23 * FemT * transpose( dF ) * FemT
            + dI1CedFe * dJem23; // value(-2./3.) * Jem23 * FemT * dot( Fe, dF);
        auto dP = dW1dI1e * d2I1CebardFe * FAinv; // + d2Wvold2F;



        // auto F = _F(dispETFESpace, disp, 0);
    	// auto gf = value(activationETFESpace, *fiberActivation);
    	// auto f0 = value(dispETFESpace, fibers);
    	// auto FAinv = _FAinv(gf, f0);
    	// auto W1 = eval(this->M_W1, F);
    	// auto FE = F * FAinv;
	// auto dP = W1 * _d2I1bardF (FE) * FAinv;
        
//	integrate ( elements ( dispETFESpace->mesh() ) ,
//		    quadRuleTetra4pt,
//		    dispETFESpace,
//		    dispETFESpace,
//		    dot ( dP , grad (phi_i) )
//		    ) >> jacobianPtr;
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

    	std::cout << "Norm of the activation = " <<  fiberActivation->norm2() << std::endl;
//    	EMAssembler::computeActiveStrainI1ResidualTerms(  disp,
//    			 	 	 	 	 	 	 	 	 	 	 	 	 dispETFESpace,
//    			 	 	 	 	 	 	 	 	 	 	 	 	 fibers,
//    			 	 	 	 	 	 	 	 	 	 	 	 	 sheets,
//    			 	 	 	 	 	 	 	 	 	 	 	 	 fiberActivation,
//    			 	 	 	 	 	 	 	 	 	 	 	 	 sheetActivation,
//    			 	 	 	 	 	 	 	 	 	 	 	 	 normalActivation,
//    			 	 	 	 	 	 	 	 	 	 	 	 	 activationETFESpace,
//    			 	 	 	 	 	 	 	 	 	 	 	 	 residualVectorPtr,
//    			 	 	 	 	 	 	 	 	 	 	 	 	 this->M_W1 );
    	using namespace ExpressionAssembly;

        MatrixSmall<3,3> Id;
        Id (0, 0) = 1.; Id (0, 1) = 0., Id (0, 2) = 0.;
        Id (1, 0) = 0.; Id (1, 1) = 1., Id (1, 2) = 0.;
        Id (2, 0) = 0.; Id (2, 1) = 0., Id (2, 2) = 1.;

        auto I = value(Id);
        // auto k  = value(1.e+5);
        auto F = I + grad(dispETFESpace, disp, 0);
        auto gf = value(activationETFESpace, *fiberActivation);
        // auto gf = value(0.);
    	auto f0 = value(dispETFESpace, fibers);
        auto _gtip1 = pow( ( 1. + ( gf ) ), 0.5 );
        auto _gm = value(-1.0) * ( gf ) / ( ( gf ) + 1.0 );
        auto _gti = pow( ( 1. + ( gf ) ), 0.5 ) - value(1.0);
    	// auto FAinv = ( _gtip1 * I + ( _gm - _gti ) * outerProduct(f0, f0) );
        auto gtip1 = pow( ( 1. + ( gf ) ), 0.5 );
		auto ggf = value(1.) - gtip1 - gf / (gf +value(1.));

        auto gn = pow( 1.+gf, -0.5 )+value(-1.);


//        auto FAinv = ( 1 - gn / (gn + value(1.) ) ) * I + ( gn / (gn + value(1.) ) - gf/(value(1.)+gf) )* outerProduct(f0, f0);
		auto FAinv = gtip1* I + ggf * outerProduct(f0, f0);

//        auto FAinv = I - gf/(value(1.)+gf)*outerProduct(f0, f0);
        // auto FAinv = I;
        auto dW1dI1e = eval(this->M_W1, F);
        // auto mu = value(10.);
        // auto dW1dI1e = value(0.5) * mu;
        auto J = det(F);
        auto Fe = F * FAinv;
        auto Je = det(Fe);
        auto FmT = minusT(F);
        auto FemT = minusT(Fe);
        // auto dWvol = value(0.5) * k * (J-value(1.)/J)* J * FmT;
        auto dI1e = value(2.0) * Fe;
        auto I1e  = dot( Fe, Fe );
        auto Jem23 = pow( Je, -2./3.);
        auto dJem23 = value(-2.0/3.0) * Jem23 * FemT;
        auto P = dW1dI1e * ( Jem23 * dI1e + I1e * dJem23 ) * FAinv; // + dWvol;
        // auto P = dW1*FAinv*( _Jm23(F) * _dI1(F) + _I1(F) * _dJm23(F) ) + dWvol;
        // auto lnJ = eval( log, J);
        // auto k = value(10.);
        // auto Wvol = k*value(0.25)*( pow(J,2) - value(1.) -lnJ );

        // // auto F = _F(dispETFESpace, disp, 0);
    	// // auto gf = value(activationETFESpace, *fiberActivation);
    	// // auto f0 = value(dispETFESpace, fibers);
    	// // auto FAinv = _FAinv(gf, f0);
    	// // auto W1 = eval(this->M_W1, F);
    	// // auto FE = F * FAinv;
		// // auto P = W1 * _dI1bar (FE) * FAinv;



//	integrate ( elements ( dispETFESpace->mesh() ) ,
//		    quadRuleTetra4pt,
//		    dispETFESpace,
//		    dot ( P, grad (phi_i) )
//		    ) >> residualVectorPtr;

    EMAssembler::computeActiveStrainI1ResidualTerms( disp,
													 dispETFESpace,
													 fibers,
													 sheets,
	 												fiberActivation,
	  												sheetActivation,
		 											normalActivation,
		 											activationETFESpace,
		 											residualVectorPtr,
													this->M_W1);
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
