/*
 * EMMaterialFunctions.hpp
 *
 *  Created on: 28/apr/2014
 *      Author: srossi
 */

#ifndef FUNCTIONSNEOHOOKEAN_HPP_
#define FUNCTIONSNEOHOOKEAN_HPP_

#include <lifev/em/solver/mechanics/EMElasticityFunctions.hpp>

//#include <lifev/em/solver/mechanics/EMETAAssembler.hpp>
#include <lifev/em/solver/mechanics/materials/functions/EMMaterialFunctions.hpp>

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
//  NEO HOOKEAN FUNCTIONS
////////////////////////////////////////////////////////////////////////
template <class Mesh>
class NeoHookean : public virtual EMMaterialFunctions<Mesh>
{
public:
    typedef EMData          data_Type;

//    virtual return_Type operator() (const MatrixSmall<3, 3>& F)
//    {
//        return 0.5 * M_mu;
//    }
//
//    //    NeoHookean() : M_mu(4960) {} // 0.496 KPa
      NeoHookean (Real mu = 4960) : M_W1 (new W1(mu)) { } // 0.496 KPa
//    NeoHookean (const NeoHookean& neoHookean)
//    {
//        M_mu = neoHookean.M_mu;
//    }
//    virtual ~NeoHookean() {}

    class W1
    {
    public:
        typedef typename MaterialFunctions::EMMaterialFunctions<Mesh>::return_Type return_Type;
        virtual return_Type operator() (const MatrixSmall<3, 3>& F)
        {
            return 0.5 * mu;
        }

        //    NeoHookean() : M_mu(4960) {} // 0.496 KPa
        W1 (Real mu = 4960) : mu (mu) {} // 0.496 KPa
        W1 (const W1& neoHookean)
        {
            mu = neoHookean.mu;
        }
        virtual ~W1() {}

        void setMu( Real mu )
        {
        	this->mu = mu;
        }

        Real mu;

    };

    virtual void computeJacobian ( const vector_Type& disp,
                                          boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > >  dispETFESpace,
                                          const vector_Type& fibers,
                                          const vector_Type& sheets,
                                          matrixPtr_Type           jacobianPtr )
    {
        EMAssembler::computeI1JacobianTerms ( disp, dispETFESpace, jacobianPtr, M_W1 );
    }

    virtual void computeResidual ( const vector_Type& disp,
                                          boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > >  dispETFESpace,
                                          const vector_Type& fibers,
                                          const vector_Type& sheets,
                                          vectorPtr_Type residualVectorPtr)
    {
        EMAssembler::computeI1ResidualTerms ( disp, dispETFESpace, residualVectorPtr, M_W1 );
    }

    void setParameters (data_Type& data)
    {
    	M_W1->setMu( data.parameter("mu") );
    }

    void showMe()
    {
    	std::cout << "Shear Modulus = " << M_W1->mu << "\n";
    }
private:
    boost::shared_ptr<NeoHookean::W1> M_W1;
};


} //EMMaterialFunctions

} //LifeV
#endif /* EMMATERIALFUNCTIONS_HPP_ */
