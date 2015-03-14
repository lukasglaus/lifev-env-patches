/*
 * EMMaterialFunctions.hpp
 *
 *  Created on: 15/mar/2015
 *      Author: tkummer
 */

#ifndef FUNCTIONSMOONEYRIVLIN_HPP_
#define FUNCTIONSMOONEYRIVLIN_HPP_

#include <lifev/em/solver/mechanics/EMElasticityFunctions.hpp>

//#include <lifev/em/solver/mechanics/EMETAAssembler.hpp>
#include <lifev/em/solver/mechanics/materials/functions/EMMaterialFunctions.hpp>

//using namespace LifeV;

//namespace MaterialFunctions
namespace LifeV
{

namespace MaterialFunctions
{



typedef VectorEpetra                            vector_Type;
typedef boost::shared_ptr<vector_Type>          vectorPtr_Type;

typedef MatrixEpetra<Real>                      matrix_Type;
typedef boost::shared_ptr<matrix_Type>          matrixPtr_Type;



////////////////////////////////////////////////////////////////////////
//  MOONEY RIVLIN FUNCTIONS
////////////////////////////////////////////////////////////////////////
template <class Mesh>
class MooneyRivlin : public virtual EMMaterialFunctions<Mesh>
{
public:
    
    typedef EMData          data_Type;

    MooneyRivlin (Real mu = 4960) : M_W2 (new W2(mu)) { } // 0.496 KPa

    class W2
    {
    public:
        typedef typename MaterialFunctions::EMMaterialFunctions<Mesh>::return_Type return_Type;
        virtual return_Type operator() (const MatrixSmall<3, 3>& F)
        {
            return 0.5 * mu;
        }

        W2 (Real mu = 4960) : mu (mu) {} // 0.496 KPa
        W2 (const W2& mooneyRivlin)
        {
            mu = mooneyRivlin.mu;
        }
        virtual ~W2() {}

        void setC2( Real mu )
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
        EMAssembler::computeI2JacobianTerms ( disp, dispETFESpace, jacobianPtr, M_W2 );
    }

    virtual void computeResidual ( const vector_Type& disp,
                                          boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > >  dispETFESpace,
                                          const vector_Type& fibers,
                                          const vector_Type& sheets,
                                          vectorPtr_Type residualVectorPtr)
    {
        EMAssembler::computeI2ResidualTerms ( disp, dispETFESpace, residualVectorPtr, M_W2 );
    }

    virtual void setParameters (data_Type& data)
    {
    	M_W2->setC2( data.solidParameter<Real>("C2") );
    }

    void showMe()
    {
    	std::cout << "Shear Modulus = " << M_W2->mu << "\n";
    }

protected:
    boost::shared_ptr<MooneyRivlin::W2> M_W2;
};




} //EMMaterialFunctions

} //LifeV
#endif /* EMMATERIALFUNCTIONS_HPP_ */
