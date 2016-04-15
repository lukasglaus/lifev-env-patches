/*
 * EMMaterialFunctions.hpp
 *
 *  Created on: 15/mar/2015
 *      Author: tkummer
 */

#ifndef FUNCTIONSMOONEYRIVLINCOMPACT_HPP_
#define FUNCTIONSMOONEYRIVLINCOMPACT_HPP_

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
class MooneyRivlinCompact : public virtual EMMaterialFunctions<Mesh>
{
public:
    
    typedef EMData          data_Type;

    MooneyRivlinCompact (Real bulk1 = 2000, Real bulk2 = 2000, Real mu1 = 4960, Real mu2 = 4960) :
        M_Vol (new Volumetric(bulk1)),
        M_dVol (new dVolumetric(bulk2)),
        M_W1 (new W1(mu1)),
        M_W2 (new W2(mu2))
    { } // 0.496 KPa

    
    class Volumetric
    {
    public:
        typedef typename MaterialFunctions::EMMaterialFunctions<Mesh>::return_Type return_Type;
        
        virtual return_Type operator() (const MatrixSmall<3, 3>& F)
        {
            Real J = Elasticity::J (F);
            //        return M_bulk * ( (J - 1) + std::log (J) / J );
            return ( M_bulk * ( J + J * std::log(J) - 1. ) ) / ( 2 * J );
        }
        
        Volumetric (Real bulk = 2000.0) : M_bulk (bulk) {}
        Volumetric (const Volumetric& v)
        {
            M_bulk = v.M_bulk;
        }
        virtual ~Volumetric() {}
        
        typedef EMData          data_Type;
        void setParameters (data_Type& data)
        {
            M_bulk = data.solidParameter<Real>("BulkModulus");
        }
        void showMe()
        {
            std::cout << "BulkModulus = " << M_bulk << "\n";
        }
    private:
        Real M_bulk;
    };
    
    class dVolumetric
    {
    public:
        typedef typename MaterialFunctions::EMMaterialFunctions<Mesh>::return_Type return_Type;
        
        virtual return_Type operator() (const MatrixSmall<3, 3>& F)
        {
            Real J = Elasticity::J (F);
            //        return M_bulk  / J / J * ( 1 + J * J - std::log (J) ) ;
            return ( M_bulk * ( J + 1. ) ) / ( 2. * J * J);
        }
        
        dVolumetric (Real bulk = 2000.0) : M_bulk (bulk) {}
        dVolumetric (const dVolumetric& v)
        {
            M_bulk = v.M_bulk;
        }
        virtual ~dVolumetric() {}
        
        typedef EMData          data_Type;
        void setParameters (data_Type& data)
        {
            M_bulk = data.solidParameter<Real>("BulkModulus");
        }
        void showMe()
        {
        }
    private:
        Real M_bulk;
    };
    
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
        W1 (const W1& mooneyRivlinCompact)
        {
            mu = mooneyRivlinCompact.mu;
        }
        virtual ~W1() {}
        
        void setMu( Real mu )
        {
            this->mu = mu;
        }
        
        Real mu;
        
    };
    
    class W2
    {
    public:
        typedef typename MaterialFunctions::EMMaterialFunctions<Mesh>::return_Type return_Type;
        virtual return_Type operator() (const MatrixSmall<3, 3>& F)
        {
            return 0.5 * mu;
        }

        W2 (Real mu = 4960) : mu (mu) {} // 0.496 KPa
        W2 (const W2& mooneyRivlinCompact)
        {
            mu = mooneyRivlinCompact.mu;
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
        EMAssembler::computePMRCJacobianTerms ( disp, dispETFESpace, jacobianPtr, M_Vol, M_dVol, M_W1, M_W2 );
    }

    virtual void computeResidual ( const vector_Type& disp,
                                          boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > >  dispETFESpace,
                                          const vector_Type& fibers,
                                          const vector_Type& sheets,
                                          vectorPtr_Type residualVectorPtr)
    {
        EMAssembler::computePMRCResidualTerms ( disp, dispETFESpace, residualVectorPtr, M_Vol, M_dVol, M_W1, M_W2 );
    }

    virtual void setParameters (data_Type& data)
    {
        M_W1->setMu( data.solidParameter<Real>("mu") );
        M_W2->setC2( data.solidParameter<Real>("C2") );
    }

    void showMe()
    {
    	std::cout << "Shear Modulus = " << M_W1->mu << "\n";
        std::cout << "Shear Modulus = " << M_W2->mu << "\n";
    }

protected:
    boost::shared_ptr<MooneyRivlinCompact::Volumetric> M_Vol;
    boost::shared_ptr<MooneyRivlinCompact::dVolumetric> M_dVol;
    boost::shared_ptr<MooneyRivlinCompact::W1> M_W1;
    boost::shared_ptr<MooneyRivlinCompact::W2> M_W2;
};




} //EMMaterialFunctions

} //LifeV
#endif /* EMMATERIALFUNCTIONSCOMPACT_HPP_ */
