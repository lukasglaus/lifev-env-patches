/*
 * EMMaterialFunctions.hpp
 *
 *  Created on: 28/apr/2014
 *      Author: srossi
 */

#ifndef FUNCTIONSVOLUMETRIC_HPP_
#define FUNCTIONSVOLUMETRIC_HPP_

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
//  VOLUMETRIC FUNCTIONS
////////////////////////////////////////////////////////////////////////
template <class Mesh>
class Volumetric : public virtual MaterialFunctions::EMMaterialFunctions<Mesh>
{
public:
    typedef typename MaterialFunctions::EMMaterialFunctions<Mesh>::return_Type return_Type;

    virtual return_Type operator() (const MatrixSmall<3, 3>& F)
    {
        Real J = Elasticity::J (F);
        return M_bulk * ( (J - 1) + std::log (J) / J );
    }

    Volumetric (Real bulk = 2000.0) : M_bulk (bulk) {}
    Volumetric (const Volumetric& v)
    {
        M_bulk = v.M_bulk;
    }
    virtual ~Volumetric() {}

    inline virtual void computeJacobian ( const vector_Type& disp,
                                          boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > >  dispETFESpace,
                                          const vector_Type& fibers,
                                          const vector_Type& sheets,
                                          matrixPtr_Type           jacobianPtr)
    {
        EMAssembler::computeVolumetricJacobianTerms (disp, dispETFESpace, jacobianPtr, this->getMe() );
    }

    inline virtual void computeResidual ( const vector_Type& disp,
                                          boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > >  dispETFESpace,
                                          const vector_Type& fibers,
                                          const vector_Type& sheets,
                                          vectorPtr_Type           residualVectorPtr)
    {
        EMAssembler::computeVolumetricResidualTerms (disp, dispETFESpace, residualVectorPtr, this->getMe() );
    }

    typedef EMData          data_Type;
    void setParameters (data_Type& data)
    {
    	M_bulk = data.parameter("BulkModulus");
    }
    void showMe()
    {
    	std::cout << "BulkModulus = " << M_bulk << "\n";
    }
private:
    Real M_bulk;
};

template <class Mesh>
class dVolumetric : public virtual EMMaterialFunctions<Mesh>
{
public:
    typedef typename MaterialFunctions::EMMaterialFunctions<Mesh>::return_Type return_Type;

    virtual return_Type operator() (const MatrixSmall<3, 3>& F)
    {
        Real J = Elasticity::J (F);
        return M_bulk  / J / J * ( 1 + J * J - std::log (J) ) ;
    }

    dVolumetric (Real bulk = 2000.0) : M_bulk (bulk) {}
    dVolumetric (const dVolumetric& v)
    {
        M_bulk = v.M_bulk;
    }
    virtual ~dVolumetric() {}

    inline virtual void computeJacobian ( const vector_Type& disp,
                                          boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > >  dispETFESpace,
                                          const vector_Type& fibers,
                                          const vector_Type& sheets,
                                          matrixPtr_Type           jacobianPtr)
    {
        EMAssembler::computeVolumetricJacobianTermsSecondDerivative (disp, dispETFESpace, jacobianPtr, this->getMe() );
    }

    typedef EMData          data_Type;
    void setParameters (data_Type& data)
    {
    	M_bulk = data.parameter("BulkModulus");
    }
    void showMe()
    {
    }
private:
    Real M_bulk;
};



} //EMMaterialFunctions

} //LifeV
#endif /* EMMATERIALFUNCTIONS_HPP_ */
