/*
 * EMMaterialFunctions.hpp
 *
 *  Created on: 28/apr/2014
 *      Author: srossi
 */

#ifndef FUNCTIONSSIMPLEACTIVESTRESS_HPP_
#define FUNCTIONSSIMPLEACTIVESTRESS_HPP_

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
//  ISOTROPIC EXPONENTIAL FUNCTIONS
////////////////////////////////////////////////////////////////////////

template <class Mesh>
class SimpleActiveStress : public virtual EMMaterialFunctions<Mesh>
{
public:
    typedef typename MaterialFunctions::EMMaterialFunctions<Mesh>::return_Type return_Type;

    virtual return_Type operator() (const MatrixSmall<3, 3>& F)
    {
        return 1.0;
    }

    virtual return_Type operator() (const Real& H)
    {
        return 0.5 * H * H * M_Tmax;
    }

    SimpleActiveStress (LifeV::Real Tmax = 49700) : M_Tmax (Tmax) {} // 0.33 KPa
    SimpleActiveStress (const SimpleActiveStress& simpleActiveStress)
    {
        M_Tmax = simpleActiveStress.M_Tmax;
    }

    virtual ~SimpleActiveStress() {}

    inline void computeJacobian ( const vector_Type& disp,
                                  boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > >  dispETFESpace,
                                  const vector_Type& fibers,
                                  const vector_Type& sheets,
                                  const vector_Type& activation,
                                  boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 1 > >  activationETFESpace,
                                  matrixPtr_Type           jacobianPtr)
    {
        //      EMAssembler::computeFiberActiveStressJacobianTerms(disp, dispETFESpace, fibers, sheets, activation, activationETFESpace,jacobianPtr, this->getMe() );
        EMAssembler::computeFiberActiveStressJacobianTerms (disp,
                                                            dispETFESpace,
                                                            fibers,
                                                            sheets,
                                                            activation,
                                                            activationETFESpace,
                                                            jacobianPtr,
                                                            this->getMe() );
    }

    inline virtual void computeResidual ( const vector_Type& disp,
                                          boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > >  dispETFESpace,
                                          const vector_Type& fibers,
                                          const vector_Type& sheets,
                                          const vector_Type& activation,
                                          boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 1 > >  activationETFESpace,
                                          vectorPtr_Type           residualVectorPtr)
    {
        EMAssembler::computeFiberActiveStressResidualTerms (disp,
                                                            dispETFESpace,
                                                            fibers,
                                                            sheets,
                                                            activation,
                                                            activationETFESpace,
                                                            residualVectorPtr,
                                                            this->getMe() );
    }

private:
    LifeV::Real M_Tmax;
};


} //EMMaterialFunctions

} //LifeV
#endif /* EMMATERIALFUNCTIONS_HPP_ */
