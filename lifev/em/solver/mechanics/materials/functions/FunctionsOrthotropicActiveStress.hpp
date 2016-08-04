/*
 * EMMaterialFunctions.hpp
 *
 *  Created on: 28/apr/2014
 *      Author: srossi
 */

#ifndef FUNCTIONSOrthotropicActiveStress_HPP_
#define FUNCTIONSOrthotropicActiveStress_HPP_

#include <lifev/em/solver/mechanics/EMElasticityFunctions.hpp>
#include <lifev/em/solver/mechanics/materials/functions/EMMaterialFunctions.hpp>

namespace LifeV
{

namespace MaterialFunctions
{

typedef VectorEpetra           vector_Type;
typedef boost::shared_ptr<vector_Type>         vectorPtr_Type;

typedef MatrixEpetra<Real>           matrix_Type;
typedef boost::shared_ptr<matrix_Type>         matrixPtr_Type;

////////////////////////////////////////////////////////////////////////
//  ORTHOTROPIC ACTIVE STRESS FUNCTIONS
////////////////////////////////////////////////////////////////////////

template <class Mesh>
class OrthotropicActiveStress : public virtual EMMaterialFunctions<Mesh>
{
public:
    typedef typename MaterialFunctions::EMMaterialFunctions<Mesh>::return_Type return_Type;

    virtual return_Type operator() (const MatrixSmall<3, 3>& F)
    {
        return 1.0;
    }

    virtual return_Type operator() (const Real& H)
    {
        return 0.5 * H * H * M_C * M_Tmax;
    }

    OrthotropicActiveStress (LifeV::Real Tmax = 49700, std::string field = "Fibers") :
        M_Tmax  (Tmax),
        M_Field (field),
        M_C     (0.0)
    {} // 0.33 KPa
    
    OrthotropicActiveStress (const OrthotropicActiveStress& OrthotropicActiveStress)
    {
        M_Tmax = OrthotropicActiveStress.M_Tmax;
    }

    virtual ~OrthotropicActiveStress() {}

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
        EMAssembler::computeOrthotropicActiveStressJacobianTerms (disp,
                                                            dispETFESpace,
                                                            fibers,
                                                            sheets,
                                                            *fiberActivation,
                                                            activationETFESpace,
                                                            jacobianPtr,
                                                            this->getMe(),
                                                            M_Field);
        
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
        EMAssembler::computeOrthotropicActiveStressResidualTerms (disp,
                                                            dispETFESpace,
                                                            fibers,
                                                            sheets,
                                                            *fiberActivation,
                                                            activationETFESpace,
                                                            residualVectorPtr,
                                                            this->getMe(),
                                                            M_Field);
    }


    typedef EMData          data_Type;
    void setParameters (data_Type& data)
    {
    	M_Tmax = data.solidParameter<Real>("MaxActiveTension");
        if ( M_Field == "Sheets")
        {
            M_C = data.solidParameter<Real>("Cs");
        }
        else if ( M_Field == "Normal" )
        {
            M_C = data.solidParameter<Real>("Cn");
        }
        else
        {
            M_C = 1.0;
        }
    }

    void showMe()
    {
    	std::cout << "Active Tension = " << M_Tmax << "\n";
    }

private:
    LifeV::Real M_Tmax;
    std::string M_Field;
    Real M_C;
};


} //EMMaterialFunctions

} //LifeV
#endif /* EMMATERIALFUNCTIONS_HPP_ */
