/*
 * RossiModel14.hpp
 *
 *  Created on: 13/mag/2014
 *      Author: srossi
 */

#ifndef ACTIVESTRAINROSSIMODEL14_HPP_
#define ACTIVESTRAINROSSIMODEL14_HPP_

#include <lifev/em/solver/activation/activeStrainModels/ActiveStrainActivation.hpp>

namespace LifeV
{

class ActiveStrainRossiModel14 : public virtual ActiveStrainActivation
{
public:

//	enum ActivationAnisotropy { TransversilyIsotropic, Orthotropic };
    //! Distributed vector // For parallel usage
    typedef VectorEpetra                                                vector_Type;

    typedef boost::shared_ptr<VectorEpetra>                             vectorPtr_Type;

    typedef ActiveStrainActivation                                      super;

    ActiveStrainRossiModel14 ( Real inverseViscosity = 0.0001,
    		                   Real coefficient = 1.,
    		                   Real threshold = 0.21  );

    virtual ~ActiveStrainRossiModel14() {}

//    void solveModel (VectorEpetra& ActiveStress, Real timeStep);

    void solveModel ( Real& timeStep);

    Real computeActiveStress(Real i4f, Real Calcium);


    Real inverseViscosity()
    {
        return M_inverseViscosity;
    }

    void setInverseViscosity (Real inverseViscosity)
    {
        M_inverseViscosity = inverseViscosity;
    }

    Real activeForceCoefficient()
    {
        return M_activeForceCoefficient;
    }

    void setActiveForceCoefficient (Real coefficient)
    {
        M_activeForceCoefficient = coefficient;
    }

    void setParameters( EMData& data);
    void setup( EMData& data, const MapEpetra& map)
    {
    	setParameters(data);
    	super::setup(data, map);
    }

    void computeI4f();

    void showMe()
    {
    	ActiveStrainActivation::showMe();
    	if(M_fiberActivationPtr->comm().MyPID() == 0)
    	{
    		std::cout << "Activation - Using ActiveStrainRossiModel14 with the following setup: " << std::endl;
    		std::cout << "Inverse of the viscosity = " << M_inverseViscosity << std::endl;
    		std::cout << "IActive Force coefficient = " << M_activeForceCoefficient << std::endl;
    		std::cout << "Threshold of the chemical species = " << M_chemicalThreshold << std::endl;
    		std::cout << "Calcium index in the electrophysiology vector = " << M_calciumIndex << std::endl;
    	}
    }

private:

    VectorSmall<3> M_PathologyCenter;
    Real PathologyRadius;
    Real PathologyStrength;
    
    Real M_inverseViscosity;
    Real M_activeForceCoefficient;
    Real M_chemicalThreshold;
    UInt M_calciumIndex;

};




inline Activation* createActiveStrainRossiModel14()
{
    return new ActiveStrainRossiModel14();
}

namespace
{
static bool registerActivation_ActiveStrainRossi14 = Activation::EMActivationFactory::instance().registerProduct ("ActiveStrainRossi14", &createActiveStrainRossiModel14 );
}



} /* namespace LifeV */

#endif /* ROSSIMODEL14_HPP_ */
