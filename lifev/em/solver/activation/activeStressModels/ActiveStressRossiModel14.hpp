/*
 * RossiModel14.hpp
 *
 *  Created on: 13/mag/2014
 *      Author: srossi
 */

#ifndef ACTIVESTRESSROSSIMODEL14_HPP_
#define ACTIVESTRESSROSSIMODEL14_HPP_

#include <lifev/em/solver/activation/activeStressModels/ActiveStressActivation.hpp>

namespace LifeV
{

class ActiveStressRossiModel14 : public virtual ActiveStressActivation
{
public:

    //! Distributed vector // For parallel usage
    typedef VectorEpetra                                                vector_Type;

    typedef boost::shared_ptr<VectorEpetra>                             vectorPtr_Type;

    typedef ActiveStressActivation                                      super;

    ActiveStressRossiModel14 (MapEpetra& map, Real beta = 2.279, Real mu = 1000, Real Tmax = 50);
    ActiveStressRossiModel14 (const MapEpetra& map, Real beta = 2.279, Real mu = 1000, Real Tmax = 50);

    virtual ~ActiveStressRossiModel14() {}

    void solveModel (VectorEpetra& potential, Real timeStep);


    inline Real coefficientBeta()
    {
        return M_coefficientBeta;
    }
    inline void setCoefficientBeta (Real beta)
    {
        M_coefficientBeta = beta;
    }

    inline Real coefficientMu()
    {
        return M_coefficientMu;
    }
    inline void setCoefficientMu (Real mu)
    {
        M_coefficientBeta = mu;
    }

    inline Real maximumActiveTenstion()
    {
        return M_maximumActiveTenstion;
    }
    inline void setMaximumActiveTenstion (Real Tmax)
    {
        M_maximumActiveTenstion = Tmax;
    }
private:

    Real M_coefficientBeta;
    Real M_coefficientMu;
    Real M_maximumActiveTenstion;

};

} /* namespace LifeV */

#endif /* ROSSIMODEL14_HPP_ */
