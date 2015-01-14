/*
 * NashPanfilovModel14.hpp
 *
 *  Created on: 13/mag/2014
 *      Author: srossi
 */

#ifndef ACTIVESTRESSNASHPANFILOVMODEL04_HPP_
#define ACTIVESTRESSNASHPANFILOVMODEL04_HPP_

#include <lifev/em/solver/activation/activeStressModels/ActiveStressActivation.hpp>

namespace LifeV
{

class ActiveStressNashPanfilovModel04 :  public virtual ActiveStressActivation
{
public:

    //! Distributed vector // For parallel usage
    typedef VectorEpetra                                                vector_Type;

    typedef boost::shared_ptr<VectorEpetra>                             vectorPtr_Type;

    typedef ActiveStressActivation                                      super;

    ActiveStressNashPanfilovModel04 (Real kTa = 49.7, Real epsilon0 = 1.0);

    virtual ~ActiveStressNashPanfilovModel04() {}

    void solveModel (Real& timeStep);

    void multiplyByEpsilon (VectorEpetra& rhs);

    inline Real kTa()
    {
        return M_kTa;
    }
    inline void setKTa (Real kTa)
    {
        M_kTa = kTa;
    }

    inline Real epsilon0()
    {
        return M_epsilon0;
    }
    inline void setEpsilon0 (Real epsilon0)
    {
        M_epsilon0 = epsilon0;
    }


    void setParameters( EMData& data);
    void setup( EMData& data, const MapEpetra& map)
    {
    	setParameters(data);
        this->M_fiberActivationPtr.reset ( new vector_Type ( map ) );
        super::setup(data, map);
    }


    void setupActivationPtrs(	vectorPtr_Type& fiberActivationPtr,
								vectorPtr_Type& sheetActivationPtr,
								vectorPtr_Type& normalActivationPtr );

    void updateActivation(	vectorPtr_Type& fiberActivationPtr,
							vectorPtr_Type& sheetActivationPtr,
							vectorPtr_Type& normalActivationPtr );

private:

    Real M_kTa;
    Real M_epsilon0;

};




inline Activation* createActiveStressNashPanfilovModel04()
{
    return new ActiveStressNashPanfilovModel04();
}

namespace
{
static bool registerActivation_ActiveStressNashPanfilov04 = Activation::EMActivationFactory::instance().registerProduct ("ActiveStressNashPanfilov04", &createActiveStressNashPanfilovModel04);
}

} /* namespace LifeV */

#endif /* NashPanfilovModel14_HPP_ */
