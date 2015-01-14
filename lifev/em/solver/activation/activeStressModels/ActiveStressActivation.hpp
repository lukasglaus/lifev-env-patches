/*
 * ActiveStressActivation.hpp
 *
 *  Created on: 15/mag/2014
 *      Author: srossi
 */

#ifndef ACTIVESTRESSACTIVATION_HPP_
#define ACTIVESTRESSACTIVATION_HPP_

#include <lifev/em/solver/activation/Activation.hpp>

namespace LifeV
{

class ActiveStressActivation : public virtual Activation
{
public:

	virtual void setup(const EMData& data, const MapEpetra& map)
	{

	}


};




} /* namespace LifeV */

#endif /* ACTIVESTRESSACTIVATION_HPP_ */
