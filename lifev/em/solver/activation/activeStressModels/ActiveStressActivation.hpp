/*
 * ActiveStressActivation.hpp
 *
 *  Created on: 15/mag/2014
 *      Author: srossi
 */

#ifndef ACTIVESTRESSACTIVATION_HPP_
#define ACTIVESTRESSACTIVATION_HPP_

#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/array/MapEpetra.hpp>


namespace LifeV {

class ActiveStressActivation {
public:
    //! Distributed vector // For parallel usage
    typedef VectorEpetra 												vector_Type;

    typedef boost::shared_ptr<VectorEpetra> 							vectorPtr_Type;

	ActiveStressActivation(MapEpetra& map) : M_activation(map){}
	ActiveStressActivation(const MapEpetra& map) : M_activation(map){}
	virtual ~ActiveStressActivation() {}

	inline VectorEpetra& activation()
	{
		return M_activation;
	}

	inline void setActivation(VectorEpetra& activation)
	{
		M_activation = activation;
	}

	inline void setActivation(vectorPtr_Type activationPtr)
	{
		M_activation = *activationPtr;
	}

protected:
	VectorEpetra M_activation;

};


} /* namespace LifeV */

#endif /* ACTIVESTRESSACTIVATION_HPP_ */
