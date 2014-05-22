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

	ActiveStressActivation(MapEpetra& map) : M_activationPtr( new vector_Type(map) ){}
	ActiveStressActivation(const MapEpetra& map) : M_activationPtr( new vector_Type(map) ){}

	ActiveStressActivation(ActiveStressActivation& activation) : M_activationPtr( new vector_Type(activation.activation()) ){}

	inline ActiveStressActivation& operator=(ActiveStressActivation& activation)
	{
		M_activationPtr.reset( new vector_Type(activation.activation()) );
		return *this;
	}

	virtual ~ActiveStressActivation() {}

	inline VectorEpetra& activation()
	{
		return *M_activationPtr;
	}

	inline vectorPtr_Type activationPtr()
	{
		return M_activationPtr;
	}

	inline void setActivation(VectorEpetra& activation)
	{
		*M_activationPtr = activation;
	}

	inline void setActivationPtr(vectorPtr_Type activationPtr)
	{
		M_activationPtr = activationPtr;
	}

	virtual void solveModel() {}

protected:
	vectorPtr_Type M_activationPtr;

};




} /* namespace LifeV */

#endif /* ACTIVESTRESSACTIVATION_HPP_ */
