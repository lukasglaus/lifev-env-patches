/*
 * ActiveStrainActivation.hpp
 *
 *  Created on: 15/mag/2014
 *      Author: srossi
 */

#ifndef ACTIVESTRAINACTIVATION_HPP_
#define ACTIVESTRAINACTIVATION_HPP_

#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/array/MapEpetra.hpp>







namespace LifeV
{


namespace ActiveStrain
{

void transversilyIsotropicActiveStrain( VectorEpetra& gf,
										VectorEpetra& gs,
										VectorEpetra& gn )
{
    gs = 1.0;
    gs /= (1.0 + gf);
    EMUtility::EpetraSqrt (gs);
    gs -= 1.0;
	gn = gs;
}

void orthotropicActiveStrain( VectorEpetra& gf,
							  VectorEpetra& gs,
							  VectorEpetra& gn,
							  Real factor)
{
     gn = factor * gf;
     gs = 1.0;
     gs /= (1.0 + gf);
     gs /= (1.0 + gn);
     gs -= 1.0;
}

}//namespace ActiveStrain


class ActiveStrainActivation
{
public:
    //! Distributed vector // For parallel usage
    typedef VectorEpetra                                                vector_Type;

    typedef boost::shared_ptr<VectorEpetra>                             vectorPtr_Type;

    ActiveStrainActivation (MapEpetra& map) :
    	                  M_gammafPtr ( new vector_Type (map) ),
    	                  M_gammasPtr ( new vector_Type (map) ),
    	                  M_gammanPtr ( new vector_Type (map) ) {}

    ActiveStrainActivation (const MapEpetra& map) :
						M_gammafPtr ( new vector_Type (map) ),
						M_gammasPtr ( new vector_Type (map) ),
						M_gammanPtr ( new vector_Type (map) ) {}

    ActiveStrainActivation (ActiveStrainActivation& activation) :
    	                M_gammafPtr ( new vector_Type (activation.gammaf() ) ),
    	                M_gammasPtr ( new vector_Type (activation.gammas() ) ),
    	                M_gammanPtr ( new vector_Type (activation.gamman() ) ) {}

    ActiveStrainActivation& operator= (ActiveStrainActivation& activation)
    {
        M_gammafPtr.reset ( new vector_Type (activation.gammaf() ) );
        M_gammasPtr.reset ( new vector_Type (activation.gammas() ) );
        M_gammanPtr.reset ( new vector_Type (activation.gamman() ) );
        return *this;
    }

    virtual ~ActiveStrainActivation() {}

    VectorEpetra& gammaf()
    {
        return *M_gammafPtr;
    }

    vectorPtr_Type gammafPtr()
    {
        return M_gammafPtr;
    }

    void setGammaf (VectorEpetra& activation)
    {
        *M_gammafPtr = activation;
    }

    void setGammafPtr (vectorPtr_Type activationPtr)
    {
        M_gammafPtr = activationPtr;
    }

    VectorEpetra& gammas()
    {
        return *M_gammasPtr;
    }

    vectorPtr_Type gammasPtr()
    {
        return M_gammasPtr;
    }

    void setGammas (VectorEpetra& activation)
    {
        *M_gammasPtr = activation;
    }

    void setGammasPtr (vectorPtr_Type activationPtr)
    {
        M_gammasPtr = activationPtr;
    }

    VectorEpetra& gamman()
    {
        return *M_gammanPtr;
    }

    vectorPtr_Type gammanPtr()
    {
        return M_gammanPtr;
    }

    void setGamman (VectorEpetra& activation)
    {
        *M_gammanPtr = activation;
    }

    void setGammanPtr (vectorPtr_Type activationPtr)
    {
        M_gammanPtr = activationPtr;
    }

    virtual void solveModel() {}

protected:
    vectorPtr_Type M_gammafPtr;
    vectorPtr_Type M_gammasPtr;
    vectorPtr_Type M_gammanPtr;

};




} /* namespace LifeV */

#endif /* ActiveStrainActivation_HPP_ */
