/*
 * ActiveStrainActivation.hpp
 *
 *  Created on: 15/mag/2014
 *      Author: srossi
 */

#ifndef ACTIVESTRAINACTIVATION_HPP_
#define ACTIVESTRAINACTIVATION_HPP_



#include <lifev/em/solver/activation/Activation.hpp>



namespace LifeV
{


//namespace ActiveStrain
//{
//
//void transversilyIsotropicActiveStrain( VectorEpetra& gf,
//										VectorEpetra& gs,
//										VectorEpetra& gn )
//{
//    gs = 1.0;
//    gs /= (1.0 + gf);
//    EMUtility::EpetraSqrt (gs);
//    gs -= 1.0;
//	gn = gs;
//}
//
//void orthotropicActiveStrain( VectorEpetra& gf,
//							  VectorEpetra& gs,
//							  VectorEpetra& gn,
//							  Real factor)
//{
//     gn = factor * gf;
//     gs = 1.0;
//     gs /= (1.0 + gf);
//     gs /= (1.0 + gn);
//     gs -= 1.0;
//}
//
//}//namespace ActiveStrain


class ActiveStrainActivation : public virtual Activation
{
public:
    //! Distributed vector // For parallel usage
    typedef VectorEpetra                                                vector_Type;

    typedef boost::shared_ptr<VectorEpetra>                             vectorPtr_Type;

    typedef Activation super;

    enum ActiveStrainType { Anisotropic, Orthotropic, TransverselyIsotropic };

    ActiveStrainActivation () :
    	M_activeStrainType(),
    	M_activeStrainOrthotropicParameter(-666.0)
        {}


    void init(const MapEpetra& map)
    {
    	switch(M_activeStrainType)
    	{
    		case TransverselyIsotropic:
    		case Orthotropic:
    			 M_fiberActivationPtr.reset( new vector_Type (map) );
    			 break;
    		case Anisotropic:
    		default:
    			M_fiberActivationPtr.reset( new vector_Type (map) );
    			M_sheetActivationPtr.reset( new vector_Type (map) );
    			M_normalActivationPtr.reset( new vector_Type (map) );
				 break;
    	}
    }

    virtual ~ActiveStrainActivation() {}

    virtual void setup(EMData& data,  const MapEpetra& map)
    {
    	M_activeStrainOrthotropicParameter = data.activationParameter<Real>("EMActiveStrainOrthotropicParameter");
    	init(map);
    }



//


    ActiveStrainType activeStrainType()
    {
    	return M_activeStrainType;
    }

    void setupActivationPtrs(	vectorPtr_Type& fiberActivationPtr,
								vectorPtr_Type& sheetActivationPtr,
								vectorPtr_Type& normalActivationPtr )
    {
     	switch(M_activeStrainType)
     	{
 			case TransverselyIsotropic:
 			case Orthotropic:
 				fiberActivationPtr  = this->M_fiberActivationPtr;
 				sheetActivationPtr.reset();
 				normalActivationPtr.reset();
 				break;
 			case Anisotropic:
 			default:
 				fiberActivationPtr = this->M_fiberActivationPtr;
 				sheetActivationPtr = this->M_sheetActivationPtr;
 				normalActivationPtr = this->M_normalActivationPtr;
 				break;
     	}
     }

    void updateActivation(	vectorPtr_Type& fiberActivationPtr,
							vectorPtr_Type& sheetActivationPtr,
							vectorPtr_Type& normalActivationPtr )
    {
    	setupActivationPtrs(fiberActivationPtr, sheetActivationPtr, normalActivationPtr);
    }


protected:

    ActiveStrainType M_activeStrainType;

//    vectorPtr_Type M_gammafPtr;
//    vectorPtr_Type M_gammasPtr;
//    vectorPtr_Type M_gammanPtr;

    Real M_activeStrainOrthotropicParameter;


};




} /* namespace LifeV */

#endif /* ActiveStrainActivation_HPP_ */
