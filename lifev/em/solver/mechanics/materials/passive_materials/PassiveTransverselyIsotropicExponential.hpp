/*
 * EMPassiveNeoHookean.hpp
 *
 *  Created on: 29/apr/2014
 *      Author: srossi
 */

#ifndef PASSIVETRANSVERSELYISOTROPICEXPONENTIAL_HPP_
#define PASSIVETRANSVERSELYISOTROPICEXPONENTIAL_HPP_


#include <lifev/em/solver/mechanics/materials/EMMaterialType.hpp>

//#include <lifev/em/solver/mechanics/materials/EMMaterialFunctions.hpp>



//class EMPassiveNeoHookean : public virtual EMPassiveMaterial
//{
//public:
//	EMPassiveNeoHookean();
//	virtual ~EMPassiveNeoHookean() {}
//};

namespace LifeV
{

template<typename Mesh>
class PassiveTransverselyIsotropicExponential : public virtual EMMaterialType<Mesh>
{
public:
	typedef EMMaterialType<Mesh> super;

	PassiveTransverselyIsotropicExponential();
	virtual ~PassiveTransverselyIsotropicExponential() {}
};

template<typename Mesh>
PassiveTransverselyIsotropicExponential<Mesh>::PassiveTransverselyIsotropicExponential() :
		super("Passive Transversely Isotropic Exponential", 6)
{
	this -> M_materialFunctionList[0].reset(new MaterialFunctions::Volumetric<Mesh>(100000.0)  );
	this -> M_materialFunctionList[1].reset(new MaterialFunctions::dVolumetric<Mesh>(100000.0) );
	this -> M_materialFunctionList[2].reset(new MaterialFunctions::IsotropicExponential<Mesh>()  );
	this -> M_materialFunctionList[3].reset(new MaterialFunctions::dIsotropicExponential<Mesh>()  );
	this -> M_materialFunctionList[4].reset(new MaterialFunctions::AnisotropicExponential<Mesh>()  );
	this -> M_materialFunctionList[5].reset(new MaterialFunctions::AnisotropicExponential<Mesh>()  );
}


template <typename MeshType>
inline EMMaterialType<MeshType>* createPassiveTransverselyIsotropicExponential()
{
    return new PassiveTransverselyIsotropicExponential<MeshType>();
}
namespace
{
static bool registerEM_passiveTIE = EMMaterialType<LifeV::RegionMesh<LinearTetra> >::EMMaterialFactory::instance().registerProduct ("PTIE", &createPassiveTransverselyIsotropicExponential<LifeV::RegionMesh<LinearTetra> > );
}

}//LifeV

#endif /* EMPASSIVENEOHOOKEAN_HPP_ */
