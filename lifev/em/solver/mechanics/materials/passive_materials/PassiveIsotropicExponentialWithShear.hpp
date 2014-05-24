/*
 * EMPassiveNeoHookean.hpp
 *
 *  Created on: 29/apr/2014
 *      Author: srossi
 */

#ifndef PassiveIsotropicExponentialWithShear_HPP_
#define PassiveIsotropicExponentialWithShear_HPP_


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
class PassiveIsotropicExponentialWithShear : public virtual EMMaterialType<Mesh>
{
public:
	typedef EMMaterialType<Mesh> super;

	PassiveIsotropicExponentialWithShear();
	virtual ~PassiveIsotropicExponentialWithShear() {}
};

template<typename Mesh>
PassiveIsotropicExponentialWithShear<Mesh>::PassiveIsotropicExponentialWithShear() :
		super("Passive Isotropic Exponential With Shear", 6)
{
	this -> M_materialFunctionList[0].reset( new MaterialFunctions::Volumetric<Mesh>(100000.0)  );
	this -> M_materialFunctionList[1].reset( new MaterialFunctions::dVolumetric<Mesh>(100000.0) );
	this -> M_materialFunctionList[2].reset( new MaterialFunctions::IsotropicExponential<Mesh>()  );
	this -> M_materialFunctionList[3].reset( new MaterialFunctions::dIsotropicExponential<Mesh>()  );
	this -> M_materialFunctionList[4].reset( new MaterialFunctions::ShearExponential<Mesh>()  );
	this -> M_materialFunctionList[5].reset( new MaterialFunctions::dShearExponential<Mesh>()  );
}


template <typename MeshType>
inline EMMaterialType<MeshType>* createPassiveIsotropicExponentialWithShear()
{
    return new PassiveIsotropicExponentialWithShear<MeshType>();
}
namespace
{
static bool registerEM_passiveIEWS = EMMaterialType<LifeV::RegionMesh<LinearTetra> >::EMMaterialFactory::instance().registerProduct ("PIEWS", &createPassiveIsotropicExponentialWithShear<LifeV::RegionMesh<LinearTetra> > );
}

}//LifeV

#endif /* EMPASSIVENEOHOOKEAN_HPP_ */
