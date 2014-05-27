/*
 * EMPassiveNeoHookean.hpp
 *
 *  Created on: 29/apr/2014
 *      Author: srossi
 */

#ifndef PASSIVENEOHOOKEAN_HPP_
#define PASSIVENEOHOOKEAN_HPP_


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
class PassiveNeoHookean : public virtual EMMaterialType<Mesh>
{
public:
	typedef EMMaterialType<Mesh> super;

	PassiveNeoHookean();
	virtual ~PassiveNeoHookean() {}

};

template<typename Mesh>
PassiveNeoHookean<Mesh>::PassiveNeoHookean() :
		super("Passive Neo-Hookean", 3)
{
	this -> M_materialFunctionList[0].reset(new MaterialFunctions::Volumetric<Mesh> (3500000.0)  );
	this -> M_materialFunctionList[1].reset(new MaterialFunctions::dVolumetric<Mesh>(3500000.0) );
	this -> M_materialFunctionList[2].reset(new MaterialFunctions::NeoHookean<Mesh> ()  );
}

template <typename MeshType>
inline EMMaterialType<MeshType>* createPassiveNeoHookean()
{
    return new PassiveNeoHookean<MeshType>();
}
namespace
{
static bool registerEM_passiveNH = EMMaterialType<LifeV::RegionMesh<LinearTetra> >::EMMaterialFactory::instance().registerProduct ("PNH", &createPassiveNeoHookean<LifeV::RegionMesh<LinearTetra> > );
}

}//LifeV

#endif /* EMPASSIVENEOHOOKEAN_HPP_ */
