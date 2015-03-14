/*
 * EMPassiveNeoHookean.hpp
 *
 *  Created on: 29/apr/2014
 *      Author: srossi
 */

#ifndef PASSIVENEOHOOKEAN_HPP_
#define PASSIVENEOHOOKEAN_HPP_


#include <lifev/em/solver/mechanics/materials/EMPassiveMaterialType.hpp>

//#include <lifev/em/solver/mechanics/materials/EMMaterialFunctions.hpp>



//class EMPassiveNeoHookean : public virtual EMPassiveMaterial
//{
//public:
//  EMPassiveNeoHookean();
//  virtual ~EMPassiveNeoHookean() {}
//};

namespace LifeV
{

template<typename Mesh>
class PassiveNeoHookean : public virtual EMPassiveMaterialType<Mesh>
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
	std::cout << "Number of Functions in List: " << this->M_materialFunctionList.size();
    this -> M_materialFunctionList[0].reset (new MaterialFunctions::Volumetric<Mesh> (3500000.0)  );
    this -> M_materialFunctionList[1].reset (new MaterialFunctions::dVolumetric<Mesh> (3500000.0) );
    this -> M_materialFunctionList[2].reset (new MaterialFunctions::NeoHookean<Mesh> ()  );
}

template <typename MeshType>
inline EMPassiveMaterialType<MeshType>* createPassiveNeoHookean()
{
    return new PassiveNeoHookean<MeshType>();
}
namespace
{
static bool registerEM_passiveNH = EMPassiveMaterialType<LifeV::RegionMesh<LinearTetra> >::EMPassiveMaterialFactory::instance().registerProduct ("PNH", &createPassiveNeoHookean<LifeV::RegionMesh<LinearTetra> > );
}

}//LifeV

#endif /* EMPASSIVENEOHOOKEAN_HPP_ */
