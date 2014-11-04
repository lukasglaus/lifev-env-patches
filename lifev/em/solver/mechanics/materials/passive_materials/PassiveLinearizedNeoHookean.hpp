/*
 * EMPassiveNeoHookean.hpp
 *
 *  Created on: 29/apr/2014
 *      Author: srossi
 */

#ifndef PASSIVELINEARIZEDNEOHOOKEAN_HPP_
#define PASSIVELINEARIZEDNEOHOOKEAN_HPP_


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
class PassiveLinearizedNeoHookean : public virtual EMPassiveMaterialType<Mesh>
{
public:
    typedef EMMaterialType<Mesh> super;

    PassiveLinearizedNeoHookean();
    virtual ~PassiveLinearizedNeoHookean() {}
};

template<typename Mesh>
PassiveLinearizedNeoHookean<Mesh>::PassiveLinearizedNeoHookean() :
    super ("Passive Linearized Neo-Hookean", 2)
{
    this -> M_materialFunctionList[0].reset (new MaterialFunctions::LinearizedNeoHookeanVolumetric<Mesh>()  );
    this -> M_materialFunctionList[1].reset (new MaterialFunctions::LinearizedNeoHookeanDeviatoric<Mesh>() );
}


template <typename MeshType>
inline EMPassiveMaterialType<MeshType>* createPassiveLinearizedNeoHookean()
{
    return new PassiveLinearizedNeoHookean<MeshType>();
}
namespace
{
static bool registerEM_passiveNHL = EMPassiveMaterialType<LifeV::RegionMesh<LinearTetra> >::EMPassiveMaterialFactory::instance().registerProduct ("PNHL", &createPassiveLinearizedNeoHookean<LifeV::RegionMesh<LinearTetra> > );
}

}//LifeV

#endif /* EMPASSIVENEOHOOKEAN_HPP_ */
