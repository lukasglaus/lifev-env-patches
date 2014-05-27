/*
 * EMPassiveNeoHookean.hpp
 *
 *  Created on: 29/apr/2014
 *      Author: srossi
 */

#ifndef PASSIVEORTHOTROPICFUNG_HPP_
#define PASSIVEORTHOTROPICFUNG_HPP_


#include <lifev/em/solver/mechanics/materials/EMMaterialType.hpp>

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
class PassiveOrthotropicFung : public virtual EMMaterialType<Mesh>
{
public:
    typedef EMMaterialType<Mesh> super;

    PassiveOrthotropicFung();
    virtual ~PassiveOrthotropicFung() {}
};

template<typename Mesh>
PassiveOrthotropicFung<Mesh>::PassiveOrthotropicFung() :
    super ("Passive Orthotropic Fung", 3)
{
    this -> M_materialFunctionList[0].reset (new MaterialFunctions::Volumetric<Mesh> (350000.0)  );
    this -> M_materialFunctionList[1].reset (new MaterialFunctions::dVolumetric<Mesh> (350000.0) );
    this -> M_materialFunctionList[2].reset (new MaterialFunctions::OrthotropicFung<Mesh>() );
}


template <typename MeshType>
inline EMMaterialType<MeshType>* createPassiveOrthotropicFung()
{
    return new PassiveOrthotropicFung<MeshType>();
}
namespace
{
static bool registerEM_passiveOF = EMMaterialType<LifeV::RegionMesh<LinearTetra> >::EMMaterialFactory::instance().registerProduct ("POF", &createPassiveOrthotropicFung<LifeV::RegionMesh<LinearTetra> > );
}

}//LifeV

#endif /* EMPASSIVENEOHOOKEAN_HPP_ */
