/*
 * EMPassiveNeoHookean.hpp
 *
 *  Created on: 29/apr/2014
 *      Author: srossi
 */

#ifndef PASSIVEISOTROPICFUNG_HPP_
#define PASSIVEISOTROPICFUNG_HPP_


#include <lifev/em/solver/mechanics/materials/EMMaterialType.hpp>


namespace LifeV
{

template<typename Mesh>
class PassiveIsotropicFung : public virtual EMMaterialType<Mesh>
{
public:
    typedef EMMaterialType<Mesh> super;

    PassiveIsotropicFung();
    virtual ~PassiveIsotropicFung() {}
};

template<typename Mesh>
PassiveIsotropicFung<Mesh>::PassiveIsotropicFung() :
    super ("Passive Isotropic Fung", 3)
{
    this -> M_materialFunctionList[0].reset (new MaterialFunctions::Volumetric<Mesh> (3500000.0)  );
    this -> M_materialFunctionList[1].reset (new MaterialFunctions::dVolumetric<Mesh>(3500000.0) );
    Real C = 100000.0;
    Real b = 1.0;//0.0
    Real bt =  1.0;
    Real bfs = 1.0;
    this -> M_materialFunctionList[2].reset (new MaterialFunctions::IsotropicFung<Mesh> (C, b) );
//    this -> M_materialFunctionList[2].reset (new MaterialFunctions::OrthotropicFung<Mesh> (C, b, bt, bfs) );
}


template <typename MeshType>
inline EMMaterialType<MeshType>* createPassiveIsotropicFung()
{
    return new PassiveIsotropicFung<MeshType>();
}
namespace
{
static bool registerEM_passiveIF = EMMaterialType<LifeV::RegionMesh<LinearTetra> >::EMMaterialFactory::instance().registerProduct ("PIF", &createPassiveIsotropicFung<LifeV::RegionMesh<LinearTetra> > );
}

}//LifeV

#endif /* EMPASSIVENEOHOOKEAN_HPP_ */
