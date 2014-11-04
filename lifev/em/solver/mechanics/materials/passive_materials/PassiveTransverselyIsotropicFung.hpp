/*
 * EMPassiveNeoHookean.hpp
 *
 *  Created on: 29/apr/2014
 *      Author: srossi
 */

#ifndef PASSIVETRANSVERSELYISOTROPICFUNG_HPP_
#define PASSIVETRANSVERSELYISOTROPICFUNG_HPP_


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
class PassiveTransverselyIsotropicFung : public virtual EMPassiveMaterialType<Mesh>
{
public:
    typedef EMMaterialType<Mesh> super;

    PassiveTransverselyIsotropicFung();
    virtual ~PassiveTransverselyIsotropicFung() {}
};

template<typename Mesh>
PassiveTransverselyIsotropicFung<Mesh>::PassiveTransverselyIsotropicFung() :
    super ("Passive Transversely Isotropic Fung", 3)
{
    this -> M_materialFunctionList[0].reset (new MaterialFunctions::Volumetric<Mesh> (3500.0)  );
    this -> M_materialFunctionList[1].reset (new MaterialFunctions::dVolumetric<Mesh> (3500.0) );
//    Real C = 8760.0;
//    Real bf = 18.48;
//    Real bt =  3.58;
//    Real bfs = 1.627;
    Real C = 20000.0;
    Real bf = 8.00;
    Real bt =  3.00;
    Real bfs = 4.00;
    this -> M_materialFunctionList[2].reset (new MaterialFunctions::OrthotropicFung<Mesh> (C, bf, bt, bt, bfs, bfs, bt) );
}


template <typename MeshType>
inline EMPassiveMaterialType<MeshType>* createPassiveTransverselyIsotropicFung()
{
    return new PassiveTransverselyIsotropicFung<MeshType>();
}
namespace
{
static bool registerEM_passiveTIF = EMPassiveMaterialType<LifeV::RegionMesh<LinearTetra> >::EMPassiveMaterialFactory::instance().registerProduct ("PTIF", &createPassiveTransverselyIsotropicFung<LifeV::RegionMesh<LinearTetra> > );
}

}//LifeV

#endif /* EMPASSIVENEOHOOKEAN_HPP_ */
