/*
 * EMPassiveNeoHookean.hpp
 *
 *  Created on: 29/apr/2014
 *      Author: srossi
 */

#ifndef PASSIVEHOLZAPFELOGDEN_HPP_
#define PASSIVEHOLZAPFELOGDEN_HPP_


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
class PassiveHolzapfelOgden : public virtual EMPassiveMaterialType<Mesh>
{
public:
    typedef EMMaterialType<Mesh> super;

    PassiveHolzapfelOgden();
    virtual ~PassiveHolzapfelOgden() {}
};

template<typename Mesh>
PassiveHolzapfelOgden<Mesh>::PassiveHolzapfelOgden() :
    super ("Passive Holzapfel Ogden", 6)
{
    //this -> M_materialFunctionList[0].reset ( new MaterialFunctions::Volumetric<Mesh> (3500000.0)  );
    //this -> M_materialFunctionList[1].reset ( new MaterialFunctions::dVolumetric<Mesh> (3500000.0) );
    //this -> M_materialFunctionList[2].reset ( new MaterialFunctions::IsotropicExponential<Mesh> (3330.0, 9.242)  );
    //this -> M_materialFunctionList[3].reset ( new MaterialFunctions::dIsotropicExponential<Mesh> (3330.0, 9.242)  );
    this -> M_materialFunctionList[0].reset ( new MaterialFunctions::AnisotropicExponential<Mesh> (185350, 15.972)  );
    this -> M_materialFunctionList[1].reset ( new MaterialFunctions::dAnisotropicExponential<Mesh> (185350, 15.972)  );
    this -> M_materialFunctionList[2].reset ( new MaterialFunctions::AnisotropicExponential<Mesh> ( 25640.0, 10.446, MaterialFunctions::AnisotropicExponential<Mesh>::Sheets )  );
    this -> M_materialFunctionList[3].reset ( new MaterialFunctions::dAnisotropicExponential<Mesh> ( 25640.0, 10.446, MaterialFunctions::AnisotropicExponential<Mesh>::Sheets )  );
    this -> M_materialFunctionList[4].reset ( new MaterialFunctions::ShearExponential<Mesh> (4170.0, 11.602)  );
    this -> M_materialFunctionList[5].reset ( new MaterialFunctions::dShearExponential<Mesh> (4170.0, 11.602)  );

    //  this -> M_materialFunctionList[6].reset(new MaterialFunctions::AnisotropicExponential<Mesh>(25640, 10.446, false)  );
    //  this -> M_materialFunctionList[7].reset(new MaterialFunctions::dAnisotropicExponential<Mesh>(25640, 10.446, false)  );
    //  this -> M_materialFunctionList[8].reset( new MaterialFunctions::ShearExponential<Mesh>(4170, 11.602)  );
    //  this -> M_materialFunctionList[9].reset( new MaterialFunctions::dShearExponential<Mesh>(4170, 11.602)  );
}


template <typename MeshType>
inline EMPassiveMaterialType<MeshType>* createPassiveHolzapfelOgden()
{
    return new PassiveHolzapfelOgden<MeshType>();
}
namespace
{
static bool registerEM_passiveHO = EMPassiveMaterialType<LifeV::RegionMesh<LinearTetra> >::EMPassiveMaterialFactory::instance().registerProduct ("PHO", &createPassiveHolzapfelOgden<LifeV::RegionMesh<LinearTetra> > );
}

}//LifeV

#endif /* EMPASSIVENEOHOOKEAN_HPP_ */
