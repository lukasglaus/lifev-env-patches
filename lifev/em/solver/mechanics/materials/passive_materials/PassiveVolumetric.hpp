/*
 * EMPassiveVolumetric.hpp
 *
 *  Created on: 29/apr/2014
 *      Author: srossi
 */

#ifndef PASSIVEVOLUMETRIC_HPP_
#define PASSIVEVOLUMETRIC_HPP_


#include <lifev/em/solver/mechanics/materials/EMPassiveMaterialType.hpp>

//#include <lifev/em/solver/mechanics/materials/EMMaterialFunctions.hpp>



//class EMPassiveVolumetric : public virtual EMPassiveMaterial
//{
//public:
//  EMPassiveVolumetric();
//  virtual ~EMPassiveVolumetric() {}
//};

namespace LifeV
{

template<typename Mesh>
class PassiveVolumetric : public virtual EMPassiveMaterialType<Mesh>
{
public:
    typedef EMMaterialType<Mesh> super;

    PassiveVolumetric();
    virtual ~PassiveVolumetric() {}

};

template<typename Mesh>
PassiveVolumetric<Mesh>::PassiveVolumetric() :
    super("Passive Volumetric", 2)
{
	std::cout << "Number of Functions in List: " << this->M_materialFunctionList.size();
    this -> M_materialFunctionList[0].reset (new MaterialFunctions::Volumetric<Mesh> (3500000.0)  );
    this -> M_materialFunctionList[1].reset (new MaterialFunctions::dVolumetric<Mesh> (3500000.0) );
}

template <typename MeshType>
inline EMPassiveMaterialType<MeshType>* createPassiveVolumetric()
{
    return new PassiveVolumetric<MeshType>();
}
namespace
{
static bool registerEM_passiveVol = EMPassiveMaterialType<LifeV::RegionMesh<LinearTetra> >::EMPassiveMaterialFactory::instance().registerProduct ("PV", &createPassiveVolumetric<LifeV::RegionMesh<LinearTetra> > );
}

}//LifeV

#endif /* EMPassiveVolumetric_HPP_ */
