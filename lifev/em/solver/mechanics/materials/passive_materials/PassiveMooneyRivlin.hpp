/*
 *  PassiveMooneyRivlin.hpp
 *
 *  Created on: 15/mar/2015
 *      Author: srossi
 */

#ifndef PASSIVEMOONEYRIVLIN_HPP_
#define PASSIVEMOONEYRIVLIN_HPP_


#include <lifev/em/solver/mechanics/materials/EMPassiveMaterialType.hpp>

//#include <lifev/em/solver/mechanics/materials/EMMaterialFunctions.hpp>


namespace LifeV
{

template<typename Mesh>
class PassiveMooneyRivlin : public virtual EMPassiveMaterialType<Mesh>
{
public:
    typedef EMMaterialType<Mesh> super;

    PassiveMooneyRivlin();
    virtual ~PassiveMooneyRivlin() {}

};

template<typename Mesh>
PassiveMooneyRivlin<Mesh>::PassiveMooneyRivlin() :
    super("Passive Mooney-Rivlin", 4)
{
	std::cout << "Number of Functions in List: " << this->M_materialFunctionList.size();
    this -> M_materialFunctionList[0].reset (new MaterialFunctions::Volumetric<Mesh> (3500000.0) );
    this -> M_materialFunctionList[1].reset (new MaterialFunctions::dVolumetric<Mesh> (3500000.0)  );
    this -> M_materialFunctionList[2].reset (new MaterialFunctions::NeoHookean<Mesh> () );
    this -> M_materialFunctionList[3].reset (new MaterialFunctions::MooneyRivlin<Mesh> () );
}

template <typename MeshType>
inline EMPassiveMaterialType<MeshType>* createPassiveMooneyRivlin()
{
    return new PassiveMooneyRivlin<MeshType>();
}
namespace
{
static bool registerEM_passiveMR = EMPassiveMaterialType<LifeV::RegionMesh<LinearTetra> >::EMPassiveMaterialFactory::instance().registerProduct ("PMR", &createPassiveMooneyRivlin<LifeV::RegionMesh<LinearTetra> > );
}

}//LifeV

#endif /* EMPASSIVEMOONEYRIVLIN_HPP_ */
