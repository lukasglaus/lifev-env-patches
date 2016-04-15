/*
 *  PassiveMooneyRivlinCompact.hpp
 *
 *  Created on: 15/mar/2015
 *      Author: srossi
 */

#ifndef PASSIVEMOONEYRIVLINCOMPACT_HPP_
#define PASSIVEMOONEYRIVLINCOMPACT_HPP_


#include <lifev/em/solver/mechanics/materials/EMPassiveMaterialType.hpp>

//#include <lifev/em/solver/mechanics/materials/EMMaterialFunctions.hpp>


namespace LifeV
{

template<typename Mesh>
class PassiveMooneyRivlinCompact : public virtual EMPassiveMaterialType<Mesh>
{
public:
    typedef EMMaterialType<Mesh> super;

    PassiveMooneyRivlinCompact();
    virtual ~PassiveMooneyRivlinCompact() {}

};

template<typename Mesh>
PassiveMooneyRivlinCompact<Mesh>::PassiveMooneyRivlinCompact() :
    super("Passive Mooney-Rivlin Compact", 1)
{
	std::cout << "Number of Functions in List: " << this->M_materialFunctionList.size();
    this -> M_materialFunctionList[0].reset (new MaterialFunctions::MooneyRivlinCompact<Mesh> (2000, 2000, 4960, 4960) );
}

template <typename MeshType>
inline EMPassiveMaterialType<MeshType>* createPassiveMooneyRivlinCompact()
{
    return new PassiveMooneyRivlinCompact<MeshType>();
}
namespace
{
static bool registerEM_passiveMRC = EMPassiveMaterialType<LifeV::RegionMesh<LinearTetra> >::EMPassiveMaterialFactory::instance().registerProduct ("PMRC", &createPassiveMooneyRivlinCompact<LifeV::RegionMesh<LinearTetra> > );
}

}//LifeV

#endif /* EMPASSIVEMOONEYRIVLINCOMPACT_HPP_ */
