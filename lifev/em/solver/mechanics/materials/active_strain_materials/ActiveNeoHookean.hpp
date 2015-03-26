/*
 * ActiveNeoHookean.hpp
 *
 *  Created on: 29/apr/2014
 *      Author: srossi
 */

#ifndef ACTIVENEOHOOKEAN_HPP_
#define ACTIVENEOHOOKEAN_HPP_


#include <lifev/em/solver/mechanics/materials/EMActiveMaterialType.hpp>


namespace LifeV
{

template<typename Mesh>
class ActiveNeoHookean : public virtual EMActiveMaterialType<Mesh>
{
public:
    typedef EMMaterialType<Mesh> super;

    ActiveNeoHookean();
    virtual ~ActiveNeoHookean() {}
};

template<typename Mesh>
ActiveNeoHookean<Mesh>::ActiveNeoHookean() :
    super ("Active Strain NeoHookean", 1)
{
    this -> M_materialFunctionList[0].reset (new MaterialFunctions::ActiveStrainNeoHookean<Mesh> (47900.0)  );
}


template <typename MeshType>
inline EMActiveMaterialType<MeshType>* createActiveNeoHookean()
{
    return new ActiveNeoHookean<MeshType>();
}
namespace
{
static bool registerEM_ActiveNeoHookean = EMActiveMaterialType<LifeV::RegionMesh<LinearTetra> >::EMActiveMaterialFactory::instance().registerProduct ("ANH", &createActiveNeoHookean<LifeV::RegionMesh<LinearTetra> > );
}

}//LifeV

#endif /* ACTIVENEOHOOKEAN_HPP_ */
