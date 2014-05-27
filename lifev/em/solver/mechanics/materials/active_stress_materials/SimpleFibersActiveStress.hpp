/*
 * EMPassiveNeoHookean.hpp
 *
 *  Created on: 29/apr/2014
 *      Author: srossi
 */

#ifndef SIMPLEFIBERSACTIVESTRESS_HPP_
#define SIMPLEFIBERSACTIVESTRESS_HPP_


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
class SimpleFibersActiveStress : public virtual EMMaterialType<Mesh>
{
public:
    typedef EMMaterialType<Mesh> super;

    SimpleFibersActiveStress();
    virtual ~SimpleFibersActiveStress() {}
};

template<typename Mesh>
SimpleFibersActiveStress<Mesh>::SimpleFibersActiveStress() :
    super ("Simple Fibers Active Stress", 1)
{
    this -> M_materialFunctionList[0].reset (new MaterialFunctions::SimpleActiveStress<Mesh> (47900.0)  );
}


template <typename MeshType>
inline EMMaterialType<MeshType>* createSimpleFibersActiveStress()
{
    return new SimpleFibersActiveStress<MeshType>();
}
namespace
{
static bool registerEM_SimpleFibersActiveStress = EMMaterialType<LifeV::RegionMesh<LinearTetra> >::EMMaterialFactory::instance().registerProduct ("SimpleActiveStress", &createSimpleFibersActiveStress<LifeV::RegionMesh<LinearTetra> > );
}

}//LifeV

#endif /* EMPASSIVENEOHOOKEAN_HPP_ */
