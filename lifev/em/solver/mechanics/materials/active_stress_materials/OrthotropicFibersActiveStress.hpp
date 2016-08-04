/*
 * EMPassiveNeoHookean.hpp
 *
 *  Created on: 29/apr/2014
 *      Author: srossi
 */

#ifndef OrthotropicFibersActiveStress_HPP_
#define OrthotropicFibersActiveStress_HPP_


#include <lifev/em/solver/mechanics/materials/EMActiveMaterialType.hpp>


namespace LifeV
{

template<typename Mesh>
class OrthotropicFibersActiveStress : public virtual EMActiveMaterialType<Mesh>
{
public:
    typedef EMMaterialType<Mesh> super;

    OrthotropicFibersActiveStress();
    virtual ~OrthotropicFibersActiveStress() {}
};

template<typename Mesh>
OrthotropicFibersActiveStress<Mesh>::OrthotropicFibersActiveStress() :
    super ("Orthotropic Active Stress", 3)
{
    this -> M_materialFunctionList[0].reset (new MaterialFunctions::OrthotropicActiveStress<Mesh> (47900.0)  );
    this -> M_materialFunctionList[1].reset (new MaterialFunctions::OrthotropicActiveStress<Mesh> (47900.0, "Sheets")  );
    this -> M_materialFunctionList[2].reset (new MaterialFunctions::OrthotropicActiveStress<Mesh> (47900.0, "Normal")  );
}


template <typename MeshType>
inline EMActiveMaterialType<MeshType>* createOrthotropicFibersActiveStress()
{
    return new OrthotropicFibersActiveStress<MeshType>();
}
namespace
{
static bool registerEM_OrthotropicFibersActiveStress = EMActiveMaterialType<LifeV::RegionMesh<LinearTetra> >::EMActiveMaterialFactory::instance().registerProduct ("OAS", &createOrthotropicFibersActiveStress<LifeV::RegionMesh<LinearTetra> > );
}

}//LifeV

#endif /* EMPASSIVENEOHOOKEAN_HPP_ */
