/*
 * ActiveIsotropicExponential.hpp
 *
 *  Created on: 02/Mar/2015
 *      Author: Luca Barbarotta
 */

#ifndef ActiveIsotropicExponential_HPP_
#define ActiveIsotropicExponential_HPP_


#include <lifev/em/solver/mechanics/materials/EMActiveMaterialType.hpp>


namespace LifeV
{

template<typename Mesh>
class ActiveIsotropicExponential : public virtual EMActiveMaterialType<Mesh>
{
public:
    typedef EMMaterialType<Mesh> super;

    ActiveIsotropicExponential();
    virtual ~ActiveIsotropicExponential() {}
};

template<typename Mesh>
ActiveIsotropicExponential<Mesh>::ActiveIsotropicExponential() :
    super ("Active Strain Isotropic Exponential", 2)
{
    this -> M_materialFunctionList[0].reset (new MaterialFunctions::ActiveStrainIsotropicExponential<Mesh> (3330.0, 9.242)  );
    this -> M_materialFunctionList[1].reset (new MaterialFunctions::dActiveStrainIsotropicExponential<Mesh> (3330.0, 9.242)  );
}


template <typename MeshType>
inline EMActiveMaterialType<MeshType>* createActiveIsotropicExponential()
{
    return new ActiveIsotropicExponential<MeshType>();
}
namespace
{
static bool registerEM_ActiveIsotropicExponential = EMActiveMaterialType<LifeV::RegionMesh<LinearTetra> >::EMActiveMaterialFactory::instance().registerProduct ("AIE", &createActiveIsotropicExponential<LifeV::RegionMesh<LinearTetra> > );
}

}//LifeV

#endif /* ACTIVEISOTROPICEXPONENTIAL_HPP_ */
