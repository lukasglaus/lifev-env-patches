/*
 * ActiveHolzapfelOgden.hpp
 *
 *  Created on: 02/Mar/2015
 *      Author: srossi, Luca Barbarotta
 */

#ifndef ACTIVE_HOLZAPFEL_OGDEN_ISOTROPIC_SPLITTING_HPP_
#define ACTIVE_HOLZAPFEL_OGDEN_ISOTROPIC_SPLITTING_HPP_


#include <lifev/em/solver/mechanics/materials/EMActiveMaterialType.hpp>

//#include <lifev/em/solver/mechanics/materials/EMMaterialFunctions.hpp>

namespace LifeV
{

    template<typename Mesh>
    class ActiveHolzapfelOgdenIsoNISplit : public virtual EMActiveMaterialType<Mesh>
    {
    public:
        typedef EMMaterialType<Mesh> super;

        ActiveHolzapfelOgdenIsoNISplit();
        virtual ~ActiveHolzapfelOgdenIsoNISplit() {}
    };

    template<typename Mesh>
    ActiveHolzapfelOgdenIsoNISplit<Mesh>::ActiveHolzapfelOgdenIsoNISplit() :
        super ("Active Strain Holzapfel Ogden Isotropic Nearly Incompressible Splitting", 8)
    {
        this -> M_materialFunctionList[0].reset ( new MaterialFunctions::ActiveStrainIsotropicExponential<Mesh> (3330.0, 9.242)  );
        this -> M_materialFunctionList[1].reset ( new MaterialFunctions::dActiveStrainIsotropicExponential<Mesh> (3330.0, 9.242)  );
        this -> M_materialFunctionList[2].reset (
            new MaterialFunctions::ActiveStrainAnisotropicExponential<Mesh> (185350,
                                                                             15.972,
                                                                             MaterialFunctions::ActiveStrainAnisotropicExponential<Mesh>::Fibers, false) );
        this -> M_materialFunctionList[3].reset (
            new MaterialFunctions::dActiveStrainAnisotropicExponential<Mesh> (185350,
                                                                              15.972,
                                                                              MaterialFunctions::ActiveStrainAnisotropicExponential<Mesh>::Fibers, false) );
        this -> M_materialFunctionList[4].reset (
            new MaterialFunctions::ActiveStrainAnisotropicExponential<Mesh> ( 25640.0,
                                                                              10.446,
                                                                              MaterialFunctions::ActiveStrainAnisotropicExponential<Mesh>::Sheets, false ) );
        this -> M_materialFunctionList[5].reset (
            new MaterialFunctions::dActiveStrainAnisotropicExponential<Mesh> ( 25640.0,
                                                                               10.446,
                                                                               MaterialFunctions::ActiveStrainAnisotropicExponential<Mesh>::Sheets, false )  );
        this -> M_materialFunctionList[6].reset ( new MaterialFunctions::ActiveStrainShearExponential<Mesh> (4170.0, 11.602, false)  );
        this -> M_materialFunctionList[7].reset ( new MaterialFunctions::dActiveStrainShearExponential<Mesh> (4170.0, 11.602, false)  );

    }


    template <typename MeshType>
    inline EMActiveMaterialType<MeshType>* createActiveHolzapfelOgdenIsoNISplit()
    {
        return new ActiveHolzapfelOgdenIsoNISplit<MeshType>();
    }
    namespace
    {
        static bool registerEM_activeHOINIS = EMActiveMaterialType<LifeV::RegionMesh<LinearTetra> >::EMActiveMaterialFactory::instance().registerProduct ("AHO-INIS", &createActiveHolzapfelOgdenIsoNISplit<LifeV::RegionMesh<LinearTetra> > );
    }

}//LifeV

#endif /* ACTIVE_HOLZAPFEL_OGDEN_ISOTROPIC_SPLITTING_HPP_ */
