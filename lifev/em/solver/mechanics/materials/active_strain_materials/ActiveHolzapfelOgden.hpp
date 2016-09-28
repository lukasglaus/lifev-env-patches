/*
 * ActiveHolzapfelOgden.hpp
 *
 *  Created on: 02/Mar/2015
 *      Author: srossi, Luca Barbarotta
 */

#ifndef ACTIVE_HOLZAPFEL_OGDEN_HPP_
#define ACTIVE_HOLZAPFEL_OGDEN_HPP_


#include <lifev/em/solver/mechanics/materials/EMActiveMaterialType.hpp>

//#include <lifev/em/solver/mechanics/materials/EMMaterialFunctions.hpp>

namespace LifeV
{

    template<typename Mesh>
    class ActiveHolzapfelOgden : public virtual EMActiveMaterialType<Mesh>
    {
    public:
        typedef EMMaterialType<Mesh> super;

        ActiveHolzapfelOgden();
        virtual ~ActiveHolzapfelOgden() {}
    };

    template<typename Mesh>
    ActiveHolzapfelOgden<Mesh>::ActiveHolzapfelOgden() :
        super ("Active Strain Holzapfel Ogden", 8)
    {
        //this -> M_materialFunctionList[8].reset ( new MaterialFunctions::Volumetric<Mesh> (3500000.0)  );
        //this -> M_materialFunctionList[9].reset ( new MaterialFunctions::dVolumetric<Mesh> (3500000.0) );
        this -> M_materialFunctionList[0].reset ( new MaterialFunctions::ActiveStrainIsotropicExponential<Mesh> (3330.0, 9.242)  );
        this -> M_materialFunctionList[1].reset ( new MaterialFunctions::dActiveStrainIsotropicExponential<Mesh> (3330.0, 9.242)  );
        this -> M_materialFunctionList[2].reset ( new MaterialFunctions::ActiveStrainAnisotropicExponential<Mesh> (185350, 15.972)  );
        this -> M_materialFunctionList[3].reset ( new MaterialFunctions::dActiveStrainAnisotropicExponential<Mesh> (185350, 15.972)  );
        this -> M_materialFunctionList[4].reset (
            new MaterialFunctions::ActiveStrainAnisotropicExponential<Mesh> ( 25640.0,
                                                                              10.446,
                                                                              MaterialFunctions::ActiveStrainAnisotropicExponential<Mesh>::Sheets ) );
        this -> M_materialFunctionList[5].reset (
            new MaterialFunctions::dActiveStrainAnisotropicExponential<Mesh> ( 25640.0,
                                                                               10.446,
                                                                               MaterialFunctions::ActiveStrainAnisotropicExponential<Mesh>::Sheets )  );
        this -> M_materialFunctionList[6].reset ( new MaterialFunctions::ActiveStrainShearExponential<Mesh> (4170.0, 11.602)  );
        this -> M_materialFunctionList[7].reset ( new MaterialFunctions::dActiveStrainShearExponential<Mesh> (4170.0, 11.602)  );

    }


    template <typename MeshType>
    inline EMActiveMaterialType<MeshType>* createActiveHolzapfelOgden()
    {
        return new ActiveHolzapfelOgden<MeshType>();
    }
    namespace
    {
        static bool registerEM_activeHO = EMActiveMaterialType<LifeV::RegionMesh<LinearTetra> >::EMActiveMaterialFactory::instance().registerProduct ("AHO", &createActiveHolzapfelOgden<LifeV::RegionMesh<LinearTetra> > );
    }

}//LifeV

#endif /* ACTIVE_HOLZAPFEL_OGDEN_HPP_ */
