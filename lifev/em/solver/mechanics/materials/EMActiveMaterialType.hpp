/*
 * EMMaterialData.hpp
 *
 *  Created on: 29/apr/2014
 *      Author: srossi
 */

#ifndef EMACTIVESTRESSTYPE_HPP_
#define EMACTIVESTRESSTYPE_HPP_

#include <lifev/em/solver/mechanics/materials/EMMaterialType.hpp>

namespace LifeV
{

template<typename Mesh>
class EMActiveMaterialType : public virtual EMMaterialType<Mesh>
{
public:

    typedef Mesh mesh_Type;

    typedef typename boost::shared_ptr<MaterialFunctions::EMMaterialFunctions<Mesh> >    materialFunctionsPtr_Type;
    typedef std::vector<materialFunctionsPtr_Type>  vectorMaterialsPtr_Type;

    typedef ETFESpace< mesh_Type, MapEpetra, 3, 1 >                        scalarETFESpace_Type;
    typedef boost::shared_ptr<ETFESpace< mesh_Type, MapEpetra, 3, 1 > >    scalarETFESpacePtr_Type;


    //template <class Mesh>
    using ETFESpacePtr_Type = boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > >;
    //template <class Mesh>
    using FESpacePtr_Type = boost::shared_ptr< FESpace< Mesh, MapEpetra >  >;

    typedef VectorEpetra           vector_Type;
    typedef boost::shared_ptr<vector_Type>         vectorPtr_Type;

    typedef MatrixEpetra<Real>           matrix_Type;
    typedef boost::shared_ptr<matrix_Type>         matrixPtr_Type;

    typedef EMMaterialType<Mesh> super;


    EMActiveMaterialType (std::string materialName = "None", UInt n = 0);
    virtual ~EMActiveMaterialType() {}


    void
    computeJacobian ( const vector_Type&       disp,
                      const vector_Type&       activation,
                      scalarETFESpacePtr_Type  aETFESpace,
                      ETFESpacePtr_Type        dispETFESpace,
                      FESpacePtr_Type          dispFESpace,
                      matrixPtr_Type           jacobianPtr);
    //



    void
    computeResidual ( const vector_Type&       disp,
                      const vector_Type&       activation,
                      scalarETFESpacePtr_Type  aETFESpace,
                      ETFESpacePtr_Type        dispETFESpace,
                      FESpacePtr_Type          dispFESpace,
                      vectorPtr_Type           residualVectorPtr);

//    inline void showMe()
//    {
//        std::cout << "Material Type: " << M_materialName << "\n";
//    }

protected:
    vectorPtr_Type          M_activeStressPtr;

};


template<typename Mesh>
EMActiveMaterialType<Mesh>::EMActiveMaterialType (std::string materialName, UInt n ) :
	super(materialName, n)
{

}



template<typename Mesh>
void
EMActiveMaterialType<Mesh>::computeJacobian ( const vector_Type&       disp,
                                              const vector_Type&       activation,
                                              scalarETFESpacePtr_Type  aETFESpace,
                                              ETFESpacePtr_Type        dispETFESpace,
                                              FESpacePtr_Type          dispFESpace,
                                              matrixPtr_Type           jacobianPtr)
{
    int n = this->M_materialFunctionList.size();
    for (int j (0); j < n; j++)
    {
    	this->M_materialFunctionList[j]->computeJacobian (disp, activation, aETFESpace, dispETFESpace, dispFESpace, jacobianPtr);
    }
}

template <typename Mesh>
void
EMActiveMaterialType<Mesh>::computeResidual ( const vector_Type&       disp,
                                              const vector_Type&       activation,
                                              scalarETFESpacePtr_Type  aETFESpace,
                                              ETFESpacePtr_Type        dispETFESpace,
                                              FESpacePtr_Type          dispFESpace,
                                              vectorPtr_Type           residualVectorPtr)
{
    int n = this->M_materialFunctionList.size();
    for (int j (0); j < n; j++)
    {
    	this->M_materialFunctionList[j]->computeResidual (disp, activation, aETFESpace, dispETFESpace, dispFESpace, residualVectorPtr);
    }
}



}//LifeV
#endif /* EMMATERIALDATA_HPP_ */
