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

    typedef FactorySingleton<Factory<EMActiveMaterialType<Mesh>, std::string> >  EMActiveMaterialFactory;


    //template <class Mesh>
    typedef boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > > ETFESpacePtr_Type;
    //template <class Mesh>
    typedef boost::shared_ptr< FESpace< Mesh, MapEpetra >  > FESpacePtr_Type;

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
                      matrixPtr_Type           jacobianPtr) {}

    void
    computeJacobian ( const vector_Type& disp,
                      ETFESpacePtr_Type  dispETFESpace,
                      const vector_Type& fibers,
                      const vector_Type& sheets,
                      const vector_Type& activation,
                      scalarETFESpacePtr_Type  activationETFESpace,
                      matrixPtr_Type     jacobianPtr);
    //



    void
    computeResidual ( const vector_Type&       disp,
                      const vector_Type&       activation,
                      scalarETFESpacePtr_Type  aETFESpace,
                      ETFESpacePtr_Type        dispETFESpace,
                      FESpacePtr_Type          dispFESpace,
                      vectorPtr_Type           residualVectorPtr) {}

    void
    computeResidual ( const vector_Type& disp,
                      ETFESpacePtr_Type  dispETFESpace,
                      const vector_Type& fibers,
                      const vector_Type& sheets,
                      const vector_Type& activation,
                      scalarETFESpacePtr_Type  activationETFESpace,
                      vectorPtr_Type     residualVectorPtr);

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
EMActiveMaterialType<Mesh>::computeJacobian ( const vector_Type& disp,
												ETFESpacePtr_Type  dispETFESpace,
												const vector_Type& fibers,
												const vector_Type& sheets,
												const vector_Type& activation,
												scalarETFESpacePtr_Type  activationETFESpace,
												matrixPtr_Type     jacobianPtr)
{
    int n = this->M_materialFunctionList.size();
    for (int j (0); j < n; j++)
    {
    	this->M_materialFunctionList[j]->computeJacobian ( disp,
    			                                           dispETFESpace,
    			                                           fibers,
    			                                           sheets,
    			                                           activation,
    			                                           activationETFESpace,
    			                                           jacobianPtr);
    }
}

template <typename Mesh>
void
EMActiveMaterialType<Mesh>::computeResidual  (  const vector_Type& disp,
												ETFESpacePtr_Type  dispETFESpace,
												const vector_Type& fibers,
												const vector_Type& sheets,
												const vector_Type& activation,
												scalarETFESpacePtr_Type  activationETFESpace,
												vectorPtr_Type     residualVectorPtr)
{
    int n = this->M_materialFunctionList.size();
    for (int j (0); j < n; j++)
    {
    	this->M_materialFunctionList[j]->computeResidual (disp,
    													  dispETFESpace,
    													  fibers,
    													  sheets,
    													  activation,
    													  activationETFESpace,
    													  residualVectorPtr );
    }
}



}//LifeV
#endif /* EMMATERIALDATA_HPP_ */
