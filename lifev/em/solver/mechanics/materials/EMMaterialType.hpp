/*
 * EMMaterialData.hpp
 *
 *  Created on: 29/apr/2014
 *      Author: srossi
 */

#ifndef EMMATERIALTYPE_HPP_
#define EMMATERIALTYPE_HPP_

#include <lifev/core/util/Factory.hpp>
#include <lifev/core/util/FactorySingleton.hpp>


#include <lifev/em/solver/mechanics/materials/functions/FunctionsList.hpp>


namespace LifeV
{
//We forward declare the following class
// as we do not need its implementation
class StructuralConstitutiveLawData;
class EMData;

///
template<typename Mesh>
class EMMaterialType
{
public:

    typedef typename boost::shared_ptr<MaterialFunctions::EMMaterialFunctions<Mesh> >    materialFunctionsPtr_Type;
    typedef std::vector<materialFunctionsPtr_Type>  vectorMaterialsPtr_Type;

//    typedef FactorySingleton<Factory<EMMaterialType<Mesh>, std::string> >  EMMaterialFactory;


    typedef Mesh mesh_Type;
    //template <class Mesh>
    typedef boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > > ETFESpacePtr_Type;
    typedef boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 1 > > scalarETFESpacePtr_Type;
    //template <class Mesh>
    typedef boost::shared_ptr< FESpace< Mesh, MapEpetra >  > FESpacePtr_Type;

    typedef VectorEpetra           vector_Type;
    typedef boost::shared_ptr<vector_Type>         vectorPtr_Type;

    typedef MatrixEpetra<Real>           matrix_Type;
    typedef boost::shared_ptr<matrix_Type>         matrixPtr_Type;

    typedef EMData          data_Type;
    typedef typename boost::shared_ptr<data_Type>  dataPtr_Type;

    EMMaterialType (std::string materialName, UInt n);
    virtual ~EMMaterialType()   {}

    inline std::string& materialName()
    {
        return M_materialName;
    }

    inline vectorMaterialsPtr_Type materialFunctionList()
    {
        return M_materialFunctionList;
    }

    virtual void
    computeJacobian ( const vector_Type& disp,
                      ETFESpacePtr_Type dispETFESpace,
                      const vector_Type& fibers,
                      const vector_Type& sheets,
                      matrixPtr_Type     jacobianPtr)
    {
    	std::cout << "\nEMMaterial: wrong call to computeJacobian";
    }

    virtual void
    computeJacobian ( const vector_Type& disp,
                      ETFESpacePtr_Type  dispETFESpace,
                      const vector_Type& fibers,
                      const vector_Type& sheets,
                      const vector_Type& fiberActivation,
                      const vector_Type& sheetActivation,
                      const vector_Type& normalActivation,
                      scalarETFESpacePtr_Type  activationETFESpace,
                      matrixPtr_Type     jacobianPtr)
    {
    	std::cout << "\nEMMaterial: wrong call to computeJacobian";
    }

    virtual void
    computeResidual ( const vector_Type& disp,
                      ETFESpacePtr_Type  dispETFESpace,
                      const vector_Type& fibers,
                      const vector_Type& sheets,
                      vectorPtr_Type     residualVectorPtr)
    {
    	std::cout << "\nEMMaterial: wrong call to computeResidual";
    }

    virtual void
    computeResidual ( const vector_Type& disp,
                      ETFESpacePtr_Type  dispETFESpace,
                      const vector_Type& fibers,
                      const vector_Type& sheets,
                      const vector_Type& fiberActivation,
                      const vector_Type& sheetActivation,
                      const vector_Type& normalActivation,
                      scalarETFESpacePtr_Type  activationETFESpace,
                      vectorPtr_Type     residualVectorPtr)
    {
    	std::cout << "\nEMMaterial: wrong call to computeResidual";
    }

    inline void showMe()
    {
        std::cout << "Material Type: " << M_materialName << "\n";
        int n = M_materialFunctionList.size();
        for (int j (0); j < n; j++)
        {
            std::cout << "\nShowing material function: " << j << std::endl;
            M_materialFunctionList[j]->showMe();
        }
        std::cout << "\nDebug Output: end EMMaterial_TYpe::showMe()\n";

    }


    void setParameters(dataPtr_Type data)
    {
    	setParameters(*data);
    }

    virtual void setParameters(data_Type& data)
    {
        int n = M_materialFunctionList.size();
        for (int j (0); j < n; j++)
        {
            M_materialFunctionList[j]->setParameters(data);
        }
    }

protected:
    std::string M_materialName;
    vectorMaterialsPtr_Type M_materialFunctionList;

};


template<typename Mesh>
EMMaterialType<Mesh>::EMMaterialType (std::string materialName, UInt n ) :
    M_materialName (materialName),
    M_materialFunctionList (n)
{
	std::cout << "\nCreating: " << materialName << " with " << n << " functions.\n";
}



}//LifeV
#endif /* EMMATERIALDATA_HPP_ */
