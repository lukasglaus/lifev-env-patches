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
    computeJacobian ( const vector_Type& disp,
                      ETFESpacePtr_Type  dispETFESpace,
                      const vector_Type& fibers,
                      const vector_Type& sheets,
                      const vectorPtr_Type& fiberActivation,
                      const vectorPtr_Type& sheetActivation,
                      const vectorPtr_Type& normalActivation,
                      scalarETFESpacePtr_Type  activationETFESpace,
                      matrixPtr_Type     jacobianPtr);
    //





    void
    computeResidual ( const vector_Type& disp,
                      ETFESpacePtr_Type  dispETFESpace,
                      const vector_Type& fibers,
                      const vector_Type& sheets,
                      const vectorPtr_Type& fiberActivation,
                      const vectorPtr_Type& sheetActivation,
                      const vectorPtr_Type& normalActivation,
                      scalarETFESpacePtr_Type  activationETFESpace,
                      vectorPtr_Type     residualVectorPtr);

//    inline void showMe()
//    {
//        std::cout << "Material Type: " << M_materialName << "\n";
//    }
//    virtual void setParameters(super::data_Type& data)
//    {
//        int n = M_materialFunctionList.size();
//        for (int j (0); j < n; j++)
//        {
//            M_materialFunctionList[j]->setParameters(data);
//        }xx
//    }


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
							                    const vectorPtr_Type& fiberActivation,
							                    const vectorPtr_Type& sheetActivation,
							                    const vectorPtr_Type& normalActivation,
							                    scalarETFESpacePtr_Type  activationETFESpace,
												matrixPtr_Type     jacobianPtr)
{
    int n = this->M_materialFunctionList.size();
    for (int j (0); j < n; j++)
    {
//    	if(sheetActivation && normalActivation)
//    	{
			this->M_materialFunctionList[j]->computeJacobian ( disp,
															   dispETFESpace,
															   fibers,
															   sheets,
															   fiberActivation,
															   sheetActivation,
															   normalActivation,
															   activationETFESpace,
															   jacobianPtr);
//    	}
//    	else
//    	{
//    		if(sheetActivation)
//    		{
//    			this->M_materialFunctionList[j]->computeJacobian ( disp,
//    															   dispETFESpace,
//    															   fibers,
//    															   sheets,
//    															   *fiberActivation,
//    															   *sheetActivation,
//    															   activationETFESpace,
//    															   jacobianPtr);
//
//    		}
//    		else if(normalActivation)
//    		{
//    			this->M_materialFunctionList[j]->computeJacobian ( disp,
//    															   dispETFESpace,
//    															   fibers,
//    															   sheets,
//    															   *fiberActivation,
//    															   *normalActivation,
//    															   activationETFESpace,
//    															   jacobianPtr);
//
//    		}
//    		else
//    		{
//    			this->M_materialFunctionList[j]->computeJacobian ( disp,
//    															   dispETFESpace,
//    															   fibers,
//    															   sheets,
//    															   *fiberActivation,
//    															   activationETFESpace,
//    															   jacobianPtr);
//    		}
//    	}
    }
}

template <typename Mesh>
void
EMActiveMaterialType<Mesh>::computeResidual  (  const vector_Type& disp,
												ETFESpacePtr_Type  dispETFESpace,
												const vector_Type& fibers,
												const vector_Type& sheets,
							                    const vectorPtr_Type& fiberActivation,
							                    const vectorPtr_Type& sheetActivation,
							                    const vectorPtr_Type& normalActivation,
							                    scalarETFESpacePtr_Type  activationETFESpace,
												vectorPtr_Type     residualVectorPtr)
{
    int n = this->M_materialFunctionList.size();
    for (int j (0); j < n; j++)
    {
//    	if(sheetActivation && normalActivation)
//    	{
			this->M_materialFunctionList[j]->computeResidual ( disp,
															   dispETFESpace,
															   fibers,
															   sheets,
															   fiberActivation,
															   sheetActivation,
															   normalActivation,
															   activationETFESpace,
															   residualVectorPtr );
//    	}
//    	else
//    	{
//    		if(sheetActivation)
//    		{
//    			this->M_materialFunctionList[j]->computeResidual ( disp,
//    															   dispETFESpace,
//    															   fibers,
//    															   sheets,
//    															   *fiberActivation,
//    															   *sheetActivation,
//    															   activationETFESpace,
//    															   residualVectorPtr );
//    		}
//    		else if (normalActivation)
//    		{
//    			this->M_materialFunctionList[j]->computeResidual ( disp,
//    															   dispETFESpace,
//    															   fibers,
//    															   sheets,
//    															   *fiberActivation,
//    															   *normalActivation,
//    															   activationETFESpace,
//    															   residualVectorPtr );
//    		}
//    		else
//    		{
//    			this->M_materialFunctionList[j]->computeResidual ( disp,
//    															   dispETFESpace,
//    															   fibers,
//    															   sheets,
//    															   *fiberActivation,
//    															   activationETFESpace,
//    															   residualVectorPtr );
//    		}
//    	}
    }
}



}//LifeV
#endif /* EMMATERIALDATA_HPP_ */
