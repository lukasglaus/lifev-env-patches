/*
 * EMMaterialData.hpp
 *
 *  Created on: 29/apr/2014
 *      Author: srossi
 */

#ifndef EMPassiveMaterialType_HPP_
#define EMPassiveMaterialType_HPP_



#include <lifev/em/solver/mechanics/materials/EMMaterialType.hpp>

namespace LifeV
{

template<typename Mesh>
class EMPassiveMaterialType : public virtual EMMaterialType<Mesh>
{
public:

    typedef typename boost::shared_ptr<MaterialFunctions::EMMaterialFunctions<Mesh> >    materialFunctionsPtr_Type;
    typedef std::vector<materialFunctionsPtr_Type>  vectorMaterialsPtr_Type;

    typedef FactorySingleton<Factory<EMPassiveMaterialType<Mesh>, std::string> >  EMPassiveMaterialFactory;

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

    typedef StructuralConstitutiveLawData          data_Type;
    typedef typename boost::shared_ptr<data_Type>  dataPtr_Type;

    typedef EMMaterialType<Mesh> super;

    EMPassiveMaterialType (std::string materialName = "None", UInt n = 0);
    virtual ~EMPassiveMaterialType()   {}

    void
    computeJacobian ( const vector_Type& disp,
                      ETFESpacePtr_Type dispETFESpace,
                      const vector_Type& fibers,
                      const vector_Type& sheets,
                      matrixPtr_Type     jacobianPtr);

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
    computeResidual ( const vector_Type& disp,
                      ETFESpacePtr_Type  dispETFESpace,
                      const vector_Type& fibers,
                      const vector_Type& sheets,
                      vectorPtr_Type           residualVectorPtr);

    void
    computeResidual ( const vector_Type& disp,
                      ETFESpacePtr_Type  dispETFESpace,
                      const vector_Type& fibers,
                      const vector_Type& sheets,
                      const vector_Type& activation,
                      scalarETFESpacePtr_Type  activationETFESpace,
                      vectorPtr_Type     residualVectorPtr);



};


template<typename Mesh>
EMPassiveMaterialType<Mesh>::EMPassiveMaterialType (std::string materialName, UInt n ) :
super(materialName, n)
{
}



template<typename Mesh>
void
EMPassiveMaterialType<Mesh>::computeJacobian ( const vector_Type& disp,
                                        ETFESpacePtr_Type  dispETFESpace,
                                        const vector_Type& fibers,
                                        const vector_Type& sheets,
                                        matrixPtr_Type           jacobianPtr)
{
    int n = this->M_materialFunctionList.size();
    for (int j (0); j < n; j++)
    {
    	this->M_materialFunctionList[j]->computeJacobian (disp, dispETFESpace, fibers, sheets, jacobianPtr);
    }
}


template<typename Mesh>
void
EMPassiveMaterialType<Mesh>::computeJacobian ( const vector_Type& disp,
                                        ETFESpacePtr_Type   dispETFESpace,
                                        const vector_Type& fibers,
                                        const vector_Type& sheets,
                                        const vector_Type& activation,
                                        scalarETFESpacePtr_Type   activationETFESpace,
                                        matrixPtr_Type           jacobianPtr)
{
    int n = this->M_materialFunctionList.size();
    for (int j (0); j < n; j++)
    	this->M_materialFunctionList[j]->computeJacobian (disp,
														dispETFESpace,
														fibers,
														sheets,
														activation,
														activationETFESpace,
														jacobianPtr);
}

template <typename Mesh>
void
EMPassiveMaterialType<Mesh>::computeResidual ( const vector_Type& disp,
                                        ETFESpacePtr_Type  dispETFESpace,
                                        const vector_Type& fibers,
                                        const vector_Type& sheets,
                                        vectorPtr_Type           residualVectorPtr)
{
    if (residualVectorPtr)
    {
        //std::cout << "EM Material Type: ResidualVectorPtr available\n";
    }
    if (dispETFESpace)
    {
        //std::cout << "EM Material Type: dispETFESpace available\n";
    }
    int n = this->M_materialFunctionList.size();
    for (int j (0); j < n; j++)
    {
    	//std::cout << "Passive residual function " << j << " = " << residualVectorPtr->norm2() << "\n";
    	this->M_materialFunctionList[j]->computeResidual (disp, dispETFESpace, fibers, sheets, residualVectorPtr);

    }
}


template <typename Mesh>
void
EMPassiveMaterialType<Mesh>::computeResidual ( const vector_Type& disp,
                                        ETFESpacePtr_Type   dispETFESpace,
                                        const vector_Type& fibers,
                                        const vector_Type& sheets,
                                        const vector_Type& activation,
                                        scalarETFESpacePtr_Type   activationETFESpace,
                                        vectorPtr_Type           residualVectorPtr)
{
    int n = this->M_materialFunctionList.size();
    for (int j (0); j < n; j++)
    	this->M_materialFunctionList[j]->computeResidual (disp,
                                                    dispETFESpace,
                                                    fibers,
                                                    sheets,
                                                    activation,
                                                    activationETFESpace,
                                                    residualVectorPtr);
}

}//LifeV
#endif /* EMMATERIALDATA_HPP_ */
