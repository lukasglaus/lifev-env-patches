/*
 * EMMaterialData.hpp
 *
 *  Created on: 29/apr/2014
 *      Author: srossi
 */

#ifndef EMMATERIALTYPE_HPP_
#define EMMATERIALTYPE_HPP_

#include <lifev/em/solver/mechanics/materials/functions/FunctionsList.hpp>

namespace LifeV
{

template<typename Mesh>
class EMMaterialType
{
public:

	typedef typename boost::shared_ptr<MaterialFunctions::EMMaterialFunctions<Mesh> >    materialFunctionsPtr_Type;
	typedef std::vector<materialFunctionsPtr_Type>  vectorMaterialsPtr_Type;

    typedef FactorySingleton<Factory<EMMaterialType<Mesh>, std::string> >  EMMaterialFactory;


	typedef Mesh mesh_Type;
	//template <class Mesh>
	using ETFESpacePtr_Type = boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 3 > >;
	using scalarETFESpacePtr_Type = boost::shared_ptr<ETFESpace<Mesh, MapEpetra, 3, 1 > >;
	//template <class Mesh>
	using FESpacePtr_Type = boost::shared_ptr< FESpace< Mesh, MapEpetra >  >;

	typedef VectorEpetra           vector_Type;
	typedef boost::shared_ptr<vector_Type>         vectorPtr_Type;

	typedef MatrixEpetra<Real>           matrix_Type;
	typedef boost::shared_ptr<matrix_Type>         matrixPtr_Type;


	EMMaterialType(std::string materialName = "None", UInt n = 0);
	virtual ~EMMaterialType()	{}


	inline std::string& materialName()
	{
		return M_materialName;
	}

	inline vectorMaterialsPtr_Type materialFunctionList()
	{
		return M_materialFunctionList;
	}

	void
	computeJacobian( const vector_Type& disp,
					 ETFESpacePtr_Type dispETFESpace,
					 const vector_Type& fibers,
					 const vector_Type& sheets,
					 matrixPtr_Type     jacobianPtr);

	void
	computeJacobian( const vector_Type& disp,
					 ETFESpacePtr_Type  dispETFESpace,
					 const vector_Type& fibers,
					 const vector_Type& sheets,
					 const vector_Type& activation,
			         scalarETFESpacePtr_Type  activationETFESpace,
			         matrixPtr_Type     jacobianPtr);
	//



	void
	computeResidual( const vector_Type& disp,
		     ETFESpacePtr_Type  dispETFESpace,
		     const vector_Type& fibers,
		     const vector_Type& sheets,
					 vectorPtr_Type           residualVectorPtr);

	void
	computeResidual( const vector_Type& disp,
					 ETFESpacePtr_Type  dispETFESpace,
					 const vector_Type& fibers,
					 const vector_Type& sheets,
					 const vector_Type& activation,
			         scalarETFESpacePtr_Type  activationETFESpace,
			         vectorPtr_Type     residualVectorPtr);

	inline void showMe()
	{
		std::cout << "Material Type: " << M_materialName << "\n";
	}

protected:
	std::string M_materialName;
	vectorMaterialsPtr_Type M_materialFunctionList;

};


template<typename Mesh>
EMMaterialType<Mesh>::EMMaterialType(std::string materialName, UInt n ) :
		M_materialName(materialName),
		M_materialFunctionList(n)
{

}



template<typename Mesh>
void
EMMaterialType<Mesh>::computeJacobian( const vector_Type& disp,
									ETFESpacePtr_Type  dispETFESpace,
									 const vector_Type& fibers,
									 const vector_Type& sheets,
									matrixPtr_Type           jacobianPtr)
{
	int n = M_materialFunctionList.size();
	for(int j(0); j < n; j++)
		M_materialFunctionList[j]->computeJacobian(disp, dispETFESpace, fibers, sheets, jacobianPtr);
}


template<typename Mesh>
void
EMMaterialType<Mesh>::computeJacobian( const vector_Type& disp,
									   ETFESpacePtr_Type   dispETFESpace,
									   const vector_Type& fibers,
									   const vector_Type& sheets,
									   const vector_Type& activation,
									   scalarETFESpacePtr_Type   activationETFESpace,
									   matrixPtr_Type           jacobianPtr)
{
	int n = M_materialFunctionList.size();
	for(int j(0); j < n; j++)
		M_materialFunctionList[j]->computeJacobian(disp,
				                                   dispETFESpace,
				                                   fibers,
				                                   sheets,
				                                   activation,
				                                   activationETFESpace,
				                                   jacobianPtr);
}

template <typename Mesh>
void
EMMaterialType<Mesh>::computeResidual( const vector_Type& disp,
										 ETFESpacePtr_Type  dispETFESpace,
										 const vector_Type& fibers,
										 const vector_Type& sheets,
										 vectorPtr_Type           residualVectorPtr)
{
	int n = M_materialFunctionList.size();
	for(int j(0); j < n; j++)
		M_materialFunctionList[j]->computeResidual(disp, dispETFESpace, fibers, sheets, residualVectorPtr);
}


template <typename Mesh>
void
EMMaterialType<Mesh>::computeResidual( const vector_Type& disp,
									   ETFESpacePtr_Type   dispETFESpace,
									   const vector_Type& fibers,
									   const vector_Type& sheets,
									   const vector_Type& activation,
									   scalarETFESpacePtr_Type   activationETFESpace,
									   vectorPtr_Type           residualVectorPtr)
{
	int n = M_materialFunctionList.size();
	for(int j(0); j < n; j++)
		M_materialFunctionList[j]->computeResidual(disp,
				                                   dispETFESpace,
				                                   fibers,
				                                   sheets,
				                                   activation,
				                                   activationETFESpace,
				                                   residualVectorPtr);
}

}//LifeV
#endif /* EMMATERIALDATA_HPP_ */
