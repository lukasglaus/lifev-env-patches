/*
 * Activation.hpp
 *
 *  Created on: 28/dic/2014
 *      Author: srossi
 */

#ifndef ACTIVATION_HPP_
#define ACTIVATION_HPP_

#include <lifev/core/util/Factory.hpp>
#include <lifev/core/util/FactorySingleton.hpp>

#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/array/MapEpetra.hpp>

#include <lifev/em/solver/EMData.hpp>

#include <lifev/core/mesh/RegionMesh.hpp>


namespace LifeV
{

class Activation
{
public:

    typedef VectorEpetra                                                vector_Type;

    typedef boost::shared_ptr<VectorEpetra>                             vectorPtr_Type;

    typedef FactorySingleton<Factory<Activation, std::string> >  EMActivationFactory;


	virtual ~Activation() {}

	virtual void setup( EMData& data,  const MapEpetra& map) = 0;
	virtual void setParameters( EMData& data ) = 0;

    virtual void setupActivationPtrs(	vectorPtr_Type& fiberActivationPtr,
    									vectorPtr_Type& sheetActivationPtr,
    									vectorPtr_Type& normalActivationPtr ) = 0;

    virtual void updateActivation(	vectorPtr_Type& fiberActivationPtr,
    								vectorPtr_Type& sheetActivationPtr,
    								vectorPtr_Type& normalActivationPtr ) = 0;

    virtual void solveModel(Real& timeStep) = 0;

    virtual void solveModelPathology ( Real& timeStep, boost::shared_ptr<RegionMesh<LinearTetra> > fullMeshPtr, const boost::shared_ptr<FESpace<RegionMesh<LinearTetra>, MapEpetra >> dFeSpace) {}


    VectorEpetra& fiberActivation()
    {
        return *M_fiberActivationPtr;
    }

    VectorEpetra& fiberActivation() const
    {
        return *M_fiberActivationPtr;
    }

    vectorPtr_Type fiberActivationPtr()
    {
        return M_fiberActivationPtr;
    }

    void setFiberActivation (VectorEpetra& activation)
    {
        *M_fiberActivationPtr = activation;
    }

    void setFiberActivationPtr (vectorPtr_Type activationPtr)
    {
        M_fiberActivationPtr = activationPtr;
    }



    VectorEpetra& sheetActivation()
    {
        return *M_sheetActivationPtr;
    }

    vectorPtr_Type sheetActivationPtr()
    {
        return M_sheetActivationPtr;
    }

    void setSheetActivation (VectorEpetra& activation)
    {
        *M_sheetActivationPtr = activation;
    }

    void setSheetActivationPtr (vectorPtr_Type activationPtr)
    {
        M_sheetActivationPtr = activationPtr;
    }


    VectorEpetra& normalActivation()
    {
        return *M_normalActivationPtr;
    }

    vectorPtr_Type normalActivationPtr()
    {
        return M_normalActivationPtr;
    }

    void setNormalActivation (VectorEpetra& activation)
    {
        *M_normalActivationPtr = activation;
    }

    void setNormalActivationPtr (vectorPtr_Type activationPtr)
    {
        M_normalActivationPtr = activationPtr;
    }


    VectorEpetra& fiber()
    {
        return *M_fiberPtr;
    }

    vectorPtr_Type fiberPtr()
    {
        return M_fiberPtr;
    }

    void setFiberPtr(vectorPtr_Type fiberPtr)
    {
    	M_fiberPtr = fiberPtr;
    }

    VectorEpetra& sheet()
    {
        return *M_sheetPtr;
    }

    vectorPtr_Type sheetPtr()
    {
        return M_sheetPtr;
    }

    void setSheetPtr(vectorPtr_Type sheetPtr)
    {
    	M_sheetPtr = sheetPtr;
    }


    VectorEpetra& I4f()
    {
        return *M_I4fPtr;
    }

    vectorPtr_Type I4fPtr()
    {
        return M_I4fPtr;
    }

    void setI4fPtr(vectorPtr_Type i4fPtr)
    {
    	M_I4fPtr = i4fPtr;
    }

    void setI4f(const vector_Type& i4f)
    {
    	*M_I4fPtr = i4f;
    }

    template < class ElectroSolver >
    void setVariablesPtr(const ElectroSolver& solver)
    {

    	int size = solver.globalSolution().size();
    	M_electroSolution.reserve(size);
    	for(int i = 0; i<size; i++)
    	{
    		M_electroSolution.push_back(  solver.globalSolution().at(i) );
    	}
    }

    virtual void showMe()
    {
    	if(M_fiberActivationPtr->comm().MyPID() == 0)
    	{
    		std::cout << "Activation - Overload the showMe() method to get information about the activation solver you are currently using!" << std::endl;
    	}
    }


protected:
    vectorPtr_Type M_fiberActivationPtr;
    vectorPtr_Type M_sheetActivationPtr;
    vectorPtr_Type M_normalActivationPtr;
    vectorPtr_Type M_I4fPtr;

    std::vector < vectorPtr_Type >  M_electroSolution;

    vectorPtr_Type M_fiberPtr;
    vectorPtr_Type M_sheetPtr;

};

} /* namespace LifeV */

#endif /* ACTIVATION_HPP_ */
