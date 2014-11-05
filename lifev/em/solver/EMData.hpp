/*
 * EMData.hpp
 *
 *  Created on: 26/ott/2014
 *      Author: srossi
 */

#ifndef EMDATA_HPP_
#define EMDATA_HPP_

#include <lifev/core/LifeV.hpp>


#include <lifev/core/filter/GetPot.hpp>
#include <Teuchos_ParameterList.hpp>


namespace LifeV {

class EMData {
public:
	virtual ~EMData() {}

	void setup(GetPot& dataFile);

	void setSolverParameterList(Teuchos::ParameterList& list)
	{
		M_solverParametersList = list;
	}

	void setMaterialParameterList(Teuchos::ParameterList& list)
	{
		M_materialParametersList = list;
	}

	void setupMaterialParamters(GetPot& dataFile, const std::string& section = "solid");

    const Teuchos::ParameterList& parametersList()
    {
    	return M_materialParametersList;
    }

    const Teuchos::ParameterList& solverParametersList()
    {
    	return M_solverParametersList;
    }

    Real parameter(std::string parameterName)
    {
    	return M_materialParametersList.get(parameterName, 0.0);
    }

    template< class Type >
    void setParameter(std::string parameterName, Type type)
    {
    	M_materialParametersList.set(parameterName, type);
    }

private:
	void setupSolver(GetPot& dataFile);
	Teuchos::ParameterList M_solverParametersList;
	Teuchos::ParameterList M_materialParametersList;
};

} /* namespace LifeV*/

#endif /* EMDATA_HPP_ */
