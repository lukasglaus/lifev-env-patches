/*
 * EMData.hpp
 *
 *  Created on: 26/ott/2014
 *      Author: srossi
 */

#ifndef EMDATA_HPP_
#define EMDATA_HPP_


#include <lifev/core/filter/GetPot.hpp>

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

private:
	void setupSolver(GetPot& dataFile);
	Teuchos::ParameterList M_solverParametersList;
	Teuchos::ParameterList M_materialParametersList;
};

} /* namespace LifeV*/

#endif /* EMDATA_HPP_ */
