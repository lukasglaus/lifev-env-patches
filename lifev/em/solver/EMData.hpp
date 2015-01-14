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

	void setSolidParameterList(Teuchos::ParameterList& list)
	{
		M_solidParametersList = list;
	}
	void setActivationParameterList(Teuchos::ParameterList& list)
	{
		M_activationParametersList = list;
	}
	void setElectroParameterList(Teuchos::ParameterList& list)
	{
		M_electroParametersList = list;
	}

	void setupSolidParameters(GetPot& dataFile, const std::string& section = "solid");
	void setupActivationParameters(GetPot& dataFile, const std::string& section = "activation");
	void setupElectrophysiologyParameters(GetPot& dataFile, const std::string& section = "electrophysiology");

    const Teuchos::ParameterList& solidParametersList()
    {
    	return M_solidParametersList;
    }
    const Teuchos::ParameterList& activationParametersList()
    {
    	return M_activationParametersList;
    }
    const Teuchos::ParameterList& electroParametersList()
    {
    	return M_electroParametersList;
    }

    const Teuchos::ParameterList& solverParametersList()
    {
    	return M_solverParametersList;
    }

    template< class Type >
    Type solidParameter(std::string parameterName)
    {
    	return M_solidParametersList.get(parameterName, Type() );
    }


    template< class Type >
    void setSolidParameter(std::string parameterName, Type type)
    {
    	M_solidParametersList.set(parameterName, type);
    }

    template< class Type >
    Type activationParameter(std::string parameterName)
    {
    	return M_activationParametersList.get(parameterName, Type() );
    }


    template< class Type >
    void setActivationParameter(std::string parameterName, Type type)
    {
    	M_activationParametersList.set(parameterName, type);
    }


    template< class Type >
    Type electroParameter(std::string parameterName)
    {
    	return M_electroParametersList.get(parameterName, Type() );
    }

    template< class Type >
    void setElectroParameter(std::string parameterName, Type type)
    {
    	M_electroParametersList.set(parameterName, type);
    }

private:
	void setupSolver(GetPot& dataFile);
	Teuchos::ParameterList M_solidParametersList;
	Teuchos::ParameterList M_activationParametersList;
	Teuchos::ParameterList M_electroParametersList;
	Teuchos::ParameterList M_solverParametersList;
};

} /* namespace LifeV*/

#endif /* EMDATA_HPP_ */
