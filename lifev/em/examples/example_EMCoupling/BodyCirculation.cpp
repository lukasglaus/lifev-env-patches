/*
 * BodyCirculation.cpp
 *
 *  Created on: Dec 9, 2014
 *      Author: kummerth
 */

#include "BodyCirculation.h"


BodyCirculation::BodyCirculation() {
	// TODO Auto-generated constructor stub
}

BodyCirculation::~BodyCirculation() {
	// TODO Auto-generated destructor stub
}

void BodyCirculation::readInput(
	const string & filename)
{
	int numLines = readInputCountLines(filename);

	d_nodes.resize(numLines, vecS(2));
	d_values.resize(numLines, vecD(6));

	string line;
	ifstream file (filename.c_str());
	if (file.is_open())
	{
	  int i = 0;
	  while ( getline (file, line) )
	  {
		  d_nodes[i][0] = readInputTypeValue(line, "N1");
		  d_nodes[i][1] = readInputTypeValue(line, "N2");

		  d_values[i][0] = atof( readInputTypeValue(line, "R").c_str() );
		  d_values[i][1] = atof( readInputTypeValue(line, "L").c_str() );
		  d_values[i][2] = atof( readInputTypeValue(line, "C").c_str() );

		  string valuesDiode = readInputTypeValue(line, "D");
		  size_t dashOne = valuesDiode.find("-");
		  size_t dashTwo = valuesDiode.find("-", dashOne+1);
		  d_values[i][3] = atof( valuesDiode.substr(0, dashOne).c_str() );
		  d_values[i][4] = atof( valuesDiode.substr(dashOne+1, dashTwo-dashOne-1).c_str() );
		  d_values[i][5] = atof( valuesDiode.substr(dashTwo+1).c_str() );

		  ++i;
	  }
	  file.close();
	}

	else cout << "Unable to open file";

	// Create unique list of elements (nodes & connections)
	list<string> elements_list;
	for (unsigned int i = 0; i < d_nodes.size(); ++i)
	{
		string N1 = d_nodes[i][0];
		string N2 = d_nodes[i][1];
		elements_list.push_back (N1);
		elements_list.push_back (N2);
		elements_list.push_back (N1 + "-" + N2);
	}
	elements_list.sort();
	elements_list.unique();
	vector<string> elements_vec (elements_list.begin(), elements_list.end());

	d_elements_vec.resize(elements_vec.size());
	d_elements_vec = elements_vec;

	return;
}// readInput

int BodyCirculation::readInputCountLines(
	const string & filename)
{
	  int numLines = 0;
	  string line;
	  ifstream file (filename.c_str());
	  if (file.is_open())
	  {
		  while ( getline (file, line) )
		  {
			  ++numLines;
		  }
		  file.close();
	  }

	  else cout << "Unable to open file";

	  return numLines;
}// readInputCountLines

string BodyCirculation::readInputTypeValue(
	const string & line,
	const string & elemTypeName)
{

	  size_t elemType = line.find(elemTypeName);
	  size_t braceOne = line.find("(", elemType);
	  size_t braceTwo = line.find(")", elemType);

	  return ( (elemType != string::npos) ? line.substr(braceOne+1, braceTwo-braceOne-1) : "0");

}// readInputTypeValue

void BodyCirculation::buildSystem(
		const double & timestep)
{
	// Build empty matrices
	d_lhs.resize(d_elements_vec.size(), d_elements_vec.size());
	d_lhs.setZero();
	d_rhs.resize(d_elements_vec.size());
	d_rhs.setZero();
	d_X.resize(d_elements_vec.size());
	d_X.setZero();
	d_X_scratch.resize(d_elements_vec.size());
	d_X_scratch.setZero();

	d_timestep = timestep;

	// Add R and L equations
	for (unsigned int i = 0; i < d_nodes.size(); ++i)
	{
		string N1 = d_nodes[i][0];
		string N2 = d_nodes[i][1];
		int ind = getVectorInd(d_elements_vec, N1 + "-" + N2);
		d_lhs(ind, getVectorInd(d_elements_vec, N1 + "-" + N2)) = d_values[i][0] + d_values[i][1]/d_timestep;
		d_lhs(ind, getVectorInd(d_elements_vec, N1)) = -1;
		d_lhs(ind, getVectorInd(d_elements_vec, N2)) = 1;
		d_rhs(ind) = d_X_scratch(ind)*d_values[i][1]/d_timestep;
	}


	// Initialize Diodes
	unsigned int valveIt = 0;
	for (unsigned int i = 0; i < d_nodes.size(); ++i)
	{
		if (d_values[i][3] != 0)
		{
			d_valveStates.resize(d_valveStates.rows() + 1, 6);

			d_valveStates(valveIt, 0) = 0; // Todo: Initial values
			d_valveStates(valveIt, 1) = 0; // Todo: Initial values

			string N1 = d_nodes[i][0];
			string N2 = d_nodes[i][1];
			int ind = getVectorInd(d_elements_vec, N1 + "-" + N2);
			d_lhs(ind, getVectorInd(d_elements_vec, N1)) = -1*(1 - d_valveStates(valveIt, 0));
			d_lhs(ind, getVectorInd(d_elements_vec, N2)) = 1*(1 - d_valveStates(valveIt, 0));

			d_valveStates(valveIt, 2) = i;
			d_valveStates(valveIt, 3) = ind;
			d_valveStates(valveIt, 4) = getVectorInd(d_elements_vec, N1);
			d_valveStates(valveIt, 5) = getVectorInd(d_elements_vec, N2);
			++ valveIt;
		}
	}


	// Add C equations
	for (unsigned int i = 0; i < d_nodes.size(); ++i)
	{
		if ( (d_values[i][0] == 0) && (d_values[i][1] == 0) )
		{
			string N1 = d_nodes[i][0];
			string N2 = d_nodes[i][1];
			int ind = getVectorInd(d_elements_vec, N1 + "-" + N2);
			d_lhs(ind, getVectorInd(d_elements_vec, N1 + "-" + N2)) = -1;
			d_lhs(ind, getVectorInd(d_elements_vec, N1)) = d_values[i][2] / d_timestep;
			d_lhs(ind, getVectorInd(d_elements_vec, N2)) = -d_values[i][2] / d_timestep;
			d_rhs(ind) = d_X_scratch( getVectorInd(d_elements_vec, N1) )*d_values[i][2] / d_timestep - d_X_scratch( getVectorInd(d_elements_vec, N1) )*d_values[i][2] / d_timestep;
		}
	}

	// Add Q equations
	for (unsigned int i = 0; i < d_elements_vec.size(); ++i)
	{
		string element = d_elements_vec[i];
		int ind = getVectorInd(d_elements_vec, element);
		vector<int> ind1 = getVectorInd(d_elements_vec, element + "-", 1);
		vector<int> ind2 = getVectorInd(d_elements_vec, "-" + element, 1);

		if ( (ind1.size() + ind2.size()) > 1 )
		{
			for (unsigned j = 0; j < ind1.size(); ++j)
			{
				d_lhs(ind, ind1[j]) = -1;
			}

			for (unsigned j = 0; j < ind2.size(); ++j)
			{
				d_lhs(ind, ind2[j]) = 1;
			}
		}
	}

	// Plot variable names
	for (unsigned int i = 0; i < d_elements_vec.size() ; ++i)
	{
		cout << d_elements_vec[i] << endl;
	}

	return;
}// buildSystem


void BodyCirculation::updateDiodes()
{
	for (unsigned int i = 0; i < d_valveStates.rows(); ++i)
	{
		const double & valN1 = d_X(d_valveStates(i, 4));
		const double & valN2 = d_X(d_valveStates(i, 5));
		const double & alpha = d_values[d_valveStates(i, 2)][3];
		const double & gammaUpp = d_values[d_valveStates(i, 2)][4];
		const double & gammaLow = d_values[d_valveStates(i, 2)][5];
		const double & prevValveState = d_valveStates(i, 1);
		const double & currValveState = d_valveStates(i, 0);

		double Np = 2*currValveState - prevValveState - (valN1 - valN2 != 0 ? pow(d_timestep, 2.0)*alpha*(valN1 - valN2)/abs(valN1 - valN2) : 0);
		double NUpp = (currValveState + gammaUpp*d_timestep) / (1 + gammaUpp*d_timestep);
		double NLow = currValveState / (1 + gammaLow*d_timestep);

		d_valveStates(i, 1) = d_valveStates(i, 0);
		d_valveStates(i, 0) = max( Np, min( NUpp, NLow ) );

		d_lhs(d_valveStates(i, 3), d_valveStates(i, 4)) = -1*(1 - d_valveStates(i, 0));
		d_lhs(d_valveStates(i, 3), d_valveStates(i, 5)) = 1*(1 - d_valveStates(i, 0));
	}
}

void BodyCirculation::updateRhs()
{
	// L
	for (unsigned int i = 0; i < d_nodes.size(); ++i)
	{
		string N1 = d_nodes[i][0];
		string N2 = d_nodes[i][1];
		int ind = getVectorInd(d_elements_vec, N1 + "-" + N2);
		d_rhs(ind) = d_X_scratch(ind)*d_values[i][1]/d_timestep;
	}

	// C
	for (unsigned int i = 0; i < d_nodes.size(); ++i)
	{
		if ( (d_values[i][0] == 0) && (d_values[i][1] == 0) )
		{
			string N1 = d_nodes[i][0];
			string N2 = d_nodes[i][1];
			int ind = getVectorInd(d_elements_vec, N1 + "-" + N2);
			d_rhs(ind) = d_X_scratch( getVectorInd(d_elements_vec, N1) )*d_values[i][2]/d_timestep - d_X_scratch( getVectorInd(d_elements_vec, N2) )*d_values[i][2]/d_timestep;
		}
	}
}

const int BodyCirculation::getFreeRowIdx()
{
	unsigned int idx = 0;
	for (unsigned int i = 0; i < d_lhs.rows(); ++i)
	{
		for (unsigned int j = 0; j < d_lhs.cols(); ++j)
		{
			if (d_lhs(i,j) != 0) break;
			if (j == (d_lhs.cols() - 1)) idx = i;
		}
	}
	return idx;
}

void BodyCirculation::addBC(
		const string & element,
		const double & value)
{
	int freeRow = getFreeRowIdx();
	int ind = getVectorInd(d_elements_vec, element);
	d_lhs(freeRow, ind) = 1;
	d_rhs(freeRow) = value;
}

void BodyCirculation::addBC(
		const string & element1,
		const string & element2,
		const double & value)
{
	int freeRow = getFreeRowIdx();
	int ind = getVectorInd(d_elements_vec, element1 + "-" + element2);
	d_lhs(freeRow, ind) = 1;
	d_rhs(freeRow) = value;
}

const int BodyCirculation::getVectorInd(
	const vector<string> elements,
	const string word)
{
	int ind;
	for (unsigned int i = 0; i < elements.size(); ++i)
	{
		if (elements[i] == word)
		{
			ind = i;
			break;
		}
	}
	return ind;
}// getListIndex

const vector<int> BodyCirculation::getVectorInd(
	const vector<string> elements,
	const string word,
	const bool partOfString)
{
	vector<int> ind;
	for (unsigned int i = 0; i < elements.size(); ++i)
	{
		string element = elements[i];
		if (element.find(word) != string::npos)
		{
			ind.push_back(i);
		}
	}
	return ind;
}// getListIndex

void BodyCirculation::solveSystem(
	const double & timestep)
{
	if (timestep != d_timestep && timestep != 0)
	{
		// rebuild system marices
		// re-apply b.c.
	}
	d_X = d_lhs.fullPivLu().solve(d_rhs);
	d_X_scratch = d_X;
	return;
}// solveSystem

void BodyCirculation::showValues()
{
	cout << "----------------------------------------------------------------------------" << endl;
	cout << d_lhs << endl << endl;
	cout << d_rhs << endl << endl;
	cout << d_X << endl << endl;;
	return;
}// getVentricularPressure
