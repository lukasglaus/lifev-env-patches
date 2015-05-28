/*
 * BodyCirculation.h
 *
 *  Created on: Dec 9, 2014
 *      Author: kummerth
 */

#ifndef BODYCIRCULATION_H_
#define BODYCIRCULATION_H_

#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <vector>
#include <string>
#include <list>
#include <algorithm>

using namespace std;
using namespace Eigen;

class BodyCirculation {
public:
	BodyCirculation();
	virtual ~BodyCirculation();

	void readInput(
			const string & filename); // read input file

	void buildSystem(
			const double & timestep); // build matrices
	
	void addBC( // Todo: time dependent B.C.
			const string & element,
			const double & value);

	void addBC( // Todo: time dependent B.C.
			const string & element1,
			const string & element2,
			const double & value);

	void setupInitCond( // Todo
			const string & element,
			const double & value);

	void setupInitCond( // Todo
			const string & element1,
			const string & element2,
			const double & value);

	void updateRhs(); // Todo

	void updateDiodes();

	void updateVariableBC(); // Todo

	void solveSystem(
			const double & timestep = 0);

	void showValues();

	void exportValues(
			const string * variables,
			const string & filename,
			const double & time); // Todo


private:
  
// Read input file
	int readInputCountLines(
		const string & filename);

	string readInputTypeValue(
		const string & line,
		const string & elemTypeName);

	const int getVectorInd(
		const vector<string> elements,
		const string word);

	const vector<int> getVectorInd(
		const vector<string> elements,
		const string word,
		const bool partOfString);

	const int getFreeRowIdx();

	typedef vector<string> vecS;
	typedef vector<vecS> matS;
	typedef vector<double> vecD;
	typedef vector<vecD> matD;

	matS d_nodes;
	matD d_values;
	vector<string> d_elements_vec;

// Build system matrices
	MatrixXd d_lhs;
	VectorXd d_rhs;
	double d_timestep;
	MatrixXd d_valveStates;

// Solve linear system
	VectorXd d_X;
	VectorXd d_X_scratch;
	Vector2d d_p_vent;

};


#endif /* BODYCIRCULATION_H_ */
