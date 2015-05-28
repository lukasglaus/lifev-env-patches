#include <iostream>
#include <Eigen/Dense>

#include "BodyCirculation.h"

using namespace std;
using namespace Eigen;

int main()
{
	BodyCirculation bc;

	bc.readInput("inputfile");

	bc.buildSystem(0.001);

	bc.addBC("lv", "sa", 0.5);
	bc.addBC("rv", 1);
	bc.addBC("sc", 5);

	for (unsigned int i = 0; i < 10; ++i)
	{
		bc.updateRhs();
		bc.updateDiodes();
		// bc.updateVariableBC();

		bc.solveSystem();
		bc.showValues();
	}
}
