/*
 * EssentialPatchBCPatchReal.h
 *
 *  Created on: May 21, 2019
 *      Author: pamstad
 */

#ifndef ESSENTIALPATCHBCPATCHREAL_HPP_
#define ESSENTIALPATCHBCPATCHREAL_HPP_

#include <stdio.h>
#include <cmath>
#include <math.h> //this is for arctan
#include <lifev/core/array/MatrixSmall.hpp>
#include <lifev/em/examples/example_EMHeart/EssentialPatchBC.hpp>

#define PI 3.14159265359

namespace LifeV
{

class EssentialPatchBCPatchReal : public EssentialPatchBC
{
public:


	typedef MatrixSmall<3,3>					matrixSmall_Type;
	//then we define a small matrix with matrixSmall_Type nameOfMatrix; we can fill entries of the matrix with nameOfMatrix (0,0) = 1.0 ; starts by zero; so max number is 2

	virtual void setup(const GetPot& dataFile, const std::string& name)
	{
		std::cout << "We are in setup function of BC specific " << std::endl;

		super::setup(dataFile, name); //here flag, name , direction and so on gets red in

		m_Phi = dataFile(("solid/boundary_conditions/" + m_Name + "/phi").c_str(), 0.0); //phi is angle for first z Axis Rotation
		m_Theta = dataFile(("solid/boundary_conditions/" + m_Name + "/theta").c_str(), 0.0); //theta is angle of rotation about the y' axis (note that it is y prime and not y!)
		//m_Psi = dataFile(("solid/boundary_conditions/" + m_Name + "/psi").c_str(), 0.0); //this is rotation about z'' (z double prime) axis

		m_a = dataFile(("solid/boundary_conditions/" + m_Name + "/a").c_str(), 0.0); //this is semimajor axis in x-direction
		m_b = dataFile(("solid/boundary_conditions/" + m_Name + "/b").c_str(), 0.0); //this is semimajor axis in y-direction
		m_c = dataFile(("solid/boundary_conditions/" + m_Name + "/c").c_str(), 0.0); //This is semimajor axis in z-direction

		m_Height = dataFile(("solid/boundary_conditions/" + m_Name + "/height").c_str(), 0.0); //this is height (measured in z-direction) of the patch
		m_Width = dataFile(("solid/boundary_conditions/" + m_Name + "/width").c_str(), 0.0); //this is width (measured in y-direction) of the patch

		m_shift = dataFile(("solid/boundary_conditions/" + m_Name + "/shift").c_str(), 0.0); //this is shift along the axis of patch; only required when no vertexEllipse is given as input


		m_maxDisplacement = dataFile ( ("solid/boundary_conditions/" + m_Name + "/displacement").c_str(), 1.0 );

		for (UInt j(0); j < 3 ; j++)
		{
			m_vertexEllipse[j] = dataFile(("solid/boundary_conditions/" + m_Name + "/vertexEllipse").c_str(), 0.0, j);
		}

		if(m_vertexEllipse[0] != 1000) //this is when you want to give vertexPoint as input and not angles
		{
			initialiseVertexEllipse();
		}

		initialiseRotationMatrices();
		calculateyMaxzMax();
		initialseMinMaxValues();
		initialiseEllipsoidMatrix();

		std::cout << "This is value of phi: " << m_Phi << std::endl;
		std::cout << "This is value of xMin: " << m_xMin << std::endl;
		std::cout << "This is value of xMax: " << m_xMax << std::endl;
		std::cout << "This is value of yMin: " << m_yMin << std::endl;
		std::cout << "This is value of yMax: " << m_yMax << std::endl;
		std::cout << "This is value of zMin: " << m_zMin << std::endl;
		std::cout << "This is value of zMax: " << m_zMax << std::endl;
	}


	virtual const bool nodeOnPatch(const Vector3D& coord, const Real& time)
	{

		bool nodeOnPatch = false;

		//we need all three coordinates of transformed points
		//so far when taking the if statement out then he finds points but way to much

		//if(coord[0] <= (m_xMax - std::cos(m_Phi*PI/180.0)*m_maxDisplacement) && coord[0] >= (m_xMin - std::cos(m_Phi*PI/180.0)*m_maxDisplacement) && coord[1] <= (m_yMax - std::sin(m_Phi*PI/180.0)*m_maxDisplacement) && coord[1] >= (m_yMin - std::sin(m_Phi*PI/180.0)*m_maxDisplacement)   && coord[2] < m_zMax && coord[2] > m_zMax )
		if(coord[0] <= (m_xMax - std::cos(m_Phi*PI/180.0)*m_maxDisplacement) && coord[0] >= (m_xMin - std::cos(m_Phi*PI/180.0)*m_maxDisplacement) && coord[1] <= (m_yMax) && coord[1] >= (m_yMin )   && coord[2] < m_zMax && coord[2] > m_zMax )
		{
			Vector3D intermediateResult;
			intermediateResult[0] = 0.0;
			intermediateResult[1] = 0.0;
			intermediateResult[2] = 0.0;

			Vector3D xVector;
			xVector[0] = coord[0] - m_xShift + std::cos(m_Phi*PI/180.0)*m_maxDisplacement; //std::cos(m_Phi*PI/180.0)*m_maxDisplacement; //I think we need + here for the displacement because we want to shift in negative x direction
			xVector[1] = coord[1] - m_yShift ; //+ std::sin(m_Phi*PI/180.0)*m_maxDisplacement;	//Also need to rethink if we want both in x and y just the same because I think we need to weight them// shifting the points we need oppsite sign than the function
			xVector[2] = coord[2];

			//here we calculate the first part of ellipse equation; [A](x-fo)

			intermediateResult = matrixVectorMultiplicator(m_Ellipsoid, xVector);
			Real ellipseEquation = xVector.dot(intermediateResult)-1.0;

			//if a point has these coordinates it is possible that it is part of the patch; but we still need to check by plugging into ellipse equation
			// this will look something like (x - shift)*[A]*(x-shift) - 1 >= 0 ; if it is larger than zero then it is outside so we want to move it
			if(ellipseEquation >= 0)
			{
					nodeOnPatch = true;
			}
			else
			{
					nodeOnPatch = false;
			}
		 }

		return nodeOnPatch;

}


	virtual const bool nodeOnPatchCurrent(const Vector3D& coord, const Real& time)
	{
		//here we plug the coordinates into the point

		bool nodeOnPatch = false;

		Vector3D intermediateResult;
		intermediateResult[0] = 0.0;
		intermediateResult[1] = 0.0;
		intermediateResult[2] = 0.0;

		Vector3D xVector;
		xVector[0] = coord[0] - m_xShift + std::cos(m_Phi*PI/180)*activationFunction(time);
		xVector[1] = coord[1] - m_yShift + std::sin(m_Phi*PI/180)*activationFunction(time);
		xVector[2] = coord[2]- m_vertexEllipse[2];

		//here we calculate the first part of ellipse equation; [A](x-fo)

		intermediateResult = matrixVectorMultiplicator(m_Ellipsoid, xVector);

		Real ellipseEquation = xVector.dot(intermediateResult)-1;

		//we need all three coordinates of transformed points
		if(coord[0] <= (m_xMax - std::cos(m_Phi*PI/180.0)*activationFunction(time)) && coord[0] >= (m_xMin - std::cos(m_Phi*PI/180.0)*activationFunction(time)) && coord[1] <= (m_yMax - std::sin(m_Phi*PI/180.0)*activationFunction(time)) && coord[1] >= (m_yMin - std::sin(m_Phi*PI/180.0)*activationFunction(time))   && coord[2] < m_zMax && coord[2] > m_zMax )
		{
			//if a point has these coordinates it is possible that it is part of the patch; but we still need to check by plugging into ellipse equation
			// this will look something like (x - shift)*[A]*(x-shift) - 1 >= 0 ; if it is larger than zero then it is outside so we want to move it
			if(ellipseEquation >= 0)
			{
				nodeOnPatch = true;
			}
			else
			{
				nodeOnPatch = false;
			}
		}

		return nodeOnPatch;

	}

	virtual void modifyPatchArea(EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver,const int& newFlag, const Real& time)
	{
		if ( solver.comm()->MyPID() == 0 ) std::cout << "WE ARE IN MODIFY PATCH AREA " << std::endl;

		auto p2FeSpace = solver.electroSolverPtr()->feSpacePtr();
		auto p2dFeSpace = solver.structuralOperatorPtr()->dispFESpacePtr();
		FESpace<RegionMesh<LinearTetra>, MapEpetra > p1FESpace (p2FeSpace->mesh(), "P1", 1, p2FeSpace->mesh()->comm());

		//create an epetra vector to set it equal to one where it is part of patch
		VectorEpetra p1ScalarFieldFaces (p1FESpace.map());
    	p1ScalarFieldFaces *= 0.0;
    	Int p1ScalarFieldFacesDof = p1ScalarFieldFaces.epetraVector().MyLength();

    	int globalIdArray[p1ScalarFieldFacesDof];
    	p1ScalarFieldFaces.blockMap().MyGlobalElements(globalIdArray);

        m_patchFlag = newFlag;
        const auto& mesh = solver.localMeshPtr();
		const auto& meshfull = solver.fullMeshPtr();
				            // Create patches by changing the markerID (flag) locally
		unsigned int numNodesOnPatch(0); //here we just initalise an unsigned integer variable
        for (int j(0); j < mesh->numBoundaryFacets(); j++) //returns number of boundary facets
        {
		      auto& face = mesh->boundaryFacet(j);
		      auto faceFlag = face.markerID();

		      int numPointsOnFace(0);

		      for (int k(0); k < 3; ++k) //k < 3 was before; this is just a test
		      {

		    	  ID pointGlobalId = face.point(k).id();
		          auto coord = face.point(k).coordinates();
		          auto pointInPatch = nodeOnPatchCurrent(coord, time);

		          if(pointInPatch == true)
		          {
		          		++numPointsOnFace;
		           		for(int n = 0; n < p1ScalarFieldFacesDof; n++)
		          		{
		           			if(pointGlobalId == globalIdArray[n])
		           			{
		                  		p1ScalarFieldFaces[pointGlobalId] = 1.0;

		           			}
		                }
		           }

		        }

		}


		            m_patchFacesLocationPtr.reset (new vector_Type (p2FeSpace->map() ));
		            *m_patchFacesLocationPtr = p2FeSpace->feToFEInterpolate(p1FESpace, p1ScalarFieldFaces);
	            //*m_patchFacesLocationPtr = p1ScalarFieldFaces;
 }


	virtual vectorPtr_Type directionalVectorField (EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver ,const boost::shared_ptr<FESpace<RegionMesh<LinearTetra>, MapEpetra >> dFeSpace, Vector3D& direction, const Real& disp, const Real& time)
	{

			auto p2PositionVector = p2PositionVectorInitial(dFeSpace, solver);

			vectorPtr_Type p2PatchDisplacement (new VectorEpetra( dFeSpace->map(), Repeated ));
	        auto nCompLocalDof = p2PatchDisplacement->epetraVector().MyLength() / 3;

	        for (int j (0); j < nCompLocalDof; ++j)
	        {
	        	                    // Get coordinates

	                 UInt iGID = p2PatchDisplacement->blockMap().GID (j);
	                 UInt jGID = p2PatchDisplacement->blockMap().GID (j + nCompLocalDof);
	                 UInt kGID = p2PatchDisplacement->blockMap().GID (j + 2 * nCompLocalDof);

	                 (*p2PatchDisplacement)[iGID] = 0.0;
	                 (*p2PatchDisplacement)[jGID] = 0.0;
	                 (*p2PatchDisplacement)[kGID] = 0.0;
	        }

	        return p2PatchDisplacement;


	}




	void initialiseVertexEllipse()
	{
		//first we calculate phi angle

		if(m_vertexEllipse[0] > 0 && m_vertexEllipse[1] > 0)
		{
			m_Phi = atan(m_vertexEllipse[1]/m_vertexEllipse[0])*180/PI;
		}

		if(m_vertexEllipse[0] < 0 && m_vertexEllipse[1] > 0)
		{
			m_Phi = atan(m_vertexEllipse[1]/m_vertexEllipse[0])*180/PI + 180;
		}

		if(m_vertexEllipse[0] > 0 && m_vertexEllipse[1] < 0)
		{
			m_Phi = atan(m_vertexEllipse[1]/m_vertexEllipse[0])*180/PI;
		}

		if(m_vertexEllipse[0] < 0 && m_vertexEllipse[1] < 0)
		{
			m_Phi = -180 + atan(m_vertexEllipse[1]/m_vertexEllipse[0])*180/PI;
		}
		if(m_vertexEllipse[0] == 0)
		{
			m_Phi = 90;
		}

		//next we calculate midpoint
		m_midPoint[0] = m_vertexEllipse[0] - m_a*std::cos(m_Phi*PI/180);
		m_midPoint[1] = m_vertexEllipse[1] - m_a*std::sin(m_Phi*PI/180);
		m_midPoint[2] = 0.0;

		std::cout << "These are coordiantes of midPoint: " << m_midPoint[0] << "      " << m_midPoint[1] << "         " << m_midPoint[2] << std::endl;

		m_shift = m_midPoint.norm();
		m_xShift = std::cos(m_Phi*PI/180)*m_shift;
		m_yShift = std::sin(m_Phi*PI/180)*m_shift;

		std::cout << "This is value of xShift: " << m_xShift << std::endl;
		std::cout << "This is value of yShift: " << m_yShift << std::endl;


	}


	void initialiseRotationMatrices()
	{
		//this is only to compute the three different Rotation Matrices and then also the totalRotation matrices
		//we start with the first Rotation

		//////FIRST ROTATION MATRIX
		m_RotationOne (0,0) = std::cos(m_Phi*PI/180.0);
		m_RotationOne (0,1) = -std::sin(m_Phi*PI/180.0);
		m_RotationOne (0,2) = 0.0;

		m_RotationOne (1,0) = std::sin(m_Phi*PI/180.0);
		m_RotationOne (1,1) = std::cos(m_Phi*PI/180.0);
		m_RotationOne (1,2) = 0.0;

		m_RotationOne (2,0) = 0.0;
		m_RotationOne (2,1) = 0.0;
		m_RotationOne (2,2) = 1.0;
		//////////////


		//std::cout << "here comes rotationMatrix: " << std::endl;
		//printMatrix(m_RotationOne);



		///////SECOND ROTATION MATRIX
		m_RotationTwo (0,0) = std::cos(m_Theta*PI/180.0);
		m_RotationTwo (0,1) = 0.0;
		m_RotationTwo (0,2) = -std::sin(m_Theta*PI/180.0);

		m_RotationTwo (1,0) = 0.0;
		m_RotationTwo (1,1) = 1.0;
		m_RotationTwo (1,2) = 0.0;

		m_RotationTwo (2,0) = std::sin(m_Theta*PI/180.0);
		m_RotationTwo (2,1) = 0.0;
		m_RotationTwo (2,2) = std::cos(m_Theta*PI/180.0);
		/////////////

		///////THIRD ROTATION MATRIX
		/*
		m_RotationThree (0,0) = std::cos(m_Psi*PI/180.0);
		m_RotationThree (0,1) = std::sin(m_Psi*PI/180.0);
		m_RotationThree (0,2) = 0.0;

		m_RotationThree (1,0) = -std::sin(m_Psi*PI/180.0);
		m_RotationThree (1,1) = std::cos(m_Psi*PI/180.0);
		m_RotationThree (1,2) = 0.0;

		m_RotationThree (2,0) = 0.0;
		m_RotationThree (2,1) = 0.0;
		m_RotationThree (2,2) = 1.0;
		*/
		/////////////

	}

	void calculateyMaxzMax()
	{
		//procedure is as described in paper; first we get the heigt of the patch by the user
		if(m_Height <= 2*m_c)
		{
			m_zMax = m_Height/2;
		}
		else
		{
			m_zMax == m_c;
			std::cout << "Your given height of the patch exceeds 2*c: check it in dataFile " << std::endl;
		}





		Real f = 1 - std::pow(m_zMax,2.0)/std::pow(m_c,2.0);
		m_yMax = m_Width/2;




		if(m_yMax > m_b*sqrt(f))
		{
			std::cout << "Your given width is too large; check the description and adopt axis (a,b,c) of ellipsoid; y_Max is set to the maximum possible value" << std::endl;
			m_yMax = m_b*sqrt(f) - 0.01; //0.01 is just to be save that we find solution later
		}

		m_xValueOne = m_a*sqrt(f - std::pow(m_yMax,2.0)/std::pow(m_b,2.0));
		m_xValueTwo = m_a*sqrt(1 - std::pow(m_yMax, 2.0)/std::pow(m_b, 2.0));
		//Now we intialise the edge Point of the Patch which we later need to transform if a given rotation is needed; again check the script of this function

		std::cout << "This is xValueOne: " << m_xValueOne << std::endl;
		std::cout << "This is xValueTwo: " << m_xValueTwo << std::endl;


		m_PointOne[0] = m_xValueOne;
		m_PointOne[1] = m_yMax;
		m_PointOne[2] = m_zMax;

		m_PointTwo[0] = m_xValueOne;
		m_PointTwo[1] = -m_yMax;
		m_PointTwo[2] = m_zMax;

		m_xMaxPoint[0] = m_xValueTwo;
		m_xMaxPoint[1] = -m_yMax;
		m_xMaxPoint[2] = 0;

		m_yMaxPoint[0] = m_xValueTwo;
		m_yMaxPoint[1] = m_yMax;
		m_yMaxPoint[2] = 0;

		m_pointSemiAxis[0] = m_a; //in MATLAB Script this is point yMaxPoint2
		m_pointSemiAxis[1] = 0.0;
		m_pointSemiAxis[2] = 0.0;

		if(m_vertexEllipse[0] == 1000)
		{
			m_xShift = std::cos(m_Phi*PI/180)*m_shift;
			m_yShift = std::sin(m_Phi*PI/180)*m_shift;
		}

	}

	void initialseMinMaxValues()
	{
		//In general, finding the min and max values of the patch area is tricky; to get a better understanding what in this function is done look at my report
		//////////////// This is range for phi from 0 to 90 degrees

		//we need to do this initalisation a little bit better
		Vector3D pointOneRotated = matrixVectorMultiplicator(m_RotationOne, m_PointOne);
		Vector3D pointTwoRotated = matrixVectorMultiplicator(m_RotationOne, m_PointTwo);
		Vector3D pointyMaxRotated = matrixVectorMultiplicator(m_RotationOne, m_yMaxPoint);
		Vector3D pointxMaxRotated = matrixVectorMultiplicator(m_RotationOne, m_xMaxPoint);
		Vector3D pointSemiAxisRotated = matrixVectorMultiplicator(m_RotationOne, m_pointSemiAxis);

		std::cout << "This is pointOneRotated: " << pointOneRotated[0] << "    " << pointOneRotated[1] << "     " << pointOneRotated[2] << std::endl;
		std::cout << "This is pointTwoRotated: " << pointTwoRotated[0] << "    " << pointTwoRotated[1] << "     " << pointTwoRotated[2] << std::endl;



		if(m_Phi < 5 && m_Phi >=0)
		{
			m_xMin = pointOneRotated[0] + m_xShift;
			m_xMax = 1.5*m_a + m_xShift;

			m_yMin = pointTwoRotated[1] + m_yShift;
			m_yMax = pointyMaxRotated[1] + m_yShift;

			m_zMin = -m_zMax + m_vertexEllipse[2];
			m_zMax = m_zMax + m_vertexEllipse[2];
		}

		if(m_Phi >= 0 && m_Phi <=84)
		{
			m_xMin = pointOneRotated[0] + m_xShift;
			m_xMax = pointxMaxRotated[0] + m_xShift;

			m_yMin = pointTwoRotated[1] + m_yShift;
			m_yMax = pointyMaxRotated[1] + m_yShift;

			m_zMin = -m_zMax + m_vertexEllipse[2];
			m_zMax = m_zMax + m_vertexEllipse[2];
		}

		if(m_Phi > 84 && m_Phi <= 90)
		{
			m_xMin = pointOneRotated[0]+ m_xShift;
			m_xMax = pointxMaxRotated[0]+ m_xShift;

			m_yMin = pointTwoRotated[1] + m_yShift;
			m_yMax = 1.5*m_a + m_yShift;

			m_zMin = -m_zMax+ m_vertexEllipse[2];
			m_zMax = m_zMax + m_vertexEllipse[2];
		}

		if(m_Phi > 90 && m_Phi < 95)
		{
			m_xMin = pointyMaxRotated[0]+ m_xShift;
			m_xMax = pointTwoRotated[0]+ m_xShift;

			m_yMin = pointOneRotated[1] + m_yShift;
			m_yMax = 1.5*m_a + m_yShift;

			m_zMin = -m_zMax + m_vertexEllipse[2];
			m_zMax = m_zMax + m_vertexEllipse[2];
		}

		if(m_Phi >= 95 && m_Phi <= 175)
		{
			m_xMin = pointyMaxRotated[0]+ m_xShift;
			m_xMax = pointTwoRotated[0]+ m_xShift;

			m_yMin = pointOneRotated[1] + m_yShift;
			m_yMax = pointxMaxRotated[1] + m_yShift; //0.2 is just to be sure

			m_zMin = -m_zMax + m_vertexEllipse[2];
			m_zMax = m_zMax + m_vertexEllipse[2];
		}

		if(m_Phi > 175 && m_Phi <= 180)
		{

			m_xMin = -m_a*1.5+ m_xShift;
			m_xMax = pointTwoRotated[0]+ m_xShift;

			m_yMin = pointOneRotated[1] + m_yShift;
			m_yMax = pointxMaxRotated[1] + m_yShift; //0.2 is just to be sure

			m_zMin = -m_zMax+ m_vertexEllipse[2];
			m_zMax = m_zMax + m_vertexEllipse[2];

		}

		if(m_Phi >= -4 && m_Phi <0)
		{
			m_xMin = pointTwoRotated[0]+ m_xShift;
			m_xMax = 1.5*m_a+ m_xShift;

			m_yMin = pointxMaxRotated[1] + m_yShift;
			m_yMax = pointOneRotated[1] + m_yShift;

			m_zMin = -m_zMax+ m_vertexEllipse[2];
			m_zMax = m_zMax + m_vertexEllipse[2];

		}

		if(m_Phi < -4 && m_Phi >= -85)
		{
			m_xMin = pointTwoRotated[0]+ m_xShift;
			m_xMax = pointyMaxRotated[0]+ m_xShift;

			m_yMin = pointxMaxRotated[1] + m_yShift;
			m_yMax = pointOneRotated[1] + m_yShift;

			m_zMin = -m_zMax+ m_vertexEllipse[2];
			m_zMax = m_zMax + m_vertexEllipse[2];

		}

		if(m_Phi >= -90 && m_Phi < -85)
		{
			m_xMin = pointTwoRotated[0]+ m_xShift;
			m_xMax = pointyMaxRotated[0]+ m_xShift;

			m_yMin = -1.5*m_a + m_yShift;
			m_yMax = pointOneRotated[1] + m_yShift;

			m_zMin = -m_zMax+ m_vertexEllipse[2];
			m_zMax = m_zMax + m_vertexEllipse[2];

		}

		if(m_Phi <= -90 && m_Phi >= -94)
		{
			m_xMin = pointxMaxRotated[0]+ m_xShift;
			m_xMax = pointOneRotated[0]+ m_xShift;

			m_yMin = -1.5*m_a + m_yShift;
			m_yMax = pointTwoRotated[1] + m_yShift;

			m_zMin = -m_zMax+ m_vertexEllipse[2];
			m_zMax = m_zMax + m_vertexEllipse[2];

		}

		if(m_Phi <= -95 && m_Phi >= -175)
		{
			m_xMin = pointxMaxRotated[0]+ m_xShift;
			m_xMax = pointOneRotated[0]+ m_xShift;

			m_yMin = pointyMaxRotated[1] + m_yShift;
			m_yMax = pointTwoRotated[1] + m_yShift;

			m_zMin = -m_zMax+ m_vertexEllipse[2];
			m_zMax = m_zMax + m_vertexEllipse[2];

		}

		if(m_Phi <= -176 && m_Phi >= -180)
		{
			m_xMin = -1.5*m_a+ m_xShift;
			m_xMax = pointOneRotated[0]+ m_xShift;

			m_yMin = pointyMaxRotated[1] + m_yShift;
			m_yMax = pointTwoRotated[1] + m_yShift;

			m_zMin = -m_zMax+ m_vertexEllipse[2];
			m_zMax = m_zMax + m_vertexEllipse[2];

		}

	}


	void printMatrix(matrixSmall_Type matrix)
	{

		for(int i=0; i < 3; i++)
		{

				std::cout << matrix(i,0) << "   " << matrix(i,1) << "  " << matrix(i,2) << std::endl;

		}

	}


	void initialiseEllipsoidMatrix()
	{

		//need to recheck this matrix; because I have changed something

		matrixSmall_Type intermediateResult;

		for(int i(0); i < 3; i++) // just set values of intermediateResult to zero
		{
			for(int j(0); j < 3; j++)
			{
				intermediateResult (i,j) = 0.0;
			}
		}

		//here we initalise the matrix for the Ellispoid
		m_A (0, 0) = 1/std::pow(m_a,2.0);
		m_A (1, 1) = 1/std::pow(m_b,2.0);
		m_A (2, 2) = 1/std::pow(m_c,2.0);

		for(int i(0); i<3;i++) // here we just set other values to zero
		{
			for(int j(0); j < 3; j++)
			{
				if (i != j)
				{
					m_A (i, j) = 0.0;
				}
			}
		}


		intermediateResult = matrixMatrixMultiplicator(m_A, m_RotationOne.transpose());

		m_Ellipsoid = matrixMatrixMultiplicator(m_RotationOne, intermediateResult);

		std::cout << "This is Ellipsoid Matrix: " << std::endl;
		printMatrix(m_Ellipsoid);


		//note this is now only with two rotations; but this needs to be rechecked because I have changed rotation1
		/*
		m_Ellipsoid (0,0) = std::pow(std::sin(m_Theta*PI/180),2.0)/std::pow(m_c,2.0) + std::pow(std::cos(m_Theta*PI/180),2.0)*(std::pow(std::cos(m_Phi*PI/180),2.0)/std::pow(m_a,2.0)+std::pow(std::sin(m_Phi*PI/180),2.0)/std::pow(m_b,2));
		m_Ellipsoid (0,1) = -std::cos(m_Theta*PI/180)*(std::cos(m_Phi*PI/180)*std::sin(m_Phi*PI/180)/std::pow(m_a,2) - std::cos(m_Phi*PI/180)*std::sin(m_Phi*PI/180)/std::pow(m_b,2));
		m_Ellipsoid (0,2) = std::cos(m_Theta*PI/180)*std::sin(m_Theta*PI/180)*(std::pow(std::cos(m_Phi*PI/180),2)/std::pow(m_a,2) + std::pow(std::sin(m_Phi*PI/180),2)/std::pow(m_b,2) - std::cos(m_Theta*PI/180)*std::sin(m_Theta*PI/180)/std::pow(m_c,2));

		m_Ellipsoid (1,0) = -std::cos(m_Theta*PI/180)*(std::cos(m_Phi*PI/180)*std::sin(m_Phi*PI/180)/std::pow(m_a,2) - std::cos(m_Phi*PI/180)*std::sin(m_Phi*PI/180)/std::pow(m_b,2));
		m_Ellipsoid (1,1) = std::pow(std::cos(m_Phi),2.0)/std::pow(m_b,2.0) + std::pow(std::sin(m_Phi*PI/180),2.0)/std::pow(m_a,2.0);
		m_Ellipsoid (1,2) = - std::sin(m_Theta*PI/180)*(std::cos(m_Phi*PI/180)*std::sin(m_Phi*PI/180)/std::pow(m_a,2) - std::cos(m_Phi*PI/180)*std::sin(m_Phi*PI/180)/std::pow(m_b,2));

		m_Ellipsoid (2,0) = std::cos(m_Theta*PI/180)*std::sin(m_Theta*PI/180)*(std::pow(std::cos(m_Phi*PI/180),2.0)/std::pow(m_a,2) + std::pow(std::sin(m_Phi), 2.0)/std::pow(m_b,2.0)) - std::cos(m_Theta*PI/180)*std::sin(m_Theta*PI/180)/std::pow(m_c,2.0);
		m_Ellipsoid (2,1) = - std::sin(m_Theta*PI/180)*(std::cos(m_Phi*PI/180)*std::sin(m_Phi*PI/180)/std::pow(m_a,2)) - std::cos(m_Phi*PI/180)*std::sin(m_Phi*PI/180)/std::pow(m_b,2.0);
		m_Ellipsoid (2,2) = std::pow(std::cos(m_Theta),2)/std::pow(m_c,2) + std::pow(std::sin(m_Theta*PI/180),2)*(std::pow(std::cos(m_Phi*PI/180),2.0)/std::pow(m_a,2) + std::pow(std::sin(m_Phi*PI/180), 2.0)/std::pow(m_b,2));
		 */
		//WHAT A PAIN

	}

	Vector3D matrixVectorMultiplicator(matrixSmall_Type matrix, Vector3D vector)
	{

		Vector3D result;
		result[0] = 0.0;
		result[1] = 0.0;
		result[2] = 0.0;

		for(int i=0; i < 3; i++)
		{
			for(int j=0; j < 3; j++)
			{
				//std::cout << "This is matrix entry: " << matrix[i][j] << std::endl;
				//std::cout << "This is vector entry: " << vector[j] << std::endl;

				result[i] += matrix[i][j]*vector[j];
			}
		}

		return result;

	}

	matrixSmall_Type matrixMatrixMultiplicator(matrixSmall_Type matrixA, matrixSmall_Type matrixB) //this function gives result of matrixA*matrixB
	{
		matrixSmall_Type result;
		for(int i(0); i < 3; i++)
		{
			for(int j(0); j < 3; j++)
			{
				result (i,j) = 0.0;
				for(int k(0); k < 3; k++)
				{
					result (i, j) += matrixA(i, k)* matrixB(k, j);
				}
			}
		}

		return result;
	}

protected:

	//we are gonna define some member variables

	Real m_Phi;
	Real m_Theta;
	Real m_Psi;

	Real m_a;
	Real m_b;
	Real m_c;

	Real m_Height;
	Real m_Width;

	//These min and max values are here for defining the size of the patch; that means width and height
	Real m_xValueOne;
	Real m_xValueTwo;
	Real m_xMin;
	Real m_xMax;
	Real m_yMin;
	Real m_yMax;
	Real m_zMin;
	Real m_zMax;

	Vector3D m_PointOne;
	Vector3D m_PointTwo;
	Vector3D m_yMaxPoint;
	Vector3D m_xMaxPoint;
	Vector3D m_pointSemiAxis;

	Real m_maxDisplacement;
	Vector3D m_vertexEllipse;
	Vector3D m_midPoint;
	Real m_shift;
	Real m_xShift;
	Real m_yShift;

	matrixSmall_Type m_RotationOne;
	matrixSmall_Type m_RotationTwo;
	matrixSmall_Type m_RotationThree;

	matrixSmall_Type m_A; //Matrix m_A is just a diagonal Matrix with the Halbachsenwerten a,b ,c squared auf der Diagonalen; also m_A = diag(1/a^2, 1/b^2, 1/c^2)
	matrixSmall_Type m_Ellipsoid;

};



REGISTER(EssentialPatchBC, EssentialPatchBCPatchReal);

}
#endif /* ESSENTIALPATCHBCPATCHREAL_HPP_ */
