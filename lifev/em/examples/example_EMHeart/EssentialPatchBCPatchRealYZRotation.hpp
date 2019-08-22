/*
 * EssentialPatchBCPatchRealYZRotation.h
 *
 *  Created on: Jun 3, 2019
 *      Author: pamstad
 */

#ifndef ESSENTIALPATCHBCPATCHREALYZROTATION_HPP_
#define ESSENTIALPATCHBCPATCHREALYZROTATION_HPP_

#include <stdio.h>
#include <typeinfo>
#include <cmath>
#include <math.h> //this is for arctan
#include <lifev/core/array/MatrixSmall.hpp>
#include <lifev/em/examples/example_EMHeart/EssentialPatchBC.hpp>

#define PI 3.14159265359

namespace LifeV
{

class EssentialPatchBCPatchRealYZRotation : public EssentialPatchBC
{


	typedef MatrixSmall<3,3>					matrixSmall_Type;
	//then we define a small matrix with matrixSmall_Type nameOfMatrix; we can fill entries of the matrix with nameOfMatrix (0,0) = 1.0 ; starts by zero; so max number is 2

	virtual void setup(const GetPot& dataFile, const std::string& name)
	{

		super::setup(dataFile, name); //here flag, name , direction and so on gets red in

		m_Phi = dataFile(("solid/boundary_conditions/" + m_Name + "/phi").c_str(), 0.0); //phi is angle for first z Axis Rotation
		m_Psi = dataFile(("solid/boundary_conditions/" + m_Name + "/psi").c_str(), 0.0); //theta is angle of rotation about the y' axis (note that it is y prime and not y!)
		//m_Psi = dataFile(("solid/boundary_conditions/" + m_Name + "/psi").c_str(), 0.0); //this is rotation about z'' (z double prime) axis

		m_a = dataFile(("solid/boundary_conditions/" + m_Name + "/a").c_str(), 0.0); //this is semimajor axis in x-direction
		m_b = dataFile(("solid/boundary_conditions/" + m_Name + "/b").c_str(), 0.0); //this is semimajor axis in y-direction
		m_c = dataFile(("solid/boundary_conditions/" + m_Name + "/c").c_str(), 0.0); //This is semimajor axis in z-direction

		m_Height = dataFile(("solid/boundary_conditions/" + m_Name + "/height").c_str(), 0.0); //this is height (measured in z-direction) of the patch
		m_Width = dataFile(("solid/boundary_conditions/" + m_Name + "/width").c_str(), 0.0); //this is width (measured in y-direction) of the patch

		m_maxDisplacement = dataFile ( ("solid/boundary_conditions/" + m_Name + "/displacement").c_str(), 1.0 );

		for (UInt j(0); j < 3 ; j++)
		{
			m_vertexEllipse[j] = dataFile(("solid/boundary_conditions/" + m_Name + "/vertexEllipse").c_str(), 0.0, j);
		}

		for (UInt j(0); j < 3 ; j++)
		{
			m_normalVector[j] = dataFile(("solid/boundary_conditions/" + m_Name + "/normalVector").c_str(), 0.0, j);
		}


		if(m_vertexEllipse[0] != 1000) //this is when you want to give vertexPoint as input and not angles
		{
			initialiseVertexEllipse();
		}

		initialiseRotationMatrices();
		initialiseHeartAxis();
		calculateyMaxzMax();
		initialiseEllipsoidBoundaries();
		//initialseMinMaxValues();
		initialiseEllipsoidMatrix();


		//std::cout << "This is value of psi: " << m_Psi << std::endl;

		
		std::cout << "This is value of phi: " << m_Phi << std::endl;
		//std::cout << "This is value of yMin: " << m_yMin << std::endl;
		//std::cout << "This is value of yMax: " << m_yMax << std::endl;
		//std::cout << "This is value of zMin: " << m_zMin << std::endl;
		//std::cout << "This is value of zMax: " << m_zMax << std::endl;
		
	}


	virtual const bool nodeOnPatch(const Vector3D& coord, const Real& time)
		{

			bool nodeOnPatch = false;
			bool coordInPatchRange = false;
			Vector3D intermediateResult;
			Vector3D xVector;


			coordInPatchRange = coordinatesInPatchRange(coord, time);

			//we need all three coordinates of transformed points
			//so far when taking the if statement out then he finds points but way to much

			if(coordInPatchRange == true)
			{

				Vector3D intermediateResult;
				intermediateResult[0] = 0.0;
				intermediateResult[1] = 0.0;
				intermediateResult[2] = 0.0;

				Vector3D xVector;
				/*	
				xVector[0] = coord[0] - m_xPointShift - m_patchDirection[0]*m_maxDisplacement; 
				xVector[1] = coord[1] - m_yPointShift;
				xVector[2] = coord[2] - m_zPointShift - m_patchDirection[2]*m_maxDisplacement;
				*/
				
				xVector[0] = coord[0] - m_xPointShift - m_patchDirection[0]*activationFunction(time); 
                                xVector[1] = coord[1] - m_yPointShift;
                                xVector[2] = coord[2] - m_zPointShift - m_patchDirection[2]*activationFunction(time);
				
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
			bool coordInPatchRange = false;
			Vector3D intermediateResult;
			Vector3D xVector;



			coordInPatchRange = coordinatesInPatchRangeCurrent(coord, time);


			//if(coord[0] <= (m_xMax - std::cos(m_Phi*PI/180.0)*activationFunction(time)) && coord[0] >= (m_xMin - std::cos(m_Phi*PI/180.0)*activationFunction(time)) && coord[1] <= (m_yMax) && coord[1] >= (m_yMin)   && coord[2] < (m_zMax- std::sin(m_Phi*PI/180.0)*activationFunction(time)) && coord[2] > (m_zMin- std::sin(m_Phi*PI/180.0)*activationFunction(time)))
			//{
			if(coordInPatchRange == true)
			{

				intermediateResult[0] = 0.0;
				intermediateResult[1] = 0.0;
				intermediateResult[2] = 0.0;

				Vector3D xVector;
				xVector[0] = coord[0] - m_xPointShift - m_patchDirection[0]*activationFunction(time);
				xVector[1] = coord[1] - m_yPointShift;
				xVector[2] = coord[2] - m_zPointShift - m_patchDirection[2]*activationFunction(time);
				/*
				xVector[0] = coord[0] - m_xPointShift + std::cos(m_Phi*PI/180.0)*activationFunction(time);
				xVector[1] = coord[1] - m_vertexEllipse[1];
				xVector[2] = coord[2] - m_zPointShift + std::sin(m_Phi*PI/180.0)*activationFunction(time);
				*/

				intermediateResult = matrixVectorMultiplicator(m_Ellipsoid, xVector);
				Real ellipseEquation = xVector.dot(intermediateResult)-1;
				//if a point has these coordinates it is possible that it is part of the patch; but we still need to check by plugging into ellipse equation
				// this will look something like (x - shift)*[A]*(x-shift) - 1 >= 0 ; if it is larger than zero then it is outside so we want to move it
				if(ellipseEquation >= 0)
				{
					if(coord[1] >= 1.1*m_yMax || coord[1] <= 0.9*m_yMin)
                                        {
                                                nodeOnPatch = nodeOnRadius(coord, time);
                                        }
					else
					{
						nodeOnPatch = true;
					}
				}
				else
				{
					nodeOnPatch = false;
				}
			}
			//}

			return nodeOnPatch;

		}

	const bool nodeOnRadius(const Vector3D& coord, const Real& time)
	{
		 bool nodeOnRadius = false;

                 Real fzeroX = m_xPointShift + m_patchDirection[0]*activationFunction(time);
                 Real fzeroY = m_yPointShift;
                 Real fzeroZ = m_zPointShift + m_patchDirection[2]*activationFunction(time);

                 Real a;
                 Real b;
                 Real c;


                 Real lambdaOne=0;
                 Real lambdaTwo = 0;
		 Real yDifference =0;
		 Real edgeSmoother = 0;

                         a = m_patchDirection[0]*(m_Ellipsoid(0, 0)*m_patchDirection[0] + m_Ellipsoid(2,0)*m_patchDirection[2]) + m_patchDirection[2]*(m_Ellipsoid(0,2)*m_patchDirection[0] + m_Ellipsoid(2, 2)*m_patchDirection[2]);

                         b = m_patchDirection[0]*(m_Ellipsoid(0, 0)*(coord[0]-fzeroX) + m_Ellipsoid(1,0)*(coord[1] - fzeroY) + m_Ellipsoid(2,0)*(coord[2] - fzeroZ)) + m_patchDirection[2]*(m_Ellipsoid(0,2)*(coord[0] - fzeroX) + m_Ellipsoid(1, 2)*(coord[1] - fzeroY) + m_Ellipsoid(2, 2)*(coord[2] - fzeroZ)) + (coord[0] - fzeroX)*(m_Ellipsoid(0,0)*m_patchDirection[0] + m_Ellipsoid(2,0)*m_patchDirection[2]) + (coord[1] - fzeroY)*(m_Ellipsoid(0,1)*m_patchDirection[0] + m_Ellipsoid(2,1)*m_patchDirection[2]) + (coord[2] - fzeroZ)*(m_Ellipsoid(0,2)*m_patchDirection[0] + m_Ellipsoid(2,2)*m_patchDirection[2]);

                        c = (coord[0] - fzeroX)*(m_Ellipsoid(0, 0)*(coord[0] - fzeroX) + m_Ellipsoid(1,0)*(coord[1] - fzeroY) + m_Ellipsoid(2, 0)*(coord[2] - fzeroZ)) + (coord[1] - fzeroY)*(m_Ellipsoid(0, 1)*(coord[0] - fzeroX) + m_Ellipsoid(1, 1)*(coord[1] - fzeroY) + m_Ellipsoid(2,1)*(coord[2] - fzeroZ)) + (coord[2] - fzeroZ)*(m_Ellipsoid(0,2)*(coord[0] - fzeroX) + m_Ellipsoid(1,2)*(coord[1] - fzeroY) + m_Ellipsoid(2,2)*(coord[2] - fzeroZ)) -1.0;

                        lambdaOne = (-b + sqrt(std::pow(b,2) - 4*a*c))/(2*a);
                        lambdaTwo = (-b - sqrt(std::pow(b,2) - 4*a*c))/(2*a);

                        bool checkRange = false;
                        checkRange = coord[1] <= 0.9*m_yMin;
                        if(checkRange == true)
                        {
                               yDifference = coord[1] - 0.9*m_yMin;
                        }
                        else
                        {
                               yDifference = coord[1] - 1.1*m_yMax;
                        }

                        edgeSmoother = 1*std::pow(yDifference, 2.0);
                 	Real displacement = std::abs(lambdaOne) < std::abs(lambdaTwo) ? std::abs(lambdaOne) : std::abs(lambdaTwo);
                        displacement = displacement - edgeSmoother;

                        if(displacement > 0.0)
                        {
                                nodeOnRadius = true;
                        }
                        else
                        {
                                nodeOnRadius = false;
                        }


                        return nodeOnRadius;
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
			
				//std::cout << "This is faceFlag in modifyPatchArea: " << faceFlag << std::endl;
		//	if(faceFlag == m_PrevFlag)
		//	{

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
							//if((coord[0] == 2.5966) && (coord[1] == -4.63969) && (coord[2] == -4.77856))
							//{
							//	p1ScalarFieldFaces[pointGlobalId] = 1.0;
							//}

			                  		p1ScalarFieldFaces[pointGlobalId] = 1.0;

			           			}
			          	      }
			           }

			        }


	        	 if(numPointsOnFace >= 2)
                  	 {

			     //std::cout << "This is faceFlag before changing: " ;
		             //face.Marker::showMe(std::cout);
		             //std::cout << "This is faceFlag before changing: " << faceFlag << std::endl;
                             face.setMarkerID(m_patchFlag);
			     auto faceFlagChanged = face.markerID();
			     //std::cout << "This is faceFlagChanged in modify patch Area: " << faceFlagChanged << std::endl;
			     //std::cout << "This is faceFlag after changing: ";
			     //face.Marker::showMe(std::cout);

                  	 }
		//	
			
		}
				
				 solver.bcInterfacePtr()->handler()->bcUpdate( *solver.structuralOperatorPtr()->dispFESpacePtr()->mesh(), solver.structuralOperatorPtr()->dispFESpacePtr()->feBd(), solver.structuralOperatorPtr()->dispFESpacePtr()->dof() ); //this is a test if we can change flags over and over

			            m_patchFacesLocationPtr.reset (new vector_Type (p2FeSpace->map() ));
			            *m_patchFacesLocationPtr = p2FeSpace->feToFEInterpolate(p1FESpace, p1ScalarFieldFaces);
		            //*m_patchFacesLocationPtr = p1ScalarFieldFaces;
	 }


	virtual vectorPtr_Type directionalVectorField (EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver ,const boost::shared_ptr<FESpace<RegionMesh<LinearTetra>, MapEpetra >> dFeSpace, Vector3D& direction, const Real& disp, const Real& time)
		{
				
				//Real dispAdder = 0.005;
				Vector3D coord;
				Vector3D intermediateResult;
				Real yDifference;
				Real edgeSmoother;
				Real displacement;
				bool coordInPatchRange = false;
				bool coordOnPatch = false;
				Real fzeroX = m_xPointShift + m_patchDirection[0]*activationFunction(time);;
				//Real fzeroY = m_vertexEllipse[1];
				Real fzeroY = m_yPointShift;
				Real fzeroZ = m_zPointShift + m_patchDirection[2]*activationFunction(time);;

				//These are coefficients of quadratic equation we need to solve for calculating the distance
				Real a;
				Real b;
				Real c;


				Real lambdaOne=0;
				Real lambdaTwo = 0;

				auto p2PositionVector = p2PositionVectorInitial(dFeSpace, solver);
			//	std::cout << "p2PositionVector is of type: " << typeid(p2PositionVector).name() << std::endl;

				vectorPtr_Type p2PatchDisplacement (new VectorEpetra( dFeSpace->map(), Repeated ));
				*p2PatchDisplacement *= 0.0;
		        auto nCompLocalDof = p2PatchDisplacement->epetraVector().MyLength() / 3;

		        for (int j (0); j < nCompLocalDof; ++j)
		        {
		        	                    // Get coordinates

		             UInt iGID = p2PatchDisplacement->blockMap().GID (j);
		             UInt jGID = p2PatchDisplacement->blockMap().GID (j + nCompLocalDof);
		             UInt kGID = p2PatchDisplacement->blockMap().GID (j + 2 * nCompLocalDof);
			
			
	                 coord[0] = p2PositionVector[iGID];
	                 coord[1] = p2PositionVector[jGID];
	                 coord[2] = p2PositionVector[kGID];
			
			
			/*
			coord[0] = (*m_currentPositionVector)[iGID];
                        coord[1] = (*m_currentPositionVector)[jGID];
                        coord[2] = (*m_currentPositionVector)[kGID];
			*/

	                 coordInPatchRange = coordinatesInPatchRangeCurrent(coord, time);

	                 if (coordInPatchRange == true)
	                 {
	                	 intermediateResult[0] = 0.0;
	                	 intermediateResult[1] = 0.0;
	                	 intermediateResult[2] = 0.0;

	                	 Vector3D xVector;
	                	 xVector[0] = coord[0] - m_xPointShift - m_patchDirection[0]*activationFunction(time);
	                	// xVector[1] = coord[1] - m_vertexEllipse[1];
	                	 xVector[1] = coord[1] - m_yPointShift;
	                	 xVector[2] = coord[2] - m_zPointShift - m_patchDirection[2]*activationFunction(time);

	                	 intermediateResult = matrixVectorMultiplicator(m_Ellipsoid, xVector);
	                	 Real ellipseEquation = xVector.dot(intermediateResult)-1;

	                	 if(ellipseEquation >= 0)
	                	 {
	                		 a = m_patchDirection[0]*(m_Ellipsoid(0, 0)*m_patchDirection[0] + m_Ellipsoid(2,0)*m_patchDirection[2]) + m_patchDirection[2]*(m_Ellipsoid(0,2)*m_patchDirection[0] + m_Ellipsoid(2, 2)*m_patchDirection[2]);
	                		 b = m_patchDirection[0]*(m_Ellipsoid(0, 0)*(coord[0]-fzeroX) + m_Ellipsoid(1,0)*(coord[1] - fzeroY) + m_Ellipsoid(2,0)*(coord[2] - fzeroZ)) + m_patchDirection[2]*(m_Ellipsoid(0,2)*(coord[0] - fzeroX) + m_Ellipsoid(1, 2)*(coord[1] - fzeroY) + m_Ellipsoid(2, 2)*(coord[2] - fzeroZ)) + (coord[0] - fzeroX)*(m_Ellipsoid(0,0)*m_patchDirection[0] + m_Ellipsoid(2,0)*m_patchDirection[2]) + (coord[1] - fzeroY)*(m_Ellipsoid(0,1)*m_patchDirection[0] + m_Ellipsoid(2,1)*m_patchDirection[2]) + (coord[2] - fzeroZ)*(m_Ellipsoid(0,2)*m_patchDirection[0] + m_Ellipsoid(2,2)*m_patchDirection[2]);
	                		 c = (coord[0] - fzeroX)*(m_Ellipsoid(0, 0)*(coord[0] - fzeroX) + m_Ellipsoid(1,0)*(coord[1] - fzeroY) + m_Ellipsoid(2, 0)*(coord[2] - fzeroZ)) + (coord[1] - fzeroY)*(m_Ellipsoid(0, 1)*(coord[0] - fzeroX) + m_Ellipsoid(1, 1)*(coord[1] - fzeroY) + m_Ellipsoid(2,1)*(coord[2] - fzeroZ)) + (coord[2] - fzeroZ)*(m_Ellipsoid(0,2)*(coord[0] - fzeroX) + m_Ellipsoid(1,2)*(coord[1] - fzeroY) + m_Ellipsoid(2,2)*(coord[2] - fzeroZ)) -1.0;

	                		 //now we can solve for lambda one and lambda two

	                		 lambdaOne = (-b + sqrt(std::pow(b,2) - 4*a*c))/(2*a);
	                		 lambdaTwo = (-b - sqrt(std::pow(b,2) - 4*a*c))/(2*a);
					
					/*
	                		 std::cout << "This is value of lambdaOne: " << lambdaOne << std::endl;
	                		 std::cout << "This is value of lambdaTwo: " << lambdaTwo << std::endl;
						
					std::cout << "This is yMax: " << m_yMax << std::endl;
					*/
					 if(coord[1] >= 1.1*m_yMax || coord[1] <= 0.9*m_yMin)
					 {
						
						bool checkRange = false;
						checkRange = coord[1] <= 0.9*m_yMin;
						if(checkRange == true)
						{
							yDifference = coord[1] - 0.9*m_yMin;
						}
						else
						{
							yDifference = coord[1] - 1.1*m_yMax;
						}
								
						edgeSmoother = 1*std::pow(yDifference, 2.0);
												
						displacement = std::abs(lambdaOne) < std::abs(lambdaTwo) ? std::abs(lambdaOne) : std::abs(lambdaTwo);
						displacement = displacement - edgeSmoother;
						if(displacement < 0.0)
						{
							displacement = 0.0;
						}
												
						//displacement = displacement/1.4;
					 }
					 else
					 {

					
	                			displacement = std::abs(lambdaOne) < std::abs(lambdaTwo) ? std::abs(lambdaOne) : std::abs(lambdaTwo);
						//displacement = displacement/1.4; 
					//displacement = 0.004;						
					 }
					//displacement = displacement/3;
					//displacement = displacement/2.5;
					//displacement = displacement/2.2;
					//displacement = displacement/2.1;
					/*
					(*p2PatchDisplacement)[iGID] = displacement*m_patchDirection[0]; //0.0
                        		(*p2PatchDisplacement)[jGID] = 0.0;
                       			(*p2PatchDisplacement)[kGID] = displacement*m_patchDirection[2]; //0.0
					*/

					//We want to write displacement to a file and see how it evolves over time
					//checkDisplacement(solver,coord, displacement, time);
					
					//displacement = 0.3;  //0.08 + m_dispAdder;
	                		 //std::cout << "This is value of displacement: " << displacement << std::endl;
					//checkDisplacement(solver,coord, displacement, time);
	                	 }
	                 }
	                 else
	                 {
	                	 displacement = 0.0;
	                 }
		                 //in the following we need to solve a quadratic equation of the form: a*lambda^2 + b*lambda + c = 0 ; therefore we define first coefficients a, b, c
                 //need to check t;hese

			
			
	                (*p2PatchDisplacement)[iGID] = displacement*m_patchDirection[0]; //0.0
	                (*p2PatchDisplacement)[jGID] = 0.0;
	                (*p2PatchDisplacement)[kGID] =  displacement*m_patchDirection[2]; //0.0
			

		        }

			
			//*m_currentPositionVector += *p2PatchDisplacement;
			//*m_currentDisplacementVector += *p2PatchDisplacement;		


		        return p2PatchDisplacement;
		        //return m_currentDisplacementVector;
		}


	virtual void initialisePositionVector(const boost::shared_ptr<FESpace<RegionMesh<LinearTetra>, MapEpetra >> dFeSpace, EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver)
		{	
			std::cout << "We are in intialisePositionVectorFunction" << std::endl;
			//We probably need to first allocate some memory
			
			m_currentPositionVector.reset (new vector_Type (dFeSpace->map(), Repeated));
			m_currentDisplacementVector.reset (new vector_Type (dFeSpace->map(), Repeated));
			*m_currentDisplacementVector *= 0.0;

			
			*m_currentPositionVector = p2PositionVectorInitial(dFeSpace, solver);

			std::cout << "When we are here we survived the memory allocation of the position Vector" << std::endl;
		}	
	

	void checkDisplacement(EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver,Vector3D coord, Real displacement, const Real time)
	{

		auto currentprocessor = solver.comm()->MyPID();
         	std::ostringstream oss;//this is to convert int to string
         	oss << currentprocessor;

         	std::string path_coordinates = "/cluster/home/pamstad/LIFE5/lifev-env/lifev-em-build/lifev/em/examples/example_EMHeart/coordinatefiles/coordinates_" + oss.str() + ".dat";

		
		std::ofstream writer_one(path_coordinates.c_str(), std::ios_base::app);
	        if(writer_one.is_open()==false)
        	{
               		 std::cout << "error occured while opening the file" << std::endl;
         	}
        	writer_one << coord[0] << "\t\t" << coord[1] << "\t\t" << coord[2] <<  "\t\t"  << displacement  << "\t\t" << time <<  std::endl;
         	writer_one.close();

		

	}

	/*
	void initalisePositionVector()
	{
		
		currentPositionVector = p2

	}
	*/
	void initialiseVertexEllipse()
		{
			//first we calculate phi angle

			if(m_vertexEllipse[0] > 0 && m_vertexEllipse[2] > 0)
			{
				m_Phi = atan(m_vertexEllipse[2]/m_vertexEllipse[0])*180/PI;
			}

			if(m_vertexEllipse[0] < 0 && m_vertexEllipse[2] > 0)
			{
				m_Phi = 180 + atan(m_vertexEllipse[2]/m_vertexEllipse[0])*180/PI;
			}

			if(m_vertexEllipse[0] > 0 && m_vertexEllipse[2] < 0)
			{
				m_Phi = atan(m_vertexEllipse[2]/m_vertexEllipse[0])*180/PI;
			}

			if(m_vertexEllipse[0] < 0 && m_vertexEllipse[2] < 0)
			{
				m_Phi = -180 + atan(m_vertexEllipse[2]/m_vertexEllipse[0])*180/PI;
			}
			if(m_vertexEllipse[0] == 0)
			{
				m_Phi = 90;
			}
			
			std::cout << "This is value of Phi: " << m_Phi << std::endl;
			
			/*
			if(m_normalVector[0] > 0 && m_normalVector[1] > 0)
			{
				m_Psi = atan(m_normalVector[1]/m_normalVector[0])*180/PI;
			}

			if(m_normalVector[0] > 0 && m_normalVector[1] < 0)
			{
				m_Psi = atan(m_normalVector[1]/m_normalVector[0])*180/PI;
			}

			if(m_normalVector[0] < 0 && m_normalVector[1] > 0)
			{
				m_Psi = 180 + atan(m_normalVector[1]/m_normalVector[0])*180/PI;
			}

			if(m_normalVector[0] < 0 && m_normalVector[1] < 0)
			{
				m_Psi = -180 + atan(m_normalVector[1]/m_normalVector[0])*180/PI;
			}
			*/

			//m_Psi = -10;
			//m_Psi = -15;
			//m_Psi = -16;
			//m_Psi = -14;
			//m_Psi = -20;
			//m_Psi = -37;
			//m_Psi = -25;
			//m_Psi = -17;
			//m_Psi = -19;
			//m_Psi = -20;
			//m_Phi = -42;
			//m_Phi = -15;
			//next we calculate midpoint
			/*
			m_midPoint[0] = m_vertexEllipse[0] - m_a*std::cos(m_Phi*PI/180);
			m_midPoint[1] = m_vertexEllipse[1] - m_a*std::sin(m_Phi*PI/180);
			m_midPoint[2] = 0.0;
			*/


			//m_shift = m_midPoint.norm();

		}


	void initialiseRotationMatrices()
		{
			//this is only to compute the three different Rotation Matrices and then also the totalRotation matrices
			//we start with the first Rotation

			//////FIRST ROTATION MATRIX
			m_RotationOne (0,0) = std::cos(m_Phi*PI/180.0);
			m_RotationOne (0,1) = 0.0;
			m_RotationOne (0,2) = -std::sin(m_Phi*PI/180.0);

			m_RotationOne (1,0) = 0.0; //std::sin(m_Phi*PI/180.0);
			m_RotationOne (1,1) = 1.0; //std::cos(m_Phi*PI/180.0);
			m_RotationOne (1,2) = 0.0;

			m_RotationOne (2,0) = std::sin(m_Phi*PI/180.0);
			m_RotationOne (2,1) = 0.0;
			m_RotationOne (2,2) = std::cos(m_Phi*PI/180.0);
			//////////////


			//std::cout << "here comes rotationMatrixOne: " << std::endl;
			//printMatrix(m_RotationOne);

			///////SECOND ROTATION MATRIX
			m_RotationTwo (0,0) = std::cos(m_Psi*PI/180.0);
			m_RotationTwo (0,1) = -std::sin(m_Psi*PI/180.0);
			m_RotationTwo (0,2) = 0.0;

			m_RotationTwo (1,0) = std::sin(m_Psi*PI/180.0);
			m_RotationTwo (1,1) = std::cos(m_Psi*PI/180.0);
			m_RotationTwo (1,2) = 0.0;

			m_RotationTwo (2,0) = 0.0;
			m_RotationTwo (2,1) = 0.0;
			m_RotationTwo (2,2) = 1.0;
			/////////////

			//std::cout << "here comes rotationMatrixTwo: " << std::endl;
			//printMatrix(m_RotationTwo);

			m_RotationSum = matrixMatrixMultiplicator(m_RotationTwo, m_RotationOne);

			//std::cout << "here comes rotationMatrixSum: " << std::endl;
			//printMatrix(m_RotationSum);


		}

	void initialiseHeartAxis()
		{
			Vector3D axisPointApex;
			Vector3D axisPointTwo;
			Vector3D axisDirectionVector;
			Vector3D apexToVertex;
			Vector3D perpendicularPoint;
			Real lambda;

			axisPointApex[0] = 1.19793;
			axisPointApex[1] = -13.5271;
			axisPointApex[2] = -1.44197;

			axisPointTwo[0] = -0.195762;
			axisPointTwo[1] = 0.776271;
			axisPointTwo[2] = 0.65295;

			axisDirectionVector = axisPointTwo - axisPointApex;

			//std::cout << "This is axisDirectionVector: " << axisDirectionVector[0] << "      " << axisDirectionVector[1] << "       " << axisDirectionVector[2] << std::endl;

			axisDirectionVector.normalize();

			//std::cout << "This is axisDirectionVector normalised: " << axisDirectionVector[0] << "      " << axisDirectionVector[1] << "       " << axisDirectionVector[2] << std::endl;

			//now we define a vector that goes from Point a apex to the vertexEllipse point;

			apexToVertex = m_vertexEllipse - axisPointApex;

			lambda = apexToVertex.dot(axisDirectionVector); //so lambda is the projection of the apexToVertex Vector on the heart axis

			perpendicularPoint = lambda*axisDirectionVector;

			//changed here the direction of the vector; now it should be right with the moving part
			m_patchDirection = perpendicularPoint - m_vertexEllipse;
			//m_patchDirection = m_vertexEllipse - perpendicularPoint;

			m_patchDirection.normalize();
			//we now have a vector that points from the defined heart axis to the vertex Point of the Ellipse; this gives us direction in which we want to move the patch
			
			std::cout << "This is PatchDirectionVector: " << m_patchDirection[0] << "      " << m_patchDirection[1] << "     " << m_patchDirection[2] << std::endl;
			//std::cout << "This is PatchDirectionVector normalised: " << m_patchDirection[0] << "      " << m_patchDirection[1] << "       " << m_patchDirection[2] << std::endl;


		}

	void calculateyMaxzMax()
		{
			Vector3D pointOne;
			Vector3D pointTwo;
			Vector3D pointThree;
			Vector3D pointFour;
			Vector3D pointFive;
			Vector3D pointSix;
			Vector3D pointSeven;
			Vector3D pointEight;
			Vector3D semiAxisA;
			Vector3D shiftVector;

			//procedure is as described in paper; first we get the heigt of the patch by the user
			if(m_Height <= 2*m_c)
			{
				m_yMax = m_Height/2.0;
			}
			else
			{
				m_yMax == m_c;
				std::cout << "Your given height of the patch exceeds 2*c: check it in dataFile " << std::endl;
			}

			Real f = 1 - std::pow(m_yMax,2.0)/std::pow(m_b,2.0);
			m_zMax = m_Width/2.0;

			if(m_zMax > m_c*sqrt(f))
			{
				//std::cout << "Your given width is too large; check the description and adopt axis (a,b,c) of ellipsoid; z_Max is set to the maximum possible value" << std::endl;
				//m_zMax = m_c*sqrt(f) - 0.01; //0.01 is just to be save that we find solution later
			}


			//m_yMin = -m_yMax + m_vertexEllipse[1];
			//m_yMax = m_yMax + m_vertexEllipse[1];


			m_xValueOneTwo = m_a*sqrt(1 - std::pow(m_zMax, 2.0)/std::pow(m_c, 2.0));//this is x1star in MATLAB file
			m_xValueThreeFour = m_a*sqrt(1 - std::pow(m_yMax, 2.0)/std::pow(m_b, 2.0));
			m_xValueRest = m_a*sqrt(1 - std::pow(m_zMax, 2.0)/std::pow(m_c, 2.0)- std::pow(m_yMax, 2.0)/std::pow(m_b, 2.0));

                        pointThree[0] = m_a;
                        pointThree[1] = 0.0;
                        pointThree[2] = 0.0;
			
                        pointFive[2] = m_zMax;
                        pointFive[1] = 0;
                        pointFive[0] = m_a*sqrt(1-std::pow(pointFive[2], 2.0)/std::pow(m_c, 2.0));

			  pointOne[0] = pointFive[0];
                        pointOne[1] = -m_Height/2;
                        pointOne[2] = m_c*sqrt(1-std::pow(pointOne[0], 2.0)/std::pow(m_a, 2.0) - std::pow(pointOne[1], 2.0)/std::pow(m_b, 2.0));
			
                        pointSeven[2] = -m_zMax;
                        pointSeven[1] = 0.0;
                        pointSeven[0] = m_a*sqrt(1-std::pow(pointSeven[2], 2.0)/std::pow(m_c, 2.0));

			pointTwo[0] = pointSeven[0];
                        pointTwo[1] = -m_Height/2;
                        pointTwo[2] = -m_c*sqrt(1-std::pow(pointTwo[0], 2.0)/std::pow(m_a, 2.0) - std::pow(pointTwo[1], 2.0)/std::pow(m_b, 2.0));


                        pointSix[0] = pointFive[0];
                        pointSix[1] = -m_Height;
                        pointSix[2] = m_c*sqrt(1-std::pow(pointSix[0], 2.0)/std::pow(m_a, 2.0) - std::pow(pointSix[1], 2.0)/std::pow(m_b, 2.0));

			std::cout << "This is argument of point six: " << 1-std::pow(pointSix[0], 2.0)/std::pow(m_a, 2.0) - std::pow(pointSix[1], 2.0)/std::pow(m_b, 2.0) << std::endl;


                        pointEight[0] = pointSeven[0];
                        pointEight[1] = -m_Height;
                        pointEight[2] = -pointSix[2];   // -m_c*sqrt(1-std::pow(pointEight[0], 2.0)/std::pow(m_a, 2.0) - std::pow(pointEight[1], 2.0)/std::pow(m_b, 2.0));
			
                        pointFour[1] = -m_Height;
                        pointFour[2] = 0.0;
			pointFour[0] = m_a*sqrt(1-std::pow(pointFour[1], 2.0)/std::pow(m_b, 2.0));

                        semiAxisA[0] = m_a;
                        semiAxisA[1] = 0.0;
                        semiAxisA[2] = 0.0;
			

			std::cout << "This is larger Radius: " << pointFive[2] - pointSeven[2] << std::endl;
			std::cout << "This is smaller Radius: " << pointSix[2] - pointEight[2] << std::endl;			

			std::cout << "This is point One: " << pointOne[0] << "      " << pointOne[1] << "     " << pointOne[2] << std::endl;
                        std::cout << "This is point Two: " << pointTwo[0] << "      " << pointTwo[1] << "     " << pointTwo[2] << std::endl;
                        std::cout << "This is point Three: " << pointThree[0] << "      " << pointThree[1] << "     " << pointThree[3] << std::endl;
			std::cout << "This is point Four: " << pointFour[0] << "      " << pointFour[1] << "     " << pointFour[2] << std::endl;
			std::cout << "This is point Five: " << pointFive[0] << "      " << pointFive[1] << "     " << pointFive[2] << std::endl;
			std::cout << "This is point Six: " << pointSix[0] << "      " << pointSix[1] << "     " << pointSix[2] << std::endl;
			std::cout << "This is point Seven: " << pointSeven[0] << "      " << pointSeven[1] << "     " << pointSeven[2] << std::endl;
			std::cout << "This is point Eight: " << pointEight[0] << "      " << pointEight[1] << "     " << pointEight[2] << std::endl;

			/*

			pointOne[0] = m_xValueOneTwo;
			pointOne[1] = 0.0;
			pointOne[2] = m_zMax;

			pointTwo[0] = m_xValueOneTwo;
			pointTwo[1] = 0.0;
			pointTwo[2] = -m_zMax;

			pointThree[0] = m_xValueThreeFour;
			pointThree[1] = m_yMax;
			pointThree[2] = 0.0;

			pointFour[0] = m_xValueThreeFour;
			pointFour[1] = -m_yMax;
			pointFour[2] = 0.0;

			pointFive[0] = m_xValueRest;
			pointFive[1] = m_yMax;
			pointFive[2] = m_zMax;

			pointSix[0] = m_xValueRest;
			pointSix[1] = -m_yMax;
			pointSix[2] = m_zMax;

			pointSeven[0] = m_xValueRest;
			pointSeven[1] = m_yMax;
			pointSeven[2] = -m_zMax;

			pointEight[0] = m_xValueRest;
			pointEight[1] = -m_yMax;
			pointEight[2] = -m_zMax;

			semiAxisA[0] = m_a;
			semiAxisA[1] = 0.0;
			semiAxisA[2] = 0.0;

			*/

			m_pointSemiAxesARYZ = matrixVectorMultiplicator(m_RotationSum, semiAxisA);

			m_xPointShift = m_vertexEllipse[0] - m_pointSemiAxesARYZ[0];
			m_yPointShift = m_vertexEllipse[1] - m_pointSemiAxesARYZ[1];
			m_zPointShift = m_vertexEllipse[2] - m_pointSemiAxesARYZ[2];

			shiftVector[0] = m_xPointShift;
			shiftVector[1] = m_yPointShift;
			shiftVector[2] = m_zPointShift;

			m_pointOneRYZ = matrixVectorMultiplicator(m_RotationSum, pointOne) + shiftVector;
			m_pointTwoRYZ = matrixVectorMultiplicator(m_RotationSum, pointTwo) + shiftVector;
			m_pointThreeRYZ = matrixVectorMultiplicator(m_RotationSum, pointThree) + shiftVector;
			m_pointFourRYZ = matrixVectorMultiplicator(m_RotationSum, pointFour) + shiftVector;
			m_pointFiveRYZ = matrixVectorMultiplicator(m_RotationSum, pointFive) + shiftVector;
			m_pointSixRYZ = matrixVectorMultiplicator(m_RotationSum, pointSix) + shiftVector;
			m_pointSevenRYZ = matrixVectorMultiplicator(m_RotationSum, pointSeven) + shiftVector;
			m_pointEightRYZ = matrixVectorMultiplicator(m_RotationSum, pointEight) + shiftVector;

			//we need to determine yMin and yMax of selected points

			Real yValueArray[8] = {m_pointOneRYZ[1], m_pointTwoRYZ[1], m_pointThreeRYZ[1], m_pointFourRYZ[1], m_pointFiveRYZ[1], m_pointSixRYZ[1], m_pointSevenRYZ[1], m_pointEightRYZ[1]};
			Real smallestY = yValueArray[0];
			Real largestY = yValueArray[0];


			for(int i(1); i < 8; i++)
			{
				if(yValueArray[i] < smallestY)
				{
					smallestY = yValueArray[i];
				}
			}

			m_yMin = smallestY;

			for(int j(1); j < 8; j++)
			{
				if(yValueArray[j] > largestY)
				{
					largestY = yValueArray[j];
				}
			}

			m_yMax = largestY;

		}

	void initialiseEllipsoidBoundaries()
		{
			//here comes the definitions of the three slopes and 3 interception values; first we need to find the normalVector from scheitelPoint of ellipse to axis which we want to shift
			Real disp = 15.0;
			Vector3D pointUpperBoundary;
			Vector3D pointUpperBoundaryRotated;
			Vector3D pointLowerBoundary;
			Vector3D pointLowerBoundaryRotated;
			Vector3D pointSix;
			Vector3D pointEight;
			Vector3D pointSixRotated;
			Vector3D pointEightRotated;
			Vector3D vectorFiveSeven; //This is vector from point 5 to 7
			Vector3D vectorFiveSix;
			Vector3D zMaxValues;
			Vector3D zMinValues;
			Real smallestValue;
			Real largestValue;


			vectorFiveSeven = m_pointSevenRYZ - m_pointFiveRYZ;
			vectorFiveSix = m_pointSixRYZ - m_pointFiveRYZ;

			zMaxValues[0] = m_pointOneRYZ[2];
			zMaxValues[1] = m_pointFiveRYZ[2];
			zMaxValues[2] = m_pointSixRYZ[2];

			zMinValues[0] = m_pointTwoRYZ[2];
			zMinValues[1] = m_pointSevenRYZ[2];
			zMinValues[2] = m_pointEightRYZ[2];

			largestValue = zMaxValues[0];
			smallestValue = zMinValues[0];

			for(int i(1); i < 3; i++)
			{
				if(zMaxValues[i] > largestValue)
				{
					largestValue = zMaxValues[i];
				}
			}

			for(int i(1); i < 3; i++)
			{
				if(zMinValues[i] < smallestValue)
				{
					smallestValue = zMinValues[i];
				}
			}


			if(largestValue == m_pointOneRYZ[2])
			{
				m_pointOneUpperBoundary[0] = m_pointOneRYZ[0];
				m_pointOneUpperBoundary[1] = m_pointOneRYZ[1];
				m_pointOneUpperBoundary[2] = m_pointOneRYZ[2];

				m_pointTwoUpperBoundary[0] = m_pointOneUpperBoundary[0] + m_patchDirection[0]*disp;
				m_pointTwoUpperBoundary[1] = m_pointOneUpperBoundary[1];
				m_pointTwoUpperBoundary[2] = m_pointOneUpperBoundary[2] + m_patchDirection[2]*disp;

				m_SlopeOne = (m_pointOneUpperBoundary[2]- m_pointTwoUpperBoundary[2])/(m_pointOneUpperBoundary[0] - m_pointTwoUpperBoundary[0]);
				m_InterceptOne = m_pointOneUpperBoundary[2] - m_SlopeOne*m_pointOneUpperBoundary[0];

			}
			if(largestValue == m_pointFiveRYZ[2])
			{
				m_pointOneUpperBoundary[0] = m_pointFiveRYZ[0];
				m_pointOneUpperBoundary[1] = m_pointFiveRYZ[1];
				m_pointOneUpperBoundary[2] = m_pointFiveRYZ[2];

				m_pointTwoUpperBoundary[0] = m_pointOneUpperBoundary[0] + m_patchDirection[0]*disp;
				m_pointTwoUpperBoundary[1] = m_pointOneUpperBoundary[1];
				m_pointTwoUpperBoundary[2] = m_pointOneUpperBoundary[2] + m_patchDirection[2]*disp;

				m_SlopeOne = (m_pointOneUpperBoundary[2]- m_pointTwoUpperBoundary[2])/(m_pointOneUpperBoundary[0] - m_pointTwoUpperBoundary[0]);
				m_InterceptOne = m_pointOneUpperBoundary[2] - m_SlopeOne*m_pointOneUpperBoundary[0];

			}
			if(largestValue == m_pointSixRYZ[2])
			{
				m_pointOneUpperBoundary[0] = m_pointSixRYZ[0];
				m_pointOneUpperBoundary[1] = m_pointSixRYZ[1];
				m_pointOneUpperBoundary[2] = m_pointSixRYZ[2];

				m_pointTwoUpperBoundary[0] = m_pointOneUpperBoundary[0] + m_patchDirection[0]*disp;
				m_pointTwoUpperBoundary[1] = m_pointOneUpperBoundary[1];
				m_pointTwoUpperBoundary[2] = m_pointOneUpperBoundary[2] + m_patchDirection[2]*disp;

				m_SlopeOne = (m_pointOneUpperBoundary[2]- m_pointTwoUpperBoundary[2])/(m_pointOneUpperBoundary[0] - m_pointTwoUpperBoundary[0]);
				m_InterceptOne = m_pointOneUpperBoundary[2] - m_SlopeOne*m_pointOneUpperBoundary[0];

			}


			if(smallestValue == m_pointTwoRYZ[2])
			{
				m_pointOneLowerBoundary[0] = m_pointTwoRYZ[0];
				m_pointOneLowerBoundary[1] = m_pointTwoRYZ[1];
				m_pointOneLowerBoundary[2] = m_pointTwoRYZ[2];

				m_pointTwoLowerBoundary[0] = m_pointOneLowerBoundary[0] + m_patchDirection[0]*disp;
				m_pointTwoLowerBoundary[1] = m_pointOneLowerBoundary[1];
				m_pointTwoLowerBoundary[2] = m_pointOneLowerBoundary[2] + m_patchDirection[2]*disp;

				m_SlopeTwo = (m_pointOneLowerBoundary[2]- m_pointTwoLowerBoundary[2])/(m_pointOneLowerBoundary[0] - m_pointTwoLowerBoundary[0]);
				m_InterceptTwo = m_pointOneLowerBoundary[2] - m_SlopeTwo*m_pointOneLowerBoundary[0];

			}

			if(smallestValue == m_pointSevenRYZ[2])
			{
				m_pointOneLowerBoundary[0] = m_pointSevenRYZ[0];
				m_pointOneLowerBoundary[1] = m_pointSevenRYZ[1];
				m_pointOneLowerBoundary[2] = m_pointSevenRYZ[2];

				m_pointTwoLowerBoundary[0] = m_pointOneLowerBoundary[0] + m_patchDirection[0]*disp;
				m_pointTwoLowerBoundary[1] = m_pointOneLowerBoundary[1];
				m_pointTwoLowerBoundary[2] = m_pointOneLowerBoundary[2] + m_patchDirection[2]*disp;

				m_SlopeTwo = (m_pointOneLowerBoundary[2]- m_pointTwoLowerBoundary[2])/(m_pointOneLowerBoundary[0] - m_pointTwoLowerBoundary[0]);
				m_InterceptTwo = m_pointOneLowerBoundary[2] - m_SlopeTwo*m_pointOneLowerBoundary[0];

			}

			if(smallestValue == m_pointEightRYZ[2])
			{
					m_pointOneLowerBoundary[0] = m_pointEightRYZ[0];
					m_pointOneLowerBoundary[1] = m_pointEightRYZ[1];
					m_pointOneLowerBoundary[2] = m_pointEightRYZ[2];

					m_pointTwoLowerBoundary[0] = m_pointOneLowerBoundary[0] + m_patchDirection[0]*disp;
					m_pointTwoLowerBoundary[1] = m_pointOneLowerBoundary[1];
					m_pointTwoLowerBoundary[2] = m_pointOneLowerBoundary[2] + m_patchDirection[2]*disp;

					m_SlopeTwo = (m_pointOneLowerBoundary[2]- m_pointTwoLowerBoundary[2])/(m_pointOneLowerBoundary[0] - m_pointTwoLowerBoundary[0]);
					m_InterceptTwo = m_pointOneLowerBoundary[2] - m_SlopeTwo*m_pointOneLowerBoundary[0];

			}

			m_boundaryNormalVector = vectorFiveSeven.cross(vectorFiveSix);
			m_dPlane = -m_boundaryNormalVector[0]*m_pointFiveRYZ[0] - m_boundaryNormalVector[1]*m_pointFiveRYZ[1] - m_boundaryNormalVector[2]*m_pointFiveRYZ[2];

			std::cout << "This is boundaryNormalVector: " << m_boundaryNormalVector[0] << "        " << m_boundaryNormalVector[1] << "        " << m_boundaryNormalVector[2] << std::endl;

		}


	bool coordinatesInPatchRange(const Vector3D& coord, const Real& time)
		{
			bool coordInPatchRange = false;

			if(m_Phi > 0 && m_Phi < 90)
			{
				if(coord[2] - m_SlopeOne*coord[0] - m_InterceptOne <= 0 && coord[2] - m_SlopeTwo*coord[0] - m_InterceptTwo >= 0 && m_boundaryNormalVector[0]*(coord[0] - m_patchDirection[0]*m_maxDisplacement) + m_boundaryNormalVector[1]*coord[1] + m_boundaryNormalVector[2]*(coord[2] - m_patchDirection[2]*m_maxDisplacement) + m_dPlane <= 0 && coord[1] <= m_yMax && coord[1] >= m_yMin)
				{
					coordInPatchRange = true;
				}
			}

			if(m_Phi > 90 && m_Phi < 180)
			{
				if(coord[2] - m_SlopeOne*coord[0] - m_InterceptOne >= 0 && coord[2] - m_SlopeTwo*coord[0] - m_InterceptTwo <= 0 && m_boundaryNormalVector[0]*(coord[0] - m_patchDirection[0]*m_maxDisplacement) + m_boundaryNormalVector[1]*coord[1] + m_boundaryNormalVector[2]*(coord[2] - m_patchDirection[2]*m_maxDisplacement) + m_dPlane <=  0 && coord[1] <= m_yMax && coord[1] >= m_yMin)
				{
					coordInPatchRange = true;
				}
			}

			if(m_Phi > -90 && m_Phi < 0)
			{
				if(coord[2] - m_SlopeOne*coord[0] - m_InterceptOne <= 0 && coord[2] - m_SlopeTwo*coord[0] - m_InterceptTwo >= 0 && m_boundaryNormalVector[0]*(coord[0] - m_patchDirection[0]*m_maxDisplacement) + m_boundaryNormalVector[1]*coord[1] + m_boundaryNormalVector[2]*(coord[2] - m_patchDirection[2]*m_maxDisplacement) + m_dPlane <=  0 && coord[1] <= m_yMax && coord[1] >= m_yMin)
				{
					coordInPatchRange = true;
				}
			}

			if(m_Phi > -180 && m_Phi < -90)
			{
				if(coord[2] - m_SlopeOne*coord[0] - m_InterceptOne >= 0 && coord[2] - m_SlopeTwo*coord[0] - m_InterceptTwo <= 0 && m_boundaryNormalVector[0]*(coord[0] - m_patchDirection[0]*m_maxDisplacement) + m_boundaryNormalVector[1]*coord[1] + m_boundaryNormalVector[2]*(coord[2] - m_patchDirection[2]*m_maxDisplacement) + m_dPlane <=  0 && coord[1] <= m_yMax && coord[1] >= m_yMin)
				{
					coordInPatchRange = true;
				}
			}

					return coordInPatchRange;
		}



		bool coordinatesInPatchRangeCurrent(const Vector3D& coord, const Real& time)
		{
			bool coordInPatchRange = false;

			if(m_Phi > 0 && m_Phi < 90)
			{
				if(coord[2] - m_SlopeOne*coord[0] - m_InterceptOne <= 0 && coord[2] - m_SlopeTwo*coord[0] - m_InterceptTwo >= 0 && m_boundaryNormalVector[0]*(coord[0] - m_patchDirection[0]*activationFunction(time)) + m_boundaryNormalVector[1]*coord[1] + m_boundaryNormalVector[2]*(coord[2] - m_patchDirection[2]*activationFunction(time)) + m_dPlane <=  0 && coord[1] <= m_yMax && coord[1] >= m_yMin)
				{
					coordInPatchRange = true;
				}
			}

			if(m_Phi > 90 && m_Phi < 180)
			{
				if(coord[2] - m_SlopeOne*coord[0] - m_InterceptOne >= 0 && coord[2] - m_SlopeTwo*coord[0] - m_InterceptTwo <= 0 && m_boundaryNormalVector[0]*(coord[0] - m_patchDirection[0]*activationFunction(time)) + m_boundaryNormalVector[1]*coord[1] + m_boundaryNormalVector[2]*(coord[2] - m_patchDirection[2]*activationFunction(time)) + m_dPlane <=  0 && coord[1] <= m_yMax && coord[1] >= m_yMin)
				{
					coordInPatchRange = true;
				}
			}

			if(m_Phi > -90 && m_Phi < 0)
			{
				if(coord[2] - m_SlopeOne*coord[0] - m_InterceptOne <= 0 && coord[2] - m_SlopeTwo*coord[0] - m_InterceptTwo >= 0 && m_boundaryNormalVector[0]*(coord[0] - m_patchDirection[0]*activationFunction(time)) + m_boundaryNormalVector[1]*coord[1] + m_boundaryNormalVector[2]*(coord[2] - m_patchDirection[2]*activationFunction(time)) + m_dPlane <=  0 && coord[1] <= m_yMax && coord[1] >= m_yMin)
				{
					coordInPatchRange = true;
				}
			}

			if(m_Phi > -180 && m_Phi < -90)
			{
				if(coord[2] - m_SlopeOne*coord[0] - m_InterceptOne >= 0 && coord[2] - m_SlopeTwo*coord[0] - m_InterceptTwo <= 0 && m_boundaryNormalVector[0]*(coord[0] - m_patchDirection[0]*activationFunction(time)) + m_boundaryNormalVector[1]*coord[1] + m_boundaryNormalVector[2]*(coord[2] - m_patchDirection[2]*activationFunction(time)) + m_dPlane <=  0 && coord[1] <= m_yMax && coord[1] >= m_yMin)
				{
					coordInPatchRange = true;
				}

			}



			return coordInPatchRange;
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


				intermediateResult = matrixMatrixMultiplicator(m_A, m_RotationSum.transpose());

				m_Ellipsoid = matrixMatrixMultiplicator(m_RotationSum, intermediateResult);

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

			Real m_Phi;
			Real m_Theta;
			Real m_Psi;

			Real m_a;
			Real m_b;
			Real m_c;

			Real m_Height;
			Real m_Width;

			Real m_yMin;
			Real m_yMax;
			Real m_zMin;
			Real m_zMax;

			Real m_xValueOneTwo;
			Real m_xValueThreeFour;
			Real m_xValueRest;

			Vector3D m_pointOneRYZ;
			Vector3D m_pointTwoRYZ;
			Vector3D m_pointThreeRYZ;
			Vector3D m_pointFourRYZ;
			Vector3D m_pointFiveRYZ;
			Vector3D m_pointSixRYZ;
			Vector3D m_pointSevenRYZ;
			Vector3D m_pointEightRYZ;
			Vector3D m_pointSemiAxesARYZ;
			Vector3D m_pointSemiAxesBRYZ;
			Vector3D m_pointSemiAxesCRYZ;

			VectorEpetra currentPositionVector;			
			Real m_dispAdder = 0.005;

			Real m_xPointShift;
			Real m_yPointShift;
			Real m_zPointShift;


			Vector3D m_pointOneUpperBoundary;
			Vector3D m_pointTwoUpperBoundary;
			Vector3D m_pointOneLowerBoundary;
			Vector3D m_pointTwoLowerBoundary;
			Real m_SlopeOne;
			Real m_InterceptOne;
			Real m_SlopeTwo;
			Real m_InterceptTwo;
			Vector3D m_boundaryNormalVector;
			Real m_dPlane;

			Real m_maxDisplacement;
			Vector3D m_vertexEllipse;
			Vector3D m_normalVector;
			Vector3D m_patchDirection;

			vectorPtr_Type m_currentPositionVector;
			vectorPtr_Type m_currentDisplacementVector;	

			matrixSmall_Type m_RotationOne;
			matrixSmall_Type m_RotationTwo;
			matrixSmall_Type m_RotationSum;

			matrixSmall_Type m_A; //Matrix m_A is just a diagonal Matrix with the Halbachsenwerten a,b ,c squared auf der Diagonalen; also m_A = diag(1/a^2, 1/b^2, 1/c^2)
			matrixSmall_Type m_Ellipsoid;

};


REGISTER(EssentialPatchBC, EssentialPatchBCPatchRealYZRotation);

}
#endif /* ESSENTIALPATCHBCPATCHREALYZROTATION_HPP_ */
