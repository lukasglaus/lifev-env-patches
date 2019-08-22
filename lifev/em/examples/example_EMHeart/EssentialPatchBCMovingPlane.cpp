/*

 * EssentialPatchBCMovingPlane.cpp
 *
 *  Created on: Apr 5, 2019
 *      Author: pamstad
 */

#include "EssentialPatchBCMovingPlane.h"
#include <stdio.h>
#include <lifev/em/examples/example_EMHeart/EssentialPatchBC.hpp>
#include <lifev/em/examples/example_EMHeart/GenericFactory.hpp>
//#include <string>
#include <cmath>
#include <vector>
#define PI 3.14159265359

namespace LifeV
{



void EssentialPatchBCMovingPlane::setup(const GetPot& dataFile, const std::string& name)
{
	super::setup(dataFile, name);
	
	for (UInt i (0); i < 3 ;++i )
	{

		starting_point[i] = dataFile( ("solid/boundary_conditions/" + m_Name + "/startingpoint").c_str(), 0.0, i );
		

	//This point is just for plane that we have a plane that is at t=0 near the heart
	}

	for (UInt j (0); j < 3; ++j)
	{

		normal_vector[j] = dataFile ( ("solid/boundary_conditions/" + m_Name + "/direction").c_str() , 0.0 , j );

	}

	m_maxDisplacement = dataFile ( ("solid/boundary_conditions/" + m_Name + "/displacement").c_str(), 1.0 );

}



const bool EssentialPatchBCMovingPlane::nodeOnPatch(const Vector3D& coord, const Real& time)
{

	bool nodeInArea = false;

		//m_maxDisplacement = dataFile ( ("solid/boundary_conditions/" + m_Name + "/displacement").c_str(), 1.0 );

		if((normal_vector[0]*coord[0] + normal_vector[1]*coord[1] + normal_vector[2]*coord[2] - normal_vector[0]*starting_point[0]- normal_vector[1]*starting_point[1]-normal_vector[2]*starting_point[2] - m_maxDisplacement) <= 0)
			{
				nodeInArea = true;
			}
			else
			{
				nodeInArea = false;
			}


			return nodeInArea;

/*
	bool nodeInArea = 0;

//as shift we had + 0.35
	if((normal_vector[0]*coord[0] + normal_vector[1]*coord[1] + normal_vector[2]*coord[2] - normal_vector[0]*starting_point[0]- normal_vector[1]*starting_point[1]-normal_vector[2]*starting_point[2] - activationFunction(time)) <= 0)
	{
		nodeInArea = true;
	}
	else
	{
		nodeInArea = false;
	}


	return nodeInArea;
*/
}

const bool EssentialPatchBCMovingPlane::nodeOnPatchCurrent(const Vector3D& coord, const Real& time)
{
	bool nodeInArea = 0;

	//as shift we had + 0.35
		if((normal_vector[0]*coord[0] + normal_vector[1]*coord[1] + normal_vector[2]*coord[2] - normal_vector[0]*starting_point[0]- normal_vector[1]*starting_point[1]-normal_vector[2]*starting_point[2] - activationFunction(time)) <= 0)
		{
			nodeInArea = true;
		}
		else
		{
			nodeInArea = false;
		}


		return nodeInArea;

}


void EssentialPatchBCMovingPlane::modifyPatchArea(EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver,const int& newFlag, const Real& time)
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

	        //std::cout << "This is patchFlag in modify Patch Area: " << m_patchFlag << std::endl;
	        const auto& mesh = solver.localMeshPtr(); // variable mesh which we use later for for loop; we assign a local Mesh pointer to it

			const auto& meshfull = solver.fullMeshPtr();
			//auto numPoints = meshfull->numPoints();

			//getPatchRegion(solver, m_patchFlag, time);

	            // Create patches by changing the markerID (flag) locally
	            unsigned int numNodesOnPatch(0); //here we just initalise an unsigned integer variable
	            for (int j(0); j < mesh->numBoundaryFacets(); j++) //returns number of boundary facets
	                    {
	                         auto& face = mesh->boundaryFacet(j);
	                         auto faceFlag = face.markerID();
	                         //std::cout << "This is face marker ID: " << face.markerID() << std::endl;
	                      //if (faceFlag == m_PrevFlag)
	                      //{
	                         int numPointsOnFace(0);

	                         for (int k(0); k < 3; ++k) //k < 3 was before; this is just a test
	                         {
	                                    //auto coord = face.point(k).coordinates();
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
	                             			//++numPointsOnFace;
	                             			p1ScalarFieldFaces[pointGlobalId] = 1.0;

	                             			}
	                             		}
	                             	}

	                         }
				/*
	                         if (numPointsOnFace >= 1) // if there are more than two points on face we execute the if statement; not completly sure here
	                         {
	                        	 	//std::cout << "We are now changing the faceID" << std::endl;
	                        	 	//std::cout << "" << std::endl;
	                        	 	face.setMarkerID(m_patchFlag);
	                        	 	//std::cout << "This is the set face flag: " ;
	                        	 	//face.Marker::showMe(std::cout);
	                                numNodesOnPatch++;
	                         }
				*/
	                       //}
	                    }


	            m_patchFacesLocationPtr.reset (new vector_Type (p2FeSpace->map() ));
	            *m_patchFacesLocationPtr = p2FeSpace->feToFEInterpolate(p1FESpace, p1ScalarFieldFaces);
	            //*m_patchFacesLocationPtr = p1ScalarFieldFaces;

 }

//this is directional vectorfield for p2 elements

vectorPtr_Type EssentialPatchBCMovingPlane::directionalVectorField(EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver,const boost::shared_ptr<FESpace<RegionMesh<LinearTetra>, MapEpetra >> dFeSpace, Vector3D& direction, const Real& disp, const Real& time)
{
	Vector3D current_point_on_plane;
	Real distance;

		//m_p2currentPositionVector = vectorPtr_Type (new VectorEpetra( dFeSpace->map(), Repeated ));
	//std::cout << "NOW WE ARE IN DIRECTIONAL VECTOR FIELD" << std::endl;


			auto p2PositionVector = p2PositionVectorInitial(dFeSpace, solver);
	/*
		if(time == 0.0)
		{
		m_p2currentPositionVector = vectorPtr_Type (new VectorEpetra( dFeSpace->map(), Repeated ));
		*m_p2currentPositionVector = p2PositionVectorInitial(dFeSpace, solver);
		}
		*/

	        vectorPtr_Type p2PatchDisplacement (new VectorEpetra( dFeSpace->map(), Repeated ));


	        auto nCompLocalDof = p2PatchDisplacement->epetraVector().MyLength() / 3;


	        /*
	        bool sameMap = m_p2currentPositionVector->blockMap().SameAs(p2PatchDisplacement->blockMap());

	        if(sameMap == true)
	        {
	        	std::cout << "The maps are the same" << std::endl;
	        }
	        else
	        {
	        	std::cout << "The maps are not the same" << std::endl;
	        }

	        bool samePoint = m_p2currentPositionVector->blockMap().PointSameAs(p2PatchDisplacement->blockMap());
	        if(samePoint == true)
	        {
	            	std::cout << "The points are the same" << std::endl;
	        }
	        else
	        {
	             	std::cout << "The points are not the same" << std::endl;
	        }
			*/

	        direction.normalize(); //need to be careful; direction and normal_vector aren't the same anymore; after that direction is the normalised normal_vector

	        if(normal_vector[2] != 0.0)
	        {
	        	//In thoughts we set cooridnates x and y equal to zero and solve for z coordinate and store it in current_point_on_plane[0]
	        	//here we just added the max value of distance vector (3.5155), let's see how it works
	        	current_point_on_plane[2] = (normal_vector[0]*starting_point[0] + normal_vector[1]*starting_point[1] + normal_vector[2]*starting_point[2] +activationFunction(time))/normal_vector[2];
	        	current_point_on_plane[1] = 0;
	        	current_point_on_plane[0] = 0;

	        	//std::cout << "This is coordinate of current point on plane" << current_point_on_plane[2] << std::endl;

 	        }
	        else if (normal_vector[2] == 0.0 && normal_vector[1] != 0.0)
	        {
		    	
	        	current_point_on_plane[1] = (normal_vector[0]*starting_point[0] + normal_vector[1]*starting_point[1] + normal_vector[2]*starting_point[2]  + activationFunction(time))/normal_vector[1];
	        	current_point_on_plane[2] = 0;
	        	current_point_on_plane[0] = 0;
	        }
	        else if (normal_vector[2] == 0.0 && normal_vector[1] == 0.0 && normal_vector[0] != 0.0)
	        {
	        	current_point_on_plane[0] = (normal_vector[0]*starting_point[0] + normal_vector[1]*starting_point[1] + normal_vector[2]*starting_point[2]  +activationFunction(time))/normal_vector[0];
	        	current_point_on_plane[1] = 0;
	        	current_point_on_plane[2] = 0;
	        }
	        else
 	        {	
	        	std::cout << "A normal  vector in the data file of (0, 0 , 0) doesn't make sense" << std::endl;
	        }


		//std::cout << "THIS IS CURRENT POINT ON PLANE "  << current_point_on_plane[0] << "       " << current_point_on_plane[1] << "         " << current_point_on_plane[2] << std::endl;
		// std::cout << current_point_on_plane[1] << std::endl;
		// std::cout << current_point_on_plane[2] << std::endl;



	        //std::cout << "This is Length of EpetraDisplacementVector in DirectionalVectorfield: " << p2PatchDisplacement->epetraVector().MyLength() << std::endl;

	        for (int j (0); j < nCompLocalDof; ++j)
	        {
	                    // Get coordinates

	                    UInt iGID = p2PatchDisplacement->blockMap().GID (j);
	                    UInt jGID = p2PatchDisplacement->blockMap().GID (j + nCompLocalDof);
	                    UInt kGID = p2PatchDisplacement->blockMap().GID (j + 2 * nCompLocalDof);

	                    /*
	        			UInt iGID = p2PositionVector.blockMap().GID (j);
	        		    UInt jGID = p2PositionVector.blockMap().GID (j + nCompLocalDof);
	        		    UInt kGID = p2PositionVector.blockMap().GID (j + 2 * nCompLocalDof);
						*/

	                    Vector3D coordinates;

	                    coordinates(0) = p2PositionVector[iGID];
	                    coordinates(1) = p2PositionVector[jGID];
	                    coordinates(2) = p2PositionVector[kGID];

	                    /*
	                    coordinates(0) = (*m_p2currentPositionVector)[iGID];
	                    coordinates(1) = (*m_p2currentPositionVector)[jGID];
	                    coordinates(2) = (*m_p2currentPositionVector)[kGID];
						*/

	                    Vector3D QP; //define here the vector that goes from Q (point on plane) to point P

	                    QP = coordinates - current_point_on_plane;

	                    //std::cout << "THESE ARE COORDINATES OF QP: " << QP[0] << "       " << QP[1] << "       " << QP[2] << std::endl;


	                 //   QP[0] = coordinates[0] - current_point_on_plane[0];
	                 //   QP[1] = coordinates[1] - current_point_on_plane[1];
	                 //   QP[2] = coordinates[2] - current_point_on_plane[2];


	                    //next we do the projection of the vector QP ond the normalvector to the plane; we look at the sign which it has in the if statement

	                    //if(QP[0]*direction[0] + QP[1]*direction[1] + QP[2]*direction[2] <=0) //we use here normalised normal vector
	                    if(QP.dot(direction) <= 0) //here i have change to > 0
						{
	                    	//if the dot product is smaller equal zero, then we want to apply a displacement; for that we calculate the distance from point P to plane and then say this is the displacement we want
			
	                    	distance = abs(QP.dot(direction));
	                    	//distance = abs(QP[0]*direction[0]+ QP[1]*direction[1]+QP[2]*direction[2]);
	                    	//std::cout << " THIS IS CALCULATED DISTANCE IN VECTORFIELD: " << distance << std::endl;


	                    	//////////////////Here we write distance to file
	                    	auto currentprocessor = dFeSpace->mesh()->comm()->MyPID();
	                    	std::ostringstream oss;//this is to convert int to string
	                    	oss << currentprocessor;

	                    	std::string path = "/cluster/home/pamstad/LIFE5/lifev-env/lifev-em-build/lifev/em/examples/example_EMHeart/distancefiles/distances_" + oss.str() + ".dat";
	                    	std::ofstream writer(path.c_str(), std::ios_base::app);
	                    	if(writer.is_open()==false)
	                    	{
	                    		std::cout << "error occured while opening the file" << std::endl;
	                    	}
	                    	writer << coordinates(0) << "\t\t" << coordinates(1) << "\t\t" << coordinates(2) << "\t\t"  << distance << std::endl;
	                    	writer.close();
	                    	//////////////////Here the writing to the file ends


	                    }
	                    else
	                    {
	                    	distance = 0.0;
	                    }

	                    Vector3D displacement_vector;
	                    displacement_vector[0] = 0.0; // distance*direction[0];
	                    displacement_vector[1] = 0.0; //  distance*direction[1];
	                    displacement_vector[2] = 0.0; // distance*direction[2];

				//std::cout << "This is displacmeent VEctor: " <<  displacement_vector(0) << "          " << displacement_vector(1) << "         " << displacement_vector(2) << std::endl;				

	                    (*p2PatchDisplacement)[iGID] = displacement_vector[0];
	                    (*p2PatchDisplacement)[jGID] = displacement_vector[1];
	                    (*p2PatchDisplacement)[kGID] = displacement_vector[2];

	                    /*
	                    (*m_p2currentPositionVector)[iGID] = (*m_p2currentPositionVector)[iGID] + (*p2PatchDisplacement)[iGID];
	                    (*m_p2currentPositionVector)[jGID] = (*m_p2currentPositionVector)[jGID] + (*p2PatchDisplacement)[jGID];
	                    (*m_p2currentPositionVector)[kGID] = (*m_p2currentPositionVector)[kGID] + (*p2PatchDisplacement)[kGID];
						*/
	                    /*
	                    auto currentprocessor = dFeSpace->mesh()->comm()->MyPID();
	                    std::ostringstream oss;//this is to convert int to string
	                    oss << currentprocessor;

	                    std::string path_two = "/cluster/home/pamstad/LIFE5/lifev-env/lifev-em-build/lifev/em/examples/example_EMHeart/currentPositionVector/currentPositionVector_" + oss.str() + ".dat";
	                    std::ofstream writer_two(path_two.c_str(), std::ios_base::app);
	                    if(writer_two.is_open()==false)
	                    {
	                          		std::cout << "error occured while opening the file" << std::endl;
	                    }
	                    writer_two << (*m_p2currentPositionVector)[iGID] << "\t\t" << (*m_p2currentPositionVector)[jGID] << "\t\t" << (*m_p2currentPositionVector)[kGID] << std::endl;
	                    writer_two.close();
						*/


	        }
	        //p2PositionVector += *p2PatchDisplacement;
	        //*m_p2currentPositionVector += *p2PatchDisplacement;
		
		if(time == 0)
		{
			*p2PatchDisplacement *= 0.0;
		}

	        return p2PatchDisplacement;

}

//THIS IS THE CORRECTED VERSION OF THE DIRECTIONAL VECTORFIELD; LETS SEE
/*
vectorPtr_Type EssentialPatchBCMovingPlane::directionalVectorField(const boost::shared_ptr<FESpace<RegionMesh<LinearTetra>, MapEpetra >> dFeSpace, Vector3D& direction, const Real& disp, const Real& time) const
{
	Vector3D current_point_on_plane;
	Real distance;

	// auto p2PositionVector = p2PositionVectorDisplaced(dFeSpace);
	        auto p2PositionVector = p2PositionVectorInitial(dFeSpace);

	        auto nCompLocalDof_two = p2PositionVector.epetraVector().MyLength() / 3;

	        //Here we allocate memory for the epetra vector
	        vectorPtr_Type p2PatchDisplacement (new VectorEpetra( dFeSpace->map(), Repeated ));

	        auto nCompLocalDof = p2PatchDisplacement->epetraVector().MyLength() / 3;

	        direction.normalize(); //need to be careful; direction and normal_vector aren't the same anymore; after that, direction is the normalised normal_vector

	        //the next steps we do is to get a point on the plane; we call that point Q; with that point we can calculate a vector QP to any point P of the heart
	        //that vector QP we can later project onto the normalised normalvector of the plane; with that we can do two things:
	        // 1) can determine if left or right to the plane ; 2) can also determine the distance from point to plane which is distance for displacement

	        if(normal_vector[2] != 0)
	        {
	        	//In thoughts we set cooridnates x and y equal to zero and solve for z coordinate and store it in current_point_on_plane[0]
	        	//here we just added the max value of distance vector (3.5155), let's see how it works
	        	current_point_on_plane[2] = (normal_vector[0]*starting_point[0] + normal_vector[1]*starting_point[1] + normal_vector[2]*starting_point[2]-0.5 +activationFunction(time))/normal_vector[2];
	        	current_point_on_plane[1] = 0;
	        	current_point_on_plane[0] = 0;

	        	//std::cout << "This is coordinate of current point on plane" << current_point_on_plane[2] << std::endl;

 	        }
	        else if (normal_vector[2] == 0 && normal_vector[1] != 0)
	        {

	        	current_point_on_plane[1] = (normal_vector[0]*starting_point[0] + normal_vector[1]*starting_point[1] + normal_vector[2]*starting_point[2] - 0.5 + activationFunction(time))/normal_vector[1];
	        	current_point_on_plane[2] = 0;
	        	current_point_on_plane[0] = 0;
	        }
	        else if (normal_vector[2] == 0 && normal_vector[1] == 0 && normal_vector[0] != 0)
	        {
	        	current_point_on_plane[0] = (normal_vector[0]*starting_point[0] + normal_vector[1]*starting_point[1] + normal_vector[2]*starting_point[2] - 0.5  +activationFunction(time))/normal_vector[0];
	        	current_point_on_plane[1] = 0;
	        	current_point_on_plane[2] = 0;
	        }
	        else
 	        {
	        	std::cout << "A normal  vector in the data file of (0, 0 , 0) doesn't make sense" << std::endl;
	        }

	        //now we have a point on our plane which is stored in the variable current_point_on_Plane

	        //In the next steps we want to get the coordinates of the points of the heart

		    for (int j (0); j < nCompLocalDof; ++j)
	        {
	                    // Get coordinates
	                    UInt iGID = p2PatchDisplacement->blockMap().GID (j);
	                    UInt jGID = p2PatchDisplacement->blockMap().GID (j + nCompLocalDof);
	                    UInt kGID = p2PatchDisplacement->blockMap().GID (j + 2 * nCompLocalDof);

	                    //I think the coordinates are correct; what will be the bigger problem is to get right GID to fill the epetra p2PatchDisplacement EpetraVector
	                    Vector3D coordinates;
	                    coordinates(0) = p2PositionVector[iGID];
	                    coordinates(1) = p2PositionVector[jGID];
	                    coordinates(2) = p2PositionVector[kGID];


	                    //std::cout << "THESE ARE Coordinates OF NODES: " << coordinates(0) << "\t" << coordinates(1) << "\t" << coordinates(2) << std::endl;
	                    //we have here now the coordinates of of localDOF j; we call this point P

	                    Vector3D QP; //define here the vector that goes from Q (point on plane) to point P


	                    QP = coordinates - current_point_on_plane;

	                 //   QP[0] = coordinates[0] - current_point_on_plane[0];
	                 //   QP[1] = coordinates[1] - current_point_on_plane[1];
	                 //   QP[2] = coordinates[2] - current_point_on_plane[2];


	                    //next we do the projection of the vector QP ond the normalvector to the plane; we look at the sign which it has in the if statement

	                    //if(QP[0]*direction[0] + QP[1]*direction[1] + QP[2]*direction[2] <=0) //we use here normalised normal vector
	                    if(QP.dot(direction) < 0) //here i have change <= to <
						{
	                    	//if the dot product is smaller equal zero, then we want to apply a displacement; for that we calculate the distance from point P to plane and then say this is the displacement we want

	                    	distance = abs(QP[0]*direction[0]+ QP[1]*direction[1]+QP[2]*direction[2]);
	                    	//std::cout << " THIS IS CALCULATED DISTANCE IN VECTORFIELD: " << distance << std::endl;


	                    	//////////////////Here we write distance to file
	                    	auto currentprocessor = dFeSpace->mesh()->comm()->MyPID();
	                    	std::ostringstream oss;//this is to convert int to string
	                    	oss << currentprocessor;

	                    	std::string path = "/cluster/home/pamstad/LIFE5/lifev-env/lifev-em-build/lifev/em/examples/example_EMHeart/distancefiles/distances_" + oss.str() + ".dat";
	                    	std::ofstream writer(path.c_str(), std::ios_base::app);
	                    	if(writer.is_open()==false)
	                    	{
	                    		std::cout << "error occured while opening the file" << std::endl;
	                    	}
	                    	writer << distance << std::endl;
	                    	writer.close();
	                    	//////////////////Here the writing to the file ends


	                    }
	                    else
	                    {
	                    	distance = 0;
	                    }

	                    Vector3D displacement_vector;
	                    displacement_vector[0] = distance*direction[0];
	                    displacement_vector[1] = distance*direction[1];
	                    displacement_vector[2] = distance*direction[2];


	                    (*p2PatchDisplacement)[iGID] = displacement_vector[0];
	                    (*p2PatchDisplacement)[jGID] = displacement_vector[1];
	                    (*p2PatchDisplacement)[kGID] = displacement_vector[2];

	        }


	        return p2PatchDisplacement;

}
*/

//REGISTER(EssentialPatchBC, EssentialPatchBCMovingPlane);

}//this is Klammer von LifeV namespace
