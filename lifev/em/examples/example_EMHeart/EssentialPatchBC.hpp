//
//  EssentialPatchBC.hpp
//  lifev-heart
//
//  Created by Thomas Kummer on 27.09.17.
//  Copyright © 2017 Thomas Kummer. All rights reserved.
//

#ifndef EssentialPatchBC_hpp
#define EssentialPatchBC_hpp

#include <stdio.h>
#include <lifev/em/examples/example_EMHeart/GenericFactory.hpp>
#include <string>
#include <fstream>
#include <iostream>
#include <iterator>
#include <algorithm>
#include "Epetra_BlockMap.h"
#include "Epetra_Vector.h"
#include "Epetra_MpiComm.h"

namespace LifeV
{

class EssentialPatchBC
{
public:
    
    typedef VectorEpetra                                    vector_Type;
    typedef boost::shared_ptr<vector_Type>                  vectorPtr_Type;
    
    typedef BCVector                                        bcVector_Type;
    typedef boost::shared_ptr<bcVector_Type>                bcVectorPtr_Type;
    
    typedef EssentialPatchBC                                super;

    
    EssentialPatchBC(){}
    ~EssentialPatchBC(){}
    //
    virtual void setup(const GetPot& dataFile, const std::string& name) //so here we read from a file and assign some values
    {
    	
        // Patch name
        m_Name = name;
		
	//std::cout << "This is variable m_Name: " << m_Name << std::endl;
		        
        // Epicardium flag
        m_PrevFlag = dataFile ( ("solid/boundary_conditions/" + m_Name + "/flag").c_str(), 0.0 ); //what does zero stand for?
        
        // Patch motion direction   //what do we mean by patch direction? with respect to which coordinate system?
        for ( UInt j (0); j < 3; ++j ) //j (0) is just another way to set j = 0
        {
            m_patchDirection[j] = dataFile ( ("solid/boundary_conditions/" + m_Name + "/direction").c_str(), 0.0, j );
        }
        m_patchDirection.normalize();
        
        // Boundary condition components //I don't get this one here
        UInt componentSize = dataFile.vector_variable_size ( ("solid/boundary_conditions/" + m_Name + "/component").c_str() );
        for ( UInt j (0); j < componentSize; ++j )
        {
            m_patchComponent.push_back( dataFile ( ("solid/boundary_conditions/" + m_Name + "/component").c_str(), 0.0, j ) );
        }
        

        // Patch peak displacement
        m_patchDisplacement = dataFile ( ("solid/boundary_conditions/" + m_Name + "/displacement").c_str(), 1.0 );
        
        // Temporal activation parameter
        m_tmax = dataFile ( "solid/patches/tmax", 0. );
        m_tduration = dataFile ( "solid/patches/tduration", 0. );
    }
    
    void createPatchArea (EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver, const int& newFlag, const Real& time)
    {
    	auto p2FeSpace = solver.electroSolverPtr()->feSpacePtr();
    	auto p2dFeSpace = solver.structuralOperatorPtr()->dispFESpacePtr();
    	FESpace<RegionMesh<LinearTetra>, MapEpetra > p1FESpace (p2FeSpace->mesh(), "P1", 1, p2FeSpace->mesh()->comm());

    	VectorEpetra p1ScalarFieldFaces (p1FESpace.map());
    	p1ScalarFieldFaces *= 0.0;
    	Int p1ScalarFieldFacesDof = p1ScalarFieldFaces.epetraVector().MyLength();
    	//int globalIdArray[p1ScalarFieldFacesDof];

    	 m_patchFlag = newFlag;
    	 const auto& mesh = solver.localMeshPtr();

    	 int globalIdArray[p1ScalarFieldFacesDof];
    	 p1ScalarFieldFaces.blockMap().MyGlobalElements(globalIdArray);


    	 for (int j(0); j < mesh->numBoundaryFacets(); j++) //returns number of boundary facets
    	 {
    	       auto& face = mesh->boundaryFacet(j);
    	       auto faceFlag = face.markerID();
		
		//std::cout << "This is faceFlag in createPatchArea before changing: " << faceFlag << std::endl;
		//std::cout << "This is value of m_PrevFlag (should be 464): " << m_PrevFlag << std::endl;


    	       if (faceFlag == m_PrevFlag)
    	       {
    	            int numPointsOnFace(0);

    	            for (int k(0); k < 3; ++k)
    	            {

    	               	ID pointGlobalId = face.point(k).id();
    	               	auto coord = face.point(k).coordinates();
    	              // 	auto pointInPatch = nodeOnPatch

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


    	 //if (numPointsOnFace > 2) // if there are more than two points on face we execute the if statement; not completly sure here

    	 	      if(numPointsOnFace >= 1)
    	          {
    	                     face.setMarkerID(m_patchFlag);
				auto faceFlagChanged = face.markerID();
				//std::cout << "This is changed faceFlag in createPatchArea: " << faceFlagChanged << std::endl;

    	          }
    	  }
    	 }

    	 m_patchFacesLocationPtr.reset(new vector_Type (p2FeSpace->map() ));
    	 *m_patchFacesLocationPtr = p2FeSpace->feToFEInterpolate(p1FESpace, p1ScalarFieldFaces);

    	 m_patchLocationPtr.reset (new vector_Type (p2FeSpace->map() ));
    }


    //this is my createPatchArea; this is version how it worked for the first time; we had it without setting the face markers but with in apply flag = 464
    /*
    void createPatchArea (EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver, const int& newFlag, const Real& time)
    {


    	auto p2FeSpace = solver.electroSolverPtr()->feSpacePtr();
    	auto p2dFeSpace = solver.structuralOperatorPtr()->dispFESpacePtr();
    	FESpace<RegionMesh<LinearTetra>, MapEpetra > p1FESpace (p2FeSpace->mesh(), "P1", 1, p2FeSpace->mesh()->comm());

    	//create an epetra vector to set it equal to one where it is part of patch
    	VectorEpetra p1ScalarFieldFaces (p1FESpace.map());

    	p1ScalarFieldFaces *= 0.0;

    	Int p1ScalarFieldFacesDof = p1ScalarFieldFaces.epetraVector().MyLength();
    	auto numPointsMesh = p2FeSpace->mesh()->numPoints();

    	int globalIdArray[p1ScalarFieldFacesDof];

        m_patchFlag = newFlag;
        const auto& mesh = solver.localMeshPtr(); // variable mesh which we use later for for loop; we assign a local Mesh pointer to it

		//getPatchRegion(solver, m_patchFlag, time);
	

        p1ScalarFieldFaces.blockMap().MyGlobalElements(globalIdArray);

        // Create patches by changing the markerID (flag) locally
        unsigned int numNodesOnPatch(0); //here we just initalise an unsigned integer variable
        for (int j(0); j < mesh->numBoundaryFacets(); j++) //returns number of boundary facets
        {
             auto& face = mesh->boundaryFacet(j);
             auto faceFlag = face.markerID();

             //std::cout << "THESE ARE FACE FLAGS " << faceFlag << std::endl;

          if (faceFlag == m_PrevFlag)
          {
             int numPointsOnFace(0);

             for (int k(0); k < 3; ++k) //k < 3 was before; this is just a test
             {
                        //auto coord = face.point(k).coordinates();
                 	ID pointGlobalId = face.point(k).id();
                 	auto coord = face.point(k).coordinates();
                 	auto pointInPatch = nodeOnPatch(coord, time);

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

             
//if (numPointsOnFace > 2) // if there are more than two points on face we execute the if statement; not completly sure here

	      if(numPointsOnFace >= 1)
             {
                    face.setMarkerID(m_patchFlag);
                    //std::cout << "This is the set Face flag: " ;
                    //face.Marker::showMe(std::cout);
                    numNodesOnPatch++;
             }

        	}
 }


            m_patchFacesLocationPtr.reset(new vector_Type (p2FeSpace->map() ));
            *m_patchFacesLocationPtr = p2FeSpace->feToFEInterpolate(p1FESpace, p1ScalarFieldFaces);
            //*m_patchFacesLocationPtr = p1ScalarFieldFaces;


            m_patchLocationPtr.reset (new vector_Type (p2FeSpace->map() ));
           // *m_patchLocationPtr = p2FeSpace->feToFEInterpolate(p1FESpace, p1ScalarField);
           //*m_patchLocationPtr = p1ScalarField;
                      // if ( solver.comm()->MyPID() == 0 ) std::cout << "\np2Vec size: " << m_patchLocationPtr->size() << " " << p1ScalarFieldDof;


  }
*/

    //HERE IS A FUNCTION WHICH WRITES COORDINATES TO FILE
    /*
    void getPatchRegion(EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver, const int& newFlag, const Real& time)
	{
    		auto p2FeSpace = solver.electroSolverPtr()->feSpacePtr();
            auto p2dFeSpace = solver.structuralOperatorPtr()->dispFESpacePtr();
            FESpace<RegionMesh<LinearTetra>, MapEpetra > p1FESpace (p2FeSpace->mesh(), "P1", 1, p2FeSpace->mesh()->comm());


            VectorEpetra p1ScalarField (p1FESpace.map());
            p1ScalarField *= 0.0;

            Int p1ScalarFieldDof = p1ScalarField.epetraVector().MyLength();
	


            for (int j (0); j < p1ScalarFieldDof; j++)
            {
                UInt iGID = p1ScalarField.blockMap().GID(j);
                ID local_id = p2FeSpace->mesh()->point(j).localId();
				ID myGID = p2FeSpace->mesh()->point(j).id();
	
                Vector3D coord = p2FeSpace->mesh()->point(local_id).coordinates();
		
                auto currentprocessor = solver.comm()->MyPID();
                std::ostringstream oss;//this is to convert int to string
                oss << currentprocessor;



                std::string path_coordinates = "/cluster/home/pamstad/LIFE5/lifev-env/lifev-em-build/lifev/em/examples/example_EMHeart/coordinatefiles/coordinates_" + oss.str() + ".dat";
                std::string path_globaliGID = "/cluster/home/pamstad/LIFE5/lifev-env/lifev-em-build/lifev/em/examples/example_EMHeart/coordinatefiles/globalGID_" + oss.str() + ".dat";
                std::string path_globalmyGID = "/cluster/home/pamstad/LIFE5/lifev-env/lifev-em-build/lifev/em/examples/example_EMHeart/coordinatefiles/globalmyGID_" + oss.str() + ".dat";
                //Here we want to write to file

                std::ofstream writer_one(path_coordinates.c_str(), std::ios_base::app);
                if(writer_one.is_open()==false)
                {
                	std::cout << "error occured while opening the file" << std::endl;
                }

                writer_one << coord[0] << "\t\t" << coord[1] << "\t\t" << coord[2] << "\t\t"  << myGID  << std::endl;
                writer_one.close();

                if ( nodeOnPatch(coord,time) )
                {

                	std::string path_coordinates_two = "/cluster/home/pamstad/LIFE5/lifev-env/lifev-em-build/lifev/em/examples/example_EMHeart/coordinatefiles/coordinates_onPatch" + oss.str() + ".dat";
                	std::ofstream writer_two(path_coordinates_two.c_str(), std::ios_base::app);
                	if(writer_two.is_open()==false)
                	{
                		std::cout << "error occured while opening the file" << std::endl;
                	}
                	writer_two << coord[0] << "\t\t" << coord[1] << "\t\t" << coord[2] << std::endl;
                	writer_two.close();

                    p1ScalarField[iGID] = 1.0;
                }

                std::ofstream writer_three(path_globaliGID.c_str(), std::ios_base::app);

                if(writer_three.is_open()==false)
                {
                	std::cout << "error occured while opening the file" << std::endl;
                }

                writer_three << iGID <<  std::endl;
                writer_three.close();

                std::ofstream writer_four(path_globalmyGID.c_str(), std::ios_base::app);

                if(writer_four.is_open()==false)
                {
                	std::cout << "error occured while opening the file" << std::endl;
                }

                writer_four << myGID <<  std::endl;
                writer_four.close();

	}
}	
*/
/*
    void getFaceCoordinates(EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver, const int& newFlag, const Real& time)
    {
    	 m_patchFlag = newFlag;
    	 const auto& mesh = solver.localMeshPtr(); // variable mesh which we use later for for loop; we assign a local Mesh pointer to it

    	 const auto& meshfull = solver.fullMeshPtr();
    	 auto numPoints = meshfull->numPoints();

    	 unsigned int numNodesOnPatch(0); //here we just initalise an unsigned integer variable
    	 for (int j(0); j < mesh->numBoundaryFacets(); j++) //returns number of boundary facets
    	 {
    	       auto& face = mesh->boundaryFacet(j); //return boundary facet at index j
    	       auto faceFlag = face.markerID(); //probably some specific value which locates the face An entity marker ID is an integral type that is used to store information about each geometric entity

    	       //if (faceFlag == m_PrevFlag)
    	       //{
    	           int numPointsOnFace(0);

    	           for (int k(0); k < 3; ++k)
    	           {
    	                auto coord = face.point(k).coordinates();
    	                ID pointGlobalId = face.point(k).id();
    	                //unsigned int ID = face.point(k).id();
    	                auto pointInPatch = nodeOnPatch(coord, time); //here we check if the point is on the patch; if yes we set it to true

    	                auto currentprocessor = solver.comm()->MyPID();
    	                std::ostringstream oss;//this is to convert int to string
    	                oss << currentprocessor;
    	                std::string path_coordinates = "/cluster/home/pamstad/LIFE5/lifev-env/lifev-em-build/lifev/em/examples/example_EMHeart/coordinatefiles/coordinates_" + oss.str() + ".dat";

    	                std::ofstream writer_one(path_coordinates.c_str(), std::ios_base::app);
    	                if(writer_one.is_open()==false)
    	                {
    	                   	std::cout << "error occured while opening the file" << std::endl;
    	                }


    	                writer_one << "These are coordinates of points with k value: " << k << ": " <<  coord[0] <<  "\t\t" << coord[1] << "\t\t" << coord[2] << "\t\t" << pointGlobalId << std::endl;
    	                if(k == 2)
    	                {
    	                	if(writer_one.is_open() == false)
    	                	{
    	                		std::cout << "Writer one is closed " << std::endl;
    	                	}
    	                	ID faceId = face.id();
    	                	writer_one << "This is ID of the face: " << faceId << "\n\n";
    	                }
    	                writer_one.close();

    	                if (pointInPatch) //if point is on patch we increase the number of points on face
    	                {
    	                	/*
    						auto currentprocessor = solver.comm()->MyPID();
    				    	std::ostringstream oss;//this is to convert int to string
    				    	oss << currentprocessor;
    				    	std::string path_coordinates = "/cluster/home/pamstad/LIFE5/lifev-env/lifev-em-build/lifev/em/examples/example_EMHeart/coordinatefiles/coordinates_" + oss.str() + ".dat";

    				    	std::ofstream writer_one(path_coordinates.c_str(), std::ios_base::app);
    				    	if(writer_one.is_open()==false)
    				    	{
    				    	   	std::cout << "error occured while opening the file" << std::endl;
    				    	}

    				    	writer_one << coord[0] << "\t\t" << coord[1] << "\t\t" << coord[2] << "\t\t"   << std::endl;
    				    	writer_one.close();

    				        ++numPointsOnFace;
    	                }
    	           }


    	           /*
    	           auto currentprocessor = solver.comm()->MyPID();
    	           std::ostringstream oss;//this is to convert int to string
    	           oss << currentprocessor;
    	           std::string path_coordinates_two = "/cluster/home/pamstad/LIFE5/lifev-env/lifev-em-build/lifev/em/examples/example_EMHeart/coordinatefiles/coordinates_" + oss.str() + ".dat";

    	           std::ofstream writer_two(path_coordinates_two.c_str(), std::ios_base::app);
    	           if(writer_two.is_open()==false)
    	           {
    	               	std::cout << "error occured while opening the file" << std::endl;
    	           }

    	           ID faceId = face.id();
    	           writer_two << "This is id of the face with the three points: " << faceId << "\n\n" ;
    	           writer_two.close();


    	           if (numPointsOnFace > 2) // if there are more than two points on face we execute the if statement; not completly sure here
    	           {
    	                face.setMarkerID(m_patchFlag);
    	                numNodesOnPatch++;
    	           }
    	       //}
    	 }

    }

*/

    /*

    void checkMesh(EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver, const int& newFlag, const Real& time)
    {
    	auto p2FeSpace = solver.electroSolverPtr()->feSpacePtr();
    	auto numPointsMesh = p2FeSpace->mesh()->numPoints();

    	for(int k=0; k < numPointsMesh; k++)
    	{

    		Vector3D coord = p2FeSpace->mesh()->point(k).coordinates();
    		ID myGID = p2FeSpace->mesh()->point(k).id();

    		auto currentprocessor = solver.comm()->MyPID();
    		std::ostringstream oss;//this is to convert int to string
    		oss << currentprocessor;

    		std::string path_coordinates = "/cluster/home/pamstad/LIFE5/lifev-env/lifev-em-build/lifev/em/examples/example_EMHeart/coordinatefiles/coordinates_" + oss.str() + ".dat";

    		std::ofstream writer_one(path_coordinates.c_str(), std::ios_base::app);
    		if(writer_one.is_open()==false)
    		{
    		   	std::cout << "error occured while opening the file" << std::endl;
    		}
    		writer_one << coord[0] << "\t\t" << coord[1] << "\t\t" << coord[2] << "\t\t"  << myGID << std::endl;
    		writer_one.close();
    	}

    }

	*/
/*
    void understand_epetraVector(EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver, const int& newFlag, const Real& time)
    {
    	auto p2FeSpace = solver.electroSolverPtr()->feSpacePtr();
    	auto p2dFeSpace = solver.structuralOperatorPtr()->dispFESpacePtr();
    	//FESpace<RegionMesh<LinearTetra>, MapEpetra > p1dFESpace (p2dFeSpace->mesh(), "P1", 1, p2FeSpace->mesh()->comm());
    	FESpace<RegionMesh<LinearTetra>, MapEpetra > p1dFESpace (p2dFeSpace->mesh(), "P1", 3, p2FeSpace->mesh()->comm());


    	// Create P1 VectorEpetra and set it equal to 1.0 in patch regions
    	VectorEpetra p1ScalarField (p1dFESpace.map(), Repeated);

    	p1ScalarField *= 0.0;


    	//if ( solver.comm()->MyPID() == 0 ) std::cout << "\np1Vec size: " << p1ScalarField.size();


    	Int p1ScalarFieldDof = p1ScalarField.epetraVector().MyLength();
    	int globalIdArray[p1ScalarFieldDof];
    	//std::cout << "This is length of epetraVector" << p1ScalarFieldDof << std::endl;

    	Int p1ScalaLocalDof = p1ScalarFieldDof/3;

    	p1ScalarField.blockMap().MyGlobalElements(globalIdArray);

    	//for(int i(0); i < p1ScalarFieldDof; i++)
    	//{
    		//std::cout << "This is global ID on calling processor: " << globalIdArray[i] << std::endl;
    	//}

    	//p1ScalarField.Print(std::cout);

    	//p1ScalarField.showMe(std::cout);

    	p1ScalarField.blockMap().Print(std::cout);

    	for(int j(0); j < p1ScalaLocalDof; j++)
    	{
    		UInt iGID = p1ScalarField.blockMap().GID(j);
    		UInt jGID = p1ScalarField.blockMap().GID(j+ p1ScalaLocalDof);
    		UInt kGID = p1ScalarField.blockMap().GID(j +2*p1ScalaLocalDof);



    		ID pointGID = p2dFeSpace->mesh()->point(j).id();



    		//std::cout << "" << std::endl;

    		//ID pointGID = p2dFeSpace->mesh()->point(j).id();

    		p1ScalarField[iGID] = p2dFeSpace->mesh()->point (j).x();
    		p1ScalarField[jGID] = p2dFeSpace->mesh()->point (j).y();
    		p1ScalarField[kGID] = p2dFeSpace->mesh()->point (j).z();


    		//std::cout << "THIS is iGID from BLockMap: " << iGID << std::endl;
    		//std::cout << "THIS is jGID from BLockMap: " << jGID << std::endl;
    		//std::cout << "THIS is kGID from BLockMap: " << kGID << std::endl;
    		//std::cout << "This is point GLObal ID: " << pointGID << std::endl;
    		//std::cout << "These are coordinates: " << p2dFeSpace->mesh()->point (j).x() << "\t\t" << p2dFeSpace->mesh()->point (j).y() << "\t\t" << p2dFeSpace->mesh()->point (j).z() << std::endl;
    		//std::cout << "" << std::endl;
    	}

    	for(int j(0); j < p1ScalaLocalDof; j++)
    	    	{
    	    		UInt iGID = p1ScalarField.blockMap().GID (j);
    	    		UInt jGID = p1ScalarField.blockMap().GID (j + p1ScalaLocalDof);
    	    		UInt kGID = p1ScalarField.blockMap().GID (j + 2 * p1ScalaLocalDof);

    	    		auto currentprocessor = p2dFeSpace->mesh()->comm()->MyPID();
    	    		std::ostringstream oss;//this is to convert int to string
    	    		oss << currentprocessor;

    	       		std::string path = "/cluster/home/pamstad/LIFE5/lifev-env/lifev-em-build/lifev/em/examples/example_EMHeart/positionvectorfiles/coordinates_" + oss.str() + ".dat";
    	    		std::ofstream writer(path.c_str(), std::ios_base::app);
    	    		if(writer.is_open()==false)
    	       		{
    	    			std::cout << "error occured while opening the file" << std::endl;
    	    		}
    	       		//writer << iGID << "\t\t"<< jGID << "\t\t" << kGID << "\n" <<p1ScalarField[iGID] << "\t\t" << p1ScalarField[jGID] << "\t\t" << p1ScalarField[kGID] << std::endl;
    	      		writer << p1ScalarField[iGID] << "\t\t" << p1ScalarField[jGID] << "\t\t" << p1ScalarField[kGID] << std::endl;
    	    		writer.close();

    	    	}


    }
*/
/*
    void getNodeOnPatchInfo(EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver,const Vector3D& coord,const Real& time )
    {
    	//This is just that we can write it to a file and see what actually happens

    	 auto currentprocessor = solver.comm()->MyPID();
    	 std::ostringstream oss;//this is to convert int to string
    	 oss << currentprocessor;

    	 std::string path_coordinates = "/cluster/home/pamstad/LIFE5/lifev-env/lifev-em-build/lifev/em/examples/example_EMHeart/coordinatefiles/coordinates_" + oss.str() + ".dat";
    	 std::string path_functionValues = "/cluster/home/pamstad/LIFE5/lifev-env/lifev-em-build/lifev/em/examples/example_EMHeart/functionvalues/functionvalues_" + oss.str() + ".dat";
    	 std::string path_activationValues = "/cluster/home/pamstad/LIFE5/lifev-env/lifev-em-build/lifev/em/examples/example_EMHeart/functionvalues/activationvalue_" + oss.str() + ".dat";


    	 //HERE WE WRITE COORDINATES TO FILE
    	 std::ofstream writer_one(path_coordinates.c_str(), std::ios_base::app);
    	 if(writer_one.is_open()==false)
    	 {
    	    	std::cout << "error occured while opening the file" << std::endl;
    	 }
    	 writer_one << coord[0] << "\t\t" << coord[1] << "\t\t" << coord[2] << std::endl;
    	 writer_one.close();

    	 //HERE WE WRITE FUNCTION VALUES TO FILE
    	 std::ofstream writer_two(path_functionValues.c_str(), std::ios_base::app);
    	 if(writer_two.is_open()==false)
    	 {
    	    	std::cout << "error occured while opening the file" << std::endl;
    	 }
    	 writer_two << normal_vector[0]*coord[0] + normal_vector[1]*coord[1] + normal_vector[2]*coord[2] - normal_vector[0]*starting_point[0]- normal_vector[1]*starting_point[1]-normal_vector[2]*starting_point[2] - activationFunction(time) << std::endl;
    	 writer_two.close();

    	 std::ofstream writer_three(path_activationValues.c_str(), std::ios_base::app);
    	 if(writer_three.is_open()==false)
    	 {
    	     	std::cout << "error occured while opening the file" << std::endl;
    	 }
    	 writer_three << activationFunction(time) << std::endl;
    	 writer_three.close();


    }
*/


    void applyBC(EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver, const GetPot& dataFile)
    {
        auto dFeSpace = solver.structuralOperatorPtr() -> dispFESpacePtr(); //not sure what here is going on
      
        m_dispPtr = solver.structuralOperatorPtr()->displacementPtr(); //m_dispPtr is a shared pointer of type vector_type and vector_type is type definition of VectorEpetra
		
        m_patchDispPtr = directionalVectorField(solver,dFeSpace, m_patchDirection, 1e-10, 0.0);
	//we could set here the m_patchDispPtr equal to zero	

	*m_patchDispPtr *= 0.0; //just that we have no displacement at the beginning	

	//std::cout << "This is m_patchComponent: " << m_patchComponent << std::endl;

        m_patchDispBCPtr = bcVectorPtr_Type( new bcVector_Type( *m_patchDispPtr, dFeSpace -> dof().numTotalDof(), 1 ) );
//	 solver.bcInterfacePtr() -> handler()->addBC (m_Name, 464,  Essential, Component, *m_patchDispBCPtr, m_patchComponent); //This was the version where it worked
        solver.bcInterfacePtr() -> handler()->addBC (m_Name, m_patchFlag,  Essential, Component, *m_patchDispBCPtr, m_patchComponent);
       //  solver.bcInterfacePtr() -> handler()->addBC (m_Name, m_patchFlag,  Essential, Full , *m_patchDispBCPtr, m_patchComponent);
    }
    
    
    void modifyPatchBC(EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver, const Real& time, int& PatchFlag)
    {


	//std::cout << "This is value of time variable: "<< time << std::endl;
	//int adder = 12;
    	//const int constantPatchFlag = PatchFlag;
    	//const int constantPatchFlag;
	//std::string patchNameAdder = std::to_string(adder); //converts double variable time to string
	//m_Name = m_Name + patchNameAdder;
	
	const int currentPatchFlag = PatchFlag;

	/*
	if(PatchFlag == 900 && time != 0)
	{
		m_flagIncreaserOne += 10;
		currentPatchFlag = m_flagIncreaserOne;
	}

	if(PatchFlag == 901 && time != 0)
	{
		m_flagIncreaserTwo += 10;
		currentPatchFlag = m_flagIncreaserTwo;
	}
	*/
	//std::cout << "This is modified PatchName: " << m_Name << std::endl;
	
	//std::cout << "This is patchFlag in modifyPatchBC which we give modifyPatchArea: " << constantPatchFlag << std::endl;
	
	
        auto dFeSpace = solver.structuralOperatorPtr()->dispFESpacePtr();
        
        modifyPatchArea(solver, currentPatchFlag, time);

        Real currentPatchDisp = activationFunction(time) + 1e-3;
        if ( 0 == solver.comm()->MyPID() ) std::cout << "\nEssentialPatchBC: " << m_Name << " displaced by " << currentPatchDisp << " cm";

        m_patchDispPtr = directionalVectorField(solver,dFeSpace, m_patchDirection, currentPatchDisp, time);

        m_patchDispBCPtr.reset( new bcVector_Type( *m_patchDispPtr, dFeSpace->dof().numTotalDof(), 1 ) );
	/*
	if (49.99  <= time && time  <= 50.02)
	{
		 solver.bcInterfacePtr() -> handler()->addBC (m_Name, currentPatchFlag,  Essential, Component, *m_patchDispBCPtr, m_patchComponent);
		 solver.bcInterfacePtr()->handler()->modifyBC(currentPatchFlag, *m_patchDispBCPtr);
		if ( 0 == solver.comm()->MyPID() ) solver.bcInterfacePtr() -> handler() -> showMe();
	}

	if (time > 51)
	{
		std::cout << "We are now modifing the BC which we inserted later" << std::endl;
		solver.bcInterfacePtr()->handler()->modifyBC(currentPatchFlag, *m_patchDispBCPtr);
	}
	*/
	//solver.bcInterfacePtr() -> handler()->addBC (m_Name, currentPatchFlag,  Essential, Component, *m_patchDispBCPtr, m_patchComponent);
	//solver.bcInterfacePtr() -> handler()->addBC (m_Name, m_patchFlag,  Essential, Component, *m_patchDispBCPtr, m_patchComponent);
	// if ( 0 == solver.comm()->MyPID() ) solver.bcInterfacePtr() -> handler() -> showMe();
  	 solver.bcInterfacePtr()->handler()->modifyBC(currentPatchFlag, *m_patchDispBCPtr); //This was the version how it worked
        //solver.bcInterfacePtr()->handler()->modifyBC(m_patchFlag, *m_patchDispBCPtr); //this is old version
       //solver.bcInterfacePtr() -> handler()->addBC (m_Name, m_patchFlag,  Essential, Component, *m_patchDispBCPtr, m_patchComponent);//idea is now that we add everytime a new BC
    }
    
    
    vector_Type patchDisplacement(EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver)
    {
        auto dFeSpace = solver.structuralOperatorPtr()->dispFESpacePtr();
        vector_Type localPatchDisplacement ( dFeSpace->map(), Repeated );
        localPatchDisplacement *= 0.0;

        auto nCompLocalDof = localPatchDisplacement.epetraVector().MyLength() / 3;

        for (int j (0); j < nCompLocalDof; ++j) //don't know what happens in this for loop
        {
            UInt iGID = m_patchDispPtr->blockMap().GID (j);
            UInt jGID = m_patchDispPtr->blockMap().GID (j + nCompLocalDof);
            UInt kGID = m_patchDispPtr->blockMap().GID (j + 2 * nCompLocalDof);

            localPatchDisplacement[iGID] = (*m_dispPtr)[iGID]; // * (*m_patchLocationPtr)[iGID];
            localPatchDisplacement[jGID] = (*m_dispPtr)[jGID]; // * (*m_patchLocationPtr)[iGID];
            localPatchDisplacement[kGID] = (*m_dispPtr)[kGID]; // * (*m_patchLocationPtr)[iGID];
        }

        return localPatchDisplacement;
    }

    vector_Type displayDirectionalVectorField(EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver, const Real& time)
    {

    	Vector3D current_point_on_plane;
    	Real distance;
    	Vector3D normalVector;
    	Vector3D startingPoint;
    	Vector3D direction = normalVector;
    	direction.normalize();


    	startingPoint[0] = -3.76487;
    	startingPoint[1] = -10.6687;
    	startingPoint[2] = -0.36572;

    	normalVector[0] = 0.665647;
    	normalVector[1] = 0.695607;
    	normalVector[2] = -0.270367;

    	//first we want to set up the initial vector
    	auto p2dFeSpace = solver.structuralOperatorPtr()->dispFESpacePtr();
    	  	const auto& meshFull = solver.fullMeshPtr();
    	    	FESpace<RegionMesh<LinearTetra>, MapEpetra > p1dFESpace (p2dFeSpace->mesh(), "P1", 3, p2dFeSpace->mesh()->comm());

    	    	    	VectorEpetra p1PositionVector (p1dFESpace.map());
    	    	    	p1PositionVector *= 0.0;

    	    	    	Int p1nLocalPoints = p1PositionVector.epetraVector().MyLength() / 3;

    	    	    	//std::cout << "This is length of positionVector: " << p1PositionVector.epetraVector().MyLength() << std::endl;

    	    	    	for (int j (0); j < p1nLocalPoints; j++)
    	    	    	{
    	    	    	     UInt iGID = p1PositionVector.blockMap().GID (j);
    	    	    	     UInt jGID = p1PositionVector.blockMap().GID (j + p1nLocalPoints);
    	    	    	     UInt kGID = p1PositionVector.blockMap().GID (j + 2 * p1nLocalPoints);


    	    	    	     //Vector3D coord = meshFull->point(iGID).coordinates();

    	    	    	     p1PositionVector[iGID] = meshFull->point(iGID).x();
    	    	    	     p1PositionVector[jGID] = meshFull->point(iGID).y();
    	    	    	     p1PositionVector[kGID] = meshFull->point(iGID).z();



    	    	    	}

    	    	    	VectorEpetra p2PositionVector ( m_dispPtr->map() );
    	    	    	p2PositionVector = p2dFeSpace->feToFEInterpolate(p1dFESpace, p1PositionVector);

    	    	    	//now we have the positionvector


    	    	    	if(time == 0.0)
    	    	    	{
    	    	    		m_currentPositionVector = vectorPtr_Type (new VectorEpetra( p1dFESpace.map(), Repeated ));
    	    	    	    //*m_currentPositionVector = p2PositionVector;
    	    	    	}


    	    	    	vectorPtr_Type p2PatchDisplacement (new VectorEpetra( p1dFESpace.map(), Repeated ));
    	    	    	auto nCompLocalDof = p2PatchDisplacement->epetraVector().MyLength() / 3;




    	    	    	if(normalVector[2] != 0.0)
    	    	    		        {
    	    	    		        	//In thoughts we set cooridnates x and y equal to zero and solve for z coordinate and store it in current_point_on_plane[0]
    	    	    		        	//here we just added the max value of distance vector (3.5155), let's see how it works
    	    	    		        	current_point_on_plane[2] = (normalVector[0]*startingPoint[0] + normalVector[1]*startingPoint[1] + normalVector[2]*startingPoint[2] +activationFunction(time))/normalVector[2];
    	    	    		        	current_point_on_plane[1] = 0;
    	    	    		        	current_point_on_plane[0] = 0;

    	    	    		        	//std::cout << "This is coordinate of current point on plane" << current_point_on_plane[2] << std::endl;

    	    	    	 	        }
    	    	    	/*
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
							*/


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

    	    	    		                    /*
    	    	    		                    coordinates(0) = (*m_currentPositionVector)[iGID];
    	    	    		                    coordinates(1) = (*m_currentPositionVector)[jGID];
    	    	    		                    coordinates(2) = (*m_currentPositionVector)[kGID];
												*/

    	    	    		                    coordinates(0) = p2PositionVector[iGID];
    	    	    		                    coordinates(1) = p2PositionVector[jGID];
    	    	    		                    coordinates(2) = p2PositionVector[kGID];

    	    	    		                    Vector3D QP; //define here the vector that goes from Q (point on plane) to point P

    	    	    		                    QP = coordinates - current_point_on_plane;

    	    	    		                    if ( solver.comm()->MyPID() == 0 )
    	    	    		                    {
    	    	    		                    	//std::cout << "These are coordinates: " << coordinates(0) << "        " << coordinates(1) << "       " << coordinates(2) << std::endl;
    	    	    		                    	//std::cout << "This is current point on plane " << current_point_on_plane(0) << "         " << current_point_on_plane(1) << "        " << current_point_on_plane(2) << std::endl;
    	    	    		                    }


    	    	    		                    if(QP.dot(direction) <= 0) //here i have change to > 0
    	    	    		                    {
    	    	    		                    			distance = 1.0;
    	    	    		                               	//distance = abs(QP.dot(direction));
    	    	    		                    }
    	    	    		                    else
    	    	    		                      {
    	    	    		                              	distance = 0.0;
    	    	    		                      }

    	    	    		                    Vector3D displacement_vector;
    	    	    		                    displacement_vector[0] = distance*direction[0];
    	    	    		                    displacement_vector[1] = distance*direction[1];
    	    	    		                    displacement_vector[2] = distance*direction[2];

    	    	    		                    if ( solver.comm()->MyPID() == 0 )
    	    	    		                    {
    	    	    		                    	//std::cout << "This is dispalcement vector: " << displacement_vector[0] << "       " << displacement_vector[1] << "          " << displacement_vector[2] << std::endl;
    	    	    		                    }

    	    	    		                    (*p2PatchDisplacement)[iGID] = displacement_vector[0];
    	    	    		                    (*p2PatchDisplacement)[jGID] = displacement_vector[1];
    	    	    		                    (*p2PatchDisplacement)[kGID] = displacement_vector[2];

    	    	    		                    	                    (*m_currentPositionVector)[iGID] = (*m_currentPositionVector)[iGID] + (*p2PatchDisplacement)[iGID];
    	    	    		                    	                    (*m_currentPositionVector)[jGID] = (*m_currentPositionVector)[jGID] + (*p2PatchDisplacement)[jGID];
    	    	    		                    	                    (*m_currentPositionVector)[kGID] = (*m_currentPositionVector)[kGID] + (*p2PatchDisplacement)[kGID];

    	    	    		        }

    	    	    	 return *p2PatchDisplacement;
    }

    vector_Type patchVectorField(EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver, const Real& time)
    {


    	auto p2dFeSpace = solver.structuralOperatorPtr()->dispFESpacePtr();
    	const auto& meshFull = solver.fullMeshPtr();
    	FESpace<RegionMesh<LinearTetra>, MapEpetra > p1dFESpace (p2dFeSpace->mesh(), "P1", 3, p2dFeSpace->mesh()->comm());

    	VectorEpetra p1VectorField (p1dFESpace.map());
    	p1VectorField *= 0.0;

    	Int p1nLocalPoints = p1VectorField.epetraVector().MyLength() / 3;

    	for (int j (0); j < p1nLocalPoints; j++)
    	{
    	     UInt iGID = p1VectorField.blockMap().GID (j);
    	     UInt jGID = p1VectorField.blockMap().GID (j + p1nLocalPoints);
    	     UInt kGID = p1VectorField.blockMap().GID (j + 2 * p1nLocalPoints);


    	     Vector3D coord = meshFull->point(iGID).coordinates();
    	     bool pointInPatch = nodeOnPatch(coord, time);

    	     if(pointInPatch == true)
    	     {
    	    	 p1VectorField[iGID] = 1.0;
    	    	 p1VectorField[jGID] = 0.0; //1.0;
    	    	 p1VectorField[kGID] = 0.0; //1.0;

    	     }



    	}

    	return p1VectorField;


   }

    vector_Type patchLocation()
    {
        return *m_patchLocationPtr;
    }
    
    vector_Type patchFacesLocation()
    {
    	return *m_patchFacesLocationPtr;
    }

protected:
    

    //I think main message of directional vector Field is that we can say which node should move in which direction by how much
    virtual vectorPtr_Type directionalVectorField (EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver,const boost::shared_ptr<FESpace<RegionMesh<LinearTetra>, MapEpetra >> dFeSpace, Vector3D& direction, const Real& disp, const Real& time)
    {

    	//all Epetra_Vector constructors require a map argument that describes the layout of elements on the parallel machine
    	//repeated just indicates that the local vectors are replicated

    	//I'm not sure yet how we can use then the variable vectorfield

    	//via the object dFeSpace we have access to the map; need to take closer look

    	//is the same as int* pointer = new int[n] ; some dzynamic array allocation
    	//the one above can also be written as int* pointer(new int[n])
    	//so basically we have a pointer with name vectorfield of type VectorEpetra (which is a class) that points to the constructor of the class

        vectorPtr_Type vectorField (new VectorEpetra( dFeSpace->map(), Repeated )); //no idea what is happening here, but I think here happens the magic
        auto nCompLocalDof = vectorField->epetraVector().MyLength() / 3; //here we want the length of the vector; don't get it

        //then we look at the length of the Vector; so we look how many elements it has in it and divided it by 3 to get the number of nodes; because I assume at every node we have three degrees of freedom; therefore we divide by 3
        direction.normalize(); //here the input direction vector gets normalized

        direction *= disp; //then the direction vector gets multiplied by the displacement
        
        for (int j (0); j < nCompLocalDof; ++j) //What happens in this for loop?
        {
        	//does GID stand for Group identification? Or I think now that it stands for global index

        	//what I don't get yet is the one with GID

            UInt iGID = vectorField->blockMap().GID (j); //UInt is just unsigned integer 32 bit
            UInt jGID = vectorField->blockMap().GID (j + nCompLocalDof);
            UInt kGID = vectorField->blockMap().GID (j + 2 * nCompLocalDof);


            //so here we dereference our pointer but somehow I don't get it
            (*vectorField)[iGID] = direction[0]; // direction is just a 3D vector with direction[0] is x coordinate of vector; are we telling here the different nodes how much they need to displace? //here I should test that with my own vector programm and see how it goes
            (*vectorField)[jGID] = direction[1];
            (*vectorField)[kGID] = direction[2];
        }

        return vectorField;
    }


    virtual void modifyPatchArea(EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver,const int& newFlag, const Real& time){}
    


    //HERE IS THE OLD p2PositionVectorInitial function which didn't work
    /*
    virtual vector_Type p2PositionVectorInitial(const boost::shared_ptr<FESpace<RegionMesh<LinearTetra>, MapEpetra >> p2dFeSpace) const
    {
        // New P1 Space  ; P1 sind lineare element und vom typ a + bx + cy
        FESpace<RegionMesh<LinearTetra> , MapEpetra > p1dFESpace ( p2dFeSpace->mesh(), "P1", 3, p2dFeSpace->mesh()->comm() );
        
        //std::cout << "THIS IS THE CURRENT COMMUNICATOR: " << p2dFeSpace->mesh()->comm()->MyPID() <<std::endl;

        // Create P1 VectorEpetra
        VectorEpetra p1PositionVector (p1dFESpace.map());
        
        // Fill P1 vector with mesh values
        Int p1nCompLocalDof = p1PositionVector.epetraVector().MyLength() / 3;


        //const char* path = "/cluster/home/pamstad/LIFE5/lifev-env/lifev-em-build/lifev/em/examples/example_EMHeart/coordinates.dat";

        for (int j (0); j < p1nCompLocalDof; j++)
        {
            UInt iGID = p1PositionVector.blockMap().GID (j);
            UInt jGID = p1PositionVector.blockMap().GID (j + p1nCompLocalDof);
            UInt kGID = p1PositionVector.blockMap().GID (j + 2 * p1nCompLocalDof);
            
            //p1PositionVector[iGID] = p2dFeSpace->mesh()->point (iGID).x();
            //p1PositionVector[jGID] = p2dFeSpace->mesh()->point (iGID).y();
            //p1PositionVector[kGID] = p2dFeSpace->mesh()->point (iGID).z();


            std::cout << "THIS IS iGID: " << iGID << std::endl;
            std::cout << "THIS IS jGID: " << jGID << std::endl;
            std::cout << "THIS IS kGID: " << kGID << std::endl;
            std::cout << "" << std::endl;
            //here we want to check if the points have an id or how do we know which point is which?
		
		//Here is should also print out the values of iGID and j to see what happens
            ID test_GID = p2dFeSpace->mesh()->point(iGID).id();
            ID global_id = p2dFeSpace->mesh()->point(j).id();
            ID local_id = p2dFeSpace->mesh()->point(j).localId();

            ID test_local_GID = p2dFeSpace->mesh()->point(iGID).localId();
            ID test_local_j = p2dFeSpace->mesh()->point(global_id).localId();



	   p1PositionVector[iGID] = p2dFeSpace->mesh()->point (iGID).x();
	   p1PositionVector[jGID] = p2dFeSpace->mesh()->point (iGID).y();
	   p1PositionVector[kGID] = p2dFeSpace->mesh()->point (iGID).z();


	   //std::cout << "THESE ARE COORDINATES OF point " << j << " with local index j: " << p2dFeSpace->mesh()->point(local_id).x()<< "       " << p2dFeSpace->mesh()->point(local_id).y() << "      " << p2dFeSpace->mesh()->point(local_id).z() << std::endl;
	   //std::cout << "THESE ARE COORDINATES OF point " << j << " with global index from j: " << p2dFeSpace->mesh()->point(global_id).x()<< "       " << p2dFeSpace->mesh()->point(global_id).y() << "      " << p2dFeSpace->mesh()->point(global_id).z()<< std::endl;
	   //std::cout << "THESE ARE COORDINATES OF point " << j << " with global index from Blockmap iGID: " << p2dFeSpace->mesh()->point(iGID).x()<< "       " << p2dFeSpace->mesh()->point(iGID).y() << "      " << p2dFeSpace->mesh()->point(iGID).z() << std::endl;
	   //std::cout << "" << std::endl;
			//p1PositionVector[iGID] = p2dFeSpace->mesh()->point (test_j).x();
            //p1PositionVector[jGID] = p2dFeSpace->mesh()->point (test_j).y();
            //p1PositionVector[kGID] = p2dFeSpace->mesh()->point (test_j).z();

            /////////////////////////////////////////////////////
            //Here we want to get coordinates of position vector

	        //std::cout << "THIS IS X-COORDINATE of p1PositionVector:  " << p1PositionVector[iGID] << std::endl;
            //std::cout << "THIS IS Y-COORDINATE p1PositionVector:  " << p1PositionVector[jGID] << std::endl;
           //std::cout << "THIS IS Z-COORDINATE p1PositionVector:  " << p1PositionVector[kGID] << std::endl;

	    //here we close the file such that no problems occur

        }



        
        // Interpolate position vector from P1-space to current space
        VectorEpetra p2PositionVector ( m_dispPtr->map() );
        p2PositionVector = p2dFeSpace->feToFEInterpolate(p1dFESpace, p1PositionVector);
        
        return p2PositionVector;
    }
    */
    
    /*
    virtual vector_Type p1PositionVectorInitial(const boost::shared_ptr<FESpace<RegionMesh<LinearTetra>, MapEpetra >> p2dFeSpace) const
        {

        	FESpace<RegionMesh<LinearTetra> , MapEpetra > p1dFESpace ( p2dFeSpace->mesh(), "P1", 3, p2dFeSpace->mesh()->comm() );

        	// Create P1 VectorEpetra
        	VectorEpetra p1PositionVector (p1dFESpace.map());

        	// MyLength returns the local vector length on the calling processor
        	Int p1nCompLocalDof = p1PositionVector.epetraVector().MyLength() / 3;

        	for (int j (0); j < p1nCompLocalDof; j++)
        	        {
        	            UInt iGID = p1PositionVector.blockMap().GID (j);
        	            UInt jGID = p1PositionVector.blockMap().GID (j + p1nCompLocalDof);
        	            UInt kGID = p1PositionVector.blockMap().GID (j + 2 * p1nCompLocalDof);

        	            //Up to this point everything is the same as in the old function; the change comes with the argument we give the point() function

        	            ID local_id = p2dFeSpace->mesh()->point(j).localId(); //Here we get the localId of the point


        	            //Here is now big difference; we access the coordinates of the point with the LOCAL ID and not the GLOBAL ID as before; before the argument of the point function was iGID
        	            p1PositionVector[iGID] = p2dFeSpace->mesh()->point (local_id).x();
        	            p1PositionVector[jGID] = p2dFeSpace->mesh()->point (local_id).y();
        	            p1PositionVector[kGID] = p2dFeSpace->mesh()->point (local_id).z();

        	        }

        		return p1PositionVector;
        }
*/
/*
    virtual vectorPtr_Type p2PositionVectorInitial(const boost::shared_ptr<FESpace<RegionMesh<LinearTetra>, MapEpetra >> p2dFeSpace, EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver) const
        {
        	//auto p2dFeSpace = solver.structuralOperatorPtr()->dispFESpacePtr();
        	const auto& meshFull = solver.fullMeshPtr();
        	FESpace<RegionMesh<LinearTetra>, MapEpetra > p1dFESpace (p2dFeSpace->mesh(), "P1", 3, p2dFeSpace->mesh()->comm());

        	 //vectorPtr_Type vectorField (new VectorEpetra( dFeSpace->map(), Repeated ));
        	    	vectorPtr_Type p1PositionVector (new vector_Type (p1dFESpace.map(), Repeated));
        	    	//p1PositionVector *= 0.0;

        	    	Int p1nLocalPoints = p1PositionVector->epetraVector().MyLength() / 3;

        	    	//std::cout << "This is length of positionVector: " << p1PositionVector.epetraVector().MyLength() << std::endl;

        	    	for (int j (0); j < p1nLocalPoints; j++)
        	    	{
        	    	     UInt iGID = p1PositionVector->blockMap().GID (j);
        	    	     UInt jGID = p1PositionVector->blockMap().GID (j + p1nLocalPoints);
        	    	     UInt kGID = p1PositionVector->blockMap().GID (j + 2 * p1nLocalPoints);


        	    	     //Vector3D coord = meshFull->point(iGID).coordinates();

        	    	     (*p1PositionVector)[iGID] = meshFull->point(iGID).x();
        	    	     (*p1PositionVector)[jGID] = meshFull->point(iGID).y();
        	    	     (*p1PositionVector)[kGID] = meshFull->point(iGID).z();



        	    	}

        	    	VectorEpetra p2PositionVector ( m_dispPtr->map() );
        	    	p2PositionVector = p2dFeSpace->feToFEInterpolate(p1dFESpace, *p1PositionVector);

        	    	//const VectorEpetra p2PositionConst (p2PositionVector);
        	    	return p2PositionVector;
        	    	//return p2PositionConst;
        }
*/

//////////////This is working positionVectorFunction

    virtual vector_Type p2PositionVectorInitial(const boost::shared_ptr<FESpace<RegionMesh<LinearTetra>, MapEpetra >> p2dFeSpace, EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver) const
    {
    	//auto p2dFeSpace = solver.structuralOperatorPtr()->dispFESpacePtr();
    	const auto& meshFull = solver.fullMeshPtr();
    	FESpace<RegionMesh<LinearTetra>, MapEpetra > p1dFESpace (p2dFeSpace->mesh(), "P1", 3, p2dFeSpace->mesh()->comm());

    	    	VectorEpetra p1PositionVector (p1dFESpace.map());
    	    	p1PositionVector *= 0.0;

    	    	Int p1nLocalPoints = p1PositionVector.epetraVector().MyLength() / 3;

    	    	//std::cout << "This is length of positionVector: " << p1PositionVector.epetraVector().MyLength() << std::endl;

    	    	for (int j (0); j < p1nLocalPoints; j++)
    	    	{
    	    	     UInt iGID = p1PositionVector.blockMap().GID (j);
    	    	     UInt jGID = p1PositionVector.blockMap().GID (j + p1nLocalPoints);
    	    	     UInt kGID = p1PositionVector.blockMap().GID (j + 2 * p1nLocalPoints);


    	    	     //Vector3D coord = meshFull->point(iGID).coordinates();

    	    	     p1PositionVector[iGID] = meshFull->point(iGID).x();
    	    	     p1PositionVector[jGID] = meshFull->point(iGID).y();
    	    	     p1PositionVector[kGID] = meshFull->point(iGID).z();



    	    	}

    	    VectorEpetra p2PositionVector ( m_dispPtr->map() );
    	    p2PositionVector = p2dFeSpace->feToFEInterpolate(p1dFESpace, p1PositionVector);

    	    	//const VectorEpetra p2PositionConst (p2PositionVector);
    	    	return p2PositionVector;
    	    	//return p2PositionConst;
    	    	//return p1PositionVector;
    }




    virtual Real activationFunction (const Real& time) const
    {
        Real timeInPeriod = fmod(time - m_tmax + 0.5*m_tduration, 800.);
        bool inPeriod ( timeInPeriod < m_tduration && timeInPeriod > 0);
        Real sinusSquared = std::pow( std::sin(timeInPeriod * PI / m_tduration) , 2 ) * m_patchDisplacement;
        return ( inPeriod ? sinusSquared : 0 );
    }
    

    virtual const bool nodeOnPatch(const Vector3D& coord,const Real& time){} ; //this is my function
    //virtual const bool nodeOnPatch(const Vector3D& coord) const = 0;

    virtual const bool nodeOnPatchCurrent(const Vector3D& coord,const Real& time) {};
	
    virtual void initialisePositionVector(const boost::shared_ptr<FESpace<RegionMesh<LinearTetra>, MapEpetra >> dFeSpace, EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver){};

    
    std::string m_Name;
    unsigned int m_PrevFlag;
    unsigned int m_patchFlag;
    unsigned int m_flagIncreaserOne = 900;
    unsigned int m_flagIncreaserTwo = 901;   

 
    Vector3D m_patchDirection;
    Vector3D normal_vector;
    Vector3D starting_point;
    Real m_patchDisplacement;
    
    std::vector<ID> m_patchComponent;

    vectorPtr_Type m_dispPtr;

    vectorPtr_Type m_patchLocationPtr;
    vectorPtr_Type m_patchFacesLocationPtr;
    vectorPtr_Type m_currentPositionVector;

    vectorPtr_Type m_patchDispPtr;
    bcVectorPtr_Type m_patchDispBCPtr;


    Real m_tmax;
    Real m_tduration;
    
};

    

class EssentialPatchBCHandler
{
public:


    EssentialPatchBCHandler(const std::string& patchListName, const GetPot& dataFile) :
        m_patchListName ("solid/boundary_conditions/" + patchListName),
        m_dataFile (dataFile),
       m_patchNumber (( m_dataFile.vector_variable_size(m_patchListName.c_str()) ))
	{}
    
    ~EssentialPatchBCHandler(){}

    void addPatchBC(EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver, const Real& time)
    {
    	//this is my function to write coordinates to file
    	//getcoordinates(solver);

    	//this is my function to write p1positionvector to file
    	//getp1PositionVector(solver);

    	//getcoordinates_createPatchArea(solver, 900, time);

        for ( UInt i (0) ; i < m_patchNumber ; ++i )
        {
            const std::string patchName = m_dataFile ( m_patchListName.c_str(), " ", i ); //type is for example EssentialPatchBCCircular
            const std::string patchType = m_dataFile ( ("solid/boundary_conditions/" + patchName + "/type").c_str(), "EssentialPatchBCMovingPlane" );

            m_patchBCPtrVec.push_back(CREATE(EssentialPatchBC, patchType));
            m_patchBCPtrVec[i]->setup(m_dataFile, patchName);
            m_patchBCPtrVec[i]->createPatchArea(solver, 900 + i, time);
            //m_patchBCPtrVec[i]->createPatchAreaOld(solver, 900 + i, time);
        }

        if ( solver.comm()->MyPID() == 0 ) std::cout << "\nEssentialPatchBCHandler: " << __FUNCTION__ << " - done" << std::endl;
    }

    void applyPatchBC(EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver)
    {

        m_patchDisplacementVecSumPtr = vectorPtr_Type (new VectorEpetra( solver.structuralOperatorPtr()->dispFESpacePtr()->map(), Repeated ));
        m_patchVecSumPtr = vectorPtr_Type (new VectorEpetra( solver.structuralOperatorPtr()->dispFESpacePtr()->map(), Repeated ));
        m_patchLocationScalarSumPtr = vectorPtr_Type (new VectorEpetra( solver.electroSolverPtr()->potentialPtr()->map(), Repeated ));
        m_patchFacesLocationScalarSumPtr = vectorPtr_Type (new VectorEpetra( solver.electroSolverPtr()->potentialPtr()->map(), Repeated ));
        m_directionVecFieldPtr =  vectorPtr_Type (new VectorEpetra( solver.structuralOperatorPtr()->dispFESpacePtr()->map(), Repeated ));

        for (auto& patch : m_patchBCPtrVec)
        {
            patch->applyBC(solver, m_dataFile);
        }

        updatePatchDisplacementSum(solver);
        updatePatchLocationSum(solver);
        updatePatchFacesLocationSum(solver);
        updatePatchVec(solver, 0.0);
        updateDispDirectionalVec(solver, 0.0);
        
        if ( solver.comm()->MyPID() == 0 ) std::cout << "\nEssentialPatchBCHandler: " << __FUNCTION__ << " - done" << std::endl;
    }

    void modifyPatchBC(EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver, const Real& time)
    {
    	int PatchFlag = 899;
	//std::cout << PatchFlag << std::endl;
        for (auto& patch : m_patchBCPtrVec)
        {
        	PatchFlag += 1;
            patch->modifyPatchBC(solver, time, PatchFlag);
        }

        updatePatchDisplacementSum(solver);
        updatePatchLocationSum(solver);
        updatePatchFacesLocationSum(solver);
        updatePatchVec(solver, time);
        updateDispDirectionalVec(solver, time);
        
        if ( solver.comm()->MyPID() == 0 ) std::cout << "\nEssentialPatchBCHandler: " << __FUNCTION__ << " - done" << std::endl;
    }


    void getcoordinates(EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver)
    {
           	auto p2dFeSpace = solver.structuralOperatorPtr() -> dispFESpacePtr();
           	FESpace<RegionMesh<LinearTetra> , MapEpetra > p1dFESpace ( p2dFeSpace->mesh(), "P1", 3, p2dFeSpace->mesh()->comm() );

           	VectorEpetra p1PositionVector (p1dFESpace.map());

           	        // Fill P1 vector with mesh values
           	        Int p1nCompLocalDof = p1PositionVector.epetraVector().MyLength() / 3;

           	        //This works to get number of processors
           	        //std::cout << "THIS IS NUMBER OF PROCESSORS: " << p2dFeSpace->mesh()->comm()->NumProc() << std::endl;
           	        	//	std::cout << "" << std::endl;

           	        //const char* path = "/cluster/home/pamstad/LIFE5/lifev-env/lifev-em-build/lifev/em/examples/example_EMHeart/coordinates.dat";

           	        for (int j (0); j < p1nCompLocalDof; j++)
           	        {
           	            UInt iGID = p1PositionVector.blockMap().GID (j);
           	            UInt jGID = p1PositionVector.blockMap().GID (j + p1nCompLocalDof);
           	            UInt kGID = p1PositionVector.blockMap().GID (j + 2 * p1nCompLocalDof);

           	            ID local_id = p2dFeSpace->mesh()->point(j).localId();

           	            auto currentprocessor = p2dFeSpace->mesh()->comm()->MyPID();
           	            std::ostringstream oss;//this is to convert int to string
           	            oss << currentprocessor;

           	            std::string path = "/cluster/home/pamstad/LIFE5/lifev-env/lifev-em-build/lifev/em/examples/example_EMHeart/coordinatefiles/coordinates_" + oss.str() + ".dat";
           	            std::ofstream writer(path.c_str(), std::ios_base::app);
           	            if(writer.is_open()==false)
           	            {
           	            	std::cout << "error occured while opening the file" << std::endl;
           	            }
           	            writer << p2dFeSpace->mesh()->point(local_id).x() << "\t\t" << p2dFeSpace->mesh()->point(local_id).y() << "\t\t" << p2dFeSpace->mesh()->point(local_id).z() << std::endl;
           	            writer.close();

           	        }
    }


    void getp1PositionVector(EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver)
    {
    	//here in the first part we just create the p1postionvector; this is identical to part in p2position vector
    	auto p2dFeSpace = solver.structuralOperatorPtr() -> dispFESpacePtr();
    	FESpace<RegionMesh<LinearTetra> , MapEpetra > p1dFESpace ( p2dFeSpace->mesh(), "P1", 3, p2dFeSpace->mesh()->comm() );

    	// Create P1 VectorEpetra
    	VectorEpetra p1PositionVector (p1dFESpace.map());

    	// MyLength returns the local vector length on the calling processor
    	Int p1nCompLocalDof = p1PositionVector.epetraVector().MyLength() / 3;

    	for (int j (0); j < p1nCompLocalDof; j++)
    	{
    		UInt iGID = p1PositionVector.blockMap().GID (j);
    		UInt jGID = p1PositionVector.blockMap().GID (j + p1nCompLocalDof);
    		UInt kGID = p1PositionVector.blockMap().GID (j + 2 * p1nCompLocalDof);

    		ID local_id = p2dFeSpace->mesh()->point(j).localId(); //Here we get the localId of the point


    		p1PositionVector[iGID] = p2dFeSpace->mesh()->point (local_id).x();
    		p1PositionVector[jGID] = p2dFeSpace->mesh()->point (local_id).y();
    		p1PositionVector[kGID] = p2dFeSpace->mesh()->point (local_id).z();

    	}

    	//Now we have created the p1positionvector and now we want to print it out --> similar to getcoordinates
    	for(int j(0); j < p1nCompLocalDof; j++)
    	{
    		UInt iGID = p1PositionVector.blockMap().GID (j);
    		UInt jGID = p1PositionVector.blockMap().GID (j + p1nCompLocalDof);
    		UInt kGID = p1PositionVector.blockMap().GID (j + 2 * p1nCompLocalDof);

    		auto currentprocessor = p2dFeSpace->mesh()->comm()->MyPID();
    		std::ostringstream oss;//this is to convert int to string
    		oss << currentprocessor;

       		std::string path = "/cluster/home/pamstad/LIFE5/lifev-env/lifev-em-build/lifev/em/examples/example_EMHeart/positionvectorfiles/coordinates_" + oss.str() + ".dat";
    		std::ofstream writer(path.c_str(), std::ios_base::app);
    		if(writer.is_open()==false)
       		{
    			std::cout << "error occured while opening the file" << std::endl;
    		}
       		writer << p1PositionVector[iGID] << "\t\t" << p1PositionVector[jGID] << "\t\t" << p1PositionVector[kGID] << std::endl;
      		writer.close();

    	}

    }


    void getcoordinates_createPatchArea(EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver, const int& newFlag, const Real& time)
    {

    	unsigned int m_patchFlag = newFlag; //m_patchFlag is again just an unsigned integer
    	//auto just means that compiler assigns datatzype by itself
    	//const auto& mesh = solver.localMeshPtr(); // variable mesh which we use later for for loop; we assign a local Mesh pointer to it
    	const auto& mesh = solver.fullMeshPtr();


     	// Create patches by changing the markerID (flag) locally
    	unsigned int numNodesOnPatch(0); //here we just initalise an unsigned integer variable
    	for (int j(0); j < mesh->numBoundaryFacets(); j++) //returns number of boundary facets
    	{
    		auto& face = mesh->boundaryFacet(j);

                   for (int k(0); k < 3; ++k)
                    {
  	                      auto coord = face.point(k).coordinates(); //These cooridnates we want to get


  	                      auto currentprocessor = solver.comm()->MyPID();
  	                      std::ostringstream oss;//this is to convert int to string
  	                      oss << currentprocessor;

  	                      std::string path = "/cluster/home/pamstad/LIFE5/lifev-env/lifev-em-build/lifev/em/examples/example_EMHeart/coordinatefiles/coordinates_" + oss.str() + ".dat";
  	                      std::ofstream writer(path.c_str(), std::ios_base::app);
  	                      if(writer.is_open()==false)
  	                      {
  	                          	std::cout << "error occured while opening the file" << std::endl;
  	                      }
  	                      writer << coord[0] << "\t\t" << coord[1] << "\t\t" << coord[2] << std::endl;
  	                      writer.close();

   				          //auto pointInPatch = nodeOnPatch(coord, time); //here we check if the point is on the patch; if yes we set it to true

                    }

    	}


    }


    //void get_nodeOnPatchInfo




    vector_Type& patchDisplacementSum()
    {
        return *m_patchDisplacementVecSumPtr;
    }

    vectorPtr_Type patchDisplacementSumPtr()
    {
        return m_patchDisplacementVecSumPtr;
    }

    vector_Type& patchLocationSum()
    {
        return *m_patchLocationScalarSumPtr;
    }

    vectorPtr_Type patchLocationSumPtr()
    {
        return m_patchLocationScalarSumPtr;
    }

    vectorPtr_Type patchFacesLocationSumPtr()
    {
    	return m_patchFacesLocationScalarSumPtr;
    }

    vectorPtr_Type patchVecSumPtr()
    {
        	return m_patchVecSumPtr;
    }

    vectorPtr_Type directionalVecSumPtr()
    {
            	return m_directionVecFieldPtr;
    }


protected:
    
    void updatePatchDisplacementSum(EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver)
    {
        *m_patchDisplacementVecSumPtr *= 0.0;

        for (auto& patch : m_patchBCPtrVec)
        {
            *m_patchDisplacementVecSumPtr += patch->patchDisplacement(solver);
        }
    }
    
    void updatePatchLocationSum(EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver)
    {
        *m_patchLocationScalarSumPtr *= 0.0;

        for (auto& patch : m_patchBCPtrVec)
        {
            *m_patchLocationScalarSumPtr += patch->patchLocation();
        }
    }
    

    void updatePatchFacesLocationSum(EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver)
    {
    	*m_patchFacesLocationScalarSumPtr *= 0.0;
    	for (auto& patch : m_patchBCPtrVec)
    	{
    	            *m_patchFacesLocationScalarSumPtr += patch->patchFacesLocation();
    	}

    }

    void updatePatchVec(EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver, const Real& time)
    {

    	*m_patchVecSumPtr *= 0.0;
    	for (auto& patch : m_patchBCPtrVec)
    	{
    	    	            *m_patchVecSumPtr += patch->patchVectorField(solver, time);
    	}
    }

    void updateDispDirectionalVec(EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver, const Real& time)
    {
    	*m_directionVecFieldPtr *=0.0;

    	*m_directionVecFieldPtr += m_patchBCPtrVec[0]->displayDirectionalVectorField(solver, time);

    }


    const std::string m_patchListName;
    const GetPot& m_dataFile;
    const int m_patchNumber;

    vectorPtr_Type m_patchLocationScalarSumPtr;
    vectorPtr_Type m_patchDisplacementVecSumPtr;
    vectorPtr_Type m_patchFacesLocationScalarSumPtr;
    vectorPtr_Type m_patchVecSumPtr;
    vectorPtr_Type m_directionVecFieldPtr;

    std::vector<EssentialPatchBC*> m_patchBCPtrVec;

};
    
    
}

#endif /* EssentialPatchBC_hpp */
