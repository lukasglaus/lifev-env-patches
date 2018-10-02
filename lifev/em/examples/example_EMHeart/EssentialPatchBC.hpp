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
    
    virtual void setup(const GetPot& dataFile, const std::string& name)
    {
        // Patch name
        m_Name = name;
        
        // Epicardium flag
        m_PrevFlag = dataFile ( ("solid/boundary_conditions/" + m_Name + "/flag").c_str(), 0 );
        
        // Patch motion direction
        for ( UInt j (0); j < 3; ++j )
        {
            m_patchDirection[j] = dataFile ( ("solid/boundary_conditions/" + m_Name + "/direction").c_str(), 0, j );
        }
        m_patchDirection.normalize();
        
        // Boundary condition components
        UInt componentSize = dataFile.vector_variable_size ( ("solid/boundary_conditions/" + m_Name + "/component").c_str() );
        for ( UInt j (0); j < componentSize; ++j )
        {
            m_patchComponent.push_back( dataFile ( ("solid/boundary_conditions/" + m_Name + "/component").c_str(), 0, j ) );
        }
        
        // Patch peak displacement
        m_patchDisplacement = dataFile ( ("solid/boundary_conditions/" + m_Name + "/displacement").c_str(), 1.0 );
        
        // Temporal activation parameter
        m_tmax = dataFile ( "solid/patches/tmax", 0. );
        m_tduration = dataFile ( "solid/patches/tduration", 0. );
    }
    
    void createPatchArea (EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver, const int& newFlag)
    {
        m_patchFlag = newFlag;
        const auto& mesh = solver.localMeshPtr();

        // Create patches by changing the markerID (flag) locally
        unsigned int numNodesOnPatch(0);
        for (int j(0); j < mesh->numBoundaryFacets(); j++)
        {
            auto& face = mesh->boundaryFacet(j);
            auto faceFlag = face.markerID();

            if (faceFlag == m_PrevFlag)
            {
                int numPointsOnFace(0);
                
                for (int k(0); k < 3; ++k)
                {
                    auto coord = face.point(k).coordinates();
                    auto pointInPatch = nodeOnPatch(coord);

                    if (pointInPatch)
                    {
                        ++numPointsOnFace;
                    }
                }
                
                if (numPointsOnFace > 2)
                {
                    face.setMarkerID(m_patchFlag);
                    numNodesOnPatch++;
                }
            }
        }
        
        if ( solver.comm()->MyPID() == 0 ) std::cout << "\nEssentialPatchBC: " << __FUNCTION__ << ": " << numNodesOnPatch << " nodes on patch";

        // Setup P1-space
        auto p2FeSpace = solver.electroSolverPtr()->feSpacePtr();
        auto p2dFeSpace = solver.structuralOperatorPtr()->dispFESpacePtr();
        FESpace<RegionMesh<LinearTetra>, MapEpetra > p1FESpace (p2FeSpace->mesh(), "P1", 1, p2FeSpace->mesh()->comm());
        
        // Create P1 VectorEpetra and set it equal to 1.0 in patch regions
        VectorEpetra p1ScalarField (p1FESpace.map());
        p1ScalarField *= 0.0;
        
        if ( solver.comm()->MyPID() == 0 ) std::cout << "\np1Vec size: " << p1ScalarField.size();
        
        Int p1ScalarFieldDof = p1ScalarField.epetraVector().MyLength() / 3;
        for (int j (0); j < p1ScalarFieldDof; j++)
        {
            UInt iGID = p1ScalarField.blockMap().GID(j);
            
            Vector3D coord = p1FESpace.mesh()->point(iGID).coordinates();
            //if ( nodeOnPatch(coord) )
            {
                p1ScalarField[iGID] = 1.0;
            }
        }
        
        // Interpolation from P1-space to P2-space
        m_patchLocationPtr.reset (new vector_Type (p2FeSpace->map() ));
        *m_patchLocationPtr = p2FeSpace->feToFEInterpolate(p1FESpace, p1ScalarField);
        
        if ( solver.comm()->MyPID() == 0 ) std::cout << "\np2Vec size: " << m_patchLocationPtr->size();

    }
    
    
    void applyBC(EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver, const GetPot& dataFile)
    {
        auto dFeSpace = solver.structuralOperatorPtr() -> dispFESpacePtr();
        m_dispPtr = solver.structuralOperatorPtr()->displacementPtr();

        m_patchDispPtr = directionalVectorField(dFeSpace, m_patchDirection, 1e-10, 0.0);

        m_patchDispBCPtr = bcVectorPtr_Type( new bcVector_Type( *m_patchDispPtr, dFeSpace -> dof().numTotalDof(), 1 ) );
        solver.bcInterfacePtr() -> handler()->addBC (m_Name, m_patchFlag,  Essential, Component, *m_patchDispBCPtr, m_patchComponent);
    }
    
    
    void modifyPatchBC(EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver, const Real& time)
    {
        auto dFeSpace = solver.structuralOperatorPtr()->dispFESpacePtr();
        
        Real currentPatchDisp = activationFunction(time) + 1e-3;
        if ( 0 == solver.comm()->MyPID() ) std::cout << "\nEssentialPatchBC: " << m_Name << " displaced by " << currentPatchDisp << " cm";

        m_patchDispPtr = directionalVectorField(dFeSpace, m_patchDirection, currentPatchDisp, time);

        m_patchDispBCPtr.reset( new bcVector_Type( *m_patchDispPtr, dFeSpace->dof().numTotalDof(), 1 ) );
        solver.bcInterfacePtr()->handler()->modifyBC(m_patchFlag, *m_patchDispBCPtr);
    }
    
    
    vector_Type patchDisplacement(EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver)
    {
        auto dFeSpace = solver.structuralOperatorPtr()->dispFESpacePtr();
        vector_Type localPatchDisplacement ( dFeSpace->map(), Repeated );
        localPatchDisplacement *= 0.0;

        auto nCompLocalDof = m_patchDispPtr->epetraVector().MyLength() / 3;

        for (int j (0); j < nCompLocalDof; ++j)
        {
            UInt iGID = m_patchDispPtr->blockMap().GID (j);
            UInt jGID = m_patchDispPtr->blockMap().GID (j + nCompLocalDof);
            UInt kGID = m_patchDispPtr->blockMap().GID (j + 2 * nCompLocalDof);

            localPatchDisplacement[iGID] = (*m_patchDispPtr)[iGID] * (*m_patchLocationPtr)[iGID];
            localPatchDisplacement[jGID] = (*m_patchDispPtr)[jGID] * (*m_patchLocationPtr)[iGID];
            localPatchDisplacement[kGID] = (*m_patchDispPtr)[kGID] * (*m_patchLocationPtr)[iGID];
        }

        return localPatchDisplacement;
    }


    vector_Type patchLocation()
    {
        return *m_patchLocationPtr;
    }

    
protected:
    
    virtual vectorPtr_Type directionalVectorField (const boost::shared_ptr<FESpace<RegionMesh<LinearTetra>, MapEpetra >> dFeSpace, Vector3D& direction, const Real& disp, const Real& time) const
    {
        vectorPtr_Type vectorField (new VectorEpetra( dFeSpace->map(), Repeated ));
        auto nCompLocalDof = vectorField->epetraVector().MyLength() / 3;

        direction.normalize();
        
        direction *= disp;
        for (int j (0); j < nCompLocalDof; ++j)
        {
            UInt iGID = vectorField->blockMap().GID (j);
            UInt jGID = vectorField->blockMap().GID (j + nCompLocalDof);
            UInt kGID = vectorField->blockMap().GID (j + 2 * nCompLocalDof);

            (*vectorField)[iGID] = direction[0];
            (*vectorField)[jGID] = direction[1];
            (*vectorField)[kGID] = direction[2];
        }

        return vectorField;
    }
    
    
    virtual vector_Type p2PositionVectorInitial(const boost::shared_ptr<FESpace<RegionMesh<LinearTetra>, MapEpetra >> p2dFeSpace) const
    {
        // New P1 Space
        FESpace<RegionMesh<LinearTetra> , MapEpetra > p1dFESpace ( p2dFeSpace->mesh(), "P1", 3, p2dFeSpace->mesh()->comm() );
        
        // Create P1 VectorEpetra
        VectorEpetra p1PositionVector (p1dFESpace.map());
        
        // Fill P1 vector with mesh values
        Int p1nCompLocalDof = p1PositionVector.epetraVector().MyLength() / 3;
        for (int j (0); j < p1nCompLocalDof; j++)
        {
            UInt iGID = p1PositionVector.blockMap().GID (j);
            UInt jGID = p1PositionVector.blockMap().GID (j + p1nCompLocalDof);
            UInt kGID = p1PositionVector.blockMap().GID (j + 2 * p1nCompLocalDof);
            
            p1PositionVector[iGID] = p2dFeSpace->mesh()->point (iGID).x();
            p1PositionVector[jGID] = p2dFeSpace->mesh()->point (iGID).y();
            p1PositionVector[kGID] = p2dFeSpace->mesh()->point (iGID).z();
        }
        
        // Interpolate position vector from P1-space to current space
        VectorEpetra p2PositionVector ( m_dispPtr->map() );
        p2PositionVector = p2dFeSpace->feToFEInterpolate(p1dFESpace, p1PositionVector);
        
        return p2PositionVector;
    }
    
    
    virtual Real activationFunction (const Real& time) const
    {
        Real timeInPeriod = fmod(time - m_tmax + 0.5*m_tduration, 800.);
        bool inPeriod ( timeInPeriod < m_tduration && timeInPeriod > 0);
        Real sinusSquared = std::pow( std::sin(timeInPeriod * PI / m_tduration) , 2 ) * m_patchDisplacement;
        return ( inPeriod ? sinusSquared : 0 );
    }
    
    virtual const bool nodeOnPatch(const Vector3D& coord) const = 0;

    
    std::string m_Name;
    unsigned int m_PrevFlag;
    unsigned int m_patchFlag;
    
    Vector3D m_patchDirection;
    Real m_patchDisplacement;
    
    std::vector<ID> m_patchComponent;

    vectorPtr_Type m_dispPtr;

    vectorPtr_Type m_patchLocationPtr;

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

    void addPatchBC(EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver)
    {
        for ( UInt i (0) ; i < m_patchNumber ; ++i )
        {
            const std::string patchName = m_dataFile ( m_patchListName.c_str(), " ", i );
            const std::string patchType = m_dataFile ( ("solid/boundary_conditions/" + patchName + "/type").c_str(), "EssentialPatchBCCircular" );
            m_patchBCPtrVec.push_back(CREATE(EssentialPatchBC, patchType));
            m_patchBCPtrVec[i]->setup(m_dataFile, patchName);
            m_patchBCPtrVec[i]->createPatchArea(solver, 900 + i);
        }
        
        if ( solver.comm()->MyPID() == 0 ) std::cout << "\nEssentialPatchBCHandler: " << __FUNCTION__ << " - done" << std::endl;
    }

    void applyPatchBC(EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver)
    {
        m_patchDisplacementVecSumPtr = vectorPtr_Type (new VectorEpetra( solver.structuralOperatorPtr()->dispFESpacePtr()->map(), Repeated ));
        m_patchLocationScalarSumPtr = vectorPtr_Type (new VectorEpetra( solver.electroSolverPtr()->potentialPtr()->map(), Repeated ));

        for (auto& patch : m_patchBCPtrVec)
        {
            patch->applyBC(solver, m_dataFile);
        }

        updatePatchDisplacementSum(solver);
        updatePatchLocationSum(solver);
        
        if ( solver.comm()->MyPID() == 0 ) std::cout << "\nEssentialPatchBCHandler: " << __FUNCTION__ << " - done" << std::endl;
    }

    void modifyPatchBC(EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver, const Real& time)
    {
        for (auto& patch : m_patchBCPtrVec)
        {
            patch->modifyPatchBC(solver, time);
        }

        updatePatchDisplacementSum(solver);
        updatePatchLocationSum(solver);
        
        if ( solver.comm()->MyPID() == 0 ) std::cout << "\nEssentialPatchBCHandler: " << __FUNCTION__ << " - done" << std::endl;
    }

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


private:
    
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
    
    const std::string m_patchListName;
    const GetPot& m_dataFile;
    const int m_patchNumber;

    vectorPtr_Type m_patchLocationScalarSumPtr;
    vectorPtr_Type m_patchDisplacementVecSumPtr;

    std::vector<EssentialPatchBC*> m_patchBCPtrVec;

};
    
    
}

#endif /* EssentialPatchBC_hpp */
