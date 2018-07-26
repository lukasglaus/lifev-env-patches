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
        
        for (auto& mesh : solver.mesh())
        {
            for (int j(0); j < mesh->numBoundaryFacets(); j++)
            {
                auto& face = mesh->boundaryFacet(j);
                auto faceFlag = face.markerID();

                if (faceFlag == m_PrevFlag)
                {
                    int numPointsInsidePatch (0);
                    
                    for (int k(0); k < 3; ++k)
                    {
                        auto coord = face.point(k).coordinates();
                        auto pointInPatch = nodeOnPatch(coord);

                        if (pointInPatch)
                        {
                            ++numPointsInsidePatch;
                        }
                    }
                    
                    if (numPointsInsidePatch > 2)
                    {
                        face.setMarkerID(m_patchFlag);
                    }
                    
                }
            }
        }
    }
    
    
    void applyBC(EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver, const GetPot& dataFile)
    {
        if ( solver.comm()->MyPID() == 0 ) std::cout << "Applying " << m_Name << " b.c" << std::endl;

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
        if ( 0 == solver.comm()->MyPID() ) std::cout << "\nPatch " << m_Name << " displacement: " << currentPatchDisp << " cm";

        m_patchDispPtr = directionalVectorField(dFeSpace, m_patchDirection, currentPatchDisp, time);

        m_patchDispBCPtr.reset( new bcVector_Type( *m_patchDispPtr, dFeSpace->dof().numTotalDof(), 1 ) );
        solver.bcInterfacePtr()->handler()->modifyBC(m_patchFlag, *m_patchDispBCPtr);
    }
    
    vector_Type& patchDisplacement()
    {
        return *m_patchDispPtr;
    }
    
    vectorPtr_Type patchDisplacementPtr()
    {
        return m_patchDispPtr;
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
        //m_patchDisplacementSumPtr.reset(new VectorEpetra( solver.structuralOperatorPtr()->dispFESpacePtr()->map(), Repeated ));

        for ( UInt i (0) ; i < m_patchNumber ; ++i )
        {
            const std::string patchName = m_dataFile ( m_patchListName.c_str(), " ", i );
            const std::string patchType = m_dataFile ( ("solid/boundary_conditions/" + patchName + "/type").c_str(), "EssentialPatchBCCircular" );
            m_patchBCPtrVec.push_back(CREATE(EssentialPatchBC, patchType));
            m_patchBCPtrVec[i]->setup(m_dataFile, patchName);
            m_patchBCPtrVec[i]->createPatchArea(solver, 900 + i);
        }
    }

    void applyPatchBC(EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver)
    {
        for (auto& patch : m_patchBCPtrVec)
        {
            patch->applyBC(solver, m_dataFile);
        }
    }

    void modifyPatchBC(EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver, const Real& time)
    {
        for (auto& patch : m_patchBCPtrVec)
        {
            patch->modifyPatchBC(solver, time);
        }
        
        updatePatchDisplacementSumPtr();
    }

    vector_Type& patchDisplacementSum()
    {
        return *m_patchDisplacementSumPtr;
    }
    
    vectorPtr_Type patchDisplacementSumPtr()
    {
        return m_patchDisplacementSumPtr;
    }


private:
    
    void updatePatchDisplacementSumPtr()
    {
        (*m_patchDisplacementSumPtr) *= 0.0;

        for (auto& patch : m_patchBCPtrVec)
        {
            (*m_patchDisplacementSumPtr) += patch->patchDisplacement();
        }
    }
    
    
    const std::string m_patchListName;
    const GetPot& m_dataFile;
    const int m_patchNumber;
    
    vectorPtr_Type m_patchDisplacementSumPtr;
    
    std::vector<EssentialPatchBC*> m_patchBCPtrVec;

};
    
    
}

#endif /* EssentialPatchBC_hpp */
