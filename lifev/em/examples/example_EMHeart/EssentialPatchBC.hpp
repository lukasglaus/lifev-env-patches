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

    
    EssentialPatchBC(){}
    ~EssentialPatchBC(){}
    
    virtual void setup(const GetPot& dataFile, const std::string& name) = 0;

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

        UInt componentSize = dataFile.vector_variable_size ( ("solid/boundary_conditions/" + m_Name + "/component").c_str() );
        std::vector<ID> patchComponent (componentSize);
        for ( UInt j (0); j < componentSize; ++j )
        {
            patchComponent[j] = dataFile ( ("solid/boundary_conditions/" + m_Name + "/component").c_str(), 0, j );
        }

        Vector3D currentPatchDirection = assistDirection();
        m_patchDispPtr = directionalVectorField(dFeSpace, currentPatchDirection, 1e-10, 0.0);

        m_patchDispBCPtr = bcVectorPtr_Type( new bcVector_Type( *m_patchDispPtr, dFeSpace -> dof().numTotalDof(), 1 ) );
        solver.bcInterfacePtr() -> handler()->addBC (m_Name, m_patchFlag,  Essential, Component, *m_patchDispBCPtr, patchComponent);
    }
    
    
    void modifyPatchBC(EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver, const Real& time)
    {
        auto dFeSpace = solver.structuralOperatorPtr() -> dispFESpacePtr();
        
        Real currentPatchDisp = activationFunction(time) + 1e-6;
        Vector3D currentPatchDirection = assistDirection();

        m_patchDispPtr = directionalVectorField(dFeSpace, currentPatchDirection, currentPatchDisp, time);
        if ( 0 == solver.comm()->MyPID() ) std::cout << "\nCurrent " << m_Name << " displacement: " << currentPatchDisp << " cm";

        m_patchDispBCPtr.reset( new bcVector_Type( *m_patchDispPtr, dFeSpace->dof().numTotalDof(), 1 ) );
        solver.bcInterfacePtr()->handler()->modifyBC(m_patchFlag, *m_patchDispBCPtr);
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
    
    virtual const bool nodeOnPatch(Vector3D& coord) const = 0;
    
    virtual Real activationFunction (const Real& time) const = 0;

    virtual Vector3D assistDirection () const = 0;

    
    std::string m_Name;
    unsigned int m_PrevFlag;
    unsigned int m_patchFlag;
    
    Vector3D m_patchDirection;

    vectorPtr_Type m_dispPtr;
    
    vectorPtr_Type m_patchDispPtr;
    bcVectorPtr_Type m_patchDispBCPtr;
    
};


//class EssentialPatchBCHandler
//{
//public:
//
//    EssentialPatchBCHandler(){}
//    ~EssentialPatchBCHandler(){}
//
//    addPatch(EssentialPatchBCHandler* patch)
//    {
//        m_patchBCVec.push_back(patch);
//    }
//
//    setup(){}
//
//    applyBC(){}
//
//    modifyPatch(){}
//
//
//private:
//
//    std::vector<EssentialPatchBC*> m_patchBCVec
//
//};
    
    
}

#endif /* EssentialPatchBC_hpp */
