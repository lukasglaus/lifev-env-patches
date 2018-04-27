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

    void createPatchArea (EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver, const int& newFlag) const
    {
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
                        std::cout << coord;
                        auto pointInPatch = determineWhetherInPatch(coord);
                        std::cout << pointInPatch;
                        if (pointInPatch)
                        {
                            std::cout << "P";

                            ++numPointsInsidePatch;
                        }
                    }
                    
                    if (numPointsInsidePatch > 2)
                    {
                        std::cout << "O";
                        face.setMarkerID(newFlag);
                    }
                    
                }
            }
        }
    }
    
//    void applyBC()
//    {
//        patchDisplacement.push_back( dataFile ( ("solid/boundary_conditions/" + m_Name + "/displacement").c_str(), 1.0 ) );
//
//        for ( UInt j (0); j < 3; ++j )
//        {
//            patchDirection[j] = dataFile ( ("solid/boundary_conditions/" + m_Name + "/direction").c_str(), 0, j );
//        }
//
//        UInt componentSize = dataFile.vector_variable_size ( ("solid/boundary_conditions/" + m_Name + "/component").c_str() );
//        std::vector<ID> patchComponent (componentSize);
//        for ( UInt j (0); j < componentSize; ++j )
//        {
//            patchComponent[j] = dataFile ( ("solid/boundary_conditions/" + m_Name + "/component").c_str(), 0, j );
//        }
//
//        patchDispPtr = heartSolver.directionalVectorField(FESpace, patchDirection[i], 1e-10);
//
//        patchDispBCPtr = bcVectorPtr_Type( new bcVector_Type( *patchDispVecPtr[i], solver.structuralOperatorPtr() -> dispFESpacePtr() -> dof().numTotalDof(), 1 ) ) );
//        solver.bcInterfacePtr() -> handler()->addBC (patchName, (900+i),  Essential, Component, *patchDispBCVecPtr[i], patchComponent);
//    }
    
    
protected:
    
    //heartSolver.directionalVectorField
    
    virtual const bool determineWhetherInPatch(Vector3D& coord) const = 0;
    
    std::string m_Name;
    unsigned int m_PrevFlag;
    
    Real patchDisplacement;
    Vector3D patchDirection;
    
    vectorPtr_Type patchDispPtr;
    bcVectorPtr_Type patchDispBCPtr;
    
};


class EssentialPatchBCCircular : public EssentialPatchBC
{
public:
    
    EssentialPatchBCCircular(){}
    ~EssentialPatchBCCircular(){}
    

    boost::shared_ptr<VectorEpetra> directionalVectorField (const boost::shared_ptr<FESpace<RegionMesh<LinearTetra>, MapEpetra >> dFeSpace, Vector3D& direction, const Real& disp) const
    {
        boost::shared_ptr<VectorEpetra> vectorField (new VectorEpetra( dFeSpace->map(), Repeated ));
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

    
    Real sinSquared (const Real& time, const Real& Tmax, const Real& tmax, const Real& tduration) const
    {
        Real timeInPeriod = fmod(time-tmax+0.5*tduration, 800.);
        bool inPeriod ( timeInPeriod < tduration && timeInPeriod > 0);
        Real sinusSquared = std::pow( std::sin(timeInPeriod * PI / tduration) , 2 ) * Tmax;
        return ( inPeriod ? sinusSquared : 0 );
    }

    virtual void setup(const GetPot& dataFile, const std::string& name)
    {
        m_Name = name;
        m_PrevFlag = dataFile ( ("solid/boundary_conditions/" + m_Name + "/flag").c_str(), 0 );
        m_Radius= dataFile ( ("solid/boundary_conditions/" + m_Name + "/radius").c_str(), 1.0 );
        
        for ( UInt j (0); j < 3; ++j )
        {
            m_Center[j] = dataFile ( ("solid/boundary_conditions/" + m_Name + "/center").c_str(), 0, j );
        }
    }
    
    
    modifyPatchBC(){};


protected:
    
    virtual const bool determineWhetherInPatch(Vector3D& coord) const
    {
        bool pointInCircle = (coord - m_Center).norm() < m_Radius;
        return pointInCircle;
    }

    std::string m_Name;
    unsigned int m_PrevFlag;
    
    Vector3D m_Center;;
    Real m_Radius;;

};
REGISTER(EssentialPatchBC, EssentialPatchBCCircular);

}

#endif /* EssentialPatchBC_hpp */
