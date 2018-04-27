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
    
    EssentialPatchBC(){}
    ~EssentialPatchBC(){}
    
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
                        auto pointInPatch = this->determineWhetherInPatch(coord);
                        
                        if (pointInPatch)
                        {
                            ++numPointsInsidePatch;
                        }
                    }
                    
                    if (numPointsInsidePatch > 2)
                    {
                        face.setMarkerID(newFlag);
                    }
                    
                }
            }
        }
    }
    
public:
    
    virtual const bool determineWhetherInPatch(const Vector3D& coord) = 0;
    
    std::string m_Name;
    unsigned int m_PrevFlag;
    
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

    virtual const bool determineWhetherInPatch(const Vector3D& coord)
    {
        return true;
    }

    setup(const GetPot& datafile, const unsigned int& i)
    {
        m_Name = dataFile ( ( "solid/boundary_conditions/listEssentialPatchBC" ), " ", i );
        m_PrevFlag = dataFile ( ("solid/boundary_conditions/" + patchName + "/flag").c_str(), 0 );
        m_Radius= dataFile ( ("solid/boundary_conditions/" + patchName + "/radius").c_str(), 1.0 );
        
        for ( UInt j (0); j < 3; ++j )
        {
            m_Center[j] = dataFile ( ("solid/boundary_conditions/" + patchName + "/center").c_str(), 0, j );
        }
    }
    
    modifyPatchBC(){};


protected:
    
    Vector3D m_Center;;
    Real m_Radius;;

};

}

#endif /* EssentialPatchBC_hpp */
