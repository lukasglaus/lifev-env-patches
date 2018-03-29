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




class EssentialPatchBCCircular
{
public:
    
    EssentialPatchBCCircular()
    {
        
    }
    ~EssentialPatchBCCircular(){}
    

    void createPatch (EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > >& solver, const Vector3D& center, const Real& radius, const int& currentFlag, const int& newFlag) const
    {
        for (auto& mesh : solver.mesh())
        {
            for (int j(0); j < mesh->numBoundaryFacets(); j++)
            {
                auto& face = mesh->boundaryFacet(j);
                auto faceFlag = face.markerID();
                
                if (faceFlag == currentFlag || faceFlag == 470 || faceFlag == 471)
                {
                    int numPointsInsidePatch (0);
                    
                    for (int k(0); k < 3; ++k)
                    {
                        auto coord = face.point(k).coordinates();
                        bool pointInPatch = (coord - center).norm() < radius;
                        
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

    
    setupPatchBC(const GetPot& datafile)
    {
        UInt nDispPatchBC = dataFile.vector_variable_size ( ( "solid/boundary_conditions/listEssentialPatchBC" ) );
        
        for ( UInt i (0) ; i < nDispPatchBC ; ++i )
        {
            std::string patchName = dataFile ( ( "solid/boundary_conditions/listEssentialPatchBC" ), " ", i );
            Real patchFlag = dataFile ( ("solid/boundary_conditions/" + patchName + "/flag").c_str(), 0 );
            
            Real patchRadius = dataFile ( ("solid/boundary_conditions/" + patchName + "/radius").c_str(), 1.0 );
            Vector3D patchCenter;
            for ( UInt j (0); j < 3; ++j )
            {
                patchCenter[j] = dataFile ( ("solid/boundary_conditions/" + patchName + "/center").c_str(), 0, j );
            }
            
            createPatch(solver, patchCenter, patchRadius, patchFlag, (900+i));
        }
    }
    
    modifyPatchBC();


protected:

    std::string m_name;
    unsigned int m_prevFaceFlag;
    unsigned int m_currentPatchFlag;
    
    // BCFunctionBase m_bcFunctionBase;
    // BCFunctionDirectional m_bcFunctionDirectional;
    
    // PatchBCFunctionBaseCreator m_patchBCFunctionBaseCreator;
    
    Vector3D m_center { 0. , 0. , 0. };
    Real m_radius { 0. };

};


#endif /* EssentialPatchBC_hpp */
