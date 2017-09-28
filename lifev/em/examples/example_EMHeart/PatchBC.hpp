//
//  PatchBC.hpp
//  lifev-heart
//
//  Created by Thomas Kummer on 27.09.17.
//  Copyright © 2017 Thomas Kummer. All rights reserved.
//

#ifndef PatchBC_hpp
#define PatchBC_hpp

#include <stdio.h>

#endif /* PatchBC_hpp */


namespace LifeV
{

class PatchBC
{
public:
    
    typedef EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > > EMSolverType;
    typedef boost::function<Real ( const Real&, const Real&, const Real&, const Real&, const ID& ) > function_Type;
    
    PatchBC(EMSolverType& solver, const std::string& bcName, const int& prevFaceFlag, const int& patchFlag) :
        m_solver (solver),
        m_bcName (bcName),
        m_prevFaceFlag (prevFaceFlag),
        m_patchFlag (patchFlag)
    {}
    
    void setup(BCFunctionBase& bcFunctionBase, const Vector3D& center, const Real& radius)
    {
        setShapeParameters(center, radius);
        setBCFunctionBase(bcFunctionBase);
        createPatchArea();
        addPatchBC();
    }
    
protected:
    
    void setShapeParameters(const Vector3D& center, const Real& radius)
    {
        m_center = center;
        m_radius = radius;
    }
    
    void setBCFunctionBase(BCFunctionBase& bcFunctionBase)
    {
        m_bcFunctionBase.setFunction(bcFunctionBase);
    }
    
    virtual void createPatchArea()
    {
        for (auto& mesh : m_solver.mesh())
        {
            for (int j(0); j < mesh->numBoundaryFacets(); j++)
            {
                auto& face = mesh->boundaryFacet(j);
                auto faceFlag = face.markerID();
                
                if (faceFlag == m_prevFaceFlag || faceFlag == 470 || faceFlag == 471)
                {
                    int numPointsInsidePatch (0);
                    
                    for (int k(0); k < 3; ++k)
                    {
                        auto coord = face.point(k).coordinates();
                        bool pointInPatch = (coord - m_center).norm() < m_radius;
                        
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

    virtual void addPatchBC() = 0;
    
    EMSolverType m_solver;
    const std::string m_bcName;
    const int m_prevFaceFlag;
    const int m_patchFlag;
    BCFunctionBase m_bcFunctionBase;
    
    Vector3D m_center { 0. , 0. , 0. };
    Real m_radius { 0. };
};

    
class PatchCircleBCEssentialNormal : public PatchBC
{
public:
    
    using PatchBC::PatchBC;
    typedef PatchBC::function_Type function_Type;
    
protected:
    
    virtual void addPatchBC()
    {
        m_solver.bcInterfacePtr() -> handler()->addBC (m_bcName, m_patchFlag,  Essential, Normal, m_bcFunctionBase);
    }
};
   

class PatchCircleBCEssentialDirectional : public PatchBC
{
public:
    
    using PatchBC::PatchBC;
    typedef PatchBC::function_Type function_Type;
    
protected:
    
    virtual void addPatchBC()
    {
        m_solver.bcInterfacePtr() -> handler()->addBC (m_bcName, m_patchFlag,  Essential, Directional, m_bcFunctionBase);
    }
};
    
    
class PatchCircleBCEssentialComponent : public PatchBC
{
public:
    
    using PatchBC::PatchBC;
    typedef PatchBC::function_Type function_Type;
    
protected:
    
    virtual void addPatchBC()
    {
        m_solver.bcInterfacePtr() -> handler()->addBC (m_bcName, m_patchFlag,  Essential, Full, m_bcFunctionBase, 3);
    }
};

    
    
}
