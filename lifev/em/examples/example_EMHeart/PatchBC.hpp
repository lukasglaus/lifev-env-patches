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
#include <lifev/em/examples/example_EMHeart/GenericFactory.hpp>


namespace LifeV
{

    
class PatchBCFunctionBaseCreator
{
public:
    
    typedef boost::function<Real ( const Real&, const Real&, const Real&, const Real&, const ID& ) > function_Type;

    PatchBCFunctionBaseCreator(){}
    
    function_Type fct()
    {
        function_Type f;
        f = boost::bind (&PatchBCFunctionBaseCreator::patchForceFunction, this, _1, _2, _3, _4, _5);
        return f;
    }
    
    Real patchForceFunction (const Real& t, const Real&  X, const Real& Y, const Real& Z, const ID& i)
    {
        return (t * 1e-6);
//        switch (i)
//        {
//            case 0:
//                return (m_direction[0] * t * 1e-8);
//                break;
//            case 1:
//                return (m_direction[1] * t * 1e-8);
//                break;
//            case 2:
//                return (m_direction[2] * t * 1e-8);
//                break;
//            default:
//                ERROR_MSG ("This entry is not allowed");
//                return 0.;
//                break;
//        }
    }
    
    function_Type dirFct()
    {
        function_Type f;
        f = boost::bind (&PatchBCFunctionBaseCreator::directionFct, this, _1, _2, _3, _4, _5);
        return f;
    }

    
    Real directionFct (const Real& t, const Real&  X, const Real& Y, const Real& Z, const ID& i)
    {
        switch (i)
        {
            case 0:
                return 1;//(m_direction[0]);
                break;
            case 1:
                return (m_direction[1]);
                break;
            case 2:
                return (m_direction[2]);
                break;
            default:
                ERROR_MSG ("This entry is not allowed");
                return 0.;
                break;
        }
    }

    
    Real sinSquared (const Real& t, const Real& Tmax, const Real& tmax, const Real& tduration)
    {
        bool time ( fmod(t-tmax+0.5*tduration, 800.) < tduration && fmod(t-tmax+0.5*tduration, 800.) > 0);
        Real force = std::pow( std::sin(fmod(t-tmax+0.5*tduration, 800.)*3.14159265359/tduration) , 2 ) * Tmax;
        return ( time ? force : 0 );
    }

    void setDirection(const Vector3D& direction)
    {
        m_direction = direction;
    }
    
    
private:
    
    Vector3D m_direction {0.,0.,0.};
    
};

    
class PatchBC
{
public:
    
    typedef EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > > EMSolverType;
    typedef boost::function<Real ( const Real&, const Real&, const Real&, const Real&, const ID& ) > function_Type;
    
//    PatchBC(){}
    
    PatchBC(EMSolverType& solver, const std::string& bcName, const int& prevFaceFlag, const int& patchFlag) :
        m_solver (solver),
        m_bcName (bcName),
        m_prevFaceFlag (prevFaceFlag),
        m_patchFlag (patchFlag)
    {}

//    void initialize(boost::shared_ptr<EMSolverType> solver, const std::string& bcName, const int& prevFaceFlag, const int& patchFlag)
//    {
//        m_solver = solver;
//        m_bcName = bcName;
//        m_prevFaceFlag = prevFaceFlag;
//        m_patchFlag = patchFlag;
//    }
    
    void setup(Vector3D& direction, const Vector3D& center, const Real& radius)
    {
        setParameters(center, radius, direction);
        setBCFunctionBase();
        setBCFunctionDirectional();
        createPatchArea();
        addPatchBC();
    }
    
protected:
    
    void setParameters(const Vector3D& center, const Real& radius, Vector3D& direction)
    {
        m_center = center;
        m_radius = radius;
        m_patchBCFunctionBaseCreator.setDirection(direction.normalized());
    }
    
    void setBCFunctionBase()
    {
        BCFunctionBase bcFB (m_patchBCFunctionBaseCreator.fct());
        m_bcFunctionBase.setFunction(m_patchBCFunctionBaseCreator.fct());
    }
    
    void setBCFunctionDirectional()
    {
        m_bcFunctionDirectional.setFunctions_Directional(m_patchBCFunctionBaseCreator.fct(), m_patchBCFunctionBaseCreator.dirFct());
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
    BCFunctionDirectional m_bcFunctionDirectional;
    
    PatchBCFunctionBaseCreator m_patchBCFunctionBaseCreator;
    
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
        m_solver.bcInterfacePtr()->handler()->addBC (m_bcName, m_patchFlag, Essential, Normal, m_bcFunctionBase);
    }
};
   
//REGISTER(PatchBC, PatchCircleBCEssentialNormal);
    
    
class PatchCircleBCEssentialDirectional : public PatchBC
{
public:
    
    using PatchBC::PatchBC;
    typedef PatchBC::function_Type function_Type;
    
protected:
    
    virtual void addPatchBC()
    {
        BCFunctionDirectional bcFunctionDirectional (m_bcFunctionBase, m_bcFunctionDirectional);
        m_solver.bcInterfacePtr()->handler()->addBC (m_bcName, m_patchFlag, Essential, Directional, bcFunctionDirectional);
    }
};
    
//REGISTER(PatchBC, PatchCircleBCEssentialDirectional);

    
class PatchCircleBCEssentialFull : public PatchBC
{
public:
    
    using PatchBC::PatchBC;
    typedef PatchBC::function_Type function_Type;
    
protected:
    
    virtual void addPatchBC()
    {
        m_solver.bcInterfacePtr()->handler()->addBC (m_bcName, m_patchFlag, Essential, Full, m_bcFunctionBase, 3);
    }
};

//REGISTER(PatchBC, PatchCircleBCEssentialFull);


class PatchCircleBCEssentialComponent : public PatchBC
{
public:
    
    using PatchBC::PatchBC;
    typedef PatchBC::function_Type function_Type;
    
protected:
    
    virtual void addPatchBC()
    {
        m_solver.bcInterfacePtr()->handler()->addBC (m_bcName, m_patchFlag, Essential, Component, m_bcFunctionBase, 0);
    }
};

//REGISTER(PatchBC, PatchCircleBCEssentialComponent);

    
class PatchCircleBCNaturalComponent : public PatchBC
{
public:
    
    using PatchBC::PatchBC;
    typedef PatchBC::function_Type function_Type;
    
protected:
    
    virtual void addPatchBC()
    {
        m_solver.bcInterfacePtr()->handler()->addBC (m_bcName, m_patchFlag, Natural, Component, m_bcFunctionBase, 0);
    }
};

//REGISTER(PatchBC, PatchCircleBCEssentialComponent);
    
}

#endif /* PatchBC_hpp */
