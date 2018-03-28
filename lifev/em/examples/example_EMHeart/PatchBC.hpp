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

Real patchForceFunction (const Real& t, const Real&  X, const Real& Y, const Real& Z, const ID& i)
{
    return (t * 1e-5);
}

    
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
    
    static Real patchDispFun (const Real& t, const Real&  X, const Real& Y, const Real& Z, const ID& i)
    {
        switch (i)
        {
            case 0:
                return (0.00005*t);
                break;
            case 1:
                return 0;
                break;
            case 2:
                return (0.00005*t);
                break;
            default:
                ERROR_MSG ("This entry is not allowed");
                return 0.;
                break;
        }
    }

    Real patchForceFunction (const Real& t, const Real&  X, const Real& Y, const Real& Z, const ID& i)
    {
        return (t * 1e-5);
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
                return 1;//(m_direction[1]);
                break;
            case 2:
                return 1;//(m_direction[2]);
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

    
    class PatchBCEssentialCircular : public PatchBC
    {
    public:
        
        
        std::vector<Real> patchDisplacement;
        std::vector<Vector3D> patchDirection;
        
        std::vector<vectorPtr_Type> patchDispVecPtr;
        std::vector<bcVectorPtr_Type> patchDispBCVecPtr;
        
        VectorEpetra dispPreload (disp);
        
        for ( UInt i (0) ; i < nDispPatchBC ; ++i )
        {
            std::string patchName = dataFile ( ( "solid/boundary_conditions/listEssentialPatchBC" ), " ", i );
            Real patchFlag = dataFile ( ("solid/boundary_conditions/" + patchName + "/flag").c_str(), 0 );
            Real patchRadius = dataFile ( ("solid/boundary_conditions/" + patchName + "/radius").c_str(), 1.0 );
            patchDisplacement.push_back( dataFile ( ("solid/boundary_conditions/" + patchName + "/displacement").c_str(), 1.0 ) );
            
            Vector3D patchCenter;
            for ( UInt j (0); j < 3; ++j )
            {
                patchCenter[j] = dataFile ( ("solid/boundary_conditions/" + patchName + "/center").c_str(), 0, j );
            }
            
            Vector3D pd;
            for ( UInt j (0); j < 3; ++j )
            {
                pd[j] = dataFile ( ("solid/boundary_conditions/" + patchName + "/direction").c_str(), 0, j );
            }
            patchDirection.push_back(pd);
            
            UInt componentSize = dataFile.vector_variable_size ( ("solid/boundary_conditions/" + patchName + "/component").c_str() );
            std::vector<ID> patchComponent (componentSize);
            for ( UInt j (0); j < componentSize; ++j )
            {
                patchComponent[j] = dataFile ( ("solid/boundary_conditions/" + patchName + "/component").c_str(), 0, j );
            }
            
            patchDispVecPtr.push_back ( heartSolver.directionalVectorField(FESpace, patchDirection[i], 1e-10) );
            //*patchDispVecPtr[i] += dispPreload;
            patchDispBCVecPtr.push_back ( bcVectorPtr_Type( new bcVector_Type( *patchDispVecPtr[i], solver.structuralOperatorPtr() -> dispFESpacePtr() -> dof().numTotalDof(), 1 ) ) );
            solver.bcInterfacePtr() -> handler()->addBC (patchName, (900+i),  Essential, Component, *patchDispBCVecPtr[i], patchComponent);
        }
        
        Real patchDispOffset = dataFile ( "solid/patches/patch_disp_offset", 0. );
        Real tmax = dataFile ( "solid/patches/tmax", 0. );
        Real tduration = dataFile ( "solid/patches/tduration", 0. );
        
        auto modifyEssentialPatchBC = [&] (const Real& time)
        {
            for ( UInt i (0) ; i < nDispPatchBC ; ++i )
            {
                Real currentPatchDisp = heartSolver.sinSquared(time, patchDisplacement[i], tmax, tduration);
                currentPatchDisp += patchDispOffset;
                
                patchDispVecPtr[i] = heartSolver.directionalVectorField(FESpace, patchDirection[i], currentPatchDisp);
                //*patchDispVecPtr[i] += dispPreload;
                if ( 0 == comm->MyPID() ) std::cout << "\nCurrent patch-" << i << " displacement: " << currentPatchDisp << " cm";
                    
                    patchDispBCVecPtr[i].reset( new bcVector_Type( *patchDispVecPtr[i], FESpace->dof().numTotalDof(), 1 ) );
                    solver.bcInterfacePtr()->handler()->modifyBC((900+i), *patchDispBCVecPtr[i]);
                    }
                    
                    if ( 0 == comm->MyPID() ) std::cout << std::endl;
                    };
        
        
    protected:
    }
    
  
                    
                    
                    
                    
                    
                    
                    
                    

//                    Real patchForceFunction (const Real& t, const Real&  X, const Real& Y, const Real& Z, const ID& i)
//                {
//                    return (t * 1e-5);
//                }
//
//
//                    class PatchBCFunctionBaseCreator
//                {
//                public:
//
//                    typedef boost::function<Real ( const Real&, const Real&, const Real&, const Real&, const ID& ) > function_Type;
//
//                    PatchBCFunctionBaseCreator(){}
//
//                    function_Type fct()
//                    {
//                        function_Type f;
//                        f = boost::bind (&PatchBCFunctionBaseCreator::patchForceFunction, this, _1, _2, _3, _4, _5);
//                        return f;
//                    }
//
//                    static Real patchDispFun (const Real& t, const Real&  X, const Real& Y, const Real& Z, const ID& i)
//                    {
//                        switch (i)
//                        {
//                            case 0:
//                                return (0.00005*t);
//                                break;
//                            case 1:
//                                return 0;
//                                break;
//                            case 2:
//                                return (0.00005*t);
//                                break;
//                            default:
//                                ERROR_MSG ("This entry is not allowed");
//                                return 0.;
//                                break;
//                        }
//                    }
//
//                    Real patchForceFunction (const Real& t, const Real&  X, const Real& Y, const Real& Z, const ID& i)
//                    {
//                        return (t * 1e-5);
//                        //        switch (i)
//                        //        {
//                        //            case 0:
//                        //                return (m_direction[0] * t * 1e-8);
//                        //                break;
//                        //            case 1:
//                        //                return (m_direction[1] * t * 1e-8);
//                        //                break;
//                        //            case 2:
//                        //                return (m_direction[2] * t * 1e-8);
//                        //                break;
//                        //            default:
//                        //                ERROR_MSG ("This entry is not allowed");
//                        //                return 0.;
//                        //                break;
//                        //        }
//                    }
//
//                    function_Type dirFct()
//                    {
//                        function_Type f;
//                        f = boost::bind (&PatchBCFunctionBaseCreator::directionFct, this, _1, _2, _3, _4, _5);
//                        return f;
//                    }
//
//
//                    Real directionFct (const Real& t, const Real&  X, const Real& Y, const Real& Z, const ID& i)
//                    {
//                        switch (i)
//                        {
//                            case 0:
//                                return 1;//(m_direction[0]);
//                                break;
//                            case 1:
//                                return 1;//(m_direction[1]);
//                                break;
//                            case 2:
//                                return 1;//(m_direction[2]);
//                                break;
//                            default:
//                                ERROR_MSG ("This entry is not allowed");
//                                return 0.;
//                                break;
//                        }
//                    }
//
//
//                    Real sinSquared (const Real& t, const Real& Tmax, const Real& tmax, const Real& tduration)
//                    {
//                        bool time ( fmod(t-tmax+0.5*tduration, 800.) < tduration && fmod(t-tmax+0.5*tduration, 800.) > 0);
//                        Real force = std::pow( std::sin(fmod(t-tmax+0.5*tduration, 800.)*3.14159265359/tduration) , 2 ) * Tmax;
//                        return ( time ? force : 0 );
//                    }
//
//                    void setDirection(const Vector3D& direction)
//                    {
//                        m_direction = direction;
//                    }
//
//
//                private:
//
//                    Vector3D m_direction {0.,0.,0.};
//
//                };
//
//
//                    class PatchBC
//                {
//                public:
//
//                    typedef EMSolver<RegionMesh<LinearTetra>, EMMonodomainSolver<RegionMesh<LinearTetra> > > EMSolverType;
//                    typedef boost::function<Real ( const Real&, const Real&, const Real&, const Real&, const ID& ) > function_Type;
//
//                    //    PatchBC(){}
//
//                    PatchBC(EMSolverType& solver, const std::string& bcName, const int& prevFaceFlag, const int& patchFlag) :
//                    m_solver (solver),
//                    m_bcName (bcName),
//                    m_prevFaceFlag (prevFaceFlag),
//                    m_patchFlag (patchFlag)
//                    {}
//
//                    //    void initialize(boost::shared_ptr<EMSolverType> solver, const std::string& bcName, const int& prevFaceFlag, const int& patchFlag)
//                    //    {
//                    //        m_solver = solver;
//                    //        m_bcName = bcName;
//                    //        m_prevFaceFlag = prevFaceFlag;
//                    //        m_patchFlag = patchFlag;
//                    //    }
//
//                    void setup(Vector3D& direction, const Vector3D& center, const Real& radius)
//                    {
//                        setParameters(center, radius, direction);
//                        setBCFunctionBase();
//                        setBCFunctionDirectional();
//                        createPatchArea();
//                        addPatchBC();
//                    }
//
//                protected:
//
//                    void setParameters(const Vector3D& center, const Real& radius, Vector3D& direction)
//                    {
//                        m_center = center;
//                        m_radius = radius;
//                        m_patchBCFunctionBaseCreator.setDirection(direction.normalized());
//                    }
//
//                    void setBCFunctionBase()
//                    {
//                        BCFunctionBase bcFB (m_patchBCFunctionBaseCreator.fct());
//                        m_bcFunctionBase.setFunction(m_patchBCFunctionBaseCreator.fct());
//                    }
//
//                    void setBCFunctionDirectional()
//                    {
//                        m_bcFunctionDirectional.setFunctions_Directional(m_patchBCFunctionBaseCreator.fct(), m_patchBCFunctionBaseCreator.dirFct());
//                    }
//
//                    virtual void createPatchArea()
//                    {
//                        for (auto& mesh : m_solver.mesh())
//                        {
//                            for (int j(0); j < mesh->numBoundaryFacets(); j++)
//                            {
//                                auto& face = mesh->boundaryFacet(j);
//                                auto faceFlag = face.markerID();
//
//                                if (faceFlag == m_prevFaceFlag || faceFlag == 470 || faceFlag == 471)
//                                {
//                                    int numPointsInsidePatch (0);
//
//                                    for (int k(0); k < 3; ++k)
//                                    {
//                                        auto coord = face.point(k).coordinates();
//                                        bool pointInPatch = (coord - m_center).norm() < m_radius;
//
//                                        if (pointInPatch)
//                                        {
//                                            ++numPointsInsidePatch;
//                                        }
//                                    }
//
//                                    if (numPointsInsidePatch > 2)
//                                    {
//                                        face.setMarkerID(m_patchFlag);
//                                    }
//                                }
//                            }
//                        }
//                    }
//
//                    virtual void addPatchBC() = 0;
//
//                    EMSolverType m_solver;
//                    const std::string m_bcName;
//                    const int m_prevFaceFlag;
//                    const int m_patchFlag;
//
//                    BCFunctionBase m_bcFunctionBase;
//                    BCFunctionDirectional m_bcFunctionDirectional;
//
//                    PatchBCFunctionBaseCreator m_patchBCFunctionBaseCreator;
//
//                    Vector3D m_center { 0. , 0. , 0. };
//                    Real m_radius { 0. };
//                };

                    
                    
//class PatchCircleBCEssentialNormal : public PatchBC
//{
//public:
//
//    using PatchBC::PatchBC;
//    typedef PatchBC::function_Type function_Type;
//
//protected:
//
//    virtual void addPatchBC()
//    {
//        m_solver.bcInterfacePtr()->handler()->addBC (m_bcName, m_patchFlag, Essential, Normal, m_bcFunctionBase);
//    }
//};
//
////REGISTER(PatchBC, PatchCircleBCEssentialNormal);
//
//
//class PatchCircleBCEssentialDirectional : public PatchBC
//{
//public:
//
//    using PatchBC::PatchBC;
//    typedef PatchBC::function_Type function_Type;
//
//protected:
//
//    virtual void addPatchBC()
//    {
//        BCFunctionDirectional bcFunctionDirectional (m_bcFunctionBase, m_bcFunctionDirectional);
//        m_solver.bcInterfacePtr()->handler()->addBC (m_bcName, m_patchFlag, Essential, Directional, bcFunctionDirectional);
//    }
//};
//
////REGISTER(PatchBC, PatchCircleBCEssentialDirectional);
//
//
//class PatchCircleBCEssentialFull : public PatchBC
//{
//public:
//
//    using PatchBC::PatchBC;
//    typedef PatchBC::function_Type function_Type;
//
//protected:
//
//    virtual void addPatchBC()
//    {
//        m_solver.bcInterfacePtr()->handler()->addBC (m_bcName, m_patchFlag, Essential, Full, m_bcFunctionBase, 3);
//    }
//};
//
////REGISTER(PatchBC, PatchCircleBCEssentialFull);
//
//
//class PatchCircleBCEssentialComponent : public PatchBC
//{
//public:
//
//    using PatchBC::PatchBC;
//    typedef PatchBC::function_Type function_Type;
//
//protected:
//
//    virtual void addPatchBC()
//    {
//        BCFunctionBase bcFB (patchForceFunction);
//
//        m_solver.bcInterfacePtr()->handler()->addBC (m_bcName, m_patchFlag, Essential, Component, bcFB, 0);
//    }
//};
//
////REGISTER(PatchBC, PatchCircleBCEssentialComponent);
//
//
//class PatchCircleBCNaturalComponent : public PatchBC
//{
//public:
//
//    using PatchBC::PatchBC;
//    typedef PatchBC::function_Type function_Type;
//
//protected:
//
//    virtual void addPatchBC()
//    {
//        m_solver.bcInterfacePtr()->handler()->addBC (m_bcName, m_patchFlag, Natural, Component, m_bcFunctionBase, 0);
//    }
//};
//
////REGISTER(PatchBC, PatchCircleBCEssentialComponent);
    
}

#endif /* PatchBC_hpp */
