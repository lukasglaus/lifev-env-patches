//
//  CirculationTimeIntegrator.hpp
//  Circulation
//
//  Created by Thomas Kummer on 16.05.15.
//  Copyright (c) 2015 Thomas Kummer. All rights reserved.
//

#include <stdio.h>
#include <vector>
#include <string>


namespace LifeV
{

class VolumeIntegrator {
public:
    
    VolumeIntegrator(const std::vector<int>& bdFlags,
                     const std::string& domain,
                     const boost::shared_ptr<RegionMesh<LinearTetra> > mesh,
                     const boost::shared_ptr <ETFESpace<RegionMesh<LinearTetra>, MapEpetra, 3, 1> > ETFESpace) :
                M_meshPtr       ( mesh ),
                M_mesh          ( *mesh ),
                M_bdFlags       ( bdFlags ),
                M_domain        ( domain ),
                M_ETFESpace     ( ETFESpace )
    {
        initialize();
    }
    
    virtual ~VolumeIntegrator() {}
    
    void initialize()
    {
        findBoundaryPoints();
        sortBoundaryPoints();
    }
    
    void findBoundaryPoints()
    {
        std::set<unsigned int> vertexIds;
        for (UInt iBFaceIn = 0; iBFaceIn < M_mesh.numBFaces(); ++iBFaceIn)
        {
            UInt markerIdIn = M_mesh.boundaryFace(iBFaceIn).markerID();
            if ( std::find(M_bdFlags.begin(), M_bdFlags.end(), markerIdIn) != M_bdFlags.end() )
            {
                for (UInt iBFaceOut = 0; iBFaceOut < M_mesh.numBFaces(); ++iBFaceOut)
                {
                    UInt markerIdOut = M_mesh.boundaryFace(iBFaceOut).markerID();
                    if ( std::find(M_bdFlags.begin(), M_bdFlags.end(), markerIdOut) == M_bdFlags.end() )
                    {
                        for (UInt iBPointIn = 0; iBPointIn < M_mesh.boundaryFace(iBFaceIn).S_numPoints; ++iBPointIn)
                        {
                            for (UInt iBPointOut = 0; iBPointOut < M_mesh.boundaryFace(iBFaceOut).S_numPoints; ++iBPointOut)
                            {
                                UInt pointIdIn = M_mesh.boundaryFace(iBFaceIn).point(iBPointIn).id();
                                UInt pointIdOut = M_mesh.boundaryFace(iBFaceOut).point(iBPointOut).id();
                                
                                if ( pointIdIn == pointIdOut )
                                {
                                    vertexIds.insert(pointIdIn);
                                }
                            }
                        }
                    }
                }
            }
        }
        M_boundaryPoints.clear();
        for (auto it = vertexIds.begin(); it != vertexIds.end(); ++it) M_boundaryPoints.push_back(*it);
    }
    
    
    void sortBoundaryPoints()
    {
        // Determine center and normal vector of point cloud in reference position
        Vector3D center0 = center();
        Vector3D normal0 = normal(center0);
        
        
        // Order points by angles
        std::map<double, unsigned int> vertexIdsOrdered;
        Vector3D axis0  = ( M_mesh.point(*M_boundaryPoints.begin()).coordinates() - center0 ).normalized();
        Vector3D axis1  = ( normal0.cross(axis0) ).normalized();
        
        for (auto it = M_boundaryPoints.begin(); it != M_boundaryPoints.end(); ++it)
        {
            double x0 ( ( M_mesh.point(*it).coordinates() - center0 ).dot(axis0) );
            double x1 ( ( M_mesh.point(*it).coordinates() - center0 ).dot(axis1) );
            
            double angle = std::atan2( x1 , x0 );
            vertexIdsOrdered.insert( std::pair<double, unsigned int> (angle, *it) );
        }

        
        // Store sorted points in input vector
        M_boundaryPoints.clear();
        for (auto it = vertexIdsOrdered.begin(); it != vertexIdsOrdered.end(); ++it) M_boundaryPoints.push_back(it->second);
        
    }

    
    const Real computeOpenEndVolume (const VectorEpetra& disp,
                                     const int direction = - 1,
                                     const unsigned int component = 0) const
    {
        Vector3D componentVector; componentVector (component) = 1;
        
        auto positionVector ( currentPositionVector(disp) );
        Vector3D centerPoint ( center(positionVector) );

        // Compute volume
        Int nLocalDof = positionVector.epetraVector().MyLength();
        Int nComponentLocalDof = nLocalDof / 3;
        auto itNext = M_boundaryPoints.begin();
        unsigned int i (0);
        Real volume (0.0);
        for (auto it = M_boundaryPoints.begin(); it != M_boundaryPoints.end(); ++it)
        {
            if ( i++ < M_boundaryPoints.size() - 1 ) std::advance(itNext, 1);
            else std::advance(itNext, - (M_boundaryPoints.size() - 1));
            
            Vector3D P1 ( positionVector[ *it ] , positionVector[ *it + nComponentLocalDof ] , positionVector[ *it + 2 * nComponentLocalDof ] );
            Vector3D P2 ( positionVector[ *itNext ] , positionVector[ *itNext + nComponentLocalDof ] , positionVector[ *itNext + 2 * nComponentLocalDof ] );
            
            Vector3D v1 = P1 - centerPoint;
            Vector3D v2 = P2 - centerPoint;
            
            Vector3D centerTriangle = ( P1 + P2 + centerPoint ) / 3;
            Vector3D normal = ( v1.cross(v2) ).normalized();
            Real area = ( v1.cross(v2) ).norm() / 2;
            
            Real areaProjected = area * normal.dot(componentVector);
            
            volume += areaProjected * centerTriangle.dot(componentVector);
        }
        
        return direction * volume;
    }
    
    
    template<class space>
    Real computeBoundaryVolume (const VectorEpetra& disp,
                                const boost::shared_ptr <space> dETFESpace,
                                int bdFlag) const
    {
        Real fluidVolume;
        
        MatrixSmall<3, 3> Id;
        Id (0, 0) = 1.; Id (0, 1) = 0.; Id (0, 2) = 0.;
        Id (1, 0) = 0.; Id (1, 1) = 1.; Id (1, 2) = 0.;
        Id (2, 0) = 0.; Id (2, 1) = 0.; Id (2, 2) = 1.;
        VectorSmall<3> E1;
        E1 (0) = 1.; E1 (1) = 0.; E1 (2) = 0.;
        
        const VectorEpetra positionVector ( currentPositionVector(disp) );
        boost::shared_ptr<VectorEpetra> intergral ( new VectorEpetra ( positionVector.map() ) );
        
        {
            using namespace ExpressionAssembly;
            
            BOOST_AUTO_TPL (I, value (Id) );
            BOOST_AUTO_TPL (vE1, value (E1) );
            BOOST_AUTO_TPL (Grad_u, grad (dETFESpace, disp, 0) );
            BOOST_AUTO_TPL (F, (Grad_u + I) );
            BOOST_AUTO_TPL (FmT, minusT (F) );
            BOOST_AUTO_TPL (J, det (F) );
            
            QuadratureBoundary myBDQR (buildTetraBDQR (quadRuleTria4pt) );
            
            *intergral *= 0.0;
            integrate (boundary (M_meshPtr, bdFlag), myBDQR, M_ETFESpace,
                       value(-1.0) * J * dot (vE1, FmT * Nface) * phi_i) >> intergral;
            
            intergral->globalAssemble();
            
            fluidVolume = positionVector.dot (*intergral);
  
            return fluidVolume;
        }
    }
    
    
    template<class space>
    const Real volume(const VectorEpetra& disp,
                      const boost::shared_ptr <space> dETFESpace,
                      boost::shared_ptr<Epetra_Comm> comm,
                      const int direction = - 1,
                      const unsigned int component = 0)
    {
        // Compute volume over boundary
        Real volumeBoundary (0);
        for ( auto& bdFlag : M_bdFlags )
        {
            std::cout << bdFlag << std::endl;
            volumeBoundary += computeBoundaryVolume(disp, dETFESpace, bdFlag);
        }
        
        // Compute volume over open-end-boundary
        Real volumeOpenEnd = computeOpenEndVolume(disp);
        
        // Compute total volume
        Real totalVolume = volumeBoundary + volumeOpenEnd;
        
        if (comm->MyPID() == 0)
        {
            std::cout << "Volume in " << M_domain << ": " << totalVolume << std::endl;
        }
        
        return totalVolume;
    }
    
    
protected:
    
    const VectorEpetra currentPositionVector (const VectorEpetra& disp) const
    {
        VectorEpetra positionVector ( disp.map() );
        Int nLocalDof = disp.epetraVector().MyLength();
        Int nComponentLocalDof = nLocalDof / 3;
        for (int k (0); k < nComponentLocalDof; k++)
        {
            UInt iGID = positionVector.blockMap().GID (k);
            UInt jGID = positionVector.blockMap().GID (k + nComponentLocalDof);
            UInt kGID = positionVector.blockMap().GID (k + 2 * nComponentLocalDof);
            
            positionVector[iGID] = M_mesh.point (iGID).x() + disp[iGID];
            positionVector[jGID] = M_mesh.point (iGID).y() + disp[jGID];
            positionVector[kGID] = M_mesh.point (iGID).z() + disp[kGID];
        }
        return positionVector;
    }
    
    
    const Vector3D center(const VectorEpetra& positionVector) const
    {
        Int nLocalDof = positionVector.epetraVector().MyLength();
        Int nComponentLocalDof = nLocalDof / 3;
        Vector3D center;
        for (auto it = M_boundaryPoints.begin(); it != M_boundaryPoints.end(); ++it)
        {
            center (0) += positionVector[*it] / M_boundaryPoints.size();
            center (1) += positionVector[*it + nComponentLocalDof] / M_boundaryPoints.size();
            center (2) += positionVector[*it + 2 * nComponentLocalDof] / M_boundaryPoints.size();
        }
        return center;
    }
    
    
    const Vector3D center() const
    {
        Vector3D center0;
        for (auto it = M_boundaryPoints.begin(); it != M_boundaryPoints.end(); ++it)
        {
            center0 += M_mesh.point(*it).coordinates() / M_boundaryPoints.size();
        }
        return center0;
    }
    
    
    const Vector3D normal(const Vector3D& center0) const
    {
        auto itNext0 = M_boundaryPoints.begin();
        unsigned int j (0);
        Vector3D normal0;
        for (auto it = M_boundaryPoints.begin(); it != M_boundaryPoints.end(); ++it)
        {
            if ( j++ < M_boundaryPoints.size() - 1 ) std::advance(itNext0, 1);
            else std::advance(itNext0, - (M_boundaryPoints.size() - 1));
            
            Vector3D P1 ( M_mesh.point(*it).coordinates() );
            Vector3D P2 ( M_mesh.point(*itNext0).coordinates() );
            
            Vector3D v1 = P1 - center0;
            Vector3D v2 = P2 - center0;
            
            normal0 += ( v1.cross(v2) ).normalized();
        }
        return normal0.normalized();
    }
    
    const boost::shared_ptr<RegionMesh<LinearTetra> > M_meshPtr;
    const RegionMesh<LinearTetra>& M_mesh;
    const boost::shared_ptr <ETFESpace<RegionMesh<LinearTetra>, MapEpetra, 3, 1> > M_ETFESpace;
    
    const std::vector<int> M_bdFlags;
    const std::string M_domain;
    std::vector<unsigned int> M_boundaryPoints;
    
};

}