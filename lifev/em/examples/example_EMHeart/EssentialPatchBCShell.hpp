//
//  EssentialPatchBCCircularSmooth.hpp
//  lifev-heart
//
//  Created by Thomas Kummer on 30.06.18.
//  Copyright © 2018 Thomas Kummer. All rights reserved.
//

#ifndef EssentialPatchBCShell_hpp
#define EssentialPatchBCShell_hpp

#include <stdio.h>
#include <lifev/em/examples/example_EMHeart/EssentialPatchBC.hpp>

#define PI 3.14159265359

namespace LifeV
{
    
class EssentialPatchBCShell : public EssentialPatchBC
{
public:
    
    EssentialPatchBCShell(){}
    ~EssentialPatchBCShell(){}
    
    virtual void setup(const GetPot& dataFile, const std::string& name)
    {        
        super::setup(dataFile, name);
	m_Name = name;
        
	m_PrevFlag = dataFile ( ("solid/boundary_conditions/" + m_Name + "/flag").c_str(), 0 );
        m_patchDisplacement = dataFile ( ("solid/boundary_conditions/" + m_Name + "/displacement").c_str(), 1.0 );
        m_height= dataFile ( ("solid/boundary_conditions/" + m_Name + "/height").c_str(), 1.0 );
	m_lowerlength= dataFile ( ("solid/boundary_conditions/" + m_Name + "/lowerlength").c_str(), 1.0 );
	m_upperlength= dataFile ( ("solid/boundary_conditions/" + m_Name + "/upperlength").c_str(), 1.0 );
	m_width= dataFile ( ("solid/boundary_conditions/" + m_Name + "/width").c_str(), 1.0 );      

        
        for ( UInt j (0); j < 3; ++j )
        {
            m_origin[j] = dataFile ( ("solid/boundary_conditions/" + m_Name + "/origin").c_str(), 0, j );
            m_z1[j] = dataFile ( ("solid/boundary_conditions/" + m_Name + "/direction").c_str(), 0, j );
        }
        
//	m_x0=[1 0 0];
//	m_y0=[0 1 0];
//	m_z0=[0 0 1];
	m_x0[0]=1;
	m_x0[1]=0;
	m_x0[2]=0;
	m_y0[0]=0;
	m_y0[1]=1;
	m_y0[2]=0;
	m_z0[0]=0;
	m_z0[1]=0;
	m_z0[2]=1;
	
	alpha=acos(m_z0.dot(m_z1)/(m_z0.norm()*m_z1.norm()));
	
	R[0][0]=cos(alpha);
	R[1][0]=0;
	R[2][0]=sin(alpha);
	R[0][1]=0;
	R[1][1]=1;
	R[2][1]=0;
	R[0][3]=-sin(alpha);
	R[1][3]=0;
	R[3][3]=cos(alpha);

	m_x1=(R*m_x0).normalized();
	m_y1=m_y0;
	m_z1=m_z1.normalized();

       	m_tmax = dataFile ( "solid/patches/tmax", 0. );
        m_tduration = dataFile ( "solid/patches/tduration", 0. );
    }
    
    
protected:
    
       
    virtual const bool nodeOnPatch(const Vector3D& coord) const
    {
	bool pointInShell;
	Vector3D Node;
	double ly;

	ly=m_lowerlength+(m_upperlength-m_lowerlength)/m_height*(coord[1]-(m_origin[1]-m_height/2));
	Node=coord-m_origin;	
	
	pointInShell=(-m_width/2<=Node.dot(m_z1) && m_width/2>=Node.dot(m_z1)
			&& -m_height/2<=Node.dot(m_y1) && m_height/2>=Node.dot(m_y1)
			&& -ly/2<=Node.dot(m_x1) && ly/2>=Node.dot(m_x1));
        return pointInShell;
    }
 
    
    
    Real m_patchDisplacement;

    Vector3D m_origin;
    Real m_height;
    Real m_lowerlength;
    Real m_upperlength;
    Real m_width;
    Vector3D m_x0; 
    Vector3D m_y0;
    Vector3D m_z0;
    Vector3D m_x1;
    Vector3D m_y1;
    Vector3D m_z1;
    MatrixSmall<3,3> R;
    double alpha;
    
    Real m_tmax;
    Real m_tduration;
       
    Vector3D m_Center;
    Real m_Radius;
    Real m_EdgeDispFactor;

};

REGISTER(EssentialPatchBC, EssentialPatchBCShell);

}

#endif /* EssentialPatchBCShell_hpp */
