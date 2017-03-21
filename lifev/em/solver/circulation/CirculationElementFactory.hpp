//
//  CirculationElementFactory.hpp
//  Circulation
//
//  Created by Thomas Kummer on 13.05.15.
//  Copyright (c) 2015 Thomas Kummer. All rights reserved.
//

#ifndef CIRCULATIONELEMENTFACTORY_HPP_
#define CIRCULATIONELEMENTFACTORY_HPP_

#include <iostream>
#include <stdio.h>
#include <string>
#include <vector>
#include <memory>
#include <cmath>


class Element {
public:
    
    Element(const std::vector<std::string>& nodes,
            const std::vector<double>& param) :
        M_nodes(nodes),
        M_param(param),
        M_init(param.back())
    {}
    
    virtual ~Element() {}
    
    const std::string node(const unsigned int idx) const {return M_nodes[idx];} // todo: return pointer to Vertex object
    const std::vector<std::string> nodes() const {return M_nodes;}
    const double init() const {return M_init;}
    virtual void initRestart(const std::vector<double>& u, const std::vector<double>& uPrev0, const std::vector<double>& uPrev1, const double& time){}

    virtual const double rhs(const std::vector<double>& u, const double& time) = 0;
    virtual const double lhs(const unsigned int& idx, const std::vector<double>& u, const double& time) = 0;     // u containing (Q p1 p2 R)
    
    
protected:
    
    const std::vector<std::string> M_nodes;
    const std::vector<double> M_param;
    const double M_init;
    
    enum M_var { Q , dQ , p1 , dp1 , p2 , dp2 };

};


class Resistance : public Element {
public:
    
    using Element::Element;
    
    virtual const double rhs(const std::vector<double>& u, const double& time)
    {
        return 0.0;
    }
    
    virtual const double lhs(const unsigned int& idx, const std::vector<double>& u, const double& time)
    {
        switch (idx) {
            case Q:
                return M_param[0];
                break;
            case p1:
                return - 1.0;
                break;
            case p2:
                return 1.0;
                break;
            default:
                return 0.0;
                break;
        };
    }
};


class Inertance : public Element {
public:
    
    using Element::Element;
    
    virtual const double rhs(const std::vector<double>& u, const double& time)
    {
        return 0.0;
    }
    
    virtual const double lhs(const unsigned int& idx, const std::vector<double>& u, const double& time)
    {
        switch (idx) {
            case dQ:
                return M_param[0];
                break;
            case p1:
                return - 1.0;
                break;
            case p2:
                return 1.0;
                break;
            default:
                return 0.0;
                break;
        };
    }
};


class Capacitance : public Element {
public:
    
    using Element::Element;
    
    virtual const double rhs(const std::vector<double>& u, const double& time)
    {
        return 0.0;
    }
    
    virtual const double lhs(const unsigned int& idx, const std::vector<double>& u, const double& time)
    {
        switch (idx) {
            case Q:
                return 1.0;
                break;
            case dp1:
                return - M_param[0];
                break;
            case dp2:
                return M_param[0];
                break;
            default:
                return 0.0;
                break;
        };
    }
};

class OpenEndCapacitance : public Element {
public:
    
    using Element::Element;
    
    virtual const double rhs(const std::vector<double>& u, const double& time)
    {
        return 0.0;
    }
    
    virtual const double lhs(const unsigned int& idx, const std::vector<double>& u, const double& time)
    {
        switch (idx) {
            case Q:
                return 1.0;
                break;
            case dp1:
                return - M_param[0];
                break;
            default:
                return 0.0;
                break;
        };
    }
};

class ResistanceInertance : public Element {
public:
    
    using Element::Element;
    
    virtual const double rhs(const std::vector<double>& u, const double& time)
    {
        return 0.0;
    }
    
    virtual const double lhs(const unsigned int& idx, const std::vector<double>& u, const double& time)
    {
        switch (idx) {
            case Q:
                return M_param[0];
                break;
            case dQ:
                return M_param[1];
                break;
            case p1:
                return - 1.0;
                break;
            case p2:
                return 1.0;
                break;
            default:
                return 0.0;
                break;
        };
    }
};


class NonConstResistance : public Element {
public:
    
    using Element::Element;
    
    virtual const double rhs(const std::vector<double>& u, const double& time)
    {
        return 0.0;
    }
    
    virtual const double lhs(const unsigned int& idx, const std::vector<double>& u, const double& time)
    {
        switch (idx) {
            case Q:
                return M_param[0] * u[0];
                break;
            case p1:
                return - 1.0;
                break;
            case p2:
                return 1.0;
                break;
            default:
                return 0.0;
                break;
        };
    }
};


class Diode : public Element {
public:
    
    // using Element::Element;
    
    Diode(const std::vector<std::string>& nodes,
          const std::vector<double>& param) :
        Element ( nodes , param ),
        M_time  ( 0.0 )
    {}
    
    virtual const double rhs(const std::vector<double>& u, const double& time)
    {
        return 0.0;
    }
    
    virtual const double lhs(const unsigned int& idx, const std::vector<double>& u, const double& time)
    {
        initializeVariables(time);
        uptdateVariables(time);
        
        const double dp ( u[1] - u[2] );
        M_D = Di( dp );

        switch (idx) {
            case Q:
                return 1.0;
                break;
            case p1:
                return - M_D / M_param[0];
                break;
            case p2:
                return   M_D / M_param[0];
                break;
            default:
                return 0.0;
                break;
        };
    }

    virtual const double Di(const double& dp) const
    {
        return std::max( std::min( Dp(dp) , Du(dp) ) , Dl(dp) );
    }
    
    virtual const double Dp(const double& dp) const
    {
        return M_param[1] * std::pow(M_timestep, 2) * (dp == 0 ? 0.0 : dp / std::abs(dp)) + 2 * M_Dprev0 - M_Dprev1;
    }
    
    virtual const double Du(const double& dp) const
    {
        return ( M_Dprev0 + M_param[3] * M_timestep ) / ( 1 + M_param[3] * M_timestep );
    }
    
    virtual const double Dl(const double& dp) const
    {
        return M_Dprev0 / ( 1 + M_param[2] * M_timestep );
    }

    virtual void uptdateVariables(const double& time)
    {
        if ( time != M_time )
        {
            M_timestep = time - M_time;
            M_time = time;
            M_Dprev1 = M_Dprev0;
            M_Dprev0 = M_D;

        }
        return;
    }
    
    virtual void initializeVariables(const double& time)
    {
        if ( M_time == 0.0 )
        {
            M_time = time;
            M_timestep = time - 0.0;
            M_D = M_param[4];
            M_Dprev0 = M_param[4];
            M_Dprev1 = M_param[4];
        }
    }
    
    virtual void initRestart(const std::vector<double>& u, const std::vector<double>& uPrev0, const std::vector<double>& uPrev1, const double& time)
    {
        M_time = time;
        M_D = M_param[0] * u[0] / ( u[1] - u[2] );
        M_Dprev0 = M_param[0] * uPrev0[0] / ( uPrev0[1] - uPrev0[2] );
        M_Dprev1 = M_param[0] * uPrev1[0] / ( uPrev1[1] - uPrev1[2] );
    }
    
    
protected:
    
    double M_time;
    double M_timestep;
    
    bool M_initialized;
    double M_D;
    double M_Dprev0;
    double M_Dprev1;
    
};


class PDiode : public Diode {
public:
    
    using Diode::Diode;
    
    virtual const double lhs(const unsigned int& idx, const std::vector<double>& u, const double& time)
    {
        initializeVariables(time);
        uptdateVariables(time);

        const double dp ( u[1] - u[2] );
        M_D = Di( dp );

        switch (idx) {
            case Q:
                return ( M_param[0] * u[1] );
                break;
            case p1:
                return - M_D;
                break;
            case p2:
                return   M_D;
                break;
            default:
                return 0.0;
                break;
        };
    }
    
    virtual void initRestart(const std::vector<double>& u, const std::vector<double>& uPrev0, const std::vector<double>& uPrev1, const double& time)
    {
        M_time = time;
        M_D = M_param[0] * u[0] * u[1] / ( u[1] - u[2] );
        M_Dprev0 = M_param[0] * uPrev0[0] * uPrev0[1] / ( uPrev0[1] - uPrev0[2] );
        M_Dprev1 = M_param[0] * uPrev1[0] * uPrev0[1] / ( uPrev1[1] - uPrev1[2] );

    }

};


class ElementFactory {
public:
    
    static std::unique_ptr<Element> create(const std::string& elementType, const std::vector<std::string>& nodes, const std::vector<double>& param)
    {
        if ( elementType == "Resistance" ) return std::unique_ptr<Element> (new Resistance(nodes, param));
        if ( elementType == "Inertance" ) return std::unique_ptr<Element> (new Inertance(nodes, param));
        if ( elementType == "Compliance" ) return std::unique_ptr<Element> (new Capacitance(nodes, param));
        if ( elementType == "Resistance-Inertance" ) return std::unique_ptr<Element> (new ResistanceInertance(nodes, param));
        if ( elementType == "Q-Resistance" ) return std::unique_ptr<Element> (new NonConstResistance(nodes, param));
        if ( elementType == "OE-Compliance" ) return std::unique_ptr<Element> (new OpenEndCapacitance(nodes, param));
        if ( elementType == "Diode" ) return std::unique_ptr<Element> (new Diode(nodes, param));
        if ( elementType == "P-Diode" ) return std::unique_ptr<Element> (new PDiode(nodes, param));

        return NULL;
    }
    
};


class Vertex {
public:

    Vertex(const std::string& vertex, const std::vector<double>& param) : M_vertex(vertex), M_init(param.back()) {}
    
    virtual ~Vertex() {}
    
    const std::string name() const {return M_vertex;}
    const double init() const {return M_init;}
    
    
protected:
    
    const std::string M_vertex;
    const double M_init;
 
};


class VertexFactory {
public:
    
    static std::unique_ptr<Vertex> create(const std::string& elementType, const std::string& vertex, const std::vector<double>& param)
    {
        if ( elementType == "Vertex" ) return std::unique_ptr<Vertex> (new Vertex(vertex, param));
        return NULL;
    }
    
};


#endif
