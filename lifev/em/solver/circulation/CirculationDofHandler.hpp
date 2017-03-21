//
//  CirculationDofHandler.cpp
//  Circulation
//
//  Created by Thomas Kummer on 15.05.15.
//  Copyright (c) 2015 Thomas Kummer. All rights reserved.
//

#ifndef CIRCULATIONDOFHANDLER_HPP_
#define CIRCULATIONDOFHANDLER_HPP_

#include <stdio.h>
#include "CirculationElementFactory.hpp"
#include "CirculationGridView.hpp"


class DofHandler {
public:
    
    DofHandler(GridView& gv) :
        M_vertices ( gv.vertices() ),
        M_elements ( gv.elements() )
    {}
    
    virtual ~DofHandler() {}
    
    const unsigned int operator()(const std::shared_ptr<Element>& element) const
    {
        for ( unsigned int i (0) ; i < sizeElements() ; ++i )
        {
            if ( M_elements[i]->node(0) == element->node(0) && M_elements[i]->node(1) == element->node(1) ) return ( sizeVertices() + i );
        }
        return size();
    }
    
    
    const unsigned int operator()(const std::vector<std::string>& nodes) const
    {
        for ( unsigned int i (0) ; i < sizeElements() ; ++i )
        {
            if ( M_elements[i]->node(0) == nodes[0] && M_elements[i]->node(1) == nodes[1] ) return ( sizeVertices() + i );
        }
        return size();
    }
    
    
    const unsigned int operator()(const std::shared_ptr<Vertex>& vertex) const
    {
        for ( unsigned int i (0) ; i < sizeVertices() ; ++i )
        {
            if ( M_vertices[i]->name() == vertex->name() ) return i;
        }
        return size();
    }
    
    
    const unsigned int operator()(const std::string& vertex) const
    {
        for ( unsigned int i (0) ; i < sizeVertices() ; ++i )
        {
            if ( M_vertices[i]->name() == vertex ) return i;
        }
        return size();
    }
    
    
    const unsigned int size() const
    {
        return sizeVertices() + sizeElements();
    }
    
    
    const unsigned int sizeVertices() const
    {
        return (unsigned int) M_vertices.size();
    }

    
    const unsigned int sizeElements() const
    {
        return (unsigned int) M_elements.size();
    }

    
private:
    
    std::vector<std::shared_ptr<Vertex> >& M_vertices;
    std::vector<std::shared_ptr<Element> >& M_elements;
    
};

#endif
