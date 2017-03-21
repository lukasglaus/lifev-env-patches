//
//  CirculationGridView.cpp
//  Circulation
//
//  Created by Thomas Kummer on 17.05.15.
//  Copyright (c) 2015 Thomas Kummer. All rights reserved.
//

#ifndef CIRCULATIONGRIDVIEW_HPP_
#define CIRCULATIONGRIDVIEW_HPP_

#include <stdio.h>
#include <memory>


class GridView {
public:
    
    GridView() {}

    //GridView(const GridView& gv) = delete; //:
//        M_vertices (gv.M_vertices),
//        M_elements (gv.M_elements)
//    {}
    
    virtual ~GridView() {}

    std::vector<std::shared_ptr<Vertex> >& vertices()
    {
        return M_vertices;
    }
    
    std::vector<std::shared_ptr<Element> >& elements()
    {
        return M_elements;
    }
    
    
private:
    
    std::vector<std::shared_ptr<Vertex> > M_vertices;
    std::vector<std::shared_ptr<Element> > M_elements;
    
};

#endif
