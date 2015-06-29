//
//  CirculationIO.hpp
//  Circulation
//
//  Created by Thomas Kummer on 17.04.15.
//  Copyright (c) 2015 Thomas Kummer. All rights reserved.
//

#ifndef CIRCULATIONIO_H_
#define CIRCULATIONIO_H_

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <sstream>

#include "CirculationElementFactory.hpp"
#include "CirculationGridView.hpp"

class CirculationIO {
public:
    
    CirculationIO(){}
    
    virtual ~CirculationIO(){}
    
    static void readGrid(const std::string& filename, GridView& gv)
    {
        std::ifstream file (filename.c_str());
        if (file.is_open())
        {
            std::string line;
            std::string elementType;

            unsigned int i = 0;
            while ( getline (file, line) )
            {
                if ( line.empty() ) continue;
                
                std::istringstream iss(line);
                std::string column;
                
                std::vector<std::string> nodes;
                std::string vertexName;
                std::vector<double> paramElem;
                std::vector<double> paramVert;

                bool properLine (false);
                
                unsigned int j = 0;
                while(iss >> column)
                {
                    if ( column == "//" ) break;
                    
                    if ( column == "#" )
                    {
                        iss >> column;
                        elementType = column;
                        break;
                    }
                    
                    if ( j == 0 ) vertexName = column;
                    else
                    {
                        properLine = true;
                        std::istringstream elementStream(column);
                        double parameter;
                        elementStream >> parameter;
                        paramVert.push_back(parameter);
                    }
                    if ( j < 2 )
                    {
                        nodes.push_back(column);
                    }
                    else
                    {
                        std::istringstream elementStream(column);
                        double parameter;
                        elementStream >> parameter;
                        paramElem.push_back(parameter);
                    }
                    
                    ++j;
                }
                
                if ( properLine )
                {
                    std::unique_ptr<Vertex> vertex( VertexFactory::create(elementType, vertexName, paramVert) );
                    std::unique_ptr<Element> element( ElementFactory::create(elementType, nodes, paramElem) );
                    if ( vertex != nullptr ) gv.vertices().push_back( std::move( vertex ) );
                    if ( element != nullptr ) gv.elements().push_back( std::move( element ) );
                }
                ++i;
            }
            
            file.close();
        }
        
        else
        {
            throw std::runtime_error( "Unable to open file" );
        }
    }
    
//    template<class type>
//    std::vector<type> entityColumn(const std::string& filename, const std::string& entity, const unsigned int& columnNumber) const
//    {
//        std::ifstream file (filename.c_str());
//        if (file.is_open())
//        {
//            std::string line;
//            bool isEntity = false;
//            std::vector<type> entityColumnVector;
//            
//            unsigned int i = 0;
//            while ( getline (file, line) )
//            {
//                if ( line.empty() ) continue;
//                
//                bool columnFound = false;
//                std::istringstream iss(line);
//                std::string column;
//
//                unsigned int j = 0;
//                while(iss >> column)
//                {
//                    if ( column == "//" ) break;
//                    
//                    if ( column == "#" )
//                    {
//                        iss >> column;
//                        if ( column == entity ) isEntity = true;
//                        else isEntity = false;
//                        break;
//                    }
//                    
//                    if ( !isEntity ) break;
//                    
//                    if ( j == columnNumber )
//                    {
//                        std::istringstream elementStream(column);
//                        type element;
//                        elementStream >> element;
//                        
//                        entityColumnVector.push_back( element );
//                        columnFound = true;
//                        break;
//                    }
//                    
//                    ++j;
//                }
//                
//                ++i;
//            }
//            
//            file.close();
//            return entityColumnVector;
//        }
//        
//        else std::cout << "Unable to open file";
//        return std::vector<type>(0);
//    }

    template<class type>
    void exportVector(const std::string& filename, const type& vec, const bool& append = true) const
    {
        if ( !append ) std::ofstream file (filename);
        std::ofstream file (filename, std::ofstream::out | std::ofstream::app);
        
        if (file.is_open())
        {
            file << std::setprecision(8) << std::setw(1) << std::scientific;
            for ( auto i : vec )
            {
                file << i << " ";
            }
            file << "\n";
            file.close();
        }
        else std::cout << "Unable to open file";
    }
    
    
private:
    
    
};


#endif /* CIRCULATIONIO_H_ */