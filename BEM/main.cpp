//
//  main.cpp
//  BEM
//
//  Created by Matthew Hassell on 5/30/17.
//  Copyright © 2017 Matthew Hassell. All rights reserved.
//

#include <iostream>
#include <vector>
#include <array>
#include "geometry.hpp"

int main(){
     
    // first, test the basic constructor
    int N = 4;
    
    std::vector<std::vector<double> > *coords = new std::vector<std::vector<double> >(N);
    std::vector<std::vector<double> >& coordsref = *coords;
    std::vector<std::vector<int> >*elts = new std::vector<std::vector<int> >(N);
    std::vector<std::vector<int> >& eref = *elts;
    
    coordsref = {{0,0},{1,0},{0.8,0.8},{0.2,1}};
    eref = {{0,1},{2,3},{1,2},{3,0}};
    
    geometry g(coordsref,eref);
    
    std::cout << "Testing the normal vector \n";
    for(size_t i = 0; i<g.normals.size(); i++ ){
        for(int j = 0; j<2; j++){
            std::cout << g.normals[i][j] << ' ';
        }
        std::cout << std::endl;
    }
    
    
    std::cout << "Testing the lengths vector \n";
    for(size_t i = 0; i<g.lengths.size(); i++ ){
        std::cout << g.lengths[i] << ' ';
    }
    std::cout << std::endl;
    
    std::cout << "Testing the prev vector \n";
    for(size_t i = 0; i<g.prev.size(); i++ ){
        std::cout << g.prev[i] << ' ';
    }
    std::cout << std::endl;
    
    std::cout << "Testing the next vector \n";
    for(size_t i = 0; i<g.next.size(); i++ ){
        std::cout << g.next[i] << ' ';
    }
    std::cout << std::endl;
    
    // delete the new'd stuff
    delete coords;
    delete elts;
    
    return 0;
    
}
