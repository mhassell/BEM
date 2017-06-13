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
     int N = 10;
     
     std::vector<std::vector<double> > *coords = new std::vector<std::vector<double> >(N);
     std::vector<std::vector<double> >& coordsref = *coords;
     std::vector<std::vector<int> >*elts = new std::vector<std::vector<int> >(N);
     std::vector<std::vector<int> >& eref = *elts;
     geometry g(coordsref,eref);
     
     
     std::cout << "Enhanced state: ";
     std::cout << g.enhanced << std::endl;
     
     std::cout << "reference size: " << xr.size() << std::endl;
     
     for(int i = 0; i<N; i++){
     std::cout << g.xcoords[i] << ' ';
     }
     std::cout << std::endl;
     
     // is this by reference or by value?
     xr[1] = 10;
     
     std::cout << "After changing the referenced vector: \n";
     for(int i = 0; i<N; i++){
     std::cout << g.xcoords[i] << ' ';
     }
     std::cout << std::endl;
     
     std::cout << "Testing the normalx vector \n";
     for(size_t i = 0; i<g.normalx.size(); i++ ){
     std::cout << g.normalx[i] << ' ';
     }
     std::cout << std::endl;
     
     std::cout << "Testing the normaly vector \n";
     for(size_t i = 0; i<g.normaly.size(); i++ ){
     std::cout << g.normaly[i] << ' ';
     }
     std::cout << std::endl;
     
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
     
     // delete shit
     delete xp;
     delete yp;
     delete elts;
    
    return 0;
    
}
