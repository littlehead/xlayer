//
//  main.cpp
//  xlayer-module
//
//  Created by Yuping Dong on 7/6/15.
//  Copyright (c) 2015 Yuping Dong. All rights reserved.
//

#include <iostream>
#include "xlayer-modulel.hpp"
#include "types.hpp"

//using namespace ns3;

int main (int argc, char * argv[])
{
    CrossLayer cr(0.0008, 4, 0.33, 0.33, 0.33);
    RETOPT *opti = new RETOPT;
    float PhyState = 133;
    float MacState = 1;
    APPST AppState = {0,4};
    
    opti = cr.GetOpt(PhyState, MacState, AppState);
    
    std::cout << "State value of state (" << PhyState << ", " << MacState << ", { " << AppState.x << ", " << AppState.y <<"}) is " << opti->value << std::endl;
    std::cout << "Optimal external action 1 is " << opti->extaction1 << std::endl;
    std::cout << "Optimal external action 2 is " << opti->extaction2 << std::endl;
    std::cout << "Optimal external action 3 is " << opti->extaction3 << std::endl;
    std::cout << "Optimal internal action 1 is " << opti->intaction1 << std::endl;
    
}
