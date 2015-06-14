//////////////////////////////////////////////////////////
//  fluid.h                                            //
//  FluidSimulation2D                                   //
//                                                      //
//  Created by nana on 6/4/15.                          //
//  Copyright (c) 2015 na. All rights reserved.         //
//////////////////////////////////////////////////////////

#ifndef FluidSimulation2D_fluid_h
#define FluidSimulation2D_fluid_h

struct Fluidcell {
    int size; // size=N
    float dt; // time interval
    float diff;
    float visc; // viscosity
    
    float *s; // force from outside
    float *density;
    
    float *u; // first coordinate of velocity
    float *v; // second coordinate of velocity
    
    float *u0;
    float *v0;
};

#endif