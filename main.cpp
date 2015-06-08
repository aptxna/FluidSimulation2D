//////////////////////////////////////////////////////////
//  main.cpp                                            //
//  FluidSimulation2D                                   //
//                                                      //
//  Created by nana on 6/4/15.                          //
//  Copyright (c) 2015 na. All rights reserved.         //
//////////////////////////////////////////////////////////

#define N 150
#define SIZE (N+2)*(N+2)
#define IX(i,j) (i+(N+2)*j)
#define SWAP(array1, array2) {float *temp=array1; array1=array2; array2=temp;}

#include <iostream>
#include <GLUT/glut.h>
#include <math.h>

static float* u;
static float* v;
float dt = .1;
float diff = .1;

/**
 * this is to implement the force term
 * assume that the source for a given frame is provided in the array source[]
 * source[] is filled in by the mouse movement which detects source of density
 */
void addSource(float* cell, float* source, float dt) {
    for (int i=0; i<SIZE; i++) {
        cell[i] += dt * source[i];
    }
}

/**
 * set the boundary condition
 */
void setBnd(int b,float* x){

}

/**
 * This is to implement the diffusion term
 * @param current: current moument cell density
 * @param previous: previous moument cell density
 * @param diff: diffusion coefficient
 * @param dt: time interval
 * use Gauss-Seidel Relaxation
 * iteration times: k=20
 */
void diffuse(int b, float* current, float* previous, float diff, float dt) {
    float a = dt * diff * N * N;
    for (int k=0; k<20; k++) {
        for (int i=1; i<=N; i++) {
            for (int j=1; j<=N; j++) {
                current[IX(i,j)] = (previous[IX(i,j)] + a*(current[IX(i-1,j)] + current[IX(i+1,j)]
                                                        + current[IX(i,j-1)] + current[IX(i,j+1)]))/(1+4*a);
            }
        }
    }
    setBnd(b,current);
}

/**
 * This is to implement the advection term
 * @param current: current moument cell density
 * @param previous: previous moument cell density
 * @param u: the first coordinate of velocity field
 * @param v: the second coordinate of velocity field
 * @param dt: time interval
 * use Bilinear interpolation
 */
void advect(int b,float* current,float* previous,float* u,float* v,float dt){
    for (int i=1; i<=N; i++) {
        for (int j=1; j<=N; j++) {
            float x = i - dt * u[IX(i,j)];
            float y = j - dt * v[IX(i,j)];
            
            // make sure not to flow out the boundary
            if (x < 0.5) x = 0.5;
            if (x > N + 0.5) x = N + 0.5;
            if (y < 0.5) y = 0.5;
            if (y > N + 0.5) y = N + 0.5;
            
            // get the nearest four neighbors
            int i0 = (int)x;
            int i1 = i0 + 1;
            int j0 = (int)y;
            int j1 = j0 + 1;
            
            // get the four coefficient in bilinear interpolation
            float s0 = i1 - x;
            float s1 = x - i0;
            float t0 = j1 - y;
            float t1 = y - j0;
            
            // bilinear interpolation
            current[IX(i,j)] = s0 * (t0 * previous[IX(i0, j0)] + t1 * previous[IX(i0, j1)]) +
                                s1 * (t0 * previous[IX(i1, j0)] + t1 * previous[IX(i1, j1)]);
        }
    }
    
    setBnd(b,current);
}

/**
 * this is the routine to update density during a step of dt
 */
void densStep(float *u, float *v, float diff, float dt) {
    
}

int main(int argc, const char * argv[]) {
    // insert code here...
    std::cout << "Hello, World!\n";
    return 0;
}
