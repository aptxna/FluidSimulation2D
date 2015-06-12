//////////////////////////////////////////////////////////
//  main.cpp                                            //
//  FluidSimulation2D                                   //
//                                                      //
//  Created by nana on 6/4/15.                          //
//  Copyright (c) 2015 na. All rights reserved.         //
//////////////////////////////////////////////////////////

#define N 192 // the resolution size of grid
#define SIZE (N+2)*(N+2) // the size of grid
#define IX(i,j) (i+(N+2)*j) // index of the float array, i and j are the horizontal and vertical component in Cartesian coordinates
#define SWAP(x0,x) {float *temp=x0; x0=x; x=temp;} // x0 is the cause, x is the effect

#include <iostream>
#include <GLUT/glut.h>
#include <math.h>

static float* u;
static float* v;
float dt = .1;
float diff = .1;

/**
 * Set the boundary condition
 * @param b: b=1 horizontal component of velocity on the vertical wall
 *           b=2 vertical component of velocity on the horizontal wall
 * @param var: density or velocity or pressure array
 */
void setBnd(int b,float* var){
    for (int i=1; i<=N; i++) {
        var[IX(0, i)] = b == 1 ? -var[IX(1, i)] : var[IX(1, i)];
        var[IX(N+1, i)] = b == 1 ? -var[IX(N, i)] : var[IX(N, i)];
        var[IX(i, 0)] = b == 2 ? -var[IX(i, 1)] : var[IX(i, 1)];
        var[IX(i, N+1)] = b == 2 ? -var[IX(i, N)] : var[IX(i, N)];
    }
    
    var[IX(0, 0)] = 0.5 * (var[IX(1, 0)] + var[IX(0, 1)]);
    var[IX(0, N+1)] = 0.5 * (var[IX(1, N+1)] + var[IX(0, N)]);
    var[IX(N+1, 0)] = 0.5 * (var[IX(N, 0)] + var[IX(N+1, 1)]);
    var[IX(N+1, N+1)] = 0.5 *(var[IX(N, N+1)] + var[IX(N+1, N)]);
}

/**
 * This is to implement the force term
 * assume that the source for a given frame is provided in the array source[]
 * source[] is filled in by the mouse movement which detects source of density
 */
void addSource(float* cell, float* source, float dt) {
    for (int i=0; i<SIZE; i++) {
        cell[i] += dt * source[i];
    }
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
                current[IX(i, j)] = (previous[IX(i, j)] + a*(current[IX(i-1, j)] + current[IX(i+1, j)]
                                                        + current[IX(i, j-1)] + current[IX(i, j+1)])) / (1 + 4 * a);
            }
        }
    }
    setBnd(b,current);
}

/**
 * This is to implement the advection term
 * @param current: current moment cell density
 * @param previous: previous moment cell density
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
            current[IX(i, j)] = s0 * (t0 * previous[IX(i0, j0)] + t1 * previous[IX(i0, j1)]) +
                                s1 * (t0 * previous[IX(i1, j0)] + t1 * previous[IX(i1, j1)]);
        }
    }
    
    setBnd(b,current);
}

/**
 * Here we have to solve a possion equation which leads to a sparse linear system for the unknow field
 * Gauss-Seidel Relaxation is efficient to solve the equation
 * @param u,v: two components of the 2-dimension velocity field
 * Laplacian p = (...)/h^2 = Div u* = (...)/2h
 * multiply the h^2 to the right side of the equation to make simpler
 */
void project(float *u, float *v) {
    float *p, *div; // p is the pressure, div is the divergence of u*
    float h = 1.0/N;
    
    // get the divergence of u*
    for (int i=1; i<=N; i++) {
        for (int j=1; j<=N; j++) {
            div[IX(i, j)] = 0.5 * h * (u[IX(i+1, j)] - u[IX(i-1, j)] + v[IX(i, j+1)] - v[IX(i, j-1)]);
            p[IX(i, j)] = 0;
        }
    }
    
    // set the divergence and pressure on boundary, should be the same with the nearest fluid cell, so we let b=0
    setBnd(0, div);
    setBnd(0, p);
    
    // Gauss-Seidel Relaxation, solve the poisson equation to get the pressure term
    for (int k=0; k<20; k++) {
        for (int i=1; i<=N; i++) {
            for (int j=1; j<=N; j++) {
                p[IX(i, j)] = (p[IX(i-1, j)] + p[IX(i+1, j)] + p[IX(i, j-1)] + p[IX(i, j+1)] - div[IX(i, j)]) / 4;
            }
        }
        setBnd(0, p);
    }
    
    // since we got the pressure term frome above
    // then replace it in the equation to solve the next time step divergence-free velocity field
    // u^{n+1} = u* - grad p
    for (int i=1; i<=N; i++) {
        for (int j=1; j<=N; j++) {
            u[IX(i, j)] -= 0.5 * (p[IX(i+1, j)] - p[IX(i-1, j)]) / h;
            v[IX(i, j)] -= 0.5 * (p[IX(i, j+1)] - p[IX(i, j-1)]) / h;
        }
    }
    
    // set the boundary condition
    setBnd(1, u);
    setBnd(2, v);
}

/**
 * Update the density during a step of dt
 * @param x0: the given force or source
 * @param x: the density of the cell during a step of dt update, the effect result of source
 * @param u: first component of velocity
 * @param v: second component of velocity
 */
void densStep(float *x, float *x0, float *u, float *v, float diff, float dt) {
    addSource(x, x0, dt);
    SWAP(x0, x); // now swap x0 and x, the effect becomes the cause
    diffuse(0, x, x0, diff, dt);
    SWAP(x0, x);
    advect(0, x, x0, u, v, dt);
}

/**
 * Update the velocity field during a step of dt
 * @param u0, v0: the given velocity field
 * @param u, v: the velocity field during a step of dt update
 * @param dt: time step interval
 * @param visc: viscosity coefficient
 */
void velStep(float *u, float *v, float *u0, float *v0, float visc, float dt) {
    addSource(u, u0, dt);
    addSource(v, v0, dt);
    
    // update the velocity field
    SWAP(u0, u);
    SWAP(v0, v);
    
    diffuse(1, u, u0, visc, dt);
    diffuse(2, v, v0, visc, dt);
    
    // make the velocity field divergence-free for next step
    project(u, v);
    
    // update for next step
    SWAP(u0, u);
    SWAP(v0, v);
    
    // advection for the velocity field itself
    advect(1, u, u0, u0, v0, dt);
    advect(2, v, v0, u0, v0, dt);
    
    // make the velocity field divergence-free for next time step
    project(u, v);
}

int main(int argc, const char * argv[]) {
    // insert code here...
    std::cout << "Hello, World!\n";
    return 0;
}
