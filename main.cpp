//////////////////////////////////////////////////////////
//  main.cpp                                            //
//  FluidSimulation2D                                   //
//                                                      //
//  Created by nana on 6/4/15.                          //
//  Copyright (c) 2015 na. All rights reserved.         //
//////////////////////////////////////////////////////////



#include <iostream>
#include <GLUT/GLUT.h>
#include <math.h>
using namespace std;

#define N 128 // the resolution size of grid
#define M_SIZE (N+2)*(N+2) // the size of grid
#define IX(i,j) ((i)+(N+2)*(j)) // index of the float array, i and j are the horizontal and vertical component in Cartesian coordinates
#define SWAP(x0,x) {float *temp=x0; x0=x; x=temp;} // x0 is the cause, x is the effect

float *u;
float *v;
float *u_init;
float *v_init;
float *fx;
float *fy;
float *dens[3];
float *dens_init[3];
float visc = 10;
float dt = 0.1;
float diff = 0.1;
float xmult=3,ymult=3;
int currentButton=0,currentColor=1;

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
    for (int i=0; i<M_SIZE; i++) {
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
float lerp(float a, float b, float t)
{
    return (1-t)*a + t*b;
}
void advect(int b,float* current,float* previous,float* u,float* v,float dt){
    float dt0 = dt*(float)N;
    
    for (int i=1; i<=N; i++) {
        for (int j=1; j<=N; j++) {
            float x = (float)i - dt0 * u[IX(i,j)];
            float y = (float)j - dt0 * v[IX(i,j)];
            
            // make sure not to flow out the boundary
            if (x < 0.5) x = 0.5;
            if (x > N + 0.5) x = N + 0.5;
            if (y < 0.5) y = 0.5;
            if (y > N + 0.5) y = N + 0.5;
            
            // get the nearest four neighbors
            int i0 = (int)x;
            // int i1 = i0 + 1;
            int j0 = (int)y;
            // int j1 = j0 + 1;
            
            // get the four coefficient in bilinear interpolation
            /*float s0 = i1 - x;
             float s1 = x - i0;
             float t0 = j1 - y;
             float t1 = y - j0;*/
            
            //float s1 = x - i0;
            //float s0 = 1 - s1;
            //float t1 = y - j0;
            //float t0 = 1 - t1;
            //
            //// bilinear interpolation
            //current[IX(i, j)] = s0 * (t0 * previous[IX(i0, j0)] + t1 * previous[IX(i0, j1)]) +
            //                    s1 * (t0 * previous[IX(i1, j0)] + t1 * previous[IX(i1, j1)]);
            
            
            float cx = x - i0;
            float cy = y - j0;
            float v00 = previous[IX(i0,j0)];
            float v01 = previous[IX(i0+1,j0)];
            float v10 = previous[IX(i0,j0+1)];
            float v11 = previous[IX(i0+1,j0+1)];
            
            
            current[IX(i,j)] = 0.99*lerp(lerp(v00,v01,cx),
                                         lerp(v10,v11,cx),
                                         cy);
            
        }
    }
    
    setBnd(b,current);
}

/**
 * Here we have to solve a possion equation which leads to a sparse linear system for the unknow field
 * Gauss-Seidel Relaxation is efficient to solve the equation
 * @param u,v: two components of the 2-dimension velocity field
 * Laplacian p = (...)/h^2 = Div u* = (...)/2h
 * multiply the h^2 to the right side of the equation
 */
void project(float *u, float *v, float *p, float *div) {
    int i, j, k;
    float h = 1.0/(float)N;
    
    // get the divergence of u*
    for (i=1; i<=N; i++) {
        for (j=1; j<=N; j++) {
            div[IX(i, j)] = 0.5 * h * (u[IX(i+1, j)] - u[IX(i-1, j)] + v[IX(i, j+1)] - v[IX(i, j-1)]);
            p[IX(i, j)] = 0;
        }
    }
    
    // set the divergence and pressure on boundary, should be the same with the nearest fluid cell, so we let b=0
    setBnd(0, div);
    setBnd(0, p);
    
    // Gauss-Seidel Relaxation, solve the poisson equation to get the pressure term
    for (k=0; k<200; k++) {
        for (i=1; i<=N; i++) {
            for (j=1; j<=N; j++) {
                p[IX(i, j)] = (p[IX(i-1, j)] + p[IX(i+1, j)] + p[IX(i, j-1)] + p[IX(i, j+1)] - div[IX(i, j)]) / 4.0;
            }
        }
        setBnd(0, p);
    }
    
    // since we got the pressure term frome above
    // then replace it in the equation to solve the next time step divergence-free velocity field
    // u^{n+1} = u* - grad p
    for (i=1; i<=N; i++) {
        for (j=1; j<=N; j++) {
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
void densStep(float* x[],float* x0[],float* u,float* v,float diff,float dt){
    int m;
    for(m=0;m<3;m++){
        
        /* diffuse(0,x[m],x0[m],diff,dt);
         SWAP(x0[m],x[m]);*/
        advect(0,x[m],x0[m],u,v,dt);
        SWAP(x0[m],x[m]);
    }
    
}

/**
 * Update the velocity field during a step of dt
 * @param u0, v0: the given velocity field
 * @param u, v: the velocity field during a step of dt update
 * @param dt: time step interval
 * @param visc: viscosity coefficient
 */
void velStep(float *u, float *v, float *u0, float *v0, float visc, float dt) {
    
    
    //// update the velocity field
    //SWAP(u0, u);
    //SWAP(v0, v);
    //
    //diffuse(1, u, u0, visc, dt);
    //diffuse(2, v, v0, visc, dt);
    //
    //// make the velocity field divergence-free for next step
    //project(u, v, u0, v0);
    
    
    memset(u, sizeof(float)*M_SIZE, 0);
    memset(v, sizeof(float)*M_SIZE, 0);
    // advection for the velocity field itself
    advect(1, u, u0, u0, v0, dt);
    advect(2, v, v0, u0, v0, dt);
    
    //addSource(u, fx, dt);
    //addSource(v, fy, dt);
    //memset(fx, sizeof(float)*M_SIZE, 0);
    //memset(fy, sizeof(float)*M_SIZE, 0);
    
    // make the velocity field divergence-free for next time step
    project(u, v, u0, v0);
    // update for next step
    //SWAP(u0, u);
    //SWAP(v0, v);
    memcpy(u0,u, sizeof(float)*M_SIZE);
    memcpy(v0,v, sizeof(float)*M_SIZE);
    
    
    
}



/**********************************************************************************************
 *                                                                                            *
 *                                                                                            *
 *                                   Display Part Starts Here                                 *
 *                                                                                            *
 *                                                                                            *
 **********************************************************************************************/

void draw_dens(void){
    velStep(u,v,u_init,v_init,visc,dt);
    densStep(dens,dens_init,u_init,v_init,diff,dt);
    int x,y,m;
    float c[3];
    glClear(GL_COLOR_BUFFER_BIT);
    
    for(y=0;y<=N;y++){
        for(x=0;x<=N;x++){
            for(m=0;m<3;m++) c[m] = dens_init[m][IX(x,y)]/255.0;
            if(c[0]==0 && c[1]==0 && c[2]==0) continue;
            glColor3f(c[0],c[1],c[2]);
            
            glRecti(x,N-y,x+1,N-y+1);
        }
    }
    glutSwapBuffers();
    glutPostRedisplay();
}

void reshape(int w,int h){
    glViewport(0,0,(GLsizei)w,(GLsizei)h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0,N,0,N,-1,1);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    
    xmult = w/(GLfloat)N;
    ymult = h/(GLfloat)N;
}

int mouseButtons;
int mouseOldX;
int mouseOldY;
/**********************************************************************************************
 *                                                                                            *
 *                                                                                            *
 *                                Interaction Part Starts Here                                *
 *                                                                                            *
 *                                                                                            *
 **********************************************************************************************/

void mouse(int button,int state,int x,int y){
    if (state == GLUT_DOWN)
    {
        mouseButtons |= 1<<button;
    }
    else if (state == GLUT_UP)
    {
        mouseButtons = 0;
    }
    
    mouseOldX = x;
    mouseOldY = y;
    glutPostRedisplay();
}

void mouseMove(int x,int y){
    float dx, dy;
    dx = (float)(x - mouseOldX);
    dy = (float)(y - mouseOldY);
    
    if (mouseButtons == 1)
    {
        float coordx = x/500.0;
        float coordy = y/500.0;
        int buffer_i = coordx * N;
        int buffer_j = coordy * N;
        printf("%d,%d\n", buffer_i,buffer_j);
        for (int ii=max(buffer_i-3,1);ii<=min(buffer_i+3,N);ii++)
            for (int jj=max(buffer_j-3,1);jj<=min(buffer_j+3,N);jj++)
            {
                dens_init[0][IX(ii,jj)] = rand()%256;
                dens_init[1][IX(ii,jj)] = rand()%256;
                dens_init[2][IX(ii,jj)] = rand()%256;
                u_init[IX(ii,jj)] += dt*0.5*dx;
                v_init[IX(ii,jj)] += dt*0.5*dy;
            }
        glutPostRedisplay();
    }
    else
    {
        memset(fx, sizeof(float)*M_SIZE, 0);
        memset(fy, sizeof(float)*M_SIZE, 0);
        glutPostRedisplay();
    }
    
    mouseOldX = x;
    mouseOldY = y;
    
}

void fluidMainLoop(void){
    velStep(u,v,u_init,v_init,visc,dt);
    densStep(dens,dens_init,u_init,v_init,diff,dt);
    
}

int main(int argc, char* argv[]){
    int i, j;
    
    u = new float[M_SIZE];
    v = new float[M_SIZE];
    fx= new float[M_SIZE];
    fy= new float[M_SIZE];
    u_init = new float[M_SIZE];
    v_init = new float[M_SIZE];
    
    for(j=0;j<3;j++){
        dens[j] = new float[M_SIZE];
        dens_init[j] = new float[M_SIZE];
    }
    
    if(dens[0]==NULL || dens[1]==NULL || dens[2]==NULL || dens_init[0]==NULL || dens_init[1]==NULL || dens_init[2]==NULL) {
        cout << "Memory Allocation Failure.";
        exit(1);
    }
    
    for(i=0;i<M_SIZE;i++){
        u[i] = 0.0;
        v[i] = 0.0;
        u_init[i] = 0.0;
        v_init[i] = 0.0;
        for(j=0;j<3;j++){
            dens[j][i] = 0.0;
            dens_init[j][i] = 0.0;
        }
    }
    
    // GLUT
    glutInit(&argc,argv);
    glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGB);
    glutInitWindowSize(500, 500);
    glutInitWindowPosition(100, 100);
    glutCreateWindow("Fluid Simulation 2D");
    glClearColor(0,0,0,0);
    glShadeModel(GL_FLAT);
    glOrtho(0,N,0,N,-1,1);
    glutDisplayFunc(draw_dens);
    glutReshapeFunc(reshape);
    glutMouseFunc(mouse);
    glutMotionFunc(mouseMove);
    glutMainLoop();
    
    
    free(u);
    free(v);
    free(u_init);
    free(v_init);
    free(dens);
    free(dens_init);
}