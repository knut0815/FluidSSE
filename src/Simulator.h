//
//  Simulator.h
//  FluidSSE
//
//  Created by Ruby on 1/30/13.
//
//

#ifndef FluidSSE_Simulator_h
#define FluidSSE_Simulator_h

#include <x86intrin.h>
#include "mm_malloc.h"
#include "Particle.h"
#include "Node.h"

#define maxParticles 1000000

class Simulator {
    __m128 dxSub[8];
    __m128 lowBound, highBound, lowBoundS, highBoundS;
public:
    int gSizeX, gSizeY, gSizeZ, gSizeY_2, gSizeZ_2, gSize;
    int nParticles;
    Particle *particles;
    Node *grid;
    
    void InitializeGrid(int gridSizeX, int gridSizeY, int gridSizeZ) {
        gSizeX = gridSizeX;
        gSizeY = gridSizeY;
        gSizeZ = gridSizeZ;
        gSizeY_2 = gSizeY*gSizeZ-2*gSizeZ;
        gSizeZ_2 = gSizeZ-2;
        gSize = gSizeX*gSizeY*gSizeZ;
        grid = (Node*)_mm_malloc(gSize*sizeof(Node), 16);
        for (int i = 0; i < gSize; i++) {
            grid[i].Clear();
        }
        particles = (Particle*)_mm_malloc(maxParticles*sizeof(Particle), 16);
        nParticles = 0;
        
        for (int x = 0; x < 2; x++) {
            for (int y = 0; y < 2; y++) {
                for (int z = 0; z < 2; z++) {
                    dxSub[x*4+y*2+z] = _mm_set_ps(-1, z, y, x);
                }
            }
        }
        
        lowBound = _mm_set_ps(0.0f, .1f, .1f, .1f);
        highBound = _mm_set_ps(0.0f, gSizeZ-1.1, gSizeY-1.1, gSizeX-1.1);
        lowBoundS = _mm_set_ps(0.0f, 5.0f, 5.0f, 5.0f);
        highBoundS = _mm_set_ps(0.0f, gSizeZ-6, gSizeY-6, gSizeX-6);
    }
    
    void AddParticles() {
        for (int i = 0; i < 100; i++) {
            for (int j = 0; j < 100; j++) {
                for (int k = 0; k < 100; k++) {
                    particles[nParticles++].Initialize(i*.5+25, j*.5+25, k*.5+25);
                }
            }
        }
    }
    
    void Update() {
        // Clear the grid
        __m128 zero = _mm_setzero_ps();
        for (int i = 0; i < gSize; i++) {
            Node &n = grid[i];
            if (n.m > 0) {
                n.m = n.invM = 0;
                _mm_store_ps(n.gx, zero);
                _mm_store_ps(n.u, zero);
                _mm_store_ps(n.ax, zero);
            }
        }
        
        // Calculate particle kernels, and add density and density gradients to the grid
        __m128 cmul = _mm_set_ps(0, 1, gSizeZ, gSizeY*gSizeZ);
        __m128 one = _mm_set_ps1(1);
        __m128 negOne = _mm_set_ps1(-1);
        for (int i = 0; i < nParticles; i++) {
            Particle &p = particles[i];
            
            __m128 px = _mm_load_ps(p.x);
            __m128 cx = _mm_floor_ps(px);
            
            float pc;
            _mm_store_ss(&pc, _mm_dp_ps(cx, cmul, 0x7f));
            p.c = (int)pc;
            
            __m128 dx = _mm_sub_ps(px, cx);
            _mm_store_ps(p.dx, dx);
            
            __m128 temp;
            __m128 u[2], v[2], w[2];
            __m128 xdup = _mm_shuffle_ps(dx, dx, 0x00);
            u[0] = _mm_move_ss(_mm_sub_ps(one, xdup), one);
            u[1] = _mm_move_ss(xdup, negOne);
            temp = _mm_move_ss(_mm_sub_ps(one, dx), one);
            v[0] = _mm_shuffle_ps(temp, temp, _MM_SHUFFLE(1, 1, 0, 1));
            w[0] = _mm_shuffle_ps(temp, temp, _MM_SHUFFLE(2, 0, 2, 2));
            temp = _mm_move_ss(dx, negOne);
            v[1] = _mm_shuffle_ps(temp, temp, _MM_SHUFFLE(1, 1, 0, 1));
            w[1] = _mm_shuffle_ps(temp, temp, _MM_SHUFFLE(2, 0, 2, 2));
            
            Node *nodePtr = &grid[p.c];
            float *phiPtr = p.phi;
            for (int x = 0; x < 2; x++, nodePtr += gSizeY_2) {
                for (int y = 0; y < 2; y++, nodePtr += gSizeZ_2) {
                    __m128 phixy = _mm_mul_ps(u[x], v[y]);
                    for (int z = 0; z < 2 ; z++, nodePtr++, phiPtr += 4) {
                        __m128 phi = _mm_mul_ps(phixy, w[z]);
                        
                        _mm_store_ps(phiPtr, phi);
                        _mm_store_ps(nodePtr->gx, _mm_add_ps(_mm_load_ps(nodePtr->gx), phi));
                    }
                }
            }
        }
        
        // Sum particle density from grid, and add pressure and elastic forces to grid
        __m128 mask3D = _mm_castsi128_ps(_mm_set_epi32(0x00000000, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF));
        for (int i = 0; i < nParticles; i++) {
            Particle &p = particles[i];
            
            __m128 density = _mm_setzero_ps();
            __m128 pdx = _mm_load_ps(p.dx);
            
            Node *nodePtr = &grid[p.c];
            float *phiPtr = &p.phi[3];
            int dxi = 0;
            for (int x = 0; x < 2; x++, nodePtr += gSizeY_2) {
                for (int y = 0; y < 2; y++, nodePtr += gSizeZ_2) {
                    for (int z = 0; z < 2 ; z++, nodePtr++, phiPtr += 4, dxi++) {
                        density = _mm_add_ps(density, _mm_mul_ps(_mm_load1_ps(phiPtr), _mm_dp_ps(_mm_load_ps(nodePtr->gx), _mm_sub_ps(pdx, dxSub[dxi]), 0xff)));
                    }
                }
            }
            
            float densityScalar;
            _mm_store_ss(&densityScalar, density);
            
            float pressure = min(.0625*(densityScalar-8), .5);
            
            __m128 px = _mm_load_ps(p.x);
            __m128 wallforce = _mm_add_ps(_mm_max_ps(zero, _mm_sub_ps(lowBoundS, px)), _mm_min_ps(zero, _mm_sub_ps(highBoundS, px)));
            
            //wallforce = _mm_setzero_ps();
            
            nodePtr = &grid[p.c];
            phiPtr = p.phi;
            for (int x = 0; x < 2; x++, nodePtr += gSizeY_2) {
                for (int y = 0; y < 2; y++, nodePtr += gSizeZ_2) {
                    for (int z = 0; z < 2 ; z++, nodePtr++, phiPtr += 4) {
                        __m128 phi = _mm_and_ps(_mm_load_ps(phiPtr), mask3D);
                        float w = *(phiPtr+3);
                        __m128 wvec = _mm_set_ps1(w);
                        
                        __m128 a = _mm_sub_ps(_mm_mul_ps(wvec, wallforce), _mm_mul_ps(phi, _mm_set1_ps(pressure)));
                        _mm_store_ps(nodePtr->ax, _mm_add_ps(_mm_load_ps(nodePtr->ax), a));
                    }
                }
            }
        }
        
        // Divide grid accelerations by mass
        for (int i = 0; i < gSize; i++) {
            Node &n = grid[i];
            n.m = n.gx[3];
            if (n.m > 0) {
                n.invM = 1/n.m;
                _mm_store_ps(n.ax, _mm_mul_ps(_mm_load_ps(n.ax), _mm_set_ps1(n.invM)));
                n.ax[1] -= .01f;
            }
        }
        
        // Accelerate particles and interpolate velocity back to grid
        for (int i = 0; i < nParticles; i++) {
            Particle &p = particles[i];
            
            __m128 ga = _mm_setzero_ps();
            Node *nodePtr = &grid[p.c];
            float *phiPtr = p.phi;
            for (int x = 0; x < 2; x++, nodePtr += gSizeY_2) {
                for (int y = 0; y < 2; y++, nodePtr += gSizeZ_2) {
                    for (int z = 0; z < 2 ; z++, nodePtr++, phiPtr += 4) {
                        //__m128 phi = _mm_and_ps(_mm_load_ps(phiPtr), mask3D);
                        float w = *(phiPtr+3);
                        __m128 wvec = _mm_set_ps1(w);
                        ga = _mm_add_ps(ga, _mm_mul_ps(wvec, _mm_load_ps(nodePtr->ax)));
                    }
                }
            }
            
            __m128 pu = _mm_add_ps(_mm_load_ps(p.u), ga);
            _mm_store_ps(p.u, pu);
            
            nodePtr = &grid[p.c];
            phiPtr = p.phi;
            for (int x = 0; x < 2; x++, nodePtr += gSizeY_2) {
                for (int y = 0; y < 2; y++, nodePtr += gSizeZ_2) {
                    for (int z = 0; z < 2 ; z++, nodePtr++, phiPtr += 4) {
                        //__m128 phi = _mm_and_ps(_mm_load_ps(phiPtr), mask3D);
                        float w = *(phiPtr+3);
                        __m128 wvec = _mm_set_ps1(w);
                        _mm_store_ps(nodePtr->u, _mm_add_ps(_mm_mul_ps(wvec, pu), _mm_load_ps(nodePtr->u)));
                    }
                }
            }
        }
        
        // Divide grid velocity by mass
        for (int i = 0; i < gSize; i++) {
            Node &n = grid[i];
            if (n.m > 0) {
                _mm_store_ps(n.u, _mm_mul_ps(_mm_load_ps(n.u), _mm_set_ps1(n.invM)));
            }
        }
        
        // Advance particles
        for (int i = 0; i < nParticles; i++) {
            Particle &p = particles[i];
            
            __m128 gu = _mm_setzero_ps();
            
            Node *nodePtr = &grid[p.c];
            float *phiPtr = p.phi;
            for (int x = 0; x < 2; x++, nodePtr += gSizeY_2) {
                for (int y = 0; y < 2; y++, nodePtr += gSizeZ_2) {
                    for (int z = 0; z < 2 ; z++, nodePtr++, phiPtr += 4) {
                        //__m128 phi = _mm_and_ps(_mm_load_ps(phiPtr), mask3D);
                        float w = *(phiPtr+3);
                        __m128 wvec = _mm_set_ps1(w);
                        gu = _mm_add_ps(gu, _mm_mul_ps(wvec, _mm_load_ps(nodePtr->u)));
                    }
                }
            }
            
            __m128 px = _mm_add_ps(_mm_load_ps(p.x), gu);
            px = _mm_min_ps(_mm_max_ps(lowBound, px), highBound);
            
            _mm_store_ps(p.x, px);
            _mm_store_ps(p.u, gu);
        }
    }
};

#endif
