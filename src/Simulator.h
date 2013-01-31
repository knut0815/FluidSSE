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
            for (int x = 0; x < 2; x++, nodePtr += gSizeY_2) {
                for (int y = 0; y < 2; y++, nodePtr += gSizeZ_2) {
                    for (int z = 0; z < 2 ; z++, nodePtr++, phiPtr += 4) {
                        density = _mm_add_ps(density, _mm_mul_ps(_mm_load1_ps(phiPtr), _mm_dp_ps(_mm_load_ps(nodePtr->gx), pdx, 0x7f)));
                    }
                }
            }
            
            float densityScalar;
            _mm_store_ss(&densityScalar, density);
            
            float pressure = min(.25*(densityScalar-2), .5);
            
            nodePtr = &grid[p.c];
            phiPtr = p.phi;
            for (int x = 0; x < 2; x++, nodePtr += gSizeY_2) {
                for (int y = 0; y < 2; y++, nodePtr += gSizeZ_2) {
                    for (int z = 0; z < 2 ; z++, nodePtr++, phiPtr += 4) {
                        //__m128 phi = _mm_and_ps(_mm_load_ps(phiPtr), mask3D);
                        float w = *(phiPtr+3);
                        __m128 wvec = _mm_set_ps1(w);
                        nodePtr->m += w;
                    }
                }
            }
        }
        
        for (int i = 0; i < nParticles; i++) {
            Particle &p = particles[i];
            
            __m128 pu = _mm_load_ps(p.u);
            
            Node *nodePtr = &grid[p.c];
            float *phiPtr = p.phi;
            for (int x = 0; x < 2; x++, nodePtr += gSizeY_2) {
                for (int y = 0; y < 2; y++, nodePtr += gSizeZ_2) {
                    for (int z = 0; z < 2 ; z++, nodePtr++, phiPtr += 4) {
                        //__m128 phi = _mm_and_ps(_mm_load_ps(phiPtr), mask3D);
                        float w = *(phiPtr+3);
                        __m128 wvec = _mm_set_ps1(w);
                        nodePtr->m += w;
                        _mm_store_ps(nodePtr->u, _mm_add_ps(_mm_mul_ps(wvec, pu), _mm_load_ps(nodePtr->u)));
                    }
                }
            }
        }
        
        for (int i = 0; i < gSize; i++) {
            Node &n = grid[i];
            if (n.m > 0) {
                n.invM = 1/n.m;
                _mm_store_ps(n.u, _mm_mul_ps(_mm_load_ps(n.u), _mm_set_ps1(n.invM)));
            }
        }
        
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
            
            _mm_store_ps(p.x, _mm_add_ps(_mm_load_ps(p.x), gu));
        }
    }
};

#endif
