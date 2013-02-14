//
//  Simulator.h
//  FluidSSE
//
//  Created by Ruby on 1/30/13.
//
//

#ifndef FluidSSE_Simulator_h
#define FluidSSE_Simulator_h

#define numMaterials 4

#include <x86intrin.h>
#include "Material.h"
#include "Particle.h"
#include "Node.h"

#define maxParticles 4000000

class Simulator : public ofThread {
    __m128 dxSub[8];
    __m128 lowBound, highBound, lowBoundS, highBoundS;
public:
    int gSizeX, gSizeY, gSizeZ, gSizeY_2, gSizeZ_2, gSize;
    int nParticles;
    Particle *particles;
    Node *grid;
    Material materials[numMaterials];
    int loadedMaterial;
    
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
        
        lowBound = _mm_set_ps(0.0f, .001f, .001f, .001f);
        highBound = _mm_set_ps(0.0f, gSizeZ-1.001, gSizeY-1.001, gSizeX-1.001);
        lowBoundS = _mm_set_ps(0.0f, 2, 2, 2);
        highBoundS = _mm_set_ps(0.0f, gSizeZ-3, gSizeY-3, gSizeX-3);
        
        //stiffness, restDensity, bulkViscosity, elasticity, viscosity, smoothing, yieldPoint, yieldRate
        materials[0].Initialize(1, 2, 1, 0, .05, 0, .5);
        materials[1].Initialize(.0625, 8, .5, .5, .1, .1, .5);
        loadedMaterial = -1;
    }
    
    void AddParticles() {
        for (int i = 0; i < 160; i++) {
            for (int j = 0; j < 120; j++) {
                for (int k = 0; k < 40; k++) {
                    particles[nParticles++].Initialize(i*.5+3, j*.5+3, k*.5+3, 0);
                }
            }
        }
        for (int l = 0; l < 14; l++) {
        for (int i = 0; i < 10; i++) {
            for (int j = 0; j < 120; j++) {
                for (int k = 0; k < 60; k++) {
                    particles[nParticles++].Initialize(i*.5+3+l*10, j*.5+3, k*.5+23, 1);
                }
            }
        }
        }
    }
    
    void threadedFunction() {
        // Clear the grid
        __m128 zero = _mm_setzero_ps();
        #pragma omp parallel for
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
        #pragma omp parallel for
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
        loadedMaterial = -1;
        __m128 stiffness, restDensity, elasticity, viscosity;
        __m128 maxPressure = _mm_set1_ps(.5);
        float elasticityS, viscosityS, bulkViscosityS;
        #pragma omp parallel for
        for (int i = 0; i < nParticles; i++) {
            Particle &p = particles[i];
            
            if (p.material != loadedMaterial) {
                loadedMaterial = p.material;
                stiffness = _mm_set1_ps(materials[p.material].stiffness);
                restDensity = _mm_set1_ps(materials[p.material].restDensity);
                elasticity = _mm_set1_ps(materials[p.material].elasticity);
                viscosity = _mm_set1_ps(materials[p.material].viscosity);
                
                elasticityS = materials[p.material].elasticity;
                viscosityS = materials[p.material].viscosity;
                bulkViscosityS = materials[p.material].bulkViscosity;
            }
            
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
            
            __m128 pressure = _mm_min_ps(_mm_mul_ps(stiffness, _mm_sub_ps(density, restDensity)), maxPressure);
            
            __m128 px = _mm_load_ps(p.x);
            __m128 wallforce = _mm_add_ps(_mm_max_ps(zero, _mm_sub_ps(lowBoundS, px)), _mm_min_ps(zero, _mm_sub_ps(highBoundS, px)));
            
            __m128 e0 = _mm_add_ps(_mm_mul_ps(elasticity, _mm_load_ps(&p.strain[0])), _mm_mul_ps(viscosity, _mm_load_ps(&p.stress[0])));
            __m128 e1 = _mm_add_ps(_mm_mul_ps(elasticity, _mm_load_ps(&p.strain[4])), _mm_mul_ps(viscosity, _mm_load_ps(&p.stress[4])));
            __m128 e2 = _mm_add_ps(_mm_mul_ps(elasticity, _mm_load_ps(&p.strain[8])), _mm_mul_ps(viscosity, _mm_load_ps(&p.stress[8])));
            
            float trace = elasticityS*(p.strain[0]+p.strain[5]+p.strain[10])/3+(viscosityS-bulkViscosityS)*(p.stress[0]+p.stress[5]+p.stress[10])/3;
            
            pressure = _mm_sub_ps(pressure, _mm_set1_ps(trace));
            
            //e0 = _mm_mul_ps(pressure, e0);
            //e1 = _mm_mul_ps(pressure, e1);
            //e2 = _mm_mul_ps(pressure, e2);
            nodePtr = &grid[p.c];
            phiPtr = p.phi;
            for (int x = 0; x < 2; x++, nodePtr += gSizeY_2) {
                for (int y = 0; y < 2; y++, nodePtr += gSizeZ_2) {
                    for (int z = 0; z < 2 ; z++, nodePtr++, phiPtr += 4) {
                        __m128 phi = _mm_and_ps(_mm_load_ps(phiPtr), mask3D);
                        float w = *(phiPtr+3);
                        __m128 wvec = _mm_set_ps1(w);
                        
                        __m128 a = _mm_sub_ps(_mm_mul_ps(wvec, wallforce), _mm_mul_ps(phi, pressure));
                        a = _mm_sub_ps(a, _mm_add_ps(_mm_add_ps(_mm_mul_ps(e0, _mm_set1_ps(phiPtr[0])), _mm_mul_ps(e1, _mm_set1_ps(phiPtr[1]))), _mm_mul_ps(e2, _mm_set1_ps(phiPtr[2]))));
                        _mm_store_ps(nodePtr->ax, _mm_add_ps(_mm_load_ps(nodePtr->ax), a));
                    }
                }
            }
        }
        
        // Divide grid accelerations by mass
        #pragma omp parallel for
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
        #pragma omp parallel for
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
            
            __m128 px = _mm_load_ps(p.x);
            __m128 pu = _mm_add_ps(_mm_load_ps(p.u), ga);
            __m128 pxn = _mm_add_ps(px, pu);
            pxn = _mm_min_ps(_mm_max_ps(pxn, lowBound), highBound);
            pu = _mm_sub_ps(pxn, px);
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
        #pragma omp parallel for
        for (int i = 0; i < gSize; i++) {
            Node &n = grid[i];
            if (n.m > 0) {
                _mm_store_ps(n.u, _mm_mul_ps(_mm_load_ps(n.u), _mm_set_ps1(n.invM)));
            }
        }
        
        // Advance particles
        __m128 Lmul0 = _mm_set_ps(0, .5, .5, 1);
        __m128 Lmul1 = _mm_set_ps(0, .5, 1, .5);
        __m128 Lmul2 = _mm_set_ps(0, 1, .5, .5);
        __m128 LTmul0 = _mm_set_ps(0, .5, .5, 0);
        __m128 LTmul1 = _mm_set_ps(0, .5, 0, .5);
        __m128 LTmul2 = _mm_set_ps(0, 0, .5, .5);
        __m128 smoothing;
        float yieldRate;
        bool elastic;
        loadedMaterial = -1;
        #pragma omp parallel for
        for (int i = 0; i < nParticles; i++) {
            Particle &p = particles[i];
            
            if (p.material != loadedMaterial) {
                loadedMaterial = p.material;
                smoothing = _mm_set1_ps(materials[p.material].smoothing);
                yieldRate = materials[p.material].yieldRate;
                elastic = materials[p.material].elasticity > 0;
            }

            
            __m128 gu = _mm_setzero_ps();
            __m128 L0, L1, L2;
            L0 = L1 = L2 = _mm_setzero_ps();
            
            Node *nodePtr = &grid[p.c];
            float *phiPtr = p.phi;
            for (int x = 0; x < 2; x++, nodePtr += gSizeY_2) {
                for (int y = 0; y < 2; y++, nodePtr += gSizeZ_2) {
                    for (int z = 0; z < 2 ; z++, nodePtr++, phiPtr += 4) {
                        __m128 phi = _mm_and_ps(_mm_load_ps(phiPtr), mask3D);
                        float w = *(phiPtr+3);
                        __m128 wvec = _mm_set_ps1(w);
                        __m128 nu = _mm_load_ps(nodePtr->u);
                        gu = _mm_add_ps(gu, _mm_mul_ps(wvec, nu));
                        L0 = _mm_add_ps(L0, _mm_mul_ps(phi, _mm_shuffle_ps(nu, nu, 0x00)));
                        L1 = _mm_add_ps(L1, _mm_mul_ps(phi, _mm_shuffle_ps(nu, nu, 0x55)));
                        L2 = _mm_add_ps(L2, _mm_mul_ps(phi, _mm_shuffle_ps(nu, nu, 0xAA)));
                    }
                }
            }
            
            // integrate stress tensor
            __m128 LT0 = _mm_shuffle_ps(L1, L2, _MM_SHUFFLE(0, 0, 0, 0));
            __m128 LT1 = _mm_shuffle_ps(L0, L2, _MM_SHUFFLE(1, 1, 1, 1));
            __m128 LT2 = _mm_shuffle_ps(L0, L1, _MM_SHUFFLE(2, 2, 2, 2));
            LT2 = _mm_shuffle_ps(LT2, LT2, _MM_SHUFFLE(0, 0, 2, 0));
            
            L0 = _mm_add_ps(_mm_mul_ps(Lmul0, L0), _mm_mul_ps(LTmul0, LT0));
            L1 = _mm_add_ps(_mm_mul_ps(Lmul1, L1), _mm_mul_ps(LTmul1, LT1));
            L2 = _mm_add_ps(_mm_mul_ps(Lmul2, L2), _mm_mul_ps(LTmul2, LT2));
            
            _mm_store_ps(&p.stress[0], L0);
            _mm_store_ps(&p.stress[4], L1);
            _mm_store_ps(&p.stress[8], L2);
            
            if (elastic) {
                L0 = _mm_add_ps(_mm_load_ps(&p.strain[0]), L0);
                L1 = _mm_add_ps(_mm_load_ps(&p.strain[4]), L1);
                L2 = _mm_add_ps(_mm_load_ps(&p.strain[8]), L2);
                
                __m128 norm = _mm_add_ps(_mm_add_ps(_mm_dp_ps(L0, L0, 0xff), _mm_dp_ps(L1, L1, 0xff)), _mm_dp_ps(L2, L2, 0xff));
                float normScalar;
                _mm_store_ss(&normScalar, norm);
                if (normScalar > 1) {
                    normScalar = sqrtf(normScalar);
                    norm = _mm_set1_ps(yieldRate*(normScalar-1)/normScalar);
                    L0 = _mm_sub_ps(L0, _mm_mul_ps(L0, norm));
                    L1 = _mm_sub_ps(L1, _mm_mul_ps(L1, norm));
                    L2 = _mm_sub_ps(L2, _mm_mul_ps(L2, norm));
                }
                
                _mm_store_ps(&p.strain[0], L0);
                _mm_store_ps(&p.strain[4], L1);
                _mm_store_ps(&p.strain[8], L2);
            }
            __m128 px = _mm_add_ps(_mm_load_ps(p.x), gu);
            px = _mm_min_ps(_mm_max_ps(lowBound, px), highBound);
            
            _mm_store_ps(p.x, px);
            __m128 pu = _mm_load_ps(p.u);
            _mm_store_ps(p.u, _mm_add_ps(pu, _mm_mul_ps(smoothing, _mm_sub_ps(gu, pu))));
        }
        static float time = 0;
        time += .01;
        float s = sinf(time);
        float c = cosf(time);
        if (nParticles < 4000000) {
            for (int i = 0; i < 8; i++) {
                for (int j = 0; j < 8; j++) {
                    for (int k = 0; k < 10; k++) {
                        for (int l = 0; l < 10; l++) {
                            float x = (i-4)*10+k*.7;
                            float y = (j-4)*10+l*.7;
                            float tx = c*x-s*y;
                            float ty = s*x+c*y;
                            particles[nParticles++].Initialize(3, 72+tx, 60+ty, 1, 0, 0, 0);
                        }
                    }
                }
            }
        }
    }
};

#endif
