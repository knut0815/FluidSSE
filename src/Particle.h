//
//  Particle.h
//  FluidSSE
//
//  Created by Ruby on 1/30/13.
//
//

#ifndef FluidSSE_Particle_h
#define FluidSSE_Particle_h

struct Particle {
    float x[4];
    float dx[4];
    float u[4];
    float phi[4*8];
    float stress[12];
    float strain[12];
    int c, material;
    int pad[2];
    void Initialize(float px, float py, float pz, int material) {
        x[0] = px;
        x[1] = py;
        x[2] = pz;
        x[3] = 0;
        __m128 zero = _mm_setzero_ps();
        _mm_store_ps(u, zero);
        for (int i = 0; i < 8; i++) {
            _mm_store_ps(&phi[i*4], zero);
        }
        this->material = material;
    }
    void Initialize(float px, float py, float pz, float ux, float uy, float uz, int material) {
        x[0] = px;
        x[1] = py;
        x[2] = pz;
        x[3] = 0;
        u[0] = ux;
        u[1] = uy;
        u[2] = uz;
        u[3] = 0;
        __m128 zero = _mm_setzero_ps();
        for (int i = 0; i < 8; i++) {
            _mm_store_ps(&phi[i*4], zero);
        }
        this->material = material;
    }
};

#endif
