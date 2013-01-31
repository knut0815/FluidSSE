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
    int c;
    int pad[3];
    void Initialize(float px, float py, float pz) {
        x[0] = px;
        x[1] = py;
        x[2] = pz;
        x[3] = 0;
        __m128 zero = _mm_setzero_ps();
        _mm_store_ps(u, zero);
        //u[0] = .01;
        for (int i = 0; i < 8; i++) {
            _mm_store_ps(&phi[i*4], zero);
        }
    }
};

#endif
