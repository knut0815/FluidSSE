//
//  Node.h
//  FluidSSE
//
//  Created by Ruby on 1/31/13.
//
//

#ifndef FluidSSE_Node_h
#define FluidSSE_Node_h

struct __attribute__ ((__packed__)) Node {
    float gx[4];
    float u[4];
    float ax[4];
    float m, invM;
    float pad[2];
    void Clear() {
        __m128 zero = _mm_setzero_ps();
        _mm_store_ps(gx, zero);
        _mm_store_ps(u, zero);
        _mm_store_ps(ax, zero);
        m = invM = 0;
    }
};

#endif
