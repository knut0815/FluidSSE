//
//  Material.h
//  FluidSSE
//
//  Created by Ruby on 1/31/13.
//
//

#ifndef FluidSSE_Material_h
#define FluidSSE_Material_h

struct Material {
    float stiffness, restDensity, elasticity, viscosity, smoothing;
    void Initialize(float stiffness, float restDensity, float elasticity, float viscosity, float smoothing) {
        this->stiffness = stiffness;
        this->restDensity = restDensity;
        this->elasticity = elasticity;
        this->viscosity = viscosity;
        this->smoothing = smoothing;
    }
};

#endif
