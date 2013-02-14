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
    float stiffness, restDensity, bulkViscosity, elasticity, viscosity, smoothing, yieldPoint, yieldRate;
    void Initialize(float stiffness, float restDensity, float bulkViscosity, float elasticity, float viscosity, float smoothing, float yieldRate) {
        this->stiffness = stiffness*.5/restDensity;
        this->restDensity = restDensity;
        this->bulkViscosity = bulkViscosity;
        this->elasticity = elasticity;
        this->viscosity = viscosity;
        this->smoothing = smoothing;
        this->yieldRate = yieldRate;
    }
};

#endif
