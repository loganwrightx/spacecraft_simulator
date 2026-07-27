#ifndef ATMOSPHERE_H
#define ATMOSPHERE_H

#include <cmath>

typedef struct {
    double rho; // air density
    double temp; // temperature
    double p; // pressure
} atmosphere_t;

atmosphere_t get_atmospheric_conditions(double altitude);

#endif // ATMOSPHERE_H
