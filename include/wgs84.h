/*
Author: Logan Wright
Description: My own implementation of the WGS84 and EGM84 models (earth nav)
*/

#ifndef WGS84_H
#define WGS84_H

#include <math.h>
#include <stdio.h>
#include <locale.h>

#define RAD_TO_DEG 180.0 / M_PI
#define DEG_TO_RAD M_PI / 180.0

#define EGM84_J_TERMS 9

// Geometric constants of earth WGS84 reference ellipsoid
extern const double a;
extern const double b;
extern const double e_squared;
extern const double e_prime_squared;
extern const double f;
extern const double GM;
extern const double omega;

typedef struct {
    double x;
    double y;
    double z;
} ecef_t;

/// @brief Latitude and longitude (radians), and altitude (meters)
typedef struct {
    double phi;
    double lambda;
    double h;
} lla_t;

/// @brief Ecef position (meters)
typedef struct {
    double x;
    double y;
    double z;
} vector3_t;

double N(double phi);

void lla_to_ecef(lla_t lla, ecef_t* ecef);

void ecef_to_lla(ecef_t ecef, lla_t* lla);

double ecef_radius(ecef_t ecef);

vector3_t gravity_zonal(ecef_t ecef);

vector3_t gravity_somigliana(lla_t lla);

double smoothstep(double x);

vector3_t effective_gravity(lla_t lla, ecef_t ecef);

void print_lla(lla_t lla);

void print_ecef(ecef_t ecef);

void print_vector3(vector3_t vector);

void print_wgs84_model_values();

#endif // WGS84_H
