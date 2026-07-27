/*
Author: Logan Wright
Description: Rigid body dynamics formulas for a rocket
*/

#pragma once

#include <cmath>
#include <iostream>
#include <iomanip>

#include <Eigen/Dense>

#include <wgs84.h>
#include <atmosphere.h>
#include <motor_mount.h>
#include <settings.h>

typedef enum {
    R3DOF = 0,
    R6DOF
} model_t;

extern model_t model_type;
extern Eigen::Matrix3d I, I_dot;
extern Eigen::Vector3d CG, CG_dot, CP;
extern double m, m_dot;
extern MotorMount motor_mount;

typedef struct {
    double x;
    double y;
    double z;
    double vx;
    double vy;
    double vz;
    double wx;
    double wy;
    double wz;
    double qw;
    double qx;
    double qy;
    double qz;
} state_vector_t;

void print_state_vector(state_vector_t r, double t);
void print_quaternion(const Eigen::Quaterniond& q);
void enforce_surface_interaction(state_vector_t* r, double damping_factor);
Eigen::Vector3d compute_surface_normal(lla_t lla);

state_vector_t scale_state_vector(state_vector_t r, double c);
state_vector_t combine_state_vectors(state_vector_t r1, state_vector_t r2, double dt_);

namespace Dynamics {

state_vector_t r3dof(state_vector_t r, double t);
state_vector_t r6dof(state_vector_t r, double t);

};
