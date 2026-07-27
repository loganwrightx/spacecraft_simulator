/*
Author: Logan Wright
*/

#ifndef SETTINGS_H
#define SETTINGS_H

#include <Eigen/Dense>

#include <motor.h>
#include <motor_mount.h>
#include <pid_controller.h>

extern const double earth_damp_coefficient;
extern Motor motor;
extern MotorMount motor_mount;
extern PIDController x_controller, y_controller;

#endif // SETTINGS_H
