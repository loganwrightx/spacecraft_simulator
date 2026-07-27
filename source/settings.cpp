/*
Author: Logan Wright
Description: Set of all global declarations needed across the project.
*/

#include <settings.h>

// Bounciness of earth's surface - must fall between [0, 1)
//  The smaller this number, the less bouncy
const double earth_damp_coefficient = 0.1;

/* UNDERSTANDING THE SIMPLE MOTOR MODEL
 *  This model is very simplified and assumes constant thrust for a defined
 *  period of time. The motor mount then uses the motor model as well as thrust
 *  vector inputs to produce a thrust force vector. Change the thrust and burn
 *  time settings of this motor to match your configuration. Right now, this is
 *  a simplified Estes E12-6 black powder motor model.
 *
 * EDIT THIS AS NEEDED
 */
Motor motor(15.2, 3.5, 0.0);
MotorMount motor_mount(0.0, 0.0, 10.0, 0.0, Eigen::Vector3d(-0.5, 0.0, 0.0), motor);

/* UNDERSTANDING THE PID CONTROLLER
 *  These controllers are meant to set commanded actuations on the thrust vector
 *  mount to steer the vehicle to the desired trajectory. Change the PID gains
 *  to see different aggressiveness in control. Right now, this is configured to
 *  use 2 separate PID controllers for each rotation axis of the gimbal.
 *
 * EDIT AS NEEDED
 */
const double kP = 69.0, kI = 0.294117647059, kD = 0.0735294117647, min_lim = -25.0, max_lim = 25.0;
PIDController x_controller(-kP, -kI, -kD, min_lim, max_lim, 0.0), y_controller(kP, kI, kD, min_lim, max_lim, 0.0);
