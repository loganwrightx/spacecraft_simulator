
#pragma once

#include <Eigen/Dense>

#include <servo.h>
#include <motor.h>

class MotorMount {
public:
    MotorMount(double center_x, double center_y, double max_rate, double latency, Eigen::Vector3d position, Motor& motor);
    ~MotorMount();

    Eigen::Vector3d GetPositionVector();

    void Ignite(double t);

    void Command(double x_angle, double y_angle, double t);

    void Step(double dt);

    double GetBurnTime();

    bool IsBurning();

    Eigen::Vector3d GetThrustVector(double t);

private:
    Servo mX, mY;
    Eigen::Vector3d sPosition;
    Motor& mMotor;
};
