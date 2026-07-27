

#include <motor_mount.h>

MotorMount::MotorMount(double center_x, double center_y, double max_rate, double latency, Eigen::Vector3d position, Motor& motor)
    : mX(center_x, max_rate, latency), mY(center_y, max_rate, latency), sPosition(position), mMotor(motor)
{
    //
}

MotorMount::~MotorMount() {}

Eigen::Vector3d MotorMount::GetPositionVector() {
    return sPosition;
}

void MotorMount::Ignite(double t) {
    mMotor.Ignite(t);
}

/// @brief Commands thrust vector mount to new angles at time t
/// @param x_angle new x angle in degrees
/// @param y_angle new y angle in degrees
/// @param t current timestamp in seconds
void MotorMount::Command(double x_angle, double y_angle, double t) {
    mX.SetTo(x_angle, t);
    mY.SetTo(y_angle, t);
}

void MotorMount::Step(double dt) {
    mX.Step(dt);
    mY.Step(dt);
}

double MotorMount::GetBurnTime() {
    return mMotor.GetBurnTime();
}

bool MotorMount::IsBurning() {
    return mMotor.IsBurning();
}

Eigen::Vector3d MotorMount::GetThrustVector(double t) {
    double angle_x = mX.GetActuation();
    double angle_y = mY.GetActuation();
    double thrust = mMotor.GetThrust(t);

    Eigen::Vector3d thrust_vector = thrust * Eigen::Vector3d(
        cos(angle_x * M_PI / 180.0) * cos(angle_y * M_PI / 180.0), // z -> x (body-up)
        sin(angle_x * M_PI / 180.0) * cos(angle_y * M_PI / 180.0), // x -> y (body-sideways)
        sin(angle_y * M_PI / 180.0) // y -> z (body-other-sideways)
    );

    return thrust_vector;
}
