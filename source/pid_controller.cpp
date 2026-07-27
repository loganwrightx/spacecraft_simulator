
#include <pid_controller.h>

PIDController::PIDController(double kp, double ki, double kd, double min, double max, double target)
    : pKp(kp), pKi(ki), pKd(kd), pMin(min), pMax(max), pTarget(target)
{
    //
}

void PIDController::Init(double initial_state, double timestamp) {
    pLastTime = timestamp;
    pE = pTarget - initial_state;
}

void PIDController::AdjustTarget(double new_target) {
    // Update the error term first
    pE += new_target - pTarget;

    // Then set new target
    pTarget = new_target;
}

void PIDController::ResetWindup() {
    // Sets integral term to zero
    pIe = 0.0;
}

void PIDController::Step(double new_state, double timestamp) {
    double dt = timestamp - pLastTime;
    double old_state = pTarget - pE;

    // Update integral of error using a simple Euler step
    pIe += pE * dt;

    // Update derivative of error using linear slope form
    if (dt > 0.0) {
        pDe = (new_state - old_state) / dt;
    }

    // Update error term by differencing new_state from the target
    pE = pTarget - new_state;

    // Then update stateful parameters
    pLastTime = timestamp;
}

double PIDController::GetSample() {
    double response = pKp * pE + pKi * pIe + pKd * pDe;

    if (response < pMin) {
        return pMin;
    } else if (response > pMax) {
        return pMax;
    } else {
        return response;
    }
}
