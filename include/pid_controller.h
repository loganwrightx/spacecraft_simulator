

#ifndef PID_CONTROLLER_H
#define PID_CONTROLLER_H

class PIDController {
public:
    PIDController(double kp, double ki, double kd, double min, double max, double target);
    ~PIDController() {}

    void Init(double initial_state, double timestamp);
    void AdjustTarget(double new_target);
    void ResetWindup();
    void Step(double new_state, double timestamp);
    double GetSample();

private:
    double pKp, pKi, pKd, pMin, pMax;
    double pTarget;
    double pIe = 0.0, pE = 0.0, pDe = 0.0, pLastTime = 0.0;
};

#endif // PID_CONTROLLER_H
