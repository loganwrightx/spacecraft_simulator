

#include <servo.h>

Servo::Servo(double center, double max_rate, double latency) {
    sAngle = center;
    sMaxRate = max_rate;
    sLatency = latency;
}

Servo::~Servo() {}

void Servo::SetTo(double new_angle, double t) {
    sAngle = new_angle;
}

double Servo::GetActuation() {
    return sAngle;
}

void Servo::Step(double dt) {
    //
}
