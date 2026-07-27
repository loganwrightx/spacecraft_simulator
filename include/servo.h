/*
Author: Logan Wright
Description: Servo model with latency, error, and rate limiting
*/

#pragma once

class Servo {
public:
    Servo() = default;
    Servo(double initial_angle, double max_rate, double latency);
    ~Servo();

    void SetTo(double new_angle, double t);

    double GetActuation();

    void Step(double dt);

private:
    double sAngle = 0.0, sMaxRate = 0.0, sLatency = 0.0;
};
