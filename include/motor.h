

#pragma once

class Motor {
public:
    Motor() = default;
    Motor(double thrust, double burn_time, double lag);
    ~Motor();

    void Ignite(double t);

    double GetThrust(double t);

    double GetBurnTime();

    bool IsBurning();

private:
    double mThrust = 0.0, mBurnTime = 0.0, mLag = 0.0, mIgnitionTime = -1.0;
    bool mIsBurning = false;
};
