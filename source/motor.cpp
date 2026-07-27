

#include <motor.h>

Motor::Motor(double thrust, double burn_time, double lag) {
    mThrust = thrust;
    mBurnTime = burn_time;
    mLag = lag;
}

Motor::~Motor() {}

void Motor::Ignite(double t) {
    mIgnitionTime = t;
}

double Motor::GetThrust(double t) {
    if (mLag <= t - mIgnitionTime && t - mIgnitionTime - mLag <= mBurnTime && mIgnitionTime >= 0.0) {
        if (!mIsBurning) mIsBurning = true;
        return mThrust;
    }
    if (mIsBurning) mIsBurning = false;
    return 0.0;
}

double Motor::GetBurnTime() {
    return mBurnTime;
}

bool Motor::IsBurning() {
    return mIsBurning;
}
