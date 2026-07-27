#ifndef INTERPOLATION_H
#define INTERPOLATION_H

#include <vector>

class InterpolatorBase {
public:
    InterpolatorBase() = default;
    ~InterpolatorBase() {};

    virtual void Load();
};

#endif // INTERPOLATION_H
