#ifndef SOLVERS_H
#define SOLVERS_H

#include <dynamics.h>

extern double dt;

namespace Solvers {

state_vector_t rk4(state_vector_t r, double t);

};

#endif // SOLVERS_H
