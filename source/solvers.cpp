

#include <solvers.h>

double dt = 0.0;

state_vector_t Solvers::rk4(state_vector_t r, double t) {
    state_vector_t k1 = Dynamics::r6dof(r, t);

    state_vector_t k1_scaled = scale_state_vector(k1, 0.5 * dt);
    state_vector_t r_plus_k1 = combine_state_vectors(r, k1_scaled, 0.5 * dt);
    state_vector_t k2 = Dynamics::r6dof(r_plus_k1, t + 0.5 * dt);

    state_vector_t k2_scaled = scale_state_vector(k2, 0.5 * dt);
    state_vector_t r_plus_k2 = combine_state_vectors(r, k2_scaled, 0.5 * dt);
    state_vector_t k3 = Dynamics::r6dof(r_plus_k2, t + 0.5 * dt);

    state_vector_t k3_scaled = scale_state_vector(k3, dt);
    state_vector_t r_plus_k3 = combine_state_vectors(r, k3_scaled, dt);
    state_vector_t k4 = Dynamics::r6dof(r_plus_k3, t + dt);

    // Get omega updates using rk4
    double dwx = dt * (k1.wx + 2.0 * k2.wx + 2.0 * k3.wx + k4.wx) / 6.0;
    double dwy = dt * (k1.wy + 2.0 * k2.wy + 2.0 * k3.wy + k4.wy) / 6.0;
    double dwz = dt * (k1.wz + 2.0 * k2.wz + 2.0 * k3.wz + k4.wz) / 6.0;

    // Compute average of current step and forward rk4 step
    double wx_avg = r.wx + 0.5 * dwx;
    double wy_avg = r.wy + 0.5 * dwy;
    double wz_avg = r.wz + 0.5 * dwz;

    // Update attitude using averaged step
    double wMag = sqrt(wx_avg * wx_avg + wy_avg * wy_avg + wz_avg * wz_avg);
    Eigen::Quaterniond dq;
    if (wMag > 1e-12) {
        // DO NOT cut angle in half here, the quaternion constructor from AngleAxis type already handles that for us
        Eigen::AngleAxisd aa(dt * wMag, Eigen::Vector3d(wx_avg / wMag, wy_avg / wMag, wz_avg / wMag));
        dq = Eigen::Quaterniond(aa);
    } else {
        dq = Eigen::Quaterniond::Identity();
    }
    // dq applies to right since in body frame
    Eigen::Quaterniond qfinal = Eigen::Quaterniond(r.qw, r.qx, r.qy, r.qz) * dq;

    // Ensure normalized qfinal
    qfinal.normalize();

    state_vector_t result = {
        .x = r.x + dt * (k1.x + 2.0 * k2.x + 2.0 * k3.x + k4.x) / 6.0,
        .y = r.y + dt * (k1.y + 2.0 * k2.y + 2.0 * k3.y + k4.y) / 6.0,
        .z = r.z + dt * (k1.z + 2.0 * k2.z + 2.0 * k3.z + k4.z) / 6.0,

        .vx = r.vx + dt * (k1.vx + 2.0 * k2.vx + 2.0 * k3.vx + k4.vx) / 6.0,
        .vy = r.vy + dt * (k1.vy + 2.0 * k2.vy + 2.0 * k3.vy + k4.vy) / 6.0,
        .vz = r.vz + dt * (k1.vz + 2.0 * k2.vz + 2.0 * k3.vz + k4.vz) / 6.0,

        /*
         This will be slightly wrong because each k_n bases off of a different
         frame, so omega and alpha come from different coordinate systems.
         However, for small rotational rates, it's a decent approximation
        */
        .wx = r.wx + dwx,
        .wy = r.wy + dwy,
        .wz = r.wz + dwz,

        .qw = qfinal.w(),
        .qx = qfinal.x(),
        .qy = qfinal.y(),
        .qz = qfinal.z()
    };

    return result;
}
