
#include <dynamics.h>

// Solver setting
model_t model_type;

// Kinematic quantities
Eigen::Matrix3d I = Eigen::Matrix3d::Zero(), I_dot = Eigen::Matrix3d::Zero();
Eigen::Vector3d CG = Eigen::Vector3d::Zero(), CG_dot = Eigen::Vector3d::Zero(),
                CP = Eigen::Vector3d::Zero();
double m = 0.0, m_dot = 0.0;

void print_state_vector(state_vector_t r, double t) {
    std::cout << std::fixed << std::setprecision(6);

    std::cout << t << ',';

    std::cout << r.x << ',';
    std::cout << r.y << ',';
    std::cout << r.z << ',';

    std::cout << r.vx << ',';
    std::cout << r.vy << ',';
    std::cout << r.vz << ',';

    std::cout << r.wx << ',';
    std::cout << r.wy << ',';
    std::cout << r.wz << ',';

    std::cout << r.qw << ',';
    std::cout << r.qx << ',';
    std::cout << r.qy << ',';
    std::cout << r.qz << '\n';
}

void print_quaternion(const Eigen::Quaterniond& q) {
    std::cout << "w=" << q.w() << ", x=" << q.x() << ", y=" << q.y() << ", z=" << q.z() << '\n';
}

/// @brief Scales everything except for quaternion values by c
/// @param r current state vector
/// @param c scalar value
/// @return scaled state vector and identity quaternion
state_vector_t scale_state_vector(state_vector_t r, double c) {
    state_vector_t result = {
        .x = c * r.x,
        .y = c * r.y,
        .z = c * r.z,

        .vx = c * r.vx,
        .vy = c * r.vy,
        .vz = c * r.vz,

        .wx = c * r.wx,
        .wy = c * r.wy,
        .wz = c * r.wz,

        .qw = 1.0, // Quaternions don't get scaled for this purpose
        .qx = 0.0, // Quaternions don't get scaled for this purpose
        .qy = 0.0, // Quaternions don't get scaled for this purpose
        .qz = 0.0  // Quaternions don't get scaled for this purpose
    };
    return result;
}

/// @brief Combines state vector r2 onto state vector r1 - element-wise for linear terms, but multiplicative for quaternions
/// @param r1 current true state
/// @param r2 boost/differential state
/// @return boosted final state of r2 applied to r1
state_vector_t combine_state_vectors(state_vector_t r1, state_vector_t r2, double dt_) {
    Eigen::Quaterniond dq;
    // Create omega vector based on r1 only, r1 + r2 omega - euler step
    double w[3] = {r1.wx, r1.wy, r1.wz};
    double wMag = sqrt(w[0] * w[0] + w[1] * w[1] + w[2] * w[2]);
    if (wMag > 1e-12) {
        // DO NOT cut angle in half here, the quaternion constructor from AngleAxis type already handles that for us
        Eigen::AngleAxisd aa(dt_ * wMag, Eigen::Vector3d(w[0] / wMag, w[1] / wMag, w[2] / wMag));
        dq = Eigen::Quaterniond(aa);
    } else {
        dq = Eigen::Quaterniond::Identity();
    }

    // dq applies to the right since in body frame
    Eigen::Quaterniond qfinal = Eigen::Quaterniond(r1.qw, r1.qx, r1.qy, r1.qz) * dq;

    state_vector_t result = {
        .x = r1.x + r2.x,
        .y = r1.y + r2.y,
        .z = r1.z + r2.z,

        .vx = r1.vx + r2.vx,
        .vy = r1.vy + r2.vy,
        .vz = r1.vz + r2.vz,

        .wx = r1.wx + r2.wx,
        .wy = r1.wy + r2.wy,
        .wz = r1.wz + r2.wz,

        .qw = qfinal.w(),
        .qx = qfinal.x(),
        .qy = qfinal.y(),
        .qz = qfinal.z()
    };
    return result;
}

/// @brief Enforces boundary condition - reflects off of earth's surface
/// @param r current state vector
/// @param damping_factor earth's momentum-conservation factor - domain of [0-1)
void enforce_surface_interaction(state_vector_t* r, double damping_factor) {
    lla_t lla_tmp;
    ecef_t ecef_tmp = {
        .x = r->x,
        .y = r->y,
        .z = r->z
    };
    ecef_to_lla(ecef_tmp, &lla_tmp);

    // Enforce solid ground conditions with momentum conservation
    if (lla_tmp.h < 0.0) {
        // Zero out angular rates
        r->wx = 0.0;
        r->wy = 0.0;
        r->wz = 0.0;

        // Rebuild the lla, then get surface-level ECEF coordinates
        lla_tmp.h = 0.0;
        lla_to_ecef(lla_tmp, &ecef_tmp);

        // Attains new position at surface-level
        r->x = ecef_tmp.x;
        r->y = ecef_tmp.y;
        r->z = ecef_tmp.z;

        /* === SURFACE DEFLECTION === */

        // Put velocity data into Eigen vector
        Eigen::Vector3d v_i(r->vx, r->vy, r->vz);

        // Compute surface normal vector for reflection
        Eigen::Vector3d surface_normal = compute_surface_normal(lla_tmp);

        // Compute parallel projection
        Eigen::Vector3d v_parallel = v_i.dot(surface_normal) * surface_normal;

        // Compute perpendicular part via vector subtraction
        Eigen::Vector3d v_perpendicular = v_i - v_parallel;

        // Reverse and damp the terms that got reflected, then add the
        // unmodified perpendicular part back in
        r->vx = (-v_parallel.x() * damping_factor + v_perpendicular.x());
        r->vy = (-v_parallel.y() * damping_factor + v_perpendicular.y());
        r->vz = (-v_parallel.z() * damping_factor + v_perpendicular.z());

        // If any v terms are really small, just zero those ones out
        r->vx = abs(r->vx) > 5e-2 ? r->vx : 0.0;
        r->vy = abs(r->vy) > 5e-2 ? r->vy : 0.0;
        r->vz = abs(r->vz) > 5e-2 ? r->vz : 0.0;
    }
}

/// @brief Computes normal vector from geodetic latitude and longitude
/// @param lla geodetic coordinates
/// @return vector representing the normal in ECEF coordinates
Eigen::Vector3d compute_surface_normal(lla_t lla) {
    Eigen::Vector3d surface_normal(
        cos(lla.phi) * cos(lla.lambda),
        cos(lla.phi) * sin(lla.lambda),
        sin(lla.phi)
    );
    return surface_normal;
}

/// @brief Computes differential state-update given current state and time
/// @param r current state
/// @param t current time
/// @return differential state-update
state_vector_t Dynamics::r3dof(state_vector_t r, double t) {
    /* Forces only since 3dof */
    Eigen::Vector3d F(0.0, 0.0, 0.0);
    state_vector_t dr = {
        .x = r.vx,
        .y = r.vy,
        .z = r.vz,
        .vx = F.x() / m,
        .vy = F.y() / m,
        .vz = F.z() / m,
        .wx = 0.0,
        .wy = 0.0,
        .wz = 0.0,
        .qw = 1.0,
        .qx = 0.0,
        .qy = 0.0,
        .qz = 0.0
    };
    return dr;
}

/// @brief Computes differential state-update given current state and time
/// @param r current state
/// @param t current time
/// @return differential state-update
state_vector_t Dynamics::r6dof(state_vector_t r, double t) {
    /* Compute time-dependent mass properties */
    double m_t;
    Eigen::Matrix3d I_t;
    Eigen::Vector3d CG_t;
    if (!motor_mount.IsBurning()) {
        m_t = m;
        I_t = I;
        CG_t = CG;
    } else if (t <= motor_mount.GetBurnTime()) {
        m_t = m + m_dot * t;
        I_t = I + I_dot * t;
        CG_t = CG + CG_dot * t;
    } else {
        m_t = m + m_dot * motor_mount.GetBurnTime();
        I_t = I + I_dot * motor_mount.GetBurnTime();
        CG_t = CG + CG_dot * motor_mount.GetBurnTime();
    }

    /* Other preliminary quantities */
    Eigen::Matrix3d R = Eigen::Quaterniond(
        r.qw,
        r.qx,
        r.qy,
        r.qz
    ).toRotationMatrix(); // rotation from body vector -> ecef vector

    ecef_t ecef = {
        .x = r.x,
        .y = r.y,
        .z = r.z
    };
    lla_t lla;
    ecef_to_lla(ecef, &lla);

    // Get atmospheric conditions
    atmosphere_t atm = get_atmospheric_conditions(lla.h);

    Eigen::Vector3d moment_arm_b = motor_mount.GetPositionVector();

    /*--------------------------------------------------------*/
    /*------------------ Forces first - ECEF -----------------*/
    /*--------------------------------------------------------*/
    Eigen::Vector3d F(0.0, 0.0, 0.0);

    // Gravity and centrifugal acceleration for rotating reference frame components
    vector3_t acc_gravity_and_centrifugal = effective_gravity(lla, ecef);

    // Then also add in the coriolis effect for full rotating reference frame dynamics
    Eigen::Vector3d F_gravity_centrifugal_coriolis(
        (acc_gravity_and_centrifugal.x + 2 * omega * r.vy) * m_t,
        (acc_gravity_and_centrifugal.y - 2 * omega * r.vx) * m_t,
        acc_gravity_and_centrifugal.z * m_t
    );

    // Then quadratic drag force opposes motion
    double v = sqrt(r.vx * r.vx + r.vy * r.vy + r.vz * r.vz);
    double reference_area = 0.2, cd = 0.5;
    Eigen::Vector3d F_drag(
        -0.5 * atm.rho * cd * reference_area * v * r.vx,
        -0.5 * atm.rho * cd * reference_area * v * r.vy,
        -0.5 * atm.rho * cd * reference_area * v * r.vz
    );

    // And thrust force is computed in the body frame
    Eigen::Vector3d F_thrust_b = motor_mount.GetThrustVector(t);

    // Finally, add all the forces in world-frame
    F += F_gravity_centrifugal_coriolis + F_drag + R * F_thrust_b;

    /* Compute linear acceleration */
    Eigen::Vector3d velocityWorld(r.vx, r.vy, r.vz);
    Eigen::Vector3d accWorld = F / m_t;

    /*------------------------------------------------------------*/
    /*--------------------- Then moments - BODY ------------------*/
    /*------------------------------------------------------------*/
    Eigen::Vector3d M(0.0, 0.0, 0.0);

    // Thrust-vectoring moment
    Eigen::Vector3d M_thrust_vector = moment_arm_b.cross(F_thrust_b);

    // Drag force moment
    Eigen::Vector3d M_drag = (CP - CG_t).cross(R.transpose() * F_drag);

    // Finally, accumulate all the moments in the body-frame
    M += M_thrust_vector + M_drag;

    /* Compute angular acceleration */
    Eigen::Vector3d omegaBody(r.wx, r.wy, r.wz);
    Eigen::Vector3d alphaBody;
    if (t <= motor_mount.GetBurnTime()) {
        alphaBody = I_t.inverse() * (M - I_dot * omegaBody - omegaBody.cross(I_t * omegaBody));
    } else {
        alphaBody = I_t.inverse() * (M - omegaBody.cross(I_t * omegaBody));
    }

    /* Compute quaternion derivative = [0.5 (0 | w) * q] */
    Eigen::Quaterniond dq = Eigen::Quaterniond(r.qw, r.qx, r.qy, r.qz) * Eigen::Quaterniond(0.0, 0.5 * r.wx, 0.5 * r.wy, 0.5 * r.wz);

    // Build the full state-vector time derivative
    state_vector_t dr = {
        .x = r.vx,
        .y = r.vy,
        .z = r.vz,
        .vx = accWorld.x(),
        .vy = accWorld.y(),
        .vz = accWorld.z(),
        .wx = alphaBody.x(),
        .wy = alphaBody.y(),
        .wz = alphaBody.z(),
        .qw = dq.w(),
        .qx = dq.x(),
        .qy = dq.y(),
        .qz = dq.z()
    };

    return dr;
}
