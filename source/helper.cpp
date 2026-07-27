

#include <helper.h>

/// @brief Builds string that can get printed to the simulation logs
/// @param r current state vector
/// @param t current time stamp
/// @return string of all the data in CSV format
std::string build_telemetry_string(state_vector_t r, double t, int precision) {
    std::stringstream ss;
    ss << std::fixed << std::setprecision(precision);

    ss << t << ',';

    ss << r.x << ',';
    ss << r.y << ',';
    ss << r.z << ',';

    ss << r.vx << ',';
    ss << r.vy << ',';
    ss << r.vz << ',';

    ss << r.wx << ',';
    ss << r.wy << ',';
    ss << r.wz << ',';

    ss << r.qw << ',';
    ss << r.qx << ',';
    ss << r.qy << ',';
    ss << r.qz << ',';

    // Adding latitude, longitude, and altitude from ECEF
    ecef_t ecef = {
        .x = r.x,
        .y = r.y,
        .z = r.z
    };
    lla_t lla;
    ecef_to_lla(ecef, &lla);

    // Then write to the log
    ss << lla.phi * RAD_TO_DEG << ',';
    ss << lla.lambda * RAD_TO_DEG << ',';
    ss << lla.h << ',';

    // And adding temperature, pressure, and air density from altitude
    atmosphere_t atm = get_atmospheric_conditions(lla.h);

    // Write to the string stream
    ss << atm.temp << ',';
    ss << atm.p << ',';
    ss << atm.rho << ',';

    // Adding dynamic pressure to the logs, too
    double q_infinity = 0.5 * (r.vx * r.vx + r.vy * r.vy + r.vz * r.vz) * atm.rho;

    ss << q_infinity;

    return ss.str();
}
