/*
Author: Logan Wright
Description: This is the main entry point of the spacecraft_simulator program.
    To configure the desired settings, take a look at `shared.cpp` for the
    declarations of the engine, motor mount, and PID controller objects - update
    them to your needs. Within `main.cpp`, you'll want to also modify the custom
    parameters that get setup in the main() function below. Follow comments
    accordingly to understand how to get your desired behavior.
*/

#include <helper.h>
#include <dynamics.h>
#include <logger.h>

// Precision of data in the simulation log file
#define LOG_PRECISION 8

// Loggers to support plotting and debugging
Logger simulationLog("sim.csv");
Logger debugLog("debug.txt");

int main(void) {
    /* UNDERSTANDING TIME-DEPENDENT MASS
     *  The mass of a rocket is not constant, and it's not even approximately
     *  constant. Set mass and m_dot to be the initial mass and rate of change
     *  of mass per unit of burn time, respectively.
     *
     * EDIT TO YOUR DESIGN PARAMETERS
     */
    m = 0.7;                                     /* kg */
    m_dot = -0.0369 / motor_mount.GetBurnTime(); /* kg/s */

    /* UNDERSTANDING TIME-DEPENDENT MASS MOMENT OF INERTIA TENSOR
     *  This usually depends linearly on the mass of the rocket, but sometimes
     *  also depends on complex geometries. This is a constant initial property
     *  in the body-frame of the rocket, so it is typically diagonalized. Feel
     *  free to set the other 6 terms of the tensor explicitly if they're
     *  non-zero. The I_dot term is the linearized time-derivative of this
     *  quantity - set that as needed as well.
     *
     * EDIT THESE FIELDS
     */
    I(0, 0) = 5. / 8. * m * 0.074 * 0.074; /* kg m^2 */
    I(1, 1) = 1. / 12. * m;                /* kg m^2 */
    I(2, 2) = 1. / 12. * m;                /* kg m^2 */

    I_dot(0, 0) = 5. / 8. * m_dot * 0.074 * 0.074; /* kg m^2/s */
    I_dot(1, 1) = 1. / 12. * m_dot;                /* kg m^2/s */
    I_dot(2, 2) = 1. / 12. * m_dot;                /* kg m^2/s */

    /* UNDERSTANDING CENTER OF GRAVITY
     *  The center of gravity (also called center of mass) is the location where
     *  all non-aerodynamic forces are applied to a system (the rocket). Refer
     *  to the coordinate-system conventions where body-up is +x, and +y/+z
     *  point outwards such that they form a right-hand coordinate system. The
     *  CG_dot term is the linear time derivative of the center of gravity.
     *
     * EDIT AS NEEDED
     */
    CG(0) = 0.0; /* m */
    CG(1) = 0.0; /* m */
    CG(2) = 0.0; /* m */

    CG_dot(0) = 0.1 / motor_mount.GetBurnTime(); /* m/s */
    CG_dot(0) = 0.0;                             /* m/s */
    CG_dot(0) = 0.0;                             /* m/s */

    /* UNDERSTANDING CENTER OF PRESSURE
     *  The center of pressure (in body-coordinates) is the average position of
     *  all surface area of the vehicle. The aerodynamic forces all interact at
     *  that precise position. Refer to the coordinate-system conventions where
     *  body-up is +x, and +y/+z point outwards such that they form a right-hand
     *  coordinate system.
     *
     * EDIT AS NEEDED
     */
    CP(0) = -0.2;
    CP(1) = 0.0;
    CP(2) = 0.0;

    // Initialize position on earth with LLA coordinates - readable GPS position
    double phi = 35.0000, lambda = 0.0000, altitude = 0.0;

    // Convert lat and lon to rads
    phi *= DEG_TO_RAD;
    lambda *= DEG_TO_RAD;

    /* UNDERSTANDING LLA AND ECEF COORDINATE SYSTEMS
     *  WGS84 supports multiple coordinate systems, but GPS most often provides
     *  latitude, longitude, and altitude for geolocation. However, ECEF is a
     *  cartesian coordinate system that makes physics much simpler and
     *  computationally cheaper. You'll see the code alternate between these
     *  coordinate systems depending on the use-case for position information.
     *
     * DO NOT EDIT
     */
    lla_t lla = {
        .phi = phi,
        .lambda = lambda,
        .h = altitude
    };
    ecef_t ecef;
    lla_to_ecef(lla, &ecef);

    // Also initialize ECEF velocity relative to earth's surface
    Eigen::Vector3d initial_velocity(
        0.0,
        0.0,
        0.0
    );

    /* UNDERSTANDING INITIAL LOCAL ORIENTATION ON EARTH'S SURFACE
     *  In WGS84, we use geodetic coordinates (latitude and longitude). They are
     *  convenient for computing local normal vectors because they can be used
     *  like spherical angular coordinates. The following is always true to
     *  represent the rotation of a vector in local coordinates to the ECEF
     *  frame.
     *
     * DO NOT EDIT
     */
    Eigen::Quaterniond earth_normal_attitude(
        cos(lla.phi / 2.) * cos(lla.lambda / 2.),
        sin(lla.phi / 2.) * sin(lla.lambda / 2.),
        -cos(lla.lambda / 2.) * sin(lla.phi / 2.),
        cos(lla.phi / 2.) * sin(lla.lambda / 2.)
    );

    /* UNDERSTANDING THE LOCAL REFERENCE FRAME CONVENTION
     *  The operator, by default, is assumed to be looking at the rocket with
     *  the rocket body +x-axis pointed straight up into the sky, rocket +y-axis
     *  pointing to the operators left, and rocket body +z-axis points towards
     *  the operator - this builds a local right-handed reference frame. For
     *  any lat/lon coordinate that is non-zero, this mapping is true:
     *   - +x = UP
     *   - +y = EAST
     *   - +z = NORTH
     *
     * EDIT ACCORDING TO RELATIVE VEHICLE ORIENTATION IN LOCAL FRAME
     */
    Eigen::AngleAxisd vehicle_aa(
        0.0,                           /* ANGLE TO ROTATE BY IN RADIANS */
        Eigen::Vector3d(0.0, 0.0, 1.0) /* LOCAL AXIS TO ROTATE AROUND   */
    );
    Eigen::Quaterniond relative_vehicle_attitude(vehicle_aa);

    /* FULL INITIAL ATTITUDE
     *  The following is the algebraic operation required to combine rotations.
     *
     * DO NOT EDIT
     */
    Eigen::Quaterniond qi = relative_vehicle_attitude * earth_normal_attitude;

    // Initial angular velocity vector in rocket's body frame
    Eigen::Vector3d omega_b(
        0.0,
        0.0,
        0.0
    );

    /* UNDERSTANDING STATE VECTORS AND THEIR ROLE
     *  The state vector data structure holds all the kinematic information
     *  needed to track the rigid body (the rocket). Once we have all the
     *  initial conditions of the rocket, we load them all into this data
     *  structure. All numerical solvers here interface with this data type.
     *
     * DO NOT EDIT
     */
    state_vector_t r = {
        // Ecef position from lla
        .x = ecef.x,
        .y = ecef.y,
        .z = ecef.z,

        // Ecef-velocity from ecef
        .vx = initial_velocity.x(),
        .vy = initial_velocity.y(),
        .vz = initial_velocity.z(),

        // Body-attached frame
        .wx = omega_b.x(),
        .wy = omega_b.y(),
        .wz = omega_b.z(),

        // Start from initial rotation which is body->ECEF
        .qw = qi.w(),
        .qx = qi.x(),
        .qy = qi.y(),
        .qz = qi.z()
    };

    /* UNDERSTANDING TIME SETTINGS
     *  Here, we set t = 0.0 initially but this should technically be a UTC
     *  true timestamp. Relative to t, you'll also need to choose a tf (final
     *  timestamp) for when to end the simulation. And finally, you'll need to
     *  set the time step-size - the smaller, the more accurate but it becomes
     *  computationally expensive quickly.
     *
     * EDIT THESE FIELDS AS NEEDED
     */
    double t = 0.0;
    const double tf = t + 20.0;
    dt = 1e-3;

#ifdef USING_TVC
    const int update_rate = 100; // Hz
    const double tvc_dt = 1. / ((double)update_rate);
    double last_tvc_update = t;

    x_controller.Init(r.qz, t);
    y_controller.Init(r.qy, t);
#endif // USING_TVC

    motor_mount.Ignite(t);

    // Write header to telemetry log
    simulationLog.WriteLine("t,x,y,z,vx,vy,vz,wx,wy,wz,qw,qx,qy,qz,lat,lon,alt,T,p,rho,q_infinity");

    while (t <= tf) {
        // Write state to simulation log file
        simulationLog.WriteLine(build_telemetry_string(r, t, LOG_PRECISION));

        // Compute next state using RK4
        r = Solvers::rk4(r, t);

        // Enforce boundary condition with earth surface
        enforce_surface_interaction(&r, earth_damp_coefficient);

#ifdef USING_TVC
        if (t - last_tvc_update >= tvc_dt) {
            // Reset the last_tvc_update timestamp
            last_tvc_update = t;

            // Update the PID controllers
            x_controller.Step(r.qz, t);
            y_controller.Step(r.qy, t);

            // Command the thrust vector mount to a new attitude
            motor_mount.Command(
                x_controller.GetSample(),
                y_controller.GetSample(),
                t
            );
        }
#endif // USING_TVC

        // Step the motor mount for time progression
        motor_mount.Step(dt);

        t += dt;
    }

    // Write last state to simulation log file
    simulationLog.WriteLine(build_telemetry_string(r, t, LOG_PRECISION));

    return 0;
}
