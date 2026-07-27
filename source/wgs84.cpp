

#include <wgs84.h>

// Geometric constants of earth WGS84 reference ellipsoid
const double a               = 6378137.0;                   // m
const double b               = 6356752.314245;              // m
const double e_squared       = 6.69437999014e-3;            // unitless
const double e_prime_squared = e_squared * a * a / (b * b); // unitless
const double f               = 1.0 / 298.257223563;         // unitless
const double GM              = 3.986004418e14;              // m^3 / s^2
const double omega           = 72.92115e-6;                 // rad / s

// Harmonic J terms for gravity
const double J[EGM84_J_TERMS] = {
    0.0,                                /* n = 0 (unused) */
    0.0,                                /* n = 1 (unused) */
    1.082629821314e-3,                  /* J2  */
    0.0,                                /* J3  */
   -2.370911200534e-6,                  /* J4  */
    0.0,                                /* J5  */
    6.0832634129e-9,                    /* J6  */
    0.0,                                /* J7  */
   -1.4268105466e-11                    /* J8  */
};

double N(double phi) {
    return a * a / sqrt(a * a * cos(phi) * cos(phi) + b * b * sin(phi) * sin(phi));
}

void lla_to_ecef(lla_t lla, ecef_t* ecef) {
    const double n = N(lla.phi);

    ecef->x = (n + lla.h) * cos(lla.phi) * cos(lla.lambda);
    ecef->y = (n + lla.h) * cos(lla.phi) * sin(lla.lambda);
    ecef->z = (b * b / (a * a) * n + lla.h) * sin(lla.phi);
}

void ecef_to_lla(ecef_t ecef, lla_t* lla) {
    // Application of ferrari's solution
    const double p = sqrt(ecef.x * ecef.x + ecef.y * ecef.y);
    const double F = 54 * b * b * ecef.z * ecef.z;
    const double G = p * p + (1 - e_squared) * ecef.z * ecef.z - e_squared * (a * a - b * b);
    const double c = e_squared * e_squared * F * p * p / (G * G * G);
    const double s = cbrt(1 + c + sqrt(c * c + 2 * c));
    const double k = s + 1 + 1.0 / s;
    const double P = F / (3 * k * k * G * G);
    const double Q = sqrt(1 + 2 * e_squared * e_squared * P);
    const double r0 = -P * e_squared * p / (1 + Q) + sqrt(a * a / 2.0 * (1 + 1.0 / Q) - P * (1 - e_squared) * ecef.z * ecef.z / (Q * (1 + Q)) - P * p * p / 2.0);
    const double U = sqrt((p - e_squared * r0) * (p - e_squared * r0) + ecef.z * ecef.z);
    const double V = sqrt((p - e_squared * r0) * (p - e_squared * r0) + (1 - e_squared) * ecef.z * ecef.z);
    const double z0 = b * b * ecef.z / (a * V);
    lla->h = U * (1 - b * b / (a * V));
    lla->phi = atan((ecef.z + e_prime_squared * z0) / p);
    lla->lambda = atan2(ecef.y, ecef.x);
}

double ecef_radius(ecef_t ecef) {
    return sqrt(ecef.x * ecef.x + ecef.y * ecef.y + ecef.z * ecef.z);
}

vector3_t gravity_zonal(ecef_t ecef)
{
    // Initialize all parameters for the series expansion
    const double r2 = ecef.x * ecef.x + ecef.y * ecef.y + ecef.z * ecef.z;
    const double r = sqrt(r2);
    const double s = ecef.z / r;
    const double u = a / r;

    double radial_factor = 1.0;
    double z_factor = 0.0;

    // Setup the legendre polynomial terms for recursion
    double Pn_2 = 1.0, Pn_1 = s;
    double dPn_2 = 0.0, dPn_1 = 1.0;
    double un = u;

    // Iterate over the desired J terms of EGM84 and use recursion to accumulate
    // the terms
    for (int n = 2; n < EGM84_J_TERMS; ++n) {
        const double Pn  = ((2.0 * n - 1.0) / n) * s * Pn_1 - ((n - 1.0) / n) * Pn_2;
        const double dPn = ((2.0 * n - 1.0) / n) * (Pn_1 + s * dPn_1) - ((n - 1.0) / n) * dPn_2;

        const double term = J[n] * un;
        radial_factor += term * ((n + 1.0) * Pn + s * dPn);
        z_factor += term * dPn;

        Pn_2 = Pn_1;
        Pn_1 = Pn;

        dPn_2 = dPn_1;
        dPn_1 = dPn;

        un *= u;
    }

    // Compute gravitational components
    const double mu_r2 = GM / r2;
    const double gr = -mu_r2 * radial_factor;
    const double gz_ex = mu_r2 * z_factor;

    // Then load into a vector
    vector3_t g = {
        .x = gr * (ecef.x / r),
        .y = gr * (ecef.y / r),
        .z = gr * (ecef.z / r) + gz_ex
    };

    // Add the centrifugal terms
    const double w2 = omega * omega;
    g.x += w2 * ecef.x;
    g.y += w2 * ecef.y;

    return g;
}

vector3_t gravity_somigliana(lla_t lla)
{
    // Save pre-computed trigonometric terms
    const double sin_phi = sin(lla.phi);
    const double cos_phi = cos(lla.phi);
    const double sin2 = sin_phi * sin_phi;

    // Somigliana of the ellipsoid
    const double ge = 9.7803253359;
    const double k = 0.00193185265241;
    const double m = omega * omega * a * a * b / GM;

    const double gamma0 = ge * (1.0 + k * sin2) / sqrt(1.0 - e_squared * sin2);

    // Accumulate standard second-order height continuation
    const double h = lla.h;
    const double gamma = gamma0 * (1.0 - (2.0 / a) * (1.0 + f + m - 2.0 * f * sin2) * h + (3.0 / (a * a)) * h * h);

    // Set direction = outward ellipsoidal normal
    const double nx = cos_phi * cos(lla.lambda);
    const double ny = cos_phi * sin(lla.lambda);
    const double nz = sin_phi;

    // Then build final low-altitude gravity vector
    vector3_t g = {
        .x = -gamma * nx,
        .y = -gamma * ny,
        .z = -gamma * nz
    };

    return g;
}

/// @brief Smooths out the transition between piece-wise functions
/// @param x transition depth
/// @return the hermite transition value for the given transition depth
double smoothstep(double x)
{
    if (x <= 0.0) return 0.0;
    if (x >= 1.0) return 1.0;
    return x * x * (3.0 - 2.0 * x);
}

/// @brief Gets effective gravity which already includes the effect of centrifugal acceleration
/// @param lla lat, lon, alt position
/// @param ecef ecef position
/// @return effective gravity vector
vector3_t effective_gravity(lla_t lla, ecef_t ecef)
{
    // Low altitude region is 0->50km
    const double h0 = 50000.0;
    // And high altitude region is 150km+
    const double h1 = 150000.0;

    // Read ellipsoidal altitude once
    const double h = lla.h;

    if (h <= h0) {
        // pure low-altitude model
        return gravity_somigliana(lla);
    } else if (h >= h1) {
        // pure high-altitude model
        return gravity_zonal(ecef);
    }

    // Blend zone - interpolation of the two acceleration vectors
    const double t = (h - h0) / (h1 - h0);
    const double s = smoothstep(t);

    const vector3_t g_low  = gravity_somigliana(lla);
    const vector3_t g_high = gravity_zonal(ecef);

    vector3_t g = {
        .x = (1.0 - s) * g_low.x + s * g_high.x,
        .y = (1.0 - s) * g_low.y + s * g_high.y,
        .z = (1.0 - s) * g_low.z + s * g_high.z
    };

    return g;
}

void print_lla(lla_t lla) {
    printf("Latitude:  %.6f deg\n", lla.phi * RAD_TO_DEG);
    printf("Longitude: %.6f deg\n", lla.lambda * RAD_TO_DEG);
    printf("Altitude:  %.3f m\n", lla.h);
}

void print_ecef(ecef_t ecef) {
    printf("X: %'.2f m\n", ecef.x);
    printf("Y: %'.2f m\n", ecef.y);
    printf("Z: %'.2f m\n", ecef.z);
}

void print_vector3(vector3_t vector) {
    printf("|v|: %'.2f\n", sqrt(vector.x * vector.x + vector.y * vector.y + vector.z * vector.z));
    printf("vx: %'.2f\n", vector.x);
    printf("vy: %'.2f\n", vector.y);
    printf("vz: %'.2f\n", vector.z);
}

void print_wgs84_model_values() {
    printf("WGS84 Model:\n");
    printf("\ta         = %'.1f m\n", a);
    printf("\tb         = %'.6f m\n", b);
    printf("\tf         = %.6e\n", f);
    printf("\te_squared = %.11e\n", e_squared);
    printf("\tGM        = %.9e m^3/s^2\n", GM);
    printf("\tomega     = %.5e rad/s\n", omega);
}
