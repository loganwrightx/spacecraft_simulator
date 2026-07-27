#include <gtest/gtest.h>

#include <wgs84.h>
#include <cmath>
#include <cstring>

// ---------------------------------------------------------------------------
// Test fixture
// ---------------------------------------------------------------------------
class Wgs84Test : public testing::Test {
protected:
    static constexpr double TOL_POS = 1e-6;   // ~1 µm
    static constexpr double TOL_ANG = 1e-9;   // radian
    static constexpr double TOL_G   = 1e-4;   // m/s²

    // Blend window used by the production model
    static constexpr double H0 =  50000.0;    // m
    static constexpr double H1 = 150000.0;    // m

    Wgs84Test() = default;
    ~Wgs84Test() override = default;

    void SetUp() override {}
    void TearDown() override {}
};

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------
TEST_F(Wgs84Test, TestWGS84Constants) {
    EXPECT_DOUBLE_EQ(a, 6378137.0);
    EXPECT_NEAR(b, 6356752.314245, 1e-6);
    EXPECT_NEAR(f, 1.0 / 298.257223563, 1e-12);
    EXPECT_NEAR(e_squared, 6.69437999014e-3, 1e-14);
    EXPECT_DOUBLE_EQ(GM, 3.986004418e14);
    EXPECT_NEAR(omega, 72.92115e-6, 1e-12);
}

// ---------------------------------------------------------------------------
// N(φ)
// ---------------------------------------------------------------------------
TEST_F(Wgs84Test, TestN_Equator) {
    EXPECT_NEAR(N(0.0), a, TOL_POS);
}

TEST_F(Wgs84Test, TestN_Pole) {
    const double expected = (a * a) / b;
    EXPECT_NEAR(N( M_PI_2), expected, TOL_POS);
    EXPECT_NEAR(N(-M_PI_2), expected, TOL_POS);
}

TEST_F(Wgs84Test, TestN_45deg) {
    const double n = N(45.0 * DEG_TO_RAD);
    EXPECT_GT(n, a);
    EXPECT_LT(n, (a * a) / b);
}

// ---------------------------------------------------------------------------
// LLA -> ECEF
// ---------------------------------------------------------------------------
TEST_F(Wgs84Test, TestLlaToEcef_EquatorZeroHeight) {
    lla_t lla = {0.0, 0.0, 0.0};
    ecef_t ecef;
    lla_to_ecef(lla, &ecef);

    EXPECT_NEAR(ecef.x, a, TOL_POS);
    EXPECT_NEAR(ecef.y, 0.0, TOL_POS);
    EXPECT_NEAR(ecef.z, 0.0, TOL_POS);
}

TEST_F(Wgs84Test, TestLlaToEcef_EquatorWithHeight) {
    lla_t lla = {0.0, 0.0, 1000.0};
    ecef_t ecef;
    lla_to_ecef(lla, &ecef);

    EXPECT_NEAR(ecef.x, a + 1000.0, TOL_POS);
    EXPECT_NEAR(ecef.y, 0.0, TOL_POS);
    EXPECT_NEAR(ecef.z, 0.0, TOL_POS);
}

TEST_F(Wgs84Test, TestLlaToEcef_NorthPole) {
    // tiny offset avoids exact singularity in inverse later
    lla_t lla = {M_PI_2 - 1e-12, 0.0, 0.0};
    ecef_t ecef;
    lla_to_ecef(lla, &ecef);

    EXPECT_NEAR(ecef.x, 0.0, 1e-3);
    EXPECT_NEAR(ecef.y, 0.0, 1e-3);
    EXPECT_NEAR(ecef.z, b, 1e-3);
}

TEST_F(Wgs84Test, TestLlaToEcef_Arbitrary) {
    lla_t lla = {40.0 * DEG_TO_RAD, -75.0 * DEG_TO_RAD, 100.0};
    ecef_t ecef;
    lla_to_ecef(lla, &ecef);

    EXPECT_GT(ecef.x, 1.0e6);
    EXPECT_LT(ecef.x, 5.0e6);
    EXPECT_LT(ecef.y, 0.0);
    EXPECT_GT(ecef.z, 3.0e6);
}

// ---------------------------------------------------------------------------
// ECEF -> LLA  (round-trips)
// ---------------------------------------------------------------------------
TEST_F(Wgs84Test, TestEcefToLla_RoundTrip_Equator) {
    lla_t orig = {0.0, 30.0 * DEG_TO_RAD, 250.0};
    ecef_t ecef;
    lla_to_ecef(orig, &ecef);

    lla_t recovered;
    ecef_to_lla(ecef, &recovered);

    EXPECT_NEAR(recovered.phi,    orig.phi,    TOL_ANG);
    EXPECT_NEAR(recovered.lambda, orig.lambda, TOL_ANG);
    EXPECT_NEAR(recovered.h,      orig.h,      TOL_POS);
}

TEST_F(Wgs84Test, TestEcefToLla_RoundTrip_MidLatitude) {
    lla_t orig = {51.477 * DEG_TO_RAD, -0.0015 * DEG_TO_RAD, 46.0};
    ecef_t ecef;
    lla_to_ecef(orig, &ecef);

    lla_t recovered;
    ecef_to_lla(ecef, &recovered);

    EXPECT_NEAR(recovered.phi,    orig.phi,    TOL_ANG);
    EXPECT_NEAR(recovered.lambda, orig.lambda, TOL_ANG);
    EXPECT_NEAR(recovered.h,      orig.h,      TOL_POS);
}

TEST_F(Wgs84Test, TestEcefToLla_RoundTrip_HighAltitude) {
    lla_t orig = {-33.9 * DEG_TO_RAD, 151.2 * DEG_TO_RAD, 400000.0};
    ecef_t ecef;
    lla_to_ecef(orig, &ecef);

    lla_t recovered;
    ecef_to_lla(ecef, &recovered);

    EXPECT_NEAR(recovered.phi,    orig.phi,    TOL_ANG);
    EXPECT_NEAR(recovered.lambda, orig.lambda, TOL_ANG);
    EXPECT_NEAR(recovered.h,      orig.h,      1e-3);
}

// ---------------------------------------------------------------------------
// ecef_radius
// ---------------------------------------------------------------------------
TEST_F(Wgs84Test, TestEcefRadius_UnitVector) {
    ecef_t v = {3.0, 4.0, 12.0};
    EXPECT_NEAR(ecef_radius(v), 13.0, 1e-12);
}

TEST_F(Wgs84Test, TestEcefRadius_Equator) {
    ecef_t ecef = {a, 0.0, 0.0};
    EXPECT_NEAR(ecef_radius(ecef), a, TOL_POS);
}

TEST_F(Wgs84Test, TestEcefRadius_Pole) {
    ecef_t ecef = {0.0, 0.0, b};
    EXPECT_NEAR(ecef_radius(ecef), b, TOL_POS);
}

// ---------------------------------------------------------------------------
// effective_gravity – pure low-altitude regime (Somigliana)
// ---------------------------------------------------------------------------
TEST_F(Wgs84Test, TestEffectiveGravity_EquatorSurface) {
    lla_t  lla  = {0.0, 0.0, 0.0};
    ecef_t ecef = {a, 0.0, 0.0};

    vector3_t g = effective_gravity(lla, ecef);
    const double mag = std::sqrt(g.x*g.x + g.y*g.y + g.z*g.z);

    // Exact WGS 84 normal gravity at equator
    EXPECT_NEAR(mag, 9.7803253359, 5e-4);

    EXPECT_LT(g.x, 0.0);
    EXPECT_NEAR(g.y, 0.0, TOL_G);
    EXPECT_NEAR(g.z, 0.0, TOL_G);
}

TEST_F(Wgs84Test, TestEffectiveGravity_NorthPole) {
    lla_t  lla  = {M_PI_2, 0.0, 0.0};
    ecef_t ecef = {0.0, 0.0, b};

    vector3_t g = effective_gravity(lla, ecef);
    const double mag = std::sqrt(g.x*g.x + g.y*g.y + g.z*g.z);

    // Exact WGS 84 normal gravity at pole
    EXPECT_NEAR(mag, 9.8321849378, 5e-4);

    EXPECT_NEAR(g.x, 0.0, TOL_G);
    EXPECT_NEAR(g.y, 0.0, TOL_G);
    EXPECT_LT(g.z, 0.0);
}

TEST_F(Wgs84Test, TestEffectiveGravity_LowAltitude) {
    // Well inside the pure Somigliana region
    lla_t  lla  = {45.0 * DEG_TO_RAD, 10.0 * DEG_TO_RAD, 20000.0};
    ecef_t ecef;
    lla_to_ecef(lla, &ecef);

    vector3_t g = effective_gravity(lla, ecef);
    const double mag = std::sqrt(g.x*g.x + g.y*g.y + g.z*g.z);

    // Surface value at 45° is ≈ 9.806; free-air reduction ≈ 6.3 mGal/km
    // -> ≈ 9.745 m/s² at 20 km.  Allow a small window.
    EXPECT_NEAR(mag, 9.745, 0.01);
    EXPECT_LT(g.z, 0.0);
}

// ---------------------------------------------------------------------------
// effective_gravity – pure high-altitude regime (zonal SH)
// ---------------------------------------------------------------------------
TEST_F(Wgs84Test, TestEffectiveGravity_HighAltitude) {
    const double r = 2.0 * a;
    lla_t  lla  = {0.0, 0.0, r - a};
    ecef_t ecef = {r, 0.0, 0.0};

    vector3_t g = effective_gravity(lla, ecef);
    const double mag = std::sqrt(g.x*g.x + g.y*g.y + g.z*g.z);

    // Zonal model at r = 2a on equator ≈ 2.3808 m/s²
    // (Newtonian 2.4496 − J₂ correction − centrifugal)
    EXPECT_NEAR(mag, 2.3808, 0.01);

    EXPECT_LT(g.x, 0.0);
    EXPECT_NEAR(g.y, 0.0, 1e-3);
    EXPECT_NEAR(g.z, 0.0, 1e-3);
}

TEST_F(Wgs84Test, TestEffectiveGravity_VeryHighAltitude) {
    // Far above the blend window – pure zonal path
    lla_t  lla  = {0.0, 0.0, 1.0e7};          // 10 000 km altitude
    ecef_t ecef;
    lla_to_ecef(lla, &ecef);

    vector3_t g = effective_gravity(lla, ecef);
    const double r   = ecef_radius(ecef);
    const double mag = std::sqrt(g.x*g.x + g.y*g.y + g.z*g.z);

    // At this distance the dominant term is still GM/r², but centrifugal
    // and residual J₂ produce a ~0.09 m/s² difference.  Check that we
    // are close and that the vector still points inward.
    EXPECT_NEAR(mag, GM/(r*r), 0.12);
    EXPECT_LT(g.x, 0.0);
    EXPECT_NEAR(g.y, 0.0, 1e-3);
    EXPECT_NEAR(g.z, 0.0, 1e-3);
}

// ---------------------------------------------------------------------------
// effective_gravity – blend zone (exercises both helpers + smoothstep)
// ---------------------------------------------------------------------------
TEST_F(Wgs84Test, TestEffectiveGravity_BlendZoneLower) {
    // Just inside the blend window
    lla_t  lla  = {0.0, 0.0, H0 + 1000.0};
    ecef_t ecef;
    lla_to_ecef(lla, &ecef);

    vector3_t g = effective_gravity(lla, ecef);
    const double mag = std::sqrt(g.x*g.x + g.y*g.y + g.z*g.z);

    EXPECT_TRUE(std::isfinite(mag));
    EXPECT_GT(mag, 9.0);
    EXPECT_LT(mag, 9.8);
}

TEST_F(Wgs84Test, TestEffectiveGravity_BlendZoneMid) {
    // Middle of the blend window
    lla_t  lla  = {30.0 * DEG_TO_RAD, 45.0 * DEG_TO_RAD, 0.5*(H0+H1)};
    ecef_t ecef;
    lla_to_ecef(lla, &ecef);

    vector3_t g = effective_gravity(lla, ecef);
    const double mag = std::sqrt(g.x*g.x + g.y*g.y + g.z*g.z);

    EXPECT_TRUE(std::isfinite(mag));
    EXPECT_GT(mag, 8.5);
    EXPECT_LT(mag, 9.8);
}

TEST_F(Wgs84Test, TestEffectiveGravity_BlendZoneUpper) {
    // Just below the end of the blend window
    lla_t  lla  = {0.0, 0.0, H1 - 1000.0};
    ecef_t ecef;
    lla_to_ecef(lla, &ecef);

    vector3_t g = effective_gravity(lla, ecef);
    const double mag = std::sqrt(g.x*g.x + g.y*g.y + g.z*g.z);

    EXPECT_TRUE(std::isfinite(mag));
    EXPECT_GT(mag, 7.0);
    EXPECT_LT(mag, 9.5);
}

TEST_F(Wgs84Test, TestEffectiveGravity_BlendBoundaries) {
    // Exact boundaries – hit the if (h <= H0) and if (h >= H1) paths
    lla_t lla_lo = {0.0, 0.0, H0};
    lla_t lla_hi = {0.0, 0.0, H1};
    ecef_t ecef_lo, ecef_hi;
    lla_to_ecef(lla_lo, &ecef_lo);
    lla_to_ecef(lla_hi, &ecef_hi);

    vector3_t g_lo = effective_gravity(lla_lo, ecef_lo);
    vector3_t g_hi = effective_gravity(lla_hi, ecef_hi);

    EXPECT_TRUE(std::isfinite(g_lo.x));
    EXPECT_TRUE(std::isfinite(g_hi.x));
}

// ---------------------------------------------------------------------------
// Symmetry (odd Jₙ ≡ 0) – works in every regime
// ---------------------------------------------------------------------------
TEST_F(Wgs84Test, TestEffectiveGravity_OddJnAreZero_Surface) {
    lla_t  lla_n = { 30.0 * DEG_TO_RAD, 0.0, 0.0};
    lla_t  lla_s = {-30.0 * DEG_TO_RAD, 0.0, 0.0};
    ecef_t ecef_n, ecef_s;
    lla_to_ecef(lla_n, &ecef_n);
    lla_to_ecef(lla_s, &ecef_s);

    vector3_t g_n = effective_gravity(lla_n, ecef_n);
    vector3_t g_s = effective_gravity(lla_s, ecef_s);

    EXPECT_NEAR(g_n.x,  g_s.x, TOL_G);
    EXPECT_NEAR(g_n.y,  g_s.y, TOL_G);
    EXPECT_NEAR(g_n.z, -g_s.z, TOL_G);
}

TEST_F(Wgs84Test, TestEffectiveGravity_OddJnAreZero_High) {
    lla_t  lla_n = { 30.0 * DEG_TO_RAD, 0.0, 300000.0};
    lla_t  lla_s = {-30.0 * DEG_TO_RAD, 0.0, 300000.0};
    ecef_t ecef_n, ecef_s;
    lla_to_ecef(lla_n, &ecef_n);
    lla_to_ecef(lla_s, &ecef_s);

    vector3_t g_n = effective_gravity(lla_n, ecef_n);
    vector3_t g_s = effective_gravity(lla_s, ecef_s);

    EXPECT_NEAR(g_n.x,  g_s.x, TOL_G);
    EXPECT_NEAR(g_n.y,  g_s.y, TOL_G);
    EXPECT_NEAR(g_n.z, -g_s.z, TOL_G);
}

// ---------------------------------------------------------------------------
// Print helpers (coverage only)
// ---------------------------------------------------------------------------
TEST_F(Wgs84Test, TestPrintFunctions) {
    lla_t lla = {45.0 * DEG_TO_RAD, -120.0 * DEG_TO_RAD, 1234.5};
    ecef_t ecef;
    lla_to_ecef(lla, &ecef);
    vector3_t g = effective_gravity(lla, ecef);

    print_wgs84_model_values();
    print_lla(lla);
    print_ecef(ecef);
    print_vector3(g);

    SUCCEED();
}
