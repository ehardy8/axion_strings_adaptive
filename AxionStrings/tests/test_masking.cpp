#include "doctest.h"

#include "../Masking.hpp"

#include <numbers>

// milestone-1.md task 1.7, T2 (conventions.md sec.13): "a static string does
// not test the mask" -- both tests below use configurations with Pi != 0,
// so neither is blind to the failure mode a static test would miss.

namespace
{
constexpr double pi = std::numbers::pi;
}

// --- T2a: direct array-level unit test -------------------------------
// "Assert that the buffer written for the FFT is zero at precisely the set
// of points where the mask predicate holds, and bit-identical to the
// unmasked field elsewhere." (conventions.md sec.13)

TEST_CASE("T2a: scheme B is exactly zero below threshold, exactly the "
          "unmasked value at or above it")
{
    const double R         = 2.0;
    MaskingParams params{};
    params.scheme    = MaskingScheme::B;
    params.threshold = 0.8;

    // A handful of synthetic (psi1, psi2, Pi1, Pi2) points spanning below,
    // at, and above the threshold, all with Pi != 0 (not the degenerate
    // static case).
    struct Point
    {
        double psi1, psi2, Pi1, Pi2;
    };
    const Point points[] = {
        {0.1 * R, 0.0, 0.3, -0.2},        // mod = 0.1: well below threshold
        {0.7999 * R, 0.0, 0.3, -0.2},     // mod just below threshold
        {0.8 * R, 0.0, 0.3, -0.2},        // mod exactly at threshold
        {0.8001 * R, 0.0, 0.3, -0.2},     // mod just above threshold
        {0.95 * R, 0.2 * R, -0.4, 0.5},   // mod > threshold
        {1.0 * R, 0.0, 0.1, 0.1},         // mod = 1 (vacuum)
    };

    for (const auto &p : points)
    {
        const double mod = std::sqrt(p.psi1 * p.psi1 + p.psi2 * p.psi2) / R;
        const double raw = p.psi1 * p.Pi2 - p.Pi1 * p.psi2;
        const double unmasked = raw / (p.psi1 * p.psi1 + p.psi2 * p.psi2);
        const double masked =
            masked_a_dot(params, p.psi1, p.psi2, p.Pi1, p.Pi2, R);

        CAPTURE(mod);
        if (mod < params.threshold)
        {
            CHECK(masked == 0.0);
            CHECK(masking_weight(params, p.psi1, p.psi2, R) == 0.0);
        }
        else
        {
            CHECK(masked == doctest::Approx(unmasked).epsilon(1.0e-12));
            CHECK(masking_weight(params, p.psi1, p.psi2, R) == 1.0);
        }
    }
}

TEST_CASE("T2a: MaskingScheme::None matches the unmasked formula everywhere,"
          " even deep inside where B would zero it")
{
    const double R = 1.5;
    MaskingParams none_params{};
    none_params.scheme = MaskingScheme::None;

    const double psi1 = 0.05 * R;
    const double psi2 = 0.0;
    const double Pi1  = 0.2;
    const double Pi2  = -0.3;

    const double raw       = psi1 * Pi2 - Pi1 * psi2;
    const double unmasked  = raw / (psi1 * psi1 + psi2 * psi2);
    const double masked_none = masked_a_dot(none_params, psi1, psi2, Pi1, Pi2, R);

    CHECK(masked_none == doctest::Approx(unmasked).epsilon(1.0e-12));
    CHECK(masking_weight(none_params, psi1, psi2, R) == 1.0);
}

TEST_CASE("T2a: scheme A is the bare numerator, undivided, no top-hat")
{
    const double R = 1.0;
    MaskingParams a_params{};
    a_params.scheme = MaskingScheme::A;

    const double psi1 = 0.05; // deep inside where B would zero it
    const double psi2 = 0.0;
    const double Pi1  = 0.2;
    const double Pi2  = -0.3;

    const double raw = psi1 * Pi2 - Pi1 * psi2;
    CHECK(masked_a_dot(a_params, psi1, psi2, Pi1, Pi2, R) ==
          doctest::Approx(raw).epsilon(1.0e-12));
    // Scheme A has no top-hat (unlike B, nothing is ever hard-zeroed), but
    // it is not a silent no-op either: conventions.md's smooth weight
    // f = (1 + r/f_a)^2, r = |phi|-v = mod-1, still suppresses this
    // deep-inside-the-core point substantially below 1.
    const double f_a      = 1.4142135623730951; // sqrt(2) v, v = 1
    const double mod      = std::sqrt(psi1 * psi1 + psi2 * psi2) / R;
    const double f         = 1.0 + (mod - 1.0) / f_a;
    CHECK(masking_weight(a_params, psi1, psi2, R) ==
          doctest::Approx(f * f).epsilon(1.0e-12));
}

// --- T2b: boosted straight string -------------------------------------
// "A string moving with known velocity has theta_dot != 0 near the core,
// so correct masking must zero a known moving cylinder." (conventions.md
// sec.13) Analytic profile psi = R tanh(m rho) e^{i theta}, rho/theta
// built from (x - beta*tau, y): both rho and theta depend on tau only
// through this translation, so Pi_i = d(psi_i)/d(tau) is computable exactly
// by the chain rule -- no numerical differentiation, no need for the real
// RHS/Background machinery, since this tests the masking module in
// isolation.

namespace
{
struct BoostedStringState
{
    double psi1, psi2, Pi1, Pi2;
};

BoostedStringState boosted_string(double x, double y, double tau, double R0,
                                  double m0, double beta)
{
    const double X   = x - beta * tau;
    const double rho = std::sqrt(X * X + y * y);
    const double theta = std::atan2(y, X);

    const double amp = R0 * std::tanh(m0 * rho);
    const double cos_theta = std::cos(theta);
    const double sin_theta = std::sin(theta);

    BoostedStringState s{};
    s.psi1 = amp * cos_theta;
    s.psi2 = amp * sin_theta;

    if (rho == 0.0)
    {
        s.Pi1 = 0.0;
        s.Pi2 = 0.0;
        return s;
    }

    // d(rho)/d(tau) = -beta * X / rho = -beta cos(theta)
    // d(theta)/d(tau) = beta * sin(theta) / rho
    const double rho_dot   = -beta * cos_theta;
    const double theta_dot = beta * sin_theta / rho;
    const double sech2     = 1.0 - std::tanh(m0 * rho) * std::tanh(m0 * rho);
    const double amp_dot   = R0 * m0 * sech2 * rho_dot;

    // psi = amp * (cos theta, sin theta); d/dtau via product rule.
    s.Pi1 = amp_dot * cos_theta - amp * sin_theta * theta_dot;
    s.Pi2 = amp_dot * sin_theta + amp * cos_theta * theta_dot;
    return s;
}
} // namespace

TEST_CASE("T2b: boosted string has Pi != 0 (not the degenerate static case)")
{
    const auto s = boosted_string(1.0, 0.3, 0.0, /*R0=*/1.0, /*m0=*/2.0,
                                  /*beta=*/0.4);
    CHECK((std::abs(s.Pi1) > 1.0e-8 || std::abs(s.Pi2) > 1.0e-8));
}

TEST_CASE("T2b: scheme B zeros exactly the moving core, and only there, at "
          "two different times")
{
    const double R0        = 1.0;
    const double m0         = 3.0;
    const double beta       = 0.4;
    MaskingParams params{};
    params.scheme    = MaskingScheme::B;
    params.threshold = 0.8;

    // A point that sits inside the core at tau=0 (core at x=0) but well
    // outside it at tau=5 (core at x=2), and vice versa for a second point
    // -- demonstrating the masked-out region tracks the moving core rather
    // than a fixed location. y != 0 for the "away from core" cases: along
    // y=0 the phase-rotation term vanishes by symmetry (psi and Pi both
    // point along theta=0 there) regardless of masking, which would test
    // that symmetry rather than the mask.
    struct Case
    {
        double x, y, tau;
        bool expect_masked;
    };
    const Case cases[] = {
        {0.0, 0.0, 0.0, true},   // at the core at tau=0 (core at x=0)
        {0.0, 0.6, 5.0, false},  // far from the core (now at x=2) at tau=5
        {2.0, 0.0, 5.0, true},   // at the core at tau=5
        {2.0, 0.6, 0.0, false}, // far from the core (still at x=0) at tau=0
    };

    for (const auto &c : cases)
    {
        const auto s = boosted_string(c.x, c.y, c.tau, R0, m0, beta);
        const double mod =
            std::sqrt(s.psi1 * s.psi1 + s.psi2 * s.psi2) / R0;
        const double masked =
            masked_a_dot(params, s.psi1, s.psi2, s.Pi1, s.Pi2, R0);

        CAPTURE(c.x);
        CAPTURE(c.tau);
        CAPTURE(mod);
        if (c.expect_masked)
        {
            CHECK(mod < params.threshold);
            CHECK(masked == 0.0);
        }
        else
        {
            CHECK(mod >= params.threshold);
            // Away from the core the unmasked phase rate is manifestly
            // non-zero (the naive rigid translation's phase has a slowly
            // decaying tail, unlike a static string's identically-zero
            // masked a_dot) -- and scheme B must reproduce it exactly,
            // not zero it.
            const double raw =
                s.psi1 * s.Pi2 - s.Pi1 * s.psi2;
            const double unmasked =
                raw / (s.psi1 * s.psi1 + s.psi2 * s.psi2);
            CHECK(std::abs(unmasked) > 1.0e-6);
            CHECK(masked == doctest::Approx(unmasked).epsilon(1.0e-12));
        }
    }
}
