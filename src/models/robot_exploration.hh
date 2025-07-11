#include <cstdint>
#include <math.h>
#include <ppl.hh>

static const string DATA_DIR = "data_exploration/";
static const nvec_t epsilon = {5e-2, 5e-2, M_PI/32};

static ninterval_t U = {interval_t(-2, 2), interval_t(-1, 1),
                        interval_t(0, 0)};


static vector<IntervalData> Omega_0 = {
    IntervalData({interval_t(-4, -1), interval_t(-2, 2), interval_t(-M_PI, M_PI)}),
    IntervalData({interval_t(1, 4), interval_t(-2, 2), interval_t(-M_PI, M_PI)}),
    IntervalData({interval_t(-2, 2), interval_t(-4, -2), interval_t(-M_PI, M_PI)}),
    IntervalData({interval_t(-2, 2), interval_t(2, 4), interval_t(-M_PI, M_PI)}),
    IntervalData({interval_t(-1, 1), interval_t(1, 2), interval_t(-M_PI, M_PI)}),
    IntervalData({interval_t(-1, 1), interval_t(-2, -1), interval_t(-M_PI, M_PI)}),
};

static vector<ninterval_t> N_0 = {
    {interval_t(-4, -2), interval_t(-4, -2), interval_t(-M_PI, M_PI)},
    {interval_t(-4, -2), interval_t(2, 4), interval_t(-M_PI, M_PI)},
    {interval_t(2, 4), interval_t(-4, -2), interval_t(-M_PI, M_PI)},
    {interval_t(2, 4), interval_t(2, 4), interval_t(-M_PI, M_PI)},
    {interval_t(-1, 1), interval_t(-1, 1), interval_t(-M_PI, M_PI)},
};


inline C_Polyhedron A(nvec_t x_m, C_Polyhedron P) {
    return P;
}

inline C_Polyhedron B(nvec_t x_m, ninterval_t U) {
    const int64_t den = 10000000;
    const int64_t B_den = den*5;
    int64_t B[NDIM][NDIM] = {
        {rat_approx(cos(x_m[2]), den), 0, 0},
        {rat_approx(sin(x_m[2]), den), 0, 0},
        {0, den, 0},
    };

    C_Polyhedron U_p = i2p(U);
    U_p.affine_image(Variable(2), Variable(1)*B[2][1], B_den);
    U_p.affine_image(Variable(1), Variable(0)*B[1][0], B_den);
    U_p.affine_image(Variable(0), Variable(0)*B[0][0], B_den);

    return U_p;
}

inline ninterval_t Phi(ninterval_t x, nvec_t x_m) {
    return {interval_t(0,0), interval_t(0, 0)};
}

inline ninterval_t Psi(ninterval_t x, nvec_t x_m, ninterval_t U) {
    return {
        0.2*U[0]*(cos(x[2]) - cos(x_m[2])),
        0.2*U[0]*(sin(x[2]) - sin(x_m[2])),
        interval_t(0, 0),
    };
}

ninterval_t Delta = {interval_t(-0.05,0.05), interval_t(-0.05,0.05), interval_t(0,0)};
