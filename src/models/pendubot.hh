#include <cstdint>
#include <math.h>

static const string DATA_DIR = "data_pendubot/";
static const nvec_t epsilon = {1e-1, 1e-1, 1e-1, 1e-1};

static interval_t U(-2, 2);
static vector<IntervalData> Omega_0 = {
    IntervalData({interval_t(-1, 1), interval_t(-1, 1), interval_t(-1, 1), interval_t(-1, 1)})
};

static vector<ninterval_t> N_0;

inline C_Polyhedron A(nvec_t x_m, C_Polyhedron P) {
    return P;
}

inline C_Polyhedron B(nvec_t x_m, interval_t U) {
    const int64_t B[NDIM] = {
        0,
        rat_approx(-10.096*x_m[2]*x_m[2] + 44.252, INT16_MAX),
        0,
        rat_approx(37.802*x_m[2]*x_m[2] - 83.912, INT16_MAX)
    };
    const int64_t B_den = INT16_MAX*10;

    const int64_t U_den = INT16_MAX;
    static int64_t Ul = rat_approx(U.lower(), U_den);
    static int64_t Uh = rat_approx(U.upper(), U_den);

    C_Polyhedron res(NDIM, EMPTY);
    res.add_generator(point(B[0]*Ul*Variable(0) + B[1]*Ul*Variable(1)
                            + B[2]*Ul*Variable(2) + B[3]*Ul*Variable(3),  B_den*U_den));
    res.add_generator(point(B[0]*Uh*Variable(0) + B[1]*Uh*Variable(1)
                            + B[2]*Uh*Variable(2) + B[3]*Uh*Variable(3),  B_den*U_den));

    return res;
}

inline ninterval_t Phi(ninterval_t x, nvec_t x_m) {
  return {
    0.1*x[1],
    0.1*(-10.656*pow(x[0], 3) + 11.531*pow(x[0], 2)*x[2] + 7.885*x[0]*pow(x[2], 2)
         + 0.797*pow(x[1], 2)*x[2] + 0.841*x[1]*x[2]*x[3] + 21.049*pow(x[2], 3)
         + 0.420*x[2]*pow(x[3], 2) + 66.523*x[0] - 24.511*x[2]),
    0.1*x[3],
    0.1*(10.996*pow(x[0], 3) - 48.915*pow(x[0], 2)*x[2] - 6.404*x[0]*pow(x[2],2)
         - 2.396*pow(x[1], 2)*x[2] - 1.594*x[1]*x[2]*x[3] - 51.909*pow(x[2], 3)
         - 0.797*x[2]*pow(x[3], 2) - 68.642*x[0] + 103.9783*x[2])
  };
}

inline ninterval_t Psi(ninterval_t x, nvec_t x_m, interval_t U) {
  return {
    interval_t(0, 0),
    U*0.1*(-10.096*(x[2]*x[2] - x_m[2]*x_m[2])),
    interval_t(0, 0),
    U*0.1*(37.802*(x[2]*x[2] - x_m[2]*x_m[2]))
  };
}

ninterval_t Delta = {interval_t(0,0), interval_t(0,0), interval_t(0,0), interval_t(0,0)};
