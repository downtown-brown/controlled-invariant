#include <cstdint>
#include <math.h>

static const string DATA_DIR = "data_cartpole/";
static const nvec_t epsilon = {0.2, 0.2, 5, 10};

static ninterval_t U = {interval_t(5, 35), interval_t(-8500, 0), interval_t(0, 0), interval_t(0, 0)};
static ninterval_t Omega_0 = {interval_t(0, 5), interval_t(0, 5), interval_t(50, 150), interval_t(0, 150)};

static vector<ninterval_t> N_0;

inline C_Polyhedron A(nvec_t x_m, C_Polyhedron P) {
    const int64_t A[NDIM][NDIM] = {
        {1000,     0,     0,    0},
        {    0, 1000,     0,    0},
        {    0,     0, 1000-308,    308},
        {    0,     0,    867, 1000-867}
    };
    const int64_t A_den = 1000;

    C_Polyhedron res = P;
    res.affine_image(Variable(0), A[0][0] * Variable(0) + A[0][2] * Variable(2), A_den);
    res.affine_image(Variable(1), A[1][1] * Variable(1) + A[1][3] * Variable(3), A_den);
    res.affine_image(Variable(2), A[2][2] * Variable(2), A_den);
    res.affine_image(Variable(3), A[3][3] * Variable(3), A_den);
    return res;
}

inline C_Polyhedron B(nvec_t x_m, ninterval_t U) {
    const int64_t den = 10000000;
    const int64_t B_den = den*10;
    int64_t B[NDIM][NDIM] = {
        {rat_approx(x_m[0] - 5.1, den), 0, 0, 0},
        {rat_approx(x_m[1], den), 0, 0, 0},
        {rat_approx(x_m[2] - 130.0 , den), 0, 0, 0},
        {0, 1000000, 0, 0},
    };

    C_Polyhedron U_p = i2p(U);
    U_p.affine_image(Variable(3), Variable(1)*B[3][1], B_den);
    U_p.affine_image(Variable(2), Variable(0)*B[2][0], B_den);
    U_p.affine_image(Variable(1), Variable(0)*B[1][0], B_den);
    U_p.affine_image(Variable(0), Variable(0)*B[0][0], B_den);

    return U_p;

}

inline ninterval_t Phi(ninterval_t x, nvec_t x_m) {
  return {
      -1.29e+12*x[0]*exp(-9760.0/(x[2] + 273.0)) - 9.04e+9*x[0]*x[0]*exp(-8560.0/(x[2] + 273.0)),
      1.29e+12*x[0]*exp(-9760.0/(x[2] + 273.0)) - 1.29e+12*x[1]*exp(-9760.0/(x[2] + 273.0)),
      1.35e+11*exp(-8560.0/(x[2] + 273.0))*x[0]*x[0] - 1.92e+12*exp(-9760.0/(x[2] + 273.0))*x[0] + 5.03e+12*x[1]*exp(-9760.0/(x[2] + 273.0)),
      interval_t(0, 0)
  };
}

inline ninterval_t Psi(ninterval_t x, nvec_t x_m, ninterval_t U) {
  return {
      -1.0*U[0]*(x[0]-x_m[0]),
      -1.0*U[0]*(x[1]-x_m[1]),
      -1.0*U[0]*(x[2]-x_m[2]),
      interval_t(0, 0)
  };
}

ninterval_t Delta = {interval_t(0,0), interval_t(0,0), interval_t(0,0), interval_t(0,0)};
