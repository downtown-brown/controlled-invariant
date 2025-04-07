#include <math.h>

static const string DATA_DIR = "data_cartpole/";
static const nvec_t epsilon = {1.5e-1, 1.5e-1, 1.5e-1, 1.5e-1};

static interval_t U(-2, 2);
static ninterval_t Omega_0 = {
    interval_t(-.5, .5), interval_t(-0.75 * M_PI, 0.75 * M_PI), interval_t(-1, 1),
    interval_t(-1.5 * M_PI, 1.5 * M_PI)
};

inline C_Polyhedron A(nvec_t x_m, C_Polyhedron P) {
    const int64_t A[NDIM][NDIM] = {
        {25000,     0,  2500,     0},
        {    0, 25000,     0,  2500},
        {    0,  9810, 25000,  -981},
        {    0, 68670,     0, 18133}
    };
    const int64_t A_den = 25000;

    C_Polyhedron res = P;
    res.affine_image(Variable(0), A[0][0] * Variable(0) + A[0][1] * Variable(1)
                     + A[0][2] * Variable(2) + A[0][3] * Variable(3), A_den);
    res.affine_image(Variable(0), A[1][0] * Variable(0) + A[1][1] * Variable(1)
                     + A[1][2] * Variable(2) + A[1][3] * Variable(3), A_den);
    res.affine_image(Variable(0), A[2][0] * Variable(0) + A[2][1] * Variable(1)
                     + A[2][2] * Variable(2) + A[2][3] * Variable(3), A_den);
    res.affine_image(Variable(0), A[3][0] * Variable(0) + A[3][1] * Variable(1)
                     + A[3][2] * Variable(2) + A[3][3] * Variable(3), A_den);
    return res;
}

inline C_Polyhedron B(nvec_t x_m, interval_t U) {
    const int64_t B[NDIM] = {0, 0, 2, 4};
    const int64_t B_den = 50;

    const int64_t U_den = INT16_MAX;
    static int64_t Ul = rat_approx(U.lower(), U_den);
    static int64_t Uh = rat_approx(U.upper(), U_den);
q
    C_Polyhedron res(NDIM, EMPTY);
    res.add_generator(point(B[0]*Ul*Variable(0) + B[1]*Ul*Variable(1)
                            + B[2]*Ul*Variable(2) + B[3]*Ul*Variable(3),  B_den*U_den));
    res.add_generator(point(B[0]*Uh*Variable(0) + B[1]*Uh*Variable(1)
                            + B[2]*Uh*Variable(2) + B[3]*Uh*Variable(3),  B_den*U_den));

    return res;
}

inline ninterval_t Phi(ninterval_t x, nvec_t x_m) {
    return {interval_t(0,0), interval_t(0,0)};
}

inline ninterval_t Psi(ninterval_t x, nvec_t x_m, interval_t U) {
  return {
    interval_t(0, 0),
    interval_t(0, 0),
    -0.4*U - 0.3924*x[1] - (0.1*(- 0.5*sin(x[1])*x[3]*x[3] + U + 9.81*cos(x[1])*sin(x[1])))/(cos(x[1])*cos(x[1]) - 3.5),
    -0.8*U - 2.7468*x[1] - (0.002*(- 100.0*cos(x[1])*sin(x[1])*x[3]*x[3] + 6867.0*sin(x[1]) + 200.0*U*cos(x[1])))/(2.0*cos(x[1])*cos(x[1]) - 7.0)
  };
}

ninterval_t Delta = {interval_t(0,0), interval_t(0,0), interval_t(0,0), interval_t(0,0)};
