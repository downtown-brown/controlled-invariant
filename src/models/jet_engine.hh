static interval_t U(-0.35, 0.35);
static ninterval_t Omega_0 = {interval_t(-0.2, 0.2),
                              interval_t(-0.2, 0.2)};

inline C_Polyhedron A(nvec_t x_m, C_Polyhedron P) {
    const int64_t A[NDIM][NDIM] = {{10, -1}, {0, 10}};
    const int64_t A_den = 10;

    C_Polyhedron res = P;
    res.affine_image(Variable(0), A[0][0]*Variable(0) + A[0][1]*Variable(1), A_den);
    res.affine_image(Variable(1), A_den*A[1][0]*Variable(0)
                     + (A[0][0]*A[1][1] - A[0][1]*A[1][0])*Variable(1), A[0][0]*A_den);
    return res;
}

inline C_Polyhedron B(nvec_t x_m, interval_t U) {
    const int64_t B[NDIM] = {0, 1};
    const int64_t B_den = 10;

    const int64_t U_den = INT16_MAX;
    static int64_t Ul = rat_approx(U.lower(), U_den);
    static int64_t Uh = rat_approx(U.upper(), U_den);

    C_Polyhedron res(2, EMPTY);
    res.add_generator(point(B[0]*Ul*Variable(0) + B[1]*Ul*Variable(1), B_den*U_den));
    res.add_generator(point(B[0]*Uh*Variable(0) + B[1]*Uh*Variable(1), B_den*U_den));

    return res;
}

inline ninterval_t Phi(ninterval_t x, nvec_t x_m) {
    return {-0.15*pow(x[0], 2) - 0.05*pow(x[0], 2), interval_t(0,0)};
}

inline ninterval_t Psi(ninterval_t x, nvec_t x_m, interval_t U) {
    return {interval_t(0,0), interval_t(0,0)};
}
