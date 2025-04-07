#define dt 0.01
#define m 0.2
#define g 9.8
#define l 0.3
#define J 0.006
#define b 0.1

static interval_t U(-0.1, 0.1);
static ninterval_t Omega_0 = {interval_t(-0.05, 0.05),
                              interval_t(-0.01, 0.01)};

C_Polyhedron A(nvec_t x_m, C_Polyhedron P) {
    const int64_t A_den = 1000000;
    int64_t A[NDIM][NDIM] = {{A_den, A_den/100},
                             {rat_approx(m*g*l/J*dt*cos(x_m[0]), A_den),
                              rat_approx((1 - dt*b/J), A_den)}};

    C_Polyhedron res = P;
    res.affine_image(Variable(0), A[0][0]*Variable(0) + A[0][1]*Variable(1), A_den);
    res.affine_image(Variable(1), A_den*A[1][0]*Variable(0)
                     + (A[0][0]*A[1][1] - A[0][1]*A[1][0])*Variable(1), A[0][0]*A_den);

    return res;
}

C_Polyhedron B(nvec_t x_m, interval_t U) {
    const int64_t B_den = 1000000;
    int64_t B[NDIM] = {0, rat_approx(dt*l/J*cos(x_m[0]), B_den)};

    const int64_t U_den = INT16_MAX;
    static int64_t Ul = rat_approx(U.lower(), U_den);
    static int64_t Uh = rat_approx(U.upper(), U_den);

    C_Polyhedron res(2, EMPTY);
    res.add_generator(point(B[0]*Ul*Variable(0) + B[1]*Ul*Variable(1), B_den*U_den));
    res.add_generator(point(B[0]*Uh*Variable(0) + B[1]*Uh*Variable(1), B_den*U_den));

    return res;
}

ninterval_t Phi(ninterval_t x, nvec_t x_m) {
    return {interval_t(0,0),
            dt*m*g*l/J*sin(x[0] - cos(x_m[0])*x_m[0])};
}

ninterval_t Psi(ninterval_t x, nvec_t x_m, interval_t U) {
    return {interval_t(0,0),
            U*dt*l/J*(cos(x[0]) - cos(x_m[0]))};
}
