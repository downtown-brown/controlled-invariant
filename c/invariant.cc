#include "invariant.hh"
#include <array>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <ppl.hh>
#include <vector>

using namespace std;
#define epsilon 1e-3

extern uint8_t stop;

uint64_t kk = 0;
I U(-0.1, 0.1);

vector<vector<C_Polyhedron>> U_approx(vector<IntervalData> Omega) {
    uint64_t num_int = 0;
    vector<vector<C_Polyhedron>> U;
    vector<IntervalData> N;
    vector<IntervalData> E;
    auto L = vector<IntervalData>(Omega);

    auto len = Omega.size();
    vector<C_Polyhedron> Omega_p(len);
    for (int i = 0; i < len; i++) {
        Omega_p[i] = Omega[i].poly;
    }
    auto Nc = convexhull(Omega_p);
    auto Nd = regiondiff(Nc, Omega_p.begin(), Omega_p.end());

    while (!L.empty()) {

        auto x = L.back();
        L.pop_back();
        num_int++;

        auto tmp = translate_into(x.P_u_over, x.P_over, Nc, Nd);
        U.push_back(tmp);
    }

    return U;
}

vector<IntervalData> I_approx(vector<IntervalData> Omega) {
    uint64_t num_int = 0;
    vector<IntervalData> S;
    vector<IntervalData> N;
    vector<IntervalData> E;
    auto L = vector<IntervalData>(Omega);

    auto len = Omega.size();
    vector<C_Polyhedron> Omega_p(len);
    for (int i = 0; i < len; i++) {
        Omega_p[i] = Omega[i].poly;
    }
    auto Nc = convexhull(Omega_p);
    auto Nd = regiondiff(Nc, Omega_p.begin(), Omega_p.end());

    while (!L.empty()) {

        auto x = L.back();
        L.pop_back();
        num_int++;

        // cout << "Considering x: " << endl;
        // print_points(x.poly);
        // print_points(x.P_over);
        // print_points(x.P_u_over);
        // print_points(Nc);
        // print_points(Nd);


        if (!intersects(x.P_over, Omega_p)) {
            //cout << "Putting in N\n";
            //print_points(x.P_over);
            //print_points(x.P_over);
            N.push_back(x);
        } else if (can_translate_into(x.P_u_over, x.P_over, Nc, Nd)) {
            //cout << "Putting in S\n";
            S.push_back(x);
        } else if (width(x.interval) < epsilon) {
            //cout << "Putting in E, width = " << width(x.interval) << endl;
            E.push_back(x);
        } else {
            //cout << "Bisecting\n";
            auto xs = bisect(x);
            L.push_back(get<0>(xs));
            L.push_back(get<1>(xs));
        }
    }

    cout << "N: " << N.size() << ", S: " << S.size() << ", E: " << E.size() << endl;
    fprint_points(N, "n" + to_string(kk) + ".txt");
    fprint_points(E, "e" + to_string(kk) + ".txt");
    fprint_points(S, "s" + to_string(kk) + ".txt");

    kk++;
    //cout << "Over\n";
    //print_over(S);

    stop = N.empty() && E.empty();

    cout << "considered " << num_int << " intervals\n";

    return S;
}


void print_points(vector<IntervalData> P) {
    for (auto p = P.begin(); p != P.end(); ++p) {
        print_points(p->poly);
    }
}

void fprint_points(vector<IntervalData> P, string fname) {
    for (auto p = P.begin(); p != P.end(); ++p) {
        fprint_points(p->poly, fname);
    }
}

void print_over(vector<IntervalData> P) {
    for (auto p = P.begin(); p != P.end(); ++p) {
        print_points(p->P_over);
    }
}

void print_u_over(vector<IntervalData> P) {
    for (auto p = P.begin(); p != P.end(); ++p) {
        print_points(p->P_u_over);
    }
}


pair<IntervalData, IntervalData> bisect(IntervalData x) {
    pair<I, I> tmp;
    tuple<I, I> l;
    tuple<I, I> r;
    I tmp2;
    uint8_t dim;
    if (width(get<0>(x.interval)) >= width(get<1>(x.interval))) {
        tmp = bisect(get<0>(x.interval));
        tmp2 = get<1>(x.interval);
        l = {get<0>(tmp), tmp2};
        r = {get<1>(tmp), tmp2};
    }
    else {
        tmp = bisect(get<1>(x.interval));
        tmp2 = get<0>(x.interval);
        l = {tmp2, get<0>(tmp)};
        r = {tmp2, get<1>(tmp)};
    }

    IntervalData x_l(l);
    x.lchild = &x_l;
    IntervalData x_r(r);
    x.lchild = &x_r;

    return {x_l, x_r};
}

double width(tuple<I, I> interval) {
    return max(width(get<0>(interval)),
               width(get<1>(interval)));
}

IntervalData::IntervalData(tuple<I, I> x) {
    interval = x;
    array<double, 2> x_m = {median(get<0>(interval)), median(get<0>(interval))};
    poly = i2p(interval);

    //cout << "Begin interval data:\n\n";
    //print_points(poly);

    //print_points(A(x_m, poly));
    //print_points(B(x_m, U));
    //print_points(i2p(Phi(interval, x_m)));
    //print_points(i2p(Psi(interval, x_m, U)));
    //print_points(i2p(Phi(interval, x_m)) + i2p(Psi(interval, x_m, U)));
    //cout << "End interval data:\n\n";

    P_u_over = A(x_m, poly) + i2p(Phi(interval, x_m)) + i2p(Psi(interval, x_m, U));
    P_over = P_u_over + B(x_m, U);

    //print_points(A(x_m, poly) + i2p(Phi(interval, x_m)));
    //print_points(P_over);
    lchild = NULL;
    rchild = NULL;
    invariant = true;
    iter = 0;
}

#define mu 0.5
array<double, 2> midpoint(tuple<I, I> interval) {
    return {median(get<0>(interval)),
            median(get<1>(interval))};
}

#define dt 0.01
#define m 0.2
#define g 9.8
#define l 0.3
#define J 0.006
#define b 0.1

C_Polyhedron A(array<double, 2> x, C_Polyhedron P) {
    auto res = C_Polyhedron(P);
    int64_t A1_num = 1;
    int64_t A1_den = 1;
    int64_t A2_num = 1;
    int64_t A2_den = 100;
    int64_t A3_num;
    int64_t A3_den;
    rat_approx(m*g*l/J*dt*cos(x[0]), INT16_MAX, &A3_num, &A3_den);
    int64_t A4_num;
    int64_t A4_den;
    rat_approx((1 - dt*b/J), INT16_MAX, &A4_num, &A4_den);

    int64_t A1 = A1_num*A2_den*A3_den*A4_den;
    int64_t A2 = A1_den*A2_num*A3_den*A4_den;
    int64_t A3 = A1_den*A2_den*A3_num*A4_den;
    int64_t A4 = A1_den*A2_den*A3_den*A4_num;
    int64_t den = A1_den*A2_den*A3_den*A4_den;

    res.affine_image(Variable(0), A1*Variable(0) + A2*Variable(1), den);

    res.affine_image(Variable(1), den*A3*Variable(0) + (A1*A4 - A2*A3)*Variable(1), A1*den);
    return res;
}

C_Polyhedron B(array<double, 2> x, I U) {
    int64_t B1_num = 0;
    int64_t B1_den = 1;
    int64_t B2_num;
    int64_t B2_den;
    rat_approx(dt*l/J*cos(x[0]), INT16_MAX, &B2_num, &B2_den);

    int64_t nl0, dl0;
    rat_approx(U.lower(), INT32_MAX, &nl0, &dl0);
    int64_t nu0, du0;
    rat_approx(U.upper(), INT32_MAX, &nu0, &du0);

    C_Polyhedron res(2, EMPTY);
    res.add_generator(point(B1_num*nl0*Variable(0) + B2_num*nl0*Variable(1), dl0*B2_den));
    res.add_generator(point(B1_num*nu0*Variable(0) + B2_num*nu0*Variable(1), du0*B2_den));

    return res;
}

tuple<I, I> Phi(tuple<I, I> X, array<double, 2> x) {
    I Phi2 = dt*m*g*l/J*sin(get<0>(X) - cos(x[0])*x[0]);
    return {I(0,0), Phi2};
}

tuple<I, I> Psi(tuple<I, I> X, array<double, 2> x, I U) {
    I Psi2 = U*dt*l/J*(cos(get<0>(X)) - cos(x[0]));
    return {I(0,0), Psi2};
}
