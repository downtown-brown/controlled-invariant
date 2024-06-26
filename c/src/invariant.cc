#include "invariant.hh"
#include <array>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <ppl.hh>
#include <vector>

using namespace std;
#define epsilon 3e-2

extern uint8_t stop;

uint64_t kk = 0;
I U(-2, 2);

vector<vector<C_Polyhedron>> U_approx(list<IntervalData> Omega) {
    uint64_t num_int = 0;
    vector<vector<C_Polyhedron>> U;
    list<IntervalData> N;
    list<IntervalData> E;
    auto L = list<IntervalData>(Omega);

    auto len = Omega.size();
    vector<C_Polyhedron> Omega_p(len);
    int i = 0;
    for (const auto& O : Omega) {
        Omega_p[i] = O.poly;
        i++;
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

list<IntervalData> I_accel(const list<IntervalData>& Omega) {
    uint64_t num_int = 0;
    list<IntervalData> S;
    list<IntervalData> N;
    list<IntervalData> E;
    deque<IntervalData> L;

    uint32_t j = 0;

    move(
         begin(Omega),
         end(Omega),
         back_inserter(L)
    );

    auto len = L.size();
    vector<C_Polyhedron> L_p(len);
    for (uint64_t i = 0; i < len; i++) {
        L_p[i] = L[i].poly;
    }
    auto Nc = convexhull(L_p);
    auto Nd = regiondiff(Nc, L_p.begin(), L_p.end());

    int jj = 0;
    while (1) {
        jj++;
        auto x = L.front();
        L.pop_front();
        num_int++;

        if (!intersects(x.P_over, Nc)) {
            N.push_back(x);
            Nd.push_back(x.poly);
            j++;
        } else if (can_translate_into(x.P_u_over, x.P_over, Nc, Nd)) {
            x.checked = j;
            L.push_back(x);
        } else if (!wider_than(x.interval)) {
            E.push_back(x);
            Nd.push_back(x.poly);
            j++;
        } else {
            auto xs = bisect(x);
            L.push_front(get<0>(xs));
            L.push_front(get<1>(xs));
        }

        bool stop = true;
        for (const auto& i : L) {
            if (i.checked != j) {
                stop = false;
                break;
            }
        }

        if (stop)
            break;

        if (jj % 250 == 0) {
            cout << "Iteration: " << jj << endl;
        }
    }


    move(
         begin(L),
         end(L),
         back_inserter(S)
    );

    cout << "N: " << N.size() << ", S: " << S.size() << ", E: " << E.size() << endl;
    fprint_points(N, DATADIR + "n" + to_string(kk) + ".txt");
    fprint_points(E, DATADIR + "e" + to_string(kk) + ".txt");
    fprint_points(S, DATADIR + "s" + to_string(kk) + ".txt");

    kk++;
    //cout << "Over\n";
    //print_over(S);

    stop = N.empty() && E.empty();

    cout << "considered " << num_int << " intervals\n";

    return S;
}

list<IntervalData> I_approx(const list<IntervalData>& Omega) {
    uint64_t num_int = 0;
    list<IntervalData> S;
    list<IntervalData> N;
    list<IntervalData> E;
    auto L = list<IntervalData>(Omega);

    auto len = Omega.size();
    vector<C_Polyhedron> Omega_p(len);
    int i = 0;
    for (const auto& O : Omega) {
        Omega_p[i] = O.poly;
        i++;
    }
    auto Nc = convexhull(Omega_p);
    auto Nd = regiondiff(Nc, Omega_p.begin(), Omega_p.end());

    while (!L.empty()) {

        auto x = L.back();
        L.pop_back();
        num_int++;

        if (!intersects(x.P_over, Nc)) {
            N.push_back(x);
            x.status = STATUS_NOT_IN;
        } else if (can_translate_into(x.P_u_over, x.P_over, Nc, Nd)) {
            S.push_back(x);
            x.status = STATUS_IN;
        } else if (!wider_than(x.interval)) {
            E.push_back(x);
            x.status = STATUS_BOUNDARY;
        } else {
            x.status = STATUS_UNDETERMINED;
            auto xs = bisect(x);
            L.push_back(get<0>(xs));
            L.push_back(get<1>(xs));
        }
    }

    //merge(N);
    merge(S);
    //merge(E);

    cout << "N: " << N.size() << ", S: " << S.size() << ", E: " << E.size() << endl;

    fprint_points(N, DATADIR + "n" + to_string(kk) + ".txt");
    fprint_points(E, DATADIR + "e" + to_string(kk) + ".txt");
    fprint_points(S, DATADIR + "s" + to_string(kk) + ".txt");

    kk++;

    stop = N.empty() && E.empty();

    cout << "considered " << num_int << " intervals\n";

    return S;
}



pair<IntervalData, IntervalData> bisect(IntervalData x) {
    double max_width = 0;
    int bisect_dim;
    for (int i = 0; i < n; i++) {
        auto curr_width = width(x.interval[i]);
        if (curr_width > max_width) {
            max_width = curr_width;
            bisect_dim = i;
        }
    }

    auto l = x.interval;
    auto r = x.interval;

    auto tmp = bisect(x.interval[bisect_dim]);

    l[bisect_dim] = get<0>(tmp);
    r[bisect_dim] = get<1>(tmp);

    IntervalData x_l(l);
    x.lchild = &x_l;
    IntervalData x_r(r);
    x.rchild = &x_r;

    return {x_l, x_r};
}

bool wider_than(nI interval) {
    for (int i = 0; i < n; i++) {
        if (width(interval[i]) > epsilon) {
            return true;
        }
    }

    return false;
}

IntervalData::IntervalData(nI x) {
    interval = x;
    auto x_m = median(interval);

    poly = i2p(interval);

    P_u_over = A(x_m, poly) + i2p(Phi(interval, x_m) + Psi(interval, x_m, U));
    P_over = P_u_over + B(x_m, U);

    lchild = NULL;
    rchild = NULL;
    status = STATUS_UNDETERMINED;
    iter = 0;
}

array<double, n> median(nI interval) {
    array<double, n> res;

    for (int i = 0; i < n; i++) {
        res[i] = median(interval[i]);
    }

    return res;
}

/*
#define mu 0.5
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

nI Phi(nI X, array<double, 2> x) {
    I Phi2 = dt*m*g*l/J*sin(get<0>(X) - cos(x[0])*x[0]);
    return {I(0,0), Phi2};
}

nI Psi(nI X, array<double, 2> x, I U) {
    I Psi2 = U*dt*l/J*(cos(get<0>(X)) - cos(x[0]));
    return {I(0,0), Psi2};
}
*/



C_Polyhedron A(array<double, 2> x, C_Polyhedron P) {
    auto res = C_Polyhedron(P);
    const int64_t A1_num = 10-1;
    const int64_t A1_den = 10;
    const int64_t A2_num = 2;
    const int64_t A2_den = 1;
    const int64_t A3_num = 10-3;
    const int64_t A3_den = 10;
    const int64_t A4_num = 4;
    const int64_t A4_den = 1;

    static int64_t A1 = A1_num*A2_den*A3_den*A4_den;
    static int64_t A2 = A1_den*A2_num*A3_den*A4_den;
    static int64_t A3 = A1_den*A2_den*A3_num*A4_den;
    static int64_t A4 = A1_den*A2_den*A3_den*A4_num;
    static int64_t den = A1_den*A2_den*A3_den*A4_den;

    res.affine_image(Variable(0), A1*Variable(0) + A2*Variable(1), den);

    res.affine_image(Variable(1), den*A3*Variable(0) + (A1*A4 - A2*A3)*Variable(1), A1*den);
    return res;
}

C_Polyhedron B(array<double, 2> x, I U) {
    const int64_t B1_num = 0;
    //const int64_t B1_den = 1;
    const int64_t B2_num = -2;
    const int64_t B2_den = 10;

    static uint8_t i = 0;
    static int64_t nl0, dl0;
    static int64_t nu0, du0;
    if (!i) {
        rat_approx(U.lower(), INT32_MAX, &nl0, &dl0);
        rat_approx(U.upper(), INT32_MAX, &nu0, &du0);
        i = 1;
    }

    C_Polyhedron res(2, EMPTY);
    res.add_generator(point(B1_num*nl0*Variable(0) + B2_num*nl0*Variable(1), dl0*B2_den));
    res.add_generator(point(B1_num*nu0*Variable(0) + B2_num*nu0*Variable(1), du0*B2_den));

    return res;
}

nI Phi(nI X, array<double, 2> x) {
    I Phi2 = -0.25*pow(x[1],3)/10;
    return {I(0,0), Phi2};
}

nI Psi(nI X, array<double, 2> x, I U) {
    return {I(0,0), I(0,0)};
}
