#include "invariant.hh"
#include <array>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <ppl.hh>
#include <vector>
#include <mutex>
#include <atomic>
#include <thread>

using namespace std;
#define epsilon 3e-2

extern uint8_t stop;

uint64_t kk = 0;
interval_t U(-2, 2);

vector<vector<C_Polyhedron>> U_approx(vector<IntervalData> Omega) {
    vector<vector<C_Polyhedron>> Uc;
    vector<IntervalData> L = Omega;

    vector<C_Polyhedron> Omega_p(Omega.size());
    int i = 0;
    for (const IntervalData& O : Omega) {
        Omega_p[i] = O.poly;
        i++;
    }

    C_Polyhedron Nc = convexhull(Omega_p);
    vector<C_Polyhedron> Nd = regiondiff(Nc, Omega_p.begin(), Omega_p.end());

    while (!L.empty()) {

        IntervalData x = L.back();
        L.pop_back();

        vector<C_Polyhedron> tmp = translate_into(x.P_u_over, x.P_over, Nc, Nd);
        Uc.push_back(tmp);
    }

    return Uc;
}


const int NTHREAD = 12;

atomic<bool> running[NTHREAD];
mutex L_mutex;
mutex S_mutex;
mutex N_mutex;
mutex E_mutex;
atomic<int> bcount = 0;
atomic<uint64_t> num_int = 0;
void I_worker(vector<IntervalData>& L,
              vector<IntervalData>& S,
              vector<IntervalData>& N,
              vector<IntervalData>& E,
              C_Polyhedron Nc,
              vector<C_Polyhedron> Nd,
              int t) {
    int num_t = 0;
    while (1) {

        L_mutex.lock();
        if (L.empty()) {
            printf("thread %d considered %d intervals\n", t, num_t);
            L_mutex.unlock();
            return;
        }
        IntervalData x = L.back();
        L.pop_back();
        L_mutex.unlock();

        num_int++;
        num_t++;

        if (!intersects(x.P_over, Nc)) {
            x.status = STATUS_NOT_IN;
            N_mutex.lock();
            N.push_back(x);
            N_mutex.unlock();
        } else if (can_translate_into(x.P_u_over, x.P_over, Nc, Nd)) {
            x.status = STATUS_IN;
            S_mutex.lock();
            S.push_back(x);
            S_mutex.unlock();
        } else if (!wider_than(x.interval)) {
            x.status = STATUS_BOUNDARY;
            E_mutex.lock();
            E.push_back(x);
            E_mutex.unlock();
        } else {
            bcount++;
            x.status = STATUS_UNDETERMINED;
            pair<IntervalData, IntervalData> xs = bisect(x);
            L_mutex.lock();
            L.push_back(get<0>(xs));
            L.push_back(get<1>(xs));
            L_mutex.unlock();
        }
    }
}

vector<IntervalData> I_approx(const vector<IntervalData>& Omega) {
    bcount = 0;
    num_int = 0;
    vector<IntervalData> S;
    vector<IntervalData> N;
    vector<IntervalData> E;
    vector<IntervalData> L = Omega;

    vector<C_Polyhedron> Omega_p(Omega.size());
    int i = 0;
    for (const IntervalData& O : Omega) {
        Omega_p[i] = O.poly;
        i++;
    }
    C_Polyhedron Nc = convexhull(Omega_p);
    vector<C_Polyhedron> Nd = regiondiff(Nc, Omega_p.begin(), Omega_p.end());

    vector<thread> thread_vec;
    for (int t = 0; t < NTHREAD; t++) {
        thread_vec.emplace_back(make_threadable(I_worker), ref(L), ref(S), ref(N), ref(E), Nc, Nd, t);
    }

    for (int t = 0; t < NTHREAD; t++) {
        thread_vec[t].join();
    }

    //merge(N);
    merge(S);
    //merge(E);

    cout << "N: " << N.size() << ", S: " << S.size() << ", E: " << E.size() << ", B: " << bcount
         << ", Total: " << N.size() + S.size() + E.size() + bcount << endl;

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
    int max_dim = 0;

    for (int i = 0; i < n; i++) {
        double curr_width = width(x.interval[i]);
        if (curr_width > max_width) {
            max_width = curr_width;
            max_dim = i;
        }
    }

    ninterval_t l = x.interval;
    ninterval_t r = x.interval;

    pair<interval_t, interval_t> tmp = bisect(x.interval[max_dim]);

    l[max_dim] = get<0>(tmp);
    r[max_dim] = get<1>(tmp);

    IntervalData x_l(l);
    x.lchild = &x_l;
    IntervalData x_r(r);
    x.rchild = &x_r;

    return {x_l, x_r};
}

bool wider_than(ninterval_t interval) {
    for (int i = 0; i < n; i++) {
        if (width(interval[i]) > epsilon) {
            return true;
        }
    }

    return false;
}

nvec_t median(ninterval_t x) {
    nvec_t res;

    for (int i = 0; i < n; i++) {
        res[i] = median(x[i]);
    }

    return res;
}

inline C_Polyhedron A(nvec_t x, C_Polyhedron P) {
    const int64_t A[n][n] = {{9, 2}, {-3, 14}};
    const int64_t A_den = 10;

    C_Polyhedron res = P;
    res.affine_image(Variable(0), A[0][0]*Variable(0) + A[0][1]*Variable(1), A_den);
    res.affine_image(Variable(1), A_den*A[1][0]*Variable(0)
                     + (A[0][0]*A[1][1] - A[0][1]*A[1][0])*Variable(1), A[0][0]*A_den);
    return res;
}

inline C_Polyhedron B(nvec_t x, interval_t U) {
    const int64_t B[n] = {5, -20};
    const int64_t B_den = 100;

    const int64_t U_den = INT16_MAX;
    static int64_t Ul = rat_approx(U.lower(), U_den);
    static int64_t Uh = rat_approx(U.upper(), U_den);

    C_Polyhedron res(2, EMPTY);
    res.add_generator(point(B[0]*Ul*Variable(0) + B[1]*Ul*Variable(1), B_den*U_den));
    res.add_generator(point(B[0]*Uh*Variable(0) + B[1]*Uh*Variable(1), B_den*U_den));

    return res;
}

inline ninterval_t Phi(ninterval_t x, nvec_t x_m) {
    return {interval_t(0,0), -0.025*pow(x[1],3)};
}

inline ninterval_t Psi(ninterval_t x, nvec_t x_m, interval_t U) {
    return {interval_t(0,0), interval_t(0,0)};
}

IntervalData::IntervalData(ninterval_t x) {
    interval = x;
    nvec_t x_m = median(interval);

    poly = i2p(interval);

    static C_Polyhedron BU = B(x_m, U);
    P_u_over = A(x_m, poly) + i2p(Phi(interval, x_m));
    P_over = P_u_over + BU;

    lchild = NULL;
    rchild = NULL;
    status = STATUS_UNDETERMINED;
    iter = 0;
}
