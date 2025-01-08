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
#define epsilon 1e-2

extern uint8_t stop;

uint64_t kk = 0;
interval_t U(-2, 2);

vector<vector<C_Polyhedron>> U_approx(vector<IntervalData> Omega) {
    uint64_t num_int = 0;
    vector<vector<C_Polyhedron>> U;
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

    while (!L.empty()) {

        IntervalData x = L.back();
        L.pop_back();
        num_int++;

        vector<C_Polyhedron> tmp = translate_into(x.P_u_over, x.P_over, Nc, Nd);
        U.push_back(tmp);
    }

    return U;
}


const int NTHREAD = 1;

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
    double curr_width;
    int max_dim = 0;

    for (int i = 0; i < n; i++) {
        curr_width = width(x.interval[i]);
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

IntervalData::IntervalData(ninterval_t x) {
    interval = x;
    nvec_t x_m = median(interval);

    poly = i2p(interval);

    P_u_over = A(x_m, poly) + i2p(Phi(interval, x_m) + Psi(interval, x_m, U));
    P_over = P_u_over + B(x_m, U);

    lchild = NULL;
    rchild = NULL;
    status = STATUS_UNDETERMINED;
    iter = 0;
}

nvec_t median(ninterval_t interval) {
    nvec_t res;

    for (int i = 0; i < n; i++) {
        res[i] = median(interval[i]);
    }

    return res;
}

inline C_Polyhedron A(nvec_t x, C_Polyhedron P) {
    C_Polyhedron res = P;
    const int64_t A1 = 10-1;
    const int64_t A2 = 2;
    const int64_t A3 = -3;
    const int64_t A4 = 10+4;

    const int64_t den = 10;

    res.affine_image(Variable(0), A1*Variable(0) + A2*Variable(1), den);

    res.affine_image(Variable(1), den*A3*Variable(0) + (A1*A4 - A2*A3)*Variable(1), A1*den);
    return res;
}

inline C_Polyhedron B(nvec_t x, interval_t U) {
    const int64_t B1_num = 0;
    const int64_t B2_num = -20;
    const int64_t den = 100;

    static uint8_t i = 0;
    static int64_t nl0, dl0;
    static int64_t nu0, du0;
    if (!i) {
        rat_approx(U.lower(), INT32_MAX, &nl0, &dl0);
        rat_approx(U.upper(), INT32_MAX, &nu0, &du0);
        i = 1;
    }

    C_Polyhedron res(2, EMPTY);
    res.add_generator(point(B1_num*nl0*Variable(0) + B2_num*nl0*Variable(1), dl0*den));
    res.add_generator(point(B1_num*nu0*Variable(0) + B2_num*nu0*Variable(1), du0*den));

    return res;
}

inline ninterval_t Phi(ninterval_t X, nvec_t x_m) {
    interval_t Phi2 = -0.025*pow(X[1],3);
    return {interval_t(0,0), Phi2};
}

inline ninterval_t Psi(ninterval_t X, nvec_t x, interval_t U) {
    return {interval_t(0,0), interval_t(0,0)};
}
