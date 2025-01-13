#include <mutex>
#include <atomic>
#include <thread>
#include <time.h>

#include "set_computations.hh"
#include "print.hh"

void tests(void);

static const string DATA_DIR = "data/";
static const double epsilon = 3e-2;

static interval_t U(-2, 2);
static ninterval_t Omega_0 = {interval_t(-6, 6),
                              interval_t(-6, 6)};

#include "models/jet_engine.hh"

IntervalData::IntervalData(ninterval_t x) {
    interval = x;
    nvec_t x_m = median(interval);

    poly = i2p(interval);

    static C_Polyhedron BU = B(x_m, U);
    P_u_over = A(x_m, poly) + i2p(Phi(interval, x_m));
    P_over = P_u_over + BU;

    status = STATUS_UNDETERMINED;
}

static const int NTHREAD = 24;

static mutex L_mutex;
static mutex S_mutex;
static atomic<uint64_t> num_int = 0;

static atomic<uint64_t> num_N = 0;
static atomic<uint64_t> num_E = 0;
static atomic<uint64_t> num_B = 0;

void I_worker(vector<IntervalData>& L,
              vector<IntervalData>& S,
              C_Polyhedron Nc,
              vector<C_Polyhedron> Nd,
              int t) {
    int num_t = 0;
    while (1) {
        L_mutex.lock();
        if (L.empty()) {
            L_mutex.unlock();
            printf("thread %d considered %d intervals\n", t, num_t);
            return;
        }
        IntervalData x = L.back();
        L.pop_back();
        L_mutex.unlock();

        num_int++;
        num_t++;

        if (!intersects(x.P_over, Nc)) {
            x.status = STATUS_NOT_IN;
            num_N++;
        } else if (can_translate_into(x.P_u_over, x.P_over, Nc, Nd)) {
            x.status = STATUS_IN;
            S_mutex.lock();
            S.push_back(x);
            S_mutex.unlock();
        } else if (!wider_than(x.interval, epsilon)) {
            x.status = STATUS_BOUNDARY;
            num_E++;
        } else {
            num_B++;
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
    static uint64_t i_approx_iter = 0;
    vector<IntervalData> L = Omega;
    vector<IntervalData> S;

    i_approx_iter++;
    num_int = num_N = num_E = num_B = 0;

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
        thread_vec.emplace_back(make_threadable(I_worker), ref(L), ref(S), Nc, Nd, t);
    }

    for (int t = 0; t < NTHREAD; t++) {
        thread_vec[t].join();
    }

    merge(S);
    fprint_points(S, DATA_DIR + "s" + to_string(i_approx_iter) + ".txt");

    cout << "Iteration " << i_approx_iter << ", considered " << num_int << " intervals." << endl;
    cout << "N: " << num_N << ", S: " << S.size() << ", E: " << num_E << ", B: " << num_B << endl;

    return S;
}

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

        Uc.push_back(translate_into(x.P_u_over, x.P_over, Nc, Nd));
    }

    return Uc;
}

int main() {
    //tests();

    vector<IntervalData> res = {Omega_0};

    do {
        res = I_approx(res);
    } while (num_N || num_E);

    cout << "Computed invariant set.\n\n";

    return 0;
}
