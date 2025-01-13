#include <mutex>
#include <atomic>
#include <thread>
#include <time.h>

#include "set_computations.hh"
#include "print.hh"

void tests(void);
vector<IntervalData> I_approx(const vector<IntervalData>& Omega);

using namespace std;

const double epsilon = 3e-2;
uint8_t stop = 0;

interval_t U(-2, 2);
ninterval_t Omega_0 = {interval_t(-6, 6), interval_t(-6, 6)};

#define DATADIR (string)"data/"

int main() {
    tests();

    auto start = clock();
    vector<IntervalData> res = {Omega_0};
    uint64_t i = 0;
    while (!stop) {
        res = I_approx(res);
        i++;
    }
    auto end = clock();
    printf("Computed %ld iters I_approx %.10f seconds\n\n", i,
           (double)(end - start)/CLOCKS_PER_SEC);

    return 0;
}

uint64_t kk = 0;

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
        } else if (!wider_than(x.interval, epsilon)) {
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

    merge(S);

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

#include "models/jet_engine.hh"

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
}
