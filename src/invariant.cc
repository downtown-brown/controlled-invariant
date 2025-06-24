#include <cstdint>
#include <mutex>
#include <atomic>
#include <ppl.hh>
#include <thread>
#include <time.h>
#include <chrono>
#include <thread>
#include <iomanip>

#include "set_computations.hh"
#include "print.hh"

void tests(void);

//#include "models/artificial_system.hh"
//#include "models/jet_engine.hh"
//#include "models/cart.hh"
//#include "models/mass_spring_damper.hh"
//#include "models/cartpole_pendulum.hh"
//#include "models/van_der_pol.hh"
//#include "models/cartpole.hh"
//#include "models/pendubot.hh"
//#include "models/pendulum_CDC24.hh"
//#include "models/robot_exploration.hh"
#include "models/dist_mass.hh"

static const int n = 3;

vector<ninterval_t> delta(n);

//ninterval_t delta_2 = Omega_0;
//ninterval_t delta_nm1 = Omega_0;

IntervalData::IntervalData(ninterval_t x, int node) : node(node) {
    interval = x;
    nvec_t x_m = median(interval);

    poly = i2p(interval);

    if (node == 0) {
        P_u_over = A(x_m, poly) + i2p(Phi_1(interval, x_m, delta[1]));
    } else if (node == n - 1) {
        P_u_over = A(x_m, poly) + i2p(Phi_n(interval, x_m, delta[n-2]));
    } else {
      P_u_over = A(x_m, poly) + i2p(Phi_i(interval, x_m,
                                          delta[node-1], delta[node+1]));
    }

    P_over = P_u_over + B(x_m, U);
}

static const unsigned int NTHREAD = 24;

static mutex L_mutex;
static mutex S_mutex;
static mutex N_mutex;
static mutex print_mutex;
static atomic<uint64_t> num_int = 0;

static atomic<uint64_t> num_N = 0;
static atomic<uint64_t> num_E = 0;
static atomic<uint64_t> num_B = 0;

static vector<vector<ninterval_t>> N(n);
static vector<vector<ninterval_t>> E(n);
static vector<vector<ninterval_t>> S(n);

static atomic<uint64_t> n_waiting;

void I_worker(vector<IntervalData> &L, vector<IntervalData> &L_next,
              const C_Polyhedron& Nc, const vector<C_Polyhedron>& Nd, int t, int node) {
    bool waiting = false;
    while (n_waiting < NTHREAD) {
        L_mutex.lock();
        if (L.empty()) {
            L_mutex.unlock();
            if (!waiting) {
                waiting = true;
                n_waiting++;
            }
            continue;
        }

        if (waiting) {
            waiting = false;
            n_waiting--;
        }

        IntervalData x = L.back();
        L.pop_back();
        L_mutex.unlock();

        num_int++;

        if (!intersects(x.P_over, Nc)) {
            num_N++;
            N_mutex.lock();
            N[node].push_back(x.interval);
            N_mutex.unlock();
        } else if (can_translate_into(x.P_u_over, x.P_over, Nc, Nd)) {
            S_mutex.lock();
            S[node].push_back(x.interval);
            L_next.push_back(x);
            S_mutex.unlock();
        } else if (!wider_than(x.interval, epsilon)) {
            num_E++;
            S_mutex.lock();
            E[node].push_back(x.interval);
            L_next.push_back(x);
            S_mutex.unlock();
        } else {
            num_B++;
            pair<IntervalData, IntervalData> xs = bisect(x);
            L_mutex.lock();
            L.push_back(get<0>(xs));
            L.push_back(get<1>(xs));
            L_mutex.unlock();
        }

        if (num_int % 1000 == 0) {
            print_mutex.lock();
            cout << "num_int: " << setw(20) << num_int << '\r' << flush;
            print_mutex.unlock();
        }
    }
}

vector<IntervalData> I_approx(const vector<IntervalData>& Omega, int node) {
    static uint64_t i_approx_iter = 0;
    vector<IntervalData> L = Omega;
    vector<IntervalData> L_next;

    i_approx_iter++;
    num_int = num_N = num_E = num_B = 0;

    S[node].insert(S[node].end(), E[node].begin(), E[node].end());

    ninterval_t hull = intervalhull(S[node]);
    C_Polyhedron Nc = i2p(hull);

    vector<C_Polyhedron> Nd;

    for (const ninterval_t& x : N[node]) {
        if (intersects(x, hull)) {
            Nd.emplace_back(i2p(x));
        }
    }

    S[node].clear();

    n_waiting = 0;
    vector<thread> threads;
    for (uint64_t t = 0; t < NTHREAD; t++) {
        threads.emplace_back(make_threadable(I_worker), ref(L), ref(L_next), Nc, Nd, t, node);
    }

    for (uint64_t t = 0; t < NTHREAD; t++) {
        threads[t].join();
    }

    merge(S[node]);
    merge(N[node]);

    fprint_points(S[node], DATA_DIR + "s" + to_string(i_approx_iter) + to_string(node) + ".txt");

    cout << "Iteration " << i_approx_iter << ", considered " << num_int << " intervals." << endl;
    cout << "N: " << num_N << ", S: " << S[node].size() << ", E: " << num_E << ", B: " << num_B << endl;

    return L_next;
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
    // tests();
    vector<vector<IntervalData>> res(n);

    for (int i = 0; i < n; i++) {
        res[i].emplace_back(Omega_0, i);
        S[i].push_back(Omega_0);
    }

    do {
        for (int i = 0; i < n; i++) {
            res[i] = I_approx(res[i], i);
            delta[i] = intervalhull(S[i]);
        }

        //delta_2 = intervalhull(S[1]);
        //delta_nm1 = intervalhull(S[0]);
    } while (num_N || num_E);

    cout << "Computed invariant set.\n\n";

    return 0;
}
