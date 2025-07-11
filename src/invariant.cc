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
//include "models/cartpole.hh"
#include "models/pendubot.hh"
//#include "models/cstr.hh"
//#include "models/pendulum_CDC24.hh"
//#include "models/robot_exploration.hh"

IntervalData::IntervalData(ninterval_t x) {
    interval = x;
    nvec_t x_m = median(interval);

    C_Polyhedron poly = i2p(interval);

    P_u_over = A(x_m, poly) + i2p(Phi(interval, x_m) + Psi(interval, x_m) + Delta);
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

static vector<ninterval_t> N;
static vector<ninterval_t> S;

static atomic<uint64_t> n_waiting;
static atomic i_total = 0;

void I_worker(vector<IntervalData> &L, vector<IntervalData> &L_next,
              const C_Polyhedron& Nc, const vector<C_Polyhedron>& Nd, int t) {
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
        i_total++;

        if (!intersects(x.P_over, Nc)) {
            num_N++;
            std::scoped_lock lock(N_mutex);
            N.push_back(x.interval);
        } else if (can_translate_into(x.P_u_over, x.P_over, Nc, Nd)) {
            std::scoped_lock lock(S_mutex);
            S.push_back(x.interval);
            L_next.push_back(x);
        } else if (!wider_than(x.interval, epsilon)) {
            num_E++;
            std::scoped_lock lock(S_mutex);
            S.push_back(x.interval);
        } else {
            num_B++;
            pair<IntervalData, IntervalData> xs = bisect(x);
            std::scoped_lock lock(L_mutex);
            L.push_back(get<0>(xs));
            L.push_back(get<1>(xs));
        }

        if (num_int % 1000 == 0) {
            std::scoped_lock lock(print_mutex);
            cout << "num_int: " << setw(20) << num_int << '\r' << flush;
        }
    }
}

vector<IntervalData> I_approx(const vector<IntervalData>& Omega) {
    static uint64_t i_approx_iter = 0;
    vector<IntervalData> L = Omega;
    vector<IntervalData> L_next;

    i_approx_iter++;
    num_int = num_N = num_E = num_B = 0;

    ninterval_t hull = intervalhull(S);
    C_Polyhedron Nc = i2p(hull);

    vector<C_Polyhedron> Nd;

    for (const ninterval_t& x : N) {
        if (intersects(x, hull)) {
            Nd.emplace_back(i2p(x));
        }
    }

    S.clear();

    n_waiting = 0;
    vector<thread> threads;
    for (uint64_t t = 0; t < NTHREAD; t++) {
        threads.emplace_back(make_threadable(I_worker), ref(L), ref(L_next), Nc, Nd, t);
    }

    for (uint64_t t = 0; t < NTHREAD; t++) {
        threads[t].join();
    }

    merge(S);
    merge(N);

    fprint_points(S, DATA_DIR + "s" + to_string(i_approx_iter) + ".txt");

    cout << "Iteration " << i_approx_iter << ", considered " << num_int << " intervals." << endl;
    cout << "N: " << num_N << ", S: " << S.size() << ", E: " << num_E << ", B: " << num_B << endl;

    return L_next;
}

int main() {
    // tests();
    vector<IntervalData> res = {Omega_0};

    S = {Omega_0};
    do {
        res = I_approx(res);
    } while (num_N || num_E);

    cout << "Computed invariant set." << endl;
    cout << i_total << " total intervals" << endl;

    return 0;
}
