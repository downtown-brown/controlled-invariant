#include <cstdint>
#include <stdint.h>
#include <ppl.hh>
#include <type_traits>
#include <algorithm>
#include <iterator>
#include <ppl.hh>
#include <stdio.h>
#include <time.h>
#include <type_traits>
#include <vector>
#include <boost/numeric/interval.hpp>

#define DATADIR (string)"data/"

const int NDIM = 2;

using namespace std;
using namespace Parma_Polyhedra_Library;

namespace bil = boost::numeric::interval_lib;
using interval_policy = bil::policies<
    bil::save_state<bil::rounded_transc_std<double>>,
    bil::checking_base<double>
    >;

using interval_t = boost::numeric::interval<double, interval_policy>;
using ninterval_t = array<interval_t, NDIM>;
using nvec_t = array<double, NDIM>;

enum InvarianceStatus {
    STATUS_IN,
    STATUS_NOT_IN,
    STATUS_BOUNDARY,
    STATUS_UNDETERMINED
};

typedef struct IntervalData {
    ninterval_t interval;
    C_Polyhedron poly;
    C_Polyhedron P_u_over;
    C_Polyhedron P_over;
    IntervalData* lchild;
    IntervalData* rchild;
    enum InvarianceStatus status;
    IntervalData(ninterval_t x);
} IntervalData;

vector<IntervalData> I_approx(const vector<IntervalData>& Omega);
vector<IntervalData> I_accel(const vector<IntervalData>& Omega);

pair<IntervalData, IntervalData> bisect(IntervalData x);

bool wider_than(ninterval_t interval);
nvec_t median(ninterval_t interval);

void print_points(const vector<IntervalData>& P);
void fprint_points(const vector<IntervalData>& P, string fname);
void print_points(const vector<C_Polyhedron>& P);
void fprint_points(const vector<C_Polyhedron>& P, string fname);
void fprint_points(const C_Polyhedron& P, string fname, bool append);
void print_points(const C_Polyhedron& P, ostream& f=cout);

vector<vector<C_Polyhedron>> U_approx(vector<IntervalData> Omega);

vector<C_Polyhedron> regiondiff(C_Polyhedron P,
                                vector<C_Polyhedron>::iterator curr,
                                vector<C_Polyhedron>::iterator end);
bool subset(C_Polyhedron P,
            vector<C_Polyhedron>::iterator curr,
            vector<C_Polyhedron>::iterator end);
bool can_translate_into(const C_Polyhedron& P,
                        const C_Polyhedron& P_over,
                        const C_Polyhedron& Nc,
                        const vector<C_Polyhedron>& Nd);
vector<C_Polyhedron> translate_into(const C_Polyhedron& P,
                                    const C_Polyhedron& P_over,
                                    const C_Polyhedron& Nc,
                                    const vector<C_Polyhedron>& Nd);
C_Polyhedron translate_into(const C_Polyhedron &C, const C_Polyhedron& N);
vector<C_Polyhedron> translate_into(const C_Polyhedron& C,
                                    const vector<C_Polyhedron>& N,
                                    const C_Polyhedron& D);
C_Polyhedron translate_touching(const C_Polyhedron& C, const C_Polyhedron& N);
vector<C_Polyhedron> translate_touching(const C_Polyhedron& C, const vector<C_Polyhedron>& N);
bool intersects(const C_Polyhedron& A, const vector<C_Polyhedron>& B);
bool intersects(const C_Polyhedron& A, const C_Polyhedron& B);

C_Polyhedron operator+(const C_Polyhedron& a, const C_Polyhedron& b);
ninterval_t operator+(const ninterval_t& A, const ninterval_t& B);

C_Polyhedron i2p(ninterval_t x_int);
int64_t rat_approx(double f, int64_t den);


C_Polyhedron convexhull(const vector<C_Polyhedron>& P_v);
C_Polyhedron intervalhull(const vector<C_Polyhedron>& P_v);

void merge(vector<IntervalData>& Omega);
