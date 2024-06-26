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

using namespace std;
using namespace boost::numeric;
using namespace interval_lib;
using namespace Parma_Polyhedra_Library;

constexpr int n = 2;

typedef interval<double, policies<save_state<rounded_transc_std<double>>,
                                  checking_base<double>>> I;

typedef array<I, n> nI;

enum InvarianceStatus {
    STATUS_IN,
    STATUS_NOT_IN,
    STATUS_BOUNDARY,
    STATUS_UNDETERMINED
};

typedef struct IntervalData {
    nI interval;
    C_Polyhedron poly;
    C_Polyhedron P_u_over;
    C_Polyhedron P_over;
    IntervalData* lchild;
    IntervalData* rchild;
    enum InvarianceStatus status;
    int64_t iter;
    uint32_t checked;
    IntervalData(nI x);
} IntervalData;

list<IntervalData> I_approx(const list<IntervalData>& Omega);
list<IntervalData> I_accel(const list<IntervalData>& Omega);

pair<IntervalData, IntervalData> bisect(IntervalData x);

bool wider_than(nI interval);
array<double, n> median(nI interval);


C_Polyhedron A(array<double, 2> x, C_Polyhedron P);

C_Polyhedron B(array<double, 2> x, I U);

nI Phi(nI X, array<double, 2> x);

nI Psi(nI X, array<double, 2> x, I u);

void print_points(const list<IntervalData>& P);
void print_over(const list<IntervalData>& P);
void print_u_over(const list<IntervalData>& P);
void fprint_points(const list<IntervalData>& P, string fname);
void print_points(const vector<C_Polyhedron>& P);
void fprint_points(const vector<C_Polyhedron>& P, string fname);
void fprint_points(const C_Polyhedron& P, string fname, bool append);
void print_points(const C_Polyhedron& P);

vector<vector<C_Polyhedron>> U_approx(list<IntervalData> Omega);

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

C_Polyhedron operator+(const C_Polyhedron& a, const C_Polyhedron& b);

C_Polyhedron i2p(nI x_int);
void rat_approx(double f, int64_t md, int64_t *num, int64_t *denom);


C_Polyhedron convexhull(const vector<C_Polyhedron>& P_v);
C_Polyhedron intervalhull(const vector<C_Polyhedron>& P_v);

list<IntervalData> merge(const IntervalData& A, const IntervalData& B);
void merge(list<IntervalData>& Omega);
