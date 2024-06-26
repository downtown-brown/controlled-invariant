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

typedef struct IntervalData {
    nI interval;
    C_Polyhedron poly;
    C_Polyhedron P_u_over;
    C_Polyhedron P_over;
    IntervalData* lchild;
    IntervalData* rchild;
    bool invariant;
    int64_t iter;
    uint32_t checked;
    IntervalData(nI x);
} IntervalData;

list<IntervalData> I_approx(list<IntervalData> Omega);
list<IntervalData> I_accel(list<IntervalData> Omega);

pair<IntervalData, IntervalData> bisect(IntervalData x);

double width(nI interval);
array<double, n> median(nI interval);


C_Polyhedron A(array<double, 2> x, C_Polyhedron P);

C_Polyhedron B(array<double, 2> x, I U);

nI Phi(nI X, array<double, 2> x);

nI Psi(nI X, array<double, 2> x, I u);

void print_points(list<IntervalData> P);
void print_over(list<IntervalData> P);
void print_u_over(list<IntervalData> P);
void fprint_points(list<IntervalData> P, string fname);

vector<vector<C_Polyhedron>> U_approx(list<IntervalData> Omega);

vector<C_Polyhedron> regiondiff(C_Polyhedron P,
                                vector<C_Polyhedron>::iterator curr,
                                vector<C_Polyhedron>::iterator end);
bool subset(C_Polyhedron P,
            vector<C_Polyhedron>::iterator curr,
            vector<C_Polyhedron>::iterator end);
bool can_translate_into(C_Polyhedron P,
                        C_Polyhedron P_over,
                        C_Polyhedron Nc,
                        vector<C_Polyhedron> Nd);
vector<C_Polyhedron> translate_into(C_Polyhedron P,
                                    C_Polyhedron P_over,
                                    C_Polyhedron Nc,
                                    vector<C_Polyhedron> Nd);
void print_points(vector<C_Polyhedron> P);
void fprint_points(vector<C_Polyhedron> P, string fname);
void fprint_points(C_Polyhedron P, string fname, bool append);
void print_points(C_Polyhedron P);
bool intersects(C_Polyhedron A, vector<C_Polyhedron> B);

C_Polyhedron translate_into(C_Polyhedron C, C_Polyhedron N);
vector<C_Polyhedron> translate_into(C_Polyhedron C, vector<C_Polyhedron> N, C_Polyhedron D);
C_Polyhedron translate_touching(C_Polyhedron C, C_Polyhedron N);
vector<C_Polyhedron> translate_touching(C_Polyhedron C, vector<C_Polyhedron> N);

C_Polyhedron operator+(C_Polyhedron a, C_Polyhedron b);

C_Polyhedron i2p(nI x_int);
void rat_approx(double f, int64_t md, int64_t *num, int64_t *denom);


C_Polyhedron convexhull(vector<C_Polyhedron> P_v);

list<IntervalData> merge(IntervalData A, IntervalData B);
void merge(list<IntervalData> &Omega);
