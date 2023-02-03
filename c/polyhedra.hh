#include <algorithm>
#include <iterator>
#include <ppl.hh>
#include <stdio.h>
#include <time.h>
#include <type_traits>
#include <vector>
#include <boost/numeric/interval.hpp>

using namespace Parma_Polyhedra_Library;
using namespace std;

using namespace boost::numeric;
using namespace interval_lib;

typedef interval<double, policies<save_state<rounded_transc_std<double>>,
                                  checking_base<double>>> I;

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
void print_points(vector<C_Polyhedron> P);
void print_points(C_Polyhedron P);
bool intersects(C_Polyhedron A, vector<C_Polyhedron> B);

C_Polyhedron translate_into(C_Polyhedron C, C_Polyhedron N);
vector<C_Polyhedron> translate_into(C_Polyhedron C, vector<C_Polyhedron> N, C_Polyhedron D);
C_Polyhedron translate_touching(C_Polyhedron C, C_Polyhedron N);
vector<C_Polyhedron> translate_touching(C_Polyhedron C, vector<C_Polyhedron> N);

C_Polyhedron operator+(C_Polyhedron a, C_Polyhedron b);

C_Polyhedron i2p(tuple<I, I> x_int);
void rat_approx(double f, int64_t md, int64_t *num, int64_t *denom);


C_Polyhedron convexhull(vector<C_Polyhedron> P_v);
