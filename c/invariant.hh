#include "polyhedra.hh"
#include <stdint.h>
#include <ppl.hh>
#include <type_traits>

using namespace std;
using namespace boost::numeric;
using namespace interval_lib;

typedef interval<double, policies<save_state<rounded_transc_std<double>>,
                                  checking_base<double>>> I;


typedef struct IntervalData {
    tuple<I, I> interval;
    C_Polyhedron poly;
    C_Polyhedron P_u_over;
    C_Polyhedron P_over;
    IntervalData* lchild;
    IntervalData* rchild;
    bool invariant;
    int64_t iter;
    IntervalData(tuple<I, I> interval);
} IntervalData;

vector<IntervalData> I_approx(vector<IntervalData> Omega);

pair<IntervalData, IntervalData> bisect(IntervalData x);

double width(tuple<I, I> interval);


C_Polyhedron A(array<double, 2> x, C_Polyhedron P);

C_Polyhedron B(array<double, 2> x, I U);

tuple<I, I> Phi(tuple<I, I> X, array<double, 2> x);

tuple<I, I> Psi(tuple<I, I> X, array<double, 2> x, I u);

void print_points(vector<IntervalData> P);
void print_over(vector<IntervalData> P);
void print_u_over(vector<IntervalData> P);
void fprint_points(vector<IntervalData> P, string fname);

vector<vector<C_Polyhedron>> U_approx(vector<IntervalData> Omega);
