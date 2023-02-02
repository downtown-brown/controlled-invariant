

#include "polyhedra.hh"
#include <cstdint>
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

void I_approx(vector<IntervalData> Omega);

pair<IntervalData, IntervalData> bisect(IntervalData x);

double width(tuple<I, I> interval);
