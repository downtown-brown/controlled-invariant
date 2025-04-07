#include <ppl.hh>
#include <boost/numeric/interval.hpp>
#include <vector>

#ifndef SET_TYPES_HH
#define SET_TYPES_HH

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

#endif
