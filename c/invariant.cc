#include "invariant.hh"
#include <array>
#include <cstdint>
#include <ppl.hh>
#include <vector>

using namespace std;
#define epsilon 1e-3

void I_approx(vector<IntervalData> Omega) {
    vector<IntervalData> S;
    vector<IntervalData> N;
    vector<IntervalData> E;
    vector<IntervalData> L;

    vector<C_Polyhedron> Omega_p;
    C_Polyhedron Nc;
    vector<C_Polyhedron> Nd;

    while (!L.empty()) {

        auto x = L.back();
        L.pop_back();

        if (!intersects(x.P_over, Omega_p)) {
            N.push_back(x);
        } else if (can_translate_into(x.P_u_over, x.P_over, Nc, Nd)) {
            S.push_back(x);
        } else if (width(x.interval) < epsilon) {
            E.push_back(x);
        } else {
            auto xs = bisect(x);
            L.push_back(get<0>(xs));
            L.push_back(get<1>(xs));
        }
    }
}


pair<IntervalData, IntervalData> bisect(IntervalData x) {
    pair<I, I> tmp;
    tuple<I, I> l;
    tuple<I, I> r;
    I tmp2;
    uint8_t dim;
    if (width(get<0>(x.interval)) >= width(get<1>(x.interval))) {
        tmp = bisect(get<0>(x.interval));
        tmp2 = get<1>(x.interval);
        l = {get<0>(tmp), tmp2};
        r = {get<1>(tmp), tmp2};
    }
    else {
        tmp = bisect(get<1>(x.interval));
        tmp2 = get<0>(x.interval);
        l = {tmp2, get<0>(tmp)};
        r = {tmp2, get<1>(tmp)};
    }

    IntervalData x_l(l);
    x.lchild = &x_l;
    IntervalData x_r(r);
    x.lchild = &x_r;

    return {x_l, x_r};
}

double width(tuple<I, I> interval) {
    return max(width(get<0>(interval)),
               width(get<1>(interval)));
}

IntervalData::IntervalData(tuple<I, I> interval) {
    interval = interval;
    poly = i2p(interval);
    P_u_over = poly;
    P_over = poly;
    lchild = NULL;
    rchild = NULL;
    invariant = true;
    iter = 0;
}

#define mu 0.5
array<double, 2> midpoint(tuple<I, I> interval) {
    return {median(get<0>(interval)),
            median(get<1>(interval))};
}

array<array<int64_t, 2>, 2> A(array<double, 2> x) {
    return {{100, 1}, {1, 100}};
}

array<int64_t, 2> B(array<double, 2> x) {
    return {100*(1-mu)*x[0], 4*(1-mu)*x[1]}
}

tuple<I, I> Phi(tuple<I, I> X, array<double, 2> x) {
    return {I(0,0), I(0,0)}
}

tuple<I, I> Psi(tuple<I, I> X, array<double, 2> x, I u) {
    return {I(0,0), I(0,0)}
}
