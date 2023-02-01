

#include "polyhedra.hh"
#include <ppl.hh>
#include <type_traits>

struct ProblemData {
    A,
    Phi,
    B,
    Psi,
    U,
    double epsilon
};


void I_approx(vector<Interval> Omega, ProblemData p) {
    C_Polyhedron c1(2, EMPTY);
    c1.affine_image(

    vector<double> L;
    while (!L.empty()) {

        auto x = L.back();
        L.pop_back();

        auto x_p = i2p(x);

        if (!intersects(image_over, Omega_p)) {
            N.push_back(x);
        } else if (can_translate_into(image_poly, image_over)) {
            S.push_back(x);
        } else if (diam(x) < p.epsilon) {
            E.push_back(x);
        } else {
            auto xs = bisect(x);
            L.push_back(xs[0]);
            L.push_back(xs[1]);
        }
    }
}
