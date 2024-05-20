#include <algorithm>
#include <cfloat>
#include <climits>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iterator>
#include <ppl.hh>
#include <stdio.h>
#include <sys/types.h>
#include <time.h>
#include <type_traits>
#include <vector>

#include <boost/numeric/interval.hpp>
#include "polyhedra.hh"


using namespace Parma_Polyhedra_Library;
using namespace std;
using namespace boost::numeric;
using namespace interval_lib;

typedef interval<double, policies<save_state<rounded_transc_std<double>>,
                                  checking_base<double>>> I;

static Variable x(0);
static Variable y(1);

vector<C_Polyhedron> regiondiff(C_Polyhedron P,
                                vector<C_Polyhedron>::iterator curr,
                                vector<C_Polyhedron>::iterator end)
{
    vector<C_Polyhedron> res;

    P = C_Polyhedron(P);

    for (;;) {

        auto tmp = C_Polyhedron(P);
        tmp.intersection_assign(*curr);

        if (tmp.affine_dimension() >= P.affine_dimension()) {
            break;
        }

        ++curr;

        if (curr == end) {
            res.push_back(P);
            return res;
        }
    }

    auto q_h = curr->constraints();

    for (auto i = q_h.begin(); i != q_h.end(); ++i) {
        auto tmp = C_Polyhedron(P);

        tmp.add_constraint(Constraint(-i->coefficient(x)*x
                                      - i->coefficient(y)*y
                                      - i->inhomogeneous_term() >= 0)); // not done

        if (tmp.affine_dimension() < P.affine_dimension()) {
            continue;
        }

        if (curr < end - 1) {
            auto tmp2 = regiondiff(tmp, curr + 1, end);
            res.insert(res.end(), tmp2.begin(), tmp2.end());
        }
        else {
            res.push_back(tmp);
        }

        P.add_constraint(*i);
    }

    return res;
}

bool subset(C_Polyhedron P,
            vector<C_Polyhedron>::iterator curr,
            vector<C_Polyhedron>::iterator end)
{
    if (P.is_empty()) {
        return true;
    }

    P = C_Polyhedron(P);

    for (;;) {

        auto tmp = C_Polyhedron(P);
        tmp.intersection_assign(*curr);
        if (tmp.affine_dimension() >= P.affine_dimension()) {
            //if (!tmp.is_empty()) {
            break;
        }

        ++curr;

        if (curr == end) {
            return false;
        }
    }

    auto q_h = curr->constraints();

    for (auto i = q_h.begin(); i != q_h.end(); ++i) {
        auto tmp = C_Polyhedron(P);

        tmp.add_constraint(Constraint(-i->coefficient(x)*x
                                      - i->coefficient(y)*y
                                      - i->inhomogeneous_term() >= 0));

        if (tmp.affine_dimension() < P.affine_dimension()) {
            continue;
        }

        if (curr == end - 1) {
            return false;
        } else if (!subset(tmp, curr + 1, end)) {
                return false;
        }

        P.add_constraint(*i);
    }

    return true;
}

C_Polyhedron translate_into(C_Polyhedron C, C_Polyhedron N) {
    C_Polyhedron res(2);
    auto N_h = N.constraints();
    auto C_p = C.generators();
    for (auto n = N_h.begin(); n != N_h.end(); ++n) {
        GMP_Integer min = 0;
        GMP_Integer mind = 1;
        GMP_Integer d;
        for (auto c = C_p.begin(); c != C_p.end(); ++c) {
            GMP_Integer tmp = n->coefficient(x)*c->coefficient(x)
                + n->coefficient(y)*c->coefficient(y);

            d = c->divisor();

            if (tmp*mind < min*d || min == 0) {
                min = tmp;
                mind = d;
            }
        }
        res.add_constraint(n->coefficient(x)*mind*x
                           + n->coefficient(y)*mind*y
                           + mind*n->inhomogeneous_term()
                           + min >= 0);
    }

    return res;
}

vector<C_Polyhedron> translate_into(C_Polyhedron C, vector<C_Polyhedron> N, C_Polyhedron D) {
    vector<C_Polyhedron> res;

    auto U1 = translate_into(C, D);
    auto U2 = translate_touching(C, N);

    return regiondiff(U1, U2.begin(), U2.end());
}

C_Polyhedron translate_touching(C_Polyhedron C, C_Polyhedron N) {
    C_Polyhedron res(2);
    auto N_h = N.constraints();
    auto C_p = C.generators();
    for (auto n = N_h.begin(); n != N_h.end(); ++n) {
        GMP_Integer max = 0;
        GMP_Integer maxd = 1;
        GMP_Integer d;
        for (auto c = C_p.begin(); c != C_p.end(); ++c) {
            GMP_Integer tmp = n->coefficient(x)*c->coefficient(x)
                + n->coefficient(y)*c->coefficient(y);

            d = c->divisor();

            if (tmp*maxd > max*d || max == 0) {
                max = tmp;
                maxd = d;
            }
        }
        res.add_constraint(n->coefficient(x)*maxd*x
                           + n->coefficient(y)*maxd*y
                           + maxd*n->inhomogeneous_term()
                           + max >= 0);
    }

    auto C_h = C.constraints();
    auto N_p = N.generators();
    for (auto c = C_h.begin(); c != C_h.end(); ++c) {
        GMP_Integer max = 0;
        GMP_Integer maxd = 1;
        GMP_Integer d;
        for (auto n = N_p.begin(); n != N_p.end(); ++n) {
            GMP_Integer tmp = c->coefficient(x)*n->coefficient(x)
                + c->coefficient(y)*n->coefficient(y);

            d = n->divisor();

            if (tmp*maxd > max*d || max == 0) {
                max = tmp;
                maxd = d;
            }
        }
        res.add_constraint(-c->coefficient(x).get_d()*maxd*x
                           - c->coefficient(y).get_d()*maxd*y
                           + maxd*c->inhomogeneous_term().get_d()
                           + max >= 0);
    }

    return res;
}

vector<C_Polyhedron> translate_touching(C_Polyhedron C, vector<C_Polyhedron> N) {
    vector<C_Polyhedron> res;

    for (auto n = N.begin(); n != N.end(); ++n) {
        res.emplace_back(translate_touching(C, *n));
    }

    return res;
}

C_Polyhedron operator+(C_Polyhedron a, C_Polyhedron b) {
    C_Polyhedron res(2, EMPTY);
    auto V_a = a.generators();
    auto V_b = b.generators();
    for (auto v_a = V_a.begin(); v_a != V_a.end(); ++v_a) {
        for (auto v_b = V_b.begin(); v_b != V_b.end(); ++v_b) {
            res.add_generator(point((v_a->coefficient(x)*v_b->divisor()
                                     + v_b->coefficient(x)*v_a->divisor())*x
                                    + (v_a->coefficient(y)*v_b->divisor()
                                       + v_b->coefficient(y)*v_a->divisor())*y,
                                    v_a->divisor()*v_b->divisor()));

        }
    }
    return C_Polyhedron(res.minimized_generators());
}

void print_points(C_Polyhedron P) {
    auto gs = P.generators();
    cout << "[";
    for (auto g = gs.begin(); g != gs.end(); ++g) {
        if (g->is_point()) {
            cout << "[" << g->coefficient(x).get_d() / g->divisor().get_d() << ","
                 << g->coefficient(y).get_d() / g->divisor().get_d() << "]; ";
        }
        else {
            g->ascii_dump();
        }
    }
    cout << "\b\b]\n\n";
}

void fprint_points(C_Polyhedron P, string fname) {
    auto gs = P.generators();
    ofstream f(fname);
    f << "[";
    for (auto g = gs.begin(); g != gs.end(); ++g) {
        if (g->is_point()) {
            f << "[" << g->coefficient(x).get_d() / g->divisor().get_d() << ","
                 << g->coefficient(y).get_d() / g->divisor().get_d() << "]; ";
        }
        else {
            g->ascii_dump();
        }
    }
    f << "]\n";
}

void print_points(vector<C_Polyhedron> P) {
    for (auto p = P.begin(); p != P.end(); ++p) {
        print_points(*p);
    }
}
void fprint_points(vector<C_Polyhedron> P, string fname) {
    for (auto p = P.begin(); p != P.end(); ++p) {
        fprint_points(*p, fname);
    }
}

bool intersects(C_Polyhedron A, vector<C_Polyhedron> B) {
    for (auto b = B.begin(); b != B.end(); ++b) {
        auto tmp = C_Polyhedron(A);
        tmp.intersection_assign(*b);

        if (!tmp.is_empty()) {
            return true;
        }
    }

    return false;
}

vector<C_Polyhedron> translate_into(C_Polyhedron P,
                                    C_Polyhedron P_over,
                                    C_Polyhedron Nc,
                                    vector<C_Polyhedron> Nd) {
    auto Ncc = C_Polyhedron(Nc);
    Ncc.intersection_assign(P_over);

    if (!Nd.empty()) {
        auto U1 = translate_into(P, Ncc);

        vector<C_Polyhedron> Ndd;
        for (auto d = Nd.begin(); d != Nd.end(); ++d) {
            auto tmp = C_Polyhedron(*d);
            tmp.intersection_assign(P_over);

            if (!tmp.is_empty()) {
                Ndd.push_back(tmp);
            }
        }

        auto U2 = translate_touching(P, Ndd);

        if (U2.empty()) {
            return {U1};
        } else {
            //cout << "U1:\n";
            //print_points(U1);
            //
            //cout << "U2:\n";
            //print_points(U2);
            //cout << "subset: " << subset(U1, U2.begin(), U2.end());
            //

            return regiondiff(U1, U2.begin(), U2.end());
        }

    } else {
        auto U1 = translate_into(P, Ncc);
        return {U1};
    }
}

bool can_translate_into(C_Polyhedron P,
                        C_Polyhedron P_over,
                        C_Polyhedron Nc,
                        vector<C_Polyhedron> Nd) {
    auto Ncc = C_Polyhedron(Nc);
    Ncc.intersection_assign(P_over);

    if (!Nd.empty()) {
        auto U1 = translate_into(P, Ncc);

        vector<C_Polyhedron> Ndd;
        for (auto d = Nd.begin(); d != Nd.end(); ++d) {
            auto tmp = C_Polyhedron(*d);
            tmp.intersection_assign(P_over);

            if (!tmp.is_empty()) {
                Ndd.push_back(tmp);
            }
        }

        auto U2 = translate_touching(P, Ndd);

        if (U2.empty()) {
            return !U1.is_empty();
        } else {
            //cout << "U1:\n";
            //print_points(U1);
            //
            //cout << "U2:\n";
            //print_points(U2);
            //cout << "subset: " << subset(U1, U2.begin(), U2.end());
            //

            return !subset(U1, U2.begin(), U2.end());
        }

    } else {
        auto U1 = translate_into(P, Ncc);
        return !U1.is_empty();
    }
}

C_Polyhedron i2p(tuple<I, I> x_int) {

    int64_t nl0, dl0;
    rat_approx(get<0>(x_int).lower(), INT16_MAX, &nl0, &dl0);
    int64_t nl1, dl1;
    rat_approx(get<1>(x_int).lower(), INT16_MAX, &nl1, &dl1);
    int64_t nu0, du0;
    rat_approx(get<0>(x_int).upper(), INT16_MAX, &nu0, &du0);
    int64_t nu1, du1;
    rat_approx(get<1>(x_int).upper(), INT16_MAX, &nu1, &du1);

    C_Polyhedron res(2, EMPTY);
    res.add_generator(point((nl0*dl1)*x + (nl1*dl0)*y, dl0*dl1));
    res.add_generator(point((nl0*du1)*x + (nu1*dl0)*y, dl0*du1));
    res.add_generator(point((nu0*dl1)*x + (nl1*du0)*y, du0*dl1));
    res.add_generator(point((nu0*du1)*x + (nu1*du0)*y, du0*du1));

    return res;
}

void rat_approx(double f, int64_t md, int64_t *num, int64_t *denom)
{
	/*  a: continued fraction coefficients. */
	int64_t a, h[3] = { 0, 1, 0 }, k[3] = { 1, 0, 0 };
	int64_t x, d, n = 1;
	int i, neg = 0;

	if (md <= 1) { *denom = 1; *num = (int64_t) f; return; }

	if (f < 0) { neg = 1; f = -f; }

	while (f != floor(f) && n < 4611686018427387904) { n <<= 1; f *= 2;}
	d = f;

	/* continued fraction and check denominator each step */
	for (i = 0; i < 64; i++) {
		a = n ? d / n : 0;
		if (i && !a) break;

		x = d; d = n; n = x % n;

		x = a;
		if (k[1] * a + k[0] >= md) {
			x = (md - k[0]) / k[1];
			if (x * 2 >= a || k[1] >= md)
				i = 65;
			else
				break;
		}

		h[2] = x * h[1] + h[0]; h[0] = h[1]; h[1] = h[2];
		k[2] = x * k[1] + k[0]; k[0] = k[1]; k[1] = k[2];
	}
	*denom = k[1];
	*num = neg ? -h[1] : h[1];
}

C_Polyhedron convexhull(vector<C_Polyhedron> P_v) {
    C_Polyhedron res(2, EMPTY);

    for (auto P = P_v.begin(); P != P_v.end(); ++P) {
        res.add_generators(P->generators());
    }

    return C_Polyhedron(res.minimized_generators());
}
