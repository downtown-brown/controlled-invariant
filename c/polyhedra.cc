#include <algorithm>
#include <climits>
#include <cstdint>
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

        if (tmp.affine_dimension() == tmp.space_dimension()) {
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

        if (tmp.affine_dimension() != tmp.space_dimension()) {
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
            return false;
        }
    }

    auto q_h = curr->constraints();

    for (auto i = q_h.begin(); i != q_h.end(); ++i) {
        auto tmp = C_Polyhedron(P);

        auto t1 = tmp.affine_dimension();
        tmp.add_constraint(Constraint(-i->coefficient(x)*x
                                      - i->coefficient(y)*y
                                      - i->inhomogeneous_term() >= 0)); // not done

        if (tmp.affine_dimension() <= t1) {
            continue;
        }

        if ((curr == end - 1) || subset(tmp, curr + 1, end)) {
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
        int min = INT_MAX;
        long mind;
        long d;
        for (auto c = C_p.begin(); c != C_p.end(); ++c) {
            int tmp = n->coefficient(x).get_si()*c->coefficient(x).get_si()
                + n->coefficient(y).get_si()*c->coefficient(y).get_si();

            d = c->divisor().get_si();

            if (min == INT_MAX || tmp*mind < min*d) {
                min = tmp;
                mind = d;
            }
        }
        res.add_constraint(n->coefficient(x).get_si()*mind*x
                           + n->coefficient(y).get_si()*mind*y
                           + mind*n->inhomogeneous_term().get_si()
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
        int max = INT_MIN;
        long maxd = 1;
        long d;
        for (auto c = C_p.begin(); c != C_p.end(); ++c) {
            int tmp = n->coefficient(x).get_si()*c->coefficient(x).get_si()
                + n->coefficient(y).get_si()*c->coefficient(y).get_si();

            d = c->divisor().get_si();

            if (max == INT_MIN || tmp*maxd > max*d) {
                max = tmp;
                maxd = d;
            }
        }
        res.add_constraint(n->coefficient(x).get_si()*maxd*x
                           + n->coefficient(y).get_si()*maxd*y
                           + maxd*n->inhomogeneous_term().get_si()
                           + max >= 0);
    }

    auto C_h = C.constraints();
    auto N_p = N.generators();
    for (auto c = C_h.begin(); c != C_h.end(); ++c) {
        int max = INT_MIN;
        long maxd = 1;
        long d;
        for (auto n = N_p.begin(); n != N_p.end(); ++n) {
            int tmp = c->coefficient(x).get_si()*n->coefficient(x).get_si()
                + c->coefficient(y).get_si()*n->coefficient(y).get_si();

            d = n->divisor().get_si();

            if (max == INT_MIN || tmp*maxd > max*d) {
                max = tmp;
                maxd = d;
            }
        }
        res.add_constraint(-c->coefficient(x).get_si()*maxd*x
                           - c->coefficient(y).get_si()*maxd*y
                           + maxd*c->inhomogeneous_term().get_si()
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
            res.add_generator(point((v_a->coefficient(x) + v_b->coefficient(x))*x
                                    + (v_a->coefficient(y) + v_b->coefficient(y))*y));

        }
    }
    return res;
}

void print_points(C_Polyhedron P) {
    auto gs = P.generators();
    cout << "[";
    for (auto g = gs.begin(); g != gs.end(); ++g) {
        cout << "[" << g->coefficient(x).get_d() / g->divisor().get_d() << ","
             << g->coefficient(y).get_d() / g->divisor().get_d() << "], ";
    }
    cout << "\b\b]\n\n";
}

void print_points(vector<C_Polyhedron> P) {
    for (auto p = P.begin(); p != P.end(); ++p) {
        print_points(*p);
    }
}

bool intersects(C_Polyhedron A, vector<C_Polyhedron> B) {
    for (auto b = B.begin(); b != B.end(); ++b) {
        A = C_Polyhedron(A);
        A.intersection_assign(*b);

        if (A.affine_dimension() > 0) {
            return true;
        }
    }

    return false;
}

bool can_translate_into(C_Polyhedron P,
                        C_Polyhedron P_over,
                        C_Polyhedron Nc,
                        vector<C_Polyhedron> Nd) {
    auto Ncc = C_Polyhedron(Nc);
    Ncc.intersection_assign(P_over);

    vector<C_Polyhedron> Ndd;
    for (auto d = Nd.begin(); d != Nd.end(); ++d) {
        auto tmp = C_Polyhedron(*d);
        tmp.intersection_assign(P_over);

        if (!tmp.is_empty()) {
            Ndd.push_back(tmp);
        }
    }

    auto U1 = translate_into(P, Ncc);
    auto U2 = translate_touching(P, Ndd);

    return subset(U1, U2.begin(), U2.end());
}

C_Polyhedron i2p(tuple<I, I> x_int) {

    int64_t nl0, dl0;
    rat_approx(get<0>(x_int).lower(), INT32_MAX, &nl0, &dl0);
    int64_t nl1, dl1;
    rat_approx(get<1>(x_int).lower(), INT32_MAX, &nl1, &dl1);
    int64_t nu0, du0;
    rat_approx(get<0>(x_int).upper(), INT32_MAX, &nu0, &du0);
    int64_t nu1, du1;
    rat_approx(get<1>(x_int).upper(), INT32_MAX, &nu1, &du1);

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

	while (f != floor(f)) { n <<= 1; f *= 2; }
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
        //auto v = P->generators();
        //for (auto p = v.begin(); p != v.end(); ++p) {
        //res.add_generator(*p);
        //}
        res.add_generators(P->generators());
    }

    return C_Polyhedron(res.minimized_generators());
}
