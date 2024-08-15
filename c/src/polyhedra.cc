#include <algorithm>
#include <cfloat>
#include <climits>
#include <cmath>
#include <cstdint>
#include <ctime>
#include <fstream>
#include <ios>
#include <iterator>
#include <ppl.hh>
#include <stdio.h>
#include <sys/types.h>
#include <time.h>
#include <type_traits>
#include <vector>

#include <boost/numeric/interval.hpp>
#include "invariant.hh"


using namespace Parma_Polyhedra_Library;
using namespace std;
using namespace boost::numeric;
using namespace interval_lib;

static Variable x(0);
static Variable y(1);

vector<C_Polyhedron> regiondiff(C_Polyhedron P,
                                vector<C_Polyhedron>::iterator curr,
                                vector<C_Polyhedron>::iterator end)
{
    vector<C_Polyhedron> res;

    while (true) {

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

    for (const auto& i : curr->constraints()) {
        auto tmp = C_Polyhedron(P);

        tmp.add_constraint(Constraint(-i.coefficient(x)*x
                                      - i.coefficient(y)*y
                                      - i.inhomogeneous_term() >= 0));

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

        P.add_constraint(i);
    }

    return res;
}

bool subset(C_Polyhedron P,
            vector<C_Polyhedron>::iterator curr,
            vector<C_Polyhedron>::iterator end)
{
    if (P.is_empty() || curr == end) {
        return true;
    }

    while (true) {

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


    for (const auto& i : curr->constraints()) {
        auto tmp = C_Polyhedron(P);

        tmp.add_constraint(Constraint(-i.coefficient(x)*x
                                      - i.coefficient(y)*y
                                      - i.inhomogeneous_term() >= 0));

        if (tmp.affine_dimension() < P.affine_dimension()) {
            continue;
        }

        if (curr == end - 1) {
            return false;
        } else if (!subset(tmp, curr + 1, end)) {
                return false;
        }

        P.add_constraint(i);
    }

    return true;
}

C_Polyhedron translate_into(const C_Polyhedron& C, const C_Polyhedron& N) {
    C_Polyhedron res(2);
    for (const auto& n : N.constraints()) {
        GMP_Integer min = 0;
        GMP_Integer mind = 1;
        GMP_Integer d;
        for (const auto& c : C.generators()) {
            GMP_Integer tmp = n.coefficient(x)*c.coefficient(x)
                + n.coefficient(y)*c.coefficient(y);

            d = c.divisor();

            if (tmp*mind < min*d || min == 0) {
                min = tmp;
                mind = d;
            }
        }
        res.add_constraint(n.coefficient(x)*mind*x
                           + n.coefficient(y)*mind*y
                           + mind*n.inhomogeneous_term()
                           + min >= 0);
    }

    return res;
}

vector<C_Polyhedron> translate_into(const C_Polyhedron& C,
                                    const vector<C_Polyhedron>& N,
                                    const C_Polyhedron& D) {

    auto U1 = translate_into(C, D);
    auto U2 = translate_touching(C, N);

    return regiondiff(U1, U2.begin(), U2.end());
}

C_Polyhedron translate_touching(const C_Polyhedron& C, const C_Polyhedron& N) {
    C_Polyhedron res(2);
    GMP_Integer max;
    GMP_Integer maxd;
    GMP_Integer d;
    GMP_Integer tmp;
    for (const auto& n : N.constraints()) {
        max = 0;
        maxd = 1;
        for (const auto& c : C.generators()) {
            tmp = n.coefficient(x)*c.coefficient(x)
                + n.coefficient(y)*c.coefficient(y);

            d = c.divisor();

            if (tmp*maxd > max*d || max == 0) {
                max = tmp;
                maxd = d;
            }
        }
        res.add_constraint(n.coefficient(x)*maxd*x
                           + n.coefficient(y)*maxd*y
                           + maxd*n.inhomogeneous_term()
                           + max >= 0);
    }

    for (const auto& c : C.constraints()) {
        max = 0;
        maxd = 1;
        for (const auto& n : N.generators()) {
            tmp = c.coefficient(x)*n.coefficient(x)
                + c.coefficient(y)*n.coefficient(y);

            d = n.divisor();

            if (tmp*maxd > max*d || max == 0) {
                max = tmp;
                maxd = d;
            }
        }
        res.add_constraint(-c.coefficient(x).get_d()*maxd*x
                           - c.coefficient(y).get_d()*maxd*y
                           + maxd*c.inhomogeneous_term().get_d()
                           + max >= 0);
    }

    return res;
}

vector<C_Polyhedron> translate_touching(const C_Polyhedron& C, const vector<C_Polyhedron>& N) {
    vector<C_Polyhedron> res;

    for (const auto& n : N) {
        res.emplace_back(translate_touching(C, n));
    }

    return res;
}

C_Polyhedron operator+(const C_Polyhedron& a, const C_Polyhedron& b) {
    C_Polyhedron res(2, EMPTY);
    for (const auto& v_a : a.generators()) {
        for (const auto& v_b : b.generators()) {
            res.add_generator(point((v_a.coefficient(x)*v_b.divisor()
                                     + v_b.coefficient(x)*v_a.divisor())*x
                                    + (v_a.coefficient(y)*v_b.divisor()
                                       + v_b.coefficient(y)*v_a.divisor())*y,
                                    v_a.divisor()*v_b.divisor()));

        }
    }
    return C_Polyhedron(res.minimized_generators());
}

bool has_separating_plane(const C_Polyhedron& A, const C_Polyhedron& B) {
    for (const auto& a : A.constraints()) {
        bool sep = true;
        for (const auto& b : B.generators()) {
            if (a.coefficient(x)*b.coefficient(x)
                + a.coefficient(y)*b.coefficient(y)
                + a.inhomogeneous_term()*b.divisor() >= 0) {
                sep = false;
                break;
            }
        }
        if (sep) {
            return true;
        }
    }
    return false;
}

bool intersects(const C_Polyhedron& A, const C_Polyhedron& B) {
    return !has_separating_plane(A, B) && !has_separating_plane(B, A);
}

bool intersects(const C_Polyhedron& A, const vector<C_Polyhedron>& B) {
    for (auto b : B) {
        if (intersects(A, b)) {
            return true;
        }
    }
    return false;
}

extern double t_ndd_tot;
extern double t_sub_tot;
vector<C_Polyhedron> translate_into(const C_Polyhedron& P,
                                    const C_Polyhedron& P_over,
                                    const C_Polyhedron& Nc,
                                    const vector<C_Polyhedron>& Nd) {
    auto Ncc = C_Polyhedron(Nc);
    Ncc.intersection_assign(P_over);

    clock_t start = clock();
    if (!Nd.empty()) {
        auto U1 = translate_into(P, Ncc);

        vector<C_Polyhedron> Ndd;
        for (const auto& d : Nd) {
            auto tmp = C_Polyhedron(d);
            tmp.intersection_assign(P_over);

            if (!tmp.is_empty()) {
                Ndd.push_back(tmp);
            }
        }

        clock_t t1 = clock();
        //printf("Computed Ndd (%ld -> %ld elements) %.10f seconds\n\n",
        //       Nd.size(), Ndd.size(), (double)(t1 - start)/CLOCKS_PER_SEC);

        t_ndd_tot += (double)(t1-start)/CLOCKS_PER_SEC;

        auto U2 = translate_touching(P, Ndd);

        if (U2.empty()) {
            return {U1};
        } else {
            return regiondiff(U1, U2.begin(), U2.end());
        }
        clock_t end = clock();
        //printf("Computed regiondiff %.10f seconds\n\n", (double)(end - t1)/CLOCKS_PER_SEC);
        t_sub_tot += (double)(end - t1)/CLOCKS_PER_SEC;

    } else {
        auto U1 = translate_into(P, Ncc);
        return {U1};
    }
    return {};
}

bool can_translate_into(const C_Polyhedron& P,
                        const C_Polyhedron& P_over,
                        const C_Polyhedron& Nc,
                        const vector<C_Polyhedron>& Nd) {
    auto Ncc = C_Polyhedron(Nc);
    Ncc.intersection_assign(P_over);

    clock_t start = clock();
    if (!Nd.empty()) {
        auto U1 = translate_into(P, Ncc);

        vector<C_Polyhedron> Ndd;
        for (const auto& d : Nd) {
            if (intersects(P_over, d)) {
                auto tmp = d;
                tmp.intersection_assign(P_over);
                Ndd.push_back(tmp);
            }
        }

        clock_t t1 = clock();
        t_ndd_tot += (double)(t1-start)/CLOCKS_PER_SEC;

        auto U2 = translate_touching(P, Ndd);

        bool res = !subset(U1, U2.begin(), U2.end()) && !U1.is_empty();
        clock_t end = clock();
        t_sub_tot += (double)(end - t1)/CLOCKS_PER_SEC;

        return res;

    } else {
        auto U1 = translate_into(P, Ncc);
        return !U1.is_empty();
    }
}

C_Polyhedron i2p(nI x_int) {

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

C_Polyhedron convexhull(const vector<C_Polyhedron>& P_v) {
    C_Polyhedron res(2, EMPTY);

    for (const auto& P : P_v) {
        res.add_generators(P.generators());
    }

    return C_Polyhedron(res.minimized_generators());
}

C_Polyhedron intervalhull(const vector<C_Polyhedron>& P_v) {
    double l0 = HUGE_VAL, l1 = HUGE_VAL, u0 = -HUGE_VAL, u1 = -HUGE_VAL, px, py;

    for (const auto& P : P_v) {
        for (const auto& p : P.generators()) {
            px = p.coefficient(x).get_d() / p.divisor().get_d();
            py = p.coefficient(y).get_d() / p.divisor().get_d();
            if (px < l0)
                l0 = px;
            if (py < l1)
                l1 = py;
            if (px > u0)
                u0 = px;
            if (py > u1)
                u1 = py;
        }
    }

    return i2p({I(l0, u0), I(l1, u1)});
}

void merge_once(list<IntervalData> &Omega) {
    for (auto Ait = Omega.begin(); Ait != Omega.end(); Ait++) {
        for (auto Bit = Omega.begin(); Bit != Omega.end(); Bit++) {
            if (Ait == Bit) {
                continue;
            }
            auto res = merge(*Ait, *Bit);
            if (res.size() == 1) {
                *Ait = res.front();
                Omega.erase(Bit);
                return;
            }
        }
    }
}

void merge(list<IntervalData> &Omega) {
    auto len = Omega.size();
    ulong len2 = 0;
    while (len != len2) {
        len2 = len;
        merge_once(Omega);
        len = Omega.size();
    }
}

list<IntervalData> merge(const IntervalData& A, const IntervalData& B) {
    bool must_be_eq = false;
    int i_n = 0;

    for (int i = 0; i < n; i++) {
        if ((A.interval[i].lower() > B.interval[i].upper()) ||
            (A.interval[i].upper() < B.interval[i].lower())) {

            return {};

        } else if ((A.interval[i].lower() != B.interval[i].lower()) ||
                   (A.interval[i].upper() != B.interval[i].upper())) {

            if (must_be_eq) {
                return {};
            } else {
                must_be_eq = true;
                i_n = i;
            }
        }
    }

    auto interval = A.interval;

    if (A.interval[i_n].lower() == B.interval[i_n].upper()) {
        interval[i_n] = I(B.interval[i_n].lower(), A.interval[i_n].upper());
    } else {
        interval[i_n] = I(A.interval[i_n].lower(), B.interval[i_n].upper());
    }

    return {IntervalData(interval)};
}

nI operator+(const nI& A, const nI& B) {
    nI res;
    for (int i = 0; i < n; i++) {
        res[i] = A[i] + B[i];
    }

    return res;
}
