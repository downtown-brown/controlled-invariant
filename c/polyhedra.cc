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

#include "polyhedra.hh"

using namespace Parma_Polyhedra_Library;
using namespace std;

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
