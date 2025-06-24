#include <optional>
#include <list>
#include <ppl.hh>
#include <ranges>

#include "set_types.hh"
#include "print.hh"
#include "set_computations.hh"

/*
  Polyhedron computations
*/

vector<C_Polyhedron> regiondiff(C_Polyhedron P,
                                vector<C_Polyhedron>::iterator curr,
                                vector<C_Polyhedron>::iterator end)
{
    vector<C_Polyhedron> res;

    /* Skip over polytopes that don't intersect with P */
    while (true) {
        C_Polyhedron tmp = C_Polyhedron(P);
        tmp.intersection_assign(*curr);

        if (tmp.affine_dimension() >= P.affine_dimension()) {
            break;
        }

        if (++curr == end) {
            res.push_back(P);
            return res;
        }
    }

    for (const Constraint& c : curr->constraints()) {

        Linear_Expression c_flip;
        for (int i = 0; i < NDIM; i++) {
            Variable v(i);
            c_flip -= c.coefficient(v)*v;
        }

        /* This condition means the flipped constraint will not affect P */
        if (P.relation_with(Constraint(c_flip - c.inhomogeneous_term() > 0)) == Poly_Con_Relation::is_disjoint()) {
            continue;
        }

        C_Polyhedron tmp = C_Polyhedron(P);
        tmp.add_constraint(c_flip - c.inhomogeneous_term() >= 0);

        if (curr < end - 1) {
            vector<C_Polyhedron> tmp2 = regiondiff(tmp, curr + 1, end);
            res.insert(res.end(), tmp2.begin(), tmp2.end());
        }
        else {
            res.push_back(tmp);
        }

        P.add_constraint(c);
    }

    return res;
}

bool subset(C_Polyhedron P,
            vector<C_Polyhedron>::iterator curr,
            vector<C_Polyhedron>::iterator end)
{
    /* Empty set is a subset of anything */
    if (P.is_empty()) {
        return true;
    }

    /* P is a single point, check if it is contained in any one polytope */
    if (P.affine_dimension() == 0) {
        while (curr != end) {
            if (curr++->contains(P)) {
                return true;
            }
        }
        return false;
    }

    /* Skip over polytopes that don't intersect with P */
    while (true) {
        /* Nothing is a subset of the empty set */
        if (curr == end) {
            return false;
        }

        if (intersects(P, *curr)) {
            break;
        }

        curr++;
    }

    for (const Constraint& c : curr->constraints()) {
        /* Flip the direction of the constraint */
        Linear_Expression c_flip;
        for (int i = 0; i < NDIM; i++) {
            Variable v(i);
            c_flip -= c.coefficient(v)*v;
        }

        /* This condition means the flipped constraint will not affect P */
        if (P.relation_with(Constraint(c_flip - c.inhomogeneous_term() > 0)) == Poly_Con_Relation::is_disjoint()) {
          continue;
        }

        C_Polyhedron tmp = P;
        tmp.add_constraint(c_flip - c.inhomogeneous_term() >= 0);

        if (!subset(tmp, curr + 1, end)) {
            return false;
        }

        P.add_constraint(c);
    }

    return true;
}

C_Polyhedron translate_into(const C_Polyhedron& C, const C_Polyhedron& N) {
    C_Polyhedron res(NDIM);
    for (const Constraint& n : N.constraints()) {
        GMP_Integer min = 0;
        GMP_Integer mind = 1;
        GMP_Integer d;
        for (const Generator& c : C.generators()) {
            GMP_Integer tmp = 0;
            for (int i = 0; i < NDIM; i++) {
                tmp += n.coefficient(Variable(i))*c.coefficient(Variable(i));
            }

            d = c.divisor();

            if (tmp*mind < min*d || min == 0) {
                min = tmp;
                mind = d;
            }
        }

        Linear_Expression con_new;
        for (int i = 0; i < NDIM; i++) {
            Variable v(i);
            con_new += n.coefficient(v)*mind*v;
        }

        res.add_constraint(con_new + mind*n.inhomogeneous_term() + min >= 0);
    }

    return res;
}

C_Polyhedron translate_touching(const C_Polyhedron& C, const C_Polyhedron& N) {
    C_Polyhedron res(NDIM);
    GMP_Integer max;
    GMP_Integer maxd;
    GMP_Integer d;
    for (const Constraint& n : N.constraints()) {
        max = 0;
        maxd = 1;
        for (const Generator& c : C.generators()) {
            GMP_Integer tmp = 0;
            for (int i = 0; i < NDIM; i++) {
                tmp += n.coefficient(Variable(i))*c.coefficient(Variable(i));
            }

            d = c.divisor();

            if (tmp*maxd > max*d || max == 0) {
                max = tmp;
                maxd = d;
            }
        }

        Linear_Expression con_new;
        for (int i = 0; i < NDIM; i++) {
            Variable v(i);
            con_new += n.coefficient(v)*maxd*v;
        }

        res.add_constraint(con_new + maxd*n.inhomogeneous_term() + max >= 0);
    }

    for (const Constraint& c : C.constraints()) {
        max = 0;
        maxd = 1;
        for (const Generator& n : N.generators()) {
            GMP_Integer tmp = 0;
            for (int i = 0; i < NDIM; i++) {
                tmp += n.coefficient(Variable(i))*c.coefficient(Variable(i));
            }

            d = n.divisor();

            if (tmp*maxd > max*d || max == 0) {
                max = tmp;
                maxd = d;
            }
        }
        Linear_Expression con_new;
        for (int i = 0; i < NDIM; i++) {
            Variable v(i);
            con_new -= c.coefficient(v)*maxd*v;
        }

        res.add_constraint(con_new + maxd*c.inhomogeneous_term() + max >= 0);
    }

    return res;
}

vector<C_Polyhedron> translate_into(const C_Polyhedron& C,
                                    const vector<C_Polyhedron>& N,
                                    const C_Polyhedron& D) {

    C_Polyhedron U1 = translate_into(C, D);
    vector<C_Polyhedron> U2 = translate_touching(C, N);

    return regiondiff(U1, U2.begin(), U2.end());
}

C_Polyhedron operator+(const C_Polyhedron &a, const C_Polyhedron &b) {
    C_Polyhedron res(NDIM, EMPTY);
    for (const Generator& v_a : a.minimized_generators()) {
        for (const Generator& v_b : b.minimized_generators()) {
            Linear_Expression tmp;
            for (int i = 0; i < NDIM; i++) {
                Variable v(i);
                tmp += (v_a.coefficient(v)*v_b.divisor() + v_b.coefficient(v)*v_a.divisor())*v;
            }
            res.add_generator(point(tmp, v_a.divisor()*v_b.divisor()));

        }
    }
    return C_Polyhedron(res.minimized_generators());
}

bool has_separating_plane(const C_Polyhedron& A, const C_Polyhedron& B) {
    for (const Constraint& a : A.constraints()) {
        if (B.relation_with(a) == Poly_Con_Relation::is_disjoint()) {
            return true;
        }
    }
    return false;
}

/*
  This function only works on full-dimensional polytopes.
*/
bool intersects(const C_Polyhedron& A, const C_Polyhedron& B) {
    return !has_separating_plane(A, B) && !has_separating_plane(B, A);
}

bool intersects(const C_Polyhedron& A, const vector<C_Polyhedron>& B) {
    for (const C_Polyhedron& b : B) {
        if (intersects(A, b)) {
            return true;
        }
    }
    return false;
}

vector<C_Polyhedron> translate_into(const C_Polyhedron& P,
                                    const C_Polyhedron& P_over,
                                    const C_Polyhedron& Nc,
                                    const vector<C_Polyhedron>& Nd) {
    C_Polyhedron Ncc = Nc;
    Ncc.intersection_assign(P_over);

    if (!Nd.empty()) {
        C_Polyhedron U1 = translate_into(P, Ncc);

        vector<C_Polyhedron> Ndd;
        for (const C_Polyhedron& d : Nd) {
            C_Polyhedron tmp = C_Polyhedron(d);
            tmp.intersection_assign(P_over);

            if (!tmp.is_empty()) {
                Ndd.push_back(tmp);
            }
        }

        vector<C_Polyhedron> U2 = translate_touching(P, Ndd);

        if (U2.empty()) {
            return {U1};
        } else {
            return regiondiff(U1, U2.begin(), U2.end());
        }
    } else {
        C_Polyhedron U1 = translate_into(P, Ncc);
        return {U1};
    }
    return {};
}

bool can_translate_into(const C_Polyhedron& P,
                        const C_Polyhedron& P_over,
                        const C_Polyhedron& Nc,
                        const vector<C_Polyhedron>& Nd) {
    C_Polyhedron Ncc = Nc;
    Ncc.intersection_assign(P_over);

    auto Ndd = Nd | views::filter([P_over](const C_Polyhedron& x){ return intersects(P_over, x); });

    C_Polyhedron U1 = translate_into(P, Ncc);
    vector<C_Polyhedron> U2 = translate_touching(P, Ndd);

    return !subset(U1, U2.begin(), U2.end());
}

int64_t rat_approx(double f, int64_t den) {
    return (int64_t) (f * (double) den);
}

C_Polyhedron i2p(ninterval_t x) {
    const int64_t den = 100000000;
    C_Polyhedron res(NDIM);

    for (int i = 0; i < NDIM; i++) {
        Variable v(i);
        res.add_constraint(den*v - (rat_approx(x[i].lower(), den) - 1) >= 0);
        res.add_constraint(den*v - (rat_approx(x[i].upper(), den) + 1) <= 0);
    }

    return res;
}

/*
  Interval functions
*/

bool intersects(const ninterval_t& A, const ninterval_t& B) {
    for (int i = 0; i < NDIM; i++) {
        if (!overlap(A[i], B[i])) {
            return false;
        }
    }

    return true;
}


bool wider_than(const ninterval_t& x, nvec_t epsilon) {
    for (int i = 0; i < NDIM; i++) {
        if (width(x[i]) > epsilon[i]) {
            return true;
        }
    }

    return false;
}

nvec_t median(const ninterval_t& x) {
    nvec_t res;

    for (int i = 0; i < NDIM; i++) {
        res[i] = median(x[i]);
    }

    return res;
}

ninterval_t intervalhull(const vector<ninterval_t>& Omega) {
    ninterval_t res = Omega[0];
    for (const ninterval_t& x: Omega) {
        for (int i = 0; i < NDIM; i++) {
            res[i] = hull(res[i], x[i]);
        }
    }

    return res;
}

C_Polyhedron convexhull(const vector<C_Polyhedron>& P_v) {
    C_Polyhedron res(NDIM, EMPTY);

    for (const C_Polyhedron& P : P_v) {
        res.add_generators(P.generators());
    }

    return C_Polyhedron(res.minimized_generators());
}

bool comp_intervals(const ninterval_t &A, const ninterval_t &B, int i) {
    for (int j = 0; j < NDIM; j++) {
        if (A[i].lower() < B[i].lower()) {
            return true;
        } else if (A[i].lower() > B[i].lower()) {
            return false;
        } else if (A[i].upper() < B[i].upper()) {
            return true;
        } else if (A[i].upper() > B[i].upper()) {
            return false;
        }
        i = (i + 1) % NDIM;
    }

    if (A[i].upper() <= B[i].lower()) {
        return true;
    } else {
        return false;
    }
}

bool merge(ninterval_t& A, const ninterval_t& B, int i_n) {
    for (int i = 0; i < NDIM; i++) {
        if (i != i_n && ((A[i].lower() != B[i].lower()) ||
                         (A[i].upper() != B[i].upper()))) {
            return false;
        }
    }

    if (A[i_n].upper() == B[i_n].lower()) {
        A[i_n] =
            interval_t(A[i_n].lower(), B[i_n].upper());
        return true;
    } else {
        return false;
    }
}

void merge_once(list<ninterval_t> &Omega) {
    for (int i = 0; i < NDIM; i++) {
        Omega.sort([i](const ninterval_t &A, const ninterval_t &B) {
            return comp_intervals(A, B, i);
        });

        auto l = Omega.begin();
        auto r = std::next(l);

        while (r != Omega.end()) {
            if (merge(*l, *r, (i + NDIM - 1) % NDIM)) {
                r = Omega.erase(r);
            } else {
                l = r++;
            }
        }
    }
}

void merge(vector<ninterval_t> &Omega) {
    int len = Omega.size();
    int len2 = 0;
    int iter = 0;
    std::list<ninterval_t> O2(Omega.begin(), Omega.end());
    while (len != len2) {
        cout << "Merge iter: " << iter << ", len: " << len << endl;
        len2 = len;
        merge_once(O2);
        len = O2.size();
        iter++;
    }
    Omega = vector<ninterval_t>(O2.begin(), O2.end());
}

ninterval_t operator+(const ninterval_t& A, const ninterval_t& B) {
    ninterval_t res;
    for (int i = 0; i < NDIM; i++) {
        res[i] = A[i] + B[i];
    }

    return res;
}

pair<IntervalData, IntervalData> bisect(IntervalData x) {
    double max_width = 0;
    int max_dim = 0;

    for (int i = 0; i < NDIM; i++) {
        double curr_width = width(x.interval[i]);
        if (curr_width > max_width) {
            max_width = curr_width;
            max_dim = i;
        }
    }

    ninterval_t l = x.interval;
    ninterval_t r = x.interval;

    pair<interval_t, interval_t> tmp = bisect(x.interval[max_dim]);

    l[max_dim] = get<0>(tmp);
    r[max_dim] = get<1>(tmp);

    IntervalData x_l(l, x.node);
    IntervalData x_r(r, x.node);

    return {x_l, x_r};
}
