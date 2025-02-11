#include <optional>

#include "set_types.hh"

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
        C_Polyhedron tmp = C_Polyhedron(P);

        Linear_Expression con_new;
        for (int i = 0; i < NDIM; i++) {
            Variable v(i);
            con_new -= c.coefficient(v)*v;
        }
        tmp.add_constraint(con_new - c.inhomogeneous_term() >= 0);

        if (tmp.affine_dimension() < P.affine_dimension()) {
            continue;
        }
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

    /* Skip over polytopes that don't intersect with P */
    while (true) {
        C_Polyhedron tmp = C_Polyhedron(P);
        tmp.intersection_assign(*curr);
        if (tmp.affine_dimension() >= P.affine_dimension()) {
            break;
        }

        if (++curr == end) {
            return false;
        }
    }


    for (const Constraint& c : curr->constraints()) {
        C_Polyhedron tmp = C_Polyhedron(P);

        Linear_Expression con_new;
        for (int i = 0; i < NDIM; i++) {
            Variable v(i);
            con_new -= c.coefficient(v)*v;
        }
        tmp.add_constraint(con_new - c.inhomogeneous_term() >= 0);

        if (tmp.affine_dimension() < P.affine_dimension()) {
            continue;
        }

        if (curr == end - 1 || !subset(tmp, curr + 1, end)) {
            return false;
        }

        P.add_constraint(c);
    }

    return true;
}

C_Polyhedron translate_into(const C_Polyhedron& C, const C_Polyhedron& N) {
    C_Polyhedron res(2);
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
    C_Polyhedron res(2);
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

vector<C_Polyhedron> translate_touching(const C_Polyhedron& C, const vector<C_Polyhedron>& N) {
    vector<C_Polyhedron> res;

    for (const C_Polyhedron& n : N) {
        res.emplace_back(translate_touching(C, n));
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

/*
   TODO: speed up using zonotopes
*/
C_Polyhedron operator+(const C_Polyhedron& a, const C_Polyhedron& b) {
    C_Polyhedron res(2, EMPTY);
    for (const Generator& v_a : a.generators()) {
        for (const Generator& v_b : b.generators()) {
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
        bool sep = true;
        for (const Generator& b : B.generators()) {
            if (Scalar_Products::sign(a, b) >= 0) {
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
    C_Polyhedron Ncc = C_Polyhedron(Nc);
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
    C_Polyhedron Ncc = C_Polyhedron(Nc);
    Ncc.intersection_assign(P_over);

    if (!Nd.empty()) {
        C_Polyhedron U1 = translate_into(P, Ncc);

        vector<C_Polyhedron> Ndd;
        for (const C_Polyhedron& d : Nd) {
            if (intersects(P_over, d)) {
                C_Polyhedron tmp = d;
                tmp.intersection_assign(P_over);
                Ndd.push_back(tmp);
            }
        }

        vector<C_Polyhedron> U2 = translate_touching(P, Ndd);
        if (U2.empty()) {
            return !U1.is_empty();
        } else {
            return !subset(U1, U2.begin(), U2.end());
        }
    } else {
        C_Polyhedron U1 = translate_into(P, Ncc);
        return !U1.is_empty();
    }
}

int64_t rat_approx(double f, int64_t den) {
    return (int64_t) (f * (double) den);
}

C_Polyhedron i2p(ninterval_t x) {
    const int64_t den = 1000000;
    C_Polyhedron res(NDIM);

    for (int i = 0; i < NDIM; i++) {
        Variable v(i);
        res.add_constraint(den*v - rat_approx(x[i].lower(), den) >= 0);
        res.add_constraint(den*v - rat_approx(x[i].upper(), den) <= 0);
    }

    return res;
}

bool wider_than(const ninterval_t& x, double epsilon) {
    for (int i = 0; i < NDIM; i++) {
        if (width(x[i]) > epsilon) {
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

C_Polyhedron convexhull(const vector<C_Polyhedron>& P_v) {
    C_Polyhedron res(NDIM, EMPTY);

    for (const C_Polyhedron& P : P_v) {
        res.add_generators(P.generators());
    }

    return C_Polyhedron(res.minimized_generators());
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

    IntervalData x_l(l);
    IntervalData x_r(r);
    x.lchild = &x_l;
    x.rchild = &x_r;

    return {x_l, x_r};
}

optional<IntervalData> merge(const IntervalData& A, const IntervalData& B) {
    bool must_be_eq = false;
    int i_n = 0;

    for (int i = 0; i < NDIM; i++) {
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

    ninterval_t interval = A.interval;

    if (A.interval[i_n].lower() == B.interval[i_n].upper()) {
        interval[i_n] = interval_t(B.interval[i_n].lower(), A.interval[i_n].upper());
    } else {
        interval[i_n] = interval_t(A.interval[i_n].lower(), B.interval[i_n].upper());
    }

    return IntervalData(interval);
}

void merge_once(vector<IntervalData> &Omega) {
    for (auto Ait = Omega.begin(); Ait != Omega.end() - 1; Ait++) {
        for (auto Bit = Ait+ 1; Bit != Omega.end(); Bit++) {
            optional<IntervalData> res = merge(*Ait, *Bit);
            if (res.has_value()) {
                *Ait = res.value();
                *Bit = Omega.back();
                Omega.pop_back();
                return;
            }
        }
    }
}

void merge(vector<IntervalData> &Omega) {
    int len = Omega.size();
    int len2 = 0;
    while (len != len2) {
        len2 = len;
        merge_once(Omega);
        len = Omega.size();
    }
}

ninterval_t operator+(const ninterval_t& A, const ninterval_t& B) {
    ninterval_t res;
    for (int i = 0; i < NDIM; i++) {
        res[i] = A[i] + B[i];
    }

    return res;
}
