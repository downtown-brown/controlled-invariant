#include <ostream>

#include "set_types.hh"
#include "print.hh"

void print_points(const C_Polyhedron& P, ostream& f) {
    f << "[";
    for (const Generator& g : P.generators()) {
        if (g.is_point()) {
            f << "[";
            for (int i = 0; i < NDIM; i++) {
                f << g.coefficient(Variable(i)).get_d() / g.divisor().get_d();
                if (i < NDIM - 1) f << ",";
            }
            f << "];";
        }
        else {
            g.ascii_dump();
        }
    }
    f << "]\n";
}

void fprint_points(const C_Polyhedron& P, string fname, bool append) {
    ostream f;
    if (append) {
        f.open(fname, ios_base::app);
    } else {
        f.open(fname);
    }

    print_points(P, f);
}

void print_points(const vector<IntervalData>& P) {
    for (const IntervalData& p : P) {
        print_points(p.poly);
    }
}

void fprint_points(const vector<IntervalData>& P, string fname) {
    bool append = false;
    for (const IntervalData& p : P) {
        fprint_points(p.poly, fname, append);
        append = true;
    }
}

void print_points(const vector<C_Polyhedron>& P) {
    for (const C_Polyhedron& p : P) {
        print_points(p);
    }
}
void fprint_points(const vector<C_Polyhedron>& P, string fname) {
    bool append = false;
    for (const C_Polyhedron& p : P) {
        fprint_points(p, fname, append);
        append = true;
    }
}
