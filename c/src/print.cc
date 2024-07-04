#include "invariant.hh"
#include <fstream>


static Variable x(0);
static Variable y(1);

void print_points(const list<IntervalData>& P) {
    for (const auto& p : P) {
        print_points(p.poly);
    }
}

void fprint_points(const list<IntervalData>& P, string fname) {
    bool append = false;
    for (const auto& p : P) {
        fprint_points(p.poly, fname, append);
        append = true;
    }
}

void print_over(const list<IntervalData>& P) {
    for (const auto& p : P) {
        print_points(p.P_over);
    }
}

void print_u_over(const list<IntervalData>& P) {
    for (const auto& p : P) {
        print_points(p.P_u_over);
    }
}

void print_points(const C_Polyhedron& P) {
    cout << "[";
    for (const auto& g : P.generators()) {
        if (g.is_point()) {
            cout << "[" << g.coefficient(x).get_d() / g.divisor().get_d() << ","
                 << g.coefficient(y).get_d() / g.divisor().get_d() << "]; ";
        }
        else {
            g.ascii_dump();
        }
    }
    cout << "\b\b]\n\n";
}

void fprint_points(const C_Polyhedron& P, string fname, bool append) {
    ofstream f;
    if (append) {
        f.open(fname, ios_base::app);
    } else {
        f.open(fname);
    }

    f << "[";
    for (const auto& g : P.generators()) {
        if (g.is_point()) {
            f << "[" << g.coefficient(x).get_d() / g.divisor().get_d() << ","
                 << g.coefficient(y).get_d() / g.divisor().get_d() << "]; ";
        }
        else {
            g.ascii_dump();
        }
    }
    f << "]\n";
}

void print_points(const vector<C_Polyhedron>& P) {
    for (const auto& p : P) {
        print_points(p);
    }
}
void fprint_points(const vector<C_Polyhedron>& P, string fname) {
    bool append = false;
    for (const auto& p : P) {
        fprint_points(p, fname, append);
        append = true;
    }
}
