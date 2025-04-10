#include <fstream>

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

void print_points(const ninterval_t& x, ostream& f) {
    f << "[";
    for (int i = 0; i < NDIM; i++) {
        f << "[" << x[i].lower() << "," << x[i].upper() << "];";
    }
    f << "]\n";
}
