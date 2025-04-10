#ifndef PRINT_HH
#define PRINT_HH

#include <fstream>
#include "set_types.hh"

void print_points(const C_Polyhedron& P, ostream& f=cout);
void print_points(const ninterval_t& P, ostream& f=cout);


template <typename T>
void print_points(const vector<T>& P) {
    for (const T& p : P) {
        print_points(p);
    }
}

template <typename T>
void fprint_points(const T& P, string fname, bool append) {
    ofstream f;
    if (append) {
        f.open(fname, ios_base::app);
    } else {
        f.open(fname);
    }

    print_points(P, f);
}

template <typename T>
void fprint_points(const vector<T>& P, string fname) {
    bool append = false;
    for (const T& p : P) {
        fprint_points(p, fname, append);
        append = true;
    }
}

#endif
