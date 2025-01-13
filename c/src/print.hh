#include "set_types.hh"

void print_points(const vector<IntervalData>& P);
void fprint_points(const vector<IntervalData>& P, string fname);
void print_points(const vector<C_Polyhedron>& P);
void fprint_points(const vector<C_Polyhedron>& P, string fname);
void fprint_points(const C_Polyhedron& P, string fname, bool append);
void print_points(const C_Polyhedron& P, ostream& f=cout);
