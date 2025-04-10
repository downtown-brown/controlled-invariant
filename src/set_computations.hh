#include "set_types.hh"

pair<IntervalData, IntervalData> bisect(IntervalData x);

bool wider_than(const ninterval_t& interval, nvec_t epsilon);
nvec_t median(const ninterval_t& interval);

vector<C_Polyhedron> regiondiff(C_Polyhedron P,
                                vector<C_Polyhedron>::iterator curr,
                                vector<C_Polyhedron>::iterator end);
bool subset(C_Polyhedron P,
            vector<C_Polyhedron>::iterator curr,
            vector<C_Polyhedron>::iterator end);
bool can_translate_into(const C_Polyhedron& P,
                        const C_Polyhedron& P_over,
                        const C_Polyhedron& Nc,
                        const vector<C_Polyhedron>& Nd);
vector<C_Polyhedron> translate_into(const C_Polyhedron& P,
                                    const C_Polyhedron& P_over,
                                    const C_Polyhedron& Nc,
                                    const vector<C_Polyhedron>& Nd);
C_Polyhedron translate_into(const C_Polyhedron &C, const C_Polyhedron& N);
vector<C_Polyhedron> translate_into(const C_Polyhedron& C,
                                    const vector<C_Polyhedron>& N,
                                    const C_Polyhedron& D);
C_Polyhedron translate_touching(const C_Polyhedron& C, const C_Polyhedron& N);
vector<C_Polyhedron> translate_touching(const C_Polyhedron& C, const vector<C_Polyhedron>& N);
bool intersects(const C_Polyhedron& A, const vector<C_Polyhedron>& B);
bool intersects(const C_Polyhedron& A, const C_Polyhedron& B);
bool intersects(const ninterval_t& A, const ninterval_t& B);

C_Polyhedron operator+(const C_Polyhedron& a, const C_Polyhedron& b);
ninterval_t operator+(const ninterval_t& A, const ninterval_t& B);

C_Polyhedron i2p(ninterval_t x_int);
int64_t rat_approx(double f, int64_t den);

C_Polyhedron convexhull(const vector<C_Polyhedron>& P_v);
ninterval_t intervalhull(const vector<IntervalData>& Omega);

void merge(vector<IntervalData>& Omega);
void merge_fast(vector<IntervalData>& Omega);

bool comp_intervals(const IntervalData &A, const IntervalData &B, int i);
