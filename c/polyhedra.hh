#include <algorithm>
#include <iterator>
#include <ppl.hh>
#include <stdio.h>
#include <time.h>
#include <type_traits>
#include <vector>

using namespace Parma_Polyhedra_Library;
using namespace std;

vector<C_Polyhedron> regiondiff(C_Polyhedron P,
                                vector<C_Polyhedron>::iterator curr,
                                vector<C_Polyhedron>::iterator end);
void print_points(vector<C_Polyhedron> P);
void print_points(C_Polyhedron P);
bool intersects(C_Polyhedron A, vector<C_Polyhedron> B);

C_Polyhedron translate_into(C_Polyhedron C, C_Polyhedron N);
vector<C_Polyhedron> translate_into(C_Polyhedron C, vector<C_Polyhedron> N, C_Polyhedron D);
C_Polyhedron translate_touching(C_Polyhedron C, C_Polyhedron N);
vector<C_Polyhedron> translate_touching(C_Polyhedron C, vector<C_Polyhedron> N);
