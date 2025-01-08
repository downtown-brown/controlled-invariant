#include <algorithm>
#include <cstdint>
#include <iterator>
#include <ostream>
#include <ppl.hh>
#include <stdio.h>
#include <time.h>
#include <type_traits>
#include <vector>
#include <boost/numeric/interval.hpp>

#include "invariant.hh"

void tests(void);

using namespace Parma_Polyhedra_Library;
using namespace std;
using namespace boost::numeric;
using namespace interval_lib;


static Variable x(0);
static Variable y(1);

uint8_t stop = 0;


int main() {
    //tests();

    IntervalData Omega({interval_t(-6, 6), interval_t(-6, 6)});

    auto start = clock();
    vector<IntervalData> res;
    uint64_t ii = 0;
    while (!stop) {
        if (ii == 0) {
            res = I_approx({Omega});
        } else {
            res = I_approx(res);
        }
        ii++;
    }
    auto end = clock();
    printf("Computed %ld iters I_approx %.10f seconds\n\n", ii, (double)(end - start)/CLOCKS_PER_SEC);

    /*
    print_points(res);

    cout << res.size() << endl;


    cout << "inputs:\n";
    auto U = U_approx(res);
    for (auto i = U.begin(); i!=U.end(); ++i) {
        cout << "{";
        print_points(*i);
        cout << "}\n";
    }
    start = clock();
    auto res2 = I_accel({Omega});
    end = clock();
    printf("Computed I_accel %.10f seconds\n\n", (double)(end - start)/CLOCKS_PER_SEC);

    cout << res.size() << endl;
    */

    return 0;

}
