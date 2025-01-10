#include <time.h>
#include <stdio.h>

#include "invariant.hh"

void tests(void);

using namespace std;

uint8_t stop = 0;

int main() {
    //tests();

    IntervalData Omega({interval_t(-6, 6), interval_t(-6, 6)});

    auto start = clock();
    vector<IntervalData> res = {Omega};
    uint64_t i = 0;
    while (!stop) {
        res = I_approx(res);
        i++;
    }
    auto end = clock();
    printf("Computed %ld iters I_approx %.10f seconds\n\n", i,
           (double)(end - start)/CLOCKS_PER_SEC);

    /*
    cout << "inputs:\n";
    auto U = U_approx(res);
    for (auto i = U.begin(); i!=U.end(); ++i) {
        cout << "{";
        print_points(*i);
        cout << "}\n";
    }
    */

    return 0;

}
