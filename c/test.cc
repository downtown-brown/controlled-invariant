/* Boost example/rational.cpp
 * example program of how to use interval< rational<> >
 *
 * Copyright 2002-2003 Guillaume Melquiond, Sylvain Pion
 *
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or
 * copy at http://www.boost.org/LICENSE_1_0.txt)
 */

// it would have been enough to only include:
//   <boost/numeric/interval.hpp>
// but it's a bit overkill to include processor intrinsics
// and transcendental functions, so we do it by ourselves
#include <cstdint>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>


#include <boost/numeric/interval.hpp>      // base class

#include <boost/numeric/interval/interval.hpp>
#include <boost/rational.hpp>
#include <iostream>

using namespace boost::numeric;
using namespace interval_lib;
typedef interval<double, policies<save_state<rounded_transc_std<double> >,
                                  checking_base<double> > > Interval;
void rat_approx(double f, int64_t md, int64_t *num, int64_t *denom);

std::ostream& operator<<(std::ostream& os, const Interval& r) {
    os << "[" << r.lower() << "," << r.upper() << "]";
    return os;
}

int main() {
    double p = 2./3.;
    double q = 4./5.;
    Interval z(4., 5.);
    Interval a(2., 3.);
    Interval b(0.9,0.95);
    a += z;
    z *= q;
    a -= p;
    a /= q;
    std::cout << z << std::endl;
    std::cout << a << std::endl;

    int64_t num, den;
    auto start = clock();
    for (int i = 0; i < 1000; i++) {
        rat_approx(sin(b).lower(), INT16_MAX, &num, &den);
    }
    auto end = clock();

    std::cout << sin(b).lower() << " = " << num << "/" << den << std::endl;
    printf("Computed in %.6f ms\n\n", (double)(end - start)/CLOCKS_PER_SEC);
}

void rat_approx(double f, int64_t md, int64_t *num, int64_t *denom)
{
	/*  a: continued fraction coefficients. */
	int64_t a, h[3] = { 0, 1, 0 }, k[3] = { 1, 0, 0 };
	int64_t x, d, n = 1;
	int i, neg = 0;

	if (md <= 1) { *denom = 1; *num = (int64_t) f; return; }

	if (f < 0) { neg = 1; f = -f; }

	while (f != floor(f)) { n <<= 1; f *= 2; }
	d = f;

	/* continued fraction and check denominator each step */
	for (i = 0; i < 64; i++) {
		a = n ? d / n : 0;
		if (i && !a) break;

		x = d; d = n; n = x % n;

		x = a;
		if (k[1] * a + k[0] >= md) {
			x = (md - k[0]) / k[1];
			if (x * 2 >= a || k[1] >= md)
				i = 65;
			else
				break;
		}

		h[2] = x * h[1] + h[0]; h[0] = h[1]; h[1] = h[2];
		k[2] = x * k[1] + k[0]; k[0] = k[1]; k[1] = k[2];
	}
	*denom = k[1];
	*num = neg ? -h[1] : h[1];
}
