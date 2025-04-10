#include <time.h>
#include <vector>
#include <algorithm>

#include "set_computations.hh"
#include "print.hh"

static Variable x(0);
static Variable y(1);

void tests(void) {

  /*
  vector<C_Polyhedron> res;

  C_Polyhedron C(2, EMPTY);
  C.add_generator(point(4*x + 4*y));
  C.add_generator(point(2*x + 4*y));
  C.add_generator(point(1*x + 5*y));

  C_Polyhedron P(2, EMPTY);
  P.add_generator(point(6*x + 3*y));
  P.add_generator(point(-5*x + 3*y));
  P.add_generator(point(-3*x + -4*y));
  P.add_generator(point(5*x + -3*y));

  C_Polyhedron P1(2, EMPTY);
  P1.add_generator(point(2*x + 2*y));
  P1.add_generator(point(1*x + 1*y));
  P1.add_generator(point(-5*x + 3*y));

  C_Polyhedron P2(2, EMPTY);
  P2.add_generator(point(2*x + 2*y));
  P2.add_generator(point(1*x + 1*y));
  P2.add_generator(point(5*x + -3*y));

  C_Polyhedron P3(2, EMPTY);
  P3.add_generator(point(-3*x + -4*y));
  P3.add_generator(point(1*x + 1*y));
  P3.add_generator(point(-5*x + 3*y));

  C_Polyhedron P4(2, EMPTY);
  P4.add_generator(point(2*x + 2*y));
  P4.add_generator(point(6*x + 3*y));
  P4.add_generator(point(5*x + -3*y));

  vector<C_Polyhedron> P_v;
  P_v.push_back(P1);
  P_v.push_back(P2);
  P_v.push_back(P3);
  P_v.push_back(P4);

  clock_t start = clock();
  for (int i = 0; i < 1000; i++) {
      res = regiondiff(P, P_v.begin(), P_v.end());
  }
  clock_t end = clock();

  print_points(res);

  printf("Computed region difference in %.15f seconds\n\n", (double)(end -
start)/CLOCKS_PER_SEC / 1000);

  start = clock();
  for (int i = 0; i < 1000; i++) {
      subset(P, P_v.begin(), P_v.end());
  }
  end = clock();

  printf("Checked subset in %.15f seconds\n\n", (double)(end -
start)/CLOCKS_PER_SEC / 1000);

  start = clock();
  bool t1, t2, t3, t4;
  for (int i = 0; i < 1000; i++) {
      t1 = intersects(P, P_v);
      t2 = intersects(P, res);
      t3 = intersects(P_v[1], res);
      t4 = intersects(C, res);
  }
  end = clock();

  printf("Checked intersection of P and P_v: %d\n", t1);
  printf("Checked intersection of P and P \\ P_v: %d\n", t2);
  printf("Checked intersection of P_v[1] and P \\ P_v]: %d\n", t3);
  printf("Checked intersection of P \\ P_v and C: %d\n", t4);
  printf("Checked intersection in %.15f seconds\n\n", (double)(end -
start)/CLOCKS_PER_SEC / 1000 / 4);

  start = clock();
  C_Polyhedron res2;
  for (int i = 0; i < 1000; i++) {
      res2 = translate_into(C, P);
  }
  end = clock();

  print_points(res2);

  printf("Translated C into P in %.15f seconds\n\n", (double)(end -
start)/CLOCKS_PER_SEC / 1000);

  start = clock();
  C_Polyhedron res3;
  for (int i = 0; i < 1000; i++) {
      res3 = translate_touching(C, P);
  }
  end = clock();

  print_points(res3);

  printf("Translated C touching P in %.15f seconds\n\n", (double)(end -
start)/CLOCKS_PER_SEC / 1000);

  start = clock();
  vector<C_Polyhedron> res4;
  for (int i = 0; i < 1000; i++) {
      res4 = translate_touching(C, res);
  }
  end = clock();

  print_points(res4);

  printf("Translated C touching P \\ P_v in %.15f seconds\n\n", (double)(end -
start)/CLOCKS_PER_SEC / 1000);

  start = clock();
  vector<C_Polyhedron> res5;
  for (int i = 0; i < 1000; i++) {
      res5 = translate_into(C, P_v, P);
  }
  end = clock();

  print_points(res5);

  printf("Translated C into P \\ P_v in %.15f seconds\n\n", (double)(end -
start)/CLOCKS_PER_SEC / 1000);

  start = clock();
  C_Polyhedron res6;
  for (int i = 0; i < 1000; i++) {
      res6 = (C + P);
  }
  end = clock();

  print_points(res6);

  printf("Added C to P in %.15f seconds\n\n", (double)(end -
start)/CLOCKS_PER_SEC / 1000);

  start = clock();
  bool res7;
  for (int i = 0; i < 1000; i++) {
      res7 = can_translate_into(C, P, P, P_v);
  }
  end = clock();

  printf("Translated C to P \\ P_v: %d in %.15f seconds\n\n",res7, (double)(end
- start)/CLOCKS_PER_SEC / 1000);

  ninterval_t x = {interval_t(0.1,0.2), interval_t(-0.1, 0.3)};

  C_Polyhedron res8;
  start = clock();
  for (int i = 0; i < 1000; i++) {
      res8 = i2p(x);
  }
  end = clock();
  print_points(res8);

  printf("i2p in %.10f seconds\n\n", (double)(end - start)/CLOCKS_PER_SEC /
1000);

  C_Polyhedron res9;
  start = clock();
  for (int i = 0; i < 1000; i++) {
      res9 = convexhull(P_v);
  }
  end = clock();
  print_points(res9);

  printf("Computed convex hull of P_v in %.10f seconds\n\n", (double)(end -
start)/CLOCKS_PER_SEC / 1000);
  */

    ninterval_t x1({interval_t(-4, 4), interval_t(-212, 2)});
    ninterval_t x2({interval_t(-4, 4), interval_t(2, 139)});
    ninterval_t x3({interval_t(-4, 4), interval_t(139, 293)});
    ninterval_t x4({interval_t(-3, 4), interval_t(2, 3)});
    ninterval_t x5({interval_t(-5, 4), interval_t(3, 4)});
    ninterval_t x6({interval_t(4, 14), interval_t(2, 3)});
    vector<ninterval_t> x10 = {x1,x2,x3,x4,x5,x6};

    clock_t start = clock();
    for (int i = 0; i < 1000; i++) {
        x10 = {x1,x2,x3,x4,x5,x6};
        merge(x10);
    }
    clock_t end = clock();
    print_points(x10);


    printf("Computed merge in %.10f seconds\n\n", (double)(end - start)/CLOCKS_PER_SEC / 1000);

    start = clock();
    for (int i = 0; i < 1000; i++) {
        x10 = {x1,x2,x3,x4,x5,x6};
        std::sort(x10.begin(), x10.end(), [](const ninterval_t &A, const ninterval_t &B) { return comp_intervals(A, B, 1); });
    }
    end = clock();

    print_points(x10);

    printf("Computed sort in %.10f seconds\n\n", (double)(end - start)/CLOCKS_PER_SEC / 1000);

    start = clock();
    for (int i = 0; i < 1000; i++) {
        x10 = {x1,x2,x3,x4,x5,x6};
        merge_fast(x10);
    }
    end = clock();

    print_points(x10);

    printf("Computed fast merge in %.10f seconds\n\n", (double)(end - start)/CLOCKS_PER_SEC / 1000);
}
