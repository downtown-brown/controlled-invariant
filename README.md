# Computing controlled invariant sets
This repository contains the code to compute the controlled invariant
sets of a discrete-time nonlinear system.

## Dependencies
### Boost
The interval analysis is performed by the boost/numeric/interval
library. It has been tested with boost 1.86 but should work with older
or newer versions. Install this from your system's package manager and
make sure the header files are in your include path.

### GNU Multiple Precision library (GMP)
The PPL library uses GMP integers for its computations. Install `gmp`
from your package manager or from the [website](gmplib.org). Again,
make sure the library and headers can be found by the compiler.

### Parma Polyhedra Library (PPL)
The polytope computations are performed using PPL. In order to take
advantage of multithreading, the `devel` branch must be used. The main
release is not thread safe; using it will usually result in a
segfault, and will definitely make the calculations wrong.

To simplify building and checking out the library, it has been
included as a submodule. It can be downloaded with

``` shell
git submodule init
git submodule update
```

#### Building PPL
See instructions inside the `ppl` folder (submodule). There are
several options you can enable or disable during the
configuration. You must use the `enable-thread-safe` option.
The following will build and install the library.

```shell
autoreconf -fi
./configure --enable-coefficients=mpz --enable-optimization=sspeed --enable-thread-safe=yes
make -j
make install      # might need sudo here
```

Sometimes the last step will hang while generating the doxygen
files. You can interrupt it at this step and the library installation
will still work.  If this succeeds you should have `ppl.hh` in one of
the directories on your include path (typically `/usr/local/include`).

## Building the program
Simply run `make` in the top level after building and installing the
dependencies. You should have an executable named `invariant`.

## Changing the system and parameters
The system models are defined in header files in the `src/models`
directory. The model must define several parameters and functions.

- `U`: The interval/polytope of acceptable control inputs
- `Omega_0`: The initial state domain
- `A(x_m, P)`: A function returning a polytope which is the image of
  `P` under the transformation matrix `A`, which is the linear part of
  the dynamics after decomposition. `A` may be parameterized
  by `x_m` which is the midpoint of the current interval.
- `B(x_m, U)`: A function returning a polytope which is the image of
  `U` under the matrix `B`, which is the linear part of the input
  dependence after decomposition. `B` may be parameterized by `x_m`
  which is the midpoint of the current interval.
- `Phi(x, x_m)`: The nonlinear function of the state after decomposition.
- `Psi(x, x_m, U)` The nonlinear function of the input after decompostion.

The functions `A` and `B` are somewhat difficult to write due to the
idiosyncrasies of the PPL interface. Future plans include a more
user-friendly interface for defining new models.
