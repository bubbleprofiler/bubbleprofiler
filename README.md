BubbleProfiler
==============

BubbleProfiler is a C++ library for calculating the bounce solution
and Euclidean action associated with a cosmological phase transition.

Quickstart
==========

Requirements
------------

The following dependencies are required to build BubbleProfiler:

 * [CMake] (version 2.8.12 or higher)
 * C++11 compatible compiler (g++ >= 4.8.5 or clang++ >= 3.3 or icpc >= 15.0)
 * [Boost] (version 1.53.0 or higher), specifically the components
   - Boost.Program\_options
   - Boost.Filesystem
   - Boost.System
 * [Eigen 3] (version 3.1.0 or higher)
 * [GiNaC] (version 1.6.2 or higher)
 * [GNU scientific library] (version 1.15 or higher)
 * [NLopt] (version 2.4.1 or higher)

The following dependencies are optional:

 * [Doxygen], for building the documentation
 * Python 2.7.x or Python 3.x, for running some examples

Building BubbleProfiler
-----------------------

BubbleProfiler uses CMake to configure its build system. To build the
library and command line executable,

1. Create a build directory in which to build the package:

        mkdir build

2. Run CMake from this directory to configure the build system:

        cd build
        cmake ..

3. Compile the core library and command line executable:

        cmake --build .

4. Optionally, build and run the unit tests:

        cmake --build . --target 'check'

5. Optionally, build the example programs:

        cmake --build . --target 'examples'

6. Optionally, build the package documentation:

        cmake --build . --target 'docs'

7. Optionally, install the library and executables:

        cmake --build . --target 'install'

After building the package, the package directory will contain the
`bin/` and `lib/` subdirectories containing the main command line
executable and the BubbleProfiler library. For more detailed instructions
and help with the build and installation process, refer to the accompanying
INSTALL.md file.

Running BubbleProfiler
----------------------

BubbleProfiler consists of two core components, a command line executable
for standalone use and a library to be used in your own programs. The
executable `run_cmd_line_potential.x` calculates the action and, optionally,
the field profiles for the bounce solution, given a potential and set of
fields. For example, the bounce action for a potential involving two scalar
fields may be calculated using

    ./bin/run_cmd_line_potential.x \
         --potential "(x^2 + y^2)*(1.8*(x - 1)^2 + 0.2*(y-1)^2 - 0.4)" \
         --field x --field y --local-minimum 0 0

A summary of the possible options for the program may be printed using
`./bin/run_cmd_line_potential.x --help`. Detailed information about
the program and its options may be found in the accompanying manual.

Package content
===============

The following subdirectories are included within the BubbleProfiler
source distribution:

 * `cmake/` contains modules and scripts used by the CMake build system
 * `docs/` contains source files required for building the documentation
 * `examples/` contains examples demonstrating the use of BubbleProfiler
 * `include/` contains the C++ header files provided by the package
 * `src/` contains the C++ source files from which the package is built
 * `test/` contains the package unit tests
 * `utils/` contains helpful utility scripts, primarily for the package
   developers

Additional information regarding the purpose and content of
individual files in the package may be found in the documentation.

License
=======

BubbleProfiler is licensed under the GNU General Public License, version 3.
See the accompanying LICENSE file for details.


[Boost]: http://www.boost.org
[CMake]: https://cmake.org
[Doxygen]: http://www.doxygen.nl
[Eigen 3]: http://eigen.tuxfamily.org
[GiNaC]: http://www.ginac.de
[GNU scientific library]: http://www.gnu.org/software/gsl
[NLopt]: https://nlopt.readthedocs.io
