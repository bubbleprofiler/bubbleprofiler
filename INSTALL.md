Building BubbleProfiler
=======================

Basic build steps
-----------------
BubbleProfiler uses CMake to generate a build system appropriate for
the target platform. On all systems, the basic build process may be
summarized as

 * Ensure all required dependencies are installed
 * Run CMake in the desired build location to generate a custom
   build system, e.g., a set of Makefiles on GNU/Linux systems
   with Make installed
 * Run said build system to build the core library and executables
 * Optionally, compile additional components and install the
   build outputs to a different location

When building BubbleProfiler, it is recommended that an
'out-of-source' build be performed, in which the build process takes
place in a different directory from the main BubbleProfiler source
directory. By default, the build directory is the directory from
which CMake is invoked. To perform an out-of-source build, first
navigate to the desired build directory, and then invoke CMake,
passing the path to the BubbleProfiler source directory (i.e.,
the directory containing this file) as an argument. For example,
on a GNU/Linux system, one would issue the shell commands

    cd /path/to/desired/build/directory
    cmake /path/to/BubbleProfiler/source/directory

At this point, CMake will attempt to detect a supported C++ compiler
and the required dependencies for the package, and generate the
build system. If this process succeeds, the main library and
command line executable may then be compiled by running the generated
build system,

    cmake --build .

The `--build` option instructs CMake to run the underlying build
system.  Alternatively, the build system may be run directly,
but the exact command to do so depends on the platform
and chosen build system. For example, on a GNU/Linux system with Make
installed, the library and command line program are built by typing

    make

within the build directory. Successful compilation of the package
will result in the BubbleProfiler library being created in the
`lib/` subdirectory of the main package directory, and the command
line executable `run_cmd_line_potential.x` will be created in the
`bin/` subdirectory of the top level directory.

After the core program components have been successfully compiled,
additional steps can be taken to test the resulting library and/or
install it in an alternative location. To check that everything
is working as expected, the unit test suite may be built and run
by building the custom `check` target, e.g.,

    cmake --build . --target 'check'

or, e.g., if the underlying build system is Make,

    make check

will execute the unit tests. Once you are satisfied that the
package is working correctly, the library and command line
executable may also be installed to a different location using
the `install` target, e.g.,

    cmake --build . --target 'install'

or directly using

    make install

when building using Make.

Additional build options
------------------------

The build process may be customized by passing options to CMake.
If required packages are installed in non-standard locations,
the paths to these packages may be passed to CMake. For example,
a custom installation of Boost may be detected by issuing the
CMake command

    cmake -DBOOST_INCLUDEDIR=/path/to/boost/headers \
           -DBOOST_LIBRARYDIR=/path/to/boost/libraries \
           -DBoost_NO_SYSTEM_PATHS:BOOL=OFF \
           /path/to/BubbleProfiler/source/directory

Please refer to the documentation of the individual CMake
modules for further details on how to provide alternative
paths to search for a given package.

In the basic build process described above, the BubbleProfiler
library and executables are immediately generated into the
`lib/` and `bin/` subdirectories of the main package
directory. A different location for these build outputs can
be set by providing the desired values of `CMAKE_RUNTIME_OUTPUT_DIRECTORY`,
`CMAKE_ARCHIVE_OUTPUT_DIRECTORY`, and `CMAKE_LIBRARY_OUTPUT_DIRECTORY`
as options to CMake.

By default, the `install` target installs the package to
the appropriate subdirectories of `/usr/local` on UNIX systems.
An alternative directory may be chosen by specifying the value
of `CMAKE_INSTALL_PREFIX`, e.g.,

    cmake -DCMAKE_INSTALL_PREFIX=/desired/install/directory /path/to/source

The default behavior is for the command line executable to beinstalled to the
location `CMAKE_INSTALL_PREFIX/bin`, the compiled library to
`CMAKE_INSTALL_PREFIX/lib`, and the package headers to
`CMAKE_INSTALL_PREFIX/include`. Finer control over the locations
in which these components are installed can be achieved by setting the
values of the variables `BUBBLEPROFILER_INSTALL_RUNTIME_DIR`,
`BUBBLEPROFILER_INSTALL_ARCHIVE_DIR`, and `BUBBLEPROFILER_INSTALL_INCLUDE_DIR`
to the desired installation locations of the executable, library, and headers
relative to the path given by `CMAKE_INSTALL_PREFIX`. For example, invoking
CMake with the variable settings

    cmake -DCMAKE_INSTALL_PREFIX="/home/<user>/BubbleProfiler" \
           -DBUBBLEPROFILER_INSTALL_RUNTIME_DIR="bin" \
           -DBUBBLEPROFILER_INSTALL_INCLUDE_DIR="include/BubbleProfiler" \
           /path/to/source

will, upon running the `install` target, install the command line
executable into `/home/<user>/BubbleProfiler/bin/`, the package
headers into `/home/<user>/BubbleProfiler/include/BubbleProfiler/`,
and the library into `/home/<user>/BubbleProfiler/lib`.

By default, the build system produces a `Release` build of BubbleProfiler.
If debugging symbols are required, a `Debug` build may be requested
by setting the value of the `CMAKE_BUILD_TYPE` variable, e.g.,

    cmake -DCMAKE_BUILD_TYPE="Debug" ..

Additionally, compiler warnings may be enabled by setting the variable
`ENABLE_COMPILER_WARNINGS` to `ON`, e.g.,

    cmake -DENABLE_COMPILER_WARNINGS:BOOL=ON ..

It is recommended that developers contributing to BubbleProfiler always
enable this option.

Building the examples
---------------------

The BubbleProfiler distribution includes several examples demonstrating
the use of the library, located in the `examples/` directory. All of the
provided examples may be built at once by building the `examples` target
in the build directory, e.g.,

    cmake --build . --target 'examples'

or, for a Make build system,

    make examples

This will produce a set of example executables in the `bin/` subdirectory
of the package directory. Alternatively, individual example programs may be
built by specifying the desired example target. The possible targets are

 * `gaussian`
 * `fubini`
 * `logarithmic`
 * `perturbative`
 * `quartic`
 * `quartic_tabulate`
 * `scale`
 * `thin`

Please refer to the source for individual examples for information on
their purpose and usage.

Building the documentation
--------------------------

Documentation for the BubbleProfiler library is generated from inline
documentation in the source code using Doxygen. To build this documentation,
first ensure that Doxygen is installed on your system and is located
by CMake when generating the build system. The HTML documentation may
then be created by building the target `docs` from the build directory,

    cmake --build . --target 'docs'

or

    make docs

when using Make. The generated documentation may be found in the
`docs/` subdirectory of the build directory.

Platform specific build instructions
====================================

Additional details and hints on installing required packages and
building BubbleProfiler for some specific operating systems and
build systems are provided below.

Building on macOS
-----------------

To build on macOS, it is recommended to use [homebrew](https://brew.sh)
to manage the installation of the package dependencies. These can be
installed via

    brew install cmake
    brew install homebrew/versions/llvm38 # this is preferred to gcc-6
    brew install install boost
    brew install ginac gsl nlopt eigen # cln should be pulled in by ginac

BubbleProfiler may then be compiled by navigating to the desired build
directory and running CMake, where we direct CMake to use the C and
C++ compilers installed using homebrew,

    cd /path/to/build/directory
    CC=clang CXX=clang++ cmake /path/to/source/directory

The main package components may then be built by typing

    make

Building on GNU/Linux
---------------------

Depending on your distribution, BubbleProfiler's dependencies may be
installable via the distribution's package manager.

On Ubuntu, the required packages may be installed using

    sudo apt install cmake cmake-extras \
                       libboost-all-dev \
                       libeigen3-dev \
                       libginac-dev \
                       libnlopt-dev \
                       libgsl-dev

You can trim `libboost-all-dev` to the specific set of headers and
libraries needed, if you need to minimize the dependency footprint.

On Fedora, the necessary packages are installed using

    sudo dnf install cmake \
                       boost-devel \
                       eigen3-devel \
                       ginac-devel \
                       NLopt-devel \
                       gsl-devel

Once all dependencies are installed, the package is built by first running
CMake in the desired build directory, and then running the generated build
system. On a system with Make installed, for example, this is achieved by

    cd /path/to/build/directory
    cmake /path/to/source/directory
    make
