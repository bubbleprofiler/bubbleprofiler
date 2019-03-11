Reporting bugs
==============

If you find a bug that is not listed below, please file a
bug report by opening an issue at the
[BubbleProfiler development repository].

To assist with finding and fixing the cause of a new bug,
please include in your bug report as much detail as
necessary about the issue and how to reproduce it. In
particular, you should include

 * Information about your system, including

   - The operating system and version (e.g., the output of
     `uname -a` or similar)
   - The version of CMake used (e.g., as shown by
     `cmake --version`)
   - The detected C++ compiler and version (e.g., as shown
     in the CMake output)
   - Any additional details that you think might be relevant,
     such as the version of any of the required packages if
     the bug appears to be related to a dependency problem

 * Information about your version of BubbleProfiler, including

   - The version or release number, or git commit hash
     if working with a cloned repository (e.g., the output of
     `git rev-parse HEAD`)
   - The configuration options used to build the package

 * Information about the bug itself, including

   - A short description of the bug, such as the erroneous
     behavior or output and, where possible, the expected
     correct result
   - How to reproduce the bug, which is most effectively
     done by including a short program or script that
     generates the problem

Known problems
==============

C++11 ABI compatibility
-----------------------

BubbleProfiler is written in the C++11 standard and this introduced ABI
changes to libraries. If you have libraries compiled with the old ABI
you may see undefined reference compile errors during the make step.

To resolve this one can either recompile the libraries with the new ABI
or pass the \_GLIBCXX\_USE\_CXX11\_ABI=0 flag during the cmake step:

    cmake -DCMAKE_CXX_FLAGS=-D_GLIBCXX_USE_CXX11_ABI=0 ..


[BubbleProfiler development repository]: https://github.com/bubbleprofiler/bubbleprofiler/issues
