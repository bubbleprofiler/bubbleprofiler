/*
 * This file is part of BubbleProfiler.
 *
 * BubbleProfiler is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * BubbleProfiler is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with BubbleProfiler.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef BUBBLEPROFILER_ERROR_HPP_INCLUDED
#define BUBBLEPROFILER_ERROR_HPP_INCLUDED

#include <stdexcept>
#include <string>

namespace BubbleProfiler {

class Error : public std::runtime_error {
public:
   using std::runtime_error::runtime_error;
};

//! Exception indicating general setup error
class Setup_error : public Error {
public:
   using Error::Error;
};

//! Exception indicating a thin-wall ansatz which we can't solve (yet!)
class Thin_wall_error : public Error {
public:
    using Error::Error;
};

//! Exception indicating failure to integrate BVP
class BVP_solver_error : public Error {
public:
   using Error::Error;

   BVP_solver_error()
      : Error("BVP_solver_error: failed to solve BVP")
      {}
};

//! Exception indicating failure of iterative procedure to converge
class No_convergence_error : public Error {
public:
   explicit No_convergence_error(int it_)
      : Error("No_convergence_error: no convergence after "
              + std::to_string(it_) + " iterations")
      , iterations(it_)
      {}
   No_convergence_error(int it_, const std::string& what_)
      : Error(what_), iterations(it_)
      {}
   No_convergence_error(int it_, const char* what_)
      : Error(what_), iterations(it_)
      {}

   int get_number_of_iterations() const { return iterations; }

private:
   int iterations{0};
};

//! Exception indicating out of bounds access error
class Out_of_bounds_error : public Error {
public:
   explicit Out_of_bounds_error(int idx_)
      : Error("Out_of_bounds_error: attempted to access out of bounds index "
              + std::to_string(idx_))
      , index(idx_)
      {}
   Out_of_bounds_error(int idx_, const std::string& what_)
      : Error(what_)
      , index(idx_)
      {}
   Out_of_bounds_error(int idx_, const char* what_)
      : Error(what_)
      , index(idx_)
      {}

   int get_index_value() const { return index; }

private:
   int index{0};
};

//! Exception indicating that memory allocation failed
class Out_of_memory_error : public Error {
public:
   using Error::Error;
};

//! Exception indicating that the optimizer failed to converge
class Optimizer_error : public Error {
public:
   using Error::Error;
};

//! Exception indicating function evaluation out of allowed domain
class Domain_error : public Error {
public:
   using Error::Error;
};

//! Exception indicating generic numerical error
class Numerical_error : public Error {
public:
   using Error::Error;
};

//! Exception indicating generic input/output error
class IO_error : public Error {
public:
   using Error::Error;
};

/*!
 * Exception indicating that a non implemented method has been called.
 * Mainly for use in defining mock classes.
 */
class Not_implemented_error : public Error {
public:
   using Error::Error;
};

} // namespace BubbleProfiler

#endif
