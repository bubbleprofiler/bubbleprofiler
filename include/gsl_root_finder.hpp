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

#ifndef BUBBLEPROFILER_GSL_ROOT_FINDER_HPP_INCLUDED
#define BUBBLEPROFILER_GSL_ROOT_FINDER_HPP_INCLUDED

#include "error.hpp"
#include "raii_guard.hpp"
#include "root_finder.hpp"

#include <Eigen/Core>
#include <gsl/gsl_multiroots.h>

#include <cmath>

namespace BubbleProfiler {

template <int N>
class GSL_root_finder : public Root_finder<Eigen::Matrix<double,N,1> > {
public:
   using Function_argument = Eigen::Matrix<double,N,1>;
   using Function_value = Eigen::Matrix<double,N,1>;

   enum class Algorithm {
      GSL_hybrids, GSL_hybrid, GSL_dnewton
   };

   virtual ~GSL_root_finder() = default;

   virtual const Function_argument& get_solution() const override {
      return solution;
   }

   void set_algorithm(Algorithm a) { solver_algorithm = a; }
   void set_relative_error(double e) { rel_err = e; }
   void set_absolute_error(double e) { abs_err = e; }
   void set_use_stepsize_convergence_check(bool flag) {
      use_stepsize_convergence = flag;
   }
   void set_max_iterations(std::size_t i) { max_iterations = i; }

protected:
   using Function = std::function<Function_value(const Function_argument&)>;

   virtual Root_finder_status solve(const Function&, const Function_argument&) override;

private:
   std::size_t n_dims;
   Algorithm solver_algorithm{Algorithm::GSL_hybrids};
   double rel_err{1.e-6};
   double abs_err{1.e-6};
   bool use_stepsize_convergence{false};
   std::size_t max_iterations{100};
   Function function{nullptr};
   Function_argument solution{};

   const gsl_multiroot_fsolver_type* get_gsl_solver_type() const;
   static int gsl_function(const gsl_vector*, void*, gsl_vector*);
};

template <int N>
int GSL_root_finder<N>::gsl_function(
   const gsl_vector* x, void* parameters, gsl_vector* f)
{
   const std::size_t length = x->size;
   for (std::size_t i = 0; i < length; ++i) {
      if (!std::isfinite(gsl_vector_get(x, i))) {
         gsl_vector_set_all(f, std::numeric_limits<double>::max());
         return GSL_EDOM;
      }
   }

   Function* func = static_cast<Function*>(parameters);
   Function_argument arg(Function_argument::Map(x->data, length));
   Function_value result(length);
   result.setConstant(std::numeric_limits<double>::max());

   int status = GSL_SUCCESS;
   try {
      result = (*func)(arg);
      bool finite_result = true;
      for (std::size_t i = 0; i < length; ++i) {
         finite_result = finite_result && std::isfinite(result(i));
      }
      status = finite_result ? GSL_SUCCESS : GSL_EDOM;
   } catch (const Error&) {
      status = GSL_EDOM;
   }

   for (std::size_t i = 0; i < length; ++i) {
      gsl_vector_set(f, i, result(i));
   }

   return status;
}

template <int N>
Root_finder_status GSL_root_finder<N>::solve(
   const Function& f, const Function_argument& initial_guess)
{
   if (!f) {
      throw Setup_error("GSL_root_finder::solve: "
                        "invalid function provided");
   }

   n_dims = initial_guess.size();
   function = f;
   solution = initial_guess;

   void* parameters = &function;
   gsl_multiroot_function func = {gsl_function, n_dims, parameters};

   gsl_multiroot_fsolver* solver
      = gsl_multiroot_fsolver_alloc(get_gsl_solver_type(), n_dims);

   if (!solver) {
      throw Out_of_memory_error("cannot allocate root finder");
   }

   const auto handler = gsl_set_error_handler_off();

   const auto solver_guard = make_raii_guard(
      [solver, handler]() {
         gsl_multiroot_fsolver_free(solver);
         gsl_set_error_handler(handler);
      });

   gsl_vector* guess = gsl_vector_alloc(n_dims);

   if (!guess) {
      throw Out_of_memory_error("cannot allocate solution guess");
   }

   const auto vector_guard = make_raii_guard(
      [guess]() {
         gsl_vector_free(guess);
      });

   for (std::size_t i = 0; i < n_dims; ++i) {
      gsl_vector_set(guess, i, initial_guess(i));
   }

   gsl_multiroot_fsolver_set(solver, &func, guess);

   std::size_t iterations = 0;
   int status;
   do {
      iterations++;

      status = gsl_multiroot_fsolver_iterate(solver);

      if (status) {
         break;
      }

      if (!use_stepsize_convergence) {
         status = gsl_multiroot_test_residual(solver->f, abs_err);
      } else {
         status = gsl_multiroot_test_delta(solver->dx, solver->x,
                                           abs_err, rel_err);
      }
   } while (status == GSL_CONTINUE && iterations < max_iterations);

   for (std::size_t i = 0; i < n_dims; ++i) {
      solution(i) = gsl_vector_get(solver->x, i);
   }

   return (status == GSL_SUCCESS ? Root_finder_status::SUCCESS
           : Root_finder_status::FAIL);
}

template <int N>
const gsl_multiroot_fsolver_type* GSL_root_finder<N>::get_gsl_solver_type() const
{
   switch (solver_algorithm) {
   case Algorithm::GSL_hybrids: return gsl_multiroot_fsolver_hybrids;
   case Algorithm::GSL_hybrid: return gsl_multiroot_fsolver_hybrid;
   case Algorithm::GSL_dnewton: return gsl_multiroot_fsolver_dnewton;
   default:
      throw Setup_error("GSL_root_finder::get_gsl_solver: "
                        "unrecognized solver type");
   }

   return nullptr;
}

} // namespace BubbleProfiler

#endif
