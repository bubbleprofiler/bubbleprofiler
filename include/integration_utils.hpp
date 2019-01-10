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

#ifndef BUBBLEPROFILER_INTEGRATION_UTILS_HPP_INCLUDED
#define BUBBLEPROFILER_INTEGRATION_UTILS_HPP_INCLUDED

#include "error.hpp"
#include "raii_guard.hpp"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

#include <tuple>
#include <iostream>
namespace BubbleProfiler {

enum class Integration_rule {
   GK15, GK21, GK31, GK41, GK51, GK61
};

int get_integration_rule_key(Integration_rule);

template <typename F>
std::tuple<int,double,double> integrate_gsl_qag(F f, double a, double b,
                                                double abs_tol, double rel_tol,
                                                std::size_t max_intervals,
                                                Integration_rule rule)
{
   gsl_integration_workspace* workspace
      = gsl_integration_workspace_alloc(max_intervals);

   if (!workspace) {
      throw Out_of_memory_error("cannot allocate integration workspace");
   }

   const auto handler = gsl_set_error_handler_off();

   const auto workspace_guard = make_raii_guard(
      [workspace, handler]() {
         gsl_integration_workspace_free(workspace);
         gsl_set_error_handler(handler);
      });

   gsl_function func;
   func.function = [](double x, void* params) {
      return (*static_cast<F*>(params))(x);
   };
   func.params = &f;

   double result = 0.;
   double error = 0.;
   const int status = gsl_integration_qag(&func, a, b, abs_tol, rel_tol,
                                          max_intervals,
                                          get_integration_rule_key(rule),
                                          workspace, &result, &error);

   return std::make_tuple(status, result, error);
}

template <typename F>
std::tuple<int,double,double> integrate_gsl_qags(F f, double a, double b,
                                                 double abs_tol, double rel_tol,
                                                 std::size_t max_intervals)
{
   gsl_integration_workspace* workspace
      = gsl_integration_workspace_alloc(max_intervals);

   if (!workspace) {
      throw Out_of_memory_error("cannot allocate integration workspace");
   }

   const auto handler = gsl_set_error_handler_off();

   const auto workspace_guard = make_raii_guard(
      [workspace, handler]() {
         gsl_integration_workspace_free(workspace);
         gsl_set_error_handler(handler);
      });

   gsl_function func;
   func.function = [](double x, void* params) {
      return (*static_cast<F*>(params))(x);
   };
   func.params = &f;

   double result = 0.;
   double error = 0.;
   const int status = gsl_integration_qags(&func, a, b, abs_tol, rel_tol,
                                           max_intervals,
                                           workspace, &result, &error);

   return std::make_tuple(status, result, error);
}

template <typename F>
std::tuple<int,double,double,int> integrate_gsl_cquad(F f, double a, double b,
                                                      double abs_tol,
                                                      double rel_tol,
                                                      std::size_t n_intervals)
{
   if (n_intervals < 3) {
      throw Setup_error("at least 3 intervals required in CQUAD routine");
   }

   gsl_integration_cquad_workspace* workspace
      = gsl_integration_cquad_workspace_alloc(n_intervals);

   if (!workspace) {
      throw Out_of_memory_error("cannot allocate integration workspace");
   }

   const auto handler = gsl_set_error_handler_off();

   const auto workspace_guard = make_raii_guard(
      [workspace, handler]() {
         gsl_integration_cquad_workspace_free(workspace);
         gsl_set_error_handler(handler);
      });

      gsl_function func;
   func.function = [](double x, void* params) {
      return (*static_cast<F*>(params))(x);
   };
   func.params = &f;

   double result = 0.;
   double error = 0.;
   std::size_t n_evals = 0;
   const int status = gsl_integration_cquad(&func, a, b, abs_tol, rel_tol,
                                            workspace, &result, &error,
                                            &n_evals);

   return std::make_tuple(status, result, error, n_evals);
}

} // namespace BubbleProfiler

#endif
