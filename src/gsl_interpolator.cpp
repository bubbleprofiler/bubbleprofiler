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

/**
   @file gsl_interpolator.cpp
   @brief contains implementation of the GSL_interpolator class
*/

#include "gsl_interpolator.hpp"
#include "error.hpp"
#include "raii_guard.hpp"

#include <gsl/gsl_errno.h>

#include <algorithm>
#include <string>

namespace BubbleProfiler {

namespace {

gsl_interp_accel* clone_gsl_interp_accel(const gsl_interp_accel* original)
{
   if (!original) {
      return nullptr;
   }

   gsl_interp_accel* cloned = gsl_interp_accel_alloc();

   if (cloned) {
      cloned->cache = original->cache;
      cloned->miss_count = original->miss_count;
      cloned->hit_count = original->hit_count;
   }

   return cloned;
}

gsl_spline* clone_gsl_spline(const gsl_spline* original)
{
   if (!original) {
      return nullptr;
   }

   const gsl_interp_type* t = original->interp->type;
   const auto size = original->size;

   gsl_spline* cloned = gsl_spline_alloc(t, size);

   if (!cloned) {
      return nullptr;
   }

   const auto old_handler = gsl_set_error_handler_off();
   const auto handler_guard = make_raii_guard([old_handler]() {
         gsl_set_error_handler(old_handler);
      });

   const int error = gsl_spline_init(cloned, original->x, original->y, size);

   if (error) {
      gsl_spline_free(cloned);
      return nullptr;
   }

   return cloned;
}

} // anonymous namespace

/**
 * The cubic spline used to obtain the interpolated function values is
 * computed using all of the given (x_old, y_old) values, of which there
 * must be at least two.  The vectors \c x_old and \c y_old must have
 * the same lengths.  Default parabolic boundary conditions, with unknown
 * derivatives, are used at the left and right ends of the domain.  The
 * obtained spline interpolant is evaluated at all points in \c x_new .
 */
Eigen::VectorXd interpolate_f_at(const Eigen::VectorXd& x_old,
                                 const Eigen::VectorXd& y_old,
                                 const Eigen::VectorXd& x_new)
{
   const std::size_t nx = x_old.size();
   const std::size_t ny = y_old.size();

   if (nx != ny) {
      throw Setup_error("initial x and y vectors have different lengths");
   }

   return interpolate_f_at(x_old, y_old, nx, x_new, x_new.size());
}

/**
 * The cubic spline used to obtain the interpolated function values is
 * computed using the first \c n_old values from \c x_old and \c y_old,
 * with \c n_old between 2 and the minimum of the number of values in
 * \c x_old and \c y_old .  The obtained spline interpolant is evaluated at
 * the first \c n_new points in \c x_new.
 */
Eigen::VectorXd interpolate_f_at(const Eigen::VectorXd& x_old,
                                 const Eigen::VectorXd& y_old, int n_old,
                                 const Eigen::VectorXd& x_new, int n_new)
{
   if (n_old < 2) {
      throw Setup_error("spline interpolation requires at least two points");
   }

   if (n_old > x_old.size() || n_old > y_old.size()) {
      throw Setup_error("requested number of points exceeds given number");
   }

   if (n_new > x_new.size()) {
      throw Setup_error(
         "interpolate_f_at: "
         "requested number of evaluation points exceeds given number");
   }

   const auto handler = gsl_set_error_handler_off();
   const auto handler_guard = make_raii_guard([handler]() {
         gsl_set_error_handler(handler);
      });

   gsl_interp_accel* acc = gsl_interp_accel_alloc();
   if (!acc) {
      throw Out_of_memory_error(
         "interpolate_f_at: "
         "unable to allocate memory for spline interpolation");
   }

   const auto accel_guard = make_raii_guard([acc]() {
         gsl_interp_accel_free(acc);
      });

   gsl_spline* spline = gsl_spline_alloc(gsl_interp_cspline, n_old);
   if (!spline) {
      throw Out_of_memory_error(
         "interpolate_f_at: "
         "unable to allocate memory for spline interpolation");
   }

   const auto spline_guard = make_raii_guard([spline]() {
         gsl_spline_free(spline);
      });

   int error = gsl_spline_init(spline, x_old.data(), y_old.data(), n_old);
   if (error) {
      throw Numerical_error("unable to initialize cubic spline");
   }

   Eigen::VectorXd y_new(Eigen::VectorXd::Zero(n_new));
   for (int i = 0; i < n_new; ++i) {
      error = gsl_spline_eval_e(spline, x_new(i), acc, &y_new(i));
      if (error) {
         throw Numerical_error("interpolate_f_at: "
                               "unable to evaluate spline interpolant");
      }
   }

   return y_new;
}

GSL_interpolator::GSL_interpolator()
   : accelerator(gsl_interp_accel_alloc())
{
   if (!accelerator) {
      throw Setup_error("GSL_interpolator::GSL_interpolator: "
                        "unable to allocate memory for interpolator");
   }
}

GSL_interpolator::GSL_interpolator(const GSL_interpolator& other)
{
   interp_type = other.interp_type;

   accelerator.reset(clone_gsl_interp_accel(other.accelerator.get()));
   if (other.accelerator && !accelerator) {
      throw Setup_error("GSL_interpolator::GSL_interpolator: "
                        "unable to allocate memory for interpolator");
   }

   interpolant.reset(clone_gsl_spline(other.interpolant.get()));
   if (other.interpolant && !interpolant) {
      throw Setup_error("GSL_interpolator::GSL_interpolator: "
                        "unable to allocate memory for interpolator");
   }
}

GSL_interpolator& GSL_interpolator::operator=(const GSL_interpolator& other)
{
   if (&other == this) {
      return *this;
   }

   interp_type = other.interp_type;

   accelerator.reset(clone_gsl_interp_accel(other.accelerator.get()));
   if (other.accelerator && !accelerator) {
      throw Setup_error("GSL_interpolator::GSL_interpolator: "
                        "unable to allocate memory for interpolator");
   }

   interpolant.reset(clone_gsl_spline(other.interpolant.get()));
   if (other.interpolant && !interpolant) {
      throw Setup_error("GSL_interpolator::GSL_interpolator: "
                        "unable to allocate memory for interpolator");
   }

   return *this;
}

/**
 * The interpolant is computed using all of the given (x, y) values.
 * The minimum required number of points depends on the choice of
 * interpolation method.  The vectors x and y must have
 * the same lengths.  The resulting interpolant can subsequently
 * be evaluated by calls to \c evaluate_at , \c evaluate_deriv_at ,
 * and \c evaluate_second_deriv_at .
 */
void GSL_interpolator::construct_interpolant(const Eigen::VectorXd& x,
                                             const Eigen::VectorXd& y)
{
   const auto nx = x.size();
   const auto ny = y.size();

   if (nx != ny) {
      throw Setup_error("GSL_interpolator::construct_interpolant: "
                        "number of x values does not match number of y values");
   }

   construct_interpolant(x, y, nx);
}

/**
 * The requested number of points must be at least the minimum required
 * by the chosen interpolation method, and less than the
 * minimum of the number of x values or number of y values.  The first
 * \c n points are then used to construct the interpolant.  The resulting
 * interpolant can subsequently be evaluated by calls to
 * \c evaluate_at , \c evaluate_deriv_at , and \c evaluate_second_deriv_at .
 */
void GSL_interpolator::construct_interpolant(const Eigen::VectorXd& x,
                                             const Eigen::VectorXd& y,
                                             int n)
{
   const gsl_interp_type* t = to_gsl_interp_type(interp_type);
   const int min_n = gsl_interp_type_min_size(t);

   if (n < min_n) {
      const std::string name(t->name);
      std::string msg("GSL_interpolator::construct_interpolant: "
                      "interpolation with interpolator '");
      msg += name;
      msg += "' requires at least " + std::to_string(min_n)
         + " points (" + std::to_string(n) + " given";
      throw Setup_error(msg);
   }

   if (n > x.size() || n > y.size()) {
      throw Setup_error("GSL_interpolator::construct_interpolant: "
                        "requested number of points exceeds given number");
   }

   reset(interp_type, n);

   const int error = gsl_spline_init(interpolant.get(), x.data(), y.data(), n);
   if (error) {
      throw Setup_error("GSL_interpolator::construct_interpolant: "
                        "unable to initialize interpolant");
   }
}

/**
 * The interpolant must previously have been computed
 * via a call to \c construct_interpolant .
 */
double GSL_interpolator::evaluate_at(double x) const
{
   check_interpolant_computed();
   check_within_bounds(x);

   const auto handler = gsl_set_error_handler_off();
   const auto handler_guard = make_raii_guard([handler]() {
         gsl_set_error_handler(handler);
      });

   double result = 0.;
   const int error = gsl_spline_eval_e(interpolant.get(), x, accelerator.get(),
                                       &result);

   if (error) {
      throw Numerical_error("GSL_interpolator::evaluate_at: "
                            "unable to evaluate interpolant ("
                            + std::to_string(error) + ")");
   }

   return result;
}

/**
 * The interpolant must previously have been computed
 * via a call to \c construct_interpolant .
 */
double GSL_interpolator::evaluate_deriv_at(double x) const
{
   check_interpolant_computed();
   check_within_bounds(x);

   const auto handler = gsl_set_error_handler_off();
   const auto handler_guard = make_raii_guard([handler]() {
         gsl_set_error_handler(handler);
      });

   double result = 0.;
   const int error = gsl_spline_eval_deriv_e(interpolant.get(), x,
                                             accelerator.get(), &result);

   if (error) {
      throw Numerical_error("GSL_interpolator::evaluate_deriv_at: "
                            "unable to evaluate interpolant");
   }

   return result;
}

/**
 * The interpolant must previously have been computed
 * via a call to \c construct_interpolant .
 */
double GSL_interpolator::evaluate_second_deriv_at(double x) const
{
   check_interpolant_computed();
   check_within_bounds(x);

   const auto handler = gsl_set_error_handler_off();
   const auto handler_guard = make_raii_guard([handler]() {
         gsl_set_error_handler(handler);
      });

   double result = 0.;
   const int error = gsl_spline_eval_deriv2_e(interpolant.get(), x,
                                              accelerator.get(), &result);

   if (error) {
      throw Numerical_error("GSL_interpolator::evaluate_second_deriv_at: "
                            "unable to evaluate interpolant");
   }

   return result;
}

const gsl_interp_type* GSL_interpolator::to_gsl_interp_type(
   Interpolation_type t) const
{
   switch (t) {
   case Interpolation_type::GSL_linear: {
      return gsl_interp_linear;
   }
   case Interpolation_type::GSL_polynomial: {
      return gsl_interp_polynomial;
   }
   case Interpolation_type::GSL_cspline: {
      return gsl_interp_cspline;
   }
   case Interpolation_type::GSL_cspline_periodic: {
      return gsl_interp_cspline_periodic;
   }
   case Interpolation_type::GSL_akima: {
      return gsl_interp_akima;
   }
   case Interpolation_type::GSL_akima_periodic: {
      return gsl_interp_akima_periodic;
   }
#ifdef GSL_MAJOR_VERSION
#if (GSL_MAJOR_VERSION >= 2)
   case Interpolation_type::GSL_steffen: {
      return gsl_interp_steffen;
   }
#endif
#endif
   default:
      throw Setup_error("GSL_interpolator::to_gsl_interp_type: "
                        "unrecognized interpolator type");
   }

   return nullptr;
}

void GSL_interpolator::reset(Interpolation_type t,
                             int n_points)
{
   gsl_interp_accel_reset(accelerator.get());

   interpolant.reset(gsl_spline_alloc(to_gsl_interp_type(t), n_points));
   if (!interpolant) {
      throw Out_of_memory_error("GSL_interpolator::construct_interpolant: "
                                "unable to allocate memory for interpolator");
   }
}

void GSL_interpolator::check_interpolant_computed() const
{
   if (!interpolant) {
      throw Setup_error("GSL_interpolator::evaluate_second_deriv_at: "
                        "interpolant not computed");
   }
}

void GSL_interpolator::check_within_bounds(double x) const
{
   const double x_min = interpolant.get()->interp->xmin;
   const double x_max = interpolant.get()->interp->xmax;

   if (x < x_min || x > x_max) {
      throw Domain_error("GSL_interpolator::check_within_bounds: "
                         "cannot interpolate outside of domain");
   }
}

} // namespace BubbleProfiler
