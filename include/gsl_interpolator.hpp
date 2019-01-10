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

#ifndef BUBBLEPROFILER_GSL_INTERPOLATOR_HPP_INCLUDED
#define BUBBLEPROFILER_GSL_INTERPOLATOR_HPP_INCLUDED

/**
 * @file gsl_interpolator.hpp
 * @brief contains the definition of the GSL_interpolator class
 */

#include <Eigen/Core>

#include <gsl/gsl_spline.h>
#include <gsl/gsl_version.h>

#include <memory>

namespace BubbleProfiler {

/**
 * @brief Computes interpolated function values at a new set of grid points
 * @param[in] x_old the initial set of x values at which the function is known
 * @param[in] y_old the initial set of y values
 * @param[in] x_new the x values at which to interpolate the function
 * @return a vector of the interpolated values at each point in \c x_new
 */
Eigen::VectorXd interpolate_f_at(const Eigen::VectorXd& x_old,
                                 const Eigen::VectorXd& y_old,
                                 const Eigen::VectorXd& x_new);

/**
 * @brief Computes interpolated function values at a new set of grid points
 * @param[in] x_old the initial set of x values at which the function is known
 * @param[in] y_old the initial set of y values
 * @param[in] n_old the number of values to use from \c x_old and \c y_old
 * @param[in] x_new the x values at which to interpolate the function
 * @param[in] n_new the number of values in \c x_new to evaluate at
 * @return a vector of the interpolated values at the values in \c x_new
 */
Eigen::VectorXd interpolate_f_at(const Eigen::VectorXd& x_old,
                                 const Eigen::VectorXd& y_old,
                                 int n_old, const Eigen::VectorXd& x_new,
                                 int n_new);

/**
 * @class GSL_interpolator
 * @brief Provides routines for carrying out 1D interpolation
 */
class GSL_interpolator {
public:
   enum class Interpolation_type {
      GSL_linear,
      GSL_polynomial,
      GSL_cspline,
      GSL_cspline_periodic,
      GSL_akima,
      GSL_akima_periodic
#ifdef GSL_MAJOR_VERSION
#if (GSL_MAJOR_VERSION >= 2)
      , GSL_steffen
#endif
#endif
   };

   GSL_interpolator();
   ~GSL_interpolator() = default;
   GSL_interpolator(const GSL_interpolator&);
   GSL_interpolator& operator=(const GSL_interpolator&);
   GSL_interpolator(GSL_interpolator&&) = default;
   GSL_interpolator& operator=(GSL_interpolator&&) = default;

   /**
    * @brief Sets the method used to perform interpolation
    * @param t the interpolation method to use
    */
   void set_interpolation_type(Interpolation_type t) { interp_type = t; }

   /**
    * @brief Computes interpolant using given x and y values
    * @param[in] x the x values to use
    * @param[in] y the function values to use
    */
   void construct_interpolant(const Eigen::VectorXd&, const Eigen::VectorXd&);

   /**
    * @brief Computes interpolant using requested number of points
    * @param[in] x the x values to use
    * @param[in] y the function values to use
    * @param[in] n the number of points to use to build the spline
    */
   void construct_interpolant(const Eigen::VectorXd&, const Eigen::VectorXd&,
                              int);

   /**
    * @brief Evaluates computed interpolant at a point
    * @param[in] x the point to evaluate at
    * @return the interpolated function value
    */
   double evaluate_at(double) const;

   /**
    * @brief Evaluates derivative of computed interpolant at a point
    * @param[in] x the point to evaluate at
    * @return the interpolated function derivative value
    */
   double evaluate_deriv_at(double) const;

   /**
    * @brief Evaluates second derivative of computed interpolant at a point
    * @param[in] x the point to evaluate at
    * @return the interpolated function second derivative value
    */
   double evaluate_second_deriv_at(double) const;

private:
   struct GSL_accel_deleter {
      void operator()(gsl_interp_accel* a) const {
         if (a) {
            gsl_interp_accel_free(a);
         }
      }
   };

   struct GSL_spline_deleter {
      void operator()(gsl_spline* s) const {
         if (s) {
            gsl_spline_free(s);
         }
      }
   };

   Interpolation_type interp_type{Interpolation_type::GSL_cspline};
   std::unique_ptr<gsl_interp_accel, GSL_accel_deleter> accelerator{nullptr};
   std::unique_ptr<gsl_spline, GSL_spline_deleter> interpolant{nullptr};

   const gsl_interp_type* to_gsl_interp_type(Interpolation_type) const;
   void reset(Interpolation_type, int);
   void check_interpolant_computed() const;
   void check_within_bounds(double) const;
};

} // namespace BubbleProfiler

#endif
