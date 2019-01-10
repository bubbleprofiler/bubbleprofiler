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

#ifndef BUBBLEPROFILER_NLOPT_OPTIMIZER_HPP_INCLUDED
#define BUBBLEPROFILER_NLOPT_OPTIMIZER_HPP_INCLUDED

#include <nlopt.hpp>

#include <Eigen/Core>

namespace BubbleProfiler {

/*!
 * Helper function to check NLopt status code
 */
bool optimization_succeeded(nlopt::result);

/*!
 * Lightweight wrapper class for the NLOpt numerical optimization library.
 */
class NLopt_optimizer {
public:
   using Index = Eigen::VectorXd::Index;
   using Function = std::function<double(const Eigen::VectorXd&)>;

   enum class Extremum_type { MIN, MAX};

   explicit NLopt_optimizer(Index);
   template <typename F>
   NLopt_optimizer(F&&, Index);

   void set_algorithm(nlopt::algorithm a) { algorithm = a; }
   void set_extremum_type(Extremum_type e) { extremum_type = e; }
   template <typename F>
   void set_function(F&& f) { function = std::forward<F>(f); }

   void set_lower_bounds(double lb);
   void set_lower_bounds(Index i, double lb);
   void set_lower_bounds(const Eigen::VectorXd& lb);
   void set_upper_bounds(double ub);
   void set_upper_bounds(Index i, double ub);
   void set_upper_bounds(const Eigen::VectorXd& ub);
   void set_xtol_rel(double xtol_rel_);
   void set_ftol_rel(double ftol_rel_);

   /*!
    * Set a timeout for the optimizer. Failure to converge before
    * timeout is considered an error, and an exception will be thrown.
    * If this parameter is set to zero (the default value), no timeout
    * is applied.
    * @param maxtime maximum runtime (seconds)
    */
   void set_max_time(double maxtime);

   Eigen::VectorXd get_extremum_location() const { return extremum; }
   double get_extremum_value() const { return extremum_value; }

   /*!
    * Run the optimizer with configured parameters.
    * @param guess Vector in which result will be stored
    * @return NLOpt status code
    */
   nlopt::result optimize(const Eigen::VectorXd& guess);

private:
   Index n_dims{};
   nlopt::algorithm algorithm{nlopt::GN_DIRECT_L};
   Extremum_type extremum_type{Extremum_type::MIN};
   Eigen::VectorXd lower_bounds{};
   Eigen::VectorXd upper_bounds{};
   double maxtime{0.};
   Eigen::VectorXd extremum{};
   double extremum_value{0.};
   Function function{nullptr};
   double xtol_rel{2. * std::numeric_limits<double>::epsilon()};
   double ftol_rel{2. * std::numeric_limits<double>::epsilon()};

   static double nlopt_function(const std::vector<double>&,
                                std::vector<double>&, void*);
};

template <typename F>
NLopt_optimizer::NLopt_optimizer(F&& f, Index n)
   : n_dims(n)
   , lower_bounds(Eigen::VectorXd::Zero(n))
   , upper_bounds(Eigen::VectorXd::Zero(n))
   , function(std::forward<F>(f))
{
}

} // namespace BubbleProfiler

#endif
