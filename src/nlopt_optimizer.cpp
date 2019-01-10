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
    @file
    @brief Wrapper for NLOpts optimizer
*/

#include "nlopt_optimizer.hpp"
#include "error.hpp"

#include <string>

namespace BubbleProfiler {

bool optimization_succeeded(nlopt::result r)
{
   switch (r) {
   case nlopt::FAILURE: return false;
   case nlopt::INVALID_ARGS: return false;
   case nlopt::OUT_OF_MEMORY: return false;
   case nlopt::ROUNDOFF_LIMITED: return false;
   case nlopt::FORCED_STOP: return false;
   case nlopt::SUCCESS: return true;
   case nlopt::STOPVAL_REACHED: return true;
   case nlopt::FTOL_REACHED: return true;
   case nlopt::XTOL_REACHED: return true;
   case nlopt::MAXEVAL_REACHED: return true;
   case nlopt::MAXTIME_REACHED: return true;
   default:
      throw Optimizer_error("unrecognized NLopt status code: "
                            + std::to_string(r));
   }
}

NLopt_optimizer::NLopt_optimizer(Index n)
   : n_dims(n)
   , lower_bounds(Eigen::VectorXd::Zero(n))
   , upper_bounds(Eigen::VectorXd::Zero(n))
{
}

void NLopt_optimizer::set_xtol_rel(double xtol_rel_)
{
   xtol_rel = xtol_rel_;
}

void NLopt_optimizer::set_ftol_rel(double ftol_rel_)
{
   ftol_rel = ftol_rel_;
}

void NLopt_optimizer::set_lower_bounds(double lb)
{
   lower_bounds = lb * Eigen::VectorXd::Ones(n_dims);
}

void NLopt_optimizer::set_lower_bounds(Index i, double lb)
{
   if (i < 0 || i >= n_dims) {
      throw Out_of_bounds_error(i, "NLopt_optimizer::set_lower_bounds: "
                                "array index out of bounds");
   }
   lower_bounds(i) = lb;
}

void NLopt_optimizer::set_lower_bounds(const Eigen::VectorXd& lb)
{
   if (lb.size() != lower_bounds.size()) {
      throw Setup_error("NLopt_optimizer::set_lower_bounds: "
                        "size of bounds vectors must match");
   }
   lower_bounds = lb;
}

void NLopt_optimizer::set_upper_bounds(double ub)
{
   upper_bounds = ub * Eigen::VectorXd::Ones(n_dims);
}

void NLopt_optimizer::set_upper_bounds(Index i, double ub)
{
   if (i < 0 || i >= n_dims) {
      throw Out_of_bounds_error(i, "NLopt_optimizer::set_upper_bounds: "
                                "array index out of bounds");
   }
   upper_bounds(i) = ub;
}

void NLopt_optimizer::set_upper_bounds(const Eigen::VectorXd& ub)
{
   if (ub.size() != upper_bounds.size()) {
      throw Setup_error("NLopt_optimizer::set_upper_bounds: "
                        "size of bounds vectors must match");
   }
   upper_bounds = ub;
}

void NLopt_optimizer::set_max_time(double maxtime_) {
   if (maxtime < 0) {
      throw Setup_error("NLopt_optimizer::set_max_time: "
                        "maximum allowed time must be non-negative");
   }
   maxtime = maxtime_;
}


nlopt::result NLopt_optimizer::optimize(const Eigen::VectorXd& guess)
{
   if (guess.size() != n_dims) {
      throw Setup_error("NLopt_optimizer::optimize: "
                        "initial guess size must match number of variables");
   }

   nlopt::opt optimizer(algorithm, n_dims);

   void* parameters = &function;
   if (extremum_type == Extremum_type::MIN) {
      optimizer.set_min_objective(nlopt_function, parameters);
   } else {
      optimizer.set_max_objective(nlopt_function, parameters);
   }

   std::vector<double> lb(lower_bounds.data(), lower_bounds.data() + n_dims);
   std::vector<double> ub(upper_bounds.data(), upper_bounds.data() + n_dims);
   optimizer.set_lower_bounds(lb);
   optimizer.set_upper_bounds(ub);

   optimizer.set_ftol_rel(ftol_rel);
   optimizer.set_xtol_rel(xtol_rel);

   // Set max timeout if appropriate
   if (maxtime > 0) {
      optimizer.set_maxtime(maxtime);
   }

   std::vector<double> x(guess.data(), guess.data() + n_dims);
   const nlopt::result status = optimizer.optimize(x, extremum_value);

   // Throw an exception if we hit the max timeout
   if(status == nlopt::MAXTIME_REACHED) {
      throw Optimizer_error("NLopt_optimizer::optimize: "
                            "Timed out before locating extremum");
   }

   extremum = Eigen::VectorXd::Map(x.data(), n_dims);

   return status;
}

double NLopt_optimizer::nlopt_function(const std::vector<double>& x,
                                       std::vector<double>& /* grad */,
                                       void* params)
{
   Function* func = static_cast<Function*>(params);
   const Eigen::VectorXd coords(Eigen::VectorXd::Map(x.data(), x.size()));
   return (*func)(coords);
}

} // namespace BubbleProfiler
