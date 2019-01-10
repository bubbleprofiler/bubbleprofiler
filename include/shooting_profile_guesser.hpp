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

#ifndef BUBBLEPROFILER_SHOOTING_PROFILE_GUESSER_HPP_INCLUDED
#define BUBBLEPROFILER_SHOOTING_PROFILE_GUESSER_HPP_INCLUDED

#include "basic_logger.hpp"
#include "profile_guesser.hpp"
#include "shooting.hpp"

#include <Eigen/Core>

namespace BubbleProfiler {

/*!
 * @class Shooting_profile_guesser
 * @brief Constructs a solution ansatz using the overshoot/undershoot method
 *
 * Given a multi-field potential and the location of the true vacuum,
 * the initial guess for the bubble profile is found by first reducing
 * the problem to an approximate 1D problem.  This is then solved using
 * the undershoot/overshoot method, and the obtained solution is converted
 * back into the original field basis.
 */
class Shooting_profile_guesser : public Profile_guesser {
public:
   virtual ~Shooting_profile_guesser() = default;

   virtual Field_profiles get_profile_guess(const Potential&,
                                            const Eigen::VectorXd&,
                                            int, double, double,
                                            double, double) override;

   void set_trial_distance(double dist) { trial_dist = dist; }
   void set_solver(const Shooting& s) { solver = s; }

   const Shooting& get_solver() const { return solver; }
   Shooting& get_solver() { return solver; }

private:
   double trial_dist{0.5};
   double alpha{0.5};
   double E{1.};
   double dist_true_vacuum{0.};
   int num_fields{0};
   Shooting solver{};
   // Change of (field space) basis matrix, takes vectors
   // in the original basis to those in the ansatz basis.
   Eigen::MatrixXd cob_matrix{};
   logging::Basic_logger logger{};

   void compute_vacuum_distance(const Potential&, const Eigen::VectorXd&);
   void calculate_potential_parameters(const Potential&,
                                       const Eigen::VectorXd&);
   double calculate_false_min_location() const;
   double calculate_barrier_location() const;
   double calculate_true_min_location() const;
   double get_large_distance_solution(double, const Field_profiles&, int) const;

   void solve_one_dimensional(int);
   Field_profiles calculate_field_profiles(int, double, double, double, double);
   Eigen::VectorXd update_spatial_grid(const Eigen::VectorXd&, double, double, int) const;
};

} // namespace BubbleProfiler

#endif
