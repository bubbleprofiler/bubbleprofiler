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

#ifndef BUBBLEPROFILER_KINK_PROFILE_GUESSER_HPP_INCLUDED
#define BUBBLEPROFILER_KINK_PROFILE_GUESSER_HPP_INCLUDED

#include "basic_logger.hpp"
#include "profile_guesser.hpp"

#include <Eigen/Core>

#include <array>

namespace BubbleProfiler {

/*!
 * @class Kink_profile_guesser
 * @brief Computes an initial guess for the bubble profile using a kink ansatz
 */
class Kink_profile_guesser : public Profile_guesser {
public:
   virtual ~Kink_profile_guesser() = default;

   /*!
    * @brief Calculate an initial guess for the bubble profile
    *
    * @param potential the potential for which the profiles are to be computed
    * @param true_vacuum the location of the true vacuum in field space
    * @param n_dimensions the number of space-time dimensions
    * @param domain_start the minimum value of the radial coordinate at which
    *                     the profiles are evaluated.  If negative, a value is
    *                     guessed based on the computed ansatz properties.
    * @param domain_end the maximum value of the radial coordinate at which
    *                   the profiles are evaluated.  If negative, a value is
    *                   guessed based on the computed ansatz properties.
    * @param initial_step_size the initial step size to be used in constructing
    *                          the discretized solution
    * @param interpolation_points_fraction fraction of grid points to use
    *                                      for interpolation
    *
    * @return an initial guess for the field profiles
    */
   virtual Field_profiles get_profile_guess(const Potential&,
                                            const Eigen::VectorXd&,
                                            int, double, double,
                                            double, double) override;

   /*!
    * @brief Sets the trial distance used in estimating the potential parameters
    *
    * The parameters defining the quartic potential from which the ansatz is
    * determined are computed by evaluating the potential at a point on the line
    * between the local and global minima.  Initially, the trial point is taken
    * to be halfway along this line.  The distance of the point along this
    * line can be changed using this method.
    *
    * @param dist the distance between the minima to use
    */
   void set_trial_distance(double dist) { trial_dist = dist; }

   /*!
    * @brief Sets the maximum allowed distance if guessing the domain start
    *
    * If the start of the integration domain is to be guessed instead of
    * being explicitly provided, this function sets the maximum distance away
    * from the origin that the returned profiles are allowed to start.  If the
    * guessed starting distance would be larger than the value of \c r, the
    * provided value of \c r is used instead.
    *
    * @param r the maximum allowed distance from the origin to start from
    */
   void set_max_domain_start(double r);

   /*!
    * @brief Sets the minimum required distance if guessing the domain end
    *
    * If the end of the integration domain is to be guessed instead of being
    * explicitly provided, this function sets the minimum distance away from
    * the origin that the returned profiles are required to reach.  If the
    * guessed domain end would otherwise be less than the value of \c r,
    * the provided value of \c r is used instead.
    *
    * @param r the minimum distance that the profiles must extend to
    */
   void set_min_domain_end(double r);

   /*!
    * @brief returns the calculated value of the parameter alpha
    * @return the value of alpha for the fitted one-dimensional potential
    */
   double get_alpha() const { return alpha; }

   /*!
    * @brief returns the calculated value of the parameter |E|
    * @return the value of |E| for the fitted one-dimensional potential
    */
   double get_aE() const { return aE; }

   /*!
    * @brief returns the value of phi0 used in the initial ansatz
    * @return the value of phi0 used in the initial ansatz
    */
   double get_phi0() const { return phi0; }

   /*!
    * @brief returns the value of delta used in the initial ansatz
    * @return the value of delta used in the initial ansatz
    */
   double get_delta() const { return delta; }

   /*!
    * @brief returns the value of lw used in the initial ansatz
    * @return the value of lw used in the initial ansatz
    */
   double get_lw() const { return lw; }

private:
   /*! Values of alpha used for interpolation in 3 spacetime dimensions */
   const static std::array<double,46> alpha_grid_3d;
   /*! Values of alpha used for thin-wall interpolation in 3 spacetime dimensions */
   const static std::array<double, 50> alpha_grid_thin_3d;
   /*! Value of phi0 for each value of alpha in 3 spacetime dimensions */
   const static std::array<double,46> phi0_grid_3d;
   /*! Value of delta for each value of alpha in 3 spacetime dimensions */
   const static std::array<double,46> delta_grid_3d;
   /*! Value of delta for each value of alpha in 3 spacetime dimensions, thin walled case*/
   const static std::array<double, 50> delta_grid_thin_3d;
   /*! Value of Lw for each value of alpha in 3 spacetime dimensions */
   const static std::array<double,46> lw_grid_3d;

   /*! Values of alpha used for interpolation in 4 spacetime dimensions */
   const static std::array<double,44> alpha_grid_4d;
    /*! Values of alpha used for thin-wall interpolation in 4 spacetime dimensions */
    const static std::array<double, 50> alpha_grid_thin_4d;
   /*! Value of phi0 for each value of alpha in 4 spacetime dimensions */
   const static std::array<double,44> phi0_grid_4d;
   /*! Value of delta for each value of alpha in 4 spacetime dimensions */
   const static std::array<double,44> delta_grid_4d;
    /*! Value of delta for each value of alpha in 4 spacetime dimensions, thin walled case*/
    const static std::array<double, 50> delta_grid_thin_4d;
   /*! Value of Lw for each value of alpha in 4 spacetime dimensions */
   const static std::array<double,44> lw_grid_4d;

   /*! Distance of trial point between minima */
   double trial_dist{0.5};
   /*! Maximum allowed starting distance from origin */
   double max_domain_start{1.e-8};
   /*! Minimum required endpoint for guessed profile */
   double min_domain_end{0.};
   /*! Value of alpha in the fitted one-dimensional quartic potential */
   double alpha{0.5};
   /*! Absolute value of E in the fitted one-dimensional quartic potential */
   double aE{1.};
   /*! Computed value of phi0 in the ansatz solution */
   double phi0{0.};
   /*! Computed value of delta in the ansatz solution */
   double delta{0.};
   /*! Computed value of Lw in the ansatz solution */
   double lw{0.};
   /*! Distance to the global minimum from the origin in field space */
   double dist_true_vacuum{0.};
   /*! Number of fields in the potential */
   int num_fields{0};

   /*!
    * Change of (field space) basis matrix, takes vectors
    * in the original basis to those in the ansatz basis.
    */
   Eigen::MatrixXd cob_matrix{};
   logging::Basic_logger logger{};

   void compute_vacuum_distance(const Potential&, const Eigen::VectorXd&);
   void calculate_potential_parameters(const Potential&,
                                       const Eigen::VectorXd&);
   void fit_ansatz_parameters(const Potential&,
                              const Eigen::VectorXd&, int);
   double evaluate_ansatz_at(double) const;
   double evaluate_ansatz_deriv_at(double) const;
   double guess_domain_start() const;
   double guess_domain_end() const;

   //! Minimum value of alpha before throwing a Thin_wall_error
   double alpha_threshold = 0.514;

   Field_profiles calculate_field_profiles(int, double, double, double, double);
};

} // namespace BubbleProfiler

#endif
