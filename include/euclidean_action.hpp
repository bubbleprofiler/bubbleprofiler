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

#ifndef BUBBLEPROFILER_EUCLIDEAN_ACTION_HPP_INCLUDED
#define BUBBLEPROFILER_EUCLIDEAN_ACTION_HPP_INCLUDED

/**
 * @file euclidean_action.hpp
 * @brief contains helper functions for calculating the Euclidean action
 */

#include "integration_utils.hpp"

namespace BubbleProfiler {

class Potential;
class Field_profiles;

/**
 * @brief Calculates the action using only the kinetic term contributions
 *
 * @param potential the potential associated with the bounce solution
 * @param profiles the field profiles representing the bounce solution
 * @param max_intervals the maximum allowed number of subintervals for
 *                      adaptive quadrature routine
 * @param rel_tol the relative error goal
 * @param abs_tol the absolute error goal
 * @param rule the quadrature rule to use
 *
 * @return numerical estimate for the Euclidean action
 */
double calculate_kinetic_action(const Potential& potential,
                                const Field_profiles& profiles,
                                std::size_t max_intervals = 1000,
                                double rel_tol = 1.e-4,
                                double abs_tol = 1.e-4,
                                Integration_rule rule = Integration_rule::GK31);

/**
 * @brief Calculates the action using only the potential term contributions
 *
 * @param potential the potential associated with the bounce solution
 * @param profiles the field profiles representing the bounce solution
 * @param max_intervals the maximum allowed number of subintervals for
 *                      adaptive quadrature routine
 * @param rel_tol the relative error goal
 * @param abs_tol the absolute error goal
 * @param rule the quadrature rule to use
 *
 * @return numerical estimate for the Euclidean action
 */
double calculate_potential_action(const Potential& potential,
                                  const Field_profiles& profiles,
                                  std::size_t max_intervals = 1000,
                                  double rel_tol = 1.e-4,
                                  double abs_tol = 1.e-4,
                                  Integration_rule rule = Integration_rule::GK31);


/**
 * @brief Calculates the action using all terms in the action integrand
 *
 * @param potential the potential associated with the bounce solution
 * @param profiles the field profiles representing the bounce solution
 * @param max_intervals the maximum allowed number of subintervals for
 *                      adaptive quadrature routine
 * @param rel_tol the relative error goal
 * @param abs_tol the absolute error goal
 * @param rule the quadrature rule to use
 *
 * @return numerical estimate for the Euclidean action
 */
double calculate_full_action(const Potential& potential,
                             const Field_profiles& profiles,
                             std::size_t max_intervals = 1000,
                             double rel_tol = 1.e-4,
                             double abs_tol = 1.e-4,
                             Integration_rule rule = Integration_rule::GK31);

/**
 * @brief Calculates the action using either the kinetic or potential terms
 *
 * @param potential the potential associated with the bounce solution
 * @param profiles the field profiles representing the bounce solution
 * @param max_intervals the maximum allowed number of subintervals for
 *                      adaptive quadrature routine
 * @param rel_tol the relative error goal
 * @param abs_tol the absolute error goal
 * @param rule the quadrature rule to use
 * @param use_kinetic flag selecting whether the kinetic or potential terms
 *                    should be used to compute the action
 *
 * @return numerical estimate for the Euclidean action
 */
double calculate_action(const Potential& potential,
                        const Field_profiles& profiles,
                        std::size_t max_intervals = 1000,
                        double rel_tol = 1.e-4,
                        double abs_tol = 1.e-4,
                        Integration_rule rule = Integration_rule::GK31,
                        bool use_kinetic = true);

} // namespace BubbleProfiler

#endif
