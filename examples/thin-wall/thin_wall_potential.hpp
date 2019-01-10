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

#ifndef BUBBLEPROFILER_THIN_WALL_POTENTIAL_HPP_INCLUDED
#define BUBBLEPROFILER_THIN_WALL_POTENTIAL_HPP_INCLUDED

#include "field_profiles.hpp"
#include "potential.hpp"

#include <Eigen/Core>

namespace BubbleProfiler {

/**
 * @class Thin_wall_potential
 * @brief Example thin wall potential as given by Coleman
 *
 * This class implements the example potential
 * \f[
 *    V(\phi) = \frac{\lambda}{8} (\phi^2 - a^2)^2
      + \frac{\epsilon}{2 a} (\phi - a) .
 * \f]
 * and provides analytic expressions for the tunnelling
 * action in the thin-wall approximation.
 */
class Thin_wall_potential : public Potential {
public:
   Thin_wall_potential() = default;
   Thin_wall_potential(double, double, double);
   virtual ~Thin_wall_potential() = default;
   Thin_wall_potential(const Thin_wall_potential&) = default;
   Thin_wall_potential(Thin_wall_potential&&) = default;
   Thin_wall_potential& operator=(
      const Thin_wall_potential&) = default;
   Thin_wall_potential& operator=(
      Thin_wall_potential&&) = default;

   virtual Thin_wall_potential* clone() const override {
      return new Thin_wall_potential(*this);
   }

   virtual double operator()(const Eigen::VectorXd&) const override;
   virtual double partial(const Eigen::VectorXd&, int) const override;
   virtual double partial(const Eigen::VectorXd&, int, int) const override;

   virtual std::size_t get_number_of_fields() const override { return 1; }

   virtual void translate_origin(const Eigen::VectorXd &translation) override {
      origin += translation(0);
   }

   virtual void apply_basis_change(const Eigen::MatrixXd& cob_matrix) override {
      scale *= cob_matrix(0,0);
   }

   virtual void add_constant_term(double _offset) override {
      offset += _offset;
   }

   /**
    * @brief evaluates the potential for the given field value
    * @param coords the value of the field to evaluate the potential at
    * @return the value of the potential
    */
   double operator()(double) const;

   /**
    * @brief evaluates the first derivative of the potential
    * @param coords the value of the field to evaluate the derivative at
    * @return the value of \f$\partial V / \partial \phi\f$
    */
   double first_deriv(double) const;

   /**
    * @brief evaluates the second derivative of the potential
    * @param coords the value of the field to evaluate the derivative at
    * @return the value of \f$\partial^ V / \partial \phi^2\f$
    */
   double second_deriv(double) const;

   /**
    * @brief returns the location of the local minimum
    * @return the location of the local minimum
    */
   double get_local_minimum_location() const;

   /**
    * @brief returns the location of the local maximum
    * @return the location of the local maximum
    */
   double get_local_maximum_location() const;

   /**
    * @brief returns the location of the global minimum
    * @return the location of the global minimum
    */
   double get_global_minimum_location() const;

   /**
    * @brief returns the approximate action in the thin-wall approximation
    * @return value of the Euclidean action in the thin-wall approximation
    */
   double get_thin_wall_action() const;

private:
   double lambda{1.};
   double a{1.};
   double epsilon{0.01};
   double origin{0.};
   double scale{1.};
   double offset{0.};
};

} // namespace BubbleProfiler

#endif
