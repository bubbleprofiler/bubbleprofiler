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

#ifndef BUBBLEPROFILER_GENERALIZED_FUBINI_POTENTIAL_HPP_INCLUDED
#define BUBBLEPROFILER_GENERALIZED_FUBINI_POTENTIAL_HPP_INCLUDED

#include "field_profiles.hpp"
#include "potential.hpp"

#include <Eigen/Core>

namespace BubbleProfiler {

class Generalized_fubini_potential : public Potential {
public:
   Generalized_fubini_potential() = default;
   Generalized_fubini_potential(double, double, double);
   virtual ~Generalized_fubini_potential() = default;
   Generalized_fubini_potential(const Generalized_fubini_potential&) = default;
   Generalized_fubini_potential(Generalized_fubini_potential&&) = default;
   Generalized_fubini_potential& operator=(
      const Generalized_fubini_potential&) = default;
   Generalized_fubini_potential& operator=(
      Generalized_fubini_potential&&) = default;

   virtual Generalized_fubini_potential* clone() const override {
      return new Generalized_fubini_potential(*this);
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

   double operator()(double) const;
   double first_deriv(double) const;
   double second_deriv(double) const;

   double get_local_minimum_location() const;
   double get_local_maximum_location() const;
   double get_bounce_solution_at(double) const;
   Field_profiles get_profile(const Eigen::VectorXd&) const;
   double get_action() const;

private:
   double u{1.};
   double v{1.};
   double m{3.};
   double origin{0.};
   double scale{1.};
   double offset{0.};
};

} // namespace BubbleProfiler

#endif
