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

#ifndef BUBBLEPROFILER_LOGARITHMIC_POTENTIAL_HPP_INCLUDED
#define BUBBLEPROFILER_LOGARITHMIC_POTENTIAL_HPP_INCLUDED

#include "field_profiles.hpp"
#include "potential.hpp"

#include <Eigen/Core>

namespace BubbleProfiler {

class Logarithmic_potential : public Potential {
public:
   Logarithmic_potential() = default;
   Logarithmic_potential(double m_, double w_)
      : m(m_), w(w_)
      {}
   virtual ~Logarithmic_potential() = default;
   Logarithmic_potential(const Logarithmic_potential&) = default;
   Logarithmic_potential(Logarithmic_potential&&) = default;
   Logarithmic_potential& operator=(const Logarithmic_potential&) = default;
   Logarithmic_potential& operator=(Logarithmic_potential&&) = default;

   virtual Logarithmic_potential* clone() const override {
      return new Logarithmic_potential(*this);
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
   double m{1.};
   double w{1.};
   double origin{0.};
   double scale{1.};
   double offset{0.};
};

} // namespace BubbleProfiler

#endif
