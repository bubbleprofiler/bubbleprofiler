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

#ifndef BUBBLEPROFILER_RESTRICTED_QUARTIC_POTENTIAL_HPP_INCLUDED
#define BUBBLEPROFILER_RESTRICTED_QUARTIC_POTENTIAL_HPP_INCLUDED

#include "potential.hpp"

#include <Eigen/Core>

namespace BubbleProfiler {

class Restricted_quartic_potential : public Potential {
public:
   Restricted_quartic_potential() = default;
   explicit Restricted_quartic_potential(double);
   Restricted_quartic_potential(double, double);
   virtual ~Restricted_quartic_potential() = default;
   Restricted_quartic_potential(const Restricted_quartic_potential&) = default;
   Restricted_quartic_potential(Restricted_quartic_potential&&) = default;
   Restricted_quartic_potential& operator=(const Restricted_quartic_potential&) = default;
   Restricted_quartic_potential& operator=(Restricted_quartic_potential&&) = default;

   virtual Restricted_quartic_potential * clone() const override {
      return new Restricted_quartic_potential(*this);
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

   double get_global_minimum_location() const;
   double get_local_minimum_location() const;
   double get_local_maximum_location() const;

private:
   double alpha{0.6};
   double E{1.};
   double origin{0.};
   double scale{1.};
   double offset{0.};
};

} // namespace BubbleProfiler

#endif
