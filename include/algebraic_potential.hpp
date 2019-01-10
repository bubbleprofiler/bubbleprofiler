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

#ifndef BUBBLEPROFILER_ALGEBRAIC_POTENTIAL_HPP_INCLUDED
#define BUBBLEPROFILER_ALGEBRAIC_POTENTIAL_HPP_INCLUDED

#include "basic_logger.hpp"
#include "potential.hpp"

#include <ginac/ginac.h>

namespace BubbleProfiler {

//! Implements a multi dimensional scalar field potential, which may be expressed as almost any algerbraic fucntion.
/*!
 * Given a vector of field names, and a GiNaCs expression specifying
 * a scalar function of these fields, this class provides a an implementation
 * of Potential suitable for calculating bubble profiles.
 *
 * For the profiler to work correctly, the potential should be a fourth degree polynomial.
 */
class Algebraic_potential : public Potential {
public:
   //! Construct from list of fields and GiNaC potential string.
   /*!
    * @param fields vector of field names. These must match the fields in expr.
    * @param expr specification of potential in GiNaCs syntax.
    *
    * Note that the potential is expected to be a quartic multivariate
    * polynomial in the specified fields.
    */
   Algebraic_potential(const std::vector<std::string>& fields,
                        const std::string& expr);

   virtual ~Algebraic_potential() = default;

   //! Get a fresh copy of this potential
   Algebraic_potential * clone() const {
      return new Algebraic_potential(*this);
   }

   virtual double operator()(const Eigen::VectorXd& coords) const override;
   virtual double partial(const Eigen::VectorXd& coords, int i) const override;
   virtual double partial(const Eigen::VectorXd& coords, int i, int j) const override;
   virtual std::size_t get_number_of_fields() const override {
      return fields.size(); }

   virtual void translate_origin(const Eigen::VectorXd&) override;
   virtual void apply_basis_change(const Eigen::MatrixXd&) override;
   virtual void add_constant_term(double) override;

   // Printable
   inline friend std::ostream& operator<<(
      std::ostream& os, const Algebraic_potential& p) {
      os << p.v;
      return os;
   }

private:
   std::vector<std::string> fields{};
   GiNaC::ex v{}; ///< GiNaC expression to hold potential
   std::vector<GiNaC::ex> first_partials{}; ///< First partials
   std::vector<std::vector<GiNaC::ex> > second_partials{}; ///< Second partials
   std::vector<GiNaC::symbol> syms{}; ///< GiNaC symbols matched to fields vec
   logging::Basic_logger logger{};

   // Evaluate a GiNaCs expression
   double eval(const GiNaC::ex&, const Eigen::VectorXd&) const;

   // Rebuild algebraic derivatives
   void build_derivatives();
};

} // namespace BubbleProfiler

#endif
