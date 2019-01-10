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

#include "algebraic_potential.hpp"
#include "error.hpp"

#include <sstream>

namespace BubbleProfiler {

// Potential constructor - from expression string
Algebraic_potential::Algebraic_potential(
   const std::vector<std::string>& fields_,
   const std::string& expr)
   : fields(fields_)
{
   // Parse the field expression & get symbol table
   GiNaC::parser reader;
   v = reader(expr);
   GiNaC::symtab table = reader.get_syms();

   // Match fields names to symbol table
   // TODO this should have error handling
   for(auto& field : fields_) {
      syms.push_back(GiNaC::ex_to<GiNaC::symbol>(table.at(field)));
   }

   build_derivatives();

   std::stringstream log_str;
   log_str << "Input potential: " << v;
   logger.log_message(logging::Log_level::Trace, log_str.str());
}

void Algebraic_potential::build_derivatives()
{
   first_partials.clear();
   second_partials.clear();

   for (const auto& sym: syms) {
      GiNaC::ex deriv = v.diff(sym);
      first_partials.push_back(deriv);

      std::vector<GiNaC::ex> row_partials;

      for (const auto& col_sym: syms) {
         GiNaC::ex deriv2 = v.diff(sym).diff(col_sym);
         row_partials.push_back(deriv2);
      }

      second_partials.push_back(row_partials);
   }
}

void Algebraic_potential::translate_origin(
   const Eigen::VectorXd& translation)
{
   const std::size_t n_coords = translation.size();
   if (n_coords != syms.size()) {
      throw Setup_error(
         "Algebraic_potential::translate_origin: "
         "dimensions of translation do not match number of fields");
   }

   // Create coordinate substitutions
   GiNaC::lst l;
   for (std::size_t i = 0; i < n_coords; ++i) {
      l.append(syms[i] == syms[i] + translation[i]);
   }

   // Apply them
   v = v.subs(l);

   // Rebuild algebraic derivatives
   build_derivatives();

   std::stringstream log_str;
   log_str << "Potential after translation: " << v;
   logger.log_message(logging::Log_level::Trace, log_str.str());
}

void Algebraic_potential::apply_basis_change(const Eigen::MatrixXd &cob_matrix) {
   Eigen::MatrixXd cob_matrix_t = cob_matrix.transpose();
   int n_coords = syms.size();

   // Create basis transformation substitution map
   GiNaC::exmap substitutions;

   for (int old_field_ix = 0; old_field_ix < n_coords; ++old_field_ix) {

      // Build the RHS of the substitution for this field.
      GiNaC::ex rhs = 0;

      for (int new_field_ix = 0; new_field_ix < n_coords; ++new_field_ix) {
         double coeff = cob_matrix(new_field_ix, old_field_ix);
         rhs += coeff * syms[new_field_ix];
      }

      // Add the substitution for this field to the transformation map
      substitutions[syms[old_field_ix]] = rhs;
   }

   std::stringstream subs_str;
   subs_str << "Subs: " << substitutions;
   logger.log_message(logging::Log_level::Trace, subs_str.str());

   // Apply all of the transformations simultaneously
   v = v.subs(substitutions);

   // Rebuild algebraic derivatives
   build_derivatives();

   std::stringstream expr_str;
   expr_str << "Potential after translation: " << v;
   logger.log_message(logging::Log_level::Trace, expr_str.str());
}

void Algebraic_potential::add_constant_term(double offset)
{
   v += offset;

   std::stringstream expr_str;
   expr_str << "Potential after adding constant term: " << v;
   logger.log_message(logging::Log_level::Trace, expr_str.str());
}

// Evaluate a GiNaC expression (i.e. potential, partial derivatives) on
// this potential.
double Algebraic_potential::eval(
   const GiNaC::ex& expr, const Eigen::VectorXd& coords) const
{
   const std::size_t n_coords = coords.size();
   if (n_coords != syms.size()) {
      throw Setup_error(
         "Algebraic_potential::eval: "
         "number of values does not match number of coordinates");
   }

   // A list of substitutions
   GiNaC::lst l;
   for (std::size_t i = 0; i < n_coords; ++i) {
      l.append(syms[i] == coords(i));
   }

   // Make subs & eval
   GiNaC::ex result = GiNaC::evalf(expr.subs(l));

   return GiNaC::ex_to<GiNaC::numeric>(result).to_double();
}

// Evaluate potential at a coordinate
double Algebraic_potential::operator()(const Eigen::VectorXd& coords) const
{
   return Algebraic_potential::eval(v, coords);
}

// Partial derivative WRT coordinate i at a point
double Algebraic_potential::partial(const Eigen::VectorXd& coords, int i) const
{
   return Algebraic_potential::eval(first_partials[i], coords);
}

// Partial derivatives WRT coordinates i, j at a a point
double Algebraic_potential::partial(const Eigen::VectorXd& coords,
                                     int i, int j) const
{
   return Algebraic_potential::eval(second_partials[i][j], coords);
}

} // namespace BubbleProfiler
