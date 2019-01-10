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

#ifndef BUBBLEPROFILER_GENERALIZED_FUBINI_OBSERVER_HPP_INCLUDED
#define BUBBLEPROFILER_GENERALIZED_FUBINI_OBSERVER_HPP_INCLUDED

#include "generalized_fubini_potential.hpp"

#include <fstream>
#include <string>

namespace BubbleProfiler {

class Field_profiles;

class Generalized_fubini_observer {
public:
   Generalized_fubini_observer(const Generalized_fubini_potential& potential_,
                               const std::string& output_file_);

   void operator()(const Field_profiles& profile,
                   const Field_profiles& perturbation);
private:
   int iteration_count{0};
   Generalized_fubini_potential potential{};
   std::string output_file{};
   std::ofstream output_stream;

   void write_file_header();
   void write_data(const Field_profiles&, const Field_profiles&,
                   const Field_profiles&);
};

} // namespace BubbleProfiler

#endif
