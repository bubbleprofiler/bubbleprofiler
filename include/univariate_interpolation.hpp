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

#ifndef BUBBLEPROFILER_UNIVARIATE_INTERPOLATION_HPP_INCLUDED
#define BUBBLEPROFILER_UNIVARIATE_INTERPOLATION_HPP_INCLUDED

#include <algorithm>

namespace BubbleProfiler {

/**
 * Assumes data sorted in order of ascending x-values
 */
template <typename Scalar, typename Data>
Scalar linear_lagrange_interpolation_at_point(Scalar x, const Data& xdata,
                                              const Data& ydata)
{
   using index_type = typename Data::size_type;

   const auto num_x_points = xdata.size();

   const Scalar xmin = xdata[0];
   const Scalar xmax = xdata[num_x_points - 1];

   index_type upper_index;
   if (x <= xmin) {
      upper_index = 1;
   } else if (x >= xmax) {
      upper_index = num_x_points - 1;
   } else {
      const auto upper_it = std::upper_bound(std::begin(xdata),
                                             std::end(xdata), x);
      upper_index = std::distance(std::begin(xdata), upper_it);
   }

   index_type lower_index = upper_index - 1;

   const Scalar x_lower = (x - xdata[upper_index])
      / (xdata[lower_index] - xdata[upper_index]);
   const Scalar x_upper = (x - xdata[lower_index])
      / (xdata[upper_index] - xdata[lower_index]);

   return ydata[lower_index] * x_lower + ydata[upper_index] * x_upper;
}

template <typename Scalar, typename Data>
Scalar quadratic_lagrange_interpolation_at_point(Scalar x, const Data& xdata,
                                                 const Data& ydata)
{
   using index_type = typename Data::size_type;

   const auto num_x_points = xdata.size();

   index_type upper_index;
   if (x <= xdata[1]) {
      upper_index = 2;
   } else if (x >= xdata[num_x_points - 1]) {
      upper_index = num_x_points - 1;
   } else {
      const auto upper_it = std::upper_bound(std::begin(xdata),
                                             std::end(xdata), x);
      upper_index = std::distance(std::begin(xdata), upper_it);
   }

   index_type center_index = upper_index - 1;
   index_type lower_index = upper_index - 2;

   const Scalar xdata_lower = xdata[lower_index];
   const Scalar xdata_center = xdata[center_index];
   const Scalar xdata_upper = xdata[upper_index];

   const Scalar x_lower = (x - xdata_center) * (x - xdata_upper)
      / ((xdata_lower - xdata_center) * (xdata_lower - xdata_upper));
   const Scalar x_center = (x - xdata_lower) * (x - xdata_upper)
      / ((xdata_center - xdata_lower) * (xdata_center - xdata_upper));
   const Scalar x_upper = (x - xdata_lower) * (x - xdata_center)
      / ((xdata_upper - xdata_lower) * (xdata_upper - xdata_center));

   return ydata[lower_index] * x_lower + ydata[center_index] * x_center
      + ydata[upper_index] * x_upper;
}

} // namespace BubbleProfiler

#endif
