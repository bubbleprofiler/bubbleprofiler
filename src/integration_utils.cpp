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

#include "integration_utils.hpp"

namespace BubbleProfiler {

int get_integration_rule_key(Integration_rule rule)
{
   switch (rule) {
   case Integration_rule::GK15: return GSL_INTEG_GAUSS15;
   case Integration_rule::GK21: return GSL_INTEG_GAUSS21;
   case Integration_rule::GK31: return GSL_INTEG_GAUSS31;
   case Integration_rule::GK41: return GSL_INTEG_GAUSS41;
   case Integration_rule::GK51: return GSL_INTEG_GAUSS51;
   case Integration_rule::GK61: return GSL_INTEG_GAUSS61;
   default:
      throw Setup_error("unrecognized integration rule");
   }

   return 0;
}

} // namespace BubbleProfiler
