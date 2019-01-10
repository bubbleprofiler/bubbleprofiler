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

/**
    @file
    @brief Numerical functions for one-dimensional shooting method and
    area of \f$n\f$-sphere
*/

#include "numeric.hpp"
#include "math_wrappers.hpp"

namespace BubbleProfiler {

double area_n_sphere(int n) {
   double v = 1.;
   double s = 2.;
   for (int i = 1; i <= n; ++i) {
      const double s_n = 2. * Pi * v;
      const double v_n = s / i;
      s = s_n;
      v = v_n;
   }
   return s;
}

double Wm1(double a) {
   /**
     @brief Approximation to negative branch of the Lambert-\f$W\f$ function.

     See Eq. 4.19 of <a href="https://link.springer.com/article/10.1007/BF02124750">this paper</a>

     This is correct up to order \f$(\log \log a / \log a)^5\f$.

     @returns Negative branch of Lambert-\f$W\f$ function
     @param a
   */
   const std::complex<double> ca(a, 0.);
   const std::complex<double> l1 = Log(ca);
   const std::complex<double> l2 = Log(l1);
   const std::complex<double> sum = l1 - l2
      + l2 / l1
      + 0.5 * l2 * (-2. + l2) / std::pow(l1, 2)
      + l2 * (6. - 9. * l2 + 2. * std::pow(l2, 2)) / (6. * std::pow(l1, 3))
      + l2 * (-12. + 36. * l2 - 22. * std::pow(l2, 2) + 3. * std::pow(l2, 2)) / (12. * std::pow(l1, 4));
   return Re(sum);
}

double asinch(double a) {
   if (a <= 1.5) {
      // Taylor expansion
      return Sqrt(6. * (a - 1.));
   } else {
      // Lambert W-function solution
      return -Wm1(-0.5 / a);
   }
}

double series(double a) {
   /**
      @brief Helper function for approximate inverse of \f$\sinc x\f$

      See in Mathematica,
      @code
      Simplify[Normal[InverseSeries[1 - Series[Sin[x]/x, {x, 0, 20}]]]/Sqrt[3/2] /. x -> x^2, x > 0]
      @endcode
      and this  <a href="https://www.dsprelated.com/showthread/comp.dsp/13099-1.php">link</a>
      for further information.
   */
   return 2. * a + 3. * std::pow(a, 3) / 10. + 321. * std::pow(a, 5) / 2800. + 3197. * std::pow(a, 7) / 56000. +
      445617. * std::pow(a, 9) / 13798400. + 1766784699. * std::pow(a, 11) / 89689600000. +
      317184685563. * std::pow(a, 13) / 25113088000000. +
      14328608561991. * std::pow(a, 15) / 1707689984000000. +
      6670995251837391. * std::pow(a, 17) / 1165411287040000000. +
      910588298588385889. * std::pow(a, 19) / 228420612259840000000.;
}

double asinc(double a) {
   return Sqrt(3. / 2.) * series(Sqrt(1. - a));
}

double approx_root_pos_4(double a) {
   if (a < 1.) {
      // Taylor expansion
      return 4. * Sqrt(a);
   } else {
      // Lambert W-function solution
      const double arg = - 2. / 3. * std::pow(2. / Pi, 1. / 3.)
         / std::pow(1. + 2. * a, 2. / 3.);
      return -1.5 * Wm1(arg);
   }
}

double approx_root_neg_4(double a) {
   return 4. * Sqrt(a);
}

}  // namespace BubbleProfiler
