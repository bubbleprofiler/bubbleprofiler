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

#ifndef BUBBLEPROFILER_NUMERIC_HPP_INCLUDED
#define BUBBLEPROFILER_NUMERIC_HPP_INCLUDED

/**
   @file
   @brief Numerical functions for one-dimensional shooting method and
   area of \f$n\f$-sphere
*/

namespace BubbleProfiler {

/**
   @returns Area of \f$n\f$-sphere from <a href="https://en.wikipedia.org/wiki/N-sphere#Recurrences">wiki</a> 
   @param n \f$n\f$
*/
double area_n_sphere(int n);
/**
   @returns Approximate inverse of the hyperbolic sinch function, \f$\sinh x / x = a\f$

   I use Taylor expansions and \f$\sinh x \approx e^x / 2\f$.

   See e.g., in Mathematica,
   @code
   Asinch[a_] := x /. FindRoot[Sinh[x] /x == a, {x, 5}]
   Taylor[a_] := Sqrt[6. (a - 1)]
   Lambert[a_] := -LambertW[-1, -0.5 / a]
   Plot[{Asinch[a], Taylor[a], Lambert[a]}, {a, 1, 10}]
   @endcode
   The Taylor expansion is more precise until about \f$a \approx 1.5\f$.
*/
double asinch(double a);
/**
   @returns Approximate inverse of the sinc function, \f$\sin x / x = a \f$

   See e.g., in Mathematica
   @code
   Asinch[a_] := InverseFunction[Sinc][a] // N // Re
   f[x_] := 2*x + 3*x^3/10 + 321*x^5/2800 + 3197*x^7/56000 +
   445617*x^9/13798400 + 1766784699*x^11/89689600000 +
   317184685563*x^13/25113088000000 +
   14328608561991*x^15/1707689984000000 +
   6670995251837391*x^17/1165411287040000000 +
   910588298588385889*x^19/228420612259840000000
   Approx[a_] := Sqrt[3./2.] * f[Sqrt[1. - a]]
   Plot[{Asinch[a], Approx[a]}, {a, 0, 1}]
   Plot[{Sinc[Approx[a]]}, {a, 0, 1}]
   @endcode
*/
double asinc(double a);
/**
   @returns Root of \f$I_1(x) / x - a - 0.5 = 0\f$

   This is used in dimension \f$d=4\f$ case. I use Taylor expansions and a
   large \f$x\f$ approximation of \f$I_1(x)\f$.

   See e.g., in Mathematica,
   @code
   Numerical[a_] := x /. FindRoot[BesselI[1, x] / x - 0.5 == a, {x, 10.}]
   Taylor[a_] := 4. Sqrt[a]
   Lambert[a_] := -(3/2) LambertW[-1, (2 (2/\[Pi])^(1/3))/(3 (-(1 + 2 a)^2)^(1/3))] // N // Re
   Plot[{Numerical[a], Taylor[a], Lambert[a]}, {a, 0, 100}]
   @endcode
*/
double approx_root_pos_4(double a);
/**
   @returns  Root of \f$ J_1(x) / x + a - 0.5 = 0 \f$

   See e.g., in Mathematica,
   @code
   Numerical[a_] := x /. FindRoot[BesselJ[1, x]/x - 0.5 == -a, {x, 2}]
   Taylor[a_] := 4. Sqrt[a]
   Plot[{Numerical[a], Taylor[a]}, {a, 0, 1}]
   @endcode
*/
double approx_root_neg_4(double a);

}  // namespace BubbleProfiler

#endif  // BUBBLEPROFILER_NUMERIC_HPP_INCLUDED
