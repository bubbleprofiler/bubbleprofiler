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

#ifndef BUBBLEPROFILER_INTEGRATION_POLICY_HPP_INCLUDED
#define BUBBLEPROFILER_INTEGRATION_POLICY_HPP_INCLUDED

#include "eigen_state_utils.hpp"

#include <boost/numeric/odeint.hpp>

namespace BubbleProfiler {

template <class Parameter,
          class State,
          class Stepper = boost::numeric::odeint::runge_kutta4<
             State,
             Parameter,
             State,
             Parameter,
             typename State_algebra_dispatcher<State>::algebra_type> >
class Fixed_step_size_integrator {
protected:
   ~Fixed_step_size_integrator() = default;

   template <class System, class Observer>
   int integrate_system(System system, State& state,
                        Parameter from, Parameter to, Parameter step_size,
                        Observer observer) const;
};

template <class Parameter, class State, class Stepper>
template <class System, class Observer>
int Fixed_step_size_integrator<Parameter, State, Stepper>::integrate_system(
   System system, State& state, Parameter from, Parameter to,
   Parameter step_size, Observer observer) const
{
   Stepper stepper;
   return boost::numeric::odeint::integrate_adaptive(
      stepper, system, state, from, to, step_size, observer);
}

template <class Parameter, class State, class Stepper>
class Controlled_step_size_integrator {
public:
   void set_absolute_error(Parameter err) { abs_err = err; }
   void set_relative_error(Parameter err) { rel_err = err; }

protected:
   ~Controlled_step_size_integrator() = default;

   template <class System, class Observer>
   int integrate_system(System system, State& state,
                        Parameter from, Parameter to, Parameter step_size,
                        Observer observer) const;

private:
   Parameter abs_err{1.e-6};
   Parameter rel_err{1.e-6};
};

template <class Parameter, class State, class Stepper>
template <class System, class Observer>
int Controlled_step_size_integrator<Parameter, State, Stepper>::integrate_system(
   System system, State& state, Parameter from, Parameter to,
   Parameter step_size, Observer observer) const
{
   auto stepper = boost::numeric::odeint::make_controlled(abs_err, rel_err,
                                                          Stepper());

   return boost::numeric::odeint::integrate_adaptive(
      stepper, system, state, from, to, step_size, observer);
}

} // namespace BubbleProfiler

#endif
