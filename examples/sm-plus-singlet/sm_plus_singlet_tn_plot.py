from __future__ import print_function

import sys
import argparse
import numpy as np
import sm_plus_singlet

def get_plot_data(min_temp, max_temp, num_points, Tc, lambda_m, lambda_s, retry):
    # List of tuples (T, action)
    data = []

    for T in np.linspace(min_temp, max_temp, num_points):
        tries = 0

        #HORRID HACK
        epsilon = 0.001
        Trun = T

        while tries < retry:
            try:
                potential_data = sm_plus_singlet.generate_potential(Trun, Tc, lambda_m, lambda_s)
                profiler_output = sm_plus_singlet.run_profiler(potential_data)
                break
            except:
                Trun = Trun + epsilon
                print("retrying at T=" + str(Trun))
                tries += 1

        if tries == retry:
            print("Error: profiler failed to converge for T = " + str(T) +
                  " after " + str(retry) + " attempts")
            sys.exit(1)

        data.append((Trun, profiler_output["action"]))

        #TEMP
        print(Trun, profiler_output["action"],
              "S_E/T=", str(profiler_output["action"] / Trun))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate a plot of S_E/T for an instance of the SM+Singlet model')

    parser.add_argument("--min-temp", type=float, default=0.0,
                        help="lower bound on temp for plot")
    parser.add_argument("--max-temp", type=float, required=True,
                        help="upper bound on temp for plot")
    parser.add_argument("--num-points", type=int, required=True,
                        help="number of points on plot")
    parser.add_argument('-Tc', '--crit-temp', type=float, required=True,
                        help='critical temperature for this potential')
    parser.add_argument('--lambda-m', type=float, required=True,
                        help='mixed quartic coupling parameter')
    parser.add_argument('--lambda-s', type=float, required=True,
                        help='singlet quartic coupling parameter')
    parser.add_argument("--retry", type=int, default=3,
                        help="retry a given point this many times before failing")

    args = parser.parse_args()

    plot_data = get_plot_data(
        args.min_temp, args.max_temp, args.num_points, args.crit_temp,
        args.lambda_m, args.lambda_s, args.retry
    )
