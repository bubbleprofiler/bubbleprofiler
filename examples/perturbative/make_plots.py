#!/usr/bin/env python

from __future__ import print_function

import argparse
import matplotlib.pyplot as plt
import numpy as np
import os

def parse_cmd_line_args():
    parser = argparse.ArgumentParser(description="Plot field profiles")
    parser.add_argument("output_dir", help="directory containing output file")
    args = parser.parse_args()
    return args.output_dir

def plot_profiles(profile_file, output_dir, prefix):
    profiles_data = np.genfromtxt(profile_file, names=True)

    fields = [f for f in profiles_data.dtype.names
              if (f != "perturbation" and f != "rho")]
    perturbations = {int(i): profiles_data[profiles_data["perturbation"] == i]
                     for i in np.unique(profiles_data["perturbation"])}

    for f in fields:
        output_file = os.path.join(output_dir, "{}_{}.png".format(prefix, f))
        fig, ax = plt.subplots()

        for p in perturbations:
            ax.plot(perturbations[p]["rho"], perturbations[p][f],
                    label="perturbation {:d}".format(p))

        ax.set_xlabel("rho")
        ax.set_ylabel(f)
        plt.legend(numpoints = 1)

        plt.savefig(output_file)
        plt.close(fig)

def plot_action(action_file, output_dir, prefix):
    action_data = np.genfromtxt(action_file, names=True)
    output_file = os.path.join(output_dir, "{}.png".format(prefix))

    fig, ax = plt.subplots()

    ax.plot(action_data["perturbation"], action_data["action"], "bo")

    ax.set_xlabel("perturbation")
    ax.set_ylabel("action")

    plt.savefig(output_file)
    plt.close(fig)


output_dir = parse_cmd_line_args()
field_profiles_file = os.path.join(output_dir, "field_profiles.txt")
corrections_file = os.path.join(output_dir, "perturbations.txt")
action_file = os.path.join(output_dir, "action.txt")

if os.path.exists(field_profiles_file):
    plot_profiles(field_profiles_file, output_dir, "field_profiles")
else:
    raise IOError("Field profiles file '" + field_profiles_file
                  + "' not found!")

if os.path.exists(corrections_file):
    plot_profiles(corrections_file, output_dir, "corrections")
else:
    raise IOError("Perturbations/corrections file '" + corrections_file
                  + "' not found!")

if os.path.exists(action_file):
    plot_action(action_file, output_dir, "action")
else:
    raise IOError("Action file '" + action_file + "' not found!")
