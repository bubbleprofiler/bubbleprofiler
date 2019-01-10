import argparse
import matplotlib.pyplot as plt
import numpy as np

def parse_cmd_line_args():
    parser = argparse.ArgumentParser(
        description="Plot action and timing values")
    parser.add_argument("data_file", help="data file to read")
    parser.add_argument("-o,--output-file", dest="output_file",
                        default="", help="name to save plot with")
    args = parser.parse_args()
    return (args.data_file, args.output_file)

def get_profilers(column_names):
    potnl_cols = ['dim', 'E', 'alpha']
    profiler_cols = [c for c in column_names if c not in potnl_cols]

    if len(profiler_cols) % 3 != 0:
        raise IOError("invalid data column format")

    return [p[:p.find("_action")] for p in profiler_cols[::3]]

def get_profiler_data(profiler, dim, E, data):
    dim_data = data["dim"]
    E_data = data["E"]

    action_data = data[profiler + "_action"]
    timing_data = data[profiler + "_timems"]
    error_data = data[profiler + "_error"]

    profiler_data = {"action": action_data[(dim_data == dim) & (E_data == E)],
                     "timing": timing_data[(dim_data == dim) & (E_data == E)],
                     "error": error_data[(dim_data == dim) & (E_data == E)]
                    }

    return profiler_data

def read_data(data_file):
    data = np.genfromtxt(data_file, names=True)
    cols = [c for c in data.dtype.fields]

    if "dim" not in cols:
        raise IOError("missing dimension column")
    dim = np.unique(data["dim"])[0]

    if "E" not in cols:
        raise IOError("missing E parameter column")
    E = np.unique(data["E"])[0]

    if "alpha" not in cols:
        raise IOError("missing alpha parameter column")
    alpha_data = data["alpha"]

    profilers = get_profilers(cols)
    profilers_data = {}
    for p in profilers:
        profilers_data[p] = get_profiler_data(p, dim, E, data)

    return (dim, E, alpha_data, profilers_data)

def plot_action(axis, alpha_values, profilers_data):
    for p in profilers_data:
        error_data = profilers_data[p]["error"]
        valid_alpha_values = alpha_values[error_data == 0]
        valid_action_data = profilers_data[p]["action"][error_data == 0]
        axis.semilogy(valid_alpha_values, valid_action_data, label=p)

    axis.set_ylabel("Action")
    axis.legend(numpoints=1)
    axis.grid()

def plot_timings(axis, alpha_values, profilers_data):
    for p in profilers_data:
        error_data = profilers_data[p]["error"]
        valid_alpha_values = alpha_values[error_data == 0]
        valid_timings_data = profilers_data[p]["timing"][error_data == 0]
        axis.semilogy(valid_alpha_values, valid_timings_data, label=p)

    axis.set_ylabel("Time/ms")
    axis.legend(numpoints=1)
    axis.grid()

def plot_results(dim, E, alpha_values, profilers_data, output_file):
    fig, ax = plt.subplots(2, 1, sharex=True)

    fig.suptitle("d = {:d}, E = {:.2f}".format(int(dim), E))

    plot_action(ax[0], alpha_values, profilers_data)
    plot_timings(ax[1], alpha_values, profilers_data)

    ax[-1].set_xlabel("alpha")

    if output_file:
        plt.savefig(output_file)
    else:
        plt.show()

    plt.close(fig)

if __name__ == "__main__":
    data_file, output_file = parse_cmd_line_args()
    dim, E, alpha_data, profilers_data = read_data(data_file)
    plot_results(dim, E, alpha_data, profilers_data, output_file)
