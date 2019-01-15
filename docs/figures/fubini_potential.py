"""
Illustration of fubini potential
================================
"""

from matplotlib import pyplot as plt
import matplotlib.ticker as plticker
import numpy as np


RC_PARAMS = {'text.latex.preamble' : [r'\usepackage{amsmath}'],
             'font.family': 'lmodern',
             'font.size': 15,
             'axes.titlesize': 10,
             'xtick.labelsize': 15,
             'ytick.labelsize': 15,
             'axes.labelsize': 15,
             'legend.fontsize': 10,
             'text.usetex': True,
             'lines.linewidth': 2,
             'axes.linewidth': 1,
             'axes.grid': True,
             'grid.linestyle': '--',
             'legend.framealpha': 1.,
             'legend.edgecolor': 'black',
             'savefig.bbox': 'tight'}


# Potential

def potential(m, u=1., v=1.):
    return lambda x: 4 * u * m**2 * (m - 1) / (2 * m + 1) * x**(2 + 1 / m) - 2 * u * v * m**2 * x**(2 + 2 / m)

# RC parameters

plt.rcParams.update(RC_PARAMS)

# Make colors

deltas = [1e-4, 1e-3, 1e-2, 0.025, 0.05, 0.075, 1e-1]

cmap = plt.get_cmap("Paired", len(deltas))

# Plot potentials

phi = np.linspace(0., 0.04, 1000)

for i, delta in enumerate(reversed(deltas)):
    p = potential(1. + delta)
    label = r"$m = 1 + {}$".format(delta)
    plt.plot(phi, -p(phi), label=label, color=cmap(i), lw=3)

# Finish up

# plt.xlim([-0.25, 1.25])
# plt.ylim([-0.05, 0.3])
# loc = plticker.MultipleLocator(base=0.25)
# plt.gca().xaxis.set_major_locator(loc)
plt.grid(True)

handles, labels = plt.gca().get_legend_handles_labels()
leg = plt.legend(handles[::-1], labels[::-1], loc='upper left')
leg.get_frame().set_alpha(0.75)

plt.xlabel(r'Field, $\phi$')
plt.ylabel(r'Upturned potential, $-V(\phi)$')
plt.title(r'$V(\phi) = 4 m^2  (m - 1) / (2  m + 1)  \phi^{(2 + 1 / m)} - 2  m^2  \phi^{(2 + 2 / m)}$')

plt.tight_layout()
plt.savefig('figures/fubini_potential.pdf')
