"""
Illustration of logarithmic potential
=====================================
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

def potential(omega, mass=1.):
    return lambda x: 0.5 * mass**2 * x**2 * (1. - np.log(x**2 / omega**2))

# RC parameters

plt.rcParams.update(RC_PARAMS)

# Make colors

omegas = np.linspace(1, 10, 10)

cmap = plt.get_cmap("Paired", len(omegas))

# Plot potentials

phi = np.linspace(0., 10., 1000)

for i, omega in enumerate(reversed(omegas)):
    p = potential(omega)
    label = r"$\omega = {}$".format(int(omega))
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
plt.title(r'$V(\phi) = \frac{1}{2} \phi^2 \left(1 - \log\left(\frac{\phi^2}{\omega^2}\right)\right)$')

plt.tight_layout()
plt.savefig('figures/logarithmic_potential.pdf')
