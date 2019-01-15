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

def potential(epsilon, lambda_=1., a=1.):
    return lambda x: lambda_ / 8 * (x**2 - a**2)**2 + 0.5 * epsilon / a * (x - a)

# RC parameters

plt.rcParams.update(RC_PARAMS)

# Make colors

epsilons = np.linspace(0.01, 0.1, 10)

cmap = plt.get_cmap("Paired", len(epsilons))

# Plot potentials

phi = np.linspace(-1.2, 1.2, 1000)

for i, epsilon in enumerate(epsilons):
    p = potential(epsilon)
    label = r"$\epsilon = {}$".format(epsilon)
    plt.plot(phi, -p(phi), label=label, color=cmap(i), lw=3)

# Finish up

# plt.xlim([-0.25, 1.25])
# plt.ylim([-0.05, 0.3])
# loc = plticker.MultipleLocator(base=0.25)
# plt.gca().xaxis.set_major_locator(loc)
plt.grid(True)

handles, labels = plt.gca().get_legend_handles_labels()
leg = plt.legend(handles[::-1], labels[::-1], loc='upper center')
leg.get_frame().set_alpha(0.75)

plt.xlabel(r'Field, $\phi$')
plt.ylabel(r'Upturned potential, $-V(\phi)$')
plt.title(r'$V(\phi) = \frac{1}{8} (\phi^2 - 1)^2 + \frac{1}{2} \epsilon (\phi - 1)$')

plt.tight_layout()
plt.savefig('figures/thin_wall_potential.pdf')
