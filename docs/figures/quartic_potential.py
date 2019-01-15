"""
Illustration of quartic potential
=================================
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

def potential(E, alpha):
    return lambda phi: E * ((-4. * alpha + 3.) * 0.5 * phi**2 - phi**3 + alpha * phi**4)

# RC parameters

plt.rcParams.update(RC_PARAMS)

# Make colors

alphas = np.arange(0.5, 0.751, 0.025)

cmap = plt.get_cmap("Paired", len(alphas))

# Plot potentials

E = 1
phi = np.linspace(-0.3, 1.3, 1000)

for i, alpha in enumerate(alphas):
    p = potential(E, alpha)
    label = r"$\alpha = {}$".format(alpha)
    plt.plot(phi, -p(phi), label=label, color=cmap(i), lw=3)

# Finish up

plt.xlim([-0.25, 1.25])
plt.ylim([-0.05, 0.3])
loc = plticker.MultipleLocator(base=0.25)
plt.gca().xaxis.set_major_locator(loc)
plt.grid(True)

handles, labels = plt.gca().get_legend_handles_labels()
leg = plt.legend(handles[::-1], labels[::-1], loc='upper left')
leg.get_frame().set_alpha(0.75)

plt.xlabel(r'Field, $\phi$')
plt.ylabel(r'Upturned potential, $-V(\phi)$')
plt.title(r'$V(\phi) = \frac{-4 \alpha + 3}{2} \phi^2 - \phi^3 + \alpha \phi^4$')

plt.tight_layout()
plt.savefig('figures/quartic_potential.pdf')
