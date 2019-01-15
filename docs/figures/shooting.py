"""
Illustration of shooting method
===============================
"""

from matplotlib import pyplot as plt

import matplotlib.patches as patches
import numpy as np


RC_PARAMS = {'text.latex.preamble' : [r'\usepackage{amsmath}'],
             'font.family': 'lmodern',
             'font.size': 30,
             'axes.titlesize': 20,
             'xtick.labelsize': 24,
             'ytick.labelsize': 24,
             'axes.labelsize': 30,
             'legend.fontsize': 20,
             'text.usetex': True,
             'lines.linewidth': 3,
             'axes.linewidth': 1,
             'axes.grid': True,
             'grid.linestyle': '--',
             'legend.framealpha': 1.,
             'legend.edgecolor': 'black',
             'savefig.bbox': 'tight'}


# Potential

def minus_potential(phi):
    E = 10.
    alpha = 0.55
    phi = 1. - phi
    return -E * ((-4. * alpha + 3.) * 0.5 * phi**2 - phi**3 + alpha * phi**4)


# Starting and resting place of field

phi_0 = 0.25
phi_tv = 0.
phi_fv = 1.
phi_b = 2. - 0.75 / 0.55
x_shift = 0.1

# Plot potential

phi = np.linspace(phi_tv - x_shift, phi_fv + x_shift, 10000)

plt.rcParams.update(RC_PARAMS)
fig, ax = plt.subplots(figsize=(16, 8))

# Overshoot region

over = (phi < phi_0) & (phi > phi_tv)
plt.plot(phi[over], minus_potential(phi[over]), color="Coral", lw=10, zorder=1, label=r'Overshoot --- go past $\phi_f$')

# Undershoot region

under = (phi > phi_0) & (phi < phi_fv)
plt.plot(phi[under], minus_potential(phi[under]), color="Maroon", lw=10, zorder=1, label=r'Undershoot --- roll into $\phi_b$')

# Solution

ax.axvline(phi_0, ls="--", label=r"Solution --- stop at $\phi_f$", color="SteelBlue", alpha=0.8, lw=5, zorder=2)

# Rest of line

left = phi < phi_tv - 0.1 * x_shift
right = phi > phi_fv + 0.1 * x_shift
plt.plot(phi[left], minus_potential(phi[left]), color="black", lw=10, alpha=0.5, zorder=1, ls="--")
plt.plot(phi[right], minus_potential(phi[right]), color="black", lw=10, alpha=0.5, zorder=1, ls="--")

# Plot and annotate arrow

y_shift = 0.05
start_coord = (phi_0 + 0.01, minus_potential(phi_0) + y_shift)
end_coord = (phi_fv, minus_potential(phi_fv) + y_shift)
p = patches.FancyArrowPatch(start_coord, end_coord, connectionstyle='arc3, rad=0.1', mutation_scale=35, fc='k', ec='k', zorder=5)
ax.add_patch(p)
plt.annotate('Field rolls and stops at false vacuum', xy=(0.5, 0.2))

# Plot rolling field

for frac in np.linspace(0., 1., 30):
    alpha = 0.9 - frac * 0.5
    angle = frac * 360
    marker = (7, 0, angle)
    phi = frac * abs(phi_fv - phi_0) + phi_0
    # plt.scatter(phi, minus_potential(phi), marker=marker, c="Crimson", edgecolors='k', s=300, zorder=3, alpha=alpha)


# Finish up

plt.xlim([phi_tv - x_shift, phi_fv + x_shift,])
plt.ylim([-0.2, 0.55])
plt.xticks([phi_tv, phi_0, phi_b, phi_fv], [r"True vacuum, $\phi_t$", r"Solution, $\phi(\rho = 0) = \phi_0$", r"Barrier, $\phi_b$", r"False vacuum, $\phi_f$"])
plt.yticks([], [])
ax.tick_params(axis='both', which='major', pad=15)
plt.grid(True)

xticks = plt.getp(ax, 'xticklines')
plt.setp(xticks, markeredgewidth=4, markersize=10)

leg = plt.legend(scatterpoints=1, loc='upper right')
leg.get_frame().set_alpha(0.75)
plt.xlabel(r'Field $\phi$')
plt.ylabel(r'Upturned potential, $-V(\phi)$')
plt.title('')

plt.tight_layout()
plt.savefig('figures/shooting.pdf')
