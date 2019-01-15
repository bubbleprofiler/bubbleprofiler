from __future__ import (division, print_function)

import matplotlib.ticker as tck
import matplotlib.pyplot as plt
import numpy as np

class SSMModel:
    def __init__(self, g_=0, gPrime_=0, h_t_=0,
                 lambda_h_=0, lambda_m_=0, lambda_s_=0,
                 vEW_=0, Tc_=0):
        self.g = g_
        self.gPrime = gPrime_
        self.h_t = h_t_
        self.lambda_h = lambda_h_
        self.lambda_m = lambda_m_
        self.lambda_s = lambda_s_
        self.vEW = vEW_
        self.Tc = Tc_

    def c_h(self):
        return ((9 * self.g**2 + 3 * self.gPrime**2
                + 2 * (6 * self.h_t + 12 * self.lambda_h + self.lambda_m))
                / 48.0)
    def c_s(self):
        return (2 * self.lambda_m + 3 * self.lambda_s) / 12.0

    def v(self, temp):
        return np.sqrt(self.vEW**2 - (self.c_h() / self.lambda_h) * temp**2)

    def m_h(self, temp):
        return np.sqrt(2.0 * self.lambda_h * self.v(temp)**2)

    def m_s_Tc(self):
        return np.sqrt(0.5 * self.v(self.Tc)**2
                       * (self.lambda_m
                          - 2 * np.sqrt(self.lambda_h * self.lambda_s)))

    def m_s(self, temp):
        m_s_Tc = self.m_s_Tc()
        c_h = self.c_h()
        c_s = self.c_s()
        return np.sqrt(m_s_Tc**2
                       + (0.5 * self.lambda_m * c_h / self.lambda_h - c_s)
                       * (self.Tc**2 - temp**2))

    def lambda_e_sq(self):
        m_s_Tc = self.m_s_Tc()
        return ((m_s_Tc**4 - self.lambda_m * self.v(self.Tc)**2 * m_s_Tc**2)
                / (self.v(self.Tc)**4))

    def w(self, temp):
        return np.sqrt((self.m_h(temp)**2 * (self.lambda_m *self.v(temp)**2
                                             - 2 * self.m_s(temp)**2)) /
                       (self.v(temp)**2 * (4 *self.lambda_e_sq()
                                           + self.lambda_m**2)))

    def Rh(self):
        return self.m_h(self.Tc)**2 / self.w(self.Tc)**2

    def Rs(self):
        return self.m_s_Tc()**2 / self.v(self.Tc)**2

    def order_param(self):
        return self.v(self.Tc) / self.Tc

    def VT(self, h, s, temp=0):
        m_h_Tc = self.m_h(self.Tc)
        v_Tc = self.v(self.Tc)
        Rs = self.Rs()
        Rh = self.Rh()
        w0 = self.w(self.Tc)
        c_h = self.c_h()
        c_s = self.c_s()

        return (0.125 * m_h_Tc**2 * v_Tc**2 * (
            4 * Rs * h**2 * s**2 / (Rh * v_Tc**2 * w0**2)
            + (h**2 / v_Tc**2 + s**2 / w0**2 - 1)**2)
                - 0.5 * (self.Tc**2 - temp**2) *
                (c_h * h**2 + c_s * s**2))

def get_benchmark_model():
    vEW = 246.22
    Mh = 125.1
    Mw = 80.385
    Mz = 91.1876
    Mt = 173.34

    g = (2.0*Mw)/vEW
    gPrime = np.sqrt((4.0*(Mz**2 - Mw**2))/(vEW**2))
    h_t = (np.sqrt(2)*Mt)/vEW
    lambda_h = (Mh**2)/(2.0*vEW**2)

    lambda_m = 1.5
    lambda_s = 0.65
    Tc = 110

    return SSMModel(g, gPrime, h_t, lambda_h, lambda_m, lambda_s, vEW, Tc)

def plot_effective_potential(model, temp, log10_vnorm=0,
                             contour_levels=None,
                             title="", outfile=""):
    vT = model.v(temp)
    wT = model.w(temp)

    n_points = 50
    h_scale = 1.4
    s_scale = 1.5

    h_range = np.linspace(0, h_scale * vT, n_points)
    s_range = np.linspace(-s_scale * wT, s_scale * wT, n_points)

    h_vals, s_vals = np.meshgrid(h_range, s_range)
    potential_vals = model.VT(h_vals, s_vals, temp)
    potential_norm = 10**log10_vnorm

    plt.rc("text", usetex=True)
    plt.rc("font", family="serif", weight="normal")
    plt.rcParams["text.latex.preamble"] = [r"\usepackage{amsmath}"]

    fig = plt.figure(figsize=(4.5,4))

    plt.gcf().subplots_adjust(bottom=0.15, left=0.2)
    ax = plt.gca()
    ax.set_axisbelow(True)
    ax.tick_params(direction="in", which="both", labelsize=10)
    ax.xaxis.set_ticks_position("both")
    ax.yaxis.set_ticks_position("both")

    if contour_levels is not None:
        cs = plt.contourf(h_vals / vT, s_vals / wT,
                          potential_vals / potential_norm,
                          contour_levels, alpha=0.9,
                          cmap="viridis")
    else:
        cs = plt.contourf(h_vals / vT, s_vals / wT,
                          potential_vals / potential_norm, 10,
                          alpha=0.9,
                          cmap="viridis")

    plt.xlim([0, h_scale])
    plt.ylim([-s_scale, s_scale])

    plt.xlabel(r"$h / v(T)$", fontsize=12)
    plt.ylabel(r"$s / w(T)$", fontsize=12)

    if title:
        plt.title(title, fontsize=12)

    cb = plt.colorbar(cs)
    cb.ax.set_ylabel(r"$V_T\;/\;10^" + str(log10_vnorm)
                     + r"\;\mathrm{GeV}^4$",
                     fontsize=12)

    if outfile:
        plt.savefig(outfile)

    plt.close(fig)

def main():
    benchmark = get_benchmark_model()
    TN = 85
    Tc = benchmark.Tc

    log10_vnorm_T0 = 7
    log10_vnorm_TN = 7
    log10_vnorm_Tc = 7

    contours_T0 = [-11, -9.5, -8, -6.5, -5, -3.5, -2, -0.5,
                   1, 2.5, 4]
    contours_TN = [-3, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1]
    contours_Tc = [-0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]

    plot_effective_potential(benchmark, 0,
                             log10_vnorm=log10_vnorm_T0,
                             contour_levels=contours_T0,
                             title=r"$T = 0$", outfile="figures/ssm_eff_pot_T0.pdf")

    plot_effective_potential(benchmark, TN,
                             log10_vnorm=log10_vnorm_TN,
                             contour_levels=contours_TN,
                             title=r"$T = T_N$", outfile="figures/ssm_eff_pot_TN.pdf")

    plot_effective_potential(benchmark, Tc,
                             log10_vnorm=log10_vnorm_Tc,
                             contour_levels=contours_Tc,
                             title=r"$T = T_c$", outfile="figures/ssm_eff_pot_Tc.pdf")

if __name__ == "__main__":
    main()
