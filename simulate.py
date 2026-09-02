"""Numerical simulations for the three examples in main.tex.

Closed-loop wave PDE--ODE cascade with
    u_tt = u_xx,  u_x(0,t)=0,  u_x(D,t) = U(t),
    dX = A X + int B1 u dx + int B2 u_t dx,
    U  = -K X - b u_t(D) - q u(D),   K = B^T P_eps (general) or Bbar^T P_eps (nilpotent).

The script produces the figures cited in the numerical-simulation section.
Only numpy and matplotlib are used.
"""
import os
import site
import sys

# The user-site Matplotlib installation is incompatible with the system's
# mpl_toolkits package. Use the matching system packages so 3-D axes work.
user_site = site.getusersitepackages()
if isinstance(user_site, str) and user_site in sys.path:
    sys.path.remove(user_site)

import matplotlib
import numpy as np
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

HERE = os.path.dirname(os.path.abspath(__file__))

# ----------------------------------------------------------------------------
# basic linear algebra helpers (no scipy)
# ----------------------------------------------------------------------------

def expm(A, tol=1e-14):
    """Matrix exponential via scaling-and-squaring + Taylor (works for any square A)."""
    A = np.asarray(A, dtype=float)
    m = A.shape[0]
    nrm = np.max(np.abs(A)) if A.size else 0.0
    s = max(0, int(np.ceil(np.log2(max(nrm, 1e-300)))))
    B = A / (2.0 ** s)
    term = np.eye(m)
    R = np.eye(m)
    for k in range(1, 120):
        term = term @ B / k
        R = R + term
        if np.max(np.abs(term)) < tol:
            break
    for _ in range(s):
        R = R @ R
    return R


def solve_are_lowgain(A, B, eps):
    """Unique pos-def solution P of  A^T P + P A - P B B^T P = -eps P.

    Solved via the equivalent Lyapunov equation for W = P^{-1}:
        W M^T + M W = B B^T,   M = A + (eps/2) I.
    """
    n = A.shape[0]
    M = A + 0.5 * eps * np.eye(n)
    # vec(W M^T) = (M kron I) vec(W),  vec(M W) = (I kron M) vec(W)
    T = np.kron(M, np.eye(n)) + np.kron(np.eye(n), M)
    rhs = (B @ B.T).reshape(-1, order="F")
    W = np.linalg.solve(T, rhs).reshape((n, n), order="F")
    W = 0.5 * (W + W.T)
    P = np.linalg.inv(W)
    P = 0.5 * (P + P.T)
    # verify the ARE residual
    res = A.T @ P + P @ A - P @ (B @ B.T) @ P + eps * P
    resid = np.max(np.abs(res))
    return P, resid


# ----------------------------------------------------------------------------
# kernel BVP (general case):  g'' + B1 = A(A g - B2), g'(0)=0, -g'(D)=(bA+qI)g(D)
# ----------------------------------------------------------------------------

def solve_kernel(A, B1fn, B2fn, b, q, D, Nk=4000):
    """Return (B = g(D), eta = g(0)) for the kernel BVP.

        g'' + B1 = A(A g - B2),   g'(0) = 0,   -g'(D) = (bA + qI) g(D).

    Solved in the state-space form  zeta' = Abar zeta + F with
    zeta = (g, g'),  Abar = [[0,I],[A^2,0]],  F = (0, -A B2 - B1).
    """
    n = A.shape[0]
    x = np.linspace(0.0, D, Nk + 1)
    h = D / Nk
    w = np.ones(Nk + 1)
    w[0] = 0.5
    w[-1] = 0.5
    w *= h
    B1 = B1fn(x)          # (Nk+1, n)
    B2 = B2fn(x)          # (Nk+1, n)
    Abar = np.zeros((2 * n, 2 * n))
    Abar[:n, n:] = np.eye(n)
    Abar[n:, :n] = A @ A
    F = np.vstack([np.zeros((n, Nk + 1)), -(A @ B2.T) - B1.T])
    PhiD = expm(Abar * D)
    integ = np.zeros(2 * n)
    for j in range(Nk + 1):
        integ = integ + w[j] * (expm(Abar * (D - x[j])) @ F[:, j])
    L = np.hstack([b * A + q * np.eye(n), np.eye(n)])
    Mmat = L @ PhiD @ np.vstack([np.eye(n), np.zeros((n, n))])
    rvec = L @ integ
    eta = np.linalg.solve(Mmat, -rvec)
    B = (PhiD @ np.concatenate([eta, np.zeros(n)]))[:n] + integ[:n]
    return B, eta


# ----------------------------------------------------------------------------
# time integrator (method of lines, RK4)
# ----------------------------------------------------------------------------

def simulate(n, A, B1mat, B2mat, K, b, q, D, X0, u0, v0, T=50.0, dt=5e-4,
             N=200, nskip=50, damping_only=False):
    """Semi-discrete wave + ODE, RK4.  B1mat/B2mat: (N+1, n) on the grid."""
    dx = D / N
    dx2 = dx * dx
    xg = np.linspace(0.0, D, N + 1)
    w = np.ones(N + 1)
    w[0] = 0.5
    w[-1] = 0.5
    w *= dx
    Np1 = N + 1
    Kk = K if (not damping_only) else np.zeros_like(K)
    bb = b
    qq = q if (not damping_only) else 0.0

    def rhs(y):
        X = y[:n]
        u = y[n:n + Np1]
        v = y[n + Np1:]
        U = -(Kk @ X) - bb * v[-1] - qq * u[-1]
        dX = A @ X
        for i in range(n):
            dX[i] += np.dot(w * B1mat[:, i], u) + np.dot(w * B2mat[:, i], v)
        du = v
        dv = np.zeros(Np1)
        dv[1:N] = (u[2:] - 2.0 * u[1:-1] + u[:-2]) / dx2
        dv[0] = 2.0 * (u[1] - u[0]) / dx2
        dv[-1] = 2.0 * (u[-2] - u[-1] + dx * U) / dx2
        return np.concatenate([dX, du, dv])

    y = np.concatenate([np.asarray(X0, float), np.asarray(u0, float),
                        np.asarray(v0, float)])
    nt = int(round(T / dt))
    rec = []
    times = []
    for step in range(nt + 1):
        t = step * dt
        if step % nskip == 0:
            X = y[:n]
            u = y[n:n + Np1]
            rec.append(np.concatenate([X, u]))
            times.append(t)
        k1 = rhs(y)
        k2 = rhs(y + 0.5 * dt * k1)
        k3 = rhs(y + 0.5 * dt * k2)
        k4 = rhs(y + dt * k3)
        y = y + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
    rec = np.array(rec)
    times = np.array(times)
    return times, rec


# ----------------------------------------------------------------------------
# figures
# ----------------------------------------------------------------------------

def plot_states(axis, times, X, title, ylabel=True):
    """Draw the two ODE coordinates with a consistent panel style."""
    for i in range(X.shape[1]):
        axis.plot(times, X[:, i], lw=1.2, label=fr"$x_{i + 1}$")
    axis.set_title(title, fontsize=9)
    axis.set_xlabel("Time (s)", fontsize=8)
    if ylabel:
        axis.set_ylabel("$x_i(t)$", fontsize=8)
    axis.grid(alpha=0.3)
    axis.tick_params(labelsize=7)
    axis.legend(fontsize=7, frameon=False, loc="best")


def plot_surface(axis, times, xg, u, title):
    """Draw a three-dimensional surface of the PDE displacement."""
    sampled_times = times[::10]
    sampled_x = xg[::4]
    time_mesh, x_mesh = np.meshgrid(sampled_times, sampled_x)
    surface = u[::10, ::4].T
    axis.plot_surface(time_mesh, x_mesh, surface, cmap="viridis", rstride=1,
                      cstride=1, linewidth=0, antialiased=True, shade=True)
    axis.set_title(title, fontsize=9, pad=4)
    axis.set_xlabel("Time (s)", fontsize=7, labelpad=0)
    axis.set_ylabel("$x$", fontsize=7, labelpad=0)
    axis.set_zlabel("$u(x,t)$", fontsize=7, labelpad=2)
    axis.tick_params(labelsize=6, pad=0)
    axis.view_init(elev=25, azim=-125)
    axis.ticklabel_format(style="sci", axis="z", scilimits=(0, 0))


def fig_example1(times, xg, X, u, path):
    """Example 1: closed-loop ODE trajectories and a 3-D wave response."""
    fig = plt.figure(figsize=(7.1, 3.15))
    grid = fig.add_gridspec(1, 2, wspace=0.25)
    plot_states(fig.add_subplot(grid[0, 0]), times, X, "(a) ODE state response")
    plot_surface(fig.add_subplot(grid[0, 1], projection="3d"), times, xg, u,
                 "(b) Wave displacement")
    fig.subplots_adjust(left=0.07, right=0.97, bottom=0.14, top=0.89)
    fig.savefig(path, dpi=180)
    plt.close(fig)


def fig_comparison(times, xg, Xcl, ucl, Xdmp, udmp, path):
    """Example 2: ODE and PDE comparison for damping-only and full feedback."""
    fig = plt.figure(figsize=(7.1, 5.7))
    grid = fig.add_gridspec(2, 2, wspace=0.25, hspace=0.35)
    plot_states(fig.add_subplot(grid[0, 0]), times, Xdmp,
                "(a) Damping-only ODE response")
    plot_surface(fig.add_subplot(grid[0, 1], projection="3d"), times, xg, udmp,
                 "(b) Damping-only wave response")
    plot_states(fig.add_subplot(grid[1, 0]), times, Xcl,
                "(c) Proposed-controller ODE response")
    plot_surface(fig.add_subplot(grid[1, 1], projection="3d"), times, xg, ucl,
                 "(d) Proposed-controller wave response")
    fig.subplots_adjust(left=0.07, right=0.97, bottom=0.06, top=0.94)
    fig.savefig(path, dpi=180)
    plt.close(fig)


def fig_epsilon(times, xg, responses, path):
    """Example 3: one ODE/PDE comparison figure for all epsilon values."""
    fig = plt.figure(figsize=(7.1, 8.7))
    grid = fig.add_gridspec(3, 2, wspace=0.25, hspace=0.36)
    for row, (eps, X, u) in enumerate(responses):
        plot_states(fig.add_subplot(grid[row, 0]), times, X,
                    fr"({chr(ord('a') + 2 * row)}) ODE response, $\varepsilon={eps:.1f}$")
        plot_surface(fig.add_subplot(grid[row, 1], projection="3d"), times, xg, u,
                     fr"({chr(ord('b') + 2 * row)}) Wave response, $\varepsilon={eps:.1f}$")
    fig.subplots_adjust(left=0.07, right=0.97, bottom=0.04, top=0.96)
    fig.savefig(path, dpi=180)
    plt.close(fig)


# ----------------------------------------------------------------------------
# Example 1 : general marginally stable oscillator
# ----------------------------------------------------------------------------
def general_case(N=200):
    """Shared model, grid, and initial data for Examples 1 and 3."""
    D, b, q, eps = 1.0, 1.0, 0.4, 0.12
    A = np.array([[0.0, -1.0], [1.0, 0.0]])

    def B1fn(x):
        g = np.exp(-30.0 * (x - 0.7) ** 2)
        return np.outer(g, np.array([1.0, 0.4]))

    def B2fn(x):
        g = np.exp(-20.0 * (x - 0.35) ** 2)
        return np.outer(g, np.array([1.0, 0.25]))

    xg = np.linspace(0.0, D, N + 1)
    X0 = np.array([0.40, -0.2301])
    u0 = 0.08 * np.cos(np.pi * xg / (2.0 * D))
    v0 = np.zeros_like(xg)
    return D, b, q, eps, A, B1fn, B2fn, xg, X0, u0, v0


def ex1():
    """General marginally stable oscillator--wave cascade."""
    D, b, q, eps, A, B1fn, B2fn, xg, X0, u0, v0 = general_case()
    B, eta = solve_kernel(A, B1fn, B2fn, b, q, D)
    rank = np.linalg.matrix_rank(np.hstack([B.reshape(-1, 1), A @ B.reshape(-1, 1)]))
    P, residual = solve_are_lowgain(A, B.reshape(-1, 1), eps)
    K = (B.reshape(1, -1) @ P).reshape(-1)
    times, rec = simulate(2, A, B1fn(xg), B2fn(xg), K, b, q, D, X0, u0, v0,
                          N=xg.size - 1)
    print("EX1 kernel: eta(g(0)) =", eta, " B = g(D) =", B, " rank[A,AB] =", rank,
          " ARE-resid =", residual)
    print("EX1 P =", P)
    print("EX1 K = B^T P =", K)
    fig_example1(times, xg, rec[:, :2], rec[:, 2:], os.path.join(HERE, "fig2.png"))


# ----------------------------------------------------------------------------
# Example 2 : nilpotent chain
# ----------------------------------------------------------------------------
def ex2():
    D, b, q, eps = 1.0, 1.0, 0.5, 0.20
    A = np.array([[0.0, 1.0], [0.0, 0.0]])
    Bbar = np.array([0.0, 2.0])  # (1/q) * int_0^D [0;1] dx  = (1/0.5)*[0;1]
    P, res2 = solve_are_lowgain(A, Bbar.reshape(-1, 1), eps)
    K = (Bbar.reshape(1, -1) @ P).reshape(-1)
    N = 200
    xg = np.linspace(0.0, D, N + 1)
    B1mat = np.outer(np.ones(N + 1), np.array([0.0, 1.0]))
    B2mat = np.zeros((N + 1, 2))
    X0 = np.array([0.30, -0.20])
    u0 = 0.04 * np.cos(np.pi * xg / (2.0 * D))
    v0 = np.zeros_like(xg)

    tcl, rcl = simulate(n=2, A=A, B1mat=B1mat, B2mat=B2mat, K=K,
                         b=b, q=q, D=D, X0=X0, u0=u0, v0=v0, N=N)
    td, rd = simulate(n=2, A=A, B1mat=B1mat, B2mat=B2mat, K=K,
                      b=b, q=q, D=D, X0=X0, u0=u0, v0=v0, N=N,
                      damping_only=True)
    print("EX2 Bbar =", Bbar, " P =", P, " ARE-resid =", res2)
    print("EX2 K = Bbar^T P =", K)
    fig_comparison(tcl, xg, rcl[:, :2], rcl[:, 2:], rd[:, :2], rd[:, 2:],
                   os.path.join(HERE, "fig3.png"))


# ----------------------------------------------------------------------------
# Example 3 : effect of the low-gain parameter
# ----------------------------------------------------------------------------
def ex3():
    D, b, q, _, A, B1fn, B2fn, xg, X0, u0, v0 = general_case()
    B, _ = solve_kernel(A, B1fn, B2fn, b, q, D)
    responses = []
    for eps in (0.2, 0.5, 0.9):
        P, residual = solve_are_lowgain(A, B.reshape(-1, 1), eps)
        K = (B.reshape(1, -1) @ P).reshape(-1)
        times, rec = simulate(n=2, A=A, B1mat=B1fn(xg), B2mat=B2fn(xg), K=K,
                              b=b, q=q, D=D, X0=X0, u0=u0, v0=v0,
                              N=xg.size - 1)
        responses.append((eps, rec[:, :2], rec[:, 2:]))
        print(f"EX3 eps = {eps:.1f}: K = {K}, ARE-resid = {residual}")
    fig_epsilon(times, xg, responses, os.path.join(HERE, "fig4.png"))


def main():
    ex1()
    ex2()
    ex3()
    print("DONE")


if __name__ == "__main__":
    main()
