"""Real simulation for main.tex Examples 1 and 2.

Closed-loop wave PDE--ODE cascade with
    u_tt = u_xx,  u_x(0,t)=0,  u_x(D,t) = U(t),
    dX = A X + int B1 u dx + int B2 u_t dx,
    U  = -K X - b u_t(D) - q u(D),   K = B^T P_eps (general) or Bbar^T P_eps (nilpotent).

Pure-numpy implementation (scipy is broken against numpy 2.2.6 on this box).
Only numpy + matplotlib are used.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

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
    uxsq = np.zeros(Np1)
    for step in range(nt + 1):
        t = step * dt
        if step % nskip == 0:
            X = y[:n]
            u = y[n:n + Np1]
            v = y[n + Np1:]
            U = -(Kk @ X) - bb * v[-1] - qq * u[-1]
            # u_x on grid: left Neumann -> 0, right -> U, interior central
            uxsq[0] = 0.0
            uxsq[-1] = U * U
            uxsq[1:N] = ((u[2:] - u[:-2]) / (2.0 * dx)) ** 2
            Iux2 = np.dot(w, uxsq)
            Iut2 = np.dot(w, v * v)
            Om = np.dot(X, X) + Iux2 + Iut2 + u[-1] ** 2
            rec.append(np.concatenate([X, u, v, [U, Om]]))
            times.append(t)
        k1 = rhs(y)
        k2 = rhs(y + 0.5 * dt * k1)
        k3 = rhs(y + 0.5 * dt * k2)
        k4 = rhs(y + dt * k3)
        y = y + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
    rec = np.array(rec)
    times = np.array(times)
    return times, rec, n, Np1, dx


# ----------------------------------------------------------------------------
# indicators
# ----------------------------------------------------------------------------

def indicators(times, rec, n, Np1, dx):
    X = rec[:, :n]
    u = rec[:, n:n + Np1]
    v = rec[:, n + Np1:2 * n + Np1]
    U = rec[:, -2]
    Om = rec[:, -1]
    Om0 = Om[0]
    ratio = Om[-1] / Om0
    KX = X.copy()  # placeholder for gain*X (computed by caller if needed)
    # least-squares decay rate over the window where ln Omega is roughly linear
    mask = (times >= 5.0) & (Om > 0)
    if mask.sum() > 2:
        tt = times[mask]
        ll = np.log(np.clip(Om[mask], 1e-300, None))
        rate = np.polyfit(tt, ll, 1)[0]
        rate = -rate
    else:
        rate = np.nan
    # crossing times
    fracs = [0.10, 0.05, 0.01]
    crosses = []
    for f in fracs:
        idx = np.where(Om <= f * Om0)[0]
        crosses.append(times[idx[0]] if idx.size else np.nan)
    return dict(ratio=ratio, rate=rate, crosses=crosses, Om0=Om0,
                X50=X[-1], uD50=u[-1, -1], uD=u[:, -1], vD=v[:, -1])


# ----------------------------------------------------------------------------
# figures
# ----------------------------------------------------------------------------

def fig_states(times, Xcl, Xdmp, labels, path):
    fig, ax = plt.subplots(1, 2, figsize=(7.0, 2.6))
    for a, X, lab in zip(ax, [Xdmp, Xcl], ["damping only", "closed loop"]):
        for i in range(X.shape[1]):
            a.plot(times, X[:, i], lw=1.2, label=f"$x_{i+1}$")
        a.set_title(lab, fontsize=9)
        a.set_xlabel("$t$", fontsize=9)
        a.grid(alpha=0.3)
        if a is ax[0]:
            a.set_ylabel("$X_i(t)$", fontsize=9)
    fig.legend(loc="upper center", ncol=2, fontsize=8, frameon=False)
    fig.suptitle(labels, fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(path, dpi=150)
    plt.close(fig)


def fig_surface(times, xg, u, path, step=4):
    """3-D surface of u(x, t)."""
    U = u[::step, ::step]
    T = times[::step]
    Xg = xg[::step]
    fig = plt.figure(figsize=(5.0, 3.0))
    ax = fig.add_subplot(111, projection="3d")
    Xm, Tm = np.meshgrid(Xg, T)
    ax.plot_surface(Tm, Xm, U, cmap="viridis", rstride=1, cstride=1,
                    linewidth=0, antialiased=True)
    ax.set_xlabel("$t$", fontsize=8)
    ax.set_ylabel("$x$", fontsize=8)
    ax.set_zlabel("$u(x,t)$", fontsize=8)
    ax.tick_params(labelsize=6)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


# ----------------------------------------------------------------------------
# Example 1 : general marginally stable oscillator
# ----------------------------------------------------------------------------
def ex1():
    D, b, q, eps = 1.0, 1.0, 0.4, 0.12
    A = np.array([[0.0, -1.0], [1.0, 0.0]])

    def B1fn(x):
        g = np.exp(-30.0 * (x - 0.7) ** 2)
        return np.outer(g, np.array([1.0, 0.4]))

    def B2fn(x):
        g = np.exp(-20.0 * (x - 0.35) ** 2)
        return np.outer(g, np.array([1.0, 0.25]))

    B, eta = solve_kernel(A, B1fn, B2fn, b, q, D)
    # controllability
    Cm = np.hstack([B.reshape(-1, 1), A @ B.reshape(-1, 1)])
    rank = np.linalg.matrix_rank(Cm)
    P, res1 = solve_are_lowgain(A, B.reshape(-1, 1), eps)
    K = (B.reshape(1, -1) @ P).reshape(-1)
    # grid for PDE
    N = 200
    xg = np.linspace(0.0, D, N + 1)
    B1mat = B1fn(xg)
    B2mat = B2fn(xg)
    X0 = np.array([0.40, -0.2301])
    u0 = 0.08 * np.cos(np.pi * xg / (2.0 * D))
    v0 = np.zeros_like(xg)

    tcl, rcl, n, Np1, dx = simulate(n=2, A=A, B1mat=B1mat, B2mat=B2mat, K=K,
                                    b=b, q=q, D=D, X0=X0, u0=u0, v0=v0)
    td, rd, _, _, _ = simulate(n=2, A=A, B1mat=B1mat, B2mat=B2mat, K=K,
                               b=b, q=q, D=D, X0=X0, u0=u0, v0=v0,
                               damping_only=True)
    ind = indicators(tcl, rcl, 2, N + 1, dx)
    indd = indicators(td, rd, 2, N + 1, dx)
    KXcl = rcl[:, :2] @ K
    maxKX = np.max(np.abs(KXcl))
    print("EX1 kernel: eta(g(0)) =", eta, " B = g(D) =", B, " rank[A,AB] =", rank,
          " ARE-resid =", res1)
    print("EX1 P =", P)
    print("EX1 K = B^T P =", K)
    print("EX1 Omega(50)/Omega(0) closed =", ind["ratio"],
          "  damping-only =", indd["ratio"])
    print("EX1 decay rate =", ind["rate"])
    print("EX1 crosses 10/5/1% =", ind["crosses"])
    print("EX1 max|KX| =", maxKX)
    print("EX1 |X(50)| =", np.linalg.norm(ind["X50"]), " |u(D,50)| =", ind["uD50"])
    print("EX1 Omega(0) =", ind["Om0"])

    fig_states(tcl, rcl[:, :2], rd[:, :2], "Example 1: ODE states",
               os.path.join(HERE, "fig2.png"))
    fig_surface(tcl, xg, rcl[:, 2:2 + N + 1], os.path.join(HERE, "fig3.png"))
    return ind, indd, K, B, P


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

    tcl, rcl, n, Np1, dx = simulate(n=2, A=A, B1mat=B1mat, B2mat=B2mat, K=K,
                                    b=b, q=q, D=D, X0=X0, u0=u0, v0=v0)
    td, rd, _, _, _ = simulate(n=2, A=A, B1mat=B1mat, B2mat=B2mat, K=K,
                               b=b, q=q, D=D, X0=X0, u0=u0, v0=v0,
                               damping_only=True)
    ind = indicators(tcl, rcl, 2, N + 1, dx)
    indd = indicators(td, rd, 2, N + 1, dx)
    KXcl = rcl[:, :2] @ K
    maxKX = np.max(np.abs(KXcl))
    print("EX2 Bbar =", Bbar, " P =", P, " ARE-resid =", res2)
    print("EX2 K = Bbar^T P =", K)
    print("EX2 Omega(50)/Omega(0) closed =", ind["ratio"],
          "  damping-only =", indd["ratio"])
    print("EX2 decay rate =", ind["rate"])
    print("EX2 crosses 10/5/1% =", ind["crosses"])
    print("EX2 max|KX| =", maxKX)
    print("EX2 |X(50)| =", np.linalg.norm(ind["X50"]), " |u(D,50)| =", ind["uD50"])
    print("EX2 Omega(0) =", ind["Om0"])

    fig_states(tcl, rcl[:, :2], rd[:, :2], "Example 2: nilpotent ODE states",
               os.path.join(HERE, "fig4.png"))
    fig_surface(tcl, xg, rcl[:, 2:2 + N + 1], os.path.join(HERE, "fig5.png"))
    return ind, indd, K, Bbar, P


if __name__ == "__main__":
    ind1, indd1, K1, B1, P1 = ex1()
    ind2, indd2, K2, Bb2, P2 = ex2()
    print("DONE")
