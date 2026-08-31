# Detailed derivation for the wave--ODE low-gain proof

This note is the calculation record for the proof in
`wave_ode_low_gain.tex`. All integrals without displayed limits are over
`[0,D]`. The derivation keeps the symbols already used in the paper:
`X`, `Z`, `u`, `u_t`, `g`, `g_0`, `\mathcal B`, `\bar{\mathcal B}`,
`E_u`, `V_1`, `V_2`, `V`, and `\Omega`. No auxiliary scalar is introduced.

## 1. Kernel identities and the auxiliary state

The general kernel equations are

\[
 g''(x)+B_1(x)=A\bigl(Ag(x)-B_2(x)\bigr),
 \qquad g_0(x)=Ag(x)-B_2(x).
\]

Therefore,

\[
 g''(x)+B_1(x)=Ag_0(x),
 \qquad B_2(x)+g_0(x)=Ag(x).
\tag{1}
\]

The auxiliary state and its difference from the original ODE state are

\[
 Z(t)=X(t)+\int_0^D g_0(x)u(x,t)\,dx
 +\int_0^D g(x)u_t(x,t)\,dx
 +b\mathcal B u(D,t),
\tag{2}
\]

\[
 \varphi(t)=Z(t)-X(t)
 =\int_0^D g_0(x)u(x,t)\,dx
 +\int_0^D g(x)u_t(x,t)\,dx
 +b\mathcal B u(D,t).
\tag{3}
\]

Differentiate (2) along a classical trajectory:

\[
\begin{aligned}
 \dot Z={}&\dot X+\int_0^Dg_0(x)u_t(x,t)\,dx
 +\int_0^Dg(x)u_{tt}(x,t)\,dx
 +b\mathcal B u_t(D,t).
\end{aligned}
\tag{4}
\]

Use the ODE and the wave equation in (4):

\[
\begin{aligned}
 \dot Z={}&AX+\int_0^D B_1(x)u(x,t)\,dx
 +\int_0^D B_2(x)u_t(x,t)\,dx\\
 &+\int_0^D g_0(x)u_t(x,t)\,dx
 +\int_0^D g(x)u_{xx}(x,t)\,dx
 +b\mathcal B u_t(D,t).
\end{aligned}
\tag{5}
\]

The wave term is integrated by parts twice. First,

\[
 \int_0^Dg u_{xx}\,dx
 =\bigl[g u_x\bigr]_0^D-\int_0^Dg'u_x\,dx.
\tag{6}
\]

For the remaining integral,

\[
 -\int_0^Dg'u_x\,dx
 =-\bigl[g'u\bigr]_0^D+\int_0^Dg''u\,dx.
\tag{7}
\]

Combining (6) and (7), and then using
`g'(0)=0`, `u_x(0,t)=0`, and `g(D)=\mathcal B`, gives

\[
\begin{aligned}
 \int_0^Dg u_{xx}\,dx
 & =g(D)u_x(D,t)-g(0)u_x(0,t)\\
 &\quad-g'(D)u(D,t)+g'(0)u(0,t)+\int_0^Dg''u\,dx\\
 & =\mathcal B u_x(D,t)-g'(D)u(D,t)+\int_0^Dg''u\,dx.
\end{aligned}
\tag{8}
\]

Insert (8) into (5), group the displacement and velocity integrals, and
apply (1):

\[
\begin{aligned}
 \dot Z={}&AX+\int_0^D\bigl(B_1+g''\bigr)u\,dx
 +\int_0^D\bigl(B_2+g_0\bigr)u_t\,dx\\
 &+\mathcal B u_x(D,t)-g'(D)u(D,t)+b\mathcal B u_t(D,t)\\
={}&AX+A\int_0^Dg_0u\,dx+A\int_0^Dg u_t\,dx\\
 &+\mathcal B u_x(D,t)-g'(D)u(D,t)+b\mathcal B u_t(D,t).
\end{aligned}
\tag{9}
\]

The right kernel boundary condition is

\[
 -g'(D)=(bA+qI)g(D)=(bA+qI)\mathcal B.
\tag{10}
\]

Using (2) and (10) in (9), the terms proportional to
`bA\mathcal B u(D,t)` cancel:

\[
\begin{aligned}
 \dot Z
 &=A\bigl(Z-b\mathcal B u(D,t)\bigr)
 +\mathcal B u_x(D,t)\\
 &\quad +(bA+qI)\mathcal B u(D,t)+b\mathcal B u_t(D,t)\\
 &=AZ+\mathcal B\bigl(u_x(D,t)+bu_t(D,t)+qu(D,t)\bigr).
\end{aligned}
\tag{11}
\]

For the general controller,

\[
 u_x(D,t)=U(t)=-bu_t(D,t)-qu(D,t)-\mathcal B^\top P_\varepsilon X(t).
\tag{12}
\]

Consequently,

\[
 \dot Z=AZ-\mathcal B\mathcal B^\top P_\varepsilon X
 =(A-\mathcal B\mathcal B^\top P_\varepsilon)Z
 +\mathcal B\mathcal B^\top P_\varepsilon\varphi.
\tag{13}
\]

The last equality uses `X=Z-\varphi`.

## 2. Trace and auxiliary-state estimates

For every `x` in `[0,D]`,

\[
 u(x,t)=u(D,t)-\int_x^D u_x(\zeta,t)\,d\zeta.
\]

The elementary square inequality and Cauchy--Schwarz give

\[
\begin{aligned}
 |u(x,t)|^2
 &\leq 2u^2(D,t)+2(D-x)\int_x^D u_x^2(\zeta,t)\,d\zeta\\
 &\leq 2u^2(D,t)+2D\int_0^D u_x^2(\zeta,t)\,d\zeta.
\end{aligned}
\]

Since `E_u` contains `q u^2(D,t)/2` and
`q>0`,

\[
 u^2(D,t)\leq \frac{2}{q}E_u(t),
 \qquad
 \int_0^D u_x^2\,dx\leq 2E_u(t).
\]

Thus

\[
 \|u(\cdot,t)\|_{L^\infty}^2
 \leq 4(D+q^{-1})E_u(t)=C_\infty E_u(t).
\tag{14}
\]

Integration over `[0,D]` gives

\[
 \|u(\cdot,t)\|^2\leq DC_\infty E_u(t)=C_0E_u(t).
\tag{15}
\]

Apply Cauchy--Schwarz separately to the three terms in (3), followed by
the three-term square inequality:

\[
\begin{aligned}
 |\varphi|^2
 &\leq 3\left|\int_0^Dg_0u\,dx\right|^2
 +3\left|\int_0^Dg u_t\,dx\right|^2
 +3b^2|\mathcal B|^2u^2(D,t)\\
 &\leq 3\|g_0\|^2\|u\|^2
 +3\|g\|^2\int_0^D u_t^2\,dx
 +3b^2|\mathcal B|^2u^2(D,t)\\
 &\leq \left(3\|g_0\|^2C_0+6\|g\|^2
 +\frac{6b^2|\mathcal B|^2}{q}\right)E_u(t).
\end{aligned}
\tag{16}
\]

Therefore the constant used in the paper is

\[
 C_\varphi=3\|g_0\|^2C_0+6\|g\|^2
 +\frac{6b^2|\mathcal B|^2}{q},
 \qquad |\varphi|^2\leq C_\varphi E_u.
\tag{17}
\]

## 3. Wave energy and multiplier derivatives

The wave energy is

\[
 E_u=\frac12\int_0^D(u_x^2+u_t^2)\,dx+\frac q2u^2(D,t).
\]

Differentiate it and use `u_{tt}=u_{xx}`:

\[
\begin{aligned}
 \dot E_u
 &=\int_0^D u_xu_{xt}\,dx+\int_0^D u_tu_{tt}\,dx
 +q u(D,t)u_t(D,t)\\
 &=\int_0^D u_xu_{xt}\,dx+\int_0^D u_tu_{xx}\,dx
 +q u(D,t)u_t(D,t).
\end{aligned}
\]

Integrate the second spatial term once:

\[
 \int_0^D u_tu_{xx}\,dx
 =\bigl[u_tu_x\bigr]_0^D-\int_0^D u_{tx}u_x\,dx.
\]

The two interior integrals cancel. Since `u_x(0,t)=0` and
`u_x(D,t)=U(t)`,

\[
 \dot E_u=U(t)u_t(D,t)+q u(D,t)u_t(D,t).
\tag{18}
\]

With the general feedback (12),

\[
 \dot E_u=-b u_t^2(D,t)
 -\bigl(\mathcal B^\top P_\varepsilon X(t)\bigr)u_t(D,t).
\tag{19}
\]

For the first multiplier, differentiate before integrating by parts:

\[
\begin{aligned}
 \frac{d}{dt}\int_0^Dxu_xu_t\,dx
 &=\int_0^Dx u_{xt}u_t\,dx+\int_0^Dx u_xu_{tt}\,dx\\
 &=\int_0^Dx u_{xt}u_t\,dx+\int_0^Dx u_xu_{xx}\,dx.
\end{aligned}
\]

The first term is

\[
 \int_0^Dx u_{xt}u_t\,dx
 =\frac12\int_0^Dx\partial_x(u_t^2)\,dx
 =\frac D2u_t^2(D,t)-\frac12\int_0^D u_t^2\,dx.
\]

The second term is

\[
 \int_0^Dx u_xu_{xx}\,dx
 =\frac12\int_0^Dx\partial_x(u_x^2)\,dx
 =\frac D2u_x^2(D,t)-\frac12\int_0^D u_x^2\,dx,
\]

because `x=0` removes the left endpoint. Hence

\[
 \frac{d}{dt}\int_0^Dxu_xu_t\,dx
 =\frac D2\bigl(U^2(t)+u_t^2(D,t)\bigr)
 -\frac12\int_0^D(u_x^2+u_t^2)\,dx.
\tag{20}
\]

For the second multiplier,

\[
\begin{aligned}
 \frac{d}{dt}\int_0^Duu_t\,dx
 &=\int_0^D u_t^2\,dx+\int_0^Duu_{tt}\,dx\\
 &=\int_0^D u_t^2\,dx+\int_0^Duu_{xx}\,dx\\
 &=\int_0^D u_t^2\,dx+\bigl[uu_x\bigr]_0^D
 -\int_0^D u_x^2\,dx\\
 &=\int_0^D(u_t^2-u_x^2)\,dx+u(D,t)U(t).
\end{aligned}
\tag{21}
\]

The left endpoint vanishes because `u_x(0,t)=0`.

## 4. Complete expansion of the wave functional

Choose

\[
 \frac{Dq}{2}<\kappa<\frac12,
 \qquad \delta_2=\kappa\delta_1.
\tag{22}
\]

The wave functional is

\[
 V_2=E_u+\delta_1\int_0^Dxu_xu_t\,dx
 +\delta_2\int_0^Duu_t\,dx.
\tag{23}
\]

Substitute (19)--(21) into the derivative of (23):

\[
\begin{aligned}
 \dot V_2={}&-b u_t^2(D,t)
 -\bigl(\mathcal B^\top P_\varepsilon X\bigr)u_t(D,t)\\
 &+\frac{\delta_1D}{2}\bigl(U^2(t)+u_t^2(D,t)\bigr)
 -\frac{\delta_1}{2}\int_0^D(u_x^2+u_t^2)\,dx\\
 &+\delta_2\int_0^D(u_t^2-u_x^2)\,dx
 +\delta_2u(D,t)U(t).
\end{aligned}
\tag{24}
\]

For the boundary square, use (12) and expand every term:

\[
\begin{aligned}
 U^2(t)={}&b^2u_t^2(D,t)+q^2u^2(D,t)
 +\bigl(\mathcal B^\top P_\varepsilon X\bigr)^2\\
 &+2bq u_t(D,t)u(D,t)\\
 &+2b u_t(D,t)\bigl(\mathcal B^\top P_\varepsilon X\bigr)\\
 &+2q u(D,t)\bigl(\mathcal B^\top P_\varepsilon X\bigr).
\end{aligned}
\tag{25}
\]

The last term in (24) is

\[
\begin{aligned}
 \delta_2u(D,t)U(t)={}&-\delta_2b u(D,t)u_t(D,t)
 -\delta_2q u^2(D,t)\\
 &-\delta_2u(D,t)\bigl(\mathcal B^\top P_\varepsilon X\bigr).
\end{aligned}
\tag{26}
\]

Collecting (24)--(26), using `\delta_2=\kappa\delta_1`, gives

\[
\begin{aligned}
 \dot V_2={}&-\frac{\delta_1(1+2\kappa)}2
       \int_0^D u_x^2\,dx
 -\frac{\delta_1(1-2\kappa)}2
       \int_0^D u_t^2\,dx\\
 &-\left[b-\frac{\delta_1D(1+b^2)}2\right]u_t^2(D,t)\\
 &-\delta_1q\left(\kappa-\frac{Dq}{2}\right)u^2(D,t)\\
 &+\delta_1b(Dq-\kappa)u(D,t)u_t(D,t)\\
 &-(1-\delta_1Db)\bigl(\mathcal B^\top P_\varepsilon X\bigr)u_t(D,t)\\
 &-\delta_1(\kappa-Dq)u(D,t)\bigl(\mathcal B^\top P_\varepsilon X\bigr)\\
 &+\frac{\delta_1D}{2}\bigl(\mathcal B^\top P_\varepsilon X\bigr)^2.
\end{aligned}
\tag{27}
\]

The signs in (27) follow directly from the four boundary contributions:
the energy identity, the square in (25), the mixed term `2bq`, and (26).

## 5. Boundary Young inequalities and the constants `c_3` and `C_2`

Choose `\delta_1>0` so that

\[
 b-\frac{\delta_1D(1+b^2)}2>0.
\tag{28}
\]

The first mixed boundary term in (27) satisfies

\[
\begin{aligned}
 \left|\delta_1b(Dq-\kappa)u(D,t)u_t(D,t)\right|
 &\leq \frac14\left[b-\frac{\delta_1D(1+b^2)}2\right]u_t^2(D,t)\\
 &\quad+\frac{\delta_1^2b^2(Dq-\kappa)^2}
 {b-\frac{\delta_1D(1+b^2)}2}u^2(D,t).
\end{aligned}
\tag{29}
\]

This is Young's inequality with the positive coefficient in (28) as the
weight; it leaves one quarter of that weighted boundary velocity square
and the displayed reciprocal-weight displacement square.

The mixed velocity--ODE term satisfies the same inequality:

\[
\begin{aligned}
 &(1-\delta_1Db)
 \left|\bigl(\mathcal B^\top P_\varepsilon X\bigr)u_t(D,t)\right|\\
 &\leq \frac14\left[b-\frac{\delta_1D(1+b^2)}2\right]u_t^2(D,t)\\
 &\quad+\frac{(1-\delta_1Db)^2}
 {b-\frac{\delta_1D(1+b^2)}2}
 \bigl(\mathcal B^\top P_\varepsilon X\bigr)^2.
\end{aligned}
\tag{30}
\]

After (29)--(30), the remaining negative coefficient of `u^2(D,t)` is
bounded below by

\[
 \delta_1q\left(\kappa-\frac{Dq}{2}\right)
 -\frac{\delta_1^2b^2(Dq-\kappa)^2}
 {b-\frac{\delta_1D(1+b^2)}2}.
\tag{31}
\]

Choose `\delta_1` small enough that (31) is positive. The last mixed
displacement--ODE term is estimated by

\[
\begin{aligned}
 &\delta_1(\kappa-Dq)
 \left|u(D,t)\mathcal B^\top P_\varepsilon X\right|\\
 &\leq \frac12\left[\delta_1q\left(\kappa-\frac{Dq}{2}\right)
 -\frac{\delta_1^2b^2(Dq-\kappa)^2}
 {b-\frac{\delta_1D(1+b^2)}2}\right]u^2(D,t)\\
 &\quad+\frac{\delta_1^2(\kappa-Dq)^2}
 {2\left[\delta_1q\left(\kappa-\frac{Dq}{2}\right)
 -\frac{\delta_1^2b^2(Dq-\kappa)^2}
 {b-\frac{\delta_1D(1+b^2)}2}\right]}
 \bigl(\mathcal B^\top P_\varepsilon X\bigr)^2.
\end{aligned}
\tag{32}
\]

The boundary velocity term left after (29)--(30) is negative. Therefore
the three negative coefficients in front of the two integrals and the
boundary displacement can be summarized by

\[
\begin{aligned}
 c_3=\min\Biggl\{&\delta_1(1+2\kappa),\quad
 \delta_1(1-2\kappa),\\
 &\frac1q\left[\delta_1q\left(\kappa-\frac{Dq}{2}\right)
 -\frac{\delta_1^2b^2(Dq-\kappa)^2}
 {b-\frac{\delta_1D(1+b^2)}2}\right]\Biggr\}>0.
\end{aligned}
\tag{33}
\]

Indeed, (33) implies

\[
\begin{aligned}
 &-\frac{\delta_1(1+2\kappa)}2\int_0^Du_x^2\,dx
 -\frac{\delta_1(1-2\kappa)}2\int_0^Du_t^2\,dx\\
 &\quad-\frac{c_3q}{2}u^2(D,t)
 \leq -c_3E_u(t).
\end{aligned}
\tag{34}
\]

The coefficient multiplying
`(\mathcal B^\top P_\varepsilon X)^2` after (30) and (32) is

\[
\begin{aligned}
 C_2={}&\frac{\delta_1D}{2}
 +\frac{(1-\delta_1Db)^2}
 {b-\frac{\delta_1D(1+b^2)}2}\\
 &+\frac{\delta_1^2(\kappa-Dq)^2}
 {2\left[\delta_1q\left(\kappa-\frac{Dq}{2}\right)
 -\frac{\delta_1^2b^2(Dq-\kappa)^2}
 {b-\frac{\delta_1D(1+b^2)}2}\right]}.
\end{aligned}
\tag{35}
\]

Equations (27)--(35) yield

\[
 \dot V_2\leq-c_3E_u+C_2
 \bigl(\mathcal B^\top P_\varepsilon X\bigr)^2.
\tag{36}
\]

## 6. Cauchy--Schwarz estimates for the ODE cross term

Since `X=Z-\varphi`,

\[
\begin{aligned}
 \bigl(\mathcal B^\top P_\varepsilon X\bigr)^2
 &\leq 2\bigl(\mathcal B^\top P_\varepsilon Z\bigr)^2
 +2\bigl(\mathcal B^\top P_\varepsilon\varphi\bigr)^2\\
 &\leq 2|\mathcal B|^2|P_\varepsilon Z|^2
 +2|\mathcal B|^2|P_\varepsilon|^2|\varphi|^2.
\end{aligned}
\tag{37}
\]

For a symmetric positive definite matrix `P_\varepsilon`,

\[
 P_\varepsilon^2\preceq |P_\varepsilon|P_\varepsilon,
 \qquad
 |P_\varepsilon Z|^2
 =Z^\top P_\varepsilon^2Z
 \leq |P_\varepsilon|Z^\top P_\varepsilon Z
 =|P_\varepsilon|V_1.
\tag{38}
\]

Using (17), (37), and (38) in (36),

\[
 \dot V_2\leq-c_3E_u+C_3|P_\varepsilon|V_1
 +C_4|P_\varepsilon|^2E_u,
\tag{39}
\]

where

\[
 C_3=2C_2|\mathcal B|^2,
 \qquad
 C_4=2C_2|\mathcal B|^2C_\varphi.
\tag{40}
\]

## 7. General Riccati functional

Let

\[
 V_1=Z^\top P_\varepsilon Z.
\]

Differentiating and using (13),

\[
\begin{aligned}
 \dot V_1={}&Z^\top(A^\top P_\varepsilon+P_\varepsilon A)Z
 -2Z^\top P_\varepsilon\mathcal B\mathcal B^\top P_\varepsilon Z\\
 &+2Z^\top P_\varepsilon\mathcal B\mathcal B^\top P_\varepsilon\varphi.
\end{aligned}
\tag{41}
\]

The algebraic Riccati equation gives

\[
 A^\top P_\varepsilon+P_\varepsilon A
 =-\varepsilon P_\varepsilon
 +P_\varepsilon\mathcal B\mathcal B^\top P_\varepsilon.
\tag{42}
\]

Substitution into (41) gives

\[
\begin{aligned}
 \dot V_1={}&-\varepsilon V_1
 -\bigl(\mathcal B^\top P_\varepsilon Z\bigr)^2\\
 &+2\bigl(\mathcal B^\top P_\varepsilon Z\bigr)
       \bigl(\mathcal B^\top P_\varepsilon\varphi\bigr).
\end{aligned}
\tag{43}
\]

The scalar inequality obtained by completing the square, applied directly
to the two displayed inner products in (43), yields

\[
 \dot V_1\leq-\varepsilon V_1
 +\bigl(\mathcal B^\top P_\varepsilon\varphi\bigr)^2.
\]

Cauchy--Schwarz and (17) then give

\[
 \dot V_1\leq-\varepsilon V_1+C_1|P_\varepsilon|^2E_u,
 \qquad
 C_1=|\mathcal B|^2C_\varphi.
\tag{44}
\]

## 8. Equivalence of `V_2` and `E_u`

The first multiplier term satisfies

\[
\begin{aligned}
 \left|\delta_1\int_0^Dxu_xu_t\,dx\right|
 &\leq \frac{\delta_1D}{2}\int_0^D(u_x^2+u_t^2)\,dx\\
 &\leq \delta_1D E_u.
\end{aligned}
\tag{45}
\]

For the second multiplier term, Cauchy--Schwarz and the two-term Young
inequality give

\[
\begin{aligned}
 \left|\delta_2\int_0^Duu_t\,dx\right|
 &\leq \frac{\delta_2}{2}
 \left(\int_0^Du^2\,dx+\int_0^Du_t^2\,dx\right)\\
 &\leq \kappa\delta_1\left(\frac{C_0}{2}+1\right)E_u.
\end{aligned}
\tag{46}
\]

Define the already used constant

\[
 K_\delta=D+\kappa\left(\frac{C_0}{2}+1\right).
\tag{47}
\]

Equations (45)--(46) imply

\[
 (1-\delta_1K_\delta)E_u
 \leq V_2\leq (1+\delta_1K_\delta)E_u.
\tag{48}
\]

For `0<\delta_1\leq(2K_\delta)^{-1}`, set

\[
 m_\delta=1-\delta_1K_\delta>0,
 \qquad M_\delta=1+\delta_1K_\delta.
\tag{49}
\]

Then `m_\delta E_u\leq V_2\leq M_\delta E_u`.

## 9. Total functional and the general small-gain bound

Set

\[
 V=V_1+\varepsilon V_2.
\tag{50}
\]

Combining (39) and (44),

\[
\begin{aligned}
 \dot V\leq{}&-\varepsilon V_1+C_1|P_\varepsilon|^2E_u
 +\varepsilon\left[-c_3E_u+C_3|P_\varepsilon|V_1
 +C_4|P_\varepsilon|^2E_u\right]\\
={}&-\varepsilon(1-C_3|P_\varepsilon|)V_1\\
 &-\left[\varepsilon c_3
 -(C_1+\varepsilon C_4)|P_\varepsilon|^2\right]E_u.
\end{aligned}
\tag{51}
\]

The Riccati estimate is `|P_\varepsilon|\leq K_P\varepsilon` for
`0<\varepsilon\leq\varepsilon_0`. The condition

\[
 \varepsilon\leq\frac{1}{2C_3K_P}
\tag{52}
\]

gives `C_3|P_\varepsilon|\leq1/2`. The two separate conditions

\[
 \varepsilon\leq\frac{c_3}{4C_1K_P^2},
 \qquad
 \varepsilon\leq\sqrt{\frac{c_3}{4C_4K_P^2}}
\tag{53}
\]

give, respectively,

\[
 C_1|P_\varepsilon|^2\leq\frac{\varepsilon c_3}{4},
 \qquad
 \varepsilon C_4|P_\varepsilon|^2\leq\frac{\varepsilon c_3}{4}.
\tag{54}
\]

If one of `C_1` or `C_4` is zero, the corresponding restriction is
omitted. Thus the explicit sufficient choice is

\[
 \varepsilon^*=\min\left\{
 \varepsilon_0,\frac{1}{2C_3K_P},
 \frac{c_3}{4C_1K_P^2},
 \sqrt{\frac{c_3}{4C_4K_P^2}}
 \right\},
\tag{55}
\]

with zero-denominator terms omitted. Equations (51)--(54) imply

\[
 \dot V\leq-\frac{\varepsilon}{2}V_1
 -\frac{\varepsilon c_3}{2}E_u.
\tag{56}
\]

Because `V_2\leq M_\delta E_u`,

\[
 V=V_1+\varepsilon V_2
 \leq V_1+\varepsilon M_\delta E_u.
\]

Therefore (56) implies

\[
 \dot V\leq-\beta V,
 \qquad
 \beta=\min\left\{\frac{\varepsilon}{2},
 \frac{c_3}{2M_\delta}\right\}>0.
\tag{57}
\]

The two terms in the minimum separately absorb `V_1` and
`\varepsilon M_\delta E_u`.

## 10. Explicit equivalence of `V` and `\Omega`

First, `Z=X+\varphi` implies

\[
 |X|^2\leq2|Z|^2+2|\varphi|^2.
\tag{58}
\]

Also,

\[
 \int_0^D(u_x^2+u_t^2)\,dx+u^2(D,t)
 \leq 2(1+q^{-1})E_u.
\tag{59}
\]

Since `V_1\geq\lambda_{\min}(P_\varepsilon)|Z|^2`,
(17), (58), and (59) give

\[
 \Omega\leq 2\lambda_{\min}(P_\varepsilon)^{-1}V_1
 +\bigl[2C_\varphi+2(1+q^{-1})\bigr]E_u.
\tag{60}
\]

From `V_1\leq V` and `E_u\leq V_2/m_\delta\leq V/(\varepsilon m_\delta)`,

\[
 \Omega\leq
 \max\left\{2\lambda_{\min}(P_\varepsilon)^{-1},
 \frac{2C_\varphi+2(1+q^{-1})}{\varepsilon m_\delta}\right\}V.
\tag{61}
\]

Consequently,

\[
 m_\varepsilon\Omega\leq V,
\quad
 m_\varepsilon=\left[
 \max\left\{2\lambda_{\min}(P_\varepsilon)^{-1},
 \frac{2C_\varphi+2(1+q^{-1})}{\varepsilon m_\delta}\right\}
 \right]^{-1}.
\tag{62}
\]

For the upper bound, use `Z=X+\varphi`, (17), and (49):

\[
\begin{aligned}
 V_1&\leq |P_\varepsilon||Z|^2\\
 &\leq2|P_\varepsilon||X|^2
 +2C_\varphi|P_\varepsilon|E_u,\\
 \varepsilon V_2&\leq\varepsilon M_\delta E_u.
\end{aligned}
\tag{63}
\]

Moreover,

\[
 E_u\leq \max\left\{\frac12,\frac q2\right\}
 \left[\int_0^D(u_x^2+u_t^2)\,dx+u^2(D,t)\right]
 \leq \max\left\{\frac12,\frac q2\right\}\Omega.
\tag{64}
\]

Combining (63)--(64),

\[
 V\leq M_\varepsilon\Omega,
\tag{65}
\]

where

\[
 M_\varepsilon=2|P_\varepsilon|
 +\left(2C_\varphi|P_\varepsilon|+\varepsilon M_\delta\right)
 \max\left\{\frac12,\frac q2\right\}.
\tag{66}
\]

## 11. Gronwall estimate and well-posedness

Integrating (57), or applying the differential form of Gronwall's
inequality, gives

\[
 V(t)\leq e^{-\beta t}V(0).
\tag{67}
\]

Use (62) and (65) at times `t` and `0`:

\[
 \Omega(t)\leq m_\varepsilon^{-1}V(t)
 \leq \frac{M_\varepsilon}{m_\varepsilon}
 e^{-\beta t}\Omega(0).
\tag{68}
\]

Thus the theorem holds with

\[
 \alpha=\frac{M_\varepsilon}{m_\varepsilon}.
\]

For classical compatible initial data, the wave equation with the
feedback boundary condition is a finite-dimensional boundary perturbation
of the standard Neumann wave generator. The distributed maps

\[
 u\mapsto\int_0^DB_1(x)u(x,t)\,dx,
 \qquad
 u_t\mapsto\int_0^DB_2(x)u_t(x,t)\,dx
\]

are bounded on `H^1(0,D)` and `L^2(0,D)`, respectively, by
Cauchy--Schwarz. The boundary trace `u(D,t)` is bounded on `H^1(0,D)`,
and the feedback is therefore a bounded boundary feedback on the energy
domain. Standard semigroup perturbation theory gives existence and
uniqueness. Energy solutions follow by density from compatible classical
data, and (68) passes to the limit.

## 12. Nilpotent specialization

For a nilpotent `A`, use the direct-cancellation kernel

\[
 g''(x)+B_1(x)=0,
 \qquad g'(0)=0,
 \qquad -g'(D)=qg(D),
 \qquad g_0(x)=-B_2(x).
\tag{69}
\]

Integrating `g''=-B_1` from `0` to `x` and using `g'(0)=0` gives

\[
 g'(x)=-\int_0^xB_1(\zeta)\,d\zeta.
\tag{70}
\]

At `x=D`, (69) and (70) give

\[
 qg(D)=\int_0^DB_1(\zeta)\,d\zeta,
 \qquad
 \bar{\mathcal B}=g(D)=q^{-1}\int_0^DB_1(\zeta)\,d\zeta.
\tag{71}
\]

Repeating (4)--(9), but now using
`g''+B_1=0` and `B_2+g_0=0`, gives

\[
 \dot Z=AX-g'(D)u(D,t)+\bar{\mathcal B}u_x(D,t)
 +b\bar{\mathcal B}u_t(D,t).
\tag{72}
\]

The boundary condition in (69) gives `-g'(D)=q\bar{\mathcal B}`; hence

\[
 \dot Z=AX+\bar{\mathcal B}\bigl(u_x(D,t)+bu_t(D,t)+qu(D,t)\bigr).
\tag{73}
\]

With

\[
 U=-bu_t(D,t)-qu(D,t)-\bar{\mathcal B}^{\top}P_\varepsilon X(t),
\]

equation (73) becomes

\[
 \dot Z=(A-\bar{\mathcal B}\bar{\mathcal B}^{\top}P_\varepsilon)X
 =(A-\bar{\mathcal B}\bar{\mathcal B}^{\top}P_\varepsilon)(Z-\varphi).
\tag{74}
\]

The estimates for `E_u`, `V_2`, `c_3`, and the norm-equivalence constants
are unchanged after replacing `\mathcal B` by `\bar{\mathcal B}`.

## 13. Nilpotent Riccati estimate and the extra `A\varphi` term

For the nilpotent case,

\[
 V_1=Z^\top P_\varepsilon Z,
 \qquad
 \dot Z=(A-\bar{\mathcal B}\bar{\mathcal B}^{\top}P_\varepsilon)(Z-\varphi).
\]

Direct differentiation gives

\[
\begin{aligned}
 \dot V_1={}&2Z^\top P_\varepsilon A Z
 -2Z^\top P_\varepsilon A\varphi\\
 &-2\bigl(\bar{\mathcal B}^{\top}P_\varepsilon Z\bigr)^2
 +2\bigl(\bar{\mathcal B}^{\top}P_\varepsilon Z\bigr)
       \bigl(\bar{\mathcal B}^{\top}P_\varepsilon\varphi\bigr).
\end{aligned}
\tag{75}
\]

The Riccati equation transforms the first term in (75):

\[
 2Z^\top P_\varepsilon A Z
 =-\varepsilon V_1
 +\bigl(\bar{\mathcal B}^{\top}P_\varepsilon Z\bigr)^2.
\tag{76}
\]

The two terms involving `\bar{\mathcal B}` satisfy

\[
\begin{aligned}
 &-\bigl(\bar{\mathcal B}^{\top}P_\varepsilon Z\bigr)^2
 +2\bigl(\bar{\mathcal B}^{\top}P_\varepsilon Z\bigr)
       \bigl(\bar{\mathcal B}^{\top}P_\varepsilon\varphi\bigr)\\
 &\qquad\leq
 \bigl(\bar{\mathcal B}^{\top}P_\varepsilon\varphi\bigr)^2.
\end{aligned}
\tag{77}
\]

For the additional term, apply the weighted Young inequality with weights
`\varepsilon/2` and `2/\varepsilon` to the two inner products generated by
`P_\varepsilon`:

\[
\begin{aligned}
 2|Z^\top P_\varepsilon A\varphi|
 &\leq \frac{\varepsilon}{2}Z^\top P_\varepsilon Z
 +\frac{2}{\varepsilon}\varphi^\top A^\top P_\varepsilon A\varphi\\
 &\leq \frac{\varepsilon}{2}V_1
 +2C_A\varepsilon|P_\varepsilon||\varphi|^2.
\end{aligned}
\tag{78}
\]

The second line uses

\[
 A^\top P_\varepsilon A\preceq C_A\varepsilon^2P_\varepsilon,
 \qquad
 \varphi^\top P_\varepsilon\varphi
 \leq |P_\varepsilon||\varphi|^2.
\]

Combining (75)--(78), Cauchy--Schwarz, and (17),

\[
\begin{aligned}
 \dot V_1\leq{}&-\frac{\varepsilon}{2}V_1
 +\left(|\bar{\mathcal B}|^2|P_\varepsilon|^2
 +2C_A\varepsilon|P_\varepsilon|\right)|\varphi|^2\\
\leq{}&-\frac{\varepsilon}{2}V_1
 +\left(|\bar{\mathcal B}|^2C_\varphi|P_\varepsilon|^2
 +2C_AC_\varphi\varepsilon|P_\varepsilon|\right)E_u.
\end{aligned}
\tag{79}
\]

## 14. Nilpotent small-gain selection

Use (39) with `\mathcal B` replaced by `\bar{\mathcal B}` and add
`\varepsilon` times that estimate to (79):

\[
\begin{aligned}
 \dot V\leq{}&-\frac{\varepsilon}{2}V_1
 +\varepsilon C_3|P_\varepsilon|V_1\\
 &+\left(2C_AC_\varphi\varepsilon|P_\varepsilon|
 +|\bar{\mathcal B}|^2C_\varphi|P_\varepsilon|^2
 +\varepsilon C_4|P_\varepsilon|^2
 -\varepsilon c_3\right)E_u.
\end{aligned}
\tag{80}
\]

The first restriction

\[
 \varepsilon\leq\frac{1}{4C_3K_P}
\tag{81}
\]

gives `\varepsilon C_3|P_\varepsilon|\leq\varepsilon/4`.
For the three positive terms in the coefficient of `E_u`, impose

\[
\varepsilon\leq\frac{c_3}{12C_AC_\varphi K_P},
\qquad
\varepsilon\leq\frac{c_3}
 {6|\bar{\mathcal B}|^2C_\varphi K_P^2},
\qquad
\varepsilon\leq\sqrt{\frac{c_3}{6C_4K_P^2}}.
\tag{82}
\]

Indeed, `|P_\varepsilon|\leq K_P\varepsilon` changes the three terms
into quantities bounded by `\varepsilon c_3/6` each. Hence the explicit
nilpotent bound is

\[
 \varepsilon_n^*=\min\left\{
 \varepsilon_0,\frac{1}{4C_3K_P},
 \frac{c_3}{12C_AC_\varphi K_P},
 \frac{c_3}{6|\bar{\mathcal B}|^2C_\varphi K_P^2},
 \sqrt{\frac{c_3}{6C_4K_P^2}}
 \right\},
\tag{83}
\]

again omitting zero-denominator terms. Equations (80)--(83) imply

\[
 \dot V\leq-\frac{\varepsilon}{4}V_1
 -\frac{\varepsilon c_3}{2}E_u.
\tag{84}
\]

Since `V\leq V_1+\varepsilon M_\delta E_u`,

\[
 \dot V\leq-\beta_nV,
 \qquad
 \beta_n=\min\left\{\frac{\varepsilon}{4},
 \frac{c_3}{2M_\delta}\right\}.
\tag{85}
\]

The norm-equivalence proof (58)--(68) is unchanged after replacing
`\mathcal B` by `\bar{\mathcal B}`. Gronwall's inequality therefore gives
the same exponential estimate for `\Omega`.

## 15. Order analysis and interpretation

The Riccati lemma supplies

\[
 |P_\varepsilon|=O(\varepsilon).
\]

Consequently,

\[
 |P_\varepsilon|^2E_u=O(\varepsilon^2)E_u,
 \qquad
 \varepsilon|P_\varepsilon|V_1=O(\varepsilon^2)V_1,
\]

while the principal dissipative terms in (56) are
`O(\varepsilon)V_1` and `O(\varepsilon)E_u`. The bounds (55) and (83)
make the higher-order terms no larger than fixed fractions of those
principal terms. This is why the result is a low-gain statement: the
proof does not require a large boundary feedback coefficient, only a
small enough `\varepsilon` for the explicit absorption inequalities.
