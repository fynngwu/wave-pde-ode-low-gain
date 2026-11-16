# Wave PDE–ODE 论文改写大纲 v8

> **定位**：以 Xu 等人在 *IEEE Transactions on Automatic Control* 发表的 *Low-Gain Compensation for PDE–ODE Cascade Systems With Distributed Diffusion and Counter Convection* 的 Section III 为行文与排版参照，保留本文 wave PDE 的数学模型和证法，而不照搬 diffusion PDE 的技术步骤。
>
> **v7 相比 v6 的变化**：**只修数学排版，不改变 v6 已确定的英文文风、控制器、两条定理和证明结论**。把状态所属空间、核函数的可积性、示例中的 $\varepsilon$ 取值列表等描述性数学条件融入句子，以行内公式排版；保留 $\mathfrak B$、$E$、$\Omega$、$Z$、$E_u$、$V_1$、$V_2$、$V$ 等定义的独立展示。主结果里的控制器、能量函数、Lyapunov 导数等在不超出 IEEE 双栏宽度的前提下优先把同类项写在同一行，只有真正过长时才换行。v6 的独立核方程编号、完整自变量、Assumption 3.1、控制器及 ARE 顺序全部保持。

> **v8 相比 v7 的变化**：仅对 Theorem 3.1 的 Proof 作三处局部调整，其他行文、符号和证明顺序不变。第一，在定义核函数后，连续展示一条导数恒等式和一条从 $0$ 到 $D$ 的逐行积分推导，然后显式得出 $g(D)=\mathfrak B$；保留每条核边值条件各自独立编号。第二，为能量等价不等式中的 $m_\varepsilon,M_\varepsilon$ 给出明确、可核验的表达式。第三，最后的 Consequently 只展示最终指数衰减不等式，不再按 $V(t)$ 和 $V(0)$ 拆成三行。遵循 Xu Section III 的写法：必要的代数计算展示清楚，已成立的推论尽量简短。

## 一、总体编排和风格

1. **保留章节顺序**：I. Introduction；II. Preliminaries；III. Main Results；IV. Numerical Simulation；V. Conclusion。必要的半群良定性证明放在 Appendix。
2. **像 Xu 论文一样写 III 节**：简短开场 → 参数 → 可计算的有效向量 → Assumption 3.1 → 第一个控制器及 ARE → $\Omega$ → Theorem 3.1 及 Proof → nilpotent 特例 → Assumption 3.2 → 第二个控制器及 ARE → Theorem 3.2 及 Proof → Remarks。不设置 General Case / Nilpotent Case / Discussion 等子小节。
3. **两套有效输入向量保持区别**：一般情形只有 $\mathfrak B$，nilpotent 特例用 $\bar B$。**不再使用 $B:=\mathfrak B$ 这一重命名**。Preliminaries 的低增益引理可以使用通用矩阵对 $(A,B)$，但主结果和例 1 中必须写 $\mathfrak B$。
4. **简单、清楚的英语**：优先使用 “All eigenvalues of $A$ are on the imaginary axis”, “the matrix pair is controllable”, “Define”, “Then”, “It follows that”；不要在主结果开头插入矩阵函数、Jordan 块、无 $A^{-1}$ 等解释。这些信息统一留到 Remark 3.4。
5. **定义独立成式，描述性条件留在行内**：$E,\mathfrak B,\bar B,\Omega,Z,\varphi,E_u,V_1,V_2,V,\kappa,\beta$ 等定义及重要衰减结论用 `equation`。但诸如 $X(t)\in\mathbb R^n$、初值 $(X(0),u(\cdot,0),u_t(\cdot,0))\in\mathbb R^n\times H^1(0,D)\times L^2(0,D)$、$g\in H^2(0,D;\mathbb R^n)$、$B_i\in L^2$、$\varepsilon\in\{\cdots\}$ 等，不单独占用一个公式，直接嵌入句子。关键参数选取条件 $b>0,0<q<1/D$ 保留单独展示，以便引用。独立编号的核边值问题仍然**逐条编号**。
6. **以单栏宽度控制换行，不是越多越好**：一个公式尽量合并两个或更多短项，例如 $U(t)=-\mathfrak B^\top P_\varepsilon X(t)-b\,u_t(D,t)-q\,u(D,t)$ 用一行；短的同类项可以同一行，完整的积分项和长的交叉项无法容纳时再用 `aligned` 分为两行。**禁止**为了强行单行使用整体缩小字号、`\resizebox` 或造成 overfull box；同时不要给每个加减项都另起一行。
7. **时空变量不省略**：在整个证明的推导中始终写 $\dot Z(t)$、$X(t)$、$Z(t)$、$u(x,t)$、$u_t(x,t)$、$u_x(x,t)$、$u(D,t)$、$u_t(D,t)$；若双栏无法容纳，可删除重复且非必要的中间恒等式，**不可通过省略自变量来挤进一栏**。时间固定时的范数可以写 $\|u(\cdot,t)\|_2$。
8. **减少显示公式造成的留白**：正文中的自然语言衔接行内数学时不停顿。连续推导用一个 `aligned` 环境组织，不为一个纯文字补充说明单开 display；每次定义仍紧跟其推导，不随意移动位置或更换编号方式。
9. **只做必要解释**：核心算法的可计算性、矩阵可逆性、Jordan 块、在线测量等放在 Remark 中；Proof 中紧邻定义给出计算和等式。无外加公式方框。

## 二、章节安排

### I. Introduction

原有六段背景、相关工作、Xu low-gain 的区别、控制器预告、贡献、论文安排保留。唯一与 III 节同步的记号修正，是将控制器预告中的 $B$ 改为 $\mathfrak B$：

\[
U(t)=-\mathfrak B^\top P_\varepsilon X(t)-b\,u_t(D,t)-q\,u(D,t).
\]

摘要中 *two controllers* 对应 Theorems 3.1 与 3.2。引言不给 $\mathfrak B$ 的复杂离线构造。

### II. Preliminaries

**A. Problem Formulation**：保留完整的 ODE、wave PDE、两个 Neumann 边界条件。$X(t)\in\mathbb R^n$、$U(t)\in\mathbb R$、$B_i\in L^2(0,D;\mathbb R^n)$ 及可积性条件均**写在行内说明文字中，不另列编号公式**。说明在线只测 $X(t)$、$u(D,t)$、$u_t(D,t)$。不提前给 $\Omega$。

**B. Technical Lemmas**：Lemma 2.1 写一般矩阵对 $(A,B)$ 的参数化 ARE，保留 $\|P_\varepsilon\|\leq K_P\varepsilon$ 和 nilpotent 情况下的 $A^\top P_\varepsilon A\preceq3n^2\varepsilon^2P_\varepsilon$；Lemma 2.2 写 Poincaré 估计。这里的 $B$ 是独立的引理参数，不是第 III 节额外定义的有效向量。

### III. Main Results

#### 开头：像 Xu 的第一段一样，自然说明两个控制器

建议正文直接使用：

> In this section, we propose two low-gain controllers for the system (·)–(·) under two different conditions on the eigenvalues of the ODE matrix $A$.

紧接着先给具体的参数条件，而不是先给 $\mathfrak B$ 后补充参数：

\[
b>0,\qquad 0<q<\frac1D.
\]

接着定义有效输入向量，只引入辅助矩阵 $E$：

\[
\mathfrak B=E^{-1}\int_0^D\cosh(Ax)
 \bigl(AB_2(x)+B_1(x)\bigr)\,\mathrm dx,
\]

\[
E=(bA+qI)\cosh(AD)+A\sinh(AD).
\]

其后**仅写**：

> The matrix $E$ is assumed to be invertible.

正文此处不要讨论 Jordan 块、矩阵函数幂级数、谱映射、矩阵指数，亦不另设关于 $E$ 的 Assumption。对 $A$ 全部特征值在虚轴上时的充分可逆性，在 **Remark 3.4** 给出一次。

#### Assumption 3.1 和第一个控制器

保持 Xu 风格，假设只用一小段：

> **Assumption 3.1:** All eigenvalues of $A$ are on the imaginary axis, and the matrix pair $(A,\mathfrak B)$ is controllable.

先给控制器：

\[
U(t)=-\mathfrak B^\top P_\varepsilon X(t)
-b\,u_t(D,t)-q\,u(D,t),
\]

然后，**如 Xu 的 controller–ARE 顺序**，立即完整展示 ARE，不让读者回翻 Lemma 2.1：

\[
A^\top P_\varepsilon+P_\varepsilon A
-P_\varepsilon\mathfrak B\mathfrak B^\top P_\varepsilon
=-\varepsilon P_\varepsilon.
\]

只写 $P_\varepsilon$ 是上述 ARE 的唯一正定解，尽可能不增添中间定义。

#### 能量、Theorem 3.1 和 Proof

在主定理前定义：

\[
\Omega(t)=|X(t)|^2+
\int_0^D\bigl(u_x^2(x,t)+u_t^2(x,t)\bigr)\,\mathrm dx+u^2(D,t).
\]

紧接着写一句 *The following theorem shows the exponential stability of the closed-loop system.*，再陈述完整的定理。

**Theorem 3.1 的陈述须保留原稿精度**：清楚逐一指出闭环由 ODE (·)、wave PDE (·)、左、右边界条件及 controller (·) 组成；在 Assumption 3.1 下存在 $\varepsilon^*>0$，对每一个 $\varepsilon\in(0,\varepsilon^*]$ 和初值 $(X(0),u(\cdot,0),u_t(\cdot,0))\in\mathbb R^n\times H^1(0,D)\times L^2(0,D)$，唯一解满足 $(X(t),u(\cdot,t),u_t(\cdot,t))\in C([0,\infty);\mathbb R^n\times H^1(0,D)\times L^2(0,D))$。这些成员关系**全部融入定理的句子，不使用 display math，不单独编号**。只有指数衰减结论独立展示：

\[
\Omega(t)\leq\alpha e^{-\beta t}\Omega(0),\qquad t\geq0.
\]

$\alpha,\beta>0$ 不依赖初值。还需保留经典解条件：$u(\cdot,0)\in H^2(0,D)$、$u_t(\cdot,0)\in H^1(0,D)$，以及满足闭环边界的兼容条件。**不能仅用 “a unique energy solution” 替代这段完整陈述。**

**Proof 建议分八个连续段落，但正文无段落小标题。**

**第 1 段：核边值问题，得出 $g(D)=\mathfrak B$，定义 $Z(t)$，求导。** 三条边值关系必须用三个独立的 `equation` 环境，不要并到一个 `aligned`：

\[
g''(x)=A^2g(x)-AB_2(x)-B_1(x),
\]

\[
g'(0)=0,
\]

\[
-g'(D)=(bA+qI)g(D).
\]

紧接着直接展示下面**两条连续的完整推导公式**，不跳过从定积分到 $Eg(D)$ 的中间等号：

\[
\begin{aligned}
&\frac{\mathrm d}{\mathrm dx}
  [\cosh(Ax)g'(x)-A\sinh(Ax)g(x)]\\
&\quad=\cosh(Ax)(g''(x)-A^2g(x))\\
&\quad=-\cosh(Ax)(AB_2(x)+B_1(x)).
\end{aligned}
\]

\[
\begin{aligned}
&\int_0^D\cosh(Ax)(AB_2(x)+B_1(x))\,\mathrm dx\\
&=-[\cosh(Ax)g'(x)-A\sinh(Ax)g(x)]_0^D\\
&=-\cosh(AD)g'(D)+A\sinh(AD)g(D)\\
&=[(bA+qI)\cosh(AD)+A\sinh(AD)]g(D)\\
&=Eg(D).
\end{aligned}
\]

随即独立展示结论，与 Main Results 开头给出的有效输入向量保持完全一致：

\[
\begin{aligned}
g(D)&=E^{-1}\int_0^D\cosh(Ax)(AB_2(x)+B_1(x))\,\mathrm dx\\
&=\mathfrak B.
\end{aligned}
\]

**排版要求**：第一、第二条各自一个 `equation` + `aligned`，保持连等号；第三条结论单独编号。用 *It follows that*, *Integrating over $[0,D]$ gives*, *Thus* 串联，**不拆成只有 $Eg(D)=\cdots$ 的跳步公式**。不插入额外解释，也不在此重复论证 $E$ 的可逆性。然后：

\[
\begin{aligned}
Z(t)={}&X(t)+\int_0^D(Ag(x)-B_2(x))u(x,t)\,\mathrm dx\\
&+\int_0^Dg(x)u_t(x,t)\,\mathrm dx+b\mathfrak B\,u(D,t).
\end{aligned}
\]

$\dot Z(t)$ 推导显示关键积分消项、分部积分、右端边界项以及控制器代入；可删去冗余的中间等号，但保留所有自变量。推得

\[
\dot Z(t)=AZ(t)-\mathfrak B\mathfrak B^\top P_\varepsilon X(t).
\]

**此时才定义**

\[
\varphi(t)=Z(t)-X(t),
\]

并得到

\[
\dot Z(t)=(A-\mathfrak B\mathfrak B^\top P_\varepsilon)Z(t)
+\mathfrak B\mathfrak B^\top P_\varepsilon\varphi(t).
\]

**第 2 段：wave 能量和变换误差。**

\[
E_u(t)=\frac12\int_0^D\bigl(u_x^2(x,t)+u_t^2(x,t)\bigr)\,\mathrm dx
+\frac q2u^2(D,t).
\]

\[
\dot E_u(t)=-b\,u_t^2(D,t)
-u_t(D,t)\mathfrak B^\top P_\varepsilon X(t).
\]

紧随其后由 Poincaré、Cauchy–Schwarz 逐项证明 $|\varphi(t)|^2\leq C_\varphi E_u(t)$，并单独展示

\[
C_\varphi=12D(D+q^{-1})\|Ag-B_2\|_2^2
+6\|g\|_2^2+\frac{6b^2|\mathfrak B|^2}{q}.
\]

**第 3 段：有限维 Riccati 能量。**

\[
V_1(t)=Z^\top(t)P_\varepsilon Z(t).
\]

借助刚刚展示过的主结果 ARE，完整推导至

\[
\dot V_1(t)\leq-\varepsilon V_1(t)
+C_\varphi|\mathfrak B|^2\|P_\varepsilon\|^2E_u(t).
\]

**第 4 段：wave multiplier。** 单独定义 $\kappa=(1+Dq)/4$，再定义

\[
\begin{aligned}
V_2(t)={}&E_u(t)+\delta\int_0^D x\,u_x(x,t)u_t(x,t)\,\mathrm dx\\
&+\delta\kappa\int_0^D u(x,t)u_t(x,t)\,\mathrm dx.
\end{aligned}
\]

证明 $E_u(t)/2\leq V_2(t)\leq3E_u(t)/2$。**两个 multiplier 求导各用一个编号公式**，再写完整 $\dot V_2(t)$。三项混合积的 Young 估计建议各用一个编号公式，三个 $\delta$ 充分条件也分别编号。最后得到

\[
\dot V_2(t)\leq-cE_u(t)
+C\bigl(\mathfrak B^\top P_\varepsilon X(t)\bigr)^2.
\]

$c,C$ 的取值各在单独的公式中展示，而不是只说 “for some positive constants”。

**第 5 段：利用 $X(t)=Z(t)-\varphi(t)$ 消掉反馈交叉项。**

\[
\begin{aligned}
(\mathfrak B^\top P_\varepsilon X(t))^2
&\leq2|\mathfrak B|^2\|P_\varepsilon\|V_1(t)\\
&\quad+2C_\varphi|\mathfrak B|^2\|P_\varepsilon\|^2E_u(t).
\end{aligned}
\]

立刻并入 $\dot V_2(t)$ 估计，避免重定义常数。

**第 6 段：最后定义总 Lyapunov 泛函。**

\[
V(t)=V_1(t)+\varepsilon V_2(t).
\]

组合 $\dot V_1(t)$ 与 $\dot V_2(t)$，最后才利用 $\|P_\varepsilon\|\leq K_P\varepsilon$，得到

\[
\dot V(t)\leq-\frac\varepsilon2V_1(t)-\frac{\varepsilon c}3V_2(t)
\leq-\beta V(t).
\]

$\beta=\min\{\varepsilon/2,c/3\}$ 独立成式。

**第 7 段：$V(t)$ 与 $\Omega(t)$ 的等价性，以及显式常数。** 先写 $|X(t)|^2\leq2|Z(t)|^2+2C_\varphi E_u(t)$、反向不等式和 $P_\varepsilon$ 的上下界，得到

\[
m_\varepsilon\Omega(t)\leq V(t)\leq M_\varepsilon\Omega(t).
\]

随后紧邻给出**能够直接从已证不等式核验的具体常数**，不再只写 “for some $m_\varepsilon,M_\varepsilon>0$”：

\[
m_\varepsilon=\min\left\{\frac{\lambda_{\min}(P_\varepsilon)}2,
\frac{\varepsilon}{4(C_\varphi+\max\{1,q^{-1}\})}\right\},
\]

\[
M_\varepsilon=\max\left\{2\|P_\varepsilon\|,
\max\{1,q\}\left(C_\varphi\|P_\varepsilon\|+\frac{3\varepsilon}{4}\right)\right\}.
\]

这些常数来自 $E_u(t)\geq\frac12\min\{1,q\}(\|u_x(\cdot,t)\|_2^2+\|u_t(\cdot,t)\|_2^2+u^2(D,t))$、$E_u(t)\leq\frac12\max\{1,q\}(\|u_x(\cdot,t)\|_2^2+\|u_t(\cdot,t)\|_2^2+u^2(D,t))$，以及 $V_2(t)$ 的 $\frac12 E_u(t)$、$\frac32 E_u(t)$ 上下界。在正式正文中只需显示 $m_\varepsilon,M_\varepsilon$，无需把这些辅助界再铺成新一串编号公式。固定 $\varepsilon>0$ 时两个常数均严格正。

**第 8 段：只写最终指数估计并结束 Proof。** *Consequently* 后直接给

\[
\Omega(t)\leq\frac{M_\varepsilon}{m_\varepsilon}
 e^{-\beta t}\Omega(0),\qquad t\geq0.
\]

不再逐行展示 $m_\varepsilon^{-1}V(t)$、$e^{-\beta t}V(0)$ 等已可立即推出的中间项。
随后简短指向 Appendix 的良定性证明；不要在证毕之前插入长段解释。

#### 第二个控制器：自然过渡，不新建小节

过渡句参考 Xu：

> Next, we consider the case where all eigenvalues of $A$ are at the origin. In this case, a simpler effective input vector can be obtained.

先给

\[
\bar B=\frac1q\int_0^D B_1(x)\,\mathrm dx.
\]

接着 **Assumption 3.2**：全部特征值在原点、$(A,\bar B)$ controllable。然后第二个 controller 和对应 ARE 各自独立编号：

\[
U(t)=-\bar B^\top P_\varepsilon X(t)-b\,u_t(D,t)-q\,u(D,t),
\]

\[
A^\top P_\varepsilon+P_\varepsilon A
-P_\varepsilon\bar B\bar B^\top P_\varepsilon=-\varepsilon P_\varepsilon.
\]

**Theorem 3.2** 同样明确所考虑的 ODE、wave、两条边界条件、控制器、初值空间、唯一解正则性、指数估计及经典解条件；直接复用 $\Omega$ 定义。初值空间和解空间仍作为**行内公式**，不要回到 v6 的独立公式写法。

**Proof**：特例核方程及两个边界条件用三个独立 `equation`：

\[
g''(x)=-B_1(x),\qquad g'(0)=0,\qquad-g'(D)=qg(D).
\]

在 LaTeX 排版中这三个要拆开，不能合并成同一个带编号的式子。得到 $g(D)=\bar B$，定义另一套 $Z(t)$；求导带全 $X(t),Z(t),u(D,t),u_t(D,t)$；定义 $\varphi(t)$ 和 $V_1(t)$，利用 nilpotent Riccati 引理控制附加的 $Z^\top(t)P_\varepsilon A\varphi(t)$。然后复用 wave multiplier 的 $V_2(t)$ 估计，单独定义 $V(t)$，证明严格衰减。**不把两种补偿核视为相同。**

#### Proof 结束后四条 Remarks

- **Remark 3.1**：在线只需 $X(t),u(D,t),u_t(D,t)$，核函数与积分型 $Z(t)$ 仅供证明和离线分析。
- **Remark 3.2**：$\|P_\varepsilon\|=O(\varepsilon)$；低增益与慢收敛的取舍；充分小增益界通常保守。
- **Remark 3.3**：Assumption 3.1 不要求 $A$ 可对角化，覆盖振荡器与积分器链；相较于允许不稳定 ODE 的 full-state forwarding/backstepping 方法，代价是谱约束。
- **Remark 3.4**：这里才说明 $\cosh(Ax)$、$\sinh(Ax)$ 的标准矩阵函数意义，并验证对 $A$ 的任意虚轴特征值 $\lambda=i\omega$，有

\[
(b\lambda+q)\cosh(\lambda D)+\lambda\sinh(\lambda D)
=(q+ib\omega)\cos(\omega D)-\omega\sin(\omega D)\ne0,
\]

从而所选参数和 Assumption 3.1 实际上足以确保 $E$ 可逆。**该证明只写在 Remark，不重复到 Main Results 开头或 Proof 的核步骤。**

### IV. Numerical Simulation

保留三组原有示例及数值，Example 1 的一般有效输入向量和控制增益统一写 $\mathfrak B$；Example 2 写 $\bar B$；Example 3 固定 plant 对比不同 $\varepsilon$，其中 $\varepsilon$ 的几个测试取值为**行内集合**，不再单独编号。沿用既有图像文件名 `fig2.png`、`fig3.png`、`fig4.png`；未提供图片时显示占位，不能编造结果或图片。

### V. Conclusion

保持单段，顺序为：研究问题 → 两个有限维低增益控制器 → 指数稳定性 → nilpotent 公式 → 仅需边界与 ODE 测量 → 谱限制及未来工作。

### Appendix. Well-Posedness

保留闭环半群生成证明，尤其是生成元定义域、边界兼容条件和由经典解延拓到能量解的论证；与两条主定理的正则性陈述一致。

## 三、改稿核对表

- [x] III 节第一句话说明提出两个控制器。
- [x] 在 $\mathfrak B$ 前给出 $b>0,0<q<1/D$。
- [x] $E$ 后仅给一句 assumed invertible，充分条件放 Remark 3.4。
- [x] Assumptions / Theorems / Remarks 显示为 `3.1` / `3.2`，不是 `III.1`。
- [x] Assumption 使用“all eigenvalues of $A$”，不用 spectrum 符号写法。
- [x] 一般结果始终使用 $\mathfrak B$，controller 后再次独立展示 ARE。
- [x] Theorem 3.1 保留完整解空间、初值空间、指数估计和经典解条件；所属空间在文字中使用行内公式。
- [x] 核微分方程和两个边界条件分别编号，第二个控制器也一样。
- [x] $\varphi, E_u, V_1, V_2, V$ 等定义单独展示。
- [x] 证明中的 $u(x,t),u_t(x,t),u_x(x,t),X(t),Z(t)$ 不省略变量。
- [x] 控制器和其他短公式优先合并成一行，长式仅在双栏宽度需要时换行。
- [x] 变量空间、核函数可积条件与 Example 3 的参数集合在正文行内表达。
- [x] 不用框住公式；不在 III 节另设小节标题。
- [x] 引言控制器以及 Example 1 的 $\mathfrak B$ 记号与主结果一致。
