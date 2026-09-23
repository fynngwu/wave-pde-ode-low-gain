最终目标设为 8 页
---






### 预计节省

**0.2~0.3 页**。

---

# 五、第三刀：Lemma 1 也不要把所有低增益性质全部写出来

你现在：

```latex
\begin{enumerate}
\item ...
\item ...
\item ...
\end{enumerate}
```

尤其是：

```latex
\frac{dP_\varepsilon}{d\varepsilon}>0
```

对后面的主证明其实完全没有作用。

你真正使用的核心是：

$$
P_\varepsilon>0,\qquad
\|P_\varepsilon\|=O(\varepsilon),
$$

以及 nilpotent 情况：

$$
A^\top P_\varepsilon A
\preceq C\varepsilon^2P_\varepsilon.
$$

所以可以压缩成：

```latex
We use the standard low-gain Riccati properties:
$P_\varepsilon>0$,
$\|P_\varepsilon\|=O(\varepsilon)$,
and, in the nilpotent case,
$A^\top P_\varepsilon A\preceq C_A\varepsilon^2P_\varepsilon$.
```

然后 cite。

这个尤其重要，因为你 Introduction 里已经引用了 low-gain 文献，没有必要在 conference paper 中再复述一整套定理。

### 预计节省

**0.2 页左右**。

---

# 六、第四刀：Forwarding transformation 部分可以砍掉约一半

你现在这一段：

```latex
Differentiating...
...
Integrating the wave term by parts twice yields...
...
Using the kernel identities...
...
Under the feedback law...
...
Since $X=Z-\varphi$...
```

其实逻辑非常好，但是太“逐行证明”。

会议论文可以只保留两个关键结果：

### 第一块

定义 \(Z\)：

$$
Z=X+\cdots
$$

### 第二块

直接给出：

$$
\dot Z
=
AZ+B(u_x(D)+bu_t(D)+qu(D)).
$$

然后一句：

> Direct differentiation, integration by parts, and the kernel identities yield (X).

没必要把 integration by parts 的完整过程写出来。

---

## 你可以删掉这个公式

```latex
\begin{equation}
\begin{aligned}
\dint{g(x)u_{xx}(x,t)\,\mathrm{d}x}
={}&...
\end{aligned}
\end{equation}
```

因为这是标准计算。

### 预计节省

**0.2~0.4 页**。

---

# 七、第五刀：Kernel offline computation 现在太详细

你当前写了：

$$
\zeta(x),\quad
\zeta'(x),\quad
\Xi,\quad
E,
$$

然后 Proposition + proof。

实际上核心贡献只是：

> **distributed kernels \(\to\) offline BVP \(\to\) effective input vector \(B=g(D)\)**

这是文章真正值得读者记住的东西。

可以压缩成一个 proposition：

$$
B = \mathcal{K}(A,B_1,B_2,b,q).
$$

甚至把 \(E,\Xi\) 保留，但删掉 Proposition 的完整证明。

直接：

> By variation of constants, the endpoint value \(B=g(D)\) is obtained from (X).

然后把 proof 变成一句。

---

# 八、Assumption 1 和 Assumption 3 有些重复，可以合并

你现在有：

### Assumption 1

* \(b>0\)
* \(0<q<1/D\)
* \(E\) invertible
* controllable
* purely imaginary spectrum

然后：

### Assumption 3

* 前三项又说一次
* \(\sigma(A)\subset\{\Re s\le0\}\)
* critical pair controllable

这种在会议论文中会让人感觉“假设系统在不停升级”。

更好的组织方式：

## Assumption 1

统一定义：

$$
b>0,\quad 0<q<1/D,\quad E \text{ invertible}.
$$

## Assumption 2

Marginal case：

$$
\sigma(A)\subset i\mathbb R,\quad (A,B)\text{ controllable}.
$$

## Assumption 3

Nonexpansive extension：

$$
\sigma(A)\subset\{\Re s\le0\},
$$

and only critical pair controllable.

这样读者会非常容易理解：

> **common structural assumptions + spectral assumptions**

---

# 九、Theorem 3 我建议大幅压缩

你现在的 Nonexpansive ODE Extension 本质是：

> Theorem 1 的 proof 经过 spectral decomposition 后重复一次。

而会议论文最忌讳：

> 同一个证明讲两遍。

你现在这部分大约能砍掉 **40%~60%**。

可以直接写：

> For the nonexpansive case, decompose \(A\) into its Hurwitz and critical spectral subspaces as in (XX). Applying the same forwarding transformation to the critical subsystem and a standard quadratic Lyapunov function to the Hurwitz subsystem yields ...

然后只保留：

$$
\dot V_c\le-\varepsilon V_c+C_c\|P_c\|^2E_u
$$

$$
\dot V_s\le-\lambda_sV_s+C_s\|P_c\|^2(|Z_c|^2+E_u)
$$

最后一句：

$$
V=V_c+V_s+\varepsilon V_2
$$

satisfies

$$
\dot V\le-\beta V.
$$

没必要再写一遍“norm-equivalence and well-posedness arguments are the same...”。

一句话解决：

> The norm equivalence and well-posedness follow exactly as in Theorem 1 and are omitted.

### 预计节省

**0.4~0.6 页**。

---

# 十、Nilpotent case 不需要独立写这么长

这其实是你论文里一个很漂亮的 special case，但现在写得有些过度。

你真正需要保留的是：

$$
g''+B_1=0
$$

$$
\bar B
=
\frac1q\int_0^D B_1(x)\,dx.
$$

这是亮点。

然后直接说：

> The closed-loop stability proof is identical to Theorem 1 except for the additional \(A^\top P_\varepsilon A\) term, which is \(O(\varepsilon^2)\) by the nilpotent Riccati estimate.

这样即可。

你现在这一段：

```latex
Using ... gives
...
Define ...
...
The wave functional ...
...
The same calculation gives ...
...
norm-equivalence...
```

大量重复 Theorem 1。

### 最好做法

保留：

1. explicit \(\bar B\)
2. controller
3. theorem
4. 一段 8~10 行的 proof sketch

就够了。

### 预计节省

**0.5 页左右**。

---

# 十一、Introduction 可以削掉至少 30%

你的 Introduction 本身写得不错，但对于 8 页会议稿略重。

目前最大问题是：

> literature review 太详细，贡献又重复叙述。

例如下面几处信息高度重叠：

> backstepping controllers tend to be algebraically involved...

随后又：

> Compared with the backstepping designs...

然后 contribution 2 再说一次：

> proposed controller is structurally simpler...

最后 remark 又说一次。

实际上可以只保留一次。

---

## 我建议 Introduction 只保留四段逻辑

### Paragraph 1

应用背景：

wave actuator + distributed coupling。

### Paragraph 2

现有方法：

backstepping exists but online complexity high；
low-gain works for parabolic case；
wave case remains open。

### Paragraph 3

本文方法：

forwarding + boundary damping + low gain。

### Paragraph 4

贡献：

只保留 **3 点**。

不要 5 点，甚至不要太长。

---


# 十三、Simulation 是你第二个最大的“浪费页数”来源

你前面说的判断其实是对的：

> **不要再专门做 \(\varepsilon=0.2,0.5,0.9\) 这种性能参数扫描作为主实验。**

这个 Example 3 很容易让审稿人问：

> Why is \(\varepsilon=0.5\) unstable if the theorem says sufficiently small \(\varepsilon\)?

这反而把注意力引到了 theorem 的 conservative bound 上。

而且你自己之前也已经发现：

> “50s 后指标”很像为了凑实验。

这里完全可以删。

---

# 十四、建议整个 Simulation 重构成“两组主实验”

我更推荐：

## Example 1：General Marginally Stable Case

目的：

> 证明 general controller 有效。

保留：

* \(A\)
* \(B_1,B_2\)
* computed \(B\)
* controller
* ODE trajectory
* PDE surface

### 一张双子图：

**(a)** \(X_1,X_2\)

**(b)** \(u(x,t)\)

足够。

---

## Example 2：Nilpotent Case + Baseline Comparison

目的：

> 证明 nilpotent controller 有效，并且说明 boundary damping alone 不够。

这一组实际上比 Example 3 更有价值。

保留：

### (a) baseline

boundary damping only

### (b) proposed

full controller

可以做成：

$$
2\times2
$$

比如：

* (a) ODE baseline
* (b) wave baseline
* (c) ODE proposed
* (d) wave proposed

这张图的信息密度很高。

---

# 十五、把 Example 3 删除

我会非常明确地建议：

> **直接删除整个 Example 3。**

不是移动，而是删。

因为它属于：

> supplementary validation

不是：

> core evidence

在 8 页限制下，完全不划算。

你可以把一句话放进 Example 1：

> In simulations, sufficiently small values of \(\varepsilon\) are used, consistent with the small-gain condition.

结束。

---

# 十六、Simulation 的方法部分也应该压缩

你现在这一大段：

```latex
Introducing v...
...
N = 200...
...
ghost points...
...
RK4...
...
fine quadrature...
...
result is insensitive...
```

太像 numerical-methods paper。

CAC 会议完全可以：

> The wave PDE is discretized by a uniform method-of-lines scheme with \(N=200\) spatial intervals and RK4 time integration. The Neumann boundary conditions are implemented using ghost points, while distributed integrals are evaluated by the trapezoidal rule. Unless otherwise stated, \(D=1\) and \(\Delta t=5\times10^{-4}\).

就够了。

---

# 十七、不要把所有数值计算结果都写成公式

这是一个很容易节省篇幅的地方。

例如：

```latex
P_{0.12}=
\begin{pmatrix}
...
\end{pmatrix}
```

然后：

```latex
B^\top P_{0.12}
=
...
```

然后：

```latex
U(t)=...
```

这些全写一遍。

对于会议论文，可以只保留最终 controller：

$$
U(t)
=
-0.1786x_1+0.3689x_2-u_t(D,t)-0.4u(D,t).
$$

\(P\) 不一定要展示。

读者真正关心的是：

> controller 被实现了什么。

而不是矩阵内部元素。

---



