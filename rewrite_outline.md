wave_ode_low_gain.tex 技术文献改写大纲（修订审核版）
0. 改写定位与边界
0.1 目标
将现有稿件改写为一篇具有 IEEE TAC 论文论证节奏的技术文献，重点迁移参考论文 Low-Gain Compensation for PDE-ODE Cascade Systems With Distributed Diffusion and Counter Convection 第 2、3、4 章的组织范式，不复制其公式或具体结论。改写对象为有限维 ODE 与一维波动 PDE 构成的串联耦合系统，核心方法保留离线核函数构造、辅助补偿变量设计和低增益边界反馈控制。
0.2 保留的研究对象
研究模型保持为标准 PDE-ODE 串联耦合形式，修正积分项书写格式，统一领域常规函数空间表述：
\[
    \dot X(t)=AX(t)+\int_0^D B_1(x)u(x,t)\,dx
     +\int_0^D B_2(x)u_t(x,t)\,dx,
    \]
\[
    u_{tt}(x,t)=u_{xx}(x,t),\quad x\in(0,D),\qquad
    u_x(0,t)=0,\quad u_x(D,t)=U(t).
    \]
其中位移耦合项系数 $B_1(x)\in L^2(0,D;\mathbb R^n)$，速度耦合项系数 $B_2(x)\in L^2(0,D;\mathbb R^n)$，全程不引入 Radon 测度概念，仅在文末备注可拓展至积分有界的广义耦合场景。现阶段不改写仿真代码、不生成新图、不处理图文件引用；第 4 章仅搭建规范化算例写作框架。
0.3 拟采用的整体章节
1. Section II — Preliminaries
- 2.1 Problem Formulation
- 2.2 Technical Lemmas
- 2.3 Offline Kernel Construction
2. Section III — Main Results
- 3.1 General Marginally Stable Case
- 3.2 Nilpotent Case
- 3.3 Boundary of the Result and Implementation Remarks

0.4 符号纪律与排版要求
后续证明必须优先复用 Preliminaries 中已经定义的原始符号
$X,Z,u,u_t,g,g_0,\mathcal B,P_\varepsilon,E_u,\Omega$ 及其已有下标，
不得为了缩短公式随意增加新的状态、标量或中间变量，尤其不得引入
$s,v,p$ 作为证明中的临时记号。若原始符号导致公式过长，应采用
`align`、`aligned` 等环境换行，不得通过新增变量规避排版问题。每次改写后
必须编译两遍，并检查日志中不存在 `Overfull`、`Underfull`、未定义引用或
重复定义警告。

1. Section II — Preliminaries
1.1 2.1 Problem Formulation
1.1.1 研究对象与状态方程
定义空间区间 $x\in[0,D]$、时间 $t\ge 0$，完整列出有限维 ODE 状态方程、一维波动 PDE 方程及两端 Neumann 边界条件。给出系统初值条件：
\[
    X(0)=X_0,\qquad u(\cdot,0)=u_0,\qquad u_t(\cdot,0)=u_1,
    \]
明确系统经典解、能量解对应的状态空间，以及初值需满足的边界相容性条件。全文统一使用 $X$ 表示 ODE 原始状态，与后续构造的辅助状态 $Z$ 严格区分，避免符号混淆。
1.1.2 符号和物理含义逐项定义
仅定义系统固有状态与模型参数，删除控制器参数 $b,q,\varepsilon$ 的前置解释（统一后置控制器章节定义），逐项说明核心符号：
- $X(t)\in\mathbb R^n$：有限维被控对象状态；
- $A\in\mathbb R^{n\times n}$：开环 ODE 系统矩阵，决定系统中性、幂零模态特性；
- $u(x,t)$：波动 PDE 状态，表征柔性介质位移场；
- $U(t)$：PDE 右端边界控制输入；
- $D>0$：波动介质空间长度；
- $B_1(x)\in L^2(0,D;\mathbb R^n)$：位移分布式耦合系数；
- $B_2(x)\in L^2(0,D;\mathbb R^n)$：速度分布式耦合系数。
定义系统能量函数与整体状态度量函数，明确控制目标：实现闭环系统全局指数稳定，即存在常数 $\alpha,\beta>0$，使得 $\Omega(t)\le \alpha e^{-\beta t}\Omega(0)$ 成立。
\[
    E_u(t)=\frac12\int_0^D(u_x^2+u_t^2)\,dx+\frac q2u^2(D,t),
    \]
\[
    \Omega(t)=|X(t)|^2+\int_0^D(u_x^2+u_t^2)\,dx+u^2(D,t).
    \]
1.1.3 工程背景与结构解释
本模型可描述有限维控制器/机械部件与柔性波动介质、传输线或声学管道之间的耦合动力学行为。一维波动方程为无固有耗散的保守系统，常规无约束松弛边界无法实现能量耗散，本文采用耦合边界条件构建闭环控制框架。该分布式耦合模型具备通用性，可兼容多种典型耦合场景，适配不同工程柔性传动系统建模需求。
1.1.4 广义拓展说明（替代原退化特例）
本文耦合系数均采用 $L^2$ 空间标准定义，模型具备良好的广义性。在后续备注中可拓展说明：本文理论框架可自然兼容积分有界的广义位移耦合场景（含单点耦合近似形式），无需额外修改证明体系，全程规避测度相关表述。删除原冗余退化特例分析，相关场景统一在仿真算例中验证。
1.2 2.2 Technical Lemmas
仅保留核心低增益 Riccati 引理，删除所有其他辅助引理，仅陈述稳定性证明必需的核心工具，不展开证明过程，明确引理适用条件、数学表达式及理论用途。
Lemma 1 — Low-gain Riccati property
在有限维系统 $(A,B)$ 可控、且矩阵 $A$ 谱集全部位于虚轴的前提下，定义低增益 Riccati 方程：
\[
    A^\top P_\varepsilon+P_\varepsilon A
     -P_\varepsilon BB^\top P_\varepsilon=-\varepsilon P_\varepsilon.
    \]
明确核心性质：方程存在唯一正定解 $P_\varepsilon>0$，且满足阶次估计 $P_\varepsilon=O(\varepsilon)$；针对幂零系统，额外满足约束 $A^\top P_\varepsilon A\preceq C_A\varepsilon^2P_\varepsilon$。本引理为有限维系统低增益稳定设计、PDE-ODE 交叉扰动项吸收的核心理论依据。
1.3 2.3 Offline Kernel Construction
精简章节标题，重构辅助变量设计核心思路，删除 Artstein 命名、拆分模块与幂零前置内容，统一整合核函数求解体系。
1.3.1 辅助系统设计核心思路
针对 PDE-ODE 强耦合带来的稳定性分析难题，本文构建辅助状态子空间，核心设计逻辑：波动 PDE 系统可通过边界反馈实现自身能量耗散稳定，在此基础上构造耦合辅助状态 $Z$，实现 PDE 与 ODE 动力学解耦。只需保证辅助状态子空间可实现自稳定，即可推导出原始有限维 ODE 状态的稳定性，彻底规避直接处理耦合交叉项的难点。
定义全局辅助状态变量：
\[
    Z=X+\int_0^Dg_0(x)u(x,t)\,dx
     +\int_0^Dg(x)u_t(x,t)\,dx+b\,g(D)u(D,t),
    \]
1.3.2 核函数控制方程与求解条件
为实现耦合动力学完全解耦，推导核函数 $g(x),g_0(x)$ 满足的边值方程组：
\[
    g''+B_1=A(Ag-B_2),\qquad g_0=Ag-B_2,
    \]
边界约束条件：
\[
    g'(0)=0,\qquad -g'(D)=(bA+qI)g(D).
    \]
将核函数方程转化为一阶线性状态系统 $Y'=\mathcal A_gY+\mathcal F(B_1,B_2)$（其中 $Y=(g,g')^\top$），通过状态转移矩阵推导核函数可解的充要条件：系统对应系数矩阵可逆，保证核函数存在唯一有界解。定义边界有效耦合向量 $\mathcal B=g(D)$，该向量可通过离线求解得到，无需在线迭代计算，大幅降低控制器实现复杂度。
Remark：核函数的唯一可解性是本文控制方案成立的前置基础假设，保证辅助状态解耦设计的有效性，是后续闭环稳定性分析的前提条件。
1. Section III — Main Results
2.1 3.1 General Marginally Stable Case
严格遵循 Assumption → ARE → Controller → Theorem → Proof → Remark 标准化论证框架，证明过程补充所有系数 $C_1,C_2,C_3,C_4$ 显式推导步骤、参数阶次分析、Lyapunov 泛函与状态度量 $\Omega(t)$ 的显式上下界约束，仅保留 Remark1、Remark4。
2.1.1 Assumption 1
逐条列明系统建模、设计参数、可控性与正则性假设，分类明确模型固有条件与设计条件：
1. 控制参数满足 $b>0$ 且 $0<q<1/D$；
2. 核函数求解系数矩阵可逆，保证辅助系统唯一可解；
3. 有效耦合通道 $(A,\mathcal B)$ 满足可控性条件；
4. 开环 ODE 矩阵 $A$ 的全部特征值位于虚轴，系统为边际稳定；
5. 系统初值满足 $X_0\in\mathbb R^n$，$u_0\in H^1(0,D)$，$u_1\in L^2(0,D)$，且满足边界相容性条件。
2.1.2 ARE 与控制器
基于上述假设，引入低增益 Riccati 方程正定解 $P_\varepsilon$：
\[
    A^\top P_\varepsilon+P_\varepsilon A
     -P_\varepsilon\mathcal B\mathcal B^\top P_\varepsilon
     =-\varepsilon P_\varepsilon.
    \]
设计 PDE 右端边界低增益反馈控制律：
\[
    U(t)=-bu_t(D,t)-qu(D,t)-\mathcal B^\top P_\varepsilon X(t).
    \]
控制律三项功能明确：边界速度阻尼项耗散波动系统能量、边界位移刚度项抑制波动常值模态、低增益状态反馈项镇定有限维边际稳定模态。
2.1.3 Theorem 1
存在常数 $\varepsilon^*>0$，当低增益参数满足 $0<\varepsilon\le\varepsilon^*$ 且系统初值满足上述正则性条件时，闭环 PDE-ODE 耦合系统存在唯一经典解，且存在正实数 $\alpha,\beta$，使得系统整体状态满足全局指数衰减估计：
\[
    \Omega(t)\le\alpha e^{-\beta t}\Omega(0),\qquad t\ge0.
    \]
明确 $\varepsilon^*$ 为交叉扰动项吸收对应的理论上界，仅为充分条件而非必要条件。
2.1.4 Proof of Theorem 1（完整系数显式推导版）
严格沿用五步论证框架，补充所有关键系数显式推导、阶次分析、泛函等价性证明，补齐 $V\sim\Omega$ 的上下界显式常数。
Step 1 — 辅助变换和动力学重写：定义辅助偏差变量 $\varphi=Z-X$，对辅助状态 $Z$ 逐阶求导，结合波动方程动力学、Neumann 边界条件、核函数边值方程，逐项消去位移与速度耦合项 $B_1,B_2$，严格推导解耦后辅助系统动力学方程：$\dot Z=AZ+\mathcal B\bigl(u_x(D)+bu_t(D)+qu(D)\bigr)$。代入边界控制律，替换原始状态 $X=Z-\varphi$，得到含波动能量扰动的低增益闭环有限维动力学模型，明确所有耦合余项表达式。
Step 2 — 构造多分量 Lyapunov 泛函：有限维子系统泛函 $V_1=Z^\top P_\varepsilon Z$；波动子系统构造带交叉乘子的 Lyapunov 泛函：
\[
    V_2=E_u+\delta_1\int_0^Dxu_xu_t\,dx
     +\delta_2\int_0^Duu_t\,dx,
    \]
固定参数关系 $\delta_2=\kappa\delta_1$、$Dq/2<\kappa<1/2$，严格推导 $V_2$ 与原始波动能量 $E_u$ 的等价性常数（显式上下界）。定义系统总 Lyapunov 泛函 $V=V_1+\varepsilon V_2$，完成泛函正定性证明。
Step 3 — 轨迹求导与逐项放缩（核心系数显式推导）：沿闭环轨迹对 $V_1,V_2$ 分别求导，结合低增益 Riccati 性质、$L^2$ 空间不等式，逐项计算导数项。精准推导所有放缩系数 $C_1,C_2,C_3,C_4$ 的显式表达式，明确系数与系统参数 $D,b,q$、耦合范数、矩阵维数的对应关系，最终得到：
\[
    \dot V_1\le -\varepsilon V_1+C_1|P_\varepsilon|^2E_u,
    \]
\[
    \dot V_2\le -c_3E_u+C_3|P_\varepsilon|V_1
     +C_4|P_\varepsilon|^2E_u.
    \]
Step 4 — 小参数阶次分析：基于 $P_\varepsilon=O(\varepsilon)$ 的核心性质，精准标注所有扰动项阶次：$|P_\varepsilon|^2E_u=O(\varepsilon^2)E_u$、$\varepsilon|P_\varepsilon|V_1=O(\varepsilon^2)V_1$。通过显式系数不等式推导，确定低增益上界 $\varepsilon^*$ 的取值依据，保证高阶小扰动项可被系统负定耗散项完全吸收，给出 $\varepsilon^*$ 的显式约束不等式。
Step 5 — 指数稳定与适定性证明（补齐 $\boldsymbol{V\sim\Omega}$ 显式界）：结合前述放缩结果，推导总泛函导数负定约束 $\dot V\le -\beta_1V_1-\varepsilon\beta_2E_u\le-\beta V$。严格推导 Lyapunov 泛函 $V$ 与系统状态度量 $\Omega(t)$ 的等价关系，给出显式常数 $m,M>0$，满足 $m\Omega(t)\le V(t)\le M\Omega(t)$。结合 Gronwall 不等式完成全局指数稳定证明，最后通过半群理论验证闭环系统解的存在唯一性。
2.1.5 Remark 1、4（保留指定内容，精简优化）
Remark 1（与 backstepping 方法结构对比）：相较于传统 PDE-ODE 串联系统的 backstepping 控制方案，本文所提低增益边界反馈策略无需全域 PDE 状态信息与在线分布式核运算，仅需有限维 ODE 状态、边界位移与速度信号，大幅降低在线硬件实现复杂度。对应的代价是本文方法对有限维系统矩阵施加虚轴谱约束与有效通道可控性假设，该类约束在部分通用 backstepping 设计中无需满足。
Remark 4（方法适用边界与退化场景）：本文定理可自然覆盖无速度耦合、纯分布位移耦合等退化场景；但方法存在明确适用边界，对于开环 ODE 矩阵含右半平面特征值、有效耦合通道不可控的系统，本文低增益控制方案无法实现闭环稳定，需改用高增益控制或 backstepping 设计方案。同时本文 $L^2$ 耦合框架可拓展至积分有界的广义耦合场景，适配工程常见的单点局部耦合工况。
2.2 3.2 Nilpotent Case
沿用通用情形标准化论证框架，删除原前置核函数简化推导，统一在本章节完成幂零系统专属推导，补充系数显式放缩与阶次分析。
2.2.1 Assumption 2
1. 开环 ODE 矩阵 $A$ 为幂零矩阵，全部特征值位于原点；
2. 简化有效耦合通道 $(A,\bar{\mathcal B})$ 满足可控性条件，其中 $\bar{\mathcal B}=q^{-1}B_1([0,D])$；
3. 控制参数与系统初值满足与一般情形一致的正则性、边界相容性条件。
2.2.2 简化 ARE 与控制器
针对幂零系统，低增益 Riccati 方程形式不变，仅有效耦合向量替换为显式简化形式：
\[
    A^\top P_\varepsilon+P_\varepsilon A
     -P_\varepsilon\bar{\mathcal B}\bar{\mathcal B}^\top P_\varepsilon
     =-\varepsilon P_\varepsilon,
    \]
对应边界控制律：
\[
    U=-bu_t(D)-qu(D)-\bar{\mathcal B}^\top P_\varepsilon X.
    \]
2.2.3 Theorem 2
存在 $\varepsilon^*>0$，当 $0<\varepsilon\le\varepsilon^*$ 时，幂零型 PDE-ODE 耦合闭环系统解存在唯一，且系统状态满足全局指数衰减估计 $\Omega(t)\le\alpha e^{-\beta t}\Omega(0)$。
2.2.4 Proof of Theorem 2
沿用一般情形五步证明框架，保留所有通用放缩与泛函推导，仅补充幂零系统专属核心推导：1. 代入显式有效耦合向量 $\bar{\mathcal B}$ 完成动力学解耦；2. 利用幂零矩阵专属 Riccati 阶次估计 $A^\top P_\varepsilon A\preceq C_A\varepsilon^2P_\varepsilon$，精准控制耦合高阶扰动项；3. 补齐幂零情形下所有系数显式表达式与 $V\sim\Omega$ 等价常数，其余波动系统稳定性、适定性论证与一般情形完全一致。
2.3 3.3 Theoretical Boundary and Implementation Remarks
集中梳理全文理论约束与工程实现边界，规避重复论述：
1. 参数条件 $q<1/D$ 为本文 Lyapunov 乘子设计对应的充分稳定条件，非闭环系统稳定的必要条件；
2. 低增益上界 $\varepsilon^*$ 由 Riccati 矩阵阶次与扰动吸收系数共同决定，仅为理论保障边界，不等同于工程最大可用增益；
3. 控制器在线实现仅需有限维状态、边界位移与速度信号，无需全域核函数与 PDE 状态信息，离线仅需完成核函数与有效耦合向量求解；
4. 本文全部理论基于 $L^2$ 空间耦合假设，可兼容积分有界广义耦合场景，但无法拓展至无界速度耦合工况；
5. 方法核心优势为在线结构简洁、离线在线分工清晰，固有局限为存在系统谱约束与边界测量要求。

证明内容双版本拆分说明（详细手推MD + 期刊精简TeX）
一、整体拆分规则（统一执行标准）
为同时满足可复现手推学习与IEEE TAC期刊投稿规范，将第三章所有稳定性证明内容拆分为两套独立内容，严格区分用途、粒度和详略：
- 版本1：详细手推 MD 文件（备查底稿）：零跳步、全展开、每一步放缩/求导/不等式引用均标注依据，所有系数 \(C_1,C_2,C_3,C_4,\beta,m,M\) 逐行推导、显式写出，适合个人复盘、审稿人追问、复现推导、学生研读。
- 版本2：正文 TeX 精简版（投稿正文）：删除琐碎中间计算、合并同类推导、只保留核心逻辑链、关键放缩不等式、最终系数结论与稳定性结果，符合顶刊正文篇幅与叙事节奏，不堆砌基础计算步骤。
二、MD详细手推版——内容详略标准（无跳步全推导）
该版本为完整推导底稿，覆盖 Theorem 1（一般边际稳定情形）、Theorem 2（幂零情形）全套五步证明，执行以下细则：
2.1 求导与展开无跳步
对 \(\dot V_1,\dot V_2\) 逐分项展开，所有乘积求导、分部积分、边界代入完整书写，不省略任意基础微积分步骤；分部积分明确写出积分上下限代入结果、余项消去过程。
2.2 不等式引用精准溯源
每一次使用 Young 不等式、Cauchy-Schwarz 不等式、Riccati 阶次性质，均标注：不等式形式、参数选取规则、放缩目的、被吸收的扰动项类型，不直接给出放缩结果。
2.3 所有常数系数完全显式推导
完整推导出所有核心常数的闭合表达式，包含：
- 放缩系数 \(C_1,C_3,C_4\)（扰动上界系数）
- 耗散系数 \(c_3\)（波动能量固有负定系数）
- 低增益上界 \(\varepsilon^*\) 的显式约束不等式
- Lyapunov 等价上下界 \(m,M\)（满足 \(m\Omega\le V\le M\Omega\)）
- 最终指数衰减率 \(\beta\) 的组合表达式
2.4 阶次分析逐行标注
对所有含 \(\varepsilon\) 扰动项，逐行标注阶次 \(O(\varepsilon),O(\varepsilon^2)\)，清晰区分可被吸收的高阶小项、主导耗散项，完整展示 \(\varepsilon^*\) 的选取逻辑与可行性证明。
2.5 幂零情形专属推导全覆盖
单独展开幂零矩阵 Riccati 估计 \(A^\top P_\varepsilon A\preceq C_A\varepsilon^2P_\varepsilon\) 的使用场景、额外耦合扰动项的处理过程、专属系数修正步骤，与一般情形推导形成完整对照。
2.6 适定性收尾完整
最后完整串联 Gronwall 不等式应用、解的存在唯一性论证、状态空间相容性核验，无逻辑断点。
三、TeX期刊正文精简版——内容取舍标准（投稿专用）
该版本为最终投稿正文，删减重复、琐碎、基础性计算，只保留学术核心内容，严格适配 IEEE TAC 行文范式：
3.1 保留内容（核心必写）
- 辅助状态变换、核心动力学解耦结论（省略逐阶求导展开）
- Lyapunov 泛函完整定义、泛函等价性核心结论
- 关键导数放缩的最终不等式形式（即 \(\dot V_1,\dot V_2\) 收尾公式）
- 低增益阶次分析核心逻辑与 \(\varepsilon^*\) 存在性结论
- \(V\sim\Omega\) 等价性结论、Gronwall 指数稳定最终结果
- 幂零情形的核心修正项与最终稳定结论
3.2 删除内容（全部省略，不占正文篇幅）
- 基础求导、分部积分的冗长展开步骤
- 不等式放缩的试探性推导、中间试算过程
- 各阶系数的原始推导草稿（仅保留最终显式结果）
- 阶次分析的重复铺垫话术
3.3 正文呈现形式
正文证明采用「逻辑步骤标题 + 关键公式 + 最终结论」的结构化写法，五步框架不变，但每一步仅输出可直接用于审稿的核心成果，干净紧凑、无冗余计算。
四、双版本对应关系与使用说明
1. 投稿交付：仅上传 TeX 精简版正文，符合期刊篇幅与审美，证明逻辑清晰、重点突出。
2. 备查交付：单独附带 MD 详细手推版，作为补充材料 / 作者底稿 / 审稿回复附件，全程无跳步、可逐行复现。
3. 系数统一性：两个版本所有常数、公式、结论完全一致，仅详略不同，不存在结论冲突、公式差异问题。
