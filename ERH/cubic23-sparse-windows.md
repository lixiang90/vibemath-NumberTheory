---
title: "从连续幂次窗口到判别式 −23 的离散素数判据"
date: "2026-11-25"
abstract: "把 Alpöge–Claude 的抽象窗口引理用于固定的三阶理想类特征。每个远端幂次窗口仅需 O((log Y)^{3/2}) 个高斯加权素数和，允许个别取样点有大误差。一次正权求积把离散的单侧平均送入原引理。"
---

# 从连续幂次窗口到判别式 −23 的离散素数判据

令 $P(X)=X^3-X-1$，$K=\mathbb Q(\alpha)$，$P(\alpha)=0$。以下只数能表示成 $r^2+rt+6t^2$ 的未分歧有理素数，每个素数计一次。这些素数恰在 $K$ 的六次闭包 $M$ 中完全分裂。前稿要求在每个充分大的尺度检验一条单侧下界。本篇证明，可以自由选择一列趋于无穷的幂次窗口，每窗取 $1+\lfloor(\log Y)^{3/2}/2\rfloor$ 个样点，将各点的短缺合计。各样点使用相互重叠的对数高斯权。[域的分解与二次型](cubic23-zeros.md) · [前稿的连续判据](splitting-prime-criterion.md)

## 每个窗口需要哪几个数

令 $k_v(t)=(2\pi v)^{-1/2}\exp(-t^2/(2v))$，其中 $v>0$。若 $u=\log y$，定义

$$
 P_1(u)=\sum_{p\ \text{完全分裂}}(\log p)\,k_1(u-\log p),\qquad
 D(u)=\frac{e^{u+1/2}}6-P_1(u).
$$

求和中始终除去 23。主项的常数来自完全分裂的密度 $1/6$ 和 $\int_{\mathbb R}e^t k_1(t)\,dt=e^{1/2}$。取 $U=\log Y\ge1$，$h=U^{-1/2}$，置

$$
 u_j=U+jh,\qquad 0\le j\le\left\lfloor\frac{U}{2h}\right\rfloor,\qquad
 H(U)=h\sum_j e^{-u_j/2}\,[D(u_j)]_+,\qquad [t]_+=\max(t,0). \tag{1}
$$

**定理。** 在已核验的 Dirichlet GRH 下，下列三项等价。

1. $\zeta_K$ 满足 ERH。
2. $\displaystyle\liminf_{U\to\infty}U^{-1}\log(1+H(U))=0$。
3. 将 (1) 中的 $[D]_+$ 全部改为 $[-D]_+$ 后，第二项成立。

三项也分别等价于 $\zeta_M$ 的 ERH。条件 2 只要求对每个 $\eta>0$ 有任意大的 $U$ 使 $H(U)\le e^{\eta U}$。窗口是 $[Y,Y^{3/2}]$，其中有 $1+\lfloor(\log Y)^{3/2}/2\rfloor$ 个取样尺度。每个素数和还可以截在 $|\log p-u_j|\le3\sqrt U$，保留同一判据。事实上，由 $\sum_{p\le x}\log p\ll x$ 将尾部分成单位对数区间，逐段用高斯衰减，在 $U\le u_j\le3U/2$ 上每点删去的总量为 $O(e^{u_j-9U/2+O(\sqrt U)})$。随 $Y$ 改变的这个素数范围说明了离散化所用的平滑程度。

## 接入原 Lean 证明的条件

记 $L=\zeta_K/\zeta$，$B=-L'/L=\sum_{n\ge2}b(n)n^{-s}$（$\Re s>1$），设

$$
 S_v(u)=\sum_{n\ge2} b(n)k_v(u-\log n).
$$

这里 $b(n)$ 为实数，在 $p^j$ 处为局部特征值的 $j$ 次幂和乘 $\log p$，绝对值至多 $2\log p$。素数一次项有 $b(p)/\log p=3\mathbf1_{p\ \text{完全分裂}}-\mathbf1_{\chi_{-23}(p)=1}$（$p\ne23$）。令 $Q_v$ 是对 $\chi_{-23}(p)=1$ 的素数所作的同权和。模 23 的 Dirichlet GRH 和素数定理给出 $Q_1(u)=e^{u+1/2}/2+O((1+u)^2 e^{u/2})$。例如，先对 $\sum_{p\le x,\,\chi(p)=1}\log p=x/2+O(\sqrt x\log^2(2x))$ 分部积分，再将权的导数对数化。高斯的各阶指数矩都有限。同样估计高次素数幂和有限个分歧项，得

$$
 -S_1(u)=3D(u)+O((1+u)^2e^{u/2})\qquad(u\ge1). \tag{2}
$$

对任意固定 $v>0$，$S_v$ 连续，且 $|S_v(u)|\le C_v e^{2u}$ 对所有实 $u$ 成立。这由 $k_v(t)\le C_v'e^{2t}$ 与 $\sum |b(n)|n^{-2}<\infty$ 得出。在 $u<0$ 它还具有高斯衰减。逐项变换给出

$$
 \int_{\mathbb R}[-S_2(u)]e^{-su}\,du=-e^{s^2}B(s),\qquad \Re s>1.
$$

减去 $u<0$ 的整函数项，即得 $g=-S_2$ 在 $[0,\infty)$ 上的 Laplace 变换。每个实部为正的 $L$ 的 $m$ 重零点 $\rho$ 对应一个留数为 $m e^{\rho^2}$ 的简单极点。此留数恒非零，$L(s)>0$ 对实 $s>0$ 成立，故实轴无需处理例外。$L$ 的零点离散，余部局部解析，恰满足原代码 `Deduction.PoleData` 的增长、隔离及留数条件。取 $g=+S_2$ 时非实留数变号，其余条件仍成立。[原始 PoleData 与 WindowHypU 的定义](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/DedekindLFunctions/LandauZeroFreeCriterion.lean) · [η 乘积给出的实轴正性](cubic23-zeros.md)

原定理 `Deduction.landau_window` 对满足这些条件的任意实函数 $g$ 适用。给定 $0<\theta<1$、$\delta>0$，若对每个 $\eta>0$ 都有任意大的 $V$，在 $[V,(1+\delta)V]$ 上存在一个非负上界 $R\ge g$，且 $\int R(u)e^{-\theta u}du\le e^{\eta V}$，则上述极点均满足 $\Re\rho\le\theta$。这是可直接引用的抽象结论。将新函数放入它的条件以及下述求积估计，均在本文中作解析证明。[原始窗口定理及其证明](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/DedekindLFunctions/LandauTheorems/landauWindow.lean)

## 为什么只控制离散点上的单侧误差就够了

对任意实数 $u,w,z$，配方得到

$$
 k_1(u-w)k_1(w-z)=k_2(u-z)k_{1/2}\big(w-(u+z)/2\big).
$$

在无限网格 $w_j=U+jh$ 上应用 Poisson 求和。$h\sum_j k_{1/2}(w_j-x)$ 与 1 的差的绝对值，对所有 $x$ 一致不超过 $2q/(1-q)$，其中 $q=\exp(-\pi^2/h^2)$。分别对每个 $z=\log n$ 应用上式并取绝对值，得

$$
 \left|S_2(u)-h\sum_{j\in\mathbb Z}k_1(u-w_j)S_1(w_j)\right|
 \le \frac{2q}{1-q}\sum_n |b(n)| k_2(u-\log n)
 \ll e^{2u-\pi^2 U}. \tag{3}
$$

在 $I=[9U/8,11U/8]$ 上，舍去 (1) 的有限网格以外的项，误差为 $O(e^{-cU^2})$，$c>0$。因所删点与 $u$ 的距离至少为 $U/8$，用 $|S_1(w_j)|\ll e^{2w_j}$ 并对高斯配方即可。这一步保留了每个网格权的非负性。因此对 $g=-S_2$ 可取 $R=[-S_2]_+$，由 (2)、(3) 得

$$
 \int_I R(u)e^{-u/2}du
 \le 3e^{1/8}H(U)+O(1+U^3). \tag{4}
$$

这里每个 $j$ 所付的因子至多为 $\int_{\mathbb R}k_1(t)e^{-t/2}dt=e^{1/8}$，(2) 的误差累计为 $O(U^3)$。取 $V=9U/8$，则 $I=[V,(1+2/9)V]$。条件 2 与 (4) 满足原定理在 $\theta=1/2$ 的全部窗口量词，于是 $L$ 的零点实部至多为 $1/2$，函数方程给出 RH。对 $+S_2$ 使用相同论证即得条件 3 的方向。卷积中的高斯 Mellin 因子在任何复点均非零，因而这个过程保留任意高度的偏线零点。[高斯求积的 Poisson 原理](https://doi.org/10.1137/130932132)

反向若 $\zeta_K$ 满足 ERH，则 $L$ 满足 RH，由 $\zeta_M=\zeta_{\mathbb Q(\sqrt{-23})}L^2$ 知闭包也满足 ERH。固定域的素理想定理给出 $\sum_{p\le x,\ \text{完全分裂}}\log p=x/6+O(\sqrt x\log^2(2x))$，因为一次未分歧素理想只来自完全分裂素数，每素数贡献六次。再对高斯权作分部积分得到 $|D(u)|=O((1+u)^2e^{u/2})$，故两个符号的网格总量都为 $O(1+U^3)$。[Grenié–Molteni 的素理想定理](https://arxiv.org/abs/1312.4463)

## 留下的算术问题

只要存在一个偏线零点，条件 2、3 就都失败。对每种符号，存在 $c>0$ 使**每个充分靠后的网格**的对应总量大于 $e^{cU}$。反过来，只须找到一列网格，使其中任一种符号的总量为 $e^{o(U)}$。这收紧了对误差发生频率和符号的要求。仓库中对二面体四次域的研究通过改变特征、对许多导子求均方得到多数域的结果。这里改变的是同一数域的取样尺度，必须为判别式 $-23$ 本身取得无界窗口序列。[已有的二面体域族平均](../zeros-exact-distribution/dihedral-quartics-central-strip.md)

对于同一个 $P_1$，另一篇给出在 $30\le u\le1048476$ 上的显式平方根界，以及素数平方项导致的有符号次项。它验明了很长的连续尺度范围。条件 2 或 3 中的无穷序列仍是从该有限结论推进到全高度 ERH 所需的算术输入。[有限尺度上的素数估计](cubic23-gaussian-bias.md)

若允许高斯方差随 $\log Y$ 增长，可以在每个整数对数尺度取一个样点，并把每段的长度放宽到任意慢地趋于无穷。所需的最小方差斜率可由零点坐标精确表示。[全尺度估计与最小斜率](cubic23-gaussian-threshold.md)
