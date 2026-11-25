---
title: "完全分裂素数的一条局部下界与数域 ERH"
date: "2026-11-25"
abstract: "从已核验的 Dirichlet GRH 出发，证明数域 ERH 等价于每个有限伽罗瓦扩张中完全分裂素数在 [y,4y] 上的一条单侧下界。对固定数域，还给出任意固定模数所能恢复的一次素理想数及其最小均方误差。"
---

# 完全分裂素数的一条局部下界与数域 ERH

设 $M/\mathbb Q$ 是次数 $n$ 的有限伽罗瓦扩张。要把 Alpöge–Claude 的定理推广到所有 Dedekind ζ 函数，只须对每个这样的 $M$ 和每个 $\varepsilon>0$ 证明

$$
 \sum_{\substack{p\ \text{在 }M\text{ 完全分裂}\\ p\ \text{未分歧}}}
 w(p/y)\log p
 \ \ge\ \frac{2\log 2}{n}\,y-C_{M,\varepsilon}y^{1/2+\varepsilon}
 \qquad (y\ge 2). \tag{1}
$$

这里的 $w$ 是下文写出的同一个非负函数，支撑在 $[1,4]$。反之，$M$ 的 ERH 蕴涵 (1)。取 $M$ 为判别式 $-23$ 三次域的六次闭包，就得到一个可直接用二元二次型描述的试验问题。伴篇给出该域的局部因子和高度 64 的零点证书。[三次域的计算与证书](cubic23-zeros.md)

## 已核验的起点与伽罗瓦闭包

Alpöge 于 2026 年 11 月 20 日公布的 Lean 题面遍历正整数模数 $N$、$N$ 的每个 Dirichlet 特征以及复数 $s$。`LFunction_ne_zero_of_half_lt_re` 给出 $\Re s>1/2$ 的非消失性，平凡特征在 $s=1$ 的极点另计。本地核读的社区维护版本为 `9c9f665ed79ce0fa8da61bd803b7a4bf1ad62ed0`，原始发布版的五个定理可与 Sutherland 的独立核验逐一对照。本文引用该形式化结果，下面的数域命题给出纸面证明。[作者的公告](https://alpo.ge/remark.pdf) · [原始 Lean 定理](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Main.lean) · [独立核验的题面与范围](https://math.mit.edu/~drew/grh/)

下文 $\operatorname{ERH}(K)$ 表示 $\zeta_K$ 在 $0<\Re s<1$ 的零点均有实部 $1/2$。若 $F/\mathbb Q$ 是有限阿贝尔扩张，Kronecker–Weber 定理和 Artin 分解给出

$$
 \zeta_F(s)=\prod_{\chi\in\widehat{\operatorname{Gal}(F/\mathbb Q)}}L(s,\chi^{*}).
$$

$\chi^{*}$ 是与 $\chi$ 对应、取自身导子的本原 Dirichlet 特征，平凡因子为 $\zeta(s)$。这条恒等式包括分歧处的局部因子。由 Dirichlet GRH 立即得到 $\operatorname{ERH}(F)$。[Milne，《Class Field Theory》，导论与第 V 章](https://www.jmilne.org/math/CourseNotes/CFT.pdf)

若 $K/\mathbb Q$ 任意，取伽罗瓦闭包 $M$。扩张 $M/K$ 是伽罗瓦扩张，Aramata–Brauer 定理给出 $\zeta_M/\zeta_K$ 为整函数。于是 $K$ 的每个临界带零点都是 $M$ 的零点，重数在上升时只会增加。对所有伽罗瓦域证明 ERH，便覆盖所有数域。对单个给定的 $K$，直接处理它的置换表示有时还能减少要控制的因子。[Aramata–Brauer 定理及其证明](https://doi.org/10.1007/s00605-025-02060-7)

还可以把全体数域的问题完全写成一维特征问题。对任一数域 $k$ 上的有限阶 Hecke 特征 $\psi$，取对应的有限阿贝尔类域 $N/k$，使 $\psi$ 经由 $\operatorname{Gal}(N/k)$ 分解。此时 $\zeta_N=\prod_\eta L_k(s,\eta)$，各因子在开临界带内全纯，每个 $L_k(s,\psi)$ 的零点都保留在乘积中。反过来，在每个数域上取平凡特征便恢复 Dedekind ζ 函数。因此全体数域的 ERH 等价于全体数域上有限阶 Hecke 特征的 RH。改变基域使理想类特征进入了题面。[Milne，Artin 互反律与一维 L 函数](https://www.jmilne.org/math/CourseNotes/CFT.pdf)

## 一个固定模数能恢复多少素理想信息

设 $G=\operatorname{Gal}(M/\mathbb Q)$，$H=\operatorname{Gal}(M/K)$。对未分歧素数 $p$，记 $r_K(p)$ 为 $K$ 中范数等于 $p$ 的素理想数。它等于 $\operatorname{Frob}_p$ 在 $G/H$ 上的不动点数。给定正整数 $q$，令 $F_q=K\cap\mathbb Q(\mu_q)$，$d_q=[F_q:\mathbb Q]$。对 $a\in(\mathbb Z/q\mathbb Z)^\times$，以 $a|_{F_q}$ 表示 $\zeta_q\mapsto\zeta_q^a$ 在 $F_q$ 上的限制。写

$$
 b_q(a)=d_q\,\mathbf 1_{a|_{F_q}=1},\qquad R=|H\backslash G/H|.
$$

**命题 1（最小均方误差）。** 对模 $q$ 的单位群上的任意固定复值函数 $u$，有精确极限

$$
\lim_{x\to\infty}\frac{1}{\pi(x)}
 \sum_{\substack{p\le x\\p\nmid qD_M}}
 |r_K(p)-u(p\bmod q)|^2
 =
 R-d_q+\frac{1}{\varphi(q)}
 \sum_{a\in(\mathbb Z/q\mathbb Z)^\times}|b_q(a)-u(a)|^2 . \tag{2}
$$

这里 $D_M$ 为域判别式，$\pi(x)$ 为有理素数计数函数。唯一的最优函数是 $b_q$。令 $F=K\cap\mathbb Q^{\rm ab}$ 为 $K$ 的最大阿贝尔子域。取 $q$ 为 $F$ 的导子的倍数，便达到所有固定模数中的最小误差 $R-[F:\mathbb Q]$。

**证明。** 置 $E=M\cap\mathbb Q(\mu_q)$，$J=\operatorname{Gal}(M/E)$。合成域的伽罗瓦群由限制到 $E$ 相同的两端自同构对组成。按自然密度形式的 Chebotarev 定理，条件 $p\equiv a\bmod q$ 下，任一共轭不变函数的均值即其在 $G$ 的对应陪集 $gJ$ 上的均值。作用在 $V=\mathbb C[G/H]$ 上时，对 $J$ 平均给出投影到 $V^J$，所以该陪集上的平均迹是 $g$ 在 $V^J$ 上的迹。$J$ 的轨道集是 $G/HJ$。由于 $G/J$ 阿贝尔，$HJ$ 正规，且 $M^{HJ}=K\cap E=F_q$。于是 $G/HJ$ 在自身上作正则作用，迹恰为 $b_q(a)$。[Milne，Chebotarev 定理 3.23 与自然密度说明 3.26(d)](https://www.jmilne.org/math/CourseNotes/CFT.pdf#page=173)

另一方面，$G$ 在 $(G/H)^2$ 上的轨道数是 $R$。按不动点计数，$r_K$ 的二次矩便是 $R$。$b_q$ 在单位群上的二次矩为 $d_q$，因为限制到 $F_q$ 的映射满射。在每个剩余类内展开平方，偏离条件均值的交叉项为零，得到 (2)。Kronecker–Weber 定理保证存在 $q$ 使 $F_q=F$。

特别地，$K/\mathbb Q$ 伽罗瓦时，最佳误差是 $|G|-|G^{\rm ab}|$。更一般地，只有阿贝尔数域能使最佳误差为零。若为零，两个类函数在每个群元素上取值相同，在单位元处就有 $[K:\mathbb Q]=[F:\mathbb Q]$。式 (2) 的极限先固定 $q$，故其结论容许其他方法选用随 $x$ 增长的模数。

## 去掉阿贝尔因子后的单侧判据

取

$$
 w(u)=
 \begin{cases}
  2-2/u & 1\le u\le 2,\\
  4/u-1 & 2<u\le 4,\\
  0 & \text{其余}.
 \end{cases}
 \qquad \int_1^4 w(u)\,du=2\log 2 . \tag{3}
$$

**定理 2。** 固定数域 $K$，仍取其最大阿贝尔子域 $F$ 和伽罗瓦闭包 $M$。在上述 Dirichlet GRH 下，$\operatorname{ERH}(K)$ 等价于对每个 $\varepsilon>0$ 存在 $C_{K,\varepsilon}$，对每个 $y\ge2$ 成立

$$
 S_K(y):=\sum_{p\nmid D_M}(r_K(p)-r_F(p))\,w(p/y)\log p
 \ \ge -C_{K,\varepsilon}y^{1/2+\varepsilon}. \tag{4}
$$

证明同样适用于 $K$ 中任意已知满足 ERH 的子域 $F$，包含 $F=\mathbb Q$。取最大阿贝尔子域的好处由 (2) 量化。

**证明。** 令 $B(s)=-\zeta'_K(s)/\zeta_K(s)+\zeta'_F(s)/\zeta_F(s)$。在 $\Re s>1$ 写作 $\sum b(m)m^{-s}$。系数支撑在素数幂上。若 $f_{\mathfrak p}$ 是 $\mathfrak p\mid p$ 的剩余次数，则 $-\zeta_K'/\zeta_K$ 在 $p^j$ 处的系数为 $(\log p)\sum_{f_{\mathfrak p}\mid j}f_{\mathfrak p}$，求和遍历满足此整除条件的素理想 $\mathfrak p\mid p$。对 $F$ 作差即得 $b(p^j)$。于是 $|b(p^j)|\le 2[K:\mathbb Q]\log p$，且未分歧处 $b(p)=(r_K(p)-r_F(p))\log p$。

写 $T(y)=\sum_m b(m)w(m/y)$。计数 $p^j\le 4y,\ j\ge2$，加回有限个分歧素数，可得 $T(y)-S_K(y)=O_K(\sqrt y\log^2(2y))$。若 $K$ 满足 ERH，对 $K,F$ 分别应用带误差的素理想定理，得到 $\sum_{m\le x}b(m)=O_K(\sqrt x\log^2(2x))$。对 (3) 分部求和即得 $|T(y)|=O_K(\sqrt y\log^2(2y))$，从而有 (4)。[Grenié 与 Molteni，ERH 下的素理想定理](https://arxiv.org/abs/1312.4463)

反向使用原 Lean 文件中的斜坡核。令 $U(y)=\sum_m b(m)(y/m-1)_+$，直接逐项相减得到 $T(y)=U(4y)-3U(2y)+2U(y)$。在 $\Re s>1$ 逐项积分，

$$
 \int_0^\infty T(y)y^{-s-1}\,dy
 =\Phi(s)B(s),\qquad
 \Phi(s)=\frac{(2^s-1)(2^s-2)}{s(s-1)}. \tag{5}
$$

$0,1$ 处取可去延拓。也可写成 $\Phi(s)=\int_1^4 w(u)u^{s-1}du$，所以 $\Phi$ 在整个实轴上严格为正。由乘积式，它在 $0<\Re s<1$ 中处处非零。尺度 $1,2,4$ 的组合还精确消除了 $U$ 中的线性项。[Alpöge–Claude 原始斜坡引理](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/DedekindLFunctions/LandauTheorems/MellinRamp.lean)

固定 $1/2<\theta<1$，用 (4) 中更小的 $\varepsilon$ 并吸收素数幂误差，可选 $C$ 使 $h(y)=T(y)+Cy^\theta\ge0$ 对一切 $y\ge1$ 成立。其 Mellin 积分 $J(s)=\int_1^\infty h(y)y^{-s-1}dy$ 在 $\Re s>1$ 等于

$$
 J(s)=\Phi(s)B(s)+\frac{C}{s-\theta}
       -\int_0^1 T(y)y^{-s-1}dy. \tag{6}
$$

末项整，因为 $T(y)$ 在充分靠近零处恒为零。若 $J$ 的收敛横坐标 $\sigma_c>\theta$，则 $\sigma_c\le1$，Landau 的正性论证要求实点 $\sigma_c$ 是奇点。具体地，若在该点解析，取稍右的 $\sigma_0$，在 $\sigma_0$ 对 $J(\sigma_0-t)$ 作 Taylor 展开，可取半径大于 $\sigma_0-\sigma_c$。各项是 $h(y)y^{-\sigma_0-1}(t\log y)^j/j!$ 的非负积分，单调收敛给出 $\sigma_c$ 左侧仍收敛，矛盾。

在 $s=1$，两条 Dedekind ζ 函数的单极点在 $B$ 中抵消。在 $\theta<s<1$ 的实轴上，$\zeta_F$ 无零点。若 $\zeta_K$ 在 $\sigma_c$ 有 $m$ 重零点，则 (6) 在此的留数为 $-m\Phi(\sigma_c)<0$，沿实轴从右趋近时 $J(s)\to-\infty$，与非负积分矛盾。其余实点 (6) 解析。因此 $\sigma_c\le\theta$，$J$ 在 $\Re s>\theta$ 全纯。由 (5)、(6) 的亚纯延拓以及 $\Phi$ 在开带内非零，$\zeta_K$ 在 $\theta<\Re s<1$ 无零点。令 $\theta\downarrow1/2$ 并用函数方程，即得 ERH。

现令 $K=M$ 并取 $F=\mathbb Q$。此时 $r_M(p)=n\,\mathbf1_{p\text{ 完全分裂}}$，$r_{\mathbb Q}(p)=1$。有理数的 RH 与 (3) 给出 $\sum_p w(p/y)\log p=2y\log 2+O(\sqrt y\log^2(2y))$。这就证明了 (1) 的等价性以及开头关于所有数域的结论。

## 从原形式化证明还需推进哪一步

原证明的 `QH.MeanSq` 对每个满足 $|\kappa(m)|\le1$ 的权成立，但其左侧还乘以 Jacobi 符号 $(m/q)$，并对辅助素数 $q$ 求均方。将 $b(p^j)/\log p$ 按域次归一化，可以放进这个有界权的位置。得到的是带辅助扭转的平均估计。最终使用的 `QH.Ind` 和 `QH.Box` 限定了本原 Dirichlet 特征及其在互素整数上的取值，`EveryCharacter` 的量词也仍是 Dirichlet 特征。判别式 $-23$ 的例子在同一剩余类中同时出现迹 $2$ 与 $-1$，式 (2) 对**任意固定模数**给出了这项差异。[原归档中的 MeanSq、Ind、Box 和 GoodQ](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/MeanSquareBounds/PrimeSetTwoBounds.lean) · [带窗口与辅助素数的均方定理](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/MeanSquareBounds/meanSqK.lean)

因此，一个具体的推广目标是在每个尺度 $y$ 控制 (4) 中**固定数域、去掉辅助扭转的 Frobenius 迹**。定理 2 允许上侧任意波动，所需的算术输入只是一致的下界。单个区间或若干未覆盖全部尺度的区间，要推出同一结论还需要相应的跨尺度论证。对第一个三次例子，这个下界有一个只用二元二次型表示素数的版本。[判别式 −23 的具体不等式](cubic23-zeros.md)
