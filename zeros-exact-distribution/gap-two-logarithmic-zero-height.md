---
title: "对数高度的 Dirichlet 零点与稀疏孪生素数的检测"
date: "2026-11-24"
abstract: "构造对数坐标中的紧支撑权。固定差 2 的完整线性显式公式只需高度 O(log X) 的零点及两个完全幂修正，即可给出与孪生素数无穷性等价的判据。对无平方邻项的小除数部分，利用 Montgomery 的零点密度定理把截高无条件降至 O(log log X)，并保留奇数因子消去和四因子项的符号。"
---

# 对数高度的 Dirichlet 零点与稀疏孪生素数的检测

对尺度 $X$，可以只用高度 $O(\log X)$、导子至多 $2X+2$ 的本原 Dirichlet 零点，连同两个可逐项计算的完全幂修正，构造孪生素数无穷性的等价判据。它允许有孪生素数的窗口任意稀疏。若窗口中根本没有一对，统计量至多是 $O(\log X)$。一对落在窗口中段的真孪生素数能产生 $\gg(\log X)^2$ 的信号。整个判据无条件。

相同的权还适用于第 24 篇的互补除数方法。用平方因子展开直接选出全部无平方邻项，可把所需的小除数模数限制在 $\sqrt X(\log X)^{64}$ 的常数倍内。在这个范围，素数在等差数列中的分布给出跨模数的平均上界，Montgomery 的零点密度定理将截高无条件降至 $O(\log\log X)$。这个较低截高只用于已经可以估计的小除数部分。完整判据和平均上界分别见 (7)、(12)。[前篇的互补除数分解](https://lixiang90.github.io/vibemath-NumberTheory/zeros-exact-distribution/articles/gap-two-complementary-divisors.html)

作为比较，Korevaar 用固定移位 Dirichlet 级数的边界性质刻画 Hardy–Littlewood 渐近，Keating–Smith 从 $\zeta$ 零点对关联出发作形式傅里叶反演。这里的判据只要求无穷性，输入是随 $X$ 增长的一族已知 $L$ 函数的单个零点，截高与误差都可在这族函数上统一控制。[Korevaar，推论 2.1](https://staff.fnwi.uva.nl/j.korevaar/DisWieIke.pdf) · [Keating–Smith，素数对与零点对关联](https://arxiv.org/abs/1903.07057)

## 1. 从原证明的 Mellin 权得到两个截高

Alpöge–Claude 归档中的引理 mellin_ramp 对 $\Re s>1$ 将 $(y/n-1)_+$ 的 Mellin 积分化为 $n^{-s}/[s(s-1)]$。对每个固定的基准尺度，把 $y$ 改成 $ye^{-v}$，再按偶的概率密度 $H(v)$ 积分，就在变换侧乘以 $\Phi(s)=\int H(v)e^{-sv}\,dv$。这里对 $0<y<\infty$ 积分。若 $H$ 支撑在 $[-h,h]$，则在 $Ye^h\le y\le Y^A e^{-h}$ 上，所有取样点仍在原幂次窗口 $[Y,Y^A]$ 内。原误差 $C y^\theta$ 经平均至多多出常数 $e^{|\theta|h}$。[原引理及其精确积分范围](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/DedekindLFunctions/LandauTheorems/MellinRamp.lean)

本文用两次不同阶数的卷积，给同一个紧支撑权提供两个高度范围的控制。设 $L=\log X$，$g=1/64$，$U_a$ 是 $[-a,a]$ 上的均匀概率密度。以下 $X$ 趋于无穷，取

$$
N_1=\max(2,\lceil200\log L\rceil),\qquad
N_2=\max(2,\lceil20L\rceil),\qquad
H=U_{g/N_1}^{*N_1}*U_{g/N_2}^{*N_2},\qquad h=2g .
\tag{1}
$$

$H$ 的支撑在 $[-h,h]$，方差是 $\sigma_H^2=g^2(N_1^{-1}+N_2^{-1})/3$。写 $S(z)=\sinh(z)/z$、$S(0)=1$，则
$\Phi(s)=\prod_{i=1}^2 S(gs/N_i)^{N_i}$。它在开带 $0<\Re s<1$ 内处处非零，因此对固定 $X$ 平滑后的 Mellin 变换仍保留这个带内原有极点及其阶数。反复卷积得到的双曲正弦权有既有文献。这里使用两个增长阶数，并将同一密度用于除数配对。[Thorner–Zaman，引理 4.2](https://arxiv.org/pdf/1604.01750#page=16)

令 $a_0=\log(5/4)$，$b_0=\log(8/5)$，$f=\mathbf1_{[a_0,b_0]}*U_g*U_g$。定义具有有限高阶光滑性的窗口权及其 Mellin 变换

$$
W_X(e^x)=\frac{(f*H)(x)}{c_X},\qquad
\widehat W_X(s)=\frac{e^{b_0s}-e^{a_0s}}{s\,c_X}S(gs)^2\Phi(s),
\qquad \widehat W_X(1)=1 .
\tag{2}
$$

最后一个条件确定 $c_X$，它有与 $X$ 无关的正上、下界。权支撑在 $1<u<2$，在 $4/3\le u\le3/2$ 上恒为 $1/c_X$，因为 $\log(16/15)>4g$。$W_X$ 及其一阶导数一致有界。令 $W_0$ 是只用 $f$ 并归一化的固定权，则 $\|W_X-W_0\|_\infty\ll\sigma_H^2$，这是由于 $f'$ 为 Lipschitz 函数，且 $H$ 的一阶矩为零。用素数对上界筛可知，将两者用于窗口孪生素数质量，差为 $O(X/\log L)$。下文的稀疏判据直接使用 (2) 的平顶，避免了这项近似误差。[固定差素数对的上界筛](https://terrytao.wordpress.com/2015/01/21/254a-notes-4-some-sieve-theory/)

**导数界。** 对 $N=N_1$ 或 $N_2$，可把 $j\le N$ 次微分各放在这一组的一个均匀因子上。单个因子的分布导数是两个点质量之差，全变差为 $N/g$。于是 $\|D^jH\|_{\rm TV}\le (N/g)^j$，这里 $j=0$ 也成立。特别地，沿 $0\le\beta\le1$，在对数坐标分部积分 $N$ 次给

$$
|\widehat W_X(\beta+i\gamma)|\ll (97N/|\gamma|)^N,\qquad |\gamma|\ge1 .
\tag{3}
$$

常数 $97$ 也覆盖第 4 节额外的半尺度因子。紧支撑消去边界项。对每个 $\chi^*$，在 $T<|\Im\rho|\le2T$ 的零点数为 $O(T\log(2qT))$，重数在内。将 (3) 对二进高度段相加，在 $T=256N$ 以后得到 $O(T\log(2qT)\,2^{-N})$ 的变换尾界。[Bennett、Martin、O’Bryant 与 Rechnitzer，Dirichlet 零点计数](https://arxiv.org/abs/2005.02989)

## 2. 完整谱和中的两个有限修正

令 $D=\lfloor2X+2\rfloor$，$a_d=-\mu(d)\log d/\varphi(d)$（$d$ 为奇数），$C(D)=\sum_{d\le D,\,2\nmid d}a_d$。对 $\chi\bmod d$ 写 $\chi^*$ 为诱导本原特征，主特征取 $\zeta$。定义

$$
Z_T(X)=\sum_{\substack{d\le D\\2\nmid d}}a_d
\sum_{\chi\bmod d}\overline{\chi(-2)}
\sum_{\substack{\rho\ {\rm of}\ L(s,\chi^*),\ 0<\Re\rho<1\\|\Im\rho|\le T}}
X^\rho\widehat W_X(\rho).
\tag{4}
$$

正负高度均计，按零点重数求和。共轭特征配对使 (4) 为实数。记 $m_{\rm odd}$ 为 $m$ 的奇数部分，$\Lambda(1)=0$，写 $J(m)=\Lambda(m_{\rm odd})$。设 $\mathcal P_X$ 是 $X<n<2X$ 中 $n$ 或 $n+2$ 是完全幂的所有整数 $n$。完全幂指 $b^k$、$b,k\ge2$，集合取并集并去重。两个修正为

$$
\begin{aligned}
P_X&=\sum_{n\in\mathcal P_X}\Lambda(n)J(n+2)W_X(n/X),\\
I_{\ge2}(X)&=\sum_{\substack{\ell\ {\rm odd\ prime},\,k\ge2\\X<\ell^k<2X}}
(\log\ell)W_X(\ell^k/X)
\sum_{\substack{e\mid\ell^k+2,\ 2\nmid e\\e\ell\le D}}
\frac{-\mu(\ell e)\log(\ell e)}{\ell-1}.
\end{aligned}\tag{5}
$$

在 (5) 中，待查的完全幂端点只有 $O(\sqrt X)$ 个，连同其邻项和有关除数均可由整数运算确定，无需输入任何真孪生素数对的位置。(4) 的模数范围仍随 $X$ 线性增长。

**定理 1.** 写 $M_X=\sum_{p,\,p+2\ {\rm prime}}(\log p)\log(p+2)\,W_X(p/X)$。无条件地有

$$
\mathcal Y(X):=C(D)X-Z_{256N_2}(X)-P_X-I_{\ge2}(X)=M_X+O(L).
\tag{6}
$$

**证明。** 对任意 $m$ 有 $J(m)=-\sum_{d\mid m,\,2\nmid d}\mu(d)\log d$。由此，对 $A_X=\sum_n\Lambda(n)J(n+2)W_X(n/X)$ 使用特征正交，余类为 $-2\bmod d$。若 $n$ 为奇数，$J(n+2)=\Lambda(n+2)$。偶数 $n$ 要有非零贡献则必为 $2$ 的幂，充分大时已在 $\mathcal P_X$ 中。因此 $A_X=M_X+P_X$ 恰好成立。

先对本原函数使用光滑显式公式。$s=1$ 的极点产生 $C(D)X$，非平凡零点产生 $-Z_\infty$。平凡零点总修正为
$V=(\log3)/2\int W_X(t/X)/(t-1)\,dt-C(D)\widehat W_X(0)$，被积函数在窗口外延为零。因为奇偶投影后各模数的相位和是 $\varphi(d)(\mathbf1_{d\mid3}\pm\mathbf1_{d\mid1})/2$，且 $\zeta(0)\ne0$，故只有这两项。$C(D)=\mathfrak S(2)+o(1)$，其中 $\mathfrak S(2)=2\prod_{\ell>2,\ \ell\ {\rm prime}}(1-1/(\ell-1)^2)$，见第 4 节的 Euler 乘积，因而 $V=O(1)$。测试函数具有多于两阶的紧支撑正则性，足以用光滑逼近应用显式公式。[Dirichlet 特征的显式公式](https://kskedlaya.org/ant/chap-errorbounds2.html) · [本原函数的函数方程和平凡零点](https://dlmf.nist.gov/25.15)

非本原修正中写 $d=\ell e$。因 $\mu(d)\ne0$，$\ell\nmid e$。仅 $\ell$ 处分量为主特征时缺少 Euler 因子。对模 $e$ 的特征正交，给 $\varphi(e)\mathbf1_{e\mid\ell^k+2}$，于是 $k\ge2$ 正好贡献 (5) 的 $I_{\ge2}$。$k=1$ 时 $\ell>X$ 强制 $e=1$，剩下的修正为
$I_1=\sum_{\ell\ {\rm prime}}(\log\ell)^2 W_X(\ell/X)/(\ell-1)=O(L)$，最后一步用素数计数的切比雪夫上界。完整等式为 $M_X=C(D)X-Z_\infty-P_X-I_{\ge2}-I_1-V$。

最后，(4) 的外层绝对权 $\sum_{d\le D}|a_d|\varphi(d)\ll XL$（仍仅计奇 $d$）。对每个零点取 $|X^\rho|\le X$。由 (3)，在 $256N_2$ 以上的总尾部为 $O(X^2L\,N_2\log(2DN_2)\,2^{-N_2})=O(X^{-8})$。舍去 $I_1,V$ 即得 (6)。其中 $I_{\ge2}$ 可以带负号，且 $P_X$ 可超过单对信号，故二者都在统计量中显式保留。

## 3. 一个允许孪生素数任意稀疏的等价判据

令 $X_j=(17/16)^j$。则

$$
\boxed{\quad \#\{p:p,\ p+2\ {\rm prime}\}=\infty
\quad\Longleftrightarrow\quad
\mathcal Y(X_j)>(\log X_j)^{3/2}\ \text{对无穷多个 }j\text{ 成立}.\quad}
\tag{7}
$$

因为 $(4/3)(17/16)<3/2$，相邻的平顶区间 $[4X_j/3,3X_j/2]$ 相互覆盖，无限延伸。任何落在平顶的一对给 $M_{X_j}$ 至少 $c(\log X_j)^2$，其中 $c>0$ 固定。它最终超过阈值和 (6) 的误差之和。反过来，超过阈值的充分晚窗口必含真孪生素数，且一对只能落在有限个这样的窗口。两向都无需假设正密度。

若目标改为正比例的窗口质量，则 $P_X\ll\sqrt X L^2$，且由 $\tau(m)\ll_\nu m^\nu$ 有 $I_{\ge2}\ll_\nu X^{1/2+\nu}$。因此 $\limsup_j M_{X_j}/X_j>0$ 等价于 $\liminf_j Z_{256N_2}(X_j)/X_j<\mathfrak S(2)$。在给定尺度判出一对与估计窗口质量，对算术修正精度的要求不同。

## 4. 平方因子直接展开后可估计的小除数零点族

令 $K=H*(x\mapsto x_+)$，并令

$$
T_H(m)=\sum_{d\mid m}\mu(d)K\left(\tfrac12\log m-\log d\right),\qquad
A_{\rm sf}=\sum_p\mu^2(p+2)(\log p)W_X(p/X)T_H(p+2).
\tag{8}
$$

$A_{\rm sf}$ 包含全部无平方邻项。取 $U=L^{32}$，按 $\mu^2(m)=\sum_{c^2\mid m}\mu(c)$ 展开。对无平方的 $c$ 和内层无平方除数 $d$，写 $v=(d,c)$、$d=bv$。有 $v\mid c$，$(b,c)=1$，同时满足 $c^2\mid m,\ bv\mid m$ 对应的模数恰是 $q=c^2b$。从 $q$ 的指数为二的素因子可唯一恢复 $c$。因此可以先合并 $v$ 的贡献，再对模数取绝对值。

令 $\mathcal B$ 为奇无平方数 $c\le U,\ b\le e^h\sqrt{2X}$、$(b,c)=1$ 的所有有序对 $(b,c)$，并写

$$
G_{b,c}(u)=W_X(u)\sum_{v\mid c}\mu(v)
             K\left(\tfrac12\log(Xu)-\log(bv)\right),\qquad
Q=e^h\sqrt{2X}\,U^2 .
\tag{9}
$$

这里用 $\sqrt p$ 代替算术式 (8) 中的 $\sqrt{p+2}$。$K$ 为 1-Lipschitz，在原无平方和中作此替换，总误差至多 $O(L X^{-1}\sum_{m\le2X+2}\tau(m))=O(L^2)$。当 $b>e^h\sqrt{2X}$ 时 (9) 恒为零。删去 $c>U$ 造成的总误差为 $O(X L^{19/2}/\sqrt U)$。事实上，含某个这样的平方因子的整数至多 $O(X/U)$ 个，逐整数的总绝对权至多 $O(L^2\tau(m)^2)$，再用 Cauchy–Schwarz 及 $\sum_{m\le 3X}\tau(m)^4\ll XL^{15}$ 即得。最后一个界可由 $\tau^4\le\tau_{16}$ 及十六重除数和得到。

**定理 2.** 对每个固定 $\epsilon>0$，有

$$
A_{\rm sf}=\mathfrak S(2)X+O_\epsilon(X/L^{2-\epsilon}).
\tag{10}
$$

**等差数列误差。** 展开后的外权是 $\mu(c)\mu(b)$，内和是 $p\equiv-2\bmod q$ 上的 $(\log p)G_{b,c}(p/X)$。在 $q\le Q_0=\sqrt X/L^B$ 用 Bombieri–Vinogradov，$B$ 选得充分大。对更大的模数，$Q_0<q\le Q\ll \sqrt X L^{64}$，且 $b=q/c^2\ge Q_0/U^2$。所以 $K$ 在非零支撑上大小为 $O(\log L)$。在 $1\le u\le2$ 上，$\log(Xu)G_{b,c}(u)$ 的导数有共同上界 $O(\tau(c)L\log L)$，小模数的相应界为 $O(\tau(c)L^2)$。$\tau(c)\ll_\nu L^\nu$ 对任意固定 $\nu>0$ 成立。

在 $Q_0<q\le Q$ 的每个二进段，对每个固定 $t=Xu$，Bombieri–Friedlander–Iwaniec 第三篇的固定剩余类定理给 $\sum_q|\pi(t;q,-2)-\pi(t)/\varphi(q)|\ll_B X(\log L)^2/L^3+O_B(X/L^{10})$。这里 $|\log(2V)/\log t-1/2|\ll_B(\log L)/L$，$V$ 为段的左端，共有 $O_B(\log L)$ 段。先对每个 $u$ 求和，再对 $u$ 作分部积分，连同普通素数定理，有

$$
\sum_{(b,c)\in\mathcal B}\left|
\sum_{p\equiv-2\bmod c^2b}(\log p)G_{b,c}(p/X)
-\frac{X}{\varphi(c^2b)}\int G_{b,c}(u)\,du
\right|\ll_\epsilon X/L^{2-\epsilon}.
\tag{11}
$$

大模数部分具体为 $O_\nu(XL^{-2+\nu}(\log L)^4)$。取 $\nu<\epsilon$。模数的唯一性保证不再付出分解重数。[Bombieri–Friedlander–Iwaniec 原定理](https://doi.org/10.1090/S0894-0347-1989-0976723-6) · [Fiorilli 所列的区间形式，定理 1.2](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/6F019883E34D9946AFD27A30A4061440/S0008414X00009342a.pdf/on_a_theorem_of_bombieri_friedlander_and_iwaniec.pdf#page=2)

**主项。** 记 $V_c=\prod_{\ell\mid c}(1-1/(\ell-1))$，乘积取素数。对 $c\le U$ 以及 $\sqrt X e^{-h}/U\le y\le\sqrt{2X}e^h$，有

$$
\sum_{\substack{b\le y\\(b,2c)=1}}\frac{\mu(b)}{\varphi(b)}\log(y/b)
=\frac{\mathfrak S(2)}{V_c}+O_A(V_c^{-1}L^{-A}).
$$

为证此式，先取 $c=1$。奇数上的 $\sum\mu(n)n^{-s}/\varphi(n)=J_0(s)/\zeta(1+s)$，其中 $J_0$ 的系数在 $\Re s=-1/2$ 绝对可和，$J_0(0)=\mathfrak S(2)$。在奇素数 $\ell$ 处，$J_0$ 的因子为 $1-\ell^{-s}/[\ell(\ell-1)(1-\ell^{-1-s})]$，在 $2$ 处为 $(1-2^{-1-s})^{-1}$。将系数卷积在 $\sqrt y$ 处分开，利用 $\sum_{n\le t}\mu(n)\ll t\exp(-c_0\sqrt{\log t})$ 分部求和，即得所述任意对数精度，也给出第 2 节所用的 $C(D)$ 极限。再剔除 $\ell\mid c$，新增乘积 $\prod_{\ell\mid c}(1-\ell^{-s}/(\ell-1))^{-1}$ 的系数非负，总质量 $1/V_c$。在 $s=-1/(10\log U)$ 处，此乘积为 $(\log U)^{O(1)}$，故超过 $\sqrt y$ 的新增卷积尾部指数小于 $-c_1 L/\log U$，多项式对数因子可以吸收。这就给出对 $c$ 的一致性。

在 $K$ 的积分中取 $y=\sqrt{Xu}\,e^{-w}/v$，$|w|\le h$。上式的主项与 $v\mid c$ 无关，所以对每个 $c>1$ 都被 $\sum_{v\mid c}\mu(v)=0$ 精确消去。误差可求和，因为奇无平方数上的 $\sum_c \tau(c)/(\varphi(c^2)V_c)\le\prod_{\ell>2}(1+2/(\ell(\ell-2)))<\infty$，此处乘积仍取素数。$c=1$ 给出 $\mathfrak S(2)X\int W_X=\mathfrak S(2)X$。结合 (11) 和平方因子截断误差，得到 (10)。

**零点形式。** 在 (4) 的单模数统计中将 $\widehat W_X$ 换成 $\widehat G_{b,c}(s)=\int G_{b,c}(u)u^{s-1}du$，去掉 $a_d$ 并以 $1/\varphi(q)$ 归一化。所得量记为

$$
B_{b,c}(T)=\frac1{\varphi(q)}\sum_{\chi\bmod q}\overline{\chi(-2)}
\sum_{\substack{\rho\ {\rm nontrivial\ of}\ L(s,\chi^*)\\|\Im\rho|\le T}}
X^\rho\widehat G_{b,c}(\rho),\qquad q=c^2b .
$$

同一组测试函数无条件地满足

$$
\boxed{\displaystyle
\sum_{(b,c)\in\mathcal B}|B_{b,c}(256N_1)|\ll_\epsilon X/L^{2-\epsilon},
\qquad 256N_1=O(\log\log X).}\tag{12}
$$

先取全部高度。由显式公式，(11) 和零点式之间的高次素数幂修正汇总为 $O_\nu(X^{1/2+\nu})$，因为每个 $\ell^k\asymp X$ 至多涉及 $\tau(\ell^k+2)$ 个不同模数。每模数的缺失因子和平凡零点至多贡献 $O(\tau(c)L\log(2q))$。共有至多 $Q$ 个模数，故这些修正均可吸收。

对截高，在 (9) 的对数坐标中分别微分 $e^{\beta x}$、$W_X(e^x)$ 和 $K((L+x)/2-\log(bv))$。后一个因子的 $j$ 阶导数绝对值至多 $O(L(N/(2g))^j)$，可把导数放在均匀因子上，剩下的正部函数在相关范围内 $\ll L$。多项式展开导数给 $(1+N/g+N/(2g))^N\le(97N)^N$。故对 $0\le\beta\le1$，
$|\widehat G_{b,c}(\beta+i\gamma)|\ll \tau(c)L(97N/|\gamma|)^N$。对截尾可粗用 $\tau(c)\le L^{32}$。若所有零点都在中心线上，用 $|X^\rho|=\sqrt X$ 一步就给 $O(X/L^4)$ 的尾界。以下去掉这个假设。

写 $V=256N_2\asymp L$。高度超过 $V$ 的总贡献无条件地 $\ll XQ L^{33}V\log(2QV)2^{-N_2}\ll X/L^4$。在中间高度 $256N_1<|\gamma|\le V$，同一核的 $N_1$ 阶界可与下列族平均配合使用。

$$
\sum_{q\le Q}\frac1{\varphi(q)}\sum_{\chi\bmod q}
     \sum_{\substack{\rho=\beta+i\gamma\ {\rm nontrivial\ of}\ L(s,\chi^*)\\|\gamma|\le V}}
      X^\beta\ll XL^{85},\qquad Q=e^h\sqrt{2X}\,L^{64}.
$$

为证此界，用 $X^\beta\le\sqrt X+L\int_{1/2}^{\max(\beta,1/2)}X^\sigma\,d\sigma$。第一项对全部零点求和，由零点计数为 $O(\sqrt X Q V\log(2QV))$。对积分，固定导子 $F\le f<2F$。一个本原特征在每个 $f\mid q$ 的模数至多出现一次，诱导它的总权至多
$\sum_{k\le Q/f}1/\varphi(fk)\ll L/\varphi(f)\ll L^2/F$。Montgomery 的密度定理给

$$
\sum_{f\le2F}\sum_{\chi\bmod f}^{*}N_\chi(\sigma,V)
\ll(F^2V)^{a(\sigma)}(\log(2FV))^{13},\qquad
a(\sigma)=\frac{3(1-\sigma)}{2-\sigma},\quad \tfrac12\le\sigma\le1 .
$$

这里 $N_\chi$ 计 $\beta\ge\sigma$、$|\gamma|\le V$ 的零点，本原特征求和用星号表示。若 $a\ge1/2$，则
$X^\sigma F^{2a-1}\ll X^{\sigma+a-1/2}L^{64(2a-1)}\le XL^{64}$，因为 $1/2\le\sigma\le1$ 上 $a\le1$，且 $3/2-\sigma-a=\sigma(\sigma-1/2)/(2-\sigma)\ge0$。若 $a<1/2$，直接用 $F^{2a-1}\le1$。再计入 $V^a\ll L$、诱导权中的 $L^2$、密度定理的 $L^{13}$、$O(L)$ 个导子段和分层积分前的 $L$，得到上述界（指数 $85$ 留有余量）。导子 $1$ 也可直接用 $\zeta$ 的零点计数处理。[Montgomery，原密度定理](https://doi.org/10.1007/BF01404638) · [Inoue，引理 1，含本原特征和及对数指数](https://www.numdam.org/item/10.5802/jtnb.1162.pdf#page=10)

在中间高度逐项使用 $\tau(c)L\,2^{-N_1}$，由模数唯一性，该段总误差 $\ll XL^{118}2^{-N_1}\ll X/L^4$，因为 $200\log 2>122$。这与高段尾界及 (11) 共同证明 (12)。与 (4) 相同，截高以内仍保留模数内部的相位 $\overline{\chi(-2)}$。

## 5. 平滑对奇偶符号和单边估计的影响

$H$ 为偶概率密度，故 $K(u)-K(-u)=u$。无平方 $m>1$ 的互补除数配对仍给 $(1-\mu(m))T_H(m)=\Lambda(m)$。因此奇数个素因子的合数严格为零，素数项为 $\frac12\log m$。在偶数层还成立以下两个符号结论。

若 $m$ 有 $r\ge2$ 个不同素因子，其对数依次为 $a_1,\ldots,a_r$，记 $J_a=\mathbf1_{[-a/2,a/2]}$。在 (8) 中依次取中心差分，利用 $f(x+a/2)-f(x-a/2)=(J_a*f')(x)$ 和 $(x_+)''=\delta_0$，得到

$$
T_H(m)=
\left.\frac{d^{\,r-2}}{dx^{\,r-2}}
(H*J_{a_1}*\cdots*J_{a_r})(x)\right|_{x=0}.
\tag{13}
$$

在 (13) 中 $m$ 无平方。$r=2$ 时 $0\le T_H(m)\le\min(a_1,a_2)$，两因子对数之差至少 $2h$ 时上界取到。$r=4$ 时卷积为偶函数，且在原点取最大值，所以二阶导数 $\le0$。这里每个因子都是关于原点非增的偶密度或区间指示函数，卷积保留这种性质，可由中心区间交长及层集积分直接验证。对 $r\ge3$，若最大素因子至少为 $e^h\sqrt m$，则含它和不含它的两组除数分别落在 $K=0$ 和 $K(u)=u$ 的区域。后一区域的 Möbius 和与对数加权 Möbius 和均为零。

还有一个对任意因子分类都适用的比较。记 $T(m)$ 为将 (8) 中 $K$ 换回 $u_+$ 的原始权，则

$$
\sum_p(\log p) W_X(p/X)|T_H(p+2)-T(p+2)|
\ll X\sigma_H^2+\frac{X(\log L)^2}{L^2}
\ll X/\log L .
\tag{14}
$$

事实上 $\Delta(u)=K(u)-u_+\ge0$，支撑在 $[-h,h]$，为一致 Lipschitz 函数，$\int\Delta=\sigma_H^2/2$。最后一式可对 $w,-w$ 配对得到三角形 $(|w|-|u|)_+/2$ 后积分。以各除数的非负 $\Delta$ 贡献控制左端，只出现 $d\asymp\sqrt X$。在这些常数个二进段应用 (11) 证明所用的固定剩余类定理，测试函数 $W_X(u)\Delta(\frac12\log(Xu+2)-\log d)$ 的一阶导数几乎处处一致有界。分布误差为 $O(X(\log L)^2/L^2)$。主项中每个 $d$ 的积分 $\ll\int\Delta$，且 $\sum_{d\asymp\sqrt X}1/\varphi(d)\ll1$。因此，即使对素数 $p$ 再作任意子集限制，这次平滑也只引入 $o(X)$ 的绝对误差。

由 (8)、(10)、(13)，$M_X$ 等于 $2\mathfrak S(2)X$ 减去所有偶数因子层的两倍带符号质量，误差 $O_\epsilon(X/L^{2-\epsilon})$。四因子项仍有利于下界。由 $T_H(ab)\le\log\min(a,b)$ 和前篇的两线性式上界筛，较小素因子不超过 $X^\eta$ 的半素数仅耗去 $O(X(\eta+L^{-1}))$，$0<\eta\le1/4$。两因子均大的半素数和偶数高层仍需联合估计。(7) 给出了对稀疏序列也有效的低高度判据，(12) 估计了其中一条算术分解的小除数部分。使 (7) 的阈值无穷多次被越过，或使上述偶数层低于 $\mathfrak S(2)X$ 一个固定比例，仍待证明。

## 6. 数值和精确算术核验

对 $3\le m\le20000$ 的全部 9999 个奇数，逐除数用整数核对 $c^2b$ 的分组系数和 $\mu^2(m)\mu(d)$ 一致。另以两阶、三阶均匀卷积组成一个可精确积分的有理数模型，在 627 组正长度上核对 (13) 的二、三、四、五因子结论。

对 (2) 的实际窗口，分别从真素数、奇除数和、本原与诱导特征的差计算。下表不作除以 $X$ 的归一化。$X=5$ 说明，在没有窗口质量时遗漏 $P_X$ 也会造出很大的正信号。$X=40$ 则检验了缺失 Euler 因子修正的负号情形。

| $X$ | $M_X$ | $P_X$ | $I_{\ge2}$ | $I_1$ |
|---:|---:|---:|---:|---:|
| 4 | 4.47365747 | 0 | 0 | 0.92502652 |
| 5 | 0 | 7.70104473 | 0 | 1.80297709 |
| 20 | 33.03514240 | 15.89967412 | 1.09264141 | 2.27963407 |
| 40 | 47.88818926 | 约 0 | −0.06434043 | 2.48966444 |
| 100 | 141.08477058 | 11.13682780 | 0.92500658 | 3.43445861 |

另取 $X=4$，直接读取 $\zeta$ 的前 2469 个正高度零点，连同其共轭计算 $\sum_n\Lambda(n)W_X(n/X)$。对应高度 $64,256,1024,3000$ 与独立算术值 $2.29900515943$ 的绝对偏差依次为 $0.259$、$0.00104$、$4.21\times10^{-9}$、$6.93\times10^{-14}$。这一项检验 (2) 的 Mellin 归一化、零点重数和正负高度求和。

最后，对 $X=4,5$ 把 (4) 中的**全部特征**纳入计算。在 $-1/2\le\Re s\le3/2,\ |\Im s|\le T$ 的矩形上对 $X^s\widehat W_X(s)L'/L(s,\chi^*)$ 作正向围道积分，按留数定理扣除 $s=0,1$ 的贡献。$X=4$ 涉及模 $3,5,7$，$X=5$ 还涉及模 $11$。表中的 $\mathcal Y_T$ 将 (6) 的截高改成给定的 $T$，加权重建值再扣去 $I_1,V$。

| $X$ | $T$ | $\mathcal Y_T$ | 加权重建值 | 直接计数的 $M_X$ | $(\log X)^{3/2}$ |
|---:|---:|---:|---:|---:|---:|
| 4 | 32.25 | 4.19139115 | 3.69528868 | 4.47365747 | 1.63223687 |
| 4 | 256.25 | 4.96719888 | 4.47109641 | 4.47365747 | 1.63223687 |
| 5 | 32.25 | 1.42274012 | 0.23729901 | 0 | 2.04179126 |
| 5 | 256.25 | 1.19206723 | 0.00662612 | 0 | 2.04179126 |

围道积分采用分段 Gauss–Legendre 求积，从每段 8 节点加密至 12 节点的最大变化小于 $5.4\times10^{-5}$。这是数值对照，(6) 对充分大 $X$ 的一致误差来自证明中的零点计数和 Mellin 尾界。复现程序使用 Python、NumPy 和 python-flint，运行时添加参数 --contour 可核对完整特征族。[复现程序](gap-two-logarithmic-zero-height/verify.py)
