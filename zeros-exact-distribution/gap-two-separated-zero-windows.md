---
title: "四个分离高度窗中的零点四次和与孪生素数"
date: "2026-11-25"
abstract: "将陈氏定理与二次特征的零点矩结合。在不损失陈氏下界的数量级下，可把半素数的较小素因子限制到严格位于指数 3/11 与 1/2 之间的区间。用四个独立指定的零点高度窗分辨这些候选。窗宽为双对数时，消去半素数后的矩给出任意固定对数精度。窗宽取对数量级时，对候选集的四次和取最近整数即得精确孪生素数计数。校正后的正下界仍未得到。"
---

# 四个分离高度窗中的零点四次和与孪生素数

固定一个素数 $q$，写 $x=q+2$。前三篇的区别在于它们能告诉我们多少关于**一批**这样的 $q$ 的信息。

| 研究 | 能控制的量 | 对孪生素数下界还缺什么 |
|---|---|---|
| 第 25 篇，全区间线性显式公式 | 高度 $O(\log X)$ 的零点和能检测孤立的一对，另有两项完全幂修正 | 需要给有符号的全区间和一个下界 |
| 第 26 篇，单个素数导子 | 高度 $O(\log\log q)$ 的一次和能逐点辨认素数邻项 | 需要证明越过阈值的导子反复出现 |
| 第 27 篇，陈氏定理与异特征零点对 | 每个大二倍区间有 $\gg X/\log^2 X$ 个导子的分数至少为 $1/2-o(1)$ | 要区分给出 $1/2$ 的半素数和给出 $1$ 的真孪生素数 |

本文选择第三条路线。它已有规模正确的候选集，可以在这个集合上研究一个带符号的矩。下文将这个矩写成四个不同本原 $L$ 函数的零点和，并估计在半素数上消去时的误差。[第 25 篇](https://lixiang90.github.io/vibemath-NumberTheory/zeros-exact-distribution/articles/gap-two-logarithmic-zero-height.html) · [第 26 篇](https://lixiang90.github.io/vibemath-NumberTheory/zeros-exact-distribution/articles/gap-two-single-prime-conductor.html) · [第 27 篇](https://lixiang90.github.io/vibemath-NumberTheory/zeros-exact-distribution/articles/gap-two-chen-zero-pairs.html)

用 $P_2(X)$ 表示 $X<q\le2X$ 中的孪生素数对 $(q,q+2)$ 的个数。用 $B_2(X)$ 表示同一区间中 $q$ 为素数、$q+2=rs$ 且 $x^{1/4}\le r\le s$ 为素数的个数，允许 $r=s$。对充分大的 $X$，这两类 $q$ 都模 $6$ 余 $5$。我们将仅用零点和阈值选出它们的并集，在其中构造两个四次和 $U_X,V_X$。对任意固定 $A>0$，可以使

$$
\frac{U_X}{V_X}-\frac12
=\frac{P_2(X)}{2P_2(X)+B_2(X)}+O_A\big((\log X)^{-A}\big).
\tag{1}
$$

$V_X\asymp X/\log^2 X$ 无条件成立。式 (1) 左端若在无穷多个区间上超过某个固定的正数，便在这些区间推出 $P_2(X)\gg X/\log^2 X$。我们证明的是余项以及分母的上下界。给左端一个独立的正下界，是孪生素数问题留给这条路线的步骤。以下结果均无须假设零点在中心线上。

对于真正妨碍这个下界的半素数，还能收窄因子的范围。存在绝对常数 $\delta>0$，使每个充分大的二倍区间仍有 $\gg X/\log^2 X$ 个素数 $q$，满足 $x=q+2$ 为素数，或 $x=rs$ 为两素数之积、$x^{3/11+\delta}\le r\le x^{1/2-\delta}$。第 5 节给出这个陈氏定理推论的筛法证明，以及对同一个零点矩怎样施加这一限制。

## 1. 保住指定代表元，同时平移零点窗口

以下取充分大的素数 $q\equiv5\pmod6$。记 $L=\log q$，$\ell=\log x$。令 $\eta$ 为模 $3$ 的非主特征，$\omega$ 取四个实特征 $1,\chi_8,\chi_{-4},\chi_{-8}$，其中 $\chi_D(n)=(D/n)$ 为 Kronecker 符号，$1$ 的导子为 $1$。对每个模 $q$ 的奇特征 $\chi_i$，令 $\psi_{i,\omega}=\chi_i\eta\omega$。奇特征共有 $M=(q-1)/2$ 个。四族的本原导子依次为 $3q,24q,12q,24q$。最后一个与第二个虽有同一导子，其模 $8$ 的本原分量不同。四个二进分量均为实特征，因此从任意两个不同族各取一个函数，它们必不相同、也不互为共轭。

取整数 $1\le N\le 3L$。令 $\operatorname{sinc}(v)=\sin(v)/v$，在 $0$ 处取 $1$，并令 $J=\int_{\mathbb R}\exp(-\cosh t)\,dt$。取

$$
\begin{aligned}
K_N(v)&=\left|\frac1N\sum_{j=0}^{N-1}e^{2ijv}\right|^2,\\
W(e^v)&=\frac{K_N(v)\operatorname{sinc}(4v)^{2N}}{J}
           \int_{\mathbb R}e^{-\cosh t-32itv}\,dt,\\
H(s)&=\int_0^\infty W(u)u^{s-1}\,du .
\end{aligned}
\tag{2}
$$

乘上 Fejér 核 $K_N$ 有两个方便之处。它把中心的有效对数宽度由 $N^{-1/2}$ 缩至 $N^{-1}$，而在变换侧只额外使用 $2(N-1)$ 的高度。确切地，令 $U$ 为 $2N$ 个独立的 $[-4,4]$ 均匀变量之和，再独立取整数 $j$，概率为 $(N-|j|)/N^2$，$|j|<N$。置 $V=U+2j$。傅里叶反演给出

$$
H(s)=\frac{2\pi}{32J}\,\mathbb E\exp[-\cos((s+iV)/32)],
\qquad |V|\le 10N-2 .
\tag{3}
$$

该式在所用条带 $|\Re s|\le 5$ 成立。$W(1)=1$，$|W(u)|\le3\min(u^{12},u^{-12})$。当 $|\log u|\ge\log(7/4)$ 时，右侧还可乘以 $4^{-N}$。这些界由 (2) 中积分沿虚方向平移 $3/8$ 和 $|\operatorname{sinc}(4v)|\le1/(4|v|)$ 得到。加上模长为一的因子 $u^{-ih}$，Mellin 变换成为 $H(s-ih)$，而算术代表 $u=1$ 的值仍恰好是 $1$。[第 26 篇第 2、4 节的傅里叶反演与移线](https://lixiang90.github.io/vibemath-NumberTheory/zeros-exact-distribution/articles/gap-two-single-prime-conductor.html)

这里可从 Alpöge–Claude 证明中抽出一条具体信息。原定义 $\mathrm{QH.MeanSq}$ 在每个容许的 $y$ 上对所有 $|\kappa(n)|\le1$ 成立。在原定理的同一个辅助素数集和幂次窗口内，把 $\kappa(n)$ 换成 $\kappa(n)(n/y)^{-ih}$（$n=0$ 任取模长至多 $1$ 的值），可以对公共实参数 $h$ 取上确界而保持 $Cy\log y\,|P|$ 的均方界。量词要求同一次求和使用同一个 $h$。此外，原 Mellin 斜坡给出 $(y/n-1)_+$ 的变换。对固定尺度 $y$ 和 $f(t)=tW(t/y)(t/y)^{-ih}$，分部积分两次得 $\int_0^\infty(t/n-1)_+f''(t)\,dt=W(n/y)(n/y)^{-ih}$。本文在各导子分别使用这样的权，跨族估计由特征正交性另行证明。[原均方条件](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/MeanSquareBounds/PrimeSetTwoBounds.lean#L25) · [原幂次窗口定理](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/MeanSquareBounds/meanSqK.lean#L47) · [原 Mellin 斜坡引理](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/DedekindLFunctions/LandauTheorems/MellinRamp.lean#L9)

给定任意实数 $h$，置

$$
\begin{aligned}
D_h &=2+\log(24q(3+|h|)),&
B_h&=\log(24q)+\log D_h,\\
R_h&=10N-2+32\log(64B_h),&
I_q(h)&=[h-R_h,h+R_h].
\end{aligned}
\tag{4}
$$

对 $1\le y\le 2q$、$a\not\equiv0,\pm1\pmod q$ 定义

$$
\begin{aligned}
S_{i,\omega}(y,h)&=
\sum_{\substack{\rho=\beta+i\gamma,\ L(\rho,\psi_{i,\omega})=0\\
          0<\beta<1,\ \gamma\in I_q(h)}} y^\rho H(\rho-ih),\\
F_\omega(a,y,h)&=-\frac1M\sum_i\overline{\chi_i(a)}\,S_{i,\omega}(y,h).
\end{aligned}
\tag{5}
$$

零点按重数计。只需求位于指定区间的零点，复共轭运算可直接作用于已求得的复数。对每个本原特征，去掉窗口外的零点给 (5) 第一式造成的绝对误差为 $O(q^{-15})$，常数与 $h$ 无关。理由是每个以 $u$ 为起点的单位高度段内有 $O(\log(24q(3+|u|)))$ 个零点，可对 Bennett 等的对称计数公式在相邻高度作差，小高度直接用 $N_\psi(3)$。在第 $j$ 个壳 $2^jR_h<|\gamma-h|\le2^{j+1}R_h$，零点数为 $O((1+2^j R_h)(D_h+\log(2+2^jR_h)))$。由 (3)，每项的 $|H|$ 至多为 $\exp[-0.49(64B_h)^{2^j}]\le\exp[-30\cdot2^jB_h]$。乘 $y^\beta\le2q$ 并求和即可，因为 $R_h\ll B_h$，$D_h\le e^{B_h}$。[Bennett、Martin、O’Bryant、Rechnitzer，定理 1.1](https://arxiv.org/pdf/2005.02989#page=2)

若 $|h|\le\exp(q^b)$，$b>0$ 固定，则 $R_h=10N-2+32\log\log q+O_b(1)$。例如可取 $h_j=j e^q$，$1\le j\le4$。当 $N\asymp\log\log q$ 时，四个正高度窗两两分离，且每个半宽只有 $O(\log\log q)$。更高处的零点数密度会随高度增长，窗宽只度量需要哪一段频谱。

## 2. 直接恢复两个互补因子

记 $\nu_\omega=\eta\omega$，$\Delta_a(n)=\mathbf1_{n\equiv a\ (q)}-\mathbf1_{n\equiv-a\ (q)}$。本原 Dirichlet 函数的显式公式给出

$$
\begin{aligned}
F_\omega(a,y,h)
={}&\sum_n\Lambda(n)\nu_\omega(n)\Delta_a(n)W(n/y)(n/y)^{-ih}\\
 &+\sum_n\frac{\Lambda(n)}n\nu_\omega(n)\Delta_{a^{-1}}(n)
                W(1/(ny))(ny)^{ih}+O(q^{-15}).
\end{aligned}
\tag{6}
$$

在同一族内导子和奇偶性固定，故函数方程中的伽马对数导数、导子项以及跨过的平凡零点贡献对每个 $i$ 相同。因 $\sum_i\overline{\chi_i(a)}=0$，这些项在移线时精确消失，即使 $|h|$ 任意大也成立。奇特征正交性在直接和上给出 $\Delta_a$，在对偶和上给出 $\Delta_{a^{-1}}$。移线可先对固定 $h,q$ 取趋于无穷的水平边，再用 (4) 的一致尾界截断。[本原 Dirichlet 函数方程](https://dlmf.nist.gov/25.15.E5) · [第 26 篇第 2 节的围道推导](https://lixiang90.github.io/vibemath-NumberTheory/zeros-exact-distribution/articles/gap-two-single-prime-conductor.html)

取 $\mathcal D_x=\{d\in\mathbb N:x^{1/4}\le d\le\sqrt x,\ (d,6)=1\}$，并置

$$
b(d)=(2-\mathbf1_{d^2=x})\frac{\ell^2}{4\log d\log(x/d)},\qquad
a_d\equiv2d^{-1}\pmod q,\qquad
g(x)=\frac{\Lambda(x)^2+\sum_{d\in\mathcal D_x,\ d\mid x}b(d)\Lambda(d)\Lambda(x/d)}{\ell^2}.
\tag{7}
$$

所有 $b(d)\le8/3$。在零点侧遍历整个 $\mathcal D_x$，无需先找出 $x$ 的因子。对两族 $\omega,\tau$（以 $\omega\tau$ 表示特征乘积），定义

$$
\begin{aligned}
Z_{\omega,\tau}(h,k)=\ell^{-2}\Big[&
(\omega\tau)(x) F_\omega(2,x,h)\overline{F_\tau(2,x,k)}\\
&+\tau(x)\sum_{d\in\mathcal D_x}b(d)(\omega\tau)(d)
 F_\omega(d,d,h)\overline{F_\tau(a_d,x/d,k)}\Big],\\
C_{\omega,\tau}(h,k)={}&
\frac{(\omega\tau)(5)}{(\log5)^2}F_\omega(5,5,h)\overline{F_\tau(5,5,k)} .
\end{aligned}
\tag{8}
$$

**引理。** 对任意两窗、上述任意两族，一致地有

$$
Z_{\omega,\tau}(h,k)=g(x)+O(4^{-N}+q^{-2}),\qquad
C_{\omega,\tau}(h,k)=1+O(q^{-11}).
\tag{9}
$$

**证明。** 对 $F_\omega(2,x,h)$，$n=x$ 的相位是 $1$，且 $\eta(x)=1$。$q-2$ 是 $3$ 的倍数，被 $\eta$ 消去。$n=2$ 和对偶和共为 $O(q^{-12})$。其余 $n=jq\pm2$ 中 $j\ge2$，其与 $x$ 的比至少 $7/4$ 且 $\gg j$，总贡献 $O(4^{-N}L)$。所以该投影为 $\omega(x)\Lambda(x)+O(4^{-N}L+q^{-12})$。

对 $F_\omega(d,d,h)$，主项是 $\nu_\omega(d)\Lambda(d)$，其余直接项为 $O((d/q)^{12}L)$，对偶项为 $O(d^{-12})$，合计 $O(q^{-3})$。另一因子 $F_\tau(a_d,x/d,k)\ll L$，因为两个剩余类中低于 $q$ 的代表至多两个，较远代表按十二次幂衰减。将前一因子替成主项的总代价为 $O(\sqrt x\,q^{-3}L)$。

后一因子的直接代表满足 $dn=jq\pm2$。$d\ge5$ 排除 $j=0$，$j=1$ 的负号被 $\eta(n)$ 消去。正号正好给出 $dn=x$，此时
$(\omega\tau)(d)\tau(x)\nu_\omega(d)\nu_\tau(n)=1$，且两个测试权所取的坐标都是 $1$。对 $j\ge2$ 先按乘积合并，用 $\sum_{d\mid m}\Lambda(d)\Lambda(m/d)\le(\log m)^2$，全部远处代表总计 $O(4^{-N}\sum_{j\ge2}j^{-12}\log^2(jq+2))=O(4^{-N}L^2)$。对偶项每个为 $O((x/d)^{-12})$。除以 $\ell^2$ 即得 $Z$ 的复数绝对误差界。

最后，$F_\omega(5,5,h)=\nu_\omega(5)\log5+O(q^{-11})$。直接和的其余代表为 $jq\pm5$，对偶和的代表满足 $5n=jq\pm1$，$j\ge1$。两者均可逐项用十二次幂界求和。代回 (8) 即得第二式。取小素数 $5$ 的作用是用同一批零点构造值为 $1+o(1)$ 的二次量。

## 3. 四个族的有符号矩

给定四个中心 $\boldsymbol h=(h_1,h_2,h_3,h_4)$，简记

$$
\begin{aligned}
Z_{12}&=Z_{1,\chi_8}(h_1,h_2),&
Z_{34}&=Z_{\chi_{-4},\chi_{-8}}(h_3,h_4),&
C_{34}&=C_{\chi_{-4},\chi_{-8}}(h_3,h_4),\\
u_q&=\Re(Z_{12}Z_{34}),&
v_q&=\Re(Z_{12}C_{34}),&
t_q&=2u_q-v_q .
\end{aligned}
\tag{10}
$$

三个量展开后均恰为四次，每一项在各指定窗中各取一个零点。各族的二进分量排除了同函数或共轭函数，即便个别函数有相同数值的零点纵坐标也不影响这一条件。由 (9)，写 $e_q=4^{-N}+q^{-2}$，有

$$
u_q=g(x)^2+O(e_q),\qquad v_q=g(x)+O(e_q),\qquad
t_q=2g(x)^2-g(x)+O(e_q).
\tag{11}
$$

$g$ 的取值决定了这三个量的意义。素数给 $g=1$。两个均至少为 $x^{1/4}$ 的素数之积给 $g=1/2$，平方的两项各给 $1/4$。其余 $x$ 有 $0\le g\le13/36$。简证如下。两个不同素数底数只容许两个完整的素数幂配对，在 $x=p^ar^b$、$ab\ge2$ 时最多给 $1/(2ab)$。有三个不同底数的整数给零。对 $p^j$、$j\ge3$，不限制因子范围的上界是 $j^{-2}+H_{j-1}/(2j)\le 1/9+1/4$，其中 $H_k=\sum_{i\le k}1/i$。两个不同素数中较小者低于 $x^{1/4}$ 时也给零。[第 27 篇第 4 节的取值表](https://lixiang90.github.io/vibemath-NumberTheory/zeros-exact-distribution/articles/gap-two-chen-zero-pairs.html)

**定理。** 令 $\mathcal Q(X)$ 为 $X<q\le2X$ 且 $q\equiv5\pmod6$ 的素数。每个 $q$ 可各自选择四个中心。取 $N\to\infty$ 且 $N\le3L$，定义

$$
\mathcal C(X)=\{q\in\mathcal Q(X):\Re Z_{12}>5/12\},\quad
U_X=\sum_{q\in\mathcal C(X)}u_q,\qquad V_X=\sum_{q\in\mathcal C(X)}v_q .
\tag{12}
$$

充分大时，$\mathcal C(X)$ 恰为 $P_2(X)$ 和 $B_2(X)$ 计数的两类素数的并集。对其中**同一批** $\gg X/\log^2X$ 个素数，无论四个中心如何选择，均有 $u_q\ge1/4-O(e_q)$。写 $E_X=\max_{q\in\mathcal Q(X)}e_q$。则

$$
\begin{aligned}
U_X&=P_2(X)+\tfrac14 B_2(X)+O(E_X|\mathcal C(X)|),\\
V_X&=P_2(X)+\tfrac12 B_2(X)+O(E_X|\mathcal C(X)|),\\
2U_X-V_X&=P_2(X)+O(E_X|\mathcal C(X)|).
\end{aligned}
\tag{13}
$$

**证明。** $13/36<5/12<1/2$ 给候选集两侧留下固定间隙。(9) 的误差一致，故当 $q$ 足够大时，无论怎样逐个选窗，候选集相同。Green–Tao 引用的定量陈氏定理保证 $(X,2X]$ 中 $\gg X/\log^2X$ 个素数 $q$ 的邻项是素数或两个均超过 $q^{3/11}$ 的素数之积。因 $q^{3/11}>x^{1/4}$，这些 $q$ 全在 $\mathcal C(X)$。另一方面，候选要求 $q(q+2)$ 没有低于 $X^{1/8}$ 的素因子。对奇素数筛去 $0,-2$ 两个剩余类，维数二上界筛给 $|\mathcal C(X)|\ll X/\log^2X$。具体可在算术大筛中用水平 $X^{1/3}$，只取 $3\le p\le X^{1/24}$。乘积 $\prod_p(1-2/p)^{-1}\gg\log^2X$，按展开项归一化后 $\log d$ 的期望是 $2\sum_p(\log p)/p=(1/12+o(1))\log X$，故 $d\le X^{1/3}$ 的分母仍保留至少一半。对候选集在 (11) 中逐项求和即得 (13)。[Green–Tao，定理 6.1 与脚注 3](https://www.numdam.org/article/JTNB_2006__18_1_147_0.pdf#page=20) · [Montgomery–Vaughan，推论 1 的算术大筛](https://personal.science.psu.edu/rcv4/personal/Publications/large_sieve.pdf#page=2)

取 $N=\max(1,\lceil A\log L/\log4\rceil)$，便有 $E_X\ll_A(\log X)^{-A}$。又 $V_X\asymp |\mathcal C(X)|\asymp X/\log^2X$，所以 (13) 相除给出 (1)。在 $|h_j|\le\exp(q^b)$ 的范围，这只须半宽 $O_{A,b}(\log\log q)$。相对于第 27 篇在所有素数导子上求和，此处先筛掉了零权与高次幂的小权邻项，且四个族的局部特征自动区分每个零点的来源。式 (10) 是齐次四次和，式 (12) 另用一次阈值选择候选导子。

还可以按两个尺度使用这个结论。

- **检测稀疏的一对。** 对任意上述 $N\to\infty$，充分大时 $t_q>1/2$ 当且仅当 $q+2$ 为素数。对合数有 $g\le1/2$，故 $2g^2-g\le0$。孪生素数无穷性等价于这个四窗不等式在无穷多个素数导子上成立，即使预先指定 $h_j=j e^q$ 也成立。
- **精确计数。** 取 $N=\lceil2L/\log4\rceil$，则 $E_X\ll X^{-2}$。对上述高度范围，半宽为 $O_b(\log q)$，而 (13) 最后一式的绝对误差为 $O(1/(X\log^2X))$。所以对每个充分大的 $X$，$2U_X-V_X$ 的最近整数正好是 $P_2(X)$。这个结论也适用于只有一对或没有孪生素数的区间。它使用多个导子中全部带标签的奇特征。

陈氏定理给出的 $u_q\ge1/4-o(1)$ 可以全部由半素数贡献。在这类导子上 $2u_q-v_q$ 趋近零。由于 (9) 对四个中心一致，单靠移动取样高度，不能从半素数产生一个固定的正盈余。若要得到新的素数下界，需要从特征相位和零点之和的符号上独立估计 (13) 左端。式 (10) 保留了相位、尺度、导子和高度标签，经典的单个 $\zeta$ 函数的两两间距统计没有给出这一四次矩。[Montgomery 的对关联原文](https://public.websites.umich.edu/~hlm/paircor1.pdf) · [Ford–Zaharescu 对不同函数零点纵坐标差的研究](https://arxiv.org/abs/1305.2520)

## 4. 只用两个导子时，如何删除同特征配对

也可以要求前两个零点同属导子 $3q$，后两个同属导子 $12q$，而四个函数仍两两不同且非共轭。这时需要控制被删除的特征对。取 $N=\lceil4\log L/\log4\rceil$。在一个固定的族内，记 $w_i(a)=\bar\chi_i(a)$、$S_i=S_{i,\omega}$、$[i]=\{\chi_i,\bar\chi_i\}$，定义

$$
\begin{aligned}
A_{ij}(a,y,h)&=-\tfrac12\big(w_i(a)-w_j(a)\big)\big(S_i(y,h)-S_j(y,h)\big),\\
\mathcal P_\omega(a,y,h;b,z,k)&=
\frac1{M^4}\sum_{\substack{i,j,r,s\\([i]\cup[j])\cap([r]\cup[s])=\varnothing}}
 A_{ij}(a,y,h)\overline{A_{rs}(b,z,k)} .
\end{aligned}
\tag{14}
$$

$M^{-2}\sum_{i,j} A_{ij}=F_\omega$。差分在每对内消去公共伽马项。式 (14) 展开时每个单项只用两个不同且非共轭的函数，不会因减去全族均值而重新引入被删零点。记 $\sigma(y,h)^2=M^{-1}\sum_i|S_i(y,h)-M^{-1}\sum_j S_j(y,h)|^2$。有

$$
\mathcal P_\omega(a,y,h;b,z,k)
=F_\omega(a,y,h)\overline{F_\omega(b,z,k)}
 + O\!\left(\frac{\sigma(y,h)\sigma(z,k)}M\right).
\tag{15}
$$

事实上，固定 $(i,j)$ 至多禁掉 $8M$ 个 $(r,s)$，反之亦然。又 $\sum_{ij}|S_i-S_j|^2=2M^2\sigma^2$，且 $|w_i-w_j|\le2$。在被禁的四指标集合上作一次 Cauchy 估计，再除以 $M^4$ 即得 (15)。

对所有高度有

$$
\sigma(y,h)^2\ll yL/N\quad(x^{1/4}\le y\le2q),\qquad \sigma(5,h)^2\ll1 .
\tag{16}
$$

对单特征应用 (6) 推导中的移线，把公共项保留到求方差时精确消去。对偶和的绝对值为 $O(y^{-12})$。直接和先限制在 $y/2\le n\le2y$，用全部特征的均方控制奇特征的均方，代价至多两倍。每个剩余类至多有五个代表，正交性将其上界化为常数倍的 $\sum_n\Lambda(n)^2|W(n/y)|^2$。在 $|v|\le\log2$，由 $K_N(v)=(\sin(Nv)/(N\sin v))^2$ 及 $K_N\le1$，有 $K_N(v)\ll(1+N|v|)^{-2}$。将 $v$ 按 $1/N$ 分段，每段所含素数由 Brun–Titchmarsh 界控制为 $O(y/(N\log y))$，因为可将每段加宽到 $\asymp y/N$ 且 $\log(y/N)\asymp\log y$。素数的平方权至多 $O(L^2)$，按 $(1+|j|)^{-4}$ 求和得 $O(yL/N)$。所有高次素数幂另贡献 $O(\sqrt y\log^3(2y))$，也在此界内。[Montgomery–Vaughan，定理 2 与式 (1.12)](https://personal.science.psu.edu/rcv4/personal/Publications/large_sieve.pdf#page=3)

在中心段外，对 $n\asymp u$ 的二倍段，每个剩余类有 $O(1+u/q)$ 个代表，且 $\sum_{n\le2u}\Lambda(n)^2\ll u\log(2u)$。在均方根范数中逐段用 (2) 的十二次幂和 $4^{-N}$ 求和，得到 $O(4^{-N}\sqrt{yL})$。加上 (5) 的截断误差即得 (16)。固定 $y=5$ 时直接和与对偶和的绝对和均有绝对常数上界。

在 $Z_{\omega,\omega}$ 和 $C_{\omega,\omega}$ 中将同族的每个乘积 $F_\omega\overline{F_\omega}$ 换成 (14)，记为 $\widetilde Z_\omega,\widetilde C_\omega$。对前一式的主项用 $y=z=x$，对因子项用 $yz=x$ 并累加至多 $\sqrt x$ 个有界系数。(15)、(16) 给出

$$
\widetilde Z_\omega(h,k)=g(x)+O((LN)^{-1}),\qquad
\widetilde C_\omega(h,k)=1+O(q^{-1}).
\tag{17}
$$

故以 $\widetilde Z_1(h_1,h_2),\widetilde Z_{\chi_{-4}}(h_3,h_4),\widetilde C_{\chi_{-4}}(h_3,h_4)$ 代入 (10) 至 (13) 均成立，误差参数改为 $E_X\ll1/(\log X\log\log X)$。在这个限制下，Fejér 核使逐导子的删除误差从第 27 篇的 $O(1/(L\sqrt{\log L}))$ 降至 $O(1/(L\log L))$，并允许两端使用分离的任意高度。第 3 节使用四个局部特征，进一步消去了这项删除代价。

## 5. 把剩下的半素数限制到内部因子区间

令 $\theta=\log r/\log(q+2)$，其中 $q+2=rs$，$r\le s$ 均为素数。对于 $J\subset[1/4,1/2]$，用 $B_J(X)$ 计数 $\theta\in J$ 的这类 $q$。若 $J$ 是 $m$ 个区间的并，总长为 $|J|$，则有绝对常数使

$$
B_J(X)\ll \left(|J|+\frac{m}{\log X}\right)\frac{X}{\log^2 X}.
\tag{18}
$$

**证明。** 先取 $J=[a,b]$。$X^a\le r\le (2X+2)^b$，固定这样的素数 $r$，置 $Y=X/r$。对可能的 $s$，$s$ 与 $rs-2$ 须同时为素数，且 $s$ 在长 $O(Y)$ 的区间内。对每个 $3\le p\le Y^{1/24}$ 的素数筛去 $s\equiv0,2r^{-1}\pmod p$。因为 $r\ge X^{1/4}$ 且两个待筛素数均大于 $Y^{1/24}$，始终恰有两个禁类。第 3 节的大筛证明以 $Y^{1/3}$ 为水平，给出 $O(Y/\log^2Y)$，常数与 $r$ 无关。此处 $\log Y\asymp\log X$。由 $\pi(t)\ll t/\log t$ 及部分求和，

$$
\sum_{\substack{X^a\le r\le(2X+2)^b\\r\ {\rm prime}}}\frac1r
\ll \frac1{\log X}+
\log\!\left(\frac{b\log(2X+2)}{a\log X}\right)
\ll b-a+\frac1{\log X}.
$$

求和即得 (18)，再对各分区间相加。这里也允许 $a=b$，包括邻项为平方的端点。

**推论。** 设 $\alpha=3/11$。存在绝对的 $c,\delta>0$，使 $J_\delta=[\alpha+\delta,\,1/2-\delta]$ 满足

$$
P_2(X)+B_{J_\delta}(X)\ \ge\ c\,X/\log^2X
\qquad (X\ \text{充分大}).
\tag{19}
$$

具体地，Green–Tao 定理 6.1 给出常数 $c_0>0$。其中半素数的 $\theta$ 大于 $\alpha\log q/\log(q+2)\ge\alpha-O(1/\log X)$。用 (18) 删去 $[\alpha-O(1/\log X),\alpha+\delta]$ 和 $[1/2-\delta,1/2]$ 的半素数，所失至多为 $K(2\delta+O(1/\log X))X/\log^2X$，$K$ 为绝对常数。先取 $\delta$ 使 $2K\delta<c_0/4$ 且 $2\delta<1/2-\alpha$，再取 $X$ 充分大即可。这里从已发表的陈氏下界推得存在性，未给 $\delta$ 指定数值或最优值。[Green–Tao，定理 6.1](https://www.numdam.org/article/JTNB_2006__18_1_147_0.pdf#page=20) · [Montgomery–Vaughan，推论 1](https://personal.science.psu.edu/rcv4/personal/Publications/large_sieve.pdf#page=2)

在 (7)、(8) 中都只保留 $\log d/\ell\in J_\delta$ 的因子，记所得算术量为 $g_J$。第 2 节的误差使用绝对上界，删除这些项后同样成立。素数仍给 $g_J=1$，在内部区间的半素数给 $g_J=1/2$。其余两个不同素数之积给零，平方只余 $\Lambda(x)^2/\ell^2=1/4$，其余合数仍有 $g_J\le13/36$。因此零点阈值 $5/12$ 恰好选出 (19) 的两类。式 (1) 和 (13) 中的 $B_2(X)$ 可全部替为 $B_{J_\delta}(X)$，且 $V_X$ 保持原数量级。这样要研究的正盈余只剩内部区间的贡献。若该区间里的半素数承担了全部陈氏质量，上述推导仍然容许 $P_2(X)=0$。

## 6. 数值复核

取 $h_1,h_2,h_3,h_4=320,960,1600,2240$，$N=\lceil4\log L/\log4\rceil$。对表中每个 $q$，四个窗连同 (4) 的全部半宽都位于正半轴、互不相交。下表列出由显式公式围道侧算出的分数。两族列按 (14) 删除配对，四族列按 (8) 选局部特征。最后一列是全部族、尺度和高度中，去掉族平均后的单特征零点和与独立算术求值的最大绝对差。

| $q$ | $q+2$ | $g(q+2)$ | 两族的 $t_q$ | 四族的 $t_q$ | 最大复数差 |
|---:|---|---:|---:|---:|---:|
| 17 | $19$ | 1 | 0.0742264 | 1.0000000 | $8.06\cdot10^{-12}$ |
| 23 | $25=5^2$ | 0.5 | $-0.0000013$ | 0.0000000 | $3.39\cdot10^{-11}$ |
| 29 | $31$ | 1 | 0.2673772 | 1.0000000 | $6.67\cdot10^{-10}$ |
| 53 | $55=5\cdot11$ | 0.5 | $-0.0011252$ | 0.0000000 | $4.05\cdot10^{-9}$ |
| 173 | $175=5^2\cdot7$ | 0.25 | $-0.1038974$ | $-0.1250000$ | $8.91\cdot10^{-7}$ |

这些小导子展示了删除特征对的代价。$q=17$ 时两族矩仍远低于其渐近值 $1$。在 $(16,32]$ 中，候选阈值保留 $17,23,29$，四族有符号和取最近整数得 $2$，对应 $(17,19)$ 和 $(29,31)$。

第 5 节的因子截取另用整数算术核对。下表用 $\delta=1/132$ 作有限样本，逐个用 $x^{37}\le r^{132}\le x^{65}$ 判定，并核对截取后的 $g_J>5/12$。样本参数的取值不用于 (19) 的存在性证明。

| 区间 $(X,2X]$ 的 $X$ | 真孪生素数对 | 较小因子至少为 $x^{1/4}$ 的半素数候选 | 内部区间留下的半素数候选 |
|---:|---:|---:|---:|
| 1 000 | 26 | 20 | 14 |
| 10 000 | 137 | 146 | 123 |
| 100 000 | 936 | 894 | 780 |

程序在 $\Re s=4.5$ 直接用 Dirichlet 级数求 $L,L'$，按剩余类的离散傅里叶变换得到各特征。截取 $K=q\lceil300(\max y)^{9/7}/q\rceil$ 项，$L,L'$ 的绝对尾误差为 $O(K^{-3.5}\log K)$，且 $|1/L|\le\zeta(4.5)$。同族减均值后，左线 $1-4.5$ 上的对数导数可由右线负共轭求得，已跨过的平凡零点按族相消。对两条无限竖线作数值积分，尾部由 (3) 衰减，在窗外省去的零点另有第 1 节的界。计算使用两层梯形网格，另以素数幂直接和、对偶和复算，并在每个窗心用 Arb 核对各本原函数的 $L'/L$。对 $q=17$ 还用四指标逐项枚举核对 (14)。这些是浮点复核，表中数据没有逐个定位零点。误差阶和候选下界的证明分别见第 2 至 4 节。[复算程序](https://github.com/lixiang90/vibemath-NumberTheory/blob/main/zeros-exact-distribution/gap-two-separated-zero-windows/verify.py)
