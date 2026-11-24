---
title: "固定差二的异特征零点对与陈氏定理的半阈值"
date: "2026-11-25"
abstract: "在单个导子 3q 上配对不同且非共轭的本原 Dirichlet 特征，只取双对数高度的零点，构造一个与固定差 2 有关的有符号二次和。归一化后，真孪生素数的值趋近 1，两个足够大的素数之积给出 1/2。陈氏定理保证每个充分大的二倍区间有正确数量级的导子达到 1/2−o(1)。上极限严格超过 1/2 等价于孪生素数无穷性。"
---

# 固定差二的异特征零点对与陈氏定理的半阈值

陈氏定理给出了许多素数 $q$，使 $q+2$ 为素数或两个大素数之积。能否把这一无条件下界读成**来自不同 $L$ 函数的零点对**共同产生的大值，并看清哪些大值对应真的孪生素数？

本文给出一个这样的统计量 $V_q$。对素数 $q\equiv5\pmod6$，它只取导子 $3q$、高度 $O(\log\log q)$ 内的非平凡零点。求和时去掉同一特征和共轭特征的所有配对。对每个充分大的 $X$，区间 $(X,2X]$ 中至少有 $cX/\log^2 X$ 个素数 $q$ 满足 $V_q\ge1/2-o(1)$。无条件地，

$$
\limsup_{\substack{q\to\infty,\ q\ {\rm prime}\\q\equiv5\pmod6}} V_q
=\begin{cases}
1,&\text{孪生素数有无穷多对},\\
1/2,&\text{孪生素数只有有限多对}.
\end{cases}
\tag{1}
$$

第 26 篇用一次零点和判别孪生素数。这里把观察点改为较小的素数端点 $q$，用二次和同时测量 $q+2$ 的素数权与互补因子的素数权。陈氏定理提供严格正下界，解析上的工作是证明这个下界在删去同特征及共轭配对后仍然存在。式 (1) 中严格超过 $1/2$ 的部分仍未证明。[前篇的单导子判据](https://lixiang90.github.io/vibemath-NumberTheory/zeros-exact-distribution/articles/gap-two-single-prime-conductor.html) · [Green–Tao，第 6 节所用的定量陈氏定理](https://www.numdam.org/article/JTNB_2006__18_1_147_0.pdf#page=20)

## 1. 二次和的定义

写 $x=q+2$、$L=\log q$、$\ell=\log x$，以下取充分大的 $q\equiv5\pmod6$。令 $\eta$ 为模 $3$ 的非平凡特征，$\mathcal C_q$ 为模 $q$ 的全部奇特征，$M=(q-1)/2$。各 $\chi\eta$ 恰为导子 $3q$ 的全部本原偶特征。取

$$
N=\left\lceil4\log L/\log4\right\rceil,\qquad
T=32\log(32\log(3q))+8N.
\tag{2}
$$

令 $J=\int_{\mathbb R}e^{-\cosh t}\,dt$，$\operatorname{sinc}v=\sin(v)/v$，在零点取 $1$。沿用前篇的实权及其 Mellin 变换

$$
W(e^v)=\frac{\operatorname{sinc}(4v)^{2N}}{J}
\int_{\mathbb R}e^{-\cosh t}e^{-32itv}\,dt,
\qquad H(s)=\int_0^\infty W(u)u^{s-1}\,du .
\tag{3}
$$

对 $\chi\in\mathcal C_q$ 定义

$$
S_\chi(y)=\sum_{\substack{L(\rho,\chi\eta)=0,\ 0<\Re\rho<1\\|\Im\rho|\le T}}
 y^\rho H(\rho),\qquad
f_\chi(a,y)=-\overline{\chi(a)}S_\chi(y)/M,\qquad F(a,y)=\sum_{\chi\in\mathcal C_q}f_\chi(a,y).
\tag{4}
$$

两个方向的高度均按重数计算。令 $\mathcal D_x$ 含所有满足 $x^{1/4}\le d\le\sqrt x$、$(d,6)=1$ 的整数。对这些整数定义

$$
b_x(d)=\bigl(2-\mathbf1_{d^2=x}\bigr)\frac{\ell^2}{4\log d\log(x/d)},
\qquad a_d\equiv2d^{-1}\pmod q .
\tag{5}
$$

$b_x(d)\le8/3$。在 (5) 中遍历整个整数区间即可，系数不依赖 $d$ 是否为素数或是否整除 $x$。所用统计量为

$$
\begin{aligned}
K_q={}&
\sum_{\substack{\chi,\lambda\in\mathcal C_q\\\lambda\notin\{\chi,\bar\chi\}}}
\left[f_\chi(2,x)f_\lambda(2,x)
 +\sum_{d\in\mathcal D_x}b_x(d)f_\chi(d,d)f_\lambda(a_d,x/d)\right],\\
V_q={}&K_q/\ell^2 .
\end{aligned}
\tag{6}
$$

它是实数，因为共轭置换保持求和范围且将每项共轭。每一项均来自两种不同且非共轭的本原 $L$ 函数。比如因子项中，一对零点 $\rho,\rho'$ 的算术权含有

$$
\overline{\lambda(2)}\,x^{\rho'}
 \sum_{d\in\mathcal D_x}b_x(d)(\bar\chi\lambda)(d)d^{\rho-\rho'} .
\tag{7}
$$

在 Dirichlet GRH 下，最后一个幂化为 $e^{i(\gamma-\gamma')\log d}$。同时仍有 $x^{i\gamma'}$ 和特征标签。Montgomery 的经典对关联研究同一个 $\zeta$ 函数的纵坐标差，式 (7) 说明固定差 $2$ 在这里还要求保留哪些相位。本文所有估计均允许零点偏离中心线。[Montgomery，原始对关联论文](https://public.websites.umich.edu/~hlm/paircor1.pdf)

## 2. 不删特征时的算术量

把 (6) 的特征对先取全，得到

$$
K_q^{\rm all}=F(2,x)^2+\sum_{d\in\mathcal D_x}b_x(d) F(d,d)F(a_d,x/d).
\tag{8}
$$

以 $\Lambda$ 表示 von Mangoldt 函数，定义

$$
G(x)=\Lambda(x)^2+\sum_{\substack{d\in\mathcal D_x\\d\mid x}}
b_x(d)\Lambda(d)\Lambda(x/d).
\tag{9}
$$

**命题 1.** 一致地有 $K_q^{\rm all}=G(x)+O(L^{-2})$。

**证明。** 对 $a\not\equiv\pm1\pmod q$，令 $\Delta_a(n)=\mathbf1_{n\equiv a\ (q)}-\mathbf1_{n\equiv-a\ (q)}$。在 $1\le y\le2q$ 上，前篇的显式公式给出

$$
F(a,y)=\sum_n\Lambda(n)\eta(n)\Delta_a(n)W(n/y)
 +\sum_n\frac{\Lambda(n)}{n}\eta(n)\Delta_{a^{-1}}(n)W(1/(ny))+O(q^{-10}).
\tag{10}
$$

这里每个奇特征的系数为 $-\bar\chi(a)/M$。其和为零，所以导子项、公共的伽马因子以及偶特征在 $s=0$ 的平凡零点同时消失。式 (3) 的 Mellin 变换等于 $(2\pi/(32J))\,\mathbb E\exp[-\cos((s+iU)/32)]$，其中 $U$ 是 $2N$ 个独立的 $[-4,4]$ 均匀变量之和。故 $|U|\le8N$，从 (2) 的高度起仍有双指数衰减。对每个特征按二倍高度段使用 $N_{\chi\eta}(t)\ll t\log(3qt)$，并用 $y^{\Re\rho}\le2q$，截断误差为 $O(q^{-10})$，对单个特征同样成立。[前篇第 3、4、6 节](https://lixiang90.github.io/vibemath-NumberTheory/zeros-exact-distribution/articles/gap-two-single-prime-conductor.html) · [Bennett 等，定理 1.1](https://arxiv.org/abs/2005.02989)

有 $W(1)=1$ 和 $|W(u)|\le3\min(u^{12},u^{-12})$。当 $|\log u|\ge\log(7/4)$ 时，$1/(4|\log u|)<1/2$，后一个界还可乘以 $4^{-N}$。先对 $F(2,x)$ 用 (10)。$n=x$ 贡献 $\Lambda(x)$，因为 $\eta(x)=1$。$n=q-2$ 被 $\eta$ 消去。$n=2$ 和对偶项均为 $O(q^{-12})$。其余代表 $kq\pm2$ 满足 $k\ge2$，与 $x$ 的比值至少为 $7/4$ 且 $\gg k$，合计 $O(4^{-N}L)$。所以 $F(2,x)=\Lambda(x)+O(L^{-3})$。

对 $F(d,d)$，除 $n=d$ 之外，所有直接代表都形如 $kq\pm d$、$k\ge1$。其贡献为 $O((d/q)^{12}L)$，对偶项为 $O(d^{-12})$。因为 $x^{1/4}\le d\le\sqrt x$，得到 $F(d,d)=\eta(d)\Lambda(d)+O(q^{-3})$。

对 $F(a_d,x/d)$，直接代表恰满足 $dn=kq\pm2$。$d\ge5$ 排除 $k=0$，$d\not\equiv\pm2\pmod q$ 保证 $a_d\not\equiv\pm1$。在 $k=1$ 时，负号给出的 $q-2$ 被 $3$ 整除，而 $3\nmid d$，故 $\eta(n)=0$。正号恰在 $d\mid x$ 时给出 $n=x/d$。其对偶项为 $O((x/d)^{-12})$。这个投影的总大小为 $O(L)$，故先替换所有 $F(d,d)$ 的误差至多为 $O(\sqrt x\,q^{-3}L)$。

对剩下的 $k\ge2$，应先按乘积 $dn=kq\pm2$ 合并。对每个整数 $h$，

$$
\sum_{d\mid h}\Lambda(d)\Lambda(h/d)\le(\log h)^2
\tag{11}
$$

由 $\sum_{d\mid h}\Lambda(d)=\log h$ 即得。由于 $n/(x/d)=(kq\pm2)/x$，全部远处代表的总贡献至多为 $O(4^{-N}\sum_{k\ge2} k^{-12}\log^2(kq+2))=O(L^{-2})$。在 $dn=x$ 处又有 $\eta(d)\eta(n)=1$。这证明了命题 1。按每个 $d$ 分开取绝对误差再乘 $\sqrt x$，会丢掉 (11) 提供的精度。

## 3. 从特征均方中节省一个趋于无穷的因子

Alpöge–Claude 归档中的 Mellin 斜坡引理给出 $(y/n-1)_+$ 的变换。两次分部积分可用这些斜坡表示 (3) 的权。归档的 $\mathrm{QH.MeanSq}$ 定义了辅助素数上的 Jacobi 和均方条件。定理 $\mathrm{meanSqK}$ 在其容许的辅助素数集和幂次窗口上给出 $Cy\log y$ 乘以族大小的上界。本节改在固定模数内平均全部奇特征，保留窗口宽度得到下式。它适用于 (4) 定义的零点和。[原 Mellin 引理](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/DedekindLFunctions/LandauTheorems/MellinRamp.lean#L9) · [原 Jacobi 和的均方条件](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/MeanSquareBounds/PrimeSetTwoBounds.lean#L25) · [原幂次窗口的均方定理](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/MeanSquareBounds/meanSqK.lean#L47)

**命题 2.** 当 $x^{1/4}\le y\le2q$ 时，一致地有

$$
\frac1M\sum_{\chi\in\mathcal C_q}|S_\chi(y)|^2
\ll \frac{yL}{\sqrt N}.
\tag{12}
$$

先在单个 $L(s,\chi\eta)$ 的显式公式中保留公共项。其直接算术和为 $-\sum_n\Lambda(n)\chi(n)\eta(n)W(n/y)$，对偶和为 $O(y^{-12})$。在 $\Re s=-1/2$ 上，函数方程余下的核为

$$
C(s)=\log(3q/\pi)+\tfrac12\bigl(\Psi(s/2)+\Psi((1-s)/2)\bigr),
\qquad \Psi=\Gamma'/\Gamma .
\tag{13}
$$

公共项等于 $(2\pi i)^{-1}\int_{(-1/2)} C(s)y^s H(s)\,ds-H(0)$。由 (10) 证明中的平移表示，$\int_{\mathbb R}|H(-1/2+it)|\,dt\ll1$，带 $\log(2+|t|)$ 的相同积分为 $O(\log(2+N))$。结合伽马对数导数界，公共项为 $O(1+L/\sqrt y)$，在 (12) 的尺度内。[本原 Dirichlet 函数的函数方程](https://dlmf.nist.gov/25.15.E5)

把直接和平方后在所有模 $q$ 特征上取平均，可用两倍这个平均控制奇特征平均。对 $y/2\le n\le2y$，每个剩余类至多出现五次，正交性与 Cauchy 不等式给出上界 $O(\sum\Lambda(n)^2|W(n/y)|^2)$。在 $|v|\le\log2$ 上存在绝对常数 $c>0$ 使 $|\operatorname{sinc}(4v)|^{2N}\le e^{-cNv^2}$，这是零附近的泰勒展开与其余紧区间上的严格小于 $1$ 共同给出的。把对数坐标按 $1/\sqrt N$ 等分。每段所含素数数目 $\ll (y/\sqrt N)/\log y$，因为段长至多 $O(y/\sqrt N)$，可扩大到这一长度并用 Brun–Titchmarsh。每个素数的平方权为 $O(L^2)$，且 $\log y\asymp L$。按上述高斯界求和给 $O(yL/\sqrt N)$。中心段的所有高次素数幂另外至多贡献 $O(\sqrt y\log^3(2y))$，也在此界内。[Montgomery–Vaughan，定理 2 与式 (1.12)](https://personal.science.psu.edu/rcv4/personal/Publications/large_sieve.pdf#page=3)

在中心段以外，把 $n$ 按 $n\asymp 2^j y$ 分段。每个剩余类的占用数为 $O(1+2^j y/q)$。用 $\sum_{n\le u}\Lambda(n)^2\ll u\log(2u)$、(3) 的幂衰减及 $4^{-N}$ 因子，在均方根范数中对各段求和，得到 $O(4^{-N}\sqrt{yL})$。权中的十二次幂足以吸收占用数，使相应的几何级数收敛。连同 (10) 的单特征截断误差，命题 2 得证。

现在删除 $\lambda=\chi$ 和 $\lambda=\bar\chi$ 的全部配对，这两个关系各是一个置换。对任意两个尺度 $y,z$，由 (12) 和 Cauchy 不等式，删去项的绝对值总和在除以 $M^2$ 后至多为 $O(\sqrt{yz}\,L/(q\sqrt N))$，这里每项可带任意模长至多 $1$ 的相位。在 (6) 第一项用 $y=z=x$。在因子项用 $yz=x$，然后求和至多 $\sqrt x$ 个有界权 $b_x(d)$。于是

$$
|K_q-K_q^{\rm all}|\ll L/\sqrt N,\qquad
V_q=\frac{G(x)}{\ell^2}+O\!\left(\frac{1}{L\sqrt{\log L}}\right).
\tag{14}
$$

这一步解释了 (3) 中增长的 $N$ 的作用。只得到 $O(L)$ 的逐导子删除误差，在 $(X,2X]$ 的素数上累加可达 $O(X)$，与所求下界同阶。式 (14) 将该误差降到了 $o(X)$。

还可预先排除更多特征配对。若删去集合在每一行和每一列均至多含 $r_q$ 项，对删去集合应用 Cauchy 不等式，每一侧的平方和至多重复 $r_q$ 次，故删除代价为 $O(r_q L/\sqrt N)$。只要集合在同时取共轭下封闭，所得和仍为实数。即使 $r_q$ 随导子增长，只要 $r_q=o(\sqrt{\log\log q})$，下节的有符号平均下界和 (1) 仍保持。

## 4. 半素数给出的精确阈值

算术量 $G(x)/\ell^2$ 的取值可逐项分类。

| $x=q+2$ 的分解 | $G(x)/\ell^2$ |
|---|---:|
| 素数 | $1$ |
| $rs$，$r,s$ 为素数且 $x^{1/4}\le r\le s$，允许 $r=s$ | $1/2$ |
| $rs$，两素数不同且较小者低于 $x^{1/4}$ | $0$ |
| $p^a r^b$，$p\ne r$，$ab\ge2$ | 两个幂均至少为 $x^{1/4}$ 时为 $1/(2ab)$，否则为 $0$ |
| $p^j$，$j\ge3$ | 至多 $13/36$ |
| 至少三个不同素因子 | $0$ |

对不同的两个素数底数，只有两个完整的素数幂才能在 (9) 的互补因子中配对，(5) 恰好抵消这两个因子的对数。对 $x=p^j$，不限制中央因子范围时，上界为 $j^{-2}+\tfrac14\sum_{i=1}^{j-1}1/(i(j-i))=j^{-2}+H_{j-1}/(2j)$，其中 $H_k=\sum_{i=1}^k 1/i$。由 $H_{j-1}\le j/2$（$j\ge2$）即得表中 $13/36$。当 $j=2$ 时，$\Lambda(x)^2$ 与唯一的平方因子项各贡献 $\ell^2/4$。这里使用 $\Lambda(x)^2$ 使平方邻项也落在 $1/2$。

**定理 1.** 对任意固定 $13/36<\theta<1/2$，记 $\mathcal Q(X)=\{X<q\le2X:q\text{ 为素数},\ q\equiv5\pmod6\}$。则

$$
\#\{q\in\mathcal Q(X): V_q>\theta\}\asymp_\theta X/(\log X)^2,\qquad
\sum_{q\in\mathcal Q(X)}K_q\asymp X .
\tag{15}
$$

**下界。** Green–Tao 的定理 6.1 给出每个充分大的 $(X,2X]$ 中 $\gg X/\log^2 X$ 个素数 $q$，使 $q+2$ 为素数，或为两个均大于 $q^{3/11}$ 的素数之积。此处用的是该文明确陈述的二倍区间版本，作者在脚注 3 说明其与 Iwaniec 原叙述的区间差异。由于 $3/11>1/4$，这些素因子均超过 $x^{1/4}$。且充分大时这些 $q$ 必有 $q\equiv5\pmod6$，否则 $3\mid q+2$。因此每个这样的 $q$ 均有 $V_q\ge1/2-O(1/(L\sqrt{\log L}))$。又 $G\ge0$，据 (14) 及 $\#\mathcal Q(X)\ll X/\log X$，所有导子的删除误差总计 $O(X/\sqrt{\log\log X})$。这就给出 (15) 中有符号和的正下界。[Green–Tao，定理 6.1 与脚注 3](https://www.numdam.org/article/JTNB_2006__18_1_147_0.pdf#page=20)

**上界。** 表中前两类都要求 $q(q+2)$ 没有小于 $X^{1/8}$ 的素因子。对两个剩余类用维数二上界筛，这类 $q$ 至多 $O(X/\log^2 X)$ 个。具体地，可用水平 $X^{1/3}$ 的算术大筛，仅筛 $3\le p\le X^{1/24}$ 的两个剩余类。分母中不截断的乘积 $\prod_p(1-2/p)^{-1}\gg\log^2 X$。将乘积展开归一化为概率，所选平方自由数 $d$ 的 $\log d$ 均值为 $2\sum_{3\le p\le X^{1/24}}\log p/p=(1/12+o(1))\log X$，所以 $d\le X^{1/3}$ 保留至少半数质量。大筛分子为 $O(X)$。由表及 (14)，对固定的上述 $\theta$，充分大时 $V_q>\theta$ 恰好选出前两类。[算术大筛的剩余类形式](https://kskedlaya.org/ant/chap-largesieve2.html)

其余有非零算术权的 $x$ 或为高次幂，或被一个至少为 $X^{1/4}$ 的高次素数幂整除。后者至多有 $O(X^{7/8}\log X)$ 个，因为 $u$ 以下的高次幂数为 $O(\sqrt u\log(2u))$，由部分求和得 $\sum_{h\ge X^{1/4},\,h\ {\rm proper\ prime\ power}}h^{-1}\ll X^{-1/8}\log X$。每个 $x$ 的 $G(x)\le\ell^2$。这给出 (15) 的另一个上界。陈氏定理的下界结合表中间隔也证明了 (1)。在每个充分大的导子上，$V_q>3/4$ 当且仅当 $q+2$ 为素数。

## 5. 消去半素数后的矩条件

设 $P_2(X)$ 计数 $X<q\le2X$ 中的真孪生素数对 $(q,q+2)$，并令 $L_X=\log X$。多项式 $2u^2-u$ 在 $u=0,1/2$ 为零，在 $u=1$ 为一。由 (14)、取值表及第 4 节的稀有高次幂估计，

$$
\frac{L_X^2}{X}\sum_{q\in\mathcal Q(X)}(2V_q^2-V_q)
=\frac{L_X^2}{X}P_2(X)
 + O\!\left(\frac{1}{\sqrt{\log\log X}}+
              X^{-1/8}L_X^3\right).
\tag{16}
$$

这里求和所有逐点归一化误差后为 $O(X/(L_X^2\sqrt{\log\log X}))$。多项式在取值范围的一个固定邻域内导数有界，故相同估计适用于替换其自变量。另由 (15) 及 (14)，$\sum V_q\asymp X/L_X^2$，可直接由 $G(x)/\ell^2$ 的非负性和同样的总误差得到。因此

$$
\limsup_{X\to\infty}\frac{\sum_{q\in\mathcal Q(X)}V_q^2}
                              {\sum_{q\in\mathcal Q(X)}V_q}>\frac12
\quad\Longleftrightarrow\quad
\limsup_{X\to\infty}\frac{L_X^2P_2(X)}{X}>0 .
\tag{17}
$$

陈氏定理保证这里的分母有正确数量级。半素数给出的比值恰在 $1/2$ 的极限边界上。要取得真孪生素数的正比例下界，还需要零点的四次平均在 (17) 中有固定正超额。只要求无穷性时，相应的逐导子条件是 (1)。这个超额没有从现有的陈氏下界或相邻零点间距统计中得到。

## 6. 复算

程序对每个导子、每个尺度及每个特征分别计算 $S_\chi(y)$。第一条路径从 Hurwitz $\zeta$ 及其导数做有限矩形围道积分，扣除 $s=0$ 的平凡零点。第二条路径用单特征的两个素数幂和，在 $\Re s=1/2$ 上独立计算伽马积分。把 (13) 从 $-1/2$ 移到 $1/2$ 时跨过的留数为 $-H(0)$，所以第二条路径已包含该修正。最后形成完整的特征对矩阵，分别求出被删除的部分与保留的部分。表中的三种和均除以 $\log^2(q+2)$。

| $q$ | $q+2$ | 完整和 | 同特征或共轭部分 | 保留的 $V_q$ |
|---:|---|---:|---:|---:|
| 29 | $31$ 素数 | 1.00000000 | 0.14418515 | 0.85581485 |
| 47 | $49=7^2$ | 0.50000000 | 0.04264171 | 0.45735829 |
| 53 | $55=5\cdot11$ | 0.50000000 | 0.03888884 | 0.46111116 |
| 173 | $175=5^2\cdot7$ | 0.25000000 | 0.01356828 | 0.23643172 |
| 383 | $385=5\cdot7\cdot11$ | 0.00000000 | 0.01583236 | $-0.01583236$ |

这里的小导子也展示了删除项残留的有限尺度偏差。程序另对 $17\le q\le10^5$ 的全部 $4\,804$ 个候选素数核对算术分类，其中 $1\,221$ 个对应真孪生素数，$1\,247$ 个对应表中的大因子半素数。数值表使用浮点围道和独立算术复算，渐近定理的误差来自第 2、3 节的解析估计。[新复算程序](https://github.com/lixiang90/vibemath-NumberTheory/blob/main/zeros-exact-distribution/gap-two-chen-zero-pairs/verify.py) · [共用的 Mellin 权与 Dirichlet 函数求值程序](https://github.com/lixiang90/vibemath-NumberTheory/blob/main/zeros-exact-distribution/gap-two-single-prime-conductor/verify.py)
