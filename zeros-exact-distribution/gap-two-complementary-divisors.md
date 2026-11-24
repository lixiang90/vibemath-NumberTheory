---
title: "固定差 2 的互补除数与偶数因子余项"
date: "2026-11-24"
abstract: "在 p+2 的平方根处配对除数，奇数个素因子的无平方合数逐点抵消，四素因子项具有有利的符号。给出小因子半素数的单边界，以及移动截点对应的 Dirichlet 零点族平均估计。Dirichlet GRH 将后一个估计所需的零点截高降低一个平方根因子。"
---

# 固定差 2 的互补除数与偶数因子余项

怎样从恢复孪生素数的位置走向估计它们的数量？在最近三篇笔记中，第 21 篇用 $\zeta$ 零点差的两个傅里叶峰检测一对素数，能识别极其稀疏的孪生素数，但须证明两个指定频率同时出现峰值。第 22 篇把 $p+2\equiv0\bmod d$ 写进 Dirichlet 特征相位 $\overline{\chi(-2)}$，让零点线性出现。第 23 篇用素数在等差数列中的分布算出了小模数主项，却留下了两组绝对质量都比 $X$ 大的有符号除数和。本文沿后两篇的路线，先在**每个具体的 $p+2$ 内**做抵消，再对 $p$ 求和。

令 $m=p+2$。当 $m$ 是无平方因子的合数时，这次配对消去了所有奇数个素因子的邻项。恰有四个素因子的邻项只会改善孪生素数下界。恰有两个素因子的邻项仍是障碍，但较小因子 $\le X^\eta$ 的部分只占 $O_W(\eta X)$。这些结论及下文的误差界均为无条件推导。截断的 von Mangoldt 除数权已有 Goldston–Yıldırım 的系统研究。这里取随邻项移动的截点 $\sqrt m$，考察该权的因子构型和固定差 $2$ 的单边误差。[Goldston–Yıldırım，截断除数和及其与素数的混合相关](https://arxiv.org/abs/math/0111212) · [Murty–Vatwani，固定差素数对的奇偶问题](https://mrammurty-ai.github.io/TwinPrimes-Parity.pdf)

## 1. 一个可以无条件估计的总量

固定 $0\le W\in C_c^\infty((0,\infty))$，支撑于 $[1,2]$，$\int_1^2W(u)\,du=1$。下文 $p,q_i,\ell$ 在素数和或素数乘积中均为素数。写 $w_p=(\log p)W(p/X)$，$L=\log X$，并设

$$
\mathfrak S(2)=2\prod_{\ell>2}\left(1-\frac1{(\ell-1)^2}\right),\qquad
M_W(X)=\sum_{p,\,p+2\ {\rm prime}}w_p\log(p+2).
$$

对 $m>1$ 定义

$$
T(m)=\sum_{d\mid m}\mu(d)\left(\log\frac{\sqrt m}{d}\right)_+,
\qquad y_+=\max(y,0).
\tag{1}
$$

记 $\mathcal A_z$ 为满足 $m=p+2$ 无平方因子且最小素因子 $P^-(m)>z$ 的素数 $p$。令 $\omega(m)$ 为不同素因子数，定义

$$
A_{\rm sf}=\sum_{p\in\mathcal A_z} w_p T(p+2),\qquad
E_j=\sum_{\substack{p\in\mathcal A_z\\ \omega(p+2)=j}} w_p T(p+2),\qquad
E_{\ge6}^{\rm even}=\sum_{j\ge6,\ 2\mid j}E_j .
$$

**命题 1.** 对以下任一种筛限，均有 $A_{\rm sf}=\mathfrak S(2)X+O_W(\mathcal E)$。

| 筛限 $z$ | $\mathcal E$ |
|---|---|
| $\exp(L^{1/8})$，与上一轮相同 | $X(\log L)^4/L^{3/2}$ |
| $L^{16}$ | $X(\log L)^8/L^2$ |

对充分大的 $X$，精确恒等式及其渐近形式为

$$
\begin{aligned}
M_W(X)&=2A_{\rm sf}-2(E_2+E_4+E_{\ge6}^{\rm even}),\\
&=2\mathfrak S(2)X-2(E_2+E_4+E_{\ge6}^{\rm even})+O_W(\mathcal E).
\end{aligned}\tag{2}
$$

式中 $2\mathfrak S(2)X$ 是整个除数权产生的项。合数项在同一个 $X$ 尺度上仍须保留。精确式来自第 2 节，总量估计见第 3 节。

## 2. 每种因子构型的符号

设 $m=q_1\cdots q_j$ 无平方因子，$q_1<\cdots<q_j$。写 $s=\frac12\log m$。用 $d\leftrightarrow m/d$ 和 $u_+-(-u)_+=u$，对每个这样的 $m$ 有

$$
(1-\mu(m))T(m)=\sum_{d\mid m}\mu(d)(s-\log d)=\Lambda(m).
\tag{3}
$$

末式用了 $\sum_{d\mid m}\mu(d)=0$ 和 $-\sum_{d\mid m}\mu(d)\log d=\Lambda(m)$。所以 $m$ 为素数时 $T(m)=\frac12\log m$，$j\ge3$ 为奇数时 $T(m)=0$。对充分大的 $X$，所有孪生素数的第二端点都超过 $z$。这给出 (2) 的第一行。

* 若 $j=2$，写 $m=ab$，$a<b$，则 $T(m)=\log a$。
* 若 $j\ge3$ 且 $q_j^2>m$，则 $T(m)=0$。截点以下的除数恰好是 $m/q_j$ 的全部除数，而这个数至少有两个素因子。对它求 $\sum\mu(d)$ 及 $\sum\mu(d)\log d$ 都得零。
* 若 $j=4$，写 $a,b,c,d$ 为四个素因子的**对数**，依次递增。此时

$$
T(m)=-\min\{2a,(a+b+c-d)_+\}\le0.
\tag{4}
$$

确实，$d\ge a+b+c$ 时用上一条。否则四个单因子都在截点以下，三个入选的双因子之对数为 $a+b,\ a+c,\ \min(a+d,b+c)$。代入 (1) 得 $d-a-\min(a+d,b+c)$，即 (4)。这说明对下界直接舍去 $-2E_4$ 是允许的。

六因子层已有两种符号。下表的 $p$ 和列出的六个因子均为素数，最后一列用整数和有理数精确计算。第 2 个例子也满足最大因子小于 $\sqrt{p+2}$。

| $p$ | $p+2$ 的素因子 | $\exp(2T(p+2))$ |
|---:|---|---:|
| $1\,956\,953$ | $5,7,11,13,17,23$ | $13731656640625/279841>1$ |
| $210\,345\,133$ | $3,5,7,11,13,14009$ | $196252081/225450225<1$ |

这些例子检验的是逐点符号。它们不估计任一增长筛限下该层的总质量。

## 3. 移动截点的平均值及误差

先容许 $\mathcal A_z$ 中的 $m$ 含平方因子，仍要求 $P^-(m)>z$。以 $P=\prod_{2<\ell\le z}\ell$ 和 $V(z)=\prod_{2<\ell\le z}(1-1/(\ell-1))$ 记预筛乘积。展开 (1) 时先限制 $(d,2P)=1$，再对 $(m,P)=1$ 容斥。在准确筛选的项上前一个限制自动成立。取 $k=\lceil20\log L\rceil$，容斥保留 $\omega(r)<k$，逐点误差绝对值至多 $\sum_{r\mid(m,P),\,\omega(r)=k}1$。模数为 $q=dr$，有 $p\equiv-2\bmod q$，且可以从 $q$ 唯一恢复 $r=(q,P)$。这个唯一性允许在分布误差中对每个模数只计一次，截断误差的第 $k$ 层另作同样估计。

以下给出移动截点所需的一致性。对 $u\in[1,2]$ 令

$$
F_d(u)=W(u)\left(\log\frac{\sqrt{Xu+2}}d\right)_+,\qquad
Q_0=\frac{\sqrt X}{L^B},\qquad H=1+B\log L+k\log z .
$$

$F_d$ 绝对连续，因子在进入截点时取值零。固定足够大的 $B$。对 $q\le Q_0$，$\log(Xu)F_d(u)$ 的导数绝对值几乎处处 $\ll_W L^2$。Bombieri–Vinogradov 和分部求和给出 $O_W(X/L^5)$ 的总误差。对 $q>Q_0$，有 $d=q/r>Q_0/z^k$，故导数绝对值的共同上界降至 $O_W(LH)$。所有这些模数都 $\le\sqrt{2X+2}\,z^k$，所跨二进段数为 $O(H)$。

写 $t=Xu$。对上述每个二进段 $(V,2V]$，Bombieri–Friedlander–Iwaniec 的固定剩余类定理给

$$
\sum_{\substack{V<q\le2V\\2\nmid q}}
\left|\pi(t;q,-2)-\frac{\pi(t)}{\varphi(q)}\right|
\ll \frac{X\{H^2+(\log L)^2\}}{L^3}+O(X/L^{10}).
\tag{5}
$$

这里 $|\log(2V)/\log t-1/2|\ll H/L$，$\sum_{V<q\le2V}1/\varphi(q)\ll1$。普通素数定理允许把 $\pi(t)$ 换成 $\operatorname{li}(t)$，对所有模数总计仍可节省任意固定的对数次幂。对每个 $u$ 先用 (5) 求和，再对 $u$ 积分，正好适用于依赖 $d$ 的测试函数。这样得到大模数部分的总误差 $O_W(XH^4/L^2)$。此处取绝对值的次序不要求把定理加强成对 $t$ 取上确界的版本。[Bombieri–Friedlander–Iwaniec，第三篇](https://doi.org/10.1090/S0894-0347-1989-0976723-6) · [Fiorilli 所列的原定理，定理 1.2](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/6F019883E34D9946AFD27A30A4061440/S0008414X00009342a.pdf/on_a_theorem_of_bombieri_friedlander_and_iwaniec.pdf#page=2)

容斥第 $k$ 层的主项及有限乘积展开的尾部绝对值合计

$$
\ll_W XL^2\,\frac{1}{k!}
       \left(\sum_{2<\ell\le z}\frac1{\ell-1}\right)^k
\ll_W XL^{-18}.
$$

这里用了 $\sum_{d\le\sqrt{2X+2}}1/\varphi(d)\ll L$，以及对两种 $z$ 均成立的 $\sum_{\ell}1/(\ell-1)\le\log L$（和限同上，$X$ 充分大）。绝对尾部也可由指数级数估计，因为第 $k$ 项以后的相邻项之比小于 $1/2$。故含平方因子在内的主项为

$$
X\int_1^2 W(u)V(z)
\sum_{\substack{d\le y\\(d,2P)=1}}
\frac{\mu(d)}{\varphi(d)}\log\frac yd\,du,\qquad
y=\sqrt{Xu+2}.
\tag{6}
$$

**(6) 的系数。** 先只剔除素数 $2$。奇数上的 Dirichlet 级数有分解

$$
\sum_{2\nmid d}\frac{\mu(d)}{\varphi(d)d^v}=\frac{J(v)}{\zeta(1+v)},\qquad
J(0)=\mathfrak S(2).
$$

$J$ 的奇素数因子为 $1-\ell^{-v}/[\ell(\ell-1)(1-\ell^{-1-v})]$，在 $2$ 处为 $(1-2^{-1-v})^{-1}$，因而其系数在 $\Re v=-1/2$ 绝对可和。由经典 Möbius 素数定理、分部求和和卷积，$\sum_{d\le y,\ 2\nmid d}\mu(d)/\varphi(d)=O_A((\log y)^{-A})$，相应的 $\sum -\mu(d)\log d/\varphi(d)=\mathfrak S(2)+O_A((\log y)^{-A})$，对每个固定 $A$ 成立。具体可在 $\sqrt y$ 分开卷积，再用 $\sum_{n\le u}\mu(n)\ll u\exp(-c\sqrt{\log u})$。

剔除 $P$ 使系数再与 $\prod_{\ell\mid P}(1-\ell^{-v}/(\ell-1))^{-1}$ 的非负系数 $g_z(b)$ 卷积。$V(z)\sum_b g_z(b)=1$。该乘积在 $v=-1/(10\log z)$ 为 $(\log z)^{O(1)}$，故 $b>\sqrt y$ 的质量至多 $(\log z)^{O(1)}\exp[-\log y/(20\log z)]$。对 $b\le\sqrt y$ 将上述对数加权估计用于 $y/b$，即得 (6) 中 $V(z)$ 乘内和为 $\mathfrak S(2)+O_A(L^{-A})$，对 $u$ 一致。这里两种 $z$ 均满足 $\log y/\log z\gg L^{7/8}$，在充分大的 $X$ 上成立。

最后，仍有 $\ell^2\mid m$、$\ell>z$ 的邻项至多 $O(X/z+\sqrt X)$ 个。逐项有 $|w_p T(m)|\ll_W L^2\tau(m)$，而 $\sum_{m\le2X+2}\tau(m)^2\ll XL^3$。Cauchy–Schwarz 给出的删除误差为 $O_W((X/\sqrt z+X^{3/4})L^{7/2})$。将两种 $z$ 依次代入此界和 $XH^4/L^2$，即得命题 1。

### 移动截点对应的零点族估计

这里可以直接看出中心线信息改善了哪一项。设 $\mathcal B$ 为上面使用的 $(d,r)$，即无平方的 $d\le\sqrt{2X+2}$、$(d,2P)=1$、$r\mid P$、$\omega(r)\le k$，写 $q=dr$、$Q=\sqrt{2X+2}\,z^k$。对 $G=F_d$ 写 $\widehat G(s)=\int_0^\infty G(u)u^{s-1}\,du$。记

$$
B_{d,r}(T)=\frac1{\varphi(q)}\sum_{\chi\bmod q}\overline{\chi(-2)}
\sum_{\substack{\rho\ {\rm nontrivial\ of}\ L(s,\chi^*)\\|\Im\rho|\le T}}
X^\rho\widehat F_d(\rho).
$$

$\chi^*$ 表示诱导本原特征，主特征使用 $\zeta$。取齐正负高度并计重数，则

$$
\sum_{(d,r)\in\mathcal B}|B_{d,r}(T)|\ll_W \frac{XH^4}{L^2},
\quad
\begin{cases}
T=QL^6 & \text{无条件},\\
T=(Q/\sqrt X)L^6 & \text{假设 Dirichlet GRH}.
\end{cases}\tag{6a}
$$

先取全部高度。由显式公式，每项与对应的素数剩余类偏差相差素数幂、缺失欧拉因子和平凡零点的修正。对不同 $(d,r)$，模数 $q$ 各不相同。高次素数幂 $\ell^j\asymp X$ 的合计绝对误差 $\ll_W L\sum_{j\ge2,\,\ell^j\asymp X}\log\ell\,\tau(\ell^j+2)\ll_{W,\nu} X^{1/2+\nu}$，其中 $\nu>0$ 任意固定。缺失欧拉因子的修正每模数 $\ll_W L\log(2q)$，平凡零点每模数 $\ll_W L$。由 $Q=X^{1/2+o(1)}$ 和证明 (5) 后对移动权所作的分部求和，得到 (6a) 的全高度版本。分部求和本身就在每个模数上取了绝对值。

在对数坐标中，$F_d$ 连续且一阶导数的全变差为 $O_W(L)$，包括 $Xu+2=d^2$ 处的一次跳变。分部积分两次给 $|\widehat F_d(\beta+i\gamma)|\ll_W L/(1+|\gamma|)^2$，$0\le\beta\le1$，足以用光滑逼近证明该权的显式公式并保证零点和绝对收敛。本原零点计数给出截高误差的总上界 $O_W(QX^b L\log(2QT)/T)$，无条件取 $b=1$，在 GRH 下取 $b=1/2$。式 (6a) 的两个截高都把此误差压至 $O_W(X/L^4)$。例如 $z=L^{16}$ 时，GRH 下 $T=\exp(O((\log\log X)^2))$，且 (6a) 的右端为 $O_W(X(\log\log X)^8/(\log X)^2)$。这里控制的是各模数内保留特征相位后的绝对值之和。[本原 Dirichlet 零点计数](https://arxiv.org/abs/2005.02989)

## 4. 含小素因子的半素数可以付出多大代价

**命题 2.** 对 $0<\eta\le1/4$，将 $E_2$ 限于较小素因子 $a\le X^\eta$，所得非负量满足

$$
E_{2,\,a\le X^\eta}\ll_W X(\eta+L^{-1}).
\tag{7}
$$

常数与 $\eta$ 无关。固定奇素数 $a$，写 $p+2=ab$。在长 $\asymp X/a$ 的 $b$ 区间上，两线性式 $b,ab-2$ 同取素数的个数 $\ll (X/a)/\log^2(X/a)$，对 $a\le X^{1/4}$ 一致。事实上，除 $2,a$ 外每个素数筛掉两个剩余类，在 $2,a$ 处各筛掉一个，局部修正一致有界。用二维 Selberg 上界筛即可，例如取筛限为 $(X/a)^{1/20}$，区间取整的误差由每个无平方模数 $d$ 至多 $2^{\omega(d)}$ 个根控制。乘 $O_W(L\log a)$ 后对 $a$ 求和，用 $\sum_{a\le v,\ a\ {\rm prime}}\log a/a=\log v+O(1)$ 得 (7)。[Selberg 上界筛及素数对应用](https://terrytao.wordpress.com/2015/01/21/254a-notes-4-some-sieve-theory/)

同一论证还能控制**每个固定的因子层**。设 $j\ge2$ 固定，$0<\eta\le1/4$。有

$$
\limsup_{X\to\infty}\frac1X
\sum_{\substack{p\in\mathcal A_z,\ \omega(p+2)=j\\ P^-(p+2)\le X^\eta}}
 w_p\,|T(p+2)|
\ll_{j,W}\eta\big(1+\log(1/\eta)\big)^{j-2}.
\tag{8}
$$

为证 (8)，将 (1) 中含有最小素因子 $a$ 的除数与不含它的除数配对。函数 $u_+$ 为 1-Lipschitz，故 $|T|\le 2^{j-1}\log a$。固定前 $j-1$ 个素因子，乘积记为 $b$，令最大素因子为 $v$。由于 $v\ge(p+2)^{1/j}$，有 $X/b\gg_j X^{1/j}$。对 $v,bv-2$ 用同一个上界筛，$b$ 的至多 $j-1$ 个局部修正的乘积 $\ll_j1$。放宽中间 $j-2$ 个素因子的排列限制，所求量除以 $X$ 至多为常数乘

$$
\frac1L\sum_{a\le X^\eta,\ a\ {\rm prime}}
\frac{\log a}{a}
\left(1+\log\frac L{\log a}\right)^{j-2}.
$$

这里用了 $\sum_{a\le q\le 2X+2,\ q\ {\rm prime}}1/q\le\log(L/\log a)+O(1)$。分部求和，令 $v=\log a/L$，上极限由 $\int_0^\eta(1+\log(1/v))^{j-2}\,dv$ 控制，得到 (8)。其常数依赖 $j$，故本文在无界多个因子层的总和中保留高层余项。

## 5. 对零点和及孪生素数下界的含义

给定固定 $\varepsilon>0$，定义 $D=2X+2$，$\widehat W(s)=\int W(u)u^{s-1}\,du$。对模 $d$ 特征 $\chi$ 记其诱导本原特征为 $\chi^*$，主特征取 $\zeta$。固定差的线性谱和为

$$
Z_\varepsilon(X)=
\sum_{\substack{d\le D\\2\nmid d}}\frac{-\mu(d)\log d}{\varphi(d)}
\sum_{\chi\bmod d}\overline{\chi(-2)}
\sum_{\substack{L(\rho,\chi^*)=0,\ 0<\Re\rho<1\\|\Im\rho|\le X^\varepsilon}}
X^\rho\widehat W(\rho).
$$

内层计正负高度和重数。光滑显式公式给 $M_W=C(D)X-Z_\varepsilon+O_{W,\varepsilon}(X^{2/3})$，其中 $C(D)=\sum_{d\le D,\,2\nmid d}-\mu(d)\log d/\varphi(d)=\mathfrak S(2)+O_A(L^{-A})$。简述其误差来源如下。在 $\Lambda(n+2)=-\sum_{d\mid n+2}\mu(d)\log d$ 中先取奇 $n$，特征正交即给相位 $\overline{\chi(-2)}$。改用全部 $\Lambda(n)$ 后，高次幂端点有 $O(\sqrt X)$ 个，至多带 $O(L^2)$ 的算术权，偶数第一端点只可能是 $2$ 的幂。诱导特征缺失的因子可写成 $d=\ell e$，$n=\ell^k$。由于 $d$ 无平方，模 $e$ 的正交性强制 $e\mid\ell^k+2$，绝对修正每个 $\ell^k$ 至多 $O(L\log\ell\,\tau(\ell^k+2)/(\ell-1))$。$k=1$ 时只可能 $d=\ell$，$k\ge2$ 时对高次幂用除数界。平凡零点按奇偶特征投影，仅 $d\mid3$ 或 $d\mid1$ 留下，另扣去 $\zeta(0)$ 非零的修正，合计 $O_W(1)$。$\widehat W$ 的任意幂次衰减和标准零点计数把 $|\Im\rho|>X^\varepsilon$ 的总贡献压至 $O_{W,\varepsilon}(1)$。这些界均在 $O_{W,\varepsilon}(X^{2/3})$ 内。[Dirichlet 特征的显式公式](https://kskedlaya.org/ant/chap-errorbounds2.html) · [Bennett 等，本原 Dirichlet 零点计数](https://arxiv.org/abs/2005.02989)

合并 (2) 可得

$$
\boxed{Z_\varepsilon(X)=2(E_2+E_4+E_{\ge6}^{\rm even})
-\mathfrak S(2)X+O_{W,\varepsilon}(\mathcal E).}
\tag{9}
$$

于是 $\limsup_{X\to\infty}M_W(X)/X>0$ 当且仅当
$\liminf_{X\to\infty}(E_2+E_4+E_{\ge6}^{\rm even})/X<\mathfrak S(2)$。要得到一个只涉及较大因子的充分条件，可取固定 $\delta>0$ 和充分小的固定 $\eta>0$，使 (7) 中 $C_W\eta<\delta/2$。若在无穷多个趋于无穷的 $X$ 上有

$$
E_{2,\,a>X^\eta}+E_{\ge6}^{\rm even}\le(\mathfrak S(2)-\delta)X,
\tag{10}
$$

则由 (2)、(4)、(7) 得这些窗口中有 $\gg_{W,\delta}X/\log^2X$ 对孪生素数。(10) 中高层只需计最大素因子小于 $\sqrt{p+2}$ 的邻项。定量计入 $E_4$ 的负贡献可进一步放宽充分条件。六因子实例说明，高层内部仍有可利用的正负抵消。

还需估计两因子都大的半素数，以及最大因子受限的偶数高层的**带符号总和**。Lichtman–Teräväinen 以及 Lichtman 对 $\sum_{p\le x}\mu(p+2)$ 的结果估计全体素数上的奇偶相关。这里的 $T(p+2)$、因子大小条件和单边阈值都须保留。Alpöge–Claude 从幂次窗口排除偏离中心线的极点，由此获得的实部信息可用于 (6a) 的截高。(9) 还含跨模数相位 $\mu(d)\overline{\chi(-2)}$。本文没有得到 (10) 的固定正余量，孪生素数无穷性仍待证明。[Lichtman–Teräväinen，移位素数上的 Möbius 相关](https://arxiv.org/abs/2610.20242) · [Lichtman，定理 1.2](https://arxiv.org/html/2611.28392v1#S1) · [Alpöge–Claude 归档中的极点导数引理](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/DedekindLFunctions/LandauTheorems/LandauPoleAsym.lean)

## 6. 有限核验

取 $z=\exp(L^{1/8})$。在 $1<u<2$ 上取 $W(u)$ 正比于 $\exp(4-1/((u-1)(2-u)))$ 并归一化为积分 $1$。左端直接筛出真孪生素数，右端枚举第二端点的除数。表中质量均已除以 $X$。

| $X$ | 孪生素数对数 | $M_W/X$ | $A_{\rm sf}/X$ | $E_2/X$ | $E_4/X$ | $E_{\ge6}^{\rm even}/X$ |
|---:|---:|---:|---:|---:|---:|---:|
| $1\,000$ | 26 | 1.2678954 | 1.1178438 | 0.4838961 | 0 | 0 |
| $10\,000$ | 137 | 1.1783113 | 1.2668744 | 0.6874090 | −0.0096903 | 0 |
| $100\,000$ | 936 | 1.3148533 | 1.2966215 | 0.6852474 | −0.0460526 | 0 |
| $1\,000\,000$ | 6702 | 1.3476002 | 1.3174400 | 0.7193715 | −0.0757316 | $5.23\times10^{-13}$ |

四个窗口中，(2) 第一行两侧除以 $X$ 后的数值差均小于 $4.5\times10^{-16}$。另对 $3\le m\le100000$ 的全部 40526 个奇无平方整数，用 $\exp(2T(m))=\prod_{d\mid m,\,d^2<m}(m/d^2)^{\mu(d)}$ 做有理数精确核验，检查 (3)、(4)、半素数公式和大素因子消去。附带的 Python 程序只用标准库。有限核验检验恒等式和符号，(5)–(8) 的一致性来自解析证明。[复现程序](gap-two-complementary-divisors/verify.py)
