---
title: "原双线性和的幂指数集中在一个大素数与一个小余因子上"
date: "2026-11-27"
abstract: "保持原式（9）的系数和截断不变。对任意固定的小正数 τ，只保留超过 T 的 5/6−τ 次幂的素数，仍得到原双线性和在每个充分远完整窗口上的精确正部幂指数。证明同时控制所有高次素数幂和与第一因子共享素因子的项，并区分已达到次幂预算的区间和只证明了较小增长指数的区间。"
---

# 原双线性和的幂指数集中在一个大素数与一个小余因子上

在判别式 $229$ 的例子中，令 $F=\mathbb Q(\alpha)$、$\alpha^3-4\alpha+1=0$、$E=F(i)$，取 $M/B=F^{\rm gal}(i)/\mathbb Q(\sqrt{229},i)$ 的三阶特征 $\psi$。写 $L_B(s,\psi)=L_+(s)L_-(s)$，$L_-$ 为 $L_+$ 的 $\chi_{-4}$ 扭转。设 $\beta_*$ 为两因子非平凡零点实部的上确界，$\delta=\tfrac32(\beta_*-1/2)$。下文的等价判据保留相对于 $\mathbb Q$ 的这两个次数二因子。[两个 Hecke 因子](cm229-heights-and-zeros.md)

仍令 $X=Y^{3/2}$，$W=\lceil Y^3\rceil$，$k_n=\lfloor W/n^{3/2}\rfloor$，$r=2^{1/\lceil\log^2(2Y)\rceil}$。按 $T=2^jY$、$T\le n<2T$ 分块，限于 $Y\le n<X$。使用原式的 $z=\max(1,\lfloor T^{1/6}/\log^3(2T)\rfloor)$ 和 $c_z(h)=\mathbf1_{h=1}-(u_{\le z}*a)(h)$，其中 $a,u$ 是 $L_+,1/L_+$ 的系数。对任意实值序列 $G$ 写

$$
 \|G\|_Y=W^{-1}\sum_n k_n|G(n)|,\qquad P(G;Y)=W^{-1}\sum_n k_n[G(n)]_+ .
$$

**定理。** 固定 $0<\tau<1/3$。从原双线性和式（9）中逐项选出

$$
 Q_\tau(n)=\frac12
 \sum_{\substack{p>T^{5/6-\tau}\ {\rm prime},\ h>z\\2n/r<hp\le2n}}
 (1+\chi_{-4}(hp))c_z(h)a(p)\log p . \tag{1}
$$

则对每个 $\varepsilon>0$，在所有充分远的完整窗口上

$$
 \|R_Y-Q_\tau\|_Y\ll_{\tau,\varepsilon}Y^{(1-\tau)\delta+\varepsilon},
 \qquad
 \lim_{Y\to\infty}\frac{\log(1+P(Q_\tau;Y))}{\log Y}=\delta. \tag{2}
$$

这里 $z<h<4T^{1/6+\tau}$，并且 $p>2\sqrt T$ 对充分大的 $T$ 成立。因此 $p\nmid h$，$p$ 是乘积的唯一最大素因子。比如取 $\tau=1/6$，只须在原系数上处理 $p>T^{2/3}$、$z<h<4T^{1/3}$。可将 $\tau$ 取为任意小的固定正数。若证明 (1) 的累计正部满足 $WY^{o(1)}$（未归一化），就得到本域的 ERH。[原式的系数与完整窗口判据](cm229-cubic-hecke-carries.md)

这里有两个不同强度的估计。伴篇对 $z<\ell\le\sqrt T/z$ 的贡献已无条件证明 $\log^C(2Y)$ 的归一化绝对值上界。继续删到 $T^{5/6-\tau}$ 时，(2) 给出的是与未知 $\delta$ 有关的指数。如果 $\delta>0$，删去部分的指数严格较小，因此原和的幂增长只能由 (1) 保留下来。[已达多对数预算的素因子区间](cm229-type-ii-ranges.md) · [全双线性和的有限高度公式](cm229-type-ii-zeros.md)

## 两个截断级数的估计

对 $\nu=+,-$，令 $u_\nu,b_\nu$ 为 $1/L_\nu,-L_\nu'/L_\nu$ 的系数。单独处理一支时写 $U_z(s)=\sum_{d\le z}u_\nu(d)d^{-s}$、$C_z(s)=1-U_z(s)L_\nu(s)$ 及 $B_{(z,V]}(s)=\sum_{z<\ell\le V}b_\nu(\ell)\ell^{-s}$。两支在末尾取半和。对固定 $1/2<\sigma<1$、任意 $\eta>0$ 和序列 $g=u_\nu$ 或 $b_\nu$，有

$$
 \left|\sum_{d\le D}g(d)d^{-\sigma-it}\right|
 \ll_{\sigma,\eta}(2D)^{(\beta_*-\sigma)_++\eta}(2+|t|)^\eta .
 \tag{3}
$$

先说明使用的竖线界。固定 $\beta_*<b_1<b_2<1$，由 Hadamard 对数导数与单位高度的零点数，在 $b_1\le\Re s\le2$ 有 $L_\nu'/L_\nu\ll_{b_1}\log(2+|\Im s|)$。从 Euler 乘积一侧沿水平段积分，得到同一条带上解析分支 $\log L_\nu(s)\ll_{b_1}\log(2+|\Im s|)$。在 $\Re s=2$，$\log L_\nu$ 有界。对 $\log L_\nu(s)e^{(s-it_0)^2}$ 作三线插值，在 $b_2+it_0$ 得到 $O((\log(2+|t_0|))^{(2-b_2)/(2-b_1)})$。幂严格小于 1，故对任意 $\kappa>0$，$1/L_\nu(s)\ll_{b_2,\kappa}(2+|\Im s|)^\kappa$ 在 $\Re s\ge b_2$ 成立。$\beta_*<1$ 可用已公布的 $7/8$ 无零半平面。[局部零点计数与对数导数](cm229-critical-logarithm.md) · [Hu 的 Hecke 无零半平面](https://hexagonmath.org/2610.00216v2)

对 $1/L_\nu(s+w)$ 或 $-L_\nu'(s+w)/L_\nu(s+w)$ 作截断 Perron 积分，$s=\sigma+it$。端点取 $\lfloor D\rfloor+1/2$，起始线 $\Re w=2$，高度 $H=(2D+|t|+2)^5$。半整数距最近整数为 $1/2$，故 $|u_\nu(d)|\le d_2(d)$、$|b_\nu(d)|\le2\Lambda(d)$ 给出截断误差 $O((2D)^3/H)$。移线至 $\Re w=(\beta_*-\sigma)_++\eta'$，其中 $\eta'>0$ 任意小。此处 $w=0$ 与所有极点均在左侧。新线的费用至多为 $C(2D)^{(\beta_*-\sigma)_++\eta'}(2H+|t|)^\kappa\log(2H)$，水平边的费用为 $O((2D)^2(2H+|t|)^\kappa/H)$。先取足够小的 $\eta',\kappa$ 再吸收对数即得 (3)。

## 原和中两个较小的贡献

以下取 $z\le V<T$，各估计的隐常数对 $V$ 一致。

**素数幂与公因子。** 在 $\ell>V$ 的部分，记 $G_{\rm pow}$ 为所有 $\ell=p^j$、$j\ge2$ 的项，$G_{\rm div}$ 为 $\ell=p$、$p\mid h$ 的项，均保留原权。则

$$
 \|G_{\rm pow}\|_Y+\|G_{\rm div}\|_Y
        \ll_\varepsilon Y^{\delta/6+\varepsilon}. \tag{4}
$$

第一部分单支的级数是 $C_z(s)B_{\rm pow}(s)$，其中 $B_{\rm pow}(s)=\sum_{p^j>V,j\ge2}b_\nu(p^j)p^{-js}$ 在 $\Re s>1/2$ 绝对且一致有界于每条固定竖线。对第二部分，令局部逆多项式为 $1-a_\nu(p)Z+d_{\nu,p}Z^2$，包含分歧处的降次情形。因 $p>V\ge z$，$p$ 不整除 $U_z$ 中的任何下标。将 $a_\nu$ 的 Euler 乘积限于 $p\mid h$ 后，得到第二部分的单支级数

$$
 -U_z(s)L_\nu(s)
 \sum_{p>V}b_\nu(p)\big(a_\nu(p)p^{-2s}-d_{\nu,p}p^{-3s}\big). \tag{5}
$$

确切地说，$\sum_{p\mid j}a_\nu(j)j^{-s}=L_\nu(s)(a_\nu(p)p^{-s}-d_{\nu,p}p^{-2s})$，再乘 $b_\nu(p)p^{-s}$ 即得 (5)。末尾的素数级数也在 $\Re s>1/2$ 绝对收敛，因 $|a_\nu(p)|\le2,\ |d_{\nu,p}|\le1$。

**中间区间。** 对所有 $z<\ell\le V$，单支的级数为 $C_z(s)B_{(z,V]}(s)$，后者的 $B$ 是 $b_\nu$ 的有限 Dirichlet 多项式。$C_z$ 中的常数 1 所乘的多项式支撑于 $q\le V<T$，在原乘积区间贡献为零。其余即 $-U_z L_\nu B_{(z,V]}$。由 (3)，两个多项式的乘积在 $\sigma+it$ 至多为

$$
 C_{\sigma,\eta}(2zV)^{(\beta_*-\sigma)_++\eta}(2+|t|)^{2\eta}. \tag{6}
$$

为把这些竖线界用于每一个二进块，取 $\sigma=1/2+\theta$，$\theta>0$ 任意小。次数二函数的凸性界为 $|L_\nu(\sigma+it)|\ll_{\theta,\kappa}(2+|t|)^{1-\sigma+\kappa}$，由反射线上的函数方程和绝对收敛线插值即得。短区间 Mellin 乘子 $(1-r^{-s})/s$ 在该线上为 $O((1+|t|)^{-1})$。因此即便再乘 $(2+|t|)^{2\eta}$，只要 $\kappa+2\eta<\theta$，所得乘积仍对 $t$ 平方可积。

具体实施时，以伴篇的光滑紧支撑核按相对宽度 $1/T$ 平滑乘积区间。上述各和的绝对系数均由 $C d_4(q)\log^2(2q)$ 控制，可由 $|c_{z,\nu}|\le\delta_1+d_4$ 和 $\sum_{\ell\mid q}\Lambda(\ell)=\log q$ 逐项验证。故除去平滑和处理整数端点的平均费用均为 $O(\log^C(2T))$。将平滑 Mellin 积分移至 $\sigma$ 并作 Plancherel，(4) 中两类项的块平均绝对值各至多为 $C_\varepsilon T^{1/2+\varepsilon} z^{\beta_*-1/2}$。对于 (6) 则至多为 $C_\varepsilon T^{1/2+\varepsilon}(zV)^{\beta_*-1/2}$。这是先取足够小的 $\theta,\kappa,\eta$ 再吸收 $T$ 的任意小正幂所得。有限支撑项可在平滑前直接扣去。[平滑、整数端点及 Plancherel 的细节](cm229-type-ii-ranges.md)

每块乘 $T^{-1/2}$，再用 $z\le T^{1/6}$、$T\le X$，得到 (4)。取 $V=\lfloor T^{5/6-\tau}\rfloor$ 时有 $zV\le T^{1-\tau}$。故原和与 (1) 的归一化绝对误差至多为 $O_\varepsilon(Y^{(1-\tau)\delta+\varepsilon})$。若 $\delta>0$，此界可取严格小于 $\delta$ 的指数，而 $P(R_Y;Y)=Y^{\delta+o(1)}$ 对每个充分远窗口成立。正部的 1-Lipschitz 性遂给出 (2) 的等式。若 $\delta=0$，两者都有任意小的正幂上界，$\log(1+P)\ge0$ 给出对应下界。[原和的正部窗口定理](cm229-cubic-hecke-carries.md)

取 $V=\lfloor\sqrt T/z\rfloor$ 时，还可逐项删去 (4) 的两类项，再结合伴篇的多对数界。留下 $\ell=p>\sqrt T/z$ 且 $p\nmid h$ 的原系数和，其归一化正部同样有精确指数 $\delta$。在此选择下，全部删去项的绝对值误差为 $O_\varepsilon(Y^{\delta/6+\varepsilon})$。

## 多窄的因子范围已经有证明

由 $c_z=u_{>z}*a$，在 (1) 再写 $h=dj$ 时，有 $d>z$ 及 $j<4T^{1/6+\tau}/z\ll T^\tau\log^3(2T)$。于是第一变量的卷积逆系数只须在 $z<d<4T^{1/6+\tau}$ 上使用。若要确定问题的幂指数，固定一个 $\tau$ 即可。

还可作一个关于存在性的推论。选定递减到 0 的一列 $\tau_i\in(0,1/3)$。由每个固定 $\tau_i$ 的极限 (2)，可递归选取 $Y_i\to\infty$，使所有 $Y\ge Y_i$ 的对数比值距 $\delta$ 小于 $1/i$。在 $[Y_i,Y_{i+1})$ 使用 $\tau_i$，便有一组随 $Y$ 缓慢收缩的范围，其中 $h=T^{1/6+o(1)}$、$p=T^{5/6+o(1)}$，且每个充分远窗口仍检出原指数。这里的上界 $p\le4T/z$ 来自 $h>z$。该对角论证没有给出可预先指定的收缩速度。目前仍缺固定 $\tau$ 的 (1) 的独立次幂上界。

## 素数对的完整窗口复算

取 $Y=4096$，遍历 $4096\le n<262144$ 的所有整数。表中 $Q_{\rm cop}$ 取 $p>\sqrt T/z,\ p\nmid h$，$Q_{1/6}$ 为 (1) 的 $p>T^{2/3}$。表值均为严格下端点，上端点为表值加 $10^{-6}$。第二行使用不含对数分母的截断，以便在同一窗口检验多个 $z>1$ 的块。

| 截断 | $P(Q_{\rm cop};Y)$ | $P(Q_{1/6};Y)$ | $\lVert R_{z,z}-Q_{1/6}\rVert_Y$ |
| --- | ---: | ---: | ---: |
| 原式 | 0.606624 | 0.532182 | 1.129625 |
| $\lfloor T^{1/6}\rfloor$ | 0.566701 | 0.390309 | 0.852363 |

程序分别按互素条件和按 (5) 的欧拉多项式计算，再与被删去的素数幂、有限区间共同还原原和。每个乘积前缀都检查了这个恒等式，移动的端点及归一化窗口量用实球算术核验。[素数对及窗口核验程序](verify_cm229_type_ii.py) · [原和的同窗数据](cm229-type-ii-ranges.md)

对这里的 $Q_{1/6}$，后续已刻画原系数可能非零的余因子，并证明可在 $T^{8/15}$ 的半径内平均而将整窗绝对误差控制在 $WY^{o(1)}$。低零点的单个模式在 $z<h\ll z\log^2(2T)$ 的余因子内已有趋于 1 的乘子。[余因子的理想范数条件与低零点](cm229-coprime-cofactors.md) · [原和的高频估计与单窗判据](cm229-coprime-smoothing.md)
