---
title: "原卷积中逐项消失的互素对与三阶 Hecke 零点的余因子范围"
date: "2026-11-27"
abstract: "原余因子系数在某个截断下可能非零，当且仅当它的下标是实二次域的整理想范数。确定仍检出完整窗口正部幂指数的两类余因子，并证明对高度不超过 log T 的每个 Hecke 零点，只用 z 到常数倍 z log² T 的余因子，其光滑 Mellin 和便以任意负幂精度取值 1。"
---

# 原卷积中逐项消失的互素对与三阶 Hecke 零点的余因子范围

取 $\alpha^3-4\alpha+1=0$，$F=\mathbb Q(\alpha)$，$E=F(i)$，$B=\mathbb Q(\sqrt{229},i)$，$M=F^{\rm gal}(i)$。循环三次扩张 $M/B$ 的非平凡 Hecke 特征 $\psi$ 给出 $L_B(s,\psi)=L_+(s)L_-(s)$。这里 $L_+$ 来自 $S_3$ 的二维标准表示，$L_-$ 是其 $\chi_{-4}$ 扭转。记 $a_\nu,u_\nu$ 为 $L_\nu,1/L_\nu$ 的 Dirichlet 系数，$\nu=+,-$。无下标时取正支，写 $\chi=\chi_{-4}$。[域塔、分歧因子与特征分解](cm229-heights-and-zeros.md)

令 $X=Y^{3/2}$，$W=\lceil Y^3\rceil$，$k_n=\lfloor W/n^{3/2}\rfloor$，$r=2^{1/\lceil\log^2(2Y)\rceil}$。整数 $Y\le n<X$ 按 $T=2^jY$、$T\le n<2T$ 分块。采用原来的 $z=\max(1,\lfloor T^{1/6}/\log^3(2T)\rfloor)$ 和 $c_{z,\nu}=u_{\nu,>z}*a_\nu=\mathbf1_{h=1}-u_{\nu,\le z}*a_\nu$。考察原式（9）中的素数对

$$
 Q_T(x)=\frac12\sum_{\substack{p>T^{2/3}\ {\rm prime},\ h>z\\2x/r<hp\le2x}}
 (1+\chi(hp))c_z(h)a(p)\log p . \tag{1}
$$

当 $x=n$，有 $h<4T^{1/3}$，故在充分大的块上 $p>h$。投影因子把 (1) 写成两支同形和的半和，在偶数上也成立。对 $p$ 在 $B$ 中完全分裂的情形，$a(p)$ 在 $M/B$ 相对分裂时取 $2$、相对惰性时取 $-1$。原式还对乘积 $hp$ 作投影，以下一直保留这个权重。[从分裂权重到原双线性和](cm229-cubic-hecke-carries.md)

## 两类余因子保留完整窗口的幂指数

对 $h>1$，写 $h=st$，其中 $s$ 是 $h$ 的全部 $q\le z$ 的素因子连同重数之积，$t=h/s$。两者互素。则

$$
 c_z(h)=a(t)\big(c_z(s)-\mathbf1_{s=1}\big). \tag{2}
$$

因为每个 $d\mid h$、$d\le z$ 必整除 $s$，所以 $a(h/d)=a(t)a(s/d)$，代入定义即得。尤其 $1<s\le z$ 时，所有这样的 $(h,p)$ 均逐项为零。$s=1$ 时 $c_z(h)=-a(h)$。在 $z<h\le z^2$ 内，若一个合数含有 $q>z$，它至多含一个这样大的素因子，且余商在 $(1,z]$。由此，在此范围内只有两类项可能留下来。

- $h=q$ 为素数。此时 $c_z(q)a(p)=-a(q)a(p)$。
- $h$ 的全部素因子均 $\le z$。此时保留原有的 $c_z(h)a(p)$。

这个划分可以直接用于未平滑的整窗量。取 $V_T=\max(T^{2/3},4T/z^2)$，令 $S_T$ 为 (1) 限于 $p>V_T$ 后上面两类项之和。置 $\beta_*=\sup_{\rho}\Re\rho$，遍历 $L_+L_-$ 的非平凡零点，$\delta=\tfrac32(\beta_*-1/2)$。写 $\|G\|_Y=W^{-1}\sum k_n|G_T(n)|$ 和 $P(G;Y)=W^{-1}\sum k_n[G_T(n)]_+$。则对每个 $\varepsilon>0$

$$
 \|Q-S\|_Y\ll_\varepsilon Y^{5\delta/6+\varepsilon},\qquad
 \lim_{Y\to\infty}\frac{\log(1+P(S;Y))}{\log Y}=\delta. \tag{3}
$$

事实上 $p>V_T$ 给出 $h<z^2$。前篇对 $z<\ell\le V$ 全部对数导数系数给出的块绝对均值，按本篇权重归一化后至多为 $T^\varepsilon(zV)^{\beta_*-1/2}$。在 $V=T^{2/3}$ 和 $V=V_T$ 两端相减，再扣掉区间内的高次素数幂。这些幂按前篇式 (4) 的两个尾部作差，其整窗费用为 $O_\varepsilon(Y^{\delta/6+\varepsilon})$。充分大时 $zV_T\ll T^{5/6}\log^3(2T)$，由 $T\le X$ 得 (3) 的第一式。前篇已证明 $Q$ 的正部幂指数为 $\delta$。若 $\delta>0$，差的指数严格较小，正部的 1-Lipschitz 性给出第二式。若 $\delta=0$，任意小的正幂上界已经足够。(3) 对两类项先相加再取正部。[大素数区间定理及中间段估计](cm229-type-ii-large-prime.md)

在 $h=dj$、$d>z$ 的展开中，$p\nmid d$ 还给出 $u(d)a(p)=-u(dp)$。因此 (1) 每个 $c_z(h)a(p)$ 也可精确写为 $-\sum_{dj=h,\ d>z}u(dp)a(j)$。在 $z<h\le2z$ 只有 $d=h$ 一项。在平方自由的 $h$ 上，$u(d)=\mu(d)a(d)$ 对 $d\mid h$ 成立，故 $c_z(h)=a(h)\sum_{d\mid h,\ d>z}\mu(d)$。这些恒等式在两支分别成立。程序从欧拉系数独立构造卷积的两个截断，对 $h\le8192$ 和七个 $z$ 检查了 (2)、首段和混合合数的消失，也对 $dp\le8192$ 的全部互素 $(d,p)$ 核验了逆系数恒等式。[整数与窗口核验程序](verify_cm229_coprime_smoothing.py)

在 $Y=4096$ 的六个块上，原截断均为 $1$，此时 $S_T$ 为空。为检验 (3) 中两类权重如何相消，另取 $z=\lfloor T^{1/6}\rfloor$，保留相同的 $V_T$、窗口和正部权重。素数余因子单独给出 $0.177860$，光滑余因子单独给出 $0.184547$，相加后取正部为 $0.266649$。三者均已除以 $W$，所列下端点的误差小于 $10^{-6}$。逐个整数先加再取正部，实测比两个正部之和小 $(0.095758,0.095759)$。这个有限检验使用另设的截断，(3) 的指数结论使用原截断。[逐整数复算](verify_cm229_coprime_smoothing.py)

## 哪些整数在所有截断下都消失

令 $k=\mathbb Q(\sqrt{229})$，$\chi_k=(229/\cdot)$ 为其本原二次特征。对每个 $h>1$，有精确刻画

$$
 \exists\,1\le z<h:\ c_z(h)\ne0
 \quad\Longleftrightarrow\quad
 v_q(h)\ \text{对每个在 }k\text{ 中惰性的 }q\text{ 均为偶数}
 \quad\Longleftrightarrow\quad h=N_{k/\mathbb Q}\mathfrak a
 \quad\text{对某个整理想 }\mathfrak a. \tag{3a}
$$

在惰性 $q$ 处，标准表示的 Frobenius 是换位，局部逆多项式为 $1-Z^2$。$a(q^e)$ 与 $u(q^e)$ 在奇数 $e$ 上同时为零，因而 $v_q(h)$ 为奇数时每个 $u(d)a(h/d)$ 都是零。这包括 $q=2$。反过来，假设所有惰性指数都是偶数。在三循环素数 $q$ 处，$a(q^e)$ 依 $e\bmod3=0,1,2$ 分别为 $1,-1,0$，且 $u(q)=1$。取 $d$ 为满足 $v_q(h)\equiv2\bmod3$ 的所有三循环素数的乘积。则 $u(d)a(h/d)\ne0$，因为完全分裂处 $a(q^e)=e+1$，分歧处 $a(229^e)=1$。若 $d=1$，得到 $c_1(h)=-a(h)\ne0$。若 $d>1$，则 $d\le\sqrt h<h$，且 $c_d(h)-c_{d-1}(h)=-u(d)a(h/d)\ne0$。最后一个等价由二次域中素理想的范数给出。在扭转支上 $c_{z,-}(h)=\chi(h)c_z(h)$，所以还需 $h$ 为奇数。

式 (2) 的整数除数截断已在 Lean 4.33.1 中核验，使用互素乘法性和 $t$ 的每个素因子均大于 $z$。同一文件对任意有限素因数分解证明了 (3a) 所需的欧拉系数支撑判据，也证明有限前缀在某处非零等价于至少一次增量非零。四种欧拉类型与本域 Frobenius 的对应、理想范数解释及下文的解析估计使用本文的数学证明。[Lean 源码](CM229CofactorSupport.lean)

令 $\mathcal N(H)$ 计 $1\le h\le H$ 中满足上述范数条件的整数。其示性函数的 Dirichlet 级数记作 $D_k(s)$。逐素数平方得到

$$
 D_k(s)^2=\zeta(s)L(s,\chi_k)(1-229^{-s})^{-1}
                 \prod_{\chi_k(q)=-1}(1-q^{-2s})^{-1}.
$$

由此

$$
 \mathcal N(H)=C_k\,\frac{H}{\sqrt{\log H}}\big(1+O(1/\log H)\big),\qquad
 C_k=\left(\frac{229\,L(1,\chi_k)}{228\pi}
                   \prod_{\chi_k(q)=-1}(1-q^{-2})^{-1}\right)^{1/2}. \tag{3b}
$$

事实上，模 $229$ 的素数定理给出范数示性函数在素数上的加权平均为 $1/2$，余项为 $O_A(H/\log^A H)$，可取任意固定 $A$。代入 Landau–Selberg–Delange 均值定理得到 (3b)，首系数为 $\lim_{s\downarrow1}(s-1)^{1/2}D_k(s)/\Gamma(1/2)$。这里计的是所有截断的支撑并集，同一个固定 $z$ 还要满足 (2)。[Granville–Koukoulopoulos，定理 1](https://dms.umontreal.ca/~koukoulo/documents/publications/LSD.pdf)

实球计算给出 $0.6886417<C_k<0.6886429$。计算 $L(1,\chi_k)$ 用实偶二次特征的有限 Fourier 公式 $-229^{-1/2}\sum_{a=1}^{228}\chi_k(a)\log(2\sin(\pi a/229))$，它由正 Gauss 和 $\sqrt{229}$ 与对数级数作 Abel 极限得到。欧拉乘积截在 $P=300000$，余下乘积在 $1$ 与 $\prod_{m>P}(1-m^{-2})^{-1}=(P+1)/P$ 之间。程序还遍历 $2\le h\le8192$ 的所有系数跳跃，恰有 $2047$ 个下标在某个截断下非零，与独立的理想范数筛选逐个相同。[支撑并集、常数及严格尾界的复算](verify_cm229_coprime_smoothing.py)

## 单个零点在原余因子上的取值

设 $J>0$ 固定。对任意充分大的 $T$，存在一个取值于 $[0,1]$ 的单调截断 $\Phi_T$，在 $[0,1]$ 等于 $1$，在 $[2,\infty)$ 等于 $0$。设 $1/2\le\sigma\le1$，$s=\sigma+i\gamma$，$|\gamma|\le T$。对常数 $C_J$ 及所有 $H\ge C_J z(\log(2T)+|\gamma|)^2$，有

$$
 \sum_{h>z}c_{z,\nu}(h)h^{-s}\Phi_T(h/H)
 =1-U_{z,\nu}(s)L_\nu(s)+O_J(T^{-J}),\qquad
 U_{z,\nu}(s)=\sum_{d\le z}u_\nu(d)d^{-s}. \tag{4}
$$

**特别地**，在 $L_\nu(\rho)=0$、$\Re\rho\ge1/2$、$|\Im\rho|\le\log(2T)$ 时，取 $H=4C_J z\log^2(2T)$，(4) 的左端为 $1+O_J(T^{-J})$。它只用到

$$
 z<h\le 2H\ll_J T^{1/6}/\log(2T).
$$

当 $2n/r<hp\le2n$ 时，相应素数满足 $T^{5/6}\log(2T)\ll_J p\ll T^{5/6}\log^3(2T)$，因而全在 (1) 内。对同一个 $\Phi_T$，一整段低零点共用同一个 $H$。如果将 $H$ 放宽至 $T^{1/3}/8$，(4) 也覆盖 $|\Im\rho|\le c_J T^{1/12}\log^{3/2}(2T)$，取足够小的 $c_J>0$ 即可。

### 截断与移线的证明

取偶整数 $A=2\lceil\log(2T)\rceil$，$m=6A+20$，$\lambda=(\log 2)/m$。令 $V$ 为 $m$ 个独立的 $[0,\lambda]$ 均匀变量之和，取 $\Phi_T(y)=\Pr(e^V\ge y)$。其 Mellin 变换在 $\Re w>0$ 为

$$
 \widehat\Phi_T(w)=\frac1w\left(\frac{e^{\lambda w}-1}{\lambda w}\right)^m. \tag{5}
$$

该式在 $w=0$ 的留数为 $1$。把 $c_{z,\nu}=\mathbf1_{h=1}-u_{\nu,\le z}*a_\nu$ 代入左端，$\Phi_T(1/H)=1$。对每个 $d\le z$，对 $L_\nu(s+w)(H/d)^w\widehat\Phi_T(w)$ 作 Mellin 反演，从 $\Re w=2$ 移至 $-A$。$L_\nu$ 整，唯一留数为 $L_\nu(s)$。下面给出移线常数对 $A$ 的控制。

导子固定，两支的无穷因子分别是 $\Gamma_{\mathbb R}(s)^2$ 和 $\Gamma_{\mathbb R}(s+1)^2$。对任意偶整数 $A\ge2$、$1/2\le\sigma\le1$、实数 $v$，函数方程给出

$$
 |L_\nu(\sigma-A+iv)|\le [C(1+A+|v|)]^{2A+6}. \tag{6}
$$

具体地，从 $A=2$ 的伽马商出发，固定带上的 Stirling 界和有界高度处的连续性给出 $O((1+|v|)^6)$。每将 $A$ 增加 2，两处分子各增加一个一次因子，两处分母的倒数也各增加一个一次因子。导子和 $\pi$ 的因子至多指数增长，而反射线上的 $L_\nu(1-\sigma+A-iv)$ 绝对一致有界。这就给出 (6)。条带内由函数方程、绝对收敛和凸性有固定于该 $A$ 的多项式界，(5) 的衰减足以作上述移线。

在 $w=-A+iv$，(5) 内的平均因子模长至多为 $\min(1,2/(\lambda|w|))$。将积分分为 $|v|\le C_0m$ 和 $|v|>C_0m$，取 $C_0>4/\log 2$。前段的 $\int dv/|w|=O(1)$。后段令 $v_0=C_0m$，利用 $1+A+|\gamma|+|v|\le(1+A+|\gamma|+v_0)|v|/v_0$，多项式次数 $2A+6$ 小于 $m-1$，可积尾部也只付常数。因此移线后误差至多为

$$
 C z^{1-\sigma}\log^2(2z)\,[C(m+|\gamma|)]^6
 \left(\frac{C z(m+|\gamma|)^2}{H}\right)^A .
$$

这里用了 $|u_\nu(d)|\le d_2(d)$ 和 $\sum_{d\le z}|u_\nu(d)|d^{A-\sigma}\le z^{A+1-\sigma}(1+\log(2z))$。前因子对 $|\gamma|\le T$ 至多是 $T$ 的固定次幂，$m\ll\log(2T)$。按 $J$ 增大 $C_J$，括号的 $A$ 次幂即可给出 (4)。

在素数主项的频率 $s=1$，同一个截断给出不同的算术信息。用 Hu 的 $b=7/8$ 无零界，前篇对逆系数的 Perron 界与分部求和给出 $\sum_{d\le D}u_\nu(d)\ll_\epsilon D^{b+\epsilon}$。所以 $\sum_{d>z}u_\nu(d)/d\ll_\epsilon z^{-1/8+\epsilon}$，而 $\sum_{d=1}^\infty u_\nu(d)/d=1/L_\nu(1)$，由收敛半平面上的解析延拓得到。将 $s=1$ 代入 (4)，取 $H=4C_J z\log^2(2T)$，即有 $\sum c_{z,\nu}(h)\Phi_T(h/H)/h\ll_{\epsilon,J}z^{-1/8+\epsilon}+T^{-J}$。这给出了主项模式的实际衰减率。同一个权重在每个低零点上保留的系数由 (4) 精确给出。[逆系数的截断估计](cm229-type-ii-large-prime.md) · [Hu 的八分之七定理](https://hexagonmath.org/2610.00216v2)

### 对素数侧的含义

若把单支素数端点的一个显式公式模式 $x^\rho/\rho$ 代入 $2n/h$ 与 $2n/(rh)$ 两端，乘以 $c_{z,\nu}(h)\Phi_T(h/H)$ 再求和，所得恰为

$$
 (2n)^\rho\,\frac{1-r^{-\rho}}{\rho}
 \sum_{h>z}c_{z,\nu}(h)h^{-\rho}\Phi_T(h/H). \tag{7}
$$

对上面的低零点，(7) 中最后的因子趋于 $1$。区间的两个端点均已超过 $T^{2/3}$，可在原素数范围内作此检验。显式公式中该零点的负号与重数另乘在 (7) 外。这个计算确定了单一模式如何通过余因子卷积。所有零点、素数幂修正与其余余因子在实际整数窗口中的叠加，还要通过绝对值或正部估计来处理。对原素数和可将其中的高频误差无条件压到目标预算。[原互素对的整窗平滑定理](cm229-coprime-smoothing.md)
