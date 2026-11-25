---
title: "每个对数平方短区间内含陈氏素数加数的哥德巴赫表示"
date: "2026-11-25"
abstract: "由普通 RH 的 Selberg 均方估计、陈氏定理和四个线性式的上界筛，证明每个充分大的对数平方长度区间都有正比例的偶数，各有正确数量级下界的两素数表示，其中指定加数 p 的邻项 p+2 至多有两个大素因子。对真正的孪生素数子集，给出与其实际二倍区间计数相联系的一致误差界。"
---

# 每个对数平方短区间内含陈氏素数加数的哥德巴赫表示

存在绝对常数 $C,c>0$，使每个充分大的整数 $M$ 之后、长度 $C\log^2 M$ 的区间中，至少有 $c\log^2 M$ 个偶数 $n$ 各有 $\gg M/\log^3 M$ 种表示 $n=p+r$。两项都是真正的素数，$p$ 与 $M$ 的大小可比，并且 $p+2$ 是素数，或是两个都大于 $p^{3/11}$ 的素数之积。这使同一个加数 $r$ 还给出 $n+2=(p+2)+r$。本文使用普通黎曼猜想（RH）的经典推论。Alpöge–Claude 归档中对应的形式化定理为 `riemannHypothesis`。[原始定理](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Main.lean) · [Lean 与 Mathlib 的验证说明](https://leanprover-community.github.io/blog/posts/riemann-hypothesis/)

Goldston 证明了普通两素数表示的短区间平均定理。关于每个短区间存在一个普通哥德巴赫数，已有对区间长度给出数值常数的定理。文库已有文章把两个加数限制在指定比例窗口内。这里对其中一个素数再加上差二的因数条件。主要困难是这些加数在长度 $h$ 的移动窗口内预期只有 $h/\log^2 M$ 个。第四节还给出适用于任意孪生素数子集的误差估计，所需输入是这批素数在一个二倍区间中的**实际个数**。[Goldston，定理与推论 1、2](https://doi.org/10.1017/S001708950000937X) · [普通哥德巴赫数的显式区间定理](https://doi.org/10.1007/s40993-026-00767-6) · [文库中的平衡表示](goldbach-short-intervals.md) · [文库中的陈氏素数分解](twin-prime-scale-averages.md)

## 1. 定理及局部限制

取充分大的整数 $X$，令 $L=\log X$、$I=(X,2X]$。记 $\mathcal C_X$ 为 $I$ 中满足上述 $p+2$ 条件的素数，$\mathcal T_X=\{p\in I:p,\ p+2\text{ 均为素数}\}$。每种计数均按指定的 $p$ 计一次。Green–Tao 所引定理 6.1 给出绝对常数 $c_0>0$ 使 $|\mathcal C_X|\ge c_0X/L^2$。该文脚注 3 单列了把 Iwaniec 的区间换成二倍区间这一点。以下不需要陈氏素数的渐近式。[Green–Tao，第 6 节及定理 6.1](https://www.numdam.org/article/JTNB_2006__18_1_147_0.pdf#page=19)

对 $\mathcal A\subseteq\mathcal C_X$，置

$$
 R_{\mathcal A}(n)=\#\{p\in\mathcal A:n-p\text{ 为素数}\}.
$$

**定理 1.** 设 $0<\delta\le1$，$|\mathcal A|\ge\delta X/L^2$。在 RH 下，存在绝对常数 $C_1,c_1,c_2>0$，使对**每个**整数 $4X\le M\le5X$ 和 $C_1\delta^{-2}L^2\le h\le\sqrt X$ 的每个整数 $h$，区间 $(M,M+2h)$ 至少有 $c_1\delta^2h$ 个偶数 $n$ 满足

$$
 R_{\mathcal A}(n)\ \ge\ c_2\delta\,X/L^3. \tag{1}
$$

取 $\mathcal A=\mathcal C_X$、$\delta=\min(c_0,1)$、$X=\lfloor M/4\rfloor$，并把 $h$ 向上取整，即得开头的结论。这里 $\gg$ 的常数与 $M$ 无关，本文未将它们数值化。充分大时 $\mathcal C_X$ 中的 $p$ 都模 $3$ 余 $2$。因此 (1) 所得偶数只可能模 $3$ 余 $0$ 或 $1$，两类合并计数。

## 2. 四个线性式给出的局部二阶矩

写 $w(b)=\prod_{\ell\mid b,\ \ell>7}(1+8/\ell)$，积只取素数，$b$ 为正整数。以下三个上界使用同一个固定维数筛估计，

$$
\begin{aligned}
 |\mathcal C_X|&\ll X/L^2,\\
 \#\{p,p+d\in\mathcal C_X\}
 &\ll \frac{X}{L^4}w(d)w(|d-2|)w(d+2)
       &&(1\le d\le\sqrt X,\ d\ne2),\\
 R_{\mathcal C_X}(n)&\ll\frac{X}{L^3}w(n)w(n+2)
       &&(4X\le n\le 5X+2\sqrt X).
\end{aligned} \tag{2}
$$

具体地，三行要筛去小素因子的线性式分别为 $v,v+2$，为 $v,v+2,v+d,v+d+2$，以及为 $v,v+2,n-v$。取 $z=X^\epsilon$，$\epsilon>0$ 为充分小的绝对常数，尤其小于 $3/11$。所有被计入的行都在筛后留下。若某个素数的根覆盖所有剩余类，筛后计数为零。否则记 $\nu(\ell)$ 为不同根数，按中国剩余定理乘性延拓至无平方因子整数。对 $k$ 个线性式，$\nu(\ell)\le k\le4$，筛余项在模 $d$ 处为 $O(\nu(d))$。上界筛基本引理取水平 $D=\lfloor X^{1/2}\rfloor$，权的绝对值至多 $1$，给出上界

$$
 \ll X\prod_{\ell\le z}(1-\nu(\ell)/\ell)
       +\sum_{d\le D}\mu^2(d)k^{\omega(d)}
 \ll X\prod_{\ell\le z}(1-\nu(\ell)/\ell)+X^{1/2}L^{k-1}.
$$

这里 $\omega(d)$ 计不同素因子。筛的维数条件由 $\nu(\ell)\le k$ 和 Mertens 乘积一致成立，有限个 $\ell\le7$ 只改变常数。对 $\ell>7$，通常 $\nu(\ell)=k$，偏离只可能发生于两个根重合时，修正比 $(1-\nu(\ell)/\ell)/(1-k/\ell)\le1+8/\ell$。四个根重合所涉及的大素数都整除 $d(d-2)(d+2)$，第三行相应为 $n(n+2)$，故得 (2)。遗漏的位移 $d=2$ 恰使两式重合，此时 $p,p+2,p+4$ 都将没有小素因子，其中必有一个被 $3$ 整除，故贡献为零。[Ford，《Sieve Methods》，定理 3.6](https://ford126.web.illinois.edu/sieve2023.pdf)

设 $A_h(y)=\#(\mathcal C_X\cap(y,y+h])$。两个长 $h$ 的区间在端点相距 $d$ 时相交长度为 $(h-d)_+$，于是

$$
 \int_{\mathbb R}A_h(y)^2\,dy
 =h|\mathcal C_X|+2\sum_{1\le d<h}(h-d)\#\{p,p+d\in\mathcal C_X\}
 \ll X\left(\frac h{L^2}+\frac{h^2}{L^4}\right). \tag{3}
$$

为验证最后一步，对每个固定 $k$，把 $w(b)^k$ 展成 $\sum_{d\mid b}a_k(d)$，其中 $a_k$ 非负、支撑于无平方因子整数，且 $a_k(\ell)=O_k(1/\ell)$。故 $\sum_{b\le t}w(b)^k\le t\prod_\ell(1+a_k(\ell)/\ell)\ll_k t$。对 (2) 第二行的三个因子应用 Hölder 不等式及 $k=3$ 的界，便得它们对 $d\ne2$ 的和为 $O(h)$。

## 3. 用 RH 估计另一个加数

令 $\vartheta(u)=\sum_{r\le u,\ r\ {\rm prime}}\log r$，$E_h(u)=\vartheta(u+h)-\vartheta(u)-h$。Selberg 在 RH 下的均方估计为 $\int_1^{6X}|E_h(u)|^2du\ll XhL^2$，这里取 $1\le h\le\sqrt X$。这也是本文调用 RH 的唯一一步。[Goldston，式 (2)、(4)](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/F73ED944781B38B244B14ADD72AD7B24/S001708950000937Xa.pdf/linniks-theorem-on-goldbach-numbers-in-short-intervals.pdf)

定义 $k_h(t)=(h-|t-h|)_+$。对每个 $p\in\mathcal C_X$，考虑归一化的搭档权重

$$
 B_p(M,h)=\frac1{h^2}\sum_{r\ {\rm prime}}\log r\, k_h(p+r-M).
$$

直接按实区间交长换序，有

$$
 B_p(M,h)-1
 =\frac1{h^2}\int_{p-h}^{p} E_h(M-y)\,dy .
$$

事实上 $r\in(M-y,M-y+h]$ 与 $y\in[p-h,p)$ 的交长正是 $k_h(p+r-M)$，而常数 $h$ 的积分为 $h^2$。在所用支撑上 $1<M-y<6X$。先取绝对值再对 $p$ 求和，利用 (3)、Cauchy 不等式和 Selberg 界，得到对 $L^2\le h\le\sqrt X$ 一致成立的估计

$$
 \boxed{\displaystyle
 \sum_{p\in\mathcal C_X}|B_p(M,h)-1|
 \le \frac1{h^2}\int A_h(y)|E_h(M-y)|\,dy
 \ll \frac X{L^2}\frac{L}{\sqrt h}.} \tag{4}
$$

这一步保留了每个 $p$ 的绝对误差。因此对任意子集 $\mathcal A$，$\sum_{p\in\mathcal A}B_p(M,h)=|\mathcal A|+O(X/(L\sqrt h))$。当 $h\ge C_1\delta^{-2}L^2$ 时，取 $C_1$ 充分大，误差至多为 $|\mathcal A|/2$。因为 $k_h\le h$ 且所计 $r\asymp X$，这给出

$$
 \sum_{M<n<M+2h}R_{\mathcal A}(n)\gg \delta\, h X/L^3 . \tag{5}
$$

为控制这些表示在个别 $n$ 上的集中，置 $z_0=\log(7X)$。对 $b\le7X$，$\ell>z_0$ 的因子对 $\log w(b)$ 贡献 $O(1/\log z_0)$，因为其个数至多 $\log(7X)/\log z_0$。将剩余乘积的四次方展开，在长度为 $2h+O(1)$ 的区间中计入每个 $d$ 的倍数，得到

$$
 \sum_{M<n<M+2h+3}w(n)^4
 \ll h\prod_\ell(1+a_4(\ell)/\ell)
        +\prod_{\ell\le z_0}(1+a_4(\ell))
 \ll h+(\log z_0)^{O(1)}\ll h .
$$

所以 (2) 第三行及 Cauchy 不等式给 $\sum R_{\mathcal A}(n)^2\ll h(X/L^3)^2$。在 (5) 中删去所有 $R_{\mathcal A}(n)<c_2\delta X/L^3$ 的项，取 $c_2$ 足够小，余下的总和仍 $\gg\delta hX/L^3$。再用一次 Cauchy 不等式就有至少 $c_1\delta^2 h$ 个目标。两个素数均大于 $2$，所以这些目标全为偶数。

## 4. 真孪生对的定量结论

把 (4) 用于 $\mathcal T_X$，即得对每个允许的 $M$ 及 $L^2\le h\le\sqrt X$ 成立的公式

$$
 \frac1{h^2}\sum_{\substack{X<p\le2X\\p,\ p+2,\ r\ {\rm prime}}}
       k_h(p+r-M)\log r
 =|\mathcal T_X|+O\left(\frac{X}{L\sqrt h}\right). \tag{6}
$$

其中 $r$ 遍历所有素数。式 (6) 的误差界无需关于孪生素数数量的假设。若在某一尺度 $|\mathcal T_X|\ge\delta X/L^2$，定理 1 就在该尺度的**每一个**目标短区间给出至少 $c_1\delta^2 h$ 个 $n$，各有至少 $c_2\delta X/L^3$ 组共用同一素数 $r$ 的表示

$$
 n=p+r,\qquad n+2=(p+2)+r,\qquad p,\ p+2,\ r\ {\rm prime}.
$$

是否有固定 $\delta>0$ 使这个输入在无穷多个二倍区间成立，仍是开放问题。对陈氏素数全集的已证下界允许所有质量来自半素数邻项。

## 5. 有限核对

下表取 $M=4X$，把 $\log r$ 换成 $1$，用整数运算比较按表示和按窗口交长计算的三角和。程序也分别按窗口和按位移核对 (3) 的精确等式，并检验 $d=2$ 和模 $3$ 的局部排除。有限例子检验端点与重数，定理中的常数来自上述估计。[复算程序](chen-summands-short-intervals/verify.py)

| $X$ | $h$ | $\#\mathcal C_X$ | $\#\mathcal T_X$ | 全部陈氏素数的三角和 | 真孪生子集的三角和 |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1,000 | 19 | 43 | 26 | 2,052 | 1,296 |
| 10,000 | 31 | 269 | 137 | 24,850 | 12,826 |
| 30,000 | 45 | 646 | 344 | 115,376 | 62,286 |
