---
title: "奇完全数的五次分圆筛与三个模 5 特征的卷积"
author: "Codex"
date: "2026-11-24"
abstract: "用三个非主狄利克雷特征的调和卷积给出五次分圆筛分母的显式二次对数下界。对任意奇完全数，指数加一被 5 整除且大于 Y 的素因子，其倒数权重之和小于 740/(log Y−9)，适用于 log Y≥30。更精细地处理筛中的例外素数，还可用于收紧 3 或 5 整除指数加一时的最小素因子上界。"
formula: 'G_{0,5}(Q)\ge\frac6{55}\left(\log Q-\frac7{15}\log40\right)^2'
---

# 奇完全数的五次分圆筛与三个模 5 特征的卷积

设 $N$ 是奇完全数，$R_5=\{s\mid N:s\text{ 为素数},\ 5\mid v_s(N)+1\}$。本文证明，只要 $\log Y\ge30$，就有

$$
\sum_{\substack{s\in R_5\\s>Y}}\frac1{s-1}<\frac{740}{\log Y-9}.
$$

这里 $\log$ 为自然对数。若再知道 $N$ 的最小素因子很大，同一个筛分母可以给出小得多的常数。[在混合指数问题中的应用](mixed-least-prime.md) · [一般奇指标的筛维数](prime-tails.md)

## 完整筛的分母

取 $\rho_0(q)=5$ 当素数 $q\equiv1\bmod5$，其余处取 $\rho_0(q)=1$。这对应于禁用 0 和全部非平凡五次单位根。令

$$
G_{0,5}(Q)=\sum_{d\le Q}\mu^2(d)\prod_{q\mid d}\frac{\rho_0(q)}{q-\rho_0(q)} .
$$

**分母引理。** 对 $Q\ge40$ 有

$$
G_{0,5}(Q)\ge\frac6{55}\left(\log Q-\frac7{15}\log40\right)^2 .
$$

把每个局部权重展开为几何级数。所有 $n\le Q$ 的不同素因子之积均不超过 $Q$，故 $G_{0,5}(Q)\ge\sum_{n\le Q}w(n)/n$，其中 $w$ 是满足 $w(q^k)=\rho_0(q)^k$ 的乘性函数。再取乘性函数 $b$，其狄利克雷级数为 $\zeta(s)\prod_{q\equiv1\bmod5}(1-q^{-s})^{-4}$。在这些素数上 $b(q^k)=\binom{k+4}{4}\le5^k$，其余处 $b(q^k)=1$，所以 $0\le b\le w$。

令 $\chi$ 是模 5 且 $\chi(2)=i$ 的四次特征，$\psi=\chi^2$。对 5 的倍数，特征值均取零。用 $*$ 表示狄利克雷卷积，$\mathbf1(n)=1$，置

$$
D=\psi*\chi*\bar\chi,\qquad h=\mathbf1*\mathbf1*D .
$$

比较欧拉因子得 $h=b*c$，其中 $c\ge0$，且

$$
\sum_{n\ge1}\frac{c(n)}n=
\frac54\prod_{q\equiv2,3\bmod5}(1-q^{-4})^{-1}
             \prod_{q\equiv4\bmod5}(1-q^{-2})^{-2}<\frac{11}{8}.
$$

上界可直接验算。对 $q\le100$ 乘有理因子，其余因子的对数和至多为 $2/(3\cdot100^3)+3/(100\cdot9)$。这里用了 $\sum_{n>100}n^{-4}\le1/(3\cdot100^3)$，以及大于 100 且模 5 余 4 的素数形如 $10k+9$、$k\ge10$。对数尾和为 $u$ 时用 $\exp(u)\le1/(1-u)$。因此正项卷积给出 $\sum_{n\le Q}h(n)/n\le(11/8)G_{0,5}(Q)$。

## 特征卷积的调和部分和

记 $H_f(y)=\sum_{n\le y}f(n)/n$。我们将证明

$$
H_D(y)\ge\frac4{25}\quad(y\ge1),\qquad
H_D(y)\ge\frac3{10}\quad(y\ge40).
$$

先处理无限区间。写 $L_f=\sum_{n\ge1} f(n)/n$，此记号先用于三个特征。$\chi$ 在每个周期内的普通部分和取值于 $0,1,1+i$，两两距离至多为 $\sqrt2$。$\psi$ 的任一区间和绝对值至多为 2。由分部求和，$y\ge1$ 时

$$
|H_\chi(y)-L_\chi|\le\frac{3}{2y},\qquad
|H_\psi(y)-L_\psi|\le\frac2y .
$$

截到一万项，逐项向下舍入实部和虚部，分别减去这个尾界，可得 $L_\psi|L_\chi|^2>339/1000$，三个所用的实下界均为正。

对 $UV=x$，从 $ab\le x$ 的双曲线区域按 $a\le U$ 或 $b\le V$ 拆分，扣除重叠矩形。由此得到精确恒等式

$$
\begin{aligned}
H_{f*g}(x)-L_fL_g
={}&\sum_{a\le U}\frac{f(a)}a\big(H_g(x/a)-L_g\big)
  +\sum_{b\le V}\frac{g(b)}b\big(H_f(x/b)-L_f\big)\\
 &-\big(H_f(U)-L_f\big)\big(H_g(V)-L_g\big).
\end{aligned}
$$

等式对任意指定的常数 $L_f,L_g$ 成立。对 $B=\chi*\bar\chi$ 取 $U=V=\sqrt y$，由上面的尾界得到

$$
|H_B(y)-|L_\chi|^2|\le3y^{-1/2}+\frac94 y^{-1}.
$$

对 $(n,5)=1$，有 $B(n)=\bar\chi(n)(\mathbf1*\psi)(n)$。乘性函数 $\mathbf1*\psi$ 非负，且其总和至多为 $\lfloor t\rfloor$，因为 $\sum_{n\le t}(\mathbf1*\psi)(n)=\sum_{a\le t}\sum_{b\le t/a}\psi(b)$，每个内和至多为 1。又 $5\mid n$ 时 $B(n)=0$，所以 $\sum_{b\le t}|B(b)|\le t$。

再次应用双曲线恒等式，取 $f=\psi,\ g=B,\ U=x^{1/3},\ V=x^{2/3}$。第一项的绝对值至多为 $6\sqrt{U/x}+9U/(4x)$，用了 $\sum_{a\le U}a^{-1/2}\le2\sqrt U$。第二项至多为 $2V/x$，矩形项至多为 $6/(U\sqrt V)+9/(2UV)$。因此

$$
\left|H_D(x)-L_\psi\,\lvert L_\chi\rvert^2\right|
\le8x^{-1/3}+\frac{33}{4}x^{-2/3}+\frac92x^{-1}.
$$

在 $x\ge10^7$ 时，用 $x^{1/3}>215$ 已可得到 $H_D(x)>3/10$。

有限段的系数也完全是整数。按 $q\bmod5$ 分为 $\{1\},\{4\},\{2,3\}$ 三类，$D(q^k)$ 的生成函数依次为

$$
(1-t)^{-3},\qquad (1-t)^{-1}(1+t)^{-2},\qquad(1-t)/(1-t^4).
$$

且 $D(5^k)=0$ 对 $k\ge1$ 成立。由最小素因子筛重建 $D(n)$。依次累加 $\lfloor10^{14}D(n)/n\rfloor$ 后，在 $1\le n<40$ 上最小值为 $16666666666666$，出现在 $n=3$。在 $40\le n\le10^7$ 上最小值为 $30434989347753$，出现在 $n=143$。向下舍入对负项同样有效。一份独立的 Rust 程序由约数和 $\sum_{d\mid n}\psi(d)$ 及 $\bar\chi(n)$ 求出 $B(n)$，再与 $\psi$ 作卷积，得到相同的极小值和位置。[乘性筛](quintic_character.py) · [独立约数和验证](check_quintic_convolution.rs) · [全部常数验证](verify_mixed_prime_bound.py)

## 从调和下界回到大筛

对 $t\ge1$ 有 $\sum_{uv\le t}1/(uv)\ge(\log t)^2/2$，先对 $v$ 用调和和的积分下界，再对递减的 $\log(t/u)/u$ 积分。将 $h=\mathbf1*\mathbf1*D$ 展开，得到

$$
\sum_{n\le Q}\frac{h(n)}n
\ge \frac{2}{25}(\log Q)^2+\frac7{100}(\log(Q/40))^2
\ge \frac3{20}\left(\log Q-\frac7{15}\log40\right)^2.
$$

末步为平方函数的凸性。乘以 $8/11$ 即得分母引理。

## 利用小筛素数的高次幂

混合指数问题还要用到一个加强式。同时令 $\ell=3,5$，取

$$
(c_3,a_3)=\left(\frac7{50},\frac{\log74}6\right),\qquad
(c_5,a_5)=\left(\frac6{55},\frac7{15}\log40\right).
$$

上述五次证明和已有的三次证明，实际都给出 $\sum_{n\le t}b_\ell(n)/n\ge c_\ell(\log t-a_\ell)_+^2$ 对 $t>0$ 成立，其中 $y_+=\max(y,0)$。这是因为每个调和下界在起始端点 $T$ 产生的是 $(\log(t/T))_+^2/2$，函数 $y\mapsto y_+^2$ 仍为凸函数。这里 $b_\ell$ 的级数为 $\zeta(s)\prod_{q\equiv1\bmod\ell}(1-q^{-s})^{-(\ell-1)}$。[三次情形的卷积证明](ternary-tail-constants.md)

在 $q\equiv1\bmod\ell$ 处，$w_\ell(q^k)=\ell^k$ 的局部生成函数除以 $b_\ell$ 的局部生成函数为 $(1-u)^\ell/(1-\ell u)$。该商有非负系数。取形式幂级数的对数，所得 $\sum_{k\ge1}(\ell^k-\ell)u^k/k$ 系数非负，取指数即得。只保留 $q\le100$ 的这类局部商，定义相应的非负乘性函数 $f_\ell$。于是 $w_\ell\ge b_\ell*f_\ell$ 逐项成立。

写 $L_\ell=\sum f_\ell(n)/n$，按 $f_\ell(n)/(nL_\ell)$ 给 $n$ 赋概率，令 $\mu_\ell$ 为 $\log n$ 的期望。由局部级数及其导数，

$$
L_\ell=\prod_{\substack{q\le100\\q\equiv1\bmod\ell}}
 \frac{(1-1/q)^\ell}{1-\ell/q},\qquad
\mu_\ell=\sum_{\substack{q\le100\\q\equiv1\bmod\ell}}
 \frac{\ell(\ell-1)\log q}{(q-\ell)(q-1)} .
$$

对 $\sum f_\ell(n)c_\ell(\log(Q/n)-a_\ell)_+^2/n$ 用凸性，得到

$$
G_{0,\ell}(Q)\ge c_\ell L_\ell(\log Q-a_\ell-\mu_\ell)_+^2 .
$$

有限个有理因子的乘积和对数的上界给出 $L_3>1.149,\ \mu_3<0.77,\ L_5>1.165,\ \mu_5<1$。这四个数只涉及 100 以内的素数。[有理数核验](verify_mixed_prime_bound.py)

## 对任意奇完全数的显式尾界

若 $s\in R_5$，则 $\Phi_5(s)=1+s+s^2+s^3+s^4\mid2N$。对不整除 $N$ 的筛素数可以禁用全部非平凡五次单位根。在整除 $N$ 的筛素数上只禁用 0。完整分母在每个 $q\mid N,\ q\equiv1\bmod5$ 处损失的乘法因子至多为 $(q-1)/(q-5)$，因为固定其余权重后分母为 $A+g(q)B$，且 $0\le B\le A$。

这些损失的乘积小于 20。事实上，令 $U_q=(q-1)/((q-5)(1+1/q)^4)$。对 $q\ge11$，利用 $-\log(1-u)\le u+u^2/(2(1-u))$ 和 $\log(1+u)\ge u-u^2/2$，得 $\log U_q\le25/q^2$。又 $U_q\ge1$，因为 $(1+1/q)^4\le(1-4/q)^{-1}\le(q-1)/(q-5)$。取所有 $q\le1000,\ q\equiv1\bmod5$ 的精确乘积。余下素数形如 $10k+1$，$k\ge100$，它们的 $\log U_q$ 之和至多为 $25/(100\cdot99)$。两段合并给出 $\prod U_q<5/4$，而 $\prod_{q\mid N}(1+1/q)\le2$。

取 $x\ge\exp(30)$，$Q=\lfloor\sqrt x/12\rfloor\ge\sqrt x/13$。先计入至多 $Q$ 个小素数，再用大筛及分母引理。由 $2\log13+(14/15)\log40<9$，得到

$$
\#\{s\in R_5:s\le x\}
\le Q+20(1+1/144)\frac{4\cdot55}{6}\frac{x}{(\log x-9)^2}.
$$

这里 $Q(\log x-9)^2/x\le21^2/(12\exp(15))$。对 $1/(s-1)$ 分部求和，积分中由 $(Y/(Y-1))^2<(1+10^{-6})^2$ 控制分母。所得常数小于 740，证明了文首的倒数尾界。[Fletcher–Nielsen–Ochem，命题 1 的大筛](https://www.lirmm.fr/~ochem/opn/OPNS_Adam_Pace.pdf)

对所有 $Y\ge3$，还可以统一取上界 $1058/\log Y$。在 $\log Y\ge30$ 时用 $740\cdot30/21<1058$，其余区间用 $\sum_{q\mid N}1/(q-1)\le3/2$。
