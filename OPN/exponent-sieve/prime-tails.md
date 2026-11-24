---
title: "一般奇完全数中，哪些指数能贡献素因子的倒数和？"
date: "2026-11-24"
author: "Codex"
abstract: "对固定奇数 h≥3，按 h 整除指数加一筛选奇完全数的素因子，证明其计数上界中对数的幂次可取 h 的正约数个数。若最小素因子 P 很大，至少有常数倍 P²/log P 个素因子的指数超过任意预设阈值。对固定的有限组奇合数，还给出避开这些整除条件的素因子数的渐近下界。"
formula: '\#\{q\le x:q\mid N,\ h\mid v_q(N)+1\}\ll_h x(\log x)^{-\tau(h)}'
---

# 一般奇完全数中，哪些指数能贡献素因子的倒数和？

设 $N$ 是任意奇完全数，即 $\sigma(N)=2N$。对素数 $q\mid N$ 写 $e_q=v_q(N)$，令 $P$ 为 $N$ 的最小素因子。固定一个奇整数 $h\ge3$，以 $\tau(h)$ 表示其正约数个数，取

$$
R_h(N)=\{q\mid N:q\ {\rm prime},\ h\mid e_q+1\}.
$$

**定理 1。** 对 $x,X\ge3$，

$$
\#\{q\in R_h(N):q\le x\}\ll_h \frac{x}{(\log x)^{\tau(h)}},\qquad
\sum_{\substack{q\in R_h(N)\\q>X}}\frac1{q-1}\ll_h\frac1{(\log X)^{\tau(h)-1}}.
$$

常数可有效计算，与 $N$ 及其他指数无关。例如，$h=3$ 时它同时控制指数 $2,8,14,\ldots$ 的所有非欧拉素因子，此时倒数和的上界可取 $188/\log X$。取 $h=9$ 或 $15$ 时，倒数尾和分别为 $O((\log X)^{-2})$ 或 $O((\log X)^{-3})$。文中 $\log$ 表示自然对数。[三次筛常数的推导](ternary-tail-constants.md)

**定理 2。** 对任意固定的有限组奇素数 $\mathcal L$，存在 $c_{\mathcal L}>0$ 和 $P_{\mathcal L}$，使得 $P>P_{\mathcal L}$ 时，至少有 $c_{\mathcal L}P^2/\log P$ 个不同素因子 $q\mid N$ 满足

$$
\gcd\left(e_q+1,\prod_{\ell\in\mathcal L}\ell\right)=1 .
$$

因此，给定任意固定的 $E$，当 $P$ 充分大时，至少有常数倍 $P^2/\log P$ 个素因子的指数大于 $E$。常数可以依赖于 $E$。

若 $\mathcal H$ 是固定的有限组奇合数，则对每个 $\varepsilon>0$，当 $P$ 充分大时，有

$$
\#\{q\mid N:q\ {\rm prime},\ h\nmid e_q+1\ \text{对所有 }h\in\mathcal H\}
\ge(1/2-\varepsilon)\,\frac{P^2}{\log P}.
$$

**定理 3。** 若至多有 $K$ 个素因子的指数加一与 $\prod_{\ell\in\mathcal L}\ell$ 互素，则对每个固定 $j$，$N$ 的第 $j$ 小素因子有一个只依赖 $\mathcal L,K,j$ 的可计算上界，前提是 $\omega(N)\ge j$。下文给出递推式。

Yamada 以及 Fletcher、Nielsen 与 Ochem 用分圆因子和筛法给出了有限组指数条件下的最小素因子上界。Zelinsky 研究了奇完全数所有素因子的倒数和与最小素因子的关系。这里把倒数和按 $h\mid e_q+1$ 分组，同一组里同时利用 $h$ 的各个约数。筛掉哪些同余类可以随 $N$ 变化，控制损失所用的丰度上界对所有 $N$ 有效。[Yamada，定理 1.1](https://arxiv.org/abs/1609.08795) · [Fletcher–Nielsen–Ochem，命题 1](https://www.lirmm.fr/~ochem/opn/OPNS_Adam_Pace.pdf) · [Zelinsky](https://arxiv.org/abs/2402.14234)

对单个奇素数 $\ell$，此前按 $\Phi_\ell(q)$ 的素因子个数分组，得到了 $P^{1/(\ell-1)}\log P$ 量级的下界，并在 $P>10^8$ 处给出了 $\ell=3,5$ 的具体数值。定理 2 将充分大的 $P$ 下的量级提高到 $P^2/\log P$，并允许同时避开有限组素数。两种方法的显式阈值不同。[按分圆值分组的结果](../exponent-congruences/least-prime.md)

## 丰度控制筛去的剩余类

写 $I(m)=\sigma(m)/m$。乘法性给出

$$
1+\sum_{q\mid N}\frac1q\le\prod_{q\mid N}\left(1+\frac1q\right)
\le \prod_{q\mid N}I(q^{e_q})=2 ,
$$

这里及以下的 $q\mid N$ 在素数上取值。若 $r\in R_h(N)$，则 $1+r+\cdots+r^{h-1}$ 整除 $\sigma(r^{e_r})$，故其每个奇素因子都整除 $N$。

用素数 $q$ 筛 $r$，先禁用剩余类 0，允许 $r=q$ 这个例外。对 $q\nmid N$ 再禁用模 $q$ 的全部非平凡 $h$ 次单位根。每个这样的根都使上述几何和为零，且 $q=2$ 没有这样的根。记禁用类个数为 $\rho(q)$。对 $q\nmid N$ 有 $\rho(q)=\gcd(h,q-1)$，因为模 $q$ 的乘法群是循环群。欧拉函数恒等式和固定模数的素数定理给出

$$
\sum_{q\le u}\frac{\gcd(h,q-1)}q
=\sum_{d\mid h}\varphi(d)\sum_{\substack{q\le u\\q\equiv1\bmod d}}\frac1q
=\tau(h)\log\log u+O_h(1).
$$

在 $q\mid N$ 处改取 $\rho(q)=1$，对上式的总影响至多为 $(h-1)\sum_{q\mid N}1/q\le h-1$。因此真实筛的误差对 $N$ 一致。固定模数的显式素数定理给出了可计算的误差常数。[Bennett–Martin–O'Bryant–Rechnitzer，定理 1.2](https://arxiv.org/abs/1802.00085)

取 $Q=\lfloor\sqrt{x}\rfloor$。先计入至多 $Q$ 个较小素数，其余素数由大筛给出上界 $(x+Q^2)/G(Q)$，其中

$$
G(Q)=\sum_{d\le Q}\mu^2(d)\prod_{q\mid d}\frac{\rho(q)}{q-\rho(q)} .
$$

为给出一致的下界，令 $u=Q^{1/(4h)}$。展开 $V(u)=\prod_{q\le u}(1-\rho(q)/q)^{-1}$，将各平方自由整数 $d$ 的项归一为概率，则 $\log d$ 的期望为 $\sum_{q\le u}\rho(q)\log q/q$。写 $m=\lfloor u\rfloor$，由 $\lfloor m/q\rfloor\ge m/(2q)$ 和 $\sum_{q\le m}\lfloor m/q\rfloor\log q\le\log(m!)$，得 $\sum_{q\le u}\log q/q\le2\log u$。因此该期望至多为 $2h\log u=\tfrac12\log Q$。至少一半权重来自 $d\le Q$，故 $G(Q)\ge V(u)/2$。又 $\log V(u)\ge\sum_{q\le u}\rho(q)/q$，从而 $G(Q)\gg_h(\log x)^{\tau(h)}$。这给出计数界。分部求和时对 $1/(t(\log t)^{\tau(h)})$ 从 $X$ 积分到无穷，即得定理 1 的尾界。[大筛不等式](https://www.lirmm.fr/~ochem/opn/OPNS_Adam_Pace.pdf)

## 最小素因子较大时必须出现许多高指数

写 $R=\bigcup_{\ell\in\mathcal L}R_\ell(N)$，$T$ 为其余素因子的集合。由定理 1 并对有限个 $\ell$ 求和，存在常数 $C=C_{\mathcal L}$ 使 $\sum_{q\in R}1/(q-1)\le C/\log P$ 对充分大的 $P$ 成立，这里可从 $P-1$ 起用尾界。因 $\log I(q^{e_q})<1/(q-1)$，有

$$
\sum_{q\in T}\frac1{q-1}>\log2-\frac C{\log P}.
$$

取 $\eta=\exp(-4C-4)$、$Y=\eta P^2$。当 $P$ 充分大时，$Y>P$。素数倒数和的 Mertens 估计给出

$$
\sum_{P\le r\le Y}\frac1{r-1}
=\log\log Y-\log\log P+O\big((\log P)^{-2}\big)
=\log2+\frac{\log\eta}{2\log P}+O_\eta\big((\log P)^{-2}\big).
$$

第一式用的是素数倒数和减去 $\log\log x$ 及其常数后的 $O(1/\log^2x)$ 误差，从 $1/r$ 换为 $1/(r-1)$ 的总误差为 $O(1/P)$。这个区间内全部素数的和因而小于 $T$ 所需的和。若 $\#T$ 至多为 $\pi(Y)-\pi(P-1)$，在所有不小于 $P$ 的素数中选取最小的这么多个就能取得最大的和，矛盾。素数定理给出 $\pi(Y)-\pi(P-1)\gg_{\mathcal L} P^2/\log P$。再取 $\mathcal L$ 为不超过 $E+1$ 的全部奇素数，每个不超过 $E$ 的正偶指数都被某个 $\ell$ 覆盖。单独扣除至多一个欧拉素因子，即得高指数的结论。[Rosser–Schoenfeld，素数倒数和估计](https://doi.org/10.1215/ijm/1255631807)

若筛去的指标 $h\in\mathcal H$ 均为合数，则 $\tau(h)\ge3$，其联合集合的倒数和为 $O_{\mathcal H}(1/\log^2 P)$。对任意固定 $0<\eta<1$ 重取 $Y=\eta P^2$。上式中的负项 $(\log\eta)/(2\log P)$ 绝对值最终超过所有 $O(1/\log^2P)$ 的误差，故余下的素因子数大于 $\pi(\eta P^2)-\pi(P-1)\sim(\eta/2)P^2/\log P$。让 $\eta$ 任意接近 1，即得所述 $1/2-\varepsilon$ 的结论。

## 部分素数幂的丰度与 2 的间距

从 $N$ 取 $t\ge1$ 个完整素数幂 $q_i^{e_i}$，每个底数至多为 $B$，并且至少遗下一个素因子。令 $M=\prod_{i=1}^t I(q_i^{e_i})<2$。则

$$
2-M\ge\frac1{(4t B^t)^{2^t}} .
$$

用 $q_i/(q_i-1)$ 暂代全部有限几何和，再逐个恢复精确值。记当前乘积为 $J$，并用 $D$ 表示其未约分的分母，初值至多为 $B^t$。只要尚有一个暂代因子，$J$ 的二进赋值至多为 0。事实上，每个暂代因子贡献至多 $-1$，每个偶指数的精确因子贡献 0，唯一可能的欧拉分量在精确取值时贡献 1，因为欧拉素数及其指数均模 4 余 1。所以此时 $J\ne2$。全部恢复后 $J=M<2$。

若 $J<2$，则 $M\le J$ 且 $2-M\ge1/D$。若 $J>2$，则 $J-2\ge1/D$。记尚用暂代值的下标集为 $\mathcal D$。有

$$
J/M=\prod_{i\in\mathcal D}(1-q_i^{-(e_i+1)})^{-1}\le\zeta(2)<2,
\qquad J-M\le J\sum_{i\in\mathcal D}q_i^{-(e_i+1)} .
$$

因 $J<4$，至少一个尚未恢复的因子满足 $q_i^{e_i+1}<4tD$。恢复它后，新分母至多为 $Dq_i^{e_i}<4tD^2$。至多恢复 $t$ 次便遇到 $J<2$，各步分母均不超过 $(4tB^t)^{2^t}$。这就证明了间距估计。

## 有限个指数例外下的递推

取 $C_{\mathcal L}$ 为定理 1 中联合集合 $R$ 的尾界常数，使该界对 $X\ge3$ 成立。较小的 $X$ 可用 $\sum_{q\mid N}1/(q-1)\le3/2$ 统一处理。设至多有 $K$ 个素因子在 $R$ 以外。取 $B_0=3$，对 $j\ge1$ 置 $t=j-1$，

$$
H_j=\big(\max(1,4t)B_{j-1}^t\big)^{2^t},\qquad
B_j=\exp\!\big(3(C_{\mathcal L}+K+1)H_j\big).
$$

设 $N$ 至少有 $j$ 个不同素因子。归纳时取前 $t$ 个完整分量，其丰度乘积记作 $M$。当 $t=0$ 时 $M=1$，否则这些底数均不超过 $B_{j-1}$。上面的间距引理给出 $\log(2/M)\ge(2-M)/2\ge1/(2H_j)$。

若第 $j$ 小素因子大于 $B_j$，剩余所有分量的对数丰度至多为 $C_{\mathcal L}/\log B_j+K/(B_j-1)<1/(2H_j)$。其余分量的丰度乘积又必须等于 $2/M$。矛盾证明第 $j$ 小素因子至多为 $B_j$。例如，若至多有 $k$ 个非欧拉素因子的指数不满足 $e_q\equiv2\bmod3$，可取 $\mathcal L=\{3\}$、$K=k+1$、$C_{\{3\}}=188$。允许的指数 $2,8,14,\ldots$ 与这些素因子的数目均无预设上限。
