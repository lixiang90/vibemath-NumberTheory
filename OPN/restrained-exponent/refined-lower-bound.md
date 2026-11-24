---
title: "受限奇完全数中 3 的指数大于 3.28 × 10²²"
author: "Codex"
date: "2026-11-24"
abstract: "设恰有一个非欧拉素因子的指数大于 2。将二步筛的第二层截点提高到 300 万，结合一个 116 项的筛多项式，证明 3 的指数大于 3.28 × 10²²。计算逐项使用等差数列素数计数的显式常数，筛多项式在全部 8192 个布尔输入上用整数核验。"
formula: 'v_3(N)>3.28\cdot10^{22}'
---

# 受限奇完全数中 3 的指数大于 3.28 × 10²²

**定理。** 设 $N$ 是奇完全数，恰有一个非欧拉素因子的指数大于 2。则 $v_3(N)>3.28\cdot10^{22}$。这里对 $v_3(N)+1$ 的素因子大小没有附加条件。

后续将本文的有限筛与平滑显式公式合用，把下界提高到 $3.8\cdot10^{26}$，并排除了两个光滑指数加一。[平滑筛证明](smoothed-exponent.md)

依 Ochem 与 Zelinsky 的结构定理，写 $N=3^a p^e\prod_{s\in S}s^2$，其中 $p$ 是欧拉素数，$p$ 及所有 $s\in S$ 均大于 $3\,000\,000$。Stormwhip 的固定版本计算用分段筛和二步因子链给出 $a>1.14\cdot10^{22}$。我们先复现其二步筛数据，再把第二层检查素数的截点由 1 万提高到 10 万和 300 万。与更细的前段、116 项筛多项式及逐模数误差合用，三个截点给出的下界分别为 $2.13,2.67,3.28$ 乘以 $10^{22}$。[结构定理](https://arxiv.org/abs/2607.19746v1) · [Stormwhip 的原计算与界](https://github.com/stormwhip/OddPerfectNumberBounds/tree/dc6e8353293d8922a8ade8ffa625f4aeb2c72ae6)

## 二步筛能扣去多少权重

令 $A=328\cdot10^{20},\ B=50\,290\,000\,000\,000,\ Z=290A$。若 $a\le A$，由 $|S|\le2a+1$ 及 $w(t)=1/(t-1)-1/(Z-1)$ 得

$$
\sum_{s\in S}\frac1{s-1}
\le\frac{2A+1}{Z-1}+\sum_{\substack{s\in S\\s\le Z}} w(s).
$$

$w$ 在 $t\le Z$ 时非负。对 $s\le B$，先计入全截点筛至 $2\cdot10^{12}$ 的本地结果，再取原作者分段表的后 4829 行。合计 $H=188\,818\,064\,757$ 个候选，其 $\lceil10^{18}/(s-1)\rceil$ 之和为 $V=84\,752\,764\,168\,054\,679$。[计数界](../ordinary-components/counting.md) · [两个区间的计算](shared-indices.md) · [前段数据](ordinary-prefix-2000000000000.json) · [分段表](https://github.com/stormwhip/OddPerfectNumberBounds/blob/dc6e8353293d8922a8ade8ffa625f4aeb2c72ae6/sieve.tsv)

在 $(3\,000\,000,10^9]$ 上，原作者二步筛的初筛与本地全截点筛相同。用本地 Rust 实现重筛得到同样的 $5\,819\,954$ 个候选与整数和 $36\,514\,781\,544\,466\,812$。对每个普通素因子候选 $s$，$\Phi_3(s)$ 除去可能的因子 3 后的任一素因子 $t$ 都必须整除 $N$。若可排除 $t$ 的普通身份及欧拉身份，即可排除 $s$。若两个不同的 $t$ 都必须是欧拉素数，也可排除 $s$。剩下仅能强迫一个欧拉素数的候选按这个素数分组，实际至多保留一组。固定源码在原截点 1 万上逐行复现作者表。将截点增至 300 万后，五个互斥类别如下。[二步筛源码](https://github.com/stormwhip/OddPerfectNumberBounds/blob/dc6e8353293d8922a8ade8ffa625f4aeb2c72ae6/chain.cpp) · [原截点表](https://github.com/stormwhip/OddPerfectNumberBounds/blob/dc6e8353293d8922a8ade8ffa625f4aeb2c72ae6/chain.tsv) · [300 万截点表](chain-3000000.tsv) · [独立初筛](sieve_ordinary_prefix.rs)

| 原表类别 | 候选数 | 整数和 |
| --- | ---: | ---: |
| backward | 1,545,845 | 6,378,044,940,351,857 |
| forward_ordinary | 1,810,488 | 12,760,183,939,302,973 |
| forward_special | 1,580,046 | 11,182,719,501,787,552 |
| two_specials | 7,012 | 27,499,217,131,414 |
| special_conditional | 236,348 | 1,606,732,934,734,099 |

最后一类有 235918 组，最大组整数和为 $332\,769\,178\,653$。先减去五类全部整数和，再加回最大组，得到保留候选的整数和上界 $V'=52\,797\,916\,403\,925\,437$。对负项 $-1/(Z-1)$，则保守地从 $H$ 减去五类的全部候选数，得到所用候选超集的个数下界 $H'=188\,812\,885\,018$。由此 $(3\,000\,000,B]$ 上的 $w$ 和至多为 $V'/10^{18}-H'/(Z-1)$。表中另有一行统计分圆商为素数且已排除的情形，与前述类别有交叠。

例如，普通素因子候选 $s=3004279$ 满足 $s^2+s+1=3t$，其中 $t=3008565104707$ 为素数。$t\equiv3\bmod4$，只能取普通身份；但 $158341\mid t^2+t+1$，而 $158341$ 是被结构定理排除的素数。它是 $(t^2+t+1)/3$ 的最小素因子。截点 10 万时这条链仍未闭合，在 300 万截点上便排除了 $s$。

## 一个可以穷举核验的筛多项式

取 $D_0=\{7,13,19,31,37,43,61,67,73,79,97,103,109\}$。对素数 $r$，定义 $x_h(r)=1$ 当且仅当 $h\mid r^2+r+1$。所有普通素因子都满足全部 $x_h=0$，因为 $h\le 3\,000\,000$ 且 $h>3$。以平方自由乘积 $d$ 编号，选取多项式

$$
G(x)=\sum_d \lambda_d \prod_{h\mid d}x_h,\qquad \lambda_d\in\{-1,1\}.
$$

附件按 $[d,\lambda_d,c_\pi]$ 给出全部 116 项，$d=1$ 表示常数项，最大 $d=13699$。对每个布尔输入直接求和，可验 $G(0)=1$、全部 $8192$ 个值均非负。这保证 $G$ 逐素数控制普通素因子的指示函数。[完整系数及逐模数误差常数](bmor-majorant.json) · [整数核验](verify_refined_exponent.py)

对 $d>1$，每个 $h\mid d$ 恰有两个非平凡三次单位根，故其单项式计数恰为模 $d$ 的 $2^{\omega(d)}$ 个既约类的素数计数之和。把常数项按模 3 的两个既约类计数；这一表示在 $(B,Z]$ 覆盖每个素数。令 $F(t)$ 为这些素数计数的带符号和。对 $t\ge B$，由 Bennett、Martin、O'Bryant 与 Rechnitzer 的定理 1.3 得

$$
|F(t)-\delta\operatorname{Li}(t)|\le C\,t/(\log t)^2,
$$

其中

$$
\begin{aligned}
\delta &=1+\sum_{d>1}\lambda_d \frac{2^{\omega(d)}}{\varphi(d)}
=\frac{2492150407}{6616209600},\\
C &=2c_\pi(3)+\sum_{d>1}2^{\omega(d)} c_\pi(d)
=\frac{2759763}{5000000}.
\end{aligned}
$$

这里逐模数使用作者给出的 $c_\pi(d)$；对 $d\le10^5$，定理的起点均不超过 $8\cdot10^9<B$。空乘积的误差相应取 $2c_\pi(3)$。[定理 1.3、式 (1.11) 及附录 A](https://arxiv.org/abs/1802.00085) · [作者的逐模数常数表](https://nt.math.ubc.ca/BeMaObRe/c-psi-theta-pi/c_all_rounded.txt)

## 对数丰度比较

记 $b=\log B,\ u=\log Z,\ T=B/(B-1)$。把 $w\,dF$ 分部求和。因 $w(Z)=0$，误差的端点项至多 $CT/b^2$，积分项至多 $CT^2(1/b-1/u)$。主项使用 $\int_B^Z dt/\log t\ge Z/u-B/b+(Z-B)/u^2$，并用 $1/(t-1)\le T/t$。加入欧拉素数的对数丰度上界 $1/2999999$，得

$$
\begin{aligned}
\log\frac{\sigma(N/3^a)}{N/3^a}\le{}&
\frac{V'}{10^{18}}+\frac{2A+1-H'}{Z-1}+\frac1{2999999}\\
&+\delta\left[T\log(u/b)-
\frac{Z/u-B/b+(Z-B)/u^2}{Z-1}\right]\\
&+C\left[\frac{T}{b^2}+T^2\left(\frac1b-\frac1u\right)\right]\\
\le{}&0.287678023229577646 .
\end{aligned}
$$

完全性要求此值严格大于 $\log(4/3)>0.287682072451780927$，因为 $\sigma(3^a)/3^a<3/2$。矛盾证明 $a>A$。在附加 $2850$-光滑条件的两族中，最小的指数加一为 $2383^7=436383475091811838100527$，仍高于这个下界，14 个候选均保留。[候选族的证明](shared-indices.md)

三个截点所得的严格比较列于下表。权重上界已经向上舍入。

| 第二层截点 | $A/10^{22}$ | 对数丰度上界 |
| ---: | ---: | ---: |
| 10,000 | 2.13 | 0.287657033235868910 |
| 100,000 | 2.67 | 0.287677876327082635 |
| 3,000,000 | 3.28 | 0.287678023229577646 |

在本目录运行首条命令即可核对最强一行的二步表、五个代表见证、全部布尔值与有理数界。第二条另核验 4943 行原始分段表及 BMOR 作者表中的每一个所用常数。用 `--chain chain.tsv` 可回到原作者的 1 万截点。用原作者的 `chain.cpp` 从头进行两次新计算所需的命令也列在下方。输出顺序及代表见证均由程序确定，线程数可调整。[10 万截点输出](chain-100000.tsv) · [300 万截点输出](chain-3000000.tsv) · [原作者源码](https://github.com/stormwhip/OddPerfectNumberBounds/blob/dc6e8353293d8922a8ade8ffa625f4aeb2c72ae6/chain.cpp)

~~~sh
python -B -O verify_refined_exponent.py
python -B -O verify_refined_exponent.py --sieve sieve.tsv --bmor c_all_rounded.txt
python -B -O verify_refined_exponent.py --chain chain-100000.tsv
g++ -std=c++17 -O3 -pthread chain.cpp -o chain
./chain 1000000000 8 3000000 100000 > chain-100000.tsv
./chain 1000000000 10 3000000 3000000 > chain-3000000.tsv
~~~
