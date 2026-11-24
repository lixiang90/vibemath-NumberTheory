---
title: "奇完全数只有 3 的偶指数例外时，超过一百万个素因子模 3 余 1"
author: "Codex"
date: "2026-11-24"
abstract: "设 3 整除奇完全数 N，除 3 与欧拉素因子外，每个素因子的指数均模 3 余 2。证明 N 有超过一百万个模 3 余 1 的不同素因子，若欧拉素因子异于 5，则超过十亿个。3 的指数相应至少为 1,000,002 和十亿。证明重新核验 Ochem–Zelinsky 的小素数排除树，并利用每个三次分圆值的素因子重数估计丰度。"
formula: '\#\{t\mid N:t\ \mathrm{prime},\ t\equiv1\bmod3\}>10^6,\qquad \pi\ne5\ \Longrightarrow\ v_3(N)\ge10^9'
---

# 奇完全数只有 3 的偶指数例外时，超过一百万个素因子模 3 余 1

设 $N$ 为奇完全数，$\sigma(N)=2N$，$\pi$ 为唯一具有奇指数的欧拉素因子。记 $e_q=v_q(N)$，$a=v_3(N)$。本文假设 $3\mid N$，并且

$$
3\mid e_q+1\qquad(q\mid N,\ q\notin\{3,\pi\},\ q\text{ 为素数}).
$$

这些偶指数可独立取 $2,8,14,\ldots$。McDaniel 的定理保证它们不可能连同 $a$ 一起都模 3 余 2，故 3 正是例外所在。Ochem 与 Zelinsky 在其他偶指数全部等于 2 的子类中证明了 $a\ge23\,000\,000\,000$。下面的界允许任意多个偶指数大于 2。[Fletcher–Nielsen–Ochem 对 McDaniel 定理的陈述与证明](https://www.lirmm.fr/~ochem/opn/OPNS_Adam_Pace.pdf) · [Ochem–Zelinsky，定理 1](https://arxiv.org/html/2607.19746v1)

**定理。** 令 $T$ 为整除 $N$ 且模 3 余 1 的素数集，$M$ 为异于 3 和 $\pi$ 的不同素因子个数。下表中的个数均按不同素数计。

| 条件 | $\lvert T\rvert$ | $a=v_3(N)$ | $M$ |
| --- | ---: | ---: | ---: |
| 所有上述 $N$ | $>10^6$ | $\ge1\,000\,002$ | $>1\,280\,000\,000$ |
| 另有 $\pi\ne5$ | $>10^9$ | $\ge10^9$ | $>4\cdot10^{15}$ |

因此，若 $a<10^9$，则 $\pi=5$，且存在素因子 $r>3\,000\,000$ 同时满足 $r\equiv1\bmod5$ 和 $e_r\equiv14\bmod30$。若只知道至多有一个偶指数不模 3 余 2，允许它出现在任何素数上，也总有 $3^{2570}\mid N$。后一推论将上表第一行与例外在 3 以外时已有的 $a\ge2570$ 合并。[此前的单例外定理](one-exception.md)

## 重新核验 300 万以内的因子树

记 $P=3\,000\,000$，$\Phi(q)=q^2+q+1$。对 $q\notin\{3,\pi\}$ 有 $\Phi(q)\mid\sigma(q^{e_q})$。$\Phi(q)$ 的每个大于 3 的素因子 $t$ 都整除 $N$，且 $q$ 模 $t$ 的乘法阶为 3，故 $t\in T$。对欧拉素因子则有 $\pi\equiv1\bmod4$，且 $\pi+1$ 的奇素因子都整除 $N$。

每个 $q\in R=\{q:q\mid N,\ q\text{ 为素数},\ q\notin\{3,\pi\}\}$ 都大于 $P$。若 $\pi\le P$，则 $\pi\equiv5\bmod12$。特别地，整个 $T$ 都在 $P$ 以上。归纳的起点为 5 只能作欧拉素因子且 $7\nmid N$，这两点已由适用于当前指数条件的丰度证书证明。Ochem–Zelinsky 的公开文件从 11 起为每个素数存了一棵整除树。其原始问题中另有 $\pi\equiv1\bmod12$ 的约束，故在这里需要逐棵检验欧拉素因子的选择。[起点的证明](one-exception.md) · [Ochem–Zelinsky 第 4 节](https://arxiv.org/html/2607.19746v1#S4) · [原始树文件](https://www.lirmm.fr/~ochem/opn/restrained/tree.txt)

按根素数 $r$ 从小到大归纳。取该树中充当边起点的所有素数为顶点。依次假定其中每个模 4 余 1 的顶点是 $\pi$，再加上 $\pi$ 在顶点集以外的情形。若 $r\equiv5\bmod12$，要证的仅是 $r$ 不能作非欧拉素因子，此时略去 $\pi=r$。从 $r$ 出发，遇到指定的 $\pi$ 就沿 $q+1$ 的边，其他顶点沿 $\Phi(q)$ 的边。若到达较小素数 $s\ge5$ 且 $s\not\equiv5\bmod12$，归纳假设给出矛盾。若 $s\equiv5\bmod12$ 且已经指定树内某一顶点为 $\pi$，也有矛盾，因为树内顶点均至少为 $r>s$。当 $\pi$ 在树外时，只用第一种终点。

验证器检查了 $216\,812$ 个连续的根素数和 $286\,751$ 种欧拉素因子选择，每种都达到上述终点。程序重算每条所列整除关系，直接筛出 $P$ 以内的素数，并对 $31\,043$ 个更大的不同顶点调用 FLINT 的素性判定。这些较大顶点也用 PARI 的 $\mathrm{isprime}$ 独立核过。所用树文件按 LF 换行的 SHA-256 为 `57b69c549df9fab23ecfda1de66d0794edc1d4a2fbaa33bb882f631351a49eb7`。[因子树验证器](verify_role_tree.py) · [作者的数据与生成程序](https://www.lirmm.fr/~ochem/opn/restrained/) · [PARI 素性判定](https://pari.math.u-bordeaux.fr/dochtml/html-stable/Arithmetic_functions.html#isprime)

## 丰度与分圆值的素因子

令 $Q=\{q\in R:q\equiv5\bmod6\}$。对 $q\in Q$，$\Phi(q)$ 的素因子全在 $T$。每个 $t\in T\setminus\{\pi\}$ 又满足 $3\mid\Phi(t)\mid\sigma(t^{e_t})$。比较 $v_3(\sigma(N))=a$ 得

$$
|T|\le a+1.
$$

当 $\pi=5$ 时，每个 $t\in T$ 都贡献一个因子 3，而 $6\mid\sigma(5^{e_5})$ 再贡献一个，故 $a\ge |T|+1$。写 $I(m)=\sigma(m)/m$。由 $I(q^{e_q})<q/(q-1)$，完全性给出

$$
\sum_{q\in R}\frac1{q-1}>
\log\frac{4(\pi-1)}{3\pi}\ \ge\ \log(16/15).
$$

若 $\pi\ne5$，小素数排除给出 $\pi\ge17$，于是最后的下界可提高到 $\log(64/51)$。

对 $i=1,5$，在 $P<s\le L$ 中保留模 6 余 $i$ 的素数 $s$，并要求 $\Phi(s)$ 没有介于 5 与 $P$ 之间的素因子。用 $C_i(L)$ 记候选数，$U_i(L)$ 记将每个 $1/(s-1)$ 向上舍入到 $10^{-18}$ 后的和。$R$ 在相应区间中的每个素因子都是候选。分段筛给出下表。

| $i$ | $L$ | $C_i(L)$ | $10^{18}U_i(L)$ |
| ---: | ---: | ---: | ---: |
| 1 | $3\cdot10^8$ | 909,158 | 14,719,406,968,152,615 |
| 5 | $10^{10}$ | 26,199,184 | 24,422,990,910,195,945 |
| 1 | $3\cdot10^{11}$ | 676,766,694 | 32,152,938,648,615,242 |
| 5 | $10^{11}$ | 235,925,513 | 29,883,045,391,089,345 |

第一份筛在奇数中用模素数的非平凡三次单位根删去两个类。第二份筛只枚举 $s\equiv i\bmod6$，用 $4\Phi(s)=(2s+1)^2+3$ 和模平方根求出要删的类。两份程序所得四行整数完全相同。[三次单位根筛](../restrained-exponent/sieve_ordinary_prefix.rs) · [独立的模平方根筛](sieve_support_classes.rs)

反设 $|T|\le K$。对一个截点 $Y>P$ 置

$$
S=U_1(Y)+\frac{K-C_1(Y)}{Y-1}+\frac1{P-1},\qquad H=S\log Y.
$$

则 $\sum_{t\in T}1/(t-1)\le S$，且 $\sum_{t\in T}(\log t)/(t-1)\le H$。对第一式，可先给至多 $K$ 个素数各分配 $1/(Y-1)$，再加上所有 $Y$ 以内候选超过该权重的部分。最后一项覆盖可能不满足候选条件的欧拉素因子。对第二式同样分配 $(\log Y)/(Y-1)$，使用 $(\log t)/(t-1)$ 递减，并在 $t\le Y$ 时以 $\log Y$ 上界分子。这就把 $R$ 中模 3 余 1 的素因子的权重以及 $Q$ 在 $X$ 以内的权重合计限制为 $S+U_5(X)$。

## 按素数幂计算中间区间的丰度

取两个不交的素数集 $A,B\subset[5,P]$。对 $p\in A\cup B$，把模 $p$ 的 0 和 $\Phi$ 的根列为禁用类。其个数 $\nu_p$ 在 $p\equiv1\bmod3$ 时为 3，其余为 1。记 $w_p=\nu_p/p$，$J_p(n)$ 为命中禁用类的指示函数。筛权

$$
\mathbf1_{n\equiv5\bmod6}\prod_{p\in A}(1-J_p)
\left(1-\sum_{p\in B}J_p+\sum_{\{p,r\}\subset B}J_pJ_r\right)
$$

处处非负，在 $Q$ 上等于 1。括号中若命中 $m$ 个素数，其值为 $1-m+\binom m2$。以下的配对求和均遍历不同素数的无序对。记

$$
\rho=\frac16\prod_{p\in A}(1-w_p)
\left(1-\sum_{p\in B}w_p+\sum_{\{p,r\}\subset B}w_pw_r\right),\qquad
D=\prod_{p\in A}(1+\nu_p)
\left(1+\sum_{p\in B}\nu_p+\sum_{\{p,r\}\subset B}\nu_p\nu_r\right).
$$

对 $t\in T$ 和 $h\ge1$，$\Phi$ 模 $t^h$ 有两个根，因为模 $t$ 的两个根都为单根。$t>P$，故中国剩余定理可将每个根与筛权中每个单项式的剩余类组合。一个等差数列在 $(u,x]$ 中的计数与其平均值相差至多 1。因此，对满足 $t^h\mid\Phi(q)$ 的 $q\in Q\cap(u,x]$，其个数至多为

$$
2\rho(x-u)/t^h+2D .
$$

把这个不等式乘以 $\log t$，对 $t,h$ 求和。$\Phi(q)<(x+1)^2$，所以每个 $t$ 至多出现 $2\log(x+1)/\log t$ 个幂次。利用 $\sum_{h\ge1}t^{-h}=1/(t-1)$ 得到

$$
\sum_{\substack{q\in Q\\u<q\le x}}\log\Phi(q)
\le 2\rho H(x-u)+4DK\log(x+1).
$$

由于 $\log\Phi(q)>2\log q$，对右端以 $f(x)=1/(2(x-1)\log x)$ 分部求和。线性主项给出 $\rho H\int_u^v dx/((x-1)\log x)$。误差项至多为 $4DK[f(u)\log(u+1)+\int_u^v f(x)/(x+1)\,dx]$。用 $\int_u^\infty dx/((x-1)(x+1))\le1/(u-1)$，在 $v>u\ge10^{10}$ 时便有

$$
\sum_{\substack{q\in Q\\u<q\le v}}\frac1{q-1}
\le \rho H\,\frac{u}{u-1}\log\frac{\log v}{\log u}
+\frac{11DK}{5(u-1)}. \tag{1}
$$

## 最远端的光滑数估计与数值比较

对 $q>V$ 使用 $1/(q-1)<2/\sqrt{\Phi(q)}$。$\Phi(q)$ 严格递增，且其素因子全在 $T$。于是

$$
\sum_{\substack{q\in Q\\q>V}}\frac1{q-1}
\le 2 V^{-1/8}\prod_{t\in T}(1-t^{-7/16})^{-1}
\le 2\exp\left(\frac K{674}-\frac{\log V}{8}\right).
$$

第一步在所有 $T$ 上光滑的正整数上求和，使用了 $\Phi(q)^{1/16}>V^{1/8}$。第二步用 $P^7>675^{16}$ 和 $-\log(1-z)\le z/(1-z)$。

当 $K=10^6$ 时取 $(X,Y,V)=(10^{10},3\cdot10^8,10^{6000})$。当 $K=10^9$ 时取 $(X,Y,V)=(10^{11},3\cdot10^{11},10^{6\,000\,000})$。下表逐段给出式 (1) 中的 $(d,i,j)$。区间左端为 $u=10^d$，右端为下一项的左端，末项止于 $V$。将从 5 起的素数递增编号，$A$ 取前 $i$ 项，$B$ 取随后 $j-i$ 项。最大的筛素数为 9973。

| $K=10^6$ 时的 $(d,i,j)$ | $K=10^9$ 时的 $(d,i,j)$ |
| --- | --- |
| $(10,0,2)$ | $(11,0,0)$ |
| $(13,2,12)$ | $(14,1,6)$ |
| $(20,12,100)$ | $(20,9,50)$ |
| $(40,50,1227)$ | $(40,40,800)$ |
| $(100,160,1227)$ | $(100,160,1227)$ |
| $(400,800,1227)$ | $(400,800,1227)$ |
| $(2000,1227,1227)$ | $(2000,1227,1227)$ |

两组参数都满足 $(\log V)/8>K/674+16$，故最远端贡献小于 $10^{-6}$。将 $S+U_5(X)$、七个区间与最远端相加，$K=10^6$ 时给出上界 $0.063739537810818335$，小于 $\log(16/15)>0.06453852113757117$。$K=10^9$ 时的上界为 $0.211112814554479005$，小于 $\log(64/51)>0.22705745063534608$。两项反设分别与完全性矛盾。再用 $\pi=5$ 时的 $a\ge |T|+1$ 和 $\pi\ne5$ 时的 $a\ge |T|-1$，得到定理中的两个指数下界。[全部有理数区间运算](verify_ternary_support.py)

若 $a<10^9$，现已知 $\pi=5$。因为 $5\mid\sigma(N)$ 且 $5\nmid\sigma(5^{e_5})$，存在某个 $r\ne5$ 使 $5\mid\sigma(r^{e_r})$。$e_r+1$ 为奇数，而模 5 的乘法群阶为 4，故 $r\equiv1\bmod5$，进而 $5\mid e_r+1$。这也排除了 $r=3$。又有 $3\mid e_r+1$，两式合并即 $e_r\equiv14\bmod30$。

## 不同素因子总数

在两种剩余类合并的筛中，取 $x=2\cdot10^{11}$，有 $C=917\,142\,659$ 个候选，其权重上界为 $U=0.062696136715648189$。每个 $q>x$ 的权重至多为 $1/x$，每个 $q\le x$ 的候选权重至少为 $1/x$。因此全部 $R$ 的权重至多为 $U+(M-C)/x$。若 $M\le1\,280\,000\,000$，该值至多为 $0.064510423420648189$，小于 $\log(16/15)$。[合并两类的筛表](../restrained-exponent/ordinary-prefix-200000000000.json)

若 $\pi\ne5$，取 $x=2\cdot10^{12}$，有限段权重至多为 $0.072401332985083531$。在 $x<q\le z=5\cdot10^{17}$ 上使用已有的素数倒数和估计。它仅要求 $q$ 为素数并且 $\Phi(q)$ 避开 $7,13,\ldots,103$ 中所用的 12 个筛素数，故同样适用于 $R$。写 $\delta=64959667/169646400,\ c=31/60$，该段上界为

$$
\frac{x}{x-1}\left[\delta\log\frac{\log z}{\log x}
+c\left(\frac1{\log x}-\frac1{\log z}
+\frac3{2\log^2x}+\frac1{2\log^2z}\right)\right].
$$

对 $q>z$ 加上 $M/z$。若 $M\le4\cdot10^{15}$，总和至多为 $0.226398627629550785$，小于 $\log(64/51)$。这证明了最后一个计数下界。[有限段筛表](../restrained-exponent/ordinary-prefix-2000000000000.json) · [中间段估计与其筛权](../restrained-exponent/smooth-indices.md) · [等差数列素数估计的原文](https://arxiv.org/abs/1802.00085)

## 复算

在仓库根目录执行以下命令。因子树使用作者的 `tree.txt`，已在 Python 3.12 与 python-flint 0.9.0 上复算。数值验证只使用 Python 标准库。下列命令在 64 位系统上编译两份 Rust 分段筛，并分别对照四行整数。省略两个筛的路径时，数值验证器直接核验表中整数所给出的全部不等式。

~~~sh
python -B -O OPN/exponent-congruences/verify_role_tree.py tree.txt
rustc -O OPN/restrained-exponent/sieve_ordinary_prefix.rs -o sieve_ordinary_prefix
rustc -O OPN/exponent-congruences/sieve_support_classes.rs -o sieve_support_classes
python -I -B -O OPN/exponent-congruences/verify_ternary_support.py --sieve ./sieve_ordinary_prefix --reference-sieve ./sieve_support_classes
~~~
