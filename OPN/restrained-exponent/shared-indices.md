---
title: "受限奇完全数的 2850-光滑指数加一只余两族"
date: "2026-11-24"
author: "Codex"
abstract: "若受限奇完全数中 3 的指数加一的全部素因子不超过 2850，则其值只能是 2383 的幂，或 2579 乘以 2383 的幂。新的乘法阶与素性证书将先前 173 个候选缩到 23 个。用已公开的分段筛数据与带权的丰度估计，可再去掉 9 个，仅余 14 个整数，最大者为 2383 的 15 次幂。"
formula: 'a+1\in\{2383^j\}\cup\{2579\cdot2383^j\}'
---

# 受限奇完全数的 2850-光滑指数加一只余两族

设 $N$ 是恰有一个非欧拉素因子的指数大于 2 的奇完全数。根据 Ochem 与 Zelinsky 的结构定理，写 $N=3^a p^e\prod_{s\in S}s^2$，其中 $p\equiv1\bmod12$ 是欧拉素因子，$p$ 与所有 $s\in S$ 都大于 $3\,000\,000$。令 $n=a+1$。[结构定理与小素数排除](https://arxiv.org/html/2607.19746v1#S4)

**定理。** 若 $n$ 的每个素因子都不超过 2850，则 $n$ 属于下列 14 个整数。

$$
\begin{aligned}
&2383^j &&(7\le j\le15),\\
&2579\cdot2383^j &&(6\le j\le10).
\end{aligned}
$$

因此，若 $2383^6\nmid n$，或 $n$ 含有 $2383,2579$ 之外的素因子，则 $n$ 的最大素因子至少为 2851。还有两项无需假定光滑性的排除，即 $2707\nmid n$ 和 $2383^{16}\nmid n$。这 14 个数是 3 的指数加一的必要候选，每个对应的奇完全数是否存在仍是开放问题。[完整升序表](shared-candidates.json) · [先前的 173 个指数候选](smooth-indices.md)

这个子类中 $n$ 的上界由 $2383^{27}$ 缩至 $2383^{15}$。由 $\omega(N)\le2n+1$ 及 Nielsen 的上界，还得到 $N<2^{\,4^{2\cdot2383^{15}+1}}$。[素因子计数](../ordinary-components/counting.md) · [Nielsen 的上界](https://mathdept.byu.edu/~pace/upper_bound_web.pdf)

## 共用因子给出的整除障碍

若 $d\mid n$ 且素数 $r\ge5$ 满足 $\operatorname{ord}_r(3)=d$，则 $r\mid\sigma(3^a)\mid2N$。若 $r$ 是普通素因子，$r^2+r+1$ 的每个素因子都整除 $N$。若为欧拉素因子，则 $r\equiv1\bmod12$，且 $r+1$ 的每个奇素因子都整除 $N$。对各个允许的身份给出一条最终强迫已排除素数的链，即证明 $d\nmid n$。

例如，$r=636986691468700189849$ 为素数，模 $r$ 的 3 的乘法阶是 2707。若 $r$ 是欧拉素数，则 $5\mid r+1$。若 $r$ 是普通素数，则 $u=1878120406567177$ 整除 $r^2+r+1$。若 $u$ 是欧拉素数，$1949\mid u+1$，而 1949 作为普通素数会强迫 19。若 $u$ 是普通素数，则 $2917\mid u^2+u+1$。19 与 2917 的排除都已在先前的证书中闭合。由此得 $2707\nmid n$。这一个障碍就排除起始表中的 32 个候选。

限制 2383 幂次的证书只需一步。令 $d=2383^{16}$，取素数

$$
r=670617484191031556660627309144686105309784705364469906348834403.
$$

有 $r=2\cdot310070321\,d+1$，$3^d\equiv1\bmod r$ 而 $3^{d/2383}\not\equiv1\bmod r$。因 $r\equiv11\bmod12$，若整除 $N$，只能是普通素因子。又 $r\equiv2\bmod7$，因而 $7\mid r^2+r+1$，这与 $7\nmid N$ 矛盾。故 $2383^{16}\nmid n$。

合并见证后，新增的互不整除的障碍为

$$
\begin{gathered}
2707,\quad1109\cdot2441,\quad1789^2,\quad1789\cdot2383,\quad
2239\cdot2521,\quad971^2\cdot2797,\\
1571\cdot1949\cdot2579,\quad2441\cdot2557^2,\quad2797^3,\quad
2383^2\cdot2797^2,\quad2383^5\cdot2797,\quad2383^{11}\cdot2579,\quad2383^{16}.
\end{gathered}
$$

将旧图的已排除节点视为可复用的结论，新证明图增加 30 个素数节点。验证器递归证明素性，在每个根检查 $3^d\equiv1\bmod r$ 及所有素数 $t\mid d$ 对应的 $3^{d/t}\not\equiv1\bmod r$，然后检验每种允许身份的整除关系。先前的 173 个数中，23 个通过这些新障碍。直接重枚举全部 2850-光滑数得到相同的集合。[新增因子链](shared-chains.json) · [素性证书](shared-primality.json) · [原图与枚举证明](smooth-indices.md)

## 用带权的分段和排除较小候选

设 $A=261\cdot10^{18}$，暂设 $a\le A$。取

$$
X=2\cdot10^{12},\qquad B=50\,290\,000\,000\,000,\qquad Z=262A .
$$

已证 $\#S\le2a+1$。从每个普通素因子的权重中拿出 $1/(Z-1)$，因 $s>Z$ 时余下权重为负，得到

$$
\sum_{s\in S}\frac1{s-1}
\le\frac{2A+1}{Z-1}+\sum_{\substack{s\in S\\s\le Z}}w(s),
\qquad w(t)=\frac1{t-1}-\frac1{Z-1}.
$$

在 $s\le B$ 的部分，$w(s)$ 非负，可以用所有分段筛候选作上界。至 $X$ 使用本地公开的全截点筛，$(X,B]$ 使用 Stormwhip 表中的连续 4829 行。两个区间的列值如下，末列是对每个候选的 $10^{18}/(s-1)$ 向上取整后求和。[分量个数的证明](../ordinary-components/counting.md) · [本地有限筛及区间数据](ordinary-prefix-2000000000000.json) · [Stormwhip 的固定版本分段表](https://github.com/stormwhip/OddPerfectNumberBounds/blob/dc6e8353293d8922a8ade8ffa625f4aeb2c72ae6/sieve.tsv)

| 区间 | 候选数 | 整数和 |
| --- | ---: | ---: |
| $(3\,000\,000,X]$ | 8,391,832,364 | 72,401,332,985,083,531 |
| $(X,B]$ | 180,426,232,393 | 12,351,431,182,971,148 |
| 合计 | 188,818,064,757 | 84,752,764,168,054,679 |

记合计为 $H$ 个候选、整数和 $V$，则有限部分至多为 $V/10^{18}-H/(Z-1)$。本次重新汇总并校验了 $(X,B]$ 的逐行数据及其来源哈希，区间内的逐个素数筛选采用原作者的计算。下面的带权方法取自 Stormwhip 的验证程序。该作者另外报告用二步筛得到 $a>1.14\cdot10^{22}$。此处只用上表的直接筛即可去掉这 9 个候选。[原作者的计算与界](https://github.com/stormwhip/OddPerfectNumberBounds/tree/dc6e8353293d8922a8ade8ffa625f4aeb2c72ae6) · [带权估计源代码](https://github.com/stormwhip/OddPerfectNumberBounds/blob/dc6e8353293d8922a8ade8ffa625f4aeb2c72ae6/verify.py)

对 $(B,Z]$，沿用先前论文中逐项核对的 97 项筛多项式。其相应的等差数列素数计数线性组合 $F(t)$ 在这一区间逐点控制普通素因子的计数增量。对所有 $t\ge B$，Bennett、Martin、O'Bryant 与 Rechnitzer 的定理 1.3 给出

$$
|F(t)-\delta\,\operatorname{Li}(t)|\le C\frac{t}{(\log t)^2},
\qquad \delta=\frac{64959667}{169646400},\quad C=\frac{434}{840}.
$$

这里每个模数均介于 3 与 10000 之间。常数项按模 3 的两个互素类计入，共涉及按系数绝对值计数的 434 个剩余类。[筛多项式与系数](smooth-indices.md) · [等差数列素数误差，定理 1.3](https://arxiv.org/abs/1802.00085)

写 $b=\log B,\ z=\log Z,\ T=B/(B-1)$。对 $w\,dF$ 分部求和，因 $w(Z)=0$，误差至多为

$$
C\left(\frac{T}{b^2}+T^2\left(\frac1b-\frac1z\right)\right).
$$

其中第一项来自左端点，积分项用了 $t/(t-1)^2\le T^2/t$。主项至多为

$$
\delta\left(T\log(z/b)-\frac{Z/z-B/b+(Z-B)/z^2}{Z-1}\right).
$$

这是因为 $1/(t-1)\le T/t$，而 $\int_B^Z dt/\log t\ge Z/z-B/b+(Z-B)/z^2$，后者由分部积分得到。合并各段并加入欧拉素数至多 $1/2999999$ 的对数丰度贡献，所得上界为 $0.287621941187407121$。由于 $\sigma(3^a)/3^a<3/2$，完全性要求这部分对数丰度严格大于 $\log(4/3)>0.287682072451780927$，矛盾。因此 $a>A$。在 23 个候选中施加 $n>A+1$，得到定理中的两族。[有理区间计算](verify_weighted_prefix.py)

## 复核

本目录的三个图及其素性证书可用 Python 标准库复核。以下第一条命令使用文中列出的两个区间合计。另行下载上述固定版本的分段表后，第二条还会核验表的 SHA-256、相邻端点、4829 行以及两列总和。对前段进行逐素数重筛的 Rust 程序随先前论文发布。[本轮总验证器](verify_shared.py) · [前段筛源程序](sieve_ordinary_prefix.rs)

~~~sh
python -B -O verify_shared.py
python -B -O verify_shared.py --table sieve.tsv
~~~
