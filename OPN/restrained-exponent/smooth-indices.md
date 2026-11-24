---
title: "受限奇完全数的指数加一至少有一个不小于 2383 的素因子"
date: "2026-11-24"
author: "Codex"
abstract: "设 a 是受限奇完全数中 3 的指数。证明 a+1 的最大素因子至少为 2383。若 857 整除 a+1，则最大素因子至少为 2851。如果 a+1 的所有素因子都不超过 2850，给出包含其所有可能取值的 173 个整数，以及这一子类的有效有限上界。"
formula: 'P^+(v_3(N)+1)\ge2383'
---

# 受限奇完全数的指数加一至少有一个不小于 2383 的素因子

设 $N$ 是恰有一个偶指数大于 2 的奇完全数，记 $n=v_3(N)+1$，$P^+(n)$ 为 $n$ 的最大素因子。Ochem 与 Zelinsky 证明

$$
N=3^a p^e\prod_{s\in S}s^2,\qquad n=a+1,\quad p\equiv1\pmod{12},\quad e\equiv1\pmod4 .
$$

其中 $a$ 为偶数，除 3 外的素因子均大于 $3\,000\,000$。称 $p$ 为欧拉素因子，$S$ 的元素为普通素因子。本文沿用这些已知结构条件。[Ochem–Zelinsky，定理 1 与第 4 节](https://arxiv.org/html/2607.19746v1#S4)

**定理。** 对每个受限奇完全数，$P^+(n)\ge2383$。若 $857\mid n$，则 $P^+(n)\ge2851$。若 $P^+(n)\le2850$，则 $n$ 必在所附的 173 个整数中，其中最小者为 $2\,106\,997\,056\,806\,890\,291$，最大者为 $M=2383^{27}$。这些整数是完全性方程的指数候选。[173 个整数的升序表](smooth-candidates.json)

这也给出整类数的有限界。已有的分量计数给出 $\omega(N)\le2a+3=2n+1$。由 Nielsen 的上界，若 $P^+(n)\le2850$，必有

$$
N<2^{\,4^{2M+1}},\qquad M=2383^{27}.
$$

即使 $a$ 的指数候选很少，允许的 $N$ 仍可能极大。[普通素因子计数](../ordinary-components/counting.md) · [Nielsen，An upper bound for odd perfect numbers](https://mathdept.byu.edu/~pace/upper_bound_web.pdf)

## 从乘法阶排除整除条件

设 $d>1$，$r\ge5$ 是满足 $\operatorname{ord}_r(3)=d$ 的素数。若 $d\mid n$，则 $r\mid(3^n-1)/2=\sigma(3^a)$，从而 $r\mid N$。普通素因子 $r$ 强迫 $r^2+r+1$ 的所有素因子整除 $N$。若 $r$ 是欧拉素因子，则还须有 $r\equiv1\pmod{12}$，且 $r+1$ 的每个奇素因子都整除 $N$。只要两种可能身份都导向 $5,7,11,13$ 中已排除的素数，便得 $d\nmid n$。找不到见证的 $d$ 可以留在候选中。

例如 $d=971\cdot991=962\,261$。素数 $r=1\,924\,523$ 满足 $3^d\equiv1\pmod r$，而 $3^{971}$ 与 $3^{991}$ 模 $r$ 分别余 $832\,938$ 与 $264\,749$。$r\equiv11\pmod{12}$ 只能作为普通素因子，但 $r\equiv3\pmod{13}$ 会强迫 $13\mid N$。故 971 与 991 不能同时整除 $n$。

影响最大素因子下界的一条短链是

$$
d=1571\cdot1789^3,\qquad
r=3\,213\,328\,532\,470\,934\,771,\qquad
r\longrightarrow487\longrightarrow7 .
$$

这里 $r$ 是素数，且 $3^d\equiv1\pmod r$。将 $d$ 分别除以 1571、1789 后，模幂的余数分别为 $1\,347\,632\,346\,044\,148\,145$ 和 $1\,092\,791\,674\,411\,150\,566$。两个箭头都来自普通分量，因为 $r\equiv11\pmod{12}$，$487\equiv7\pmod{12}$。整除性可由 $r\equiv254\pmod{487}$、$487\equiv4\pmod7$ 直接复算。

另有两个较大的见证 $r=(3^\ell-1)/2$，对应 $\ell=1091,1367$，位数分别为 521 与 652。它们的普通分支都强迫 7，其欧拉分支也由证书闭合。$r-1$ 的完整因数分解可用来递归证明 $r$ 的素性。寻找见证时使用过 Cunningham 表、Brent 因数表、FactorDB 以及形如 $2kd+1$ 的试除。复核只读取附上的整数。[旧指数证书](index-factors.md) · [Brent 因数表](https://maths-people.anu.edu.au/~brent/factors.html) · [FactorDB](https://factordb.com/)

## 完整枚举光滑指数

因为 $n$ 为奇数，可先排除 $2\mid n$。在原来排除的 152 个奇素数和 $857^2$ 的基础上，再加入 737 个整除障碍。对每个见证，验证器证明素性，逐个检查 $3^d\equiv1\pmod r$ 和每个 $\ell\mid d$ 的 $3^{d/\ell}\not\equiv1\pmod r$，其中 $\ell$ 遍历素因子。它还检查对每种允许身份均有整除边，以及图中每条路径都能闭合。新增图有 1,111 个内部素数，1,608 条身份分支，另附 274 个椭圆曲线证明步骤和 1,390 个 Lucas 素性证明。[新增因子链](smooth-chains.json) · [新增素性证书](smooth-primality.json) · [原因子链](index-chains.json) · [原素性证书](primality.json)

小于等于 2850 的 413 个素数中，有 24 个留下来。依次为

$$
\begin{gathered}
857,\ 971,\ 991,\ 1109,\ 1163,\ 1571,\ 1597,\ 1627,\ 1747,\ 1787,\ 1789,\ 1949,\\
2239,\ 2333,\ 2383,\ 2441,\ 2521,\ 2557,\ 2579,\ 2591,\ 2677,\ 2683,\ 2707,\ 2797 .
\end{gathered}
$$

对每个留下的素数 $\ell$，证书至少排除了一个 $\ell^b$，故只需枚举 $v_\ell(n)$ 从 0 到 $b-1$。按素因子从小到大递归，每添加一次 $\ell$，检查最大素因子为 $\ell$ 的所有混合障碍。首次遇到整除关系后，更高次幂也都被排除。这枚举了全部满足所列条件的光滑整数。

| 只使用整除障碍时的集合 | 核验结果 |
| --- | ---: |
| 所有整数，包括 1 | 1,305 个 |
| 最大整数 | $2383^{27}$ |
| 857 的倍数 | 42 个 |
| 857 的最大倍数 | $857\cdot2383^2\cdot2591=12\,609\,455\,101\,543$ |
| 最大素因子小于 2383 时的最大整数 | $971^5\cdot1109=957\,255\,115\,116\,280\,759$ |

下一节给出 $a>1.5\cdot10^{18}$ 的复算。删去较小整数后，正好留下 173 个，并且都含有不小于 2383 的素因子，没有一个是 857 的倍数。$P^+(n)>2850$ 的情形本来就满足两个最大素因子界。由此得到定理。尤其是处理 857 的倍数时，无须分解 409 位合数 $(3^{857}-1)/2$。

## 对指数下界作独立数值核验

Stormwhip 给出了比本节所需更强的指数下界。这里采用该作者的 97 项筛多项式，自行重算有限区间，并用显式等差数列素数定理估计后段。设 $I(t)=\sigma(t)/t$。因 $I(3^a)<3/2$、$\log I(p^e)<1/(p-1)$、$\log I(s^2)<1/(s-1)$，完全性要求

$$
\log(4/3)<\frac1{2\,999\,999}+\sum_{s\in S}\frac1{s-1}.
$$

若 $a\le A$，分量计数给出 $\#S\le2A+1$。以下取 $A=1.5\cdot10^{18}$，$X=2\cdot10^{12}$，$Z=250A$。[Stormwhip 的原验证程序](https://github.com/stormwhip/OddPerfectNumberBounds/blob/main/verify.py) · [分量计数](../ordinary-components/counting.md)

在 $3\,000\,000<s\le X$ 上筛出素数，再删去 $\Phi(s)=s^2+s+1$ 含有 $[5,3\,000\,000]$ 中任一素因子的 $s$。若 $q\ne3$ 且 $q\mid\Phi(s)$，则 $q\equiv1\pmod3$，只需筛两个三阶单位根。对所有留下的 $s$ 累加整数 $\lceil10^{18}/(s-1)\rceil$。至 $2\cdot10^{12}$ 共留下 $8\,391\,832\,364$ 个素数，整数和为 $72\,401\,332\,985\,083\,531$。筛至 $2\cdot10^{11}$ 时得到 $917\,142\,659$ 个素数和整数和 $62\,696\,136\,715\,648\,189$，与 Stormwhip 的前六个已发表区间汇总相符。[分段筛源程序](sieve_ordinary_prefix.rs) · [较小区间的本地重算](ordinary-prefix-200000000000.json) · [完整区间的本地重算](ordinary-prefix-2000000000000.json) · [原分段记录](https://github.com/stormwhip/OddPerfectNumberBounds/blob/main/sieve.tsv)

在 $X<s\le Z$ 上，取

$$
\mathcal P=\{7,13,19,31,37,43,61,67,73,79,97,103\},\qquad J_q(s)=\boldsymbol{1}_{q\mid\Phi(s)}.
$$

取多项式 $G=1-\sum J_q+\sum_{q<q'}J_qJ_{q'}-\sum J_7J_iJ_j$，末和遍历 $7<i<j$、$7ij\le10000$ 的 $i,j\in\mathcal P$。逐一检查 $2^{12}$ 种布尔值，得 $G\ge0$ 且 $G(0,\ldots,0)=1$。每个普通素因子均使 $G=1$。

令 $\Theta_G(t)=\sum_{r\le t}G(r)\log r-\boldsymbol{1}_{t\ge3}\log3$，其中 $r$ 遍历素数。这样，常数项恰可用模 3 的两个互素剩余类计算，且 $(X,Z]$ 上的每个素数权重仍为 $G(r)\log r$。含有 $k\ge1$ 个指标的单项式对应模 $\prod q$ 的 $2^k$ 个互素剩余类。Bennett、Martin、O'Bryant 与 Rechnitzer 的定理 1.2 对 $3\le q\le10000$、$t\ge8\cdot10^9$ 给出每个互素类的 $\theta$ 误差小于 $t/(840\log t)$。主项及误差系数分别为

$$
\delta=\frac{64\,959\,667}{169\,646\,400},\qquad C=\frac{434}{840}=\frac{31}{60}.
$$

故上述按剩余类求和的函数满足 $|\Theta_G(t)-\delta t|\le C t/\log t$。误差系数中的 434 已计入常数项的两个类。[显式素数估计，定理 1.2](https://arxiv.org/abs/1802.00085)

把该估计对权 $w(t)=1/(t\log t)$ 分部求和，并使用 $1/(s-1)\le X/((X-1)s)$。两个端点的误差至多为 $C/\log^2 X+C/\log^2 Z$，积分误差至多为 $C\int_X^Z(1/(t\log^2t)+1/(t\log^3t))\,dt$。于是

$$
\sum_{\substack{s\in S\\X<s\le Z}}\frac1{s-1}
\le\frac X{X-1}\left[\delta\log\frac{\log Z}{\log X}
+C\left(\frac1{\log X}-\frac1{\log Z}+
\frac{3}{2\log^2 X}+\frac1{2\log^2 Z}\right)\right].
$$

最后 $s>Z$ 的贡献至多为 $(2A+1)/(Z-1)$。有限区间、两个尾段和欧拉素因子的贡献相加，上界不超过 $0.285767329932847170$，而 $\log(4/3)>0.287682$，与完全性矛盾。验证器用有理区间计算每个对数，所列比较均已向外舍入。[有理区间验证程序](verify_ordinary_bound.py)

对同一终点 $2\cdot10^{12}$，汇总 Stormwhip 表内连续的区间，得到整数和 $73\,802\,690\,054\,260\,520$。其中四段只筛到 $100\,000$，所以该汇总也是普通素因子倒数和的上界。换用这份独立数据，上段的最终上界为 $0.287168687002024159$，仍小于 $\log(4/3)$。在固定区间的七百亿个普通素因子结论中，换用这份汇总也仍满足所需的严格不等式。[原始分段表](https://github.com/stormwhip/OddPerfectNumberBounds/blob/main/sieve.tsv) · [固定区间的证明](../ordinary-components/prime-support.md)

若另外使用 Stormwhip 报告的 $a>1.14\cdot10^{22}$，候选表还剩 109 项。这项附加计数用到该作者规模更大的计算。上面的 173 项及两个最大素因子界只需本节复算的区间。[作者对较强下界的说明](https://github.com/stormwhip/OddPerfectNumberBounds)

## 离线复核

将本目录的 Python、Rust 源程序和 JSON 放在一起。在 64 位系统上使用支持整数平方根的 Rust，以优化模式编译分段筛，输入为排除上限、区间终点和工作线程数。将程序输出与附件中两个同名区间的 JSON 对照。指数验证只需要 Python 标准库，运行时重新验证旧图与新图并比较完整候选表。

~~~sh
rustc -O sieve_ordinary_prefix.rs -o sieve_ordinary_prefix
./sieve_ordinary_prefix 3000000 200000000000 4
./sieve_ordinary_prefix 3000000 2000000000000 4
python -B -O verify_ordinary_bound.py
python -B -O verify_smooth.py
~~~

[总验证程序](verify_smooth.py) · [光滑数枚举](smooth_products.py) · [原图与素性验证程序](verify_exponent.py)

后续计算把这 173 个指数候选压缩成两族，并给出了每一族的幂次范围。[两族指数及新增证书](shared-indices.md)
