---
title: "受限奇完全数的素因子集合与第 j 小普通素因子"
date: "2026-11-24"
author: "Codex"
abstract: "证明受限奇完全数在 2·10¹² 与 exp(10¹²) 之间必须有超过七百亿个不同的普通素因子。另构造一个固定素数集，包含所有受限奇完全数中异于 3 的素因子，计数函数为 O(x/(log x)²)。给出第 j 小普通素因子的显式递推上界，其中 j 可任意固定。"
formula: '\#\{r\in U:r\le x\}\ll x/(\log x)^2'
---

# 受限奇完全数的素因子集合与第 j 小普通素因子

若奇完全数恰有一个偶指数大于 2，称它为受限奇完全数。由 Ochem 与 Zelinsky，任何这样的数均可写成

$$
N=3^a p^e\prod_{s\in S}s^2,\qquad p\equiv1\pmod{12},\qquad e\equiv1\pmod4,
$$

其中各底数是不同的素数，$a$ 是偶数，且 $7\nmid N$。称 $p$ 为欧拉素因子，$S$ 的元素为普通素因子。以下所有对数均为自然对数。[Ochem–Zelinsky，定理 1、引理 6 及第 4 节](https://arxiv.org/html/2607.19746v1#S4)

**定理 1。** 存在一个固定的、可判定的素数集 $U$，使每个受限奇完全数的所有异于 3 的素因子都在 $U$ 内，且

$$
U(x)=\#\{r\in U:r\le x\}\ll \frac{x}{(\log x)^2},\qquad
\sum_{\substack{r\in U\\r>X}}\frac1{r-1}\ll\frac1{\log X}\quad (X\ge3).
$$

**定理 2。** 给定正整数 $j$，存在可显式递推的常数 $B_j$，使得每个至少有 $j$ 个普通素因子的受限奇完全数，其第 $j$ 小普通素因子至多为 $B_j$。此界与 $a,e$ 和 $\#S$ 无关。文末给出一组对所有 $j$ 同时适用的递推式。

同样，对至少有 $j$ 个模 3 余 1 的普通素因子的 $N$，这一子集中第 $j$ 小的素因子也有只依赖于 $j$ 的上界。

**定理 3。** 每个受限奇完全数在区间 $2\cdot10^{12}<s\le\exp(10^{12})$ 中有超过 $70\,000\,000\,000$ 个不同的普通素因子。这里使用显式的大筛常数与有限区间上的整数求和。

Yamada 在允许有限个例外指数时给出了最小素因子上界。Fletcher、Nielsen 与 Ochem 将大筛用于奇完全数的分圆因子及倒数和估计。这里最小素因子已经固定为 3。以下构造利用 $p+1$ 与 $s^2+s+1$ 的相互限制，从一个适用于整族的素数集推导任意固定次序的界。[Yamada](https://arxiv.org/abs/math/0511410) · [Fletcher–Nielsen–Ochem，命题 1](https://www.lirmm.fr/~ochem/opn/OPNS_Adam_Pace.pdf)

## 从两个模 84 剩余类开始

写 $\Phi(t)=t^2+t+1$。令 $F$ 为模 84 余 67 或 79 的素数集。每个 $q\in F$ 都模 12 余 7，所以只能是普通素因子。可是 $q$ 模 7 余 4 或 2，均给出 $7\mid\Phi(q)$。若 $q^2\parallel N$，则 $\Phi(q)=\sigma(q^2)\mid2N$，与 $7\nmid N$ 矛盾。因此 $F$ 中没有素数能整除任何受限奇完全数。

以下逐层定义四个集合，每次只需要判断整数的素因子。由此也给出了 $U$ 的成员判定方法。

1. $A$ 由所有满足 $s>3$ 且 $\Phi(s)$ 没有 $F$ 中素因子的素数 $s$ 组成。每个普通素因子都在 $A$ 中。
2. $E$ 由所有满足 $p\equiv1\pmod{12}$ 且 $(p+1)/2$ 的素因子全在 $A$ 中的素数 $p$ 组成。因 $e$ 是奇数，$p+1\mid\sigma(p^e)$。又 $(p+1)/2$ 为奇数、与 $3p$ 互素，故其中的素因子都是普通素因子。每个欧拉素因子都在 $E$ 中。
3. 令 $T=A\cup E$。取 $A'$ 为所有满足 $s>3$ 且 $\Phi(s)$ 的每个素因子都在 $\{3\}\cup T$ 中的素数 $s$。普通素因子的约数和必须整除 $2N$，故每个普通素因子也在 $A'$ 中。
4. 取 $U=A'\cup E$。

## 三次应用同一大筛估计

用如下上界形式。设某个素数集在筛素数 $q$ 处避开 $\rho(q)$ 个剩余类，其中包括 0，素数本身等于 $q$ 时允许例外。假定 $1\le\rho(q)\le3$、$\rho(q)<q$，且

$$
\sum_{q\le z}\frac{\rho(q)}{q}=\kappa\log\log z+O(1).
$$

则这个素数集在 $x$ 以下的元素个数为 $O(x/(\log x)^\kappa)$。可以忽略有限个筛素数。为说明所需的一致性，令 $Q=\lfloor\sqrt{x}\rfloor$，先计入至多 $Q$ 个小素数。对其余元素应用 Fletcher–Nielsen–Ochem 命题 1，所得上界为 $(x+Q^2)/G(Q)$，其中

$$
G(Q)=\sum_{d\le Q}\mu^2(d)\prod_{q\mid d}\frac{\rho(q)}{q-\rho(q)} .
$$

令 $u=Q^{1/8}$，展开 $\prod_{q\le u}(1-\rho(q)/q)^{-1}$。将各个平方自由整数项的权重归一为概率，则 $\log d$ 的期望等于 $\sum_{q\le u}\rho(q)\log q/q\le(3+o(1))\log u$。这里用了素数的 Mertens 估计。由 Markov 不等式，对充分大的 $Q$，至少半数权重来自 $d\le Q$。而整个乘积取对数后等于 $\sum_{q\le u}\rho(q)/q+O(1)$。故 $G(Q)\gg(\log x)^\kappa$。以下各次选取的余数类都由固定集合给出，所得常数对所有 $N$ 同时有效。[大筛不等式原式及乘积估计](https://www.lirmm.fr/~ochem/opn/OPNS_Adam_Pace.pdf)

对 $A$，每个筛素数 $q$ 禁用 0，遇到 $q\in F$ 再禁用 $\Phi$ 的两个根。因为 $q\equiv1\pmod3$，这两个非零根恰是三阶单位根。模 84 的两个互素剩余类给出 $\sum_{q\in F,\,q\le z}1/q=(1/12)\log\log z+O(1)$。这是固定模数的素数等差数列定理经分部求和的推论。故 $\kappa=1+2/12=7/6$，从而 $A(x)\ll x/(\log x)^{7/6}$。分部求和即得 $\sum_{q\in A}1/q<\infty$。[固定模数的素数估计](https://arxiv.org/abs/1802.00085)

对 $E$，每个奇素数 $q\notin A$ 除 0 外还禁用 $-1$。因 $\sum_{q\in A}1/q$ 收敛，筛维数为 2，故 $E(x)\ll x/(\log x)^2$。特别地，$\sum_{q\in T}1/q$ 也收敛。

最后，对 $A'$，在 $q\equiv1\pmod3$ 且 $q\notin T$ 时除 0 外禁用 $\Phi$ 的两个根。此处额外禁用类贡献 $2\bigl(\tfrac12\log\log z+O(1)\bigr)$，因为 $T$ 中素数的倒数和有限。于是 $A'(x)\ll x/(\log x)^2$。与 $E(x)$ 的界相加即得 $U(x)$ 的估计。再从 $X$ 到无穷分部求和，利用 $\int_X^\infty dt/(t\log^2t)=1/\log X$，即得定理 1。

## 固定次序的素因子不能趋于无穷

反设有一列受限奇完全数，其第 $j$ 小普通素因子趋于无穷。如果其中 $a$ 有界，由 Ochem–Zelinsky 引理 8 的 $\omega(N)\le2a+4$ 以及 Nielsen 的 $N<2^{4^{\omega(N)}}$，这样的 $N$ 只有有限多个。故可取子列使 $a\to\infty$。[Ochem–Zelinsky，引理 8](https://arxiv.org/html/2607.19746v1#S4) · [Nielsen 的上界](https://mathdept.byu.edu/~pace/upper_bound_web.pdf)

在前 $j-1$ 个位置上依次取子列，使每个位置的素数恒定或趋于无穷。令 $S_0$ 为恒定下来的普通素因子集合。剩余普通素因子中的最小者趋于无穷。写 $I(t)=\sigma(t)/t$。因为 $0<\log I(s^2)<1/(s-1)$，定理 1 保证所有剩余普通分量的 $I(s^2)$ 之积趋于 1，即使其个数无界也成立。

欧拉素数可再取子列使其趋于无穷或固定。前一情形 $1<I(p^e)<p/(p-1)$，所以极限为 1。后一情形中，$e$ 可固定或趋于无穷，极限依次为 $I(p^e)$ 或 $p/(p-1)$。记相应的极限为 $L$。由 $I(N)=2$ 取实数极限得有理数等式

$$
2=\frac32 L\prod_{s\in S_0}\frac{s^2+s+1}{s^2}.
$$

有限乘积的分子、分母均为奇数。三种情形下 $v_2(L)$ 分别为 0、1、至多 $-2$。其中固定 $p,e$ 时，$\sigma(p^e)$ 是 $e+1\equiv2\pmod4$ 个模 4 余 1 的数之和，故恰含一个因子 2。右侧的二进赋值因而分别为 $-1$、0、至多 $-3$，同左侧的 1 矛盾。这里是在得到有限个有理数的等式后取赋值。定理 2 得证。

## 模 3 余 1 的普通素因子

若 $\#S\ge2j$，取最小的 $2j$ 个普通素因子。按定理 2，它们都不超过 $B_{2j}$。在这组素数中，至多一个 $s$ 的 $\Phi(s)$ 缺少模 3 余 1 的普通素因子。这来自普通分量计数中三个例外值 $p,3p,3p^2$ 的互斥性。对其余每个 $s$，选取一个这样的素因子 $q\mid\Phi(s)$。同一 $q$ 最多被选两次，因为 $q$ 在 $2N=\sigma(N)$ 中的指数恰为 2。于是至少有 $\lceil(2j-1)/2\rceil=j$ 个不同的 $q$ 不超过 $B_{2j}^2+B_{2j}+1$。[三种例外与重数计数](counting.md)

若 $\#S<2j$，则 $\omega(N)\le2j+1$，Nielsen 的界给出 $N<2^{4^{2j+1}}$。取这个界与 $B_{2j}^2+B_{2j}+1$ 的较大者，便得到推论。

## 一个必须含有七百亿个普通素因子的固定区间

令 $m=5460=12\cdot5\cdot7\cdot13$。用 $R$ 表示模 $m$ 的互素剩余类中同时满足下列条件的类 $b$。

$$
b\equiv1\pmod3,\qquad \gcd(\Phi(b),91)>1,\qquad
\bigl(b\not\equiv1\pmod{12}\ \text{或}\ \gcd(b+1,455)>1\bigr).
$$

若素数 $q\bmod m\in R$，则 $7$ 或 $13$ 整除 $\Phi(q)$，排除了普通身份。最后一个条件排除了欧拉身份，因为欧拉素数必须模 12 余 1，且 $q+1$ 的任何奇素因子都必须整除 $N$，而 $5,7,13$ 均已排除。枚举给出 $\varphi(m)=1152$、$\#R=172$。记 $\mathcal F$ 为这些类中的素数，$\delta=172/1152$，$\alpha=2\delta=43/144$。

设 $L=10^{10}$。对 $t\ge L$，将 Bennett–Martin–O'Bryant–Rechnitzer 定理 1.2 应用于上述 172 个剩余类，得到

$$
\theta_{\mathcal F}(t)\ge\delta t-\frac{172t}{840\log t}>0.
$$

对 $1/(t\log t)$ 分部求和，保留从 $L$ 到 $u$ 的正积分并用 $-(1/(t\log t))'\ge1/(t^2\log t)$，即得

$$
\sum_{\substack{q\in\mathcal F\\q\le u}}\frac1q
\ge\delta\log\frac{\log u}{\log L}-\frac{172}{840\log L}\qquad (u\ge L).
$$

定义 $A_*$ 为使 $\Phi(s)$ 不含 $\mathcal F$ 中的素因子，也不含 $q\le10^4$、$q\equiv1\pmod3$ 的素因子的素数 $s>3$。每个普通素因子都属于 $A_*$。用 $\rho(q)=3$ 筛上述两类 $q$，其余素数取 $\rho(q)=1$。写 $V(u)=\prod_{q\le u}(1-\rho(q)/q)^{-1}$。先提取所有素数的乘积 $\prod_{q\le u}(1-1/q)^{-1}\ge\sum_{n\le u}1/n\ge\log u$。每个另筛两个根的 $q$ 还贡献 $(q-1)/(q-3)$，其对数大于 $2/q$。不属于 $\mathcal F$ 的 428 个额外小素数满足

$$
\prod_{\substack{q\le10^4,\ q\equiv1\pmod3\\q\notin\mathcal F,\ q\ {\rm prime}}}\frac{q-1}{q-3}>4 .
$$

这个比较可直接交叉相乘。由 $20<\log L<24$ 可得 $(\log L)^\alpha\exp(344/(840\log L))<3$。代回上式便有 $V(u)>\tfrac43(\log u)^{1+\alpha}$。[显式素数定理](https://arxiv.org/abs/1802.00085) · [整数常数复核](verify_support_bounds.py)

为给出带常数的计数界，取 $\log x\ge1000$、$Q=\lfloor\sqrt x\rfloor$、$u=Q^{1/12}$。有 $\log u\ge(\log x)/25$，故 $u\ge L$。令 $k=\lfloor u\rfloor$。因 $\lfloor k/q\rfloor\ge k/(2q)$，且 $\sum_{q\le k}\lfloor k/q\rfloor\log q\le\log(k!)\le k\log k$，可得 $\sum_{q\le u}(\log q)/q\le2\log u$。故大筛乘积中 $\log d$ 的期望至多为 $6\log u=\tfrac12\log Q$，从而 $G(Q)\ge V(u)/2$。用 $25^{1+\alpha}<66$ 与 $\sqrt x\le x/(\log x)^{1+\alpha}$ 得

$$
A_*(x)\le Q+\frac{x+Q^2}{G(Q)}
\le199\,\frac{x}{(\log x)^{1+\alpha}} .
$$

现在取 $Y=\exp(10^{12})$。分部求和并用 $1/(s-1)\le1.01/s$，得到

$$
\sum_{\substack{s\in A_*\\s>Y}}\frac1{s-1}
\le\frac{1.01\cdot199}{\alpha(\log Y)^\alpha}
<\frac{675}{3800}.
$$

最后一步用到 $1.01\cdot199/\alpha<675$ 和 $(\log Y)^\alpha=10^{43/12}>3800$。在 $X=2\cdot10^{12}$ 以内，分段筛给出所有可能普通素因子的倒数和 $\sum1/(s-1)\le72\,401\,332\,985\,083\,531/10^{18}$。若 $(X,Y]$ 内至多有 $70\,000\,000\,000$ 个普通素因子，它们的和至多为 $70\,000\,000\,000/X$。加上欧拉分量及上面的尾界，便有

$$
\frac1{2\,999\,999}+\sum_{s\in S}\frac1{s-1}
\le \frac1{2\,999\,999}
+\frac{72\,401\,332\,985\,083\,531}{10^{18}}
+\frac{70\,000\,000\,000}{2\cdot10^{12}}
+\frac{675}{3800}
<\frac{31}{108}<\log(4/3).
$$

最右侧的有理数是 $\log(1+1/3)$ 交错级数的前四项。另一方面，$I(3^a)<3/2$、$\log I(p^e)<1/(p-1)$、$\log I(s^2)<1/(s-1)$，且 $p>3\,000\,000$，所以完全性要求左侧大于 $\log(4/3)$。这证明了定理 3。[有限区间的计算与复算命令](../restrained-exponent/smooth-indices.md) · [筛的整数结果](../restrained-exponent/ordinary-prefix-2000000000000.json)

## 对任意次序计算上界

取 $B_0=1$。定理 2 中可对每个 $j\ge1$ 递推取

$$
D_j=3(B_{j-1}^2+B_{j-1}+1)^{j-1},\qquad
B_j=\exp\!\bigl((21600D_j^2)^4\bigr).
$$

定理 3 对前七百亿个位置给出的 $\exp(10^{12})$ 更小。一般递推使用上一节的尾界，即对 $\log X\ge1000$ 有 $\sum_{s\in A_*,\,s>X}1/(s-1)\le675/(\log X)^{43/144}$。

先证一个有理数间距估计。设 $T$ 为最简分母为 $b$ 的正有理数，$v_2(T)=2$。对任意欧拉素数及指数 $p\equiv e\equiv1\pmod4$，有

$$
\big|T-I(p^e)\big|\ge\frac1{4b^2}.
$$

若 $p>2b+1$，用 $|T-1|\ge1/b$ 与 $I(p^e)-1<1/(p-1)$ 即得。若 $p\le2b+1$ 且 $p^e\ge2b$，则 $p/(p-1)$ 的二进赋值为负，故 $|T-p/(p-1)|\ge1/(b(p-1))$。又 $p/(p-1)-I(p^e)=1/((p-1)p^e)$，至多为前一间距下界的一半。因此 $|T-I(p^e)|\ge1/(2b(p-1))\ge1/(4b^2)$。最后，若 $p^e<2b$，则 $v_2(I(p^e))=1$，交叉相乘即得 $|T-I(p^e)|\ge1/(bp^e)>1/(2b^2)$。

现在归纳证明递推界。将普通素因子排为 $s_1<s_2<\cdots$。取 $M=\prod_{i<j}I(s_i^2)$、$T=4/(3M)$，并令 $D=D_j$。每个 $I(s_i^2)$ 是奇数之比，因此 $v_2(T)=2$。由归纳假设 $s_i\le B_{j-1}$，$T$ 的最简分母至多为 $3\prod_{i<j}(s_i^2+s_i+1)\le D$。在 $j=1$ 时取空乘积为 1。

若 $3^{a+1}<32D^2$，由 $\omega(N)\le2a+4$ 和 Nielsen 的界有

$$
\log N<4^{2a+4}=16(2^{a+1})^4<16(32D^2)^4<\log B_j.
$$

若 $3^{a+1}\ge32D^2$，令 $z=3^{-(a+1)}$。反设 $s_j>B_j$。所有剩余普通分量的对数丰度之和 $H=\sum_{i\ge j}\log I(s_i^2)$ 由上述尾界满足 $H\le1/(32D^2)$，因为 $43/144>1/4$。又 $z\le1/(32D^2)$。完全性方程给出 $I(p^e)=T\exp(-H)/(1-z)$。用 $T\le4/3$、$z\le1/3$ 及 $1-\exp(-H)\le H$，得到

$$
|T-I(p^e)|\le\frac{T}{1-z}(z+H)\le\frac1{8D^2}.
$$

前面的有理数间距给出同一差的下界 $1/(4D^2)$。矛盾证明 $s_j\le B_j$，归纳完成。

另一个利用 $v_q(2N)=2$ 的筛把七百亿个普通素因子的区间上端降到 $10^{75}$，并给出更大区间中的计数。[显式三次分圆筛](../exponent-sieve/ternary-tail-constants.md)
