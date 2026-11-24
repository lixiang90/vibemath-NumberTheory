---
title: "用辅助二次 L 函数的条件协方差分离孪生素数"
date: "2026-11-25"
abstract: "在 Alpöge–Claude 的第二阶段候选素数族中，每个导子只取一个二次 L 函数。按物理端点乘积的二次符号分组后，两个不相交收缩零点窗的协方差之差检出固定差二的真孪生行。三阶混合累积量则测出半素数邻项中两个因子的相对大小。利用陈氏定理，在原来的四符号筛选和几何好窗口上，二阶与校准后的三阶质量之和有固定正下界。每个二倍区间只用三个固定窗心即可给出半素数质量的上界，相对余量至多 0.4%。"
---

# 用辅助二次 L 函数的条件协方差分离孪生素数

最近三篇用不同办法处理同一难点。第 30 篇以多区间波动标定陈氏候选的噪声。第 31 篇算出了差二的主弧主项，剩下一个补弧系数。第 32 篇保留同一行两个端点的符号，使交叉项可以单独提取。对于固定差二，我选择第三条思路，并把符号换成原辅助素数在整数 $n(n+2)$ 上的真实二次特征值。它由 $n,d$ 直接计算，能在知道端点素性以前指定。零点数据只需每个候选导子 $d$ 的一个 $L(s,\chi_d)$。[多区间波动](https://lixiang90.github.io/vibemath-NumberTheory/zeros-exact-distribution/articles/gap-two-window-quadratic-variation.html) · [主弧与补弧](https://lixiang90.github.io/vibemath-NumberTheory/zeros-exact-distribution/articles/gap-two-auxiliary-conductor-bands.html) · [前篇的平移符号](https://lixiang90.github.io/vibemath-NumberTheory/zeros-exact-distribution/articles/gap-two-translated-quadratic-characters.html)

计算表明，若 $n,n+2$ 都是素数，按这个符号分出的两族首先在**二阶**零点统计上出现差异。若一端是素数、另一端是两个不同素数的乘积，首个差异出现在**三阶**。陈氏定理给这两类的总和一个真正的正下界，原来的四道符号条件和几何选窗可以同时保留。要由此推出孪生素数的正质量，还缺三阶部分的严格上界。以下使用归档的 Dirichlet-GRH 为输入，其余结论为本文的纸面推导。[归档定理](https://github.com/leanprover-community/riemann-hypothesis#theorem-statements) · [Green–Tao 引用的陈氏定理，定理 6.1](https://arxiv.org/pdf/math/0405581#page=19)

## 1. 两种不同的辅助素数平均

沿用 ζ 构造。令 $P\to\infty$，$H=P^J$，$0<\epsilon<1/52$、$J$ 固定，$z=P^{\epsilon/5}$，$z'=z(\log P)^2$，$X=e^{z'}$。先取 $e^{z/2}<Q\le e^z$、$Q\equiv3\bmod8$ 的素数。给定 $Q$，令 $\mathcal D_Q$ 为 $(\sqrt X,X]$ 内满足下列条件的全部素数，

$$
d\equiv7\pmod {16},\qquad d\equiv-1\pmod {Q\prod_{3\le p\le z}p}.
\tag{1}
$$

乘积只取奇素数。以 $\mathbb E_d$ 记 $\mathcal D_Q$ 上的均匀平均，$\chi_d(m)=(m/d)$。原证明从这些候选取见证 $Q,d_0$，继而选出完整好窗口 $[Y,Y^2]\subset[P,H]$。计算窗内零点统计时，**先固定这个见证及 $Y$，再以整个 $\mathcal D_Q$ 为观测族**。四道符号条件始终用实际见证 $Q,d_0$。后面证明很多见证同时有效时，另对原来的两阶段抽样取概率。[原候选的同余](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/MeanSquareBounds/PrimeSetTwoBounds.lean#L37) · [取对及选窗顺序](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/CharacterBoxes/WinPairs.lean#L32)

写 $M=16Q\prod_{3\le p\le z}p$。对所有 $p\le z$ 都有 $\chi_d(p)=1$。对于奇 $p$，互反律的符号与 $(-1/p)$ 抵消。对 $2$ 用 $d\equiv7\bmod8$。设 $R_0>1$ 为与 $M$ 互素的奇平方自由整数。在 (1) 的一个可逆类中用特征正交及二次互反律，得到

$$
\big|\mathbb E_d\chi_d(R_0)\big|
 \ll \frac{\varphi(M)\log X\,\log^2(XMR_0)}{\sqrt X}.
\tag{2}
$$

展开可逆类后，每项为模 $M$ 的特征乘 Jacobi 特征 $(\cdot/R_0)$。后者在与 $M$ 互素的平方自由模上非主，故没有主项。GRH 素数特征和界及分部求和给分子 $O(\sqrt X\log^2(XMR_0))$。用同一估计分离主特征，得族的大小渐近于 $\operatorname{Li}(X)/\varphi(M)$，这里 $\log M=O(z)=o(z')$。对第一阶段的 $Q$ 也有相同论证，把 $(X,M)$ 换成 $(e^z,8)$ 即可。原归档的扭转引理控制固定范围并容许例外。式 (2) 在最终 GRH 之上使用，尤其允许下节的素数标签远超 $H$。[原扭转估计及其量词](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/MeanSquareBounds/twistK.lean#L145) · [Lee 的 GRH 特征素数和估计](https://arxiv.org/abs/2301.13457)

## 2. 两个可从真实零点计算的线性统计

取 $B=z/8$，$K=\lceil(\log P)^2\rceil$，$b=B/(4K)$，$R=2/b$。对数变量上的核及其 Fourier 变换（相位 $e^{-itu}$）为

$$
g=\mathbf1_{[-B,B]}*
  \left(\frac{\mathbf1_{[-b,b]}}{2b}\right)^{*2K},
\qquad
\widehat g(t)=\frac{2\sin(Bt)}{t}
\left(\frac{\sin(bt)}{bt}\right)^{2K}.
\tag{3}
$$

连续延拓定义零处的值。$0\le g\le1$，在 $|u|\le B/2$ 上等于 $1$，支撑于 $|u|\le3B/2$。以 $1/2+i\gamma_d$ 表示 $L(s,\chi_d)$ 的非平凡零点，按重数计。对 $|h|\le1$ 定义

$$
S_d(h)=\sum_{|\gamma_d-h|\le R}\widehat g(\gamma_d-h),\qquad
T_d(h)=-\frac12\left(S_d(h)-\log d-\mathbb E_d[S_d(h)-\log d]\right).
\tag{4}
$$

这里只要导子和一个半宽 $R=O((\log P)^2/z)$ 的零点窗。定义适用于实数 $x>1$ 的系数

$$
a_h(x)=\frac{\log x}{\sqrt x}
\sum_{j\ge0}x^{-j}g((2j+1)\log x)\cos((2j+1)h\log x).
\tag{5}
$$

**线性化。** 对任意固定 $A>0$，一致有

$$
T_d(h)=\sum_{z<p\le e^{3B/2}}a_h(p)\chi_d(p)+O_A(P^{-A}).
\tag{6}
$$

对平移的 $\widehat g$ 用 Weil 显式公式，导子项是 $g(0)\log d=\log d$，无穷处的其余项在全族中相同，因为特征均为奇特征。对实二次特征，两个素数幂项相加给 $-2\cos(h\log p^m)\chi_d(p)^m g(\log p^m)\log p/p^{m/2}$。算术支撑的上界 $e^{3B/2}<Q<\sqrt X$，所以所有偶次幂的贡献相同。小素数的奇次幂也相同。其余奇次幂按素数归并即 (5)。由 (2) 及 $\sum_p|a_h(p)|\le \exp(3B/4+O(\log B))$，在全族中心化时减去的剩余均值小于任意 $P^{-A}$。[Hughes–Rudnick，第 2 节的显式公式](https://arxiv.org/pdf/math/0208230#page=4)

为从全零点和截到 (4)，在 $2^jR<|\gamma_d-h|\le2^{j+1}R$ 上用 $|\widehat g|\le 2B(2^{j+1})^{-2K}$。该段的零点数为 $O((1+2^jR)\log(d(3+2^jR)))$，由单位长度零点计数界得到。求和给截断误差 $O(B\log(3X)\,4^{-K})=O_A(P^{-A})$，包括紧邻实轴的第一段。[Bennett、Martin、O’Bryant、Rechnitzer 的 Dirichlet 零点计数](https://arxiv.org/pdf/2005.02989#page=2)

## 3. 首次能区分两组零点的矩的阶数

设 $N$ 平方自由，素因子均在 $(z,H+2]$，$N\le(H+2)^2$，$r=\omega(N)\ge1$。用下标 $+$、$-$ 表示分别条件于 $\chi_d(N)=+1,-1$。两群的比例均为 $1/2+O_A(P^{-A})$。写 $\kappa^\pm(h_1,\ldots,h_m)$ 为 $T_d(h_1),\ldots,T_d(h_m)$ 在相应条件族的联合累积量。二阶就是条件协方差，三阶是三个条件中心化变量乘积的期望。置

$$
D(h_1,\ldots,h_m;N)=
\tfrac12\bigl(\kappa^+(h_1,\ldots,h_m)-\kappa^-(h_1,\ldots,h_m)\bigr).
$$

**定理 1。** 对固定 $m<r$，有 $D=O_A(P^{-A})$。在 $m=r$ 处，设 $p_1,\ldots,p_r$ 为 $N$ 的素因子，则

$$
D(h_1,\ldots,h_r;N)=
\sum_{\pi\in S_r}\prod_{j=1}^r a_{h_j}(p_{\pi(j)})+O_A(P^{-A}).
\tag{7}
$$

估计对这些 $N$ 和 $|h_j|\le1$ 一致。特别地，全取 $h_j=0$ 时，首个差异的主项严格正。在本范围内 $r\le 2\log(H+2)/\log z=O_{\epsilon,J}(1)$。

**证明。** 先把 (6) 中自由素数处的符号换成独立均匀正负号 $\varepsilon_p$。一个固定次数的单项式乘上可选的 $\chi_d(N)$ 后，奇重素数的乘积若为 $1$，它在两模型中均恒为 $1$。否则 (2) 给期望误差 $\exp(-(1/2-o(1))z')$，因为涉及的素数均小于 $\sqrt X$，且该乘积的对数为 $O_m(B+\log H)$。展开 $m$ 个线性和，绝对系数和至多 $\exp(O_m(B))$，误差仍小于任意负幂。条件期望可用 $1\pm\chi_d(N)$ 及组的概率写出，同样得到比较。独立和的每个固定 $j$ 阶绝对矩为 $O_j((1+B)^j)$，这由 $\sum_p a_h(p)^2\ll B^2$ 和独立符号的指数矩界推出。在实际族中比较偶次矩可得相同控制，条件化至多花一个常数。因此 (6) 的逐导子截断误差乘其他因子时，只需付出多项式大小的矩，仍可取到任意负幂精度。

在独立模型中，两个条件矩的半差等于原联合单项式乘 $\prod_{p\mid N}\varepsilon_p$ 的期望。要有贡献，$N$ 的每个素因子都须出现奇数次。低于 $r$ 阶时做不到，恰在 $r$ 阶时每个只能出现一次，遂得 (7) 的排列和。转成联合累积量时，其余分拆都只含更低阶的矩，它们在两组相同。$a_0(p)\ge\log p/\sqrt p$，因而全零窗心的首项至少为 $r!\prod_{p\mid N}\log p/\sqrt N$，误差能分辨单个这样的 $N$。

## 4. 固定差二的协方差与原窗口质量

保留 $n,n+2>1$ 且 $N=n(n+2)$ 平方自由、最小素因子大于 $z$ 的奇数行。它们两端互素，所以 $r=2$ 恰好表示两端都是素数。大于 $P$ 的所有孪生行都保留。对于 $U<n\le2U$、$P\le U$、$2U\le H$，取 $h=\pi/\log U$。设

$$
w(n)=\frac{\log n\log(n+2)}{n},\qquad
K_h(n)=-a_0(n)a_h(n+2)-a_h(n)a_0(n+2).
$$

由 (5)，对 $z<x\le H+2$ 有 $a_h(x)=(\log x/\sqrt x)(\cos(h\log x)+O(1/x))$。这对合数自变量也成立。两个余弦均为 $-1+O(1/\log^2U)$，故 $K_h(n)=(2+o(1))w(n)>0$，并且 $h>2R$。定理 1 给出

$$
-D(0,h;n(n+2))=
K_h(n)\,\mathbf1_{\mathbb P}(n)\mathbf1_{\mathbb P}(n+2)+O_A(P^{-A}).
\tag{8}
$$

左边用的是同一个 $L$ 函数两段互不相交的零点。展开零点对时没有同一零点的对角贡献。只对端点**乘积**的符号分组很重要。分别固定两个素数端点的符号会同时固定相应的两个线性项，条件中心化就会把它们消去。这里无需先判断 $n$ 或 $n+2$ 的素性。

在原好窗口内记 $I_Y=(2Y,Y^{3/2}/2]$，$L_Y=\sum_{Y\le j\le Y^{3/2}/2}1/j$（对整数求和），$h_Y=4\pi/(5\log Y)$。这次两端的余弦均至多为 $-\cos(\pi/5)+o(1)$，从而整个幂次窗口有 $K_{h_Y}(n)\asymp w(n)$，仍有 $h_Y>2R$。令 $I_{Q,d_0}(n)$ 指示两端在 $Q$ 处均为剩余，在 $d_0$ 处均为非剩余。以星号表示本节保留的行，定义

$$
\mathcal V(Y)= -\frac1{L_Y}\sum_{n\in I_Y}^{*}
 I_{Q,d_0}(n)\frac{w(n)}{K_{h_Y}(n)}D(0,h_Y;n(n+2)).
\tag{9}
$$

于是 $\mathcal V(Y)=J_{Q,d_0}(Y)+O_A(P^{-A})$，这里 $J_{Q,d_0}(Y)$ 是通过四道符号条件的真实孪生行的 $w(n)/L_Y$ 之和。相加时至多损失 $O(H)$，可预先增大误差指数。由于 $r\ge2$，两组的一阶矩之差为 $O_A(P^{-A})$，全族的均值为零。故 (9) 也可把 $D$ 换成 $\mathbb E_d[\chi_d(n(n+2))T_d(0)T_d(h_Y)]$，具有同一精度。第 17、18 篇的 $\Gamma_1$ 此时等于 $\tfrac13\limsup\mathcal V(Y)$，取同样的原好窗口见证。[第 18 篇的 $J$、$\Gamma_1$ 和四个符号](https://github.com/lixiang90/vibemath-NumberTheory/blob/main/Goldbach_TwinPrimes/twin-prime-auxiliary-signs.md)

## 5. 三阶零点统计能测出半素数的因子大小

为在原好窗口上使用已知的陈氏下界，把保留行稍加限制。记 $\mathcal R_Y$ 为 $I_Y$ 中 $N=n(n+2)$ 平方自由且每个素因子都大于 $n^{3/11}$ 的行。真孪生行仍全部保留。记 $A(Y)$ 为其中 $\omega(N)\in\{2,3\}$ 的行通过 $I_{Q,d_0}$ 后的 $w(n)/L_Y$ 总和，允许素数端点在任一侧。它由真孪生行和恰有一端为两个不同素数乘积的行组成。

写 $D_{3,k}(n)=D(0,0,h_{k,n};N)$，其中 $h_{k,n}=2\pi k/\log n$。在 $\omega(N)=2$ 时，条件独立符号模型在同时翻转所有符号下不变，所以三阶联合矩及累积量均为零。在 $\omega(N)\ge4$ 时用定理 1。在这两类上 $D_{3,k}=O_A(P^{-A})$。若 $\omega(N)=3$，将素数端记为 $p$，另一端记为 $m=st$、$s<t$，并令 $\theta=\log s/\log m$。取 $k=0$ 时，(7) 给

$$
\frac{D_{3,0}(n)}{\log n}
=w(n)\big(6\theta(1-\theta)+O(n^{-3/11})\big).
\tag{10}
$$

因为 $\sqrt{p s t}=\sqrt{n(n+2)}$，且 $\log p,\log m$ 与 $\log n$ 都相差 $O(1/n)$，两种左右次序都适用。在 $\mathcal R_Y$ 中 $\theta\ge3/11-O(1/(n\log n))$。因而半素数行的三阶权相对 $w(n)$ 落在 $[144/121-o(1),3/2+o(1)]$。特别地，同一统计赋予两个因子接近平方根的半素数更大的权。

可以用零点本身消去这一权重差异。对有界 $h_{k,n}$，排列和中第三个位置可以放三个因子之一，每项重复两次。因此在 $\omega(N)=3$ 的行上

$$
\frac{D_{3,k}(n)}{D_{3,0}(n)}
=\frac13\left(1+2\cos(2\pi k\theta)\right)
+O\left(n^{-3/11}+\frac{|k|}{n\log n}\right).
\tag{11}
$$

$p$ 的归一化对数与 $1$ 相差 $O(1/(n\log n))$，另两个分别与 $\theta,1-\theta$ 相差同阶。其余奇素数幂在每个因子的系数中给出 $O(n^{-3/11})$ 的相对尺度误差。$D_{3,0}\gg w(n)\log n$ 使 (11) 的除法合法。单取 $k=1$，右边在 $0<\theta\le1/2$ 上单调，已经保留了因子比例的信息。

选一个光滑、偶、周期为 $1$ 的实函数 $F$，使 $F(0)=0$，在 $[2/11,9/11]$ 上等于 $1/(4u(1-u))$。例如乘一个在距整数 $1/11$ 内等于零、距整数至少 $2/11$ 时等于一的光滑截断。写 $F(u)=\sum_{k\ge0}c_k\cos(2\pi ku)$，$K_0=\lfloor\sqrt{\log P}\rfloor$。仅需 $0\le k\le K_0$，所有第三个窗心均为 $O(1/\sqrt{\log P})$。同一行相邻窗心距离 $2\pi/\log n>2R$。定义

$$
\mathcal H(Y)=\frac1{L_Y}\sum_{n\in\mathcal R_Y}
\frac{I_{Q,d_0}(n)}{\log n}\sum_{k=0}^{K_0}c_kD_{3,k}(n).
\tag{12}
$$

**定理 2。** 对每个固定 $L>0$，任意原好见证及其窗口均有

$$
\mathcal H(Y)=A(Y)-J_{Q,d_0}(Y)+O_L((\log P)^{-L}).
\tag{13}
$$

证明时把 (11) 对 $k$ 加权。因 $F$ 光滑，$\sum |c_k|(1+k)<\infty$，在 $K_0$ 截断的误差为任意对数负幂。全级数的乘数是 $(F(0)+2F(\theta))/3=1/(6\theta(1-\theta))$，恰好约去 (10) 的因子。对素因子数不等于三的行，用三阶差的负幂界。最后对行求和，即使用 $\sum_{n\in I_Y}w(n)/L_Y\ll_J(\log P)^2$ 的粗界也够，只要对 Fourier 尾预留两个对数幂。式 (12) 的计算只使用粗筛、平方自由性、Jacobi 符号及零点。它也给出物理窗口内半素数与真孪生贡献的同一套对数权。

### 每个二倍区间只用三个固定窗心的上界

对于改进孪生素数下界，尤其需要防止低估半素数的质量。这里有一个只用三项、相对余量小于 $0.4\%$ 的上界。在每个标准二倍区间 $(U,2U]$ 上，令 $\widetilde D_{3,k}(n)=D(0,0,2\pi k/\log U;N)$。在 (12) 中以

$$
\frac{752\widetilde D_{3,0}(n)+465\widetilde D_{3,1}(n)
+72\widetilde D_{3,2}(n)}{1000}
$$

代替整段 Fourier 和，记所得统计为 $\mathcal H_{\rm up}(Y)$。端部不完整的二倍区间也用其标准左端。则

$$
A(Y)-J_{Q,d_0}(Y)-O_A(P^{-A})
\ \le\ \mathcal H_{\rm up}(Y)\
\le\ \frac{251}{250}\bigl(A(Y)-J_{Q,d_0}(Y)\bigr)+O_A(P^{-A}).
\tag{13a}
$$

三个窗心为 $0,2\pi/\log U,4\pi/\log U$，相邻零点窗不交。零点数据可以先于具体行一次取齐。证明只用下面一个实变量不等式。令

$$
W(v)=6v(1-v)\left(\frac{752}{1000}
+\frac{465}{3000}(1+2\cos(2\pi v))
+\frac{72}{3000}(1+2\cos(4\pi v))\right).
$$

在 $3/11\le v\le1/2$ 上有 $1.00044<W(v)<1.00372$。这个数值不等式用了可复核的区间证明。以 $v_j=3/11+j\cdot5/(22\cdot4096)$ 等分闭区间，对每个完整区间 $[v_j,v_{j+1}]$ 用 90 比特 Arb 球算术代入上式，对两个有理界逐段作严格比较。全部 4,096 段均通过。由 (10)、(11)，若先用 $\log n$ 作为频率分母，每个三因子行的相对权为 $W(\theta)+o(1)$。换成 $\log U$ 只再改动 $O(1/\log U)$，因为 $|\log n/\log U-1|\le\log2/\log U$。当 $\theta$ 比 $3/11$ 低 $O(1/(n\log n))$ 时由连续性处理。上述两个严格余量最终吸收这些变化，逐个三因子行都有 $1\le\text{相对权}\le251/250$。其余行的三阶差仍为任意负幂，总误差给出 (13a)。[包含完整区间证明的程序](https://github.com/lixiang90/vibemath-NumberTheory/blob/main/zeros-exact-distribution/gap-two-quadratic-conductor-covariance/verify.py)

## 6. 原来的好窗口上已经有哪一个正下界

**定理 3。** 存在固定 $c>0$。每个充分大的 $P$，有正比例的原两阶段候选对可作好见证。对其中这些见证的**每个**原好窗口，同时有

$$
A(Y)\ge c,\qquad
\mathcal V(Y)+\mathcal H(Y)=A(Y)+O_L((\log P)^{-L}).
\tag{14}
$$

这里的正比例对两阶段概率而言，先均匀取 $Q$，再均匀取 $d_0\in\mathcal D_Q$。还可以在预先列出的有限个窗口中选出一段，使上述结论对正比例见证在这段共同成立。

**证明。** 对每个完整的标准二倍区间 $(U,2U]\subset[P,H]$，令 $\mathcal C_U$ 是其中左端为素数、右端为素数或两个大于 $n^{3/11}$ 的不同素数之积的行。Green–Tao 定理 6.1 给 $\#\mathcal C_U\gg U/\log^2U$。剔除右端为平方的行至多花 $O(\sqrt U)$。固定一行时，两个端点有非空、互不相交的素因子集。先在第二阶段用 (2)，再在第一阶段用它的模 $8$ 版本，可知四道符号同时通过的概率为 $1/16+O_A(P^{-A})$。这是两个端点在每个辅助素数处的两道条件。若两行的所有素因子互不相交，它们同时通过的概率同理为 $1/256+O_A(P^{-A})$。各误差对所用行一致，所需物理标签乘积仅是 $P$ 的固定次幂。[二倍区间的陈氏下界及因子阈值](https://arxiv.org/pdf/math/0405581#page=19)

共享因子的两行须单独计数。每行至多涉及三个素数 $s>U^{3/11}$，共享 $s$ 的另一行须满足 $m\equiv0,-2\bmod s$，合计每行 $O(U^{1-3/11})$ 个可能。于是

$$
\mathbb E_{Q,d_0}\left|
\frac{\log^2 U}{U}\sum_{n\in\mathcal C_U}
 \left(I_{Q,d_0}(n)-\frac1{16}\right)
\right|^2
\ll (\log U)^4 U^{-3/11}+O_A(P^{-A}).
\tag{15}
$$

对 $O_J(\log P)$ 个标准二倍区间取并，以概率 $1-o(1)$，每个区间筛后都至少留下 $c_1 U/\log^2U$ 行，$c_1>0$ 固定。由原来的两阶段保留比例，可作见证的解析好对有固定正比例。与这个概率 $1-o(1)$ 的事件相交，仍有固定正比例。对每个这类好对，原证明可以再取几何好窗口。每段 $(2Y,Y^{3/2}/2]$ 含 $\gg\log Y$ 个完整标准二倍区间，每段通过的陈氏行贡献 $\gg1$ 的 $w(n)$ 总和，除以 $L_Y\asymp\log Y$ 即得 (14)。若需要同一段供多对使用，第 18 篇用原有限例外扭转指标列出固定数目的 $[Y_i,Y_i^2]$，每个解析好对至少有一段合格。有限覆盖给出共同段上的正比例。该证明先在全部原候选上做 (15)，再相交好事件。式 (7) 的观测平均仍取全 $\mathcal D_Q$。[两阶段保留比例](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/MeanSquareBounds/PairShare.lean) · [第 18 篇第 4 节的共同窗口覆盖](https://github.com/lixiang90/vibemath-NumberTheory/blob/main/Goldbach_TwinPrimes/twin-prime-auxiliary-signs.md)

因此可以精确叙述下一项需要的估计。在 (14) 给出的窗口中，若存在固定 $\eta>0$，使无穷多个窗口满足

$$
\mathcal H(Y)\le (1-\eta)A(Y),
\tag{16}
$$

则 (13) 给 $\limsup J_{Q,d_0}(Y)\ge\eta c$，从而 $\Gamma_1\ge\eta c/3$。由 (13a)，一个只涉三个共享窗心的充分条件是 $\mathcal H_{\rm up}(Y)\le(1-\eta)c$，其中 $c$ 取 (14) 已保证的下界。反过来，若这些窗上的真孪生质量趋零，(13) 迫使经因子校准的三阶统计独自达到 $A(Y)+o(1)\ge c+o(1)$。目前的平均特征和界能证明 (14)，还不能在无穷多个窗口上证明 (16)。正质量要求也强于单纯的孪生素数无穷性。

## 7. 同一族不附加算术条件时的零点分布

选 $P\le U\le H/2$、$h=\pi/\log U$。在整个 $\mathcal D_Q$ 上，两段零点的联合极限为

$$
\left(\frac{T_d(0)}{B},\frac{T_d(h)}{B}\right)
\ \Longrightarrow\
\mathcal N(0,1/2)\otimes\mathcal N(0,1/4).
\tag{17}
$$

两个正态变量独立，记号的第二个参数是方差。中心窗还有一个可定量的非高斯项。令 $\varepsilon$ 为一个均匀正负号，对固定 $j\ge2$，

$$
\kappa_{2j}(T_d(0))\sim
\frac{\kappa_{2j}(\varepsilon)}{j-1}
\frac{(\log z)^{2j-1}}{z^{j-1}},\qquad
\kappa_4(T_d(0))\sim-\frac{2(\log z)^3}{z}.
\tag{18}
$$

所有固定奇阶累积量均为 $O_A(P^{-A})$。式 (17) 是两个收缩区间上**平滑零点和**在不同导子间的分布。式 (8) 则描述加入一个算术条件后，这两个和的协方差在远小于主方差的尺度如何改变。素数导子二次特征族的低位零点单层密度已有 Andrade–Baluyot 的结果，这里的族还带随 $P$ 增长的强制同余，观测处于含许多平均间距的尺度。[Andrade–Baluyot 的单层密度](https://arxiv.org/abs/1802.03413) · [Hughes–Rudnick 的低位零点矩方法](https://arxiv.org/abs/math/0208230)

**证明。** 由第 3 节的矩比较，可在独立模型中计算任意固定阶。令 $x=\log p$。协方差及两个方差分别以

$$
\int_0^{3B/2} x g(x)^2 f(x)\,dx,\qquad
f(x)=\cos(hx),\quad 1,\quad \cos^2(hx)
$$

为主项，误差 $o(B^2)$。说明这个近似的一致性如下。略去 $p\le z$ 花费 $O((\log z)^2)$，奇幂修正由 $\sum_p(\log p)^2/p^2$ 控制。在 GRH 的素数定理下，分部求和给 $\sum_{p\le e^x}(\log p)^2/p=x^2/2+O(1)$。例如对 $\vartheta(v)=v+O(\sqrt v\log^2(2v))$ 以 $(\log v)/v$ 作 Stieltjes 积分，误差积分绝对收敛。函数 $g$ 在正轴单调递减（它是对称单峰密度与中心区间的卷积），故换为上述连续积分的额外误差为 $O(1+|h|B)$。

把核中的每个均匀卷积视为独立变量，则这 $2K$ 个变量的和除以 $B$ 后方差为 $1/(24K)$。故 $g(Bu)\to\mathbf1_{(-1,1)}(u)$ 在 $u\ne\pm1$ 处成立，支配收敛允许在积分中换成区间，代价 $o(B^2)$，对 $|f|\le1$ 一致。因 $hB\to\infty$，$\int_0^B x\cos(hx)\,dx=o(B^2)$，换 $2h$ 同样成立。三个极限系数分别为 $0,1/2,1/4$。任意固定线性组合在独立模型中是中心独立和，单个标准化系数最大值为 $o(1)$，各高阶标准化累积量趋零。矩比较和正态分布的矩确定性给 (17)。

对中心窗，独立和的 $2j$ 阶累积量等于 $\kappa_{2j}(\varepsilon)\sum_{p>z} a_0(p)^{2j}$（支撑外系数为零）。主要素数在 $z$ 以上且远低于 $e^{B/2}$。在平台内 $a_0(p)=(\log p/\sqrt p)(1+O(1/p))$，平台以外对 $j\ge2$ 是指数小的尾。素数定理及分部求和给 $\sum_{p>z}(\log p)^{2j}/p^j\sim(\log z)^{2j-1}/((j-1)z^{j-1})$。式 (6) 的精度可取到小于这个主项的任意固定幂尺度，得 (18)。

## 8. 实际零点与独立算术复算

在有限实验中取 $d\equiv23\bmod48$，对应 $d\equiv7\bmod16,\ d\equiv-1\bmod3$。取有限核 $B=3.45,b=0.75,K=2$，半宽 $6.2$，两个窗心为 $0$ 和 $h=13\pi/\log\sqrt{17\cdot19}=14.13747391\ldots$。采用奇数倍使有限范围内两窗已经分离，同时 $(17,19)$ 两端的余弦均为负。数值只检验同余、核和公式。渐近定理的 $P,Q,Y$ 与几何见证由前面各节处理。

对 $800<d\le4200$ 的全部 27 个此类导子，独立定位 $|t|\le20.76748$ 的 1,416 个零点。计算用奇二次特征的 theta 级数 $\Theta_d(x)=\sum_{m\ge1}m\chi_d(m)e^{-\pi m^2x/d}$。由其模变换，临界线上完备 $L$ 函数等于 $2\int_0^\infty\Theta_d(e^u)e^{3u/4}\cos(tu/2)\,du$。实零点变号定位，并与两种水平边网格的辐角计数逐导子比对。在三个导子的共九个临界线外测试点另与 Arb 的 Dirichlet $L$ 值对照，最大绝对差为 $3.42\times10^{-10}$。这些是浮点核对。

直接计算 (4)，再与**同一有限族中心化后**的算术和逐导子比较，两窗最大绝对差分别为 $0.000494,0.000262$。在 $N=17\cdot19$ 的两组中，负协方差半差从零点得到 $0.14245046$，逐素数计算为 $0.14251374$。独立正负号的渐近主项为 $0.35312774$。小族中真实的特征相关可以显著影响单行信号。

再取 $\sqrt{12000000}<d\le12000000$ 的全部 49,308 个同余素导子，用算术式计算有限族的二、三阶半差。第一组的两个端点都是素数，独立模型的三阶参考值为零。第三组两个端点都为合数，二、三阶参考值同时为零。

| $n,\ n+2$ | 端点分解 | 负协方差半差 | 独立参考值 | 中心窗三阶半差 | 独立参考值 |
| --- | --- | ---: | ---: | ---: | ---: |
| $17,19$ | $17,\ 19$ | 0.35673 | 0.35313 | −0.18286 | 0 |
| $53,55$ | $53,\ 5\cdot11$ | −0.00595 | 0 | 0.37063 | 0.42074 |
| $143,145$ | $11\cdot13,\ 5\cdot29$ | −0.00358 | 0 | −0.15378 | 0 |

程序另在八个标签的 256 种符号赋值上穷举，核对三行的条件协方差及三阶累积量分解。用 512 个余弦谐波复核 (12) 的因子校准，201 个因子比例上的最大归一化误差为 $3.4\times10^{-12}$ 以内。在 $10^6$ 以上两种左右次序各取 40 个实际半素数行，用一次素数幂主项复核 $\log n$ 和 $\log(n+2)$ 的尺度差，最大相对误差分别为 $8.6\times10^{-7}$ 和 $1.1\times10^{-6}$ 以内。复算结果给出 (8)、(11) 和 (13) 的符号与权重检查，式 (16) 的无穷多个窗口上界仍待证明。[复算程序](https://github.com/lixiang90/vibemath-NumberTheory/blob/main/zeros-exact-distribution/gap-two-quadratic-conductor-covariance/verify.py)
