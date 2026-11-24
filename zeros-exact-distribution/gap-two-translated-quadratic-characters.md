---
title: "在一个固定补弧频率上分离孪生素数的零点统计"
date: "2026-11-25"
abstract: "固定 Alpöge–Claude 已选的辅助素数对和物理窗口。在同一个收缩零点窗中，用平移的二次特征给两个端点配对。在固定无理频率上，四个有界余弦平均的对数比隔离出真实孪生质量，陈氏候选的相应统计在每个大二倍区间有正下界，并给出线性零点投影的负四阶累积量及其量级。反向赋号可以使差二质量为零，同时保留主弧主项以及其他各差在给定精度下的相关量。"
---

# 在一个固定补弧频率上分离孪生素数的零点统计

第 31 篇从同一导子的一小段零点同时读出所有物理端点，并算出主弧上固定差二的贡献 $\mathfrak S(2)U+O(U/\log^2U)$。剩下的问题是补弧上的抵消能否留下固定余量。本篇先做一个检验。在真正的素数上加非负权，可以让所有孪生行消失，同时几乎保持整个主弧积分和其他每个非零差的相关和。这说明还需要对哪种误差作出改进。[前篇的主弧与原好窗口公式](https://lixiang90.github.io/vibemath-NumberTheory/zeros-exact-distribution/articles/gap-two-auxiliary-conductor-bands.html)

再给出一项保留差二的新观测。固定一对已经通过原选窗条件的辅助素数 $Q,D$。用 $Q$ 的零点测量 $6m-1,6m+1$ 两端，对两端都乘上同一个平移符号 $\big(\frac{t+m}{D}\big)$。在固定频率 $\alpha_0=\sqrt2-1$ 上，式 (7) 的四个有界平均能分离孪生质量。一个常数级的质量差可以用同一辅助对的 $O(\log P)$ 个平移来分辨。把右端允许的类型从素数扩展到陈氏半素数，便能在每个大二倍区间证明这个观测有正下界。真素数右端所需的严格下界仍然开放。

## 1. 固定原窗口与零点数据

沿用原 ζ 构造。$P$ 是基准，$H=P^J$，$0<\epsilon<1/52$ 和 $J$ 固定。令 $z=P^{\epsilon/5}$，$z'=z\log^2P$。第一素数满足 $z/2<\log Q\le z$，$Q\equiv3\bmod8$。第二素数满足 $z'/2<\log D\le z'$，$D\equiv7\bmod16$，且模 $Q$ 及模每个不超过 $z$ 的奇素数余 $-1$。先取原好对和所需的 $[Y,Y^2]\subset[P,H]$，全篇的有限域变量 $t$ 随后引入。它是本篇新增的观测变量，几何证书始终使用原来的 $Q,D,Y$。[两阶段素数条件](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/MeanSquareBounds/PrimeSetTwoBounds.lean#L37) · [先取素数对再选完整窗口](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/CharacterBoxes/WinPairs.lean#L32)

用 $b_n$ 表示前篇式 (5) 的素数指标的实部。具体地，$\chi$ 遍历模 $Q$ 的奇特征，$\tfrac12+i\gamma_\chi$ 为非平凡零点，按重数计。以 $S_\chi=\sum_{|\gamma_\chi-h|\le R}\widehat g(\gamma_\chi-h)$ 表示其窗内线性统计，

$$
Z(r)=-\frac{2}{Q-1}\sum_{\chi(-1)=-1}\overline{\chi(r)}S_\chi,\qquad
b_n=\frac1{\log n}\Re\sum_{r^j=n}\mu(j)\sqrt r\,r^{ih} Z(r).
\tag{1}
$$

根取 $r\ge2,j\ge1$。前篇以 Dirichlet-GRH 为输入，给出具体的核 $g$ 和半宽 $R=O(\log P/\log Q)=o(1)$。对任意预定 $A$ 均可选核，使 $b_n=\mathbf1_{\mathbb P}(n)+O(P^{-A})$ 对 $2\le n\le H+2$ 同时成立。下文可按所需精度增大 $A$，每项构造使用同一个窗心和同一批零点。本文后续推导没有作 Lean 形式化。[窗核、尾估计与整数根反演](https://lixiang90.github.io/vibemath-NumberTheory/zeros-exact-distribution/articles/gap-two-auxiliary-conductor-bands.html) · [归档的 Dirichlet-GRH 陈述](https://github.com/leanprover-community/riemann-hypothesis#theorem-statements)

记 $\ell_m=6m-1,\ r_m=6m+1$。对左端属于给定物理区间的行，先筛去 $\ell_m r_m$ 含有不超过 $z$ 的素因子的行，保留指标集 $\mathcal J$。所有大于 $3$ 的孪生素数对都在相应的保留行中。对 $\ell_m\in(U,2U]$、$P\le U\le H$ 有

$$
\#\mathcal J\ll \frac{U}{(\log z)^2}\ll_{\epsilon,J}\frac{U}{(\log U)^2}.
\tag{2}
$$

这是二维上界筛。在变量 $m$ 上，每个 $5\le p\le z$ 恰有两个禁类。算术大筛的分母至少为常数乘 $\sum_{d\le z,(d,6)=1}\mu(d)^2 2^{\omega(d)}/d\gg\log^2 z$，分子为 $O(U+z^2)$。分母的对数平方也可由其 Dirichlet 级数除以 $\zeta(1+s)^2$ 得到，余下 Euler 乘积在 $s=0$ 绝对收敛且为正。[Montgomery–Vaughan，推论 1](https://personal.science.psu.edu/rcv4/personal/Publications/large_sieve.pdf#page=2)

允许确定的行权 $w_m\ge0$，要求 $\sum_{\mathcal J}w_m\le C_0$，$C_0$ 与 $P$ 无关。权和区间都可由原好对决定。在二倍区间可取 $w_m=(\log U)^2/U$，或取 $w_m=\log\ell_m\log r_m/U$。在原窗口上取

$$
2Y<\ell_m\le Z_Y=Y^{3/2}/2,\qquad
w_m=\frac{\log\ell_m\log r_m}{\ell_m L_Y}I_{Q,D}(\ell_m),\qquad
L_Y=\sum_{Y\le n\le Z_Y}\frac1n .
\tag{3}
$$

此处 $n$ 取整数，$I_{Q,D}$ 要求两端在 $Q$ 处都是二次剩余、在 $D$ 处都是二次非剩余。由 (2) 分二倍段相加，(3) 的总权也有界。若不作四符号筛选，取 $I_{Q,D}=1$。

## 2. 一项可以逐个平移复核的恒等式

记 $\sigma_m(t)=\big(\frac{t+m}{D}\big)$，零处取通常的值 $0$，$\mathbb E_t$ 表示对 $\mathbb F_D$ 的均匀平均。不同 $m\in\mathcal J$ 在模 $D$ 下互异。对任意复系数 $x_m,y_m$ 有

$$
\mathbb E_t\!\left[\sum_m x_m\sigma_m(t)\,
           \overline{\sum_m y_m\sigma_m(t)}\right]
+\frac1D\left(\sum_m x_m\right)\overline{\left(\sum_m y_m\right)}
=\sum_m x_m\overline{y_m}.
\tag{4}
$$

对角平均为 $1-1/D$。异标 $m,n$ 的乘积和为 $-1$，因为 $t\ne-n$ 时 $(t+m)/(t+n)$ 遍历除 $1$ 外的全域，而 $\sigma_n(t)^2=1$。这同时给出非对角修正的正号。取 $x_m=\sqrt{w_m}\,b_{\ell_m}e(\ell_m\alpha)$，$y_m=\sqrt{w_m}\,b_{r_m}e(\ell_m\alpha)$，其中 $e(u)=e^{2\pi iu}$。右端已经将右端点的相位反移了 $2$。于是任意 $\alpha$ 上都只剩 $\sum w_m b_{\ell_m}b_{r_m}$，以任意负幂精度等于真孪生质量。

特别可固定 $\alpha=\alpha_0$。对于一切有理 $a/q$，$\|\alpha_0-a/q\|_{\mathbb R/\mathbb Z}\gg q^{-2}$。将 $a/q$ 平移到距 $\alpha_0$ 最近的代表后，对其代入整数二次多项式 $x^2+2x-1$ 即得。对 $q\le\sqrt U/(\log U)^8$，这个距离超过前篇主弧的半宽 $(\log U)^2/U$。同一个 $\alpha_0$ 因此属于每个充分大二倍段的补弧。这里的单频观测同时用到了行标签、相位配准和平移平均。

## 3. 用有界的零点观测分离同一行的两端

用 $u_{m,1}=b_{\ell_m},u_{m,2}=b_{r_m}$。固定 $\alpha$，在实部和虚部两个坐标取

$$
c_{m,1}=\cos(2\pi\ell_m\alpha),\qquad c_{m,2}=\sin(2\pi\ell_m\alpha),\qquad
s_{m,j}=\operatorname{sgn}(c_{m,j})\arctan\sqrt{\tanh(w_mc_{m,j}^2)}.
\tag{5}
$$

$s_{m,j}$ 满足 $\cos(2s_{m,j})=\exp(-2w_mc_{m,j}^2)$ 和 $s_{m,j}^2\le w_mc_{m,j}^2$，这些系数可以先于素性判定算出。设

$$
L_j(t)=\sum_m s_{m,j}u_{m,1}\sigma_m(t),\qquad
R_j(t)=\sum_m s_{m,j}u_{m,2}\sigma_m(t),\qquad
C_j^\pm=\mathbb E_t\cos(L_j(t)\pm R_j(t)).
\tag{6}
$$

**定理 1。** 对每一原好对，每个这样的区间、权和频率，一致有

$$
\mathcal K(w,\alpha)=\frac12\sum_{j=1}^2\log\frac{C_j^-}{C_j^+}
=\sum_{m\in\mathcal J}w_m\,\mathbf1_{\mathbb P}(\ell_m)\mathbf1_{\mathbb P}(r_m)
   +O_A(P^{-A}).
\tag{7}
$$

四个单次观测都在 $[-1,1]$。四个平均在充分大时均至少为 $\tfrac12e^{-2C_0}$。故稀疏的正质量也可以同零质量区分，误差中没有一个随半素数数目增长的方差项。

**证明。** 先将两个 $u$ 换为精确的素数指标，再将各 $\sigma_m$ 换为独立均匀正负号。此时 $C_j^\pm$ 成为各 $\cos(s_{m,j}(u_{m,1}\pm u_{m,2}))$ 的乘积。只有一端是素数时，该行分子分母相同。两端都是素数时，该行对半对数比贡献 $-\tfrac12\log\cos(2s_{m,j})=w_mc_{m,j}^2$。两个坐标相加得到 $w_m$。各乘积至少为 $\prod_m\cos(2s_{m,j})\ge e^{-2C_0}$。

为把独立符号换回真正的平移，考虑一个实线性和 $S=\sum_\nu a_\nu\sigma_\nu(t)$，其中标签模 $D$ 互异。次数 $r\ge1$ 的任一单项式，其期望与独立符号相差 $O(r/\sqrt D)$。将奇数重标签作成平方自由多项式 $\prod(X+\nu)$ 后用 Weil 界，偶数重标签在自己的零处造成至多 $r/D$ 的误差。全部重数均为偶数时只需后一项。因此 $r$ 阶矩的比较误差至多 $O(rD^{-1/2}\|a\|_1^r)$。这是 Davenport–Erdős 平移矩方法在带权行上的应用。[Weil，On Some Exponential Sums](https://www.pnas.org/doi/10.1073/pnas.34.5.204) · [Davenport–Erdős 的平移特征矩](https://combinatorica.hu/~p_erdos/1952-11.pdf) · [Perret-Gentil 的特征函数比较](https://arxiv.org/abs/1609.09334)

截断 $e^{iS}$ 的 Taylor 式到 $2k-1$ 次，实变量上的余项由 $|S|^{2k}/(2k)!$ 控制。独立符号的矩母函数至多为 $\exp(t^2\|a\|_2^2/2)$，故矩至多为 $(Ck\|a\|_2^2)^k$。对余项也用上述矩比较，就得

$$
\left|\mathbb E_t e^{iS}-\prod_\nu\cos a_\nu\right|
\ll \left(\frac{C\|a\|_2^2}{k}\right)^k
+\frac{k^2}{\sqrt D}(1+\|a\|_1)^{2k}.
\tag{8}
$$

取 $k=\lceil(\log P)^3\rceil$。只要 $\|a\|_1\le P^{O(1)}$、$\|a\|_2^2=O((\log P)^2)$，右端就比 $\exp(-B(\log P)^2)$ 还小，$B$ 可任意固定。在 (6) 中平方范数至多 $4C_0$。把素数指标换回线性零点和时，$\sum|s_{m,j}|\le\sqrt{H C_0}$，用 (1) 和余弦的 Lipschitz 界可把附加误差压到 $O_A(P^{-A})$。对有正下界的平均取实对数即证 (7)。

这同时给出零点投影的分布。令 $\mathbf X=(L_1,R_1,L_2,R_2)$，取有界的固定实向量 $v$，写 $v\cdot\mathbf X=\sum_m d_m\sigma_m(t)$。若 $w_*=\max_m w_m=o(1)$，则

$$
\log\mathbb E_t\cos(v\cdot\mathbf X)
=-\tfrac12\sum_m d_m^2-\tfrac1{12}\sum_m d_m^4+O_v(w_*^2)+O_A(P^{-A}).
\tag{8a}
$$

事实上，$\sum d_m^2=O_v(1)$，$\max|d_m|=O_v(\sqrt{w_*})$，在 (8) 后展开 $\log\cos d_m$ 即得。正弦平均同样由 (8) 控制到任意负幂。因而在行系数给出的协方差矩阵收敛的每条子列上，四个零点投影联合趋于相应的中心高斯分布，允许退化。下一节在原素数抽样上也能应用这个结论。

## 4. 原来的第二阶段选素数也能给出这种矩比较

原归档的 $\mathrm{twistK}$ 在固定幂范围内给平方自由标签的负幂误差，允许一个例外。**在归档的最终 Dirichlet-GRH 结论之上**，对这里的第二阶段候选族可以去掉这个例外，并让次数随 $P$ 增长。以下是返回原选素数步骤所得的推论。[原引理及其例外的量词](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/MeanSquareBounds/twistK.lean#L145)

给固定的 $Q$，设 $X=e^{z'}$、$M=16Q\prod_{3\le p\le z}p$，其中 $p$ 遍历奇素数。第二阶段候选 $d$ 是 $(\sqrt X,X]$ 中单个可逆剩余类 $a\bmod M$ 内的素数。给不超过 $P^B$ 个互异奇素数标签 $r_\nu\le P^B$，均与 $M$ 互素，其中 $B$ 固定。对于每个非空子集 $S$ 一致有

$$
\left|\mathbb E_{d\mid Q}\prod_{\nu\in S}\left(\frac{r_\nu}{d}\right)\right|
\le \exp\big(-(1/2-o(1))z'\big).
\tag{9}
$$

写 $N=\prod_{\nu\in S}r_\nu$。二次互反律的符号在固定剩余类内不变。用模 $M$ 的全部特征展开这个类，得到特征 $\xi(\cdot)(\cdot/N)$，每个都非主，因为 $(M,N)=1$ 且 $N>1$ 平方自由。Dirichlet-GRH 的素数特征和界及分部求和给出分子 $O(\sqrt X\log^2(XMN))$，展开类时的 $1/\varphi(M)$ 抵消了特征个数。分母 $\gg X/(\varphi(M)\log X)$，这是对主特征分出主项、其余特征用同一界的结果。最后 $\log M=O(z)$，$\log N=O_B(P^B\log P)$，故 $\log(\varphi(M)\log X\log^2(XMN))=o(z')$。即使标签乘积超过任一固定的 $P$ 幂次，(9) 仍成立。这里没有在平方根以上水平使用素数在每个细剩余类中的相对误差。[Lee，GRH 下的特征素数和与等差数列素数定理](https://arxiv.org/abs/2301.13457)

例如可用 $H+2$ 以后的不同素数给各行编号，所需标签都在 $O(H^2)$ 内。若各行系数在抽 $d$ 前给定，用 $\big(\frac{r_\nu}{d}\big)$ 代替平移符号也有 (7)，只须在 (8) 的单项式比较中使用 (9)。式 (9) 取整个第二阶段候选族的平均。若已经看过 $D$ 才选取符号门、几何窗口或系数，式 (6) 的模 $D$ 平移仍然逐对有效。这样保留了原共同窗口的选择次序。

## 5. 多少个实际平移足以分辨信号

完整模 $D$ 平均可能非常大。取整数步长 $B_0=\lceil H\rceil+3$，对一个均匀原点 $t_0\bmod D$ 使用 $t_0+iB_0,\ 0\le i<N$。设 $N\le P^b$，这里及下文的指数界固定。将 (6) 的四个完整平均改为这 $N$ 个值的经验平均。对于预先列出的至多 $P^{b_1}$ 份观测，给定 $A>0$，以至少 $1-P^{-A}$ 的概率同时有

$$
\max_{j,\pm,\,\text{所列观测}}
|\widehat C_j^\pm-C_j^\pm|
\ll_{A,b,b_1,C_0,J}\sqrt{\log P/N}.
\tag{10}
$$

这些观测可来自已经确定的同一原好对上的不同窗口与频率。可在看到平移后从清单中选一份。若 $N\ge C\log P$ 且 $C$ 足够大，四个经验平均都为正，半对数比有同阶误差。于是分辨任一固定的正质量差只需 $N=O_{\delta,A}(\log P)$，常数也依赖固定的原窗参数。若要从这个抽样统计分辨仅有一对的二倍区间，取 $w=(\log U)^2/U$ 时，充分的样本数为 $C_Aw^{-2}\log P$。单对检测和正密度检测使用不同的精度。

**证明。** 不同副本的标签集 $m+iB_0$ 互不相交，因为每组标签直径小于 $B_0$，而 $NB_0+H<D$ 最终成立。在独立符号模型中，各副本独立，有界经验平均减其期望后的 $2r$ 阶矩至多为 $(Cr/N)^r$，可由有界变量的指数矩界直接积分得到。对真实平移比较该矩时，把至多 $2r$ 个余弦的乘积展开成带正负相位的指数函数的平均。合并相同标签后，系数的平方范数是 $O(C_0r^2)$，一范数为 $P^{O(1)}$。取 $r=\lceil c\log P\rceil$，(8) 对每项都适用。中心化并展开平均后绝对系数和至多 $2^{2r}$，矩误差仍小于任意预定的 $\exp(-B(\log P)^2)$。因 $N\le P^b$，这个误差在 Markov 不等式除以 $(C'\sqrt{\log P/N})^{2r}$ 后仍可忽略。对清单取并即得 (10)。副本的间距是这里实际用到的条件。

## 6. 陈氏候选已经有正下界，真孪生行需要进一步分离

在 (6) 的右端，还可以放入下列由同窗零点构成的次数至多为二的量，

$$
b^{\rm Ch}_{r_m}=b_{r_m}
+\sum_{\substack{d\mid r_m\\\ell_m^{3/11}<d\le\sqrt{r_m}}} b_d\,b_{r_m/d}.
\tag{11}
$$

每个真实的半素数至多有一个这样的无序素因子分解，平方也只计一次。以足够精度选 (1) 的核后，(11) 一致近似于一个 $0,1$ 指标，对 $d$ 求和的误差可由 $O(H)$ 项的平凡界吸收。令 $\mathcal K_{\rm Ch}$ 使用 (11)，$\mathcal K_{\rm semi}$ 只使用其因子和，$\mathcal K_{\rm twin}$ 使用 $b_{r_m}$。对固定的 $\alpha_0$ 和每个上述窗口，(7) 的逐行证明给出

$$
\mathcal K_{\rm Ch}=\mathcal K_{\rm semi}+\mathcal K_{\rm twin}+O_A(P^{-A}).
\tag{12}
$$

取二倍段和 $w_m=(\log U)^2/U$。陈氏定理的粗素因子版本给每个充分大的 $(U,2U]$ 至少 $cU/\log^2U$ 个左端素数。因为 $z<U^{3/11}$，这些都通过 (2) 的小素数筛，并属于模 $6$ 余 $5$ 的左端。故有固定的 $c'>0$ 使 $\mathcal K_{\rm Ch}(w,\alpha_0)\ge c'$。在未加四符号门的原幂次窗采用 (3) 的对数权，同样由完整二倍段求和得到固定正下界。至于加上四符号门的真孪生行，原归一化给出 $\Gamma_1=\tfrac13\limsup\mathcal K_{\rm twin}(w,\alpha_0)$。式 (12) 中两项间的固定正差仍须从真实素数获得。[Green–Tao，定理 6.1 与脚注 3](https://arxiv.org/pdf/math/0405581#page=19) · [原共同窗口的四符号门及斜坡因子](https://github.com/lixiang90/vibemath-NumberTheory/blob/main/Goldbach_TwinPrimes/twin-prime-auxiliary-signs.md)

**一个可证的四阶修正。** 仍取这一二倍段权和 $\alpha_0$。在左端的两个线性零点投影 $L_1,L_2$ 中，选方差较大的一个 $L_j$。每个充分大的二倍段都有

$$
\operatorname{Var}_t L_j\asymp1,\qquad
-\kappa_4(L_j)\asymp\frac{(\log U)^2}{U},\qquad
\kappa_4(L_j)=\mathbb E_tL_j^4-3(\mathbb E_tL_j^2)^2 .
\tag{12a}
$$

两个投影均值严格为零。写 $d_m=s_{m,j}b_{\ell_m}$，二阶和四阶的逐单项式比较给 $\operatorname{Var} L_j=\sum d_m^2+O_A(P^{-A})$ 及 $\kappa_4(L_j)=-2\sum d_m^4+O_A(P^{-A})$。由 (5)，每行 $\sum_j s_{m,j}^2=w_m+O(w_m^2)$。陈氏下界使两坐标的方差之和有固定正下界，而 (2) 给出固定上界。对所选 $j$，Cauchy–Schwarz 给 $\sum d_m^4\ge(\sum d_m^2)^2/\#\mathcal J\gg(\log U)^2/U$，其上界来自 $\max d_m^2\ll(\log U)^2/U$。这证明了 (12a)，其中随机变量是同一批实际零点的平移投影。

## 7. 使补弧恰好抵消主弧的非负素数权

把上述构造用作检验，回到 $(U,2U]$ 内的真实素数权 $a_n=\log n\,\mathbf1_{\mathbb P}(n)$。用 (1) 舍入 $b_n$ 可得到相同的权。取所有与此区间相交且通过 (2) 的行。在每行左端，把 $a_{\ell_m}$ 乘 $1+\sigma_m(t)$。在右端分别乘 $1+\sigma_m(t)$ 或 $1-\sigma_m(t)$，得到两组非负权 $a_n^+,a_n^-$。其余端点保留 $a_n$。记 $F^\pm(\alpha)=\sum a_n^\pm e(n\alpha)$，$F=\sum a_ne(n\alpha)$。写 $C_d(a)=\sum_n a_na_{n+d}$，区间外的权取零。

**定理 2。** 给定 $A>0$，在每对原辅助素数上，除去比例 $O_A(P^{-A})$ 的原点 $t$，下列结论同时成立。常数可依赖固定的原参数。

$$
\begin{aligned}
&C_2(a^-)=0,\qquad
 C_2(a^+)=2C_2(a)+O_A(\sqrt U(\log U)^{3/2}),\\
&C_d(a^\pm)=C_d(a)+O_A(\sqrt U(\log U)^{3/2})\qquad(1\le d<U,\ d\ne2),\\
&\sup_{\alpha\in\mathbb R/\mathbb Z}|F^\pm(\alpha)-F(\alpha)|
 \ll_A\sqrt{U\log U}.
\end{aligned}
\tag{13}
$$

证明时先去掉任一 $\sigma_m(t)=0$ 的原点，比例至多 $H/D$。每个真正孪生行的异号乘积恒为零，同号乘积为 $2+2\sigma_m$。对其他差，两个端点一旦都有标记就属于不同的行。在独立符号模型中，它们的权积期望仍为原权积。对于固定 $d$，改变一行符号最多影响四个端点乘积，每项的变化为 $O(\log^2U)$。由 (2)，有界差分不等式的平方尺度为 $O(U\log^2U)$。例如依次揭示各行符号，用条件期望差的指数矩界相乘即得这一尺度的次高斯尾。对 Fourier 和的实部或虚部，单行改变至多 $O(\log U)$，平方尺度为 $O(U)$。

取 $2r$ 阶矩，$r=\lceil c_A\log P\rceil$。相关和关于符号的次数至多为二，其中心化矩是次数至多 $4r$ 的多项式，绝对系数和至多 $(C U\log^2 U)^{2r}$。第 3 节的逐单项式估计乘上此数仍小于任意固定负幂，甚至小于 $\exp(-B(\log P)^2)$。线性 Fourier 和同理。故独立模型的矩尾可转回真平移。对全部 $d$ 以及间隔至多 $U^{-3}$ 的圆周网格取并，再用导数界 $O(U^2\log U)$ 填满网格，即证 (13)。两组权在 $d=0$ 处的独立模型期望也相同，均为 $C_0(a)$ 加上已标记端点的 $V=\sum a_n^2=O(U)$。

对单端特征估计还可再作核对。固定 $B>0$，用同一高概率原点可同时得到

$$
\sup_{\substack{q\le U^B,\ \chi\bmod q,\ |u|\le U^B\\ I\subset(U,2U]\ {\rm interval}}}
\left|\sum_{n\in I}(a_n^\pm-a_n)\chi(n)n^{iu}\right|
\ll_{A,B}\sqrt{U\log U}.
\tag{13a}
$$

每行对固定测试的改变仍为 $O(\log U)$。所有模数的特征总数、整数端点区间数、以及间隔 $U^{-3}$ 的高度网格点数都是 $U$ 的固定幂以内。对它们取并后，用关于 $u$ 的导数界 $O(U\log^2U)$ 即得 (13a)。特别在 $u=0$，改权造成的扰动低于一般 Dirichlet-GRH 特征素数定理的 $O(\sqrt U\log^2(qU))$ 误差尺度。

令 $\mathfrak M(U)$ 为前篇分母不超过 $\sqrt U/(\log U)^8$、半宽为 $(\log U)^2/U$ 的主弧，$\mathfrak m(U)$ 为补弧。弧长总和 $\ll(\log U)^{-14}$，所以 (13) 给

$$
\int_{\mathfrak M}|F^\pm-F|^2\ll_A U(\log U)^{-13},\qquad
\left|\int_{\mathfrak M} (|F^\pm|^2-|F|^2)\cos(4\pi\alpha)\,d\alpha\right|
\ll_A U(\log U)^{-6}.
\tag{14}
$$

第二式用的是前篇 $\int_{\mathfrak M}|F|^2\ll U\log U$。因此对异号权有

$$
\frac1U\int_{\mathfrak M}|F^-|^2\cos(4\pi\alpha)\,d\alpha
=\mathfrak S(2)+O((\log U)^{-2}),\qquad
\frac1U\int_{\mathfrak m}|F^-|^2\cos(4\pi\alpha)\,d\alpha
=-\mathfrak S(2)+O((\log U)^{-2}).
\tag{15}
$$

两个积分精确互为相反数，因为全圆只读出 $C_2(a^-)=0$。可见，非负性、上述主弧精度以及同时掌握其他各差到 (13) 的精度，可以和完全抵消共存。这是对所列估计的检验，原来的素数指标还满足其自身的乘法和素性约束。[主弧常数与补弧单边条件](https://lixiang90.github.io/vibemath-NumberTheory/zeros-exact-distribution/articles/gap-two-auxiliary-conductor-bands.html)

**两组零点统计之间的精确差别。** 写 $B^\pm(\alpha)=\sum_{m\in\mathcal J}(a_{\ell_m}e(\ell_m\alpha)\pm a_{r_m}e(r_m\alpha))$，继续令区间外的 $a_n=0$。用 (4) 计算 $F^\pm-F$ 的平方均值，且 $\mathbb E_t(F^\pm-F)=0$，得到对所有 $\alpha$ 的有限恒等式

$$
\mathbb E_t|F^\pm(\alpha)|^2
=|F(\alpha)|^2+V\ \pm 2 C_2(a)\cos(4\pi\alpha)-\frac{|B^\pm(\alpha)|^2}{D}.
\tag{16}
$$

特别地，令

$$
\Delta(\alpha)=\mathbb E_t\big(|F^+(\alpha)|^2-|F^-(\alpha)|^2\big)
+\frac{|B^+(\alpha)|^2-|B^-(\alpha)|^2}{D}.
$$

取单个 $\alpha_0$ 时 $\cos(4\pi\alpha_0)>0$，由 (16) 和前篇的主弧估计，通常的素数指数和 $F$ 满足

$$
\frac1U\int_{\mathfrak m(U)}|F(\alpha)|^2\cos(4\pi\alpha)\,d\alpha
=-\mathfrak S(2)+\frac{\Delta(\alpha_0)}{4U\cos(4\pi\alpha_0)}
+O((\log U)^{-2}).
\tag{17}
$$

右边只需在同一指定频率比较两个带标记的二阶矩，且两组的共同背景 $V$ 已精确消去。各指数和也可直接代入未舍入的线性零点和 (1)，只多出任意负幂的误差。由 (17)，若 $\limsup\Delta(\alpha_0)/U>0$，则有无穷多个二倍区间含 $\gg U/\log^2U$ 对孪生素数。Weil 界控制不同标签的混合矩，而这个严格不等式取决于同一标签的两个真实端点。

## 8. 两种独立复算

取 $Q=107,D=449399$，均为素数，$Q\equiv3\bmod8$，$D\equiv7\bmod16$，且 $D\equiv-1\bmod(3\cdot5\cdot7 Q)$。以下只使用这些有限同余，原幂次窗口与导子的渐近尺度由前面的证明处理。

**先用实际零点。** 对模 $107$ 的全部 $53$ 个奇特征，在 $|\gamma|\le51.37$ 内数值定位 $5018$ 个零点。以 Hardy 函数变号求根，并逐特征用两种水平边网格的辐角计数交叉核对。取窗心 $1.25$、半宽 $48$，用前篇窗核的有限参数 $B=2.2,\ b=0.14,\ k=2$，对 $\ell_m,r_m=5,7,11,13$ 的剩余类作投影。算术支撑在 $n\le15$，程序逐标签核对只有直接项并除以已知的核值。四个素数指标的最大误差为 $0.000613$。对这两个完整的行取 $w=0.2$，只预筛 $2,3$。在 $\alpha_0$ 上，式 (4) 除以 $w$ 得 $2.00051391$，式 (7) 除以 $w$ 得 $2.00057734$，两者舍入都给出 $(5,7),(11,13)$ 共两对。这是浮点核对。

**再用完整有限域与独立素数表。** 对下表各区间要求两端都在 $(U,2U]$，筛至 $7$，取共同权 $w=0.2(\log U)^2/U$。对模 $D$ 的全部平移用循环卷积算出四个真实平均。与精确素数枚举以及独立正负号的余弦乘积对照，后者的半对数比除以 $w$ 在各行均与整数计数符合至 $10^{-9}$。

| $U$ | 保留行 | 左端是素数 | 右端是素数 | 真孪生数 | 真平移的 $\mathcal K/w$ | 能量差给出的 $C_2(a)/U$ |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 512 | 36 | 22 | 22 | 12 | 11.98698460 | 1.016098807 |
| 2048 | 147 | 84 | 77 | 45 | 45.01865698 | 1.413120794 |
| 8192 | 585 | 271 | 270 | 113 | 113.05956570 | 1.213779314 |

末列按照 (16) 直接对两份复指数和取模平方、作差，包含 $1/D$ 的有限修正。它和从实际孪生逐对相加所得的对数权在所示位数完全一致。在中间一组，陈氏和纯半素数观测除以 $w$ 分别给 $82.08112632,37.00839326$，对应整数为 $82,37$。在左端 $[2411,2537]$ 的子窗，两边分别有 $6,3$ 个素数，真孪生数为零，观测值为 $0.00000586$。使用 $4,16,32$ 个等间距平移时，遍历全部原点所得四个经验均值中的最大均方根误差分别是 $0.180467,0.090254,0.063650$。在左端两个坐标中，较大方差为 $0.240612311$，对应第四累积量为 $-0.002084004$，独立符号的值为 $-0.002033105$。

最后用 $t=1234$ 生成异号权，以便在物理区间上独立检查完全抵消。有限主弧取 $q\le\lfloor\sqrt U/5\rfloor$、半宽 $2/U$。程序逐个分母检查 $\alpha_0$ 属于补弧，再以全体差的自相关和每段弧的 Fourier 积分计算系数。

| $U$ | 原主弧 $/U$ | 原补弧 $/U$ | 异号主弧 $/U$ | 异号补弧 $/U$ |
| ---: | ---: | ---: | ---: | ---: |
| 512 | 1.573760 | −0.557661 | 1.425802 | −1.425802 |
| 2048 | 1.309892 | 0.103229 | 1.325588 | −1.325588 |
| 8192 | 1.224446 | −0.010667 | 1.171094 | −1.171094 |

有限数据复核了余弦比、归一化、差二相位以及完全抵消的符号。要求无穷多个原好窗口上有固定正余量，仍须证明 (7) 中真孪生质量或 (17) 中真实能量差的下界。[复算程序及运行方式](https://github.com/lixiang90/vibemath-NumberTheory/blob/main/zeros-exact-distribution/gap-two-translated-quadratic-characters/verify.py)

