---
title: "从同一导子高度 3 以内的零点恢复哥德巴赫与孪生素数计数"
date: "2026-11-25"
abstract: "使用 Dirichlet GRH，给出零点截高与辅助素数导子之间的显式折衷。对充分大的 X，在 2X² 与 4X² 之间选一个素数导子，其全部奇特征在高度 3 内的带标签零点同时确定截至 X 的每个哥德巴赫表示数和孪生素数计数。进一步增大导子可使绝对截高趋于零。"
---

# 从同一导子高度 3 以内的零点恢复哥德巴赫与孪生素数计数

设 $X$ 是充分大的整数，$2X^2\le Q\le4X^2$ 是任一素数。在 Dirichlet GRH 下，只用导子 $Q$ 的奇 Dirichlet 特征在 $|\Im\rho|\le3$ 内的非平凡零点，保留特征标签和重数，就能同时判定 $2\le a\le X$ 中哪些整数为素数。因此能精确恢复每个偶数 $n\le X$ 的哥德巴赫表示数和每个 $x\le X$ 的孪生素数对数。每个特征只需求一次同样加权的零点和。若取 $\log(Q/X)\asymp(\log X)^2$，截高可降至 $O(1/\log X)$。

文库已有方法用模 $n$ 的零点对检验指定窗口的哥德巴赫表示，另一种方法用一个候选素数端点作导子，在双对数高度检验其邻项。本文改变的是导子和截高的取舍。增大导子在对数坐标中腾出了一段没有其他同余代表元的区间，可以用它平滑截断。代价是要保留较大导子的全部 $(Q-1)/2$ 个奇特征。[固定偶数与低高度零点对](../Goldbach_TwinPrimes/goldbach-low-zeros.md) · [以孪生候选端点作导子的方法](gap-two-single-prime-conductor.md) · [Hughes–Rudnick 关于 Dirichlet 族低零点线性统计的工作](https://arxiv.org/abs/math/0208230)

## 1. 一个可计算的逐点误差界

以下只需整数 $X\ge3$ 和素数 $Q>2X+1$。记

$$
 L=\log X,\qquad H=\frac12\log\frac{Q-X-1}{X},
 \qquad W=L+2H=\log(Q-X-1).
$$

给定整数 $m\ge3$，置 $\operatorname{sinc}(v)=\sin(v)/v$，$\operatorname{sinc}(0)=1$，并定义

$$
 \Phi(t)=W e^{itL/2}\operatorname{sinc}(Wt/2)
                    \operatorname{sinc}(Ht/m)^m.
$$

令 $\chi$ 遍历模 $Q$ 的所有奇特征，记 $1/2+i\gamma_\chi$ 为 $L(s,\chi)$ 的非平凡零点。取正负两个方向的纵坐标、按重数求和，包含中心零点。定义

$$
 S_\chi(T)=\sum_{|\gamma_\chi|\le T}\Phi(\gamma_\chi),\qquad
 D_T(a)=-\frac{2\sqrt a}{Q-1}\sum_{\chi(-1)=-1}
                      \overline{\chi(a)}\,S_\chi(T). \tag{1}
$$

**定理 1.** 假设 Dirichlet GRH。对 $T\ge em/H$，同一组数 $S_\chi(T)$ 对全部 $2\le a\le X$ 满足

$$
 |D_T(a)-\Lambda(a)|\le
 B(X,Q,m,T)=40\sqrt X\,W(1+T)\log(4Q(2+T))\,e^{-m}. \tag{2}
$$

$D_T(a)$ 为实数。这里 $\Lambda(p^j)=\log p$，其余处为零。只要 $B<(\log2)/4$，就有

$$
 a\text{ 为素数}\quad\Longleftrightarrow\quad
 D_T(a)>\tfrac34\log a. \tag{3}
$$

因为对任何高次素数幂 $\Lambda(a)\le(\log a)/2$，这个判据也覆盖平方和更高次幂。在固定 $X,Q$ 下不断增大 $m$ 并取 $T=em/H$，即可用 (2) 检查所需精度。具体数值输入若带有每个 $S_\chi$ 的误差上界 $\epsilon_\chi$，在 $B$ 上再加 $2\sqrt X\sum_\chi\epsilon_\chi/(Q-1)$ 即可。

## 2. 先投影到剩余类，再使用函数方程

Alpöge–Claude 归档中的 `mellin_ramp` 把斜坡权 $(y/n-1)_+$ 变换为 $n^{-s}/(s(s-1))$。本文保留该权中的 $\Lambda(n)\chi(n)$ 系数。用 $U_v$ 表示区间 $[-v,v]$ 上的均匀概率密度，取

$$
 g=\mathbf1_{[-H,L+H]}*U_{H/m}^{*m}.
$$

它支撑于 $[-2H,L+2H]$，在 $[0,L]$ 上恒等于 $1$，且 $\int_{\mathbb R} g(u)e^{itu}du=\Phi(t)$。令 $F(y)=\sqrt y\,g(\log y)$。对每个 $n\ge1$ 分部积分给出 $\int_1^\infty(y/n-1)_+ F''(y)\,dy=F(n)/n=n^{-1/2}g(\log n)$。这说明本测试可在斜坡的尺度变量上施加一个明确的权。下文直接用带特征的显式公式评价这个测试。纵坐标表述使用原归档的 `LFunction_eq_zero_re`。[原斜坡引理](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/DedekindLFunctions/LandauTheorems/MellinRamp.lean) · [原归档的 Dirichlet 定理](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Main.lean) · [五个定理的独立核验](https://math.mit.edu/~drew/grh/)

原归档在 `Level.Gm` 中对各特征的斜坡和取实部，再按 `Level.wt` 的实权相加，从这个标量的单侧增长界排除零点。具体整数 $a$ 的剩余类投影需要随 $a$ 改变的复系数 $\overline{\chi(a)}$。因此我们在已经取得 GRH 后回到各个带标签的显式公式，先施加这个系数，再对特征求和。一个导子由此服务于所有 $a$。[原稿的权、标量和及增长条件](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Formalism/SharpRestriction/RamificationFrobenius/TotientGaloisClasses/ExplicitFormulaErrorOfHostConditions.lean#L30-L40) · [文库中对循环层平移的分析](../Elliott–Halberstam/theta-chains-and-translated-cyclic-tests.md)

对单位剩余类 $a$ 记 $\Delta_a(n)=\mathbf1_{n\equiv a\bmod Q}-\mathbf1_{n\equiv-a\bmod Q}$。奇特征的归一化正交关系为

$$
 \frac2{Q-1}\sum_{\chi(-1)=-1}\overline{\chi(a)}\chi(n)
   =\Delta_a(n).
$$

用 $\overline{\chi(n)}$ 代替 $\chi(n)$ 得 $\Delta_{a^{-1}}(n)$。取 $n=1$ 得公共项的系数为零，因为 $2\le a\le X<Q/2$。模素数 $Q$ 的奇特征都是本原非主特征。于是完整零点和满足

$$
 \frac2{Q-1}\sum_{\chi(-1)=-1}\overline{\chi(a)}
       \sum_{\gamma_\chi}\Phi(\gamma_\chi)
 =-\sum_{n\ge2}\frac{\Lambda(n)}{\sqrt n}
 \bigl[\Delta_a(n)g(\log n)+\Delta_{a^{-1}}(n)g(-\log n)\bigr]. \tag{4}
$$

为核对符号与倒数类，取 $G(s)=\int g(u)e^{(s-1/2)u}du$，在 $\Re s=3/2$ 对 $G(s)L'(s,\chi)/L(s,\chi)$ 作 Mellin 反演。右线给 $-\sum_n\Lambda(n)\chi(n)g(\log n)/\sqrt n$。左移至 $\Re s=-1/2$ 时，奇特征的平凡零点都在两线之外。函数方程把左线上的对数导数换成 $-L'(1-s,\bar\chi)/L(1-s,\bar\chi)$ 加上只依赖 $Q$ 及奇偶性的项。后者在特征投影中消失。右线减左线，就得到 (4) 中带相同负号的两项。$G$ 在这条闭带沿竖向衰减为 $O_{X,Q,m}((1+|\Im s|)^{-m-1})$。对固定 $Q$ 可在趋于无穷的高度同时避开这有限族的零点，标准对数导数展开和单位高度的零点计数使横边积分趋零。由 GRH，$G(\rho)=\Phi(\gamma_\chi)$。共轭特征的零点纵坐标互为相反数，故 (1) 为实数。[Dirichlet 函数方程](https://dlmf.nist.gov/25.15.E5) · [Bennett、Martin、O’Bryant、Rechnitzer 的零点计数](https://arxiv.org/abs/2005.02989)

式 (4) 的第一项要求 $n\le Xe^{2H}=Q-X-1$。在这个范围内，$\pm a$ 两类的唯一正代表元是 $a$ 本身，且 $g(\log a)=1$。第二项要求 $n\le e^{2H}$。此时 $2\le an\le Q-X-1<Q-1$，从而 $an$ 不可能模 $Q$ 余 $\pm1$。所以右侧恰为 $-\Lambda(a)/\sqrt a$。排除倒数类这一步使整个 $[2,X]$ 可以共用一个测试函数。

## 3. 截高为何可以固定为 3

记 $N_\chi(v)$ 为 $|\gamma_\chi|\le v$ 的零点数。对 $v\ge0$ 可以统一使用

$$
 N_\chi(v)\le10(1+v)\log(2Q(2+v)).
$$

对 $v\ge1$ 这是 Bennett 等人定理 1.1 的宽松形式。对 $0\le v<1$ 用 $N_\chi(v)\le N_\chi(1)$，直接将该定理在 $1$ 处的主项及误差代入，也小于 $10\log(4Q)$。[定理 1.1 的原文](https://arxiv.org/pdf/2005.02989#page=2)

若 $2^jT<|\gamma|\le2^{j+1}T$，$T\ge em/H$，则 $|\Phi(\gamma)|\le W e^{-m}2^{-jm}$。相应零点数不超过 $20(1+T)2^j(j+1)\log(4Q(2+T))$。因为 $m\ge3$，有 $\sum_{j\ge0}(j+1)2^{-j(m-1)}<2$。在 (1) 中特征平均的绝对权重之和为 $1$（去掉外面的 $\sqrt a$），乘回 $\sqrt a\le\sqrt X$ 就得 (2)。

**推论 2.** 固定 $\eta>0$，在 $[2X^{1+\eta},4X^{1+\eta}]$ 中选素数 $Q$。对任一固定的 $T>e/\eta$，当 $X$ 充分大时，(3) 同时恢复所有 $a\le X$ 的素性（从 $a=2$ 起）。特别地，$\eta=1$ 可取 **$T=3$**。

事实上该范围内 $H\ge(\eta/2)L$ 对充分大的 $X$ 成立，$W=O_\eta(L)$。取 $m=\lfloor TH/e\rfloor$，则

$$
 B\ll_{\eta,T} L^2 X^{1/2-T\eta/(2e)}=o(1).
$$

素数 $Q$ 的存在性可由 Bertrand 定理得到（端点取整也可直接用素数定理）。若改在 $[2X e^{L^2},4X e^{L^2}]$ 中取素数，则 $H\ge L^2/2$、$W=O(L^2)$。取 $T=3/L$ 和同样的 $m$，此时 $B\ll L^4 X^{1/2-3/(2e)}=o(1)$。这给出了趋于零的**绝对高度**。以平均零点间距为单位衡量，$T\log Q$ 在这两个例子中都与 $\log X$ 同阶。

## 4. 两种计数与一个有限证书

通过 (3) 定义 $b(a)=\mathbf1_{D_T(a)>(3/4)\log a}$。对 $4\le n\le X$ 的偶数，以及 $2\le x\le X$，有精确整数公式

$$
 G(n)=\sum_{2\le a\le\lfloor n/2\rfloor}b(a)b(n-a),
 \qquad
 \pi_2^{\rm end}(x)=\sum_{2\le a\le x-2}b(a)b(a+2). \tag{5}
$$

$G$ 计无序表示，$\pi_2^{\rm end}$ 要求较大的素数端点不超过 $x$。相同的 $S_\chi(T)$ 用于 (5) 的每一个 $n,x$。要由这些公式证明两个完整猜想，仍须证明每个偶数的 $G(n)>0$ 和 $\pi_2^{\rm end}(x)$ 的无界性。另文给出在每个对数平方短区间内实际可证的表示下界。[带差二条件的短区间下界](../Goldbach_TwinPrimes/chen-summands-short-intervals.md)

有限例子取 $X=13,\ Q=53,\ m=17,\ T=88$。此时 $H=(\log3)/2$，满足定理的截高条件，(2) 给出 $B<0.019189<(\log2)/8$。随文程序从 Dirichlet $L$ 函数计算零点，用 Arb 实球包围全部 $26$ 个奇特征的 $4094$ 个带标签零点，再把这些球代入 (1)。对 $2\le a\le13$，直接比较得到 $\max_a|D_{88}(a)-\Lambda(a)|<1.14\times10^{-8}$。$b(a)=1$ 的整数为 $2,3,5,7,11,13$。$n=4,6,8,10,12$ 的无序表示数依次为 $1,1,1,2,1$，截至 $13$ 共有三对孪生素数。高度 $3$ 是推论 2 的渐近结论，本例以 (2) 单独核验参数。[复算程序（python-flint、SciPy）](shared-low-zeros-prime-pairs/verify.py)

程序对每个特征的 Hardy $Z$ 函数逐段找出变号区间。为认证零点无遗漏，记 $\alpha_\chi(T)$ 为从 $2+iT$ 沿水平线连续延拓得到的 $\arg L(1/2+iT,\chi)$，起点取主辐角。设 $\theta(T)=\Im\log\Gamma(3/4+iT/2)+(T/2)\log(Q/\pi)$，其中对数伽马在正实轴上为实数。由辐角原理，

$$
 N_\chi(T)=\frac{2\theta(T)+\alpha_\chi(T)+\alpha_{\bar\chi}(T)}{\pi}.
$$

具体取完成函数从 $1/2-iT$ 经 $2-iT,2+iT$ 到 $1/2+iT$ 的路径，再以函数方程关于临界线反射闭合。全周辐角变化是右侧路径的两倍。$\Re s=2$ 上 $|L(s,\chi)-1|<1$，两条水平边以区间算术细分，直到每段的函数像相对于其起点值落在右半平面，即可相加连续辐角的增量。将所得整数与 $\chi,\bar\chi$ 的正高度变号区间总数比较，并检验中心点与边界。两者相等时，每个分离区间恰含一个单零点。程序只用浮点求根选择短区间，端点变号、总数、最终阈值均用球算术核对。
