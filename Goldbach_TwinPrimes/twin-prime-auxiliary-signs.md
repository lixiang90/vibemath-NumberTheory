---
title: "固定差 2 的四个辅助二次符号与共同好窗口"
date: "2026-11-24"
abstract: "在两阶段选素数的原概率及例外指标界下，对全部子区间同时控制筛选误差。证明真实孪生素数质量的 1/16 转移和斜坡因子 1/3，明确固定正下界仍需的算术输入。"
---

# 固定差 2 的四个辅助二次符号与共同好窗口

对 ζ 目标，Alpöge–Claude 的构造依次选择两个有关联的辅助素数 $q,d$。本文问在两个真实素数 $p,p+2$ 处同时取四个指定二次符号，会保留多少质量。以原来的候选素数为概率空间，证明在全部子区间上一致的绝对二阶误差界，并将其限制到共享一扇完整上同调窗口的好对。四符号保留原质量的 $1/16$，归一化误差趋零。借助一个本原偶特征的有限反演，斜坡二次统计还带有因子 $1/3$。

论证使用以下三项外部输入，包括两阶段好对的比例、各步至多一个例外平方自由指标的一阶特征估计、以及共同好窗所容许的例外扭转指标个数。下面给出各自的数值范围与出处。(4)、(9)、(10)、(12)、(13) 是在这些输入之上的推导。

## 一　保留原来的两个选择步骤

固定 ζ 目标及原构造容许的 $0<\epsilon<1/52$。记基准尺度为 $P$，外端为 $H=P^J$，$J$ 只依赖固定数据。设 $z=P^{\epsilon/5}$，$z'=z(\log P)^2$。在以下两个有限集合上依次取均匀概率。

1. 在 $e^{z/2}<q\le e^z$、$q\equiv3\pmod8$ 的素数中抽取 $q$。
2. 给定 $q$，在 $e^{z'/2}<d\le e^{z'}$ 的素数中抽取 $d$，要求 $d\equiv7\pmod {16}$，且 $d\equiv-1$ 模 $q$ 及模每个奇素数 $\ell\le z$。这些条件给出模 $m_q=16q\prod_{3\le\ell\le z,\ \ell\ {\rm prime}}\ell$ 的一个可逆类。

记这个两步概率为 $\mathbb P_P$。各 $q$ 所带的 $d$ 的个数可以不同。对充分大的 $P$，上述每个候选集合均非空。Lean 的 $\mathrm{pair\_share}$ 给出更具体的余量。对任意固定 $0<\lambda<1$，至少 $1-\lambda$ 比例的 $q$ 通过第一步，在每个这样的 $q$ 下，至少 $1-\lambda$ 比例的 $d$ 通过第二步。取基模参数 $8$、指数 $\tau=\epsilon/2$、上限 $H$，原证明把这些条件逐对转成 $\mathrm{PairGood}(1,\epsilon,P,H,q,d)$，包括中心化特征和、在 $1$ 的对数导数及规定的无零矩形。记满足这个结论的事件为 $\mathcal G_P$。取 $\lambda=1/2$ 即有 $\mathbb P_P(\mathcal G_P)\ge1/4$。[候选集合的定义](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/MeanSquareBounds/PrimeSetTwoBounds.lean) · [两阶段保留比例](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/MeanSquareBounds/PairShare.lean) · [从两步条件到 PairGood](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/MeanSquareBounds/pairEx.lean)

符号函数 $I_{q,d}$ 来自 ζ 目标的可见类。其两坐标在未分歧素数 $p$ 处取 $((p/q),(p/d))$，原构造指定的坐标分别是 $+1,-1$。[特征坐标](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/CharacterBoxes/FieldChi.lean) · [原构造的核](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/CharacterBoxes/KerChiAny.lean)

原选择接着为同一对 $q,d$ 找到 $P\le W,\ W^2\le H$，在整个 $[W,W^2]$ 上取得上同调估计。ζ 情形的物理缩放因子为 $1$，故可取 $Y=W$，物理窗口为 $[Y,Y^{3/2}]$。固定窗内的素数求和用

$$
 Z_Y=Y^{3/2}/2,\qquad L_Y=\sum_{Y\le X\le Z_Y,\ X\in\mathbb N}X^{-1}
 =\tfrac12\log Y-\log2+O(Y^{-1}).
$$

对一对真素数 $p,p+2$，设 $w(p)=\log p\log(p+2)/p$，$\chi_r(n)=(n/r)$。记 $I_{q,d}(p)$ 为 $\chi_q(p)=\chi_q(p+2)=+1$ 且 $\chi_d(p)=\chi_d(p+2)=-1$ 的指示函数。定义

$$
 M(Y)=\frac1{L_Y}\sum_{\substack{2Y<p\le Z_Y\\p,\ p+2\ {\rm prime}}}w(p),
 \qquad
 J_{q,d}(Y)=\frac1{L_Y}\sum_{\substack{2Y<p\le Z_Y\\p,\ p+2\ {\rm prime}}}
 w(p)I_{q,d}(p). \tag{1}
$$

以下始终把 $p,p+2$ 固定为实际的素数。[原始窗口的量词](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/CharacterBoxes/ExistsPageN.lean) · [先选一对素数再取整段窗口的 WinPairs](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/CharacterBoxes/WinPairs.lean)

## 二　一个例外二次特征至多需要删一对

辅助素数选择中有一个可用于二点统计的原始估计。取第一步的全体候选，或固定 $q$ 后取第二步的全体候选，以 $r$ 表示正在抽取的素数。相应的模数 $m$ 为 $8$ 或 $m_q$。对任意固定 $D>0$，除至多一个例外的平方自由整数 $s_*$ 外，有

$$
 \left|\mathbb E_r(r/s)\right|\le P^{-D}
 \quad
 (1<s\le P^{800J},\ s\text{ 为奇平方自由数},\ (s,m)=1). \tag{2}
$$

这里 $(r/s)$ 是 Jacobi 符号。Rubinstein 的 ζ 讲义第三部分 §III.1 步骤 1 的式 (III.9) 取其基准为 $P$、$A=J$，即给此范围。其证明对模 $ms$ 的各类计素数，字符正交消去主项，Page 唯一性留下至多一个例外。两族均满足原可容许性条件，特别是 $\log m_q\le4z$，$m_q$ 的奇素因子或至多为 $z$，或为远大于任意固定次幂 $P^k$ 的 $q$。第二步的 $s_*$ 允许依赖 $q$。下文只用 $D=3$ 和 $s\le(H+2)^4$。归档中基础引理采用固定截断和 $-3$ 次幂，任意固定 $J$ 的截断范围用上述 (III.9)。[Rubinstein，§III.1 步骤 1、6、7](https://www.math.uwaterloo.ca/~mrubinst/alpoge_claude_grh/zeta-paper_v2.pdf#page=63) · [归档中的特征和 twist_bound](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Formalism/GL3MeanSquares/GenusMeanAsymptotics/JacobiCharacterSums/MSqTwist.lean)

取 $\mathcal E_P=\{p\in[P,H]:p,\ p+2\text{ 都是素数}\}$。当 $P>5$，不同的这类素数对没有公共端点，否则存在三个连续的奇数均为大于 $3$ 的素数。所有端点都大于 $z$、小于候选辅助素数，与 $m$ 互素。非空端点子集的乘积是大于 $1$ 的奇平方自由数。二次互反律把其在 $r$ 处的符号变成 (2) 中的符号乘一个固定正负号，因为每族中的 $r\bmod4$ 固定。

若某个端点乘积可能等于 $s_*$，固定 $s_*$ 的一个素因子，暂时删除含该素因子的那一对孪生素数。它只属于一对。此后用余下端点写出的任何乘积均不能等于 $s_*$。若原本就无法写出 $s_*$，一对也无需删。对任意子区间，恢复这一对至多花费

$$
 w_{\max}\ll_J(\log P)^2/P. \tag{3}
$$

这个论证每个候选族单独使用，因而容许第二步随 $q$ 改变删除对象。它给出真实素数对总质量的估计，个别预先指定的素数对仍可能涉及例外。

## 三　对所有子区间同时估计四个符号

**命题。** 令

$$
 K_P(q,d)=\sup_{I\subset[P,H]\,\text{为区间}}
 \left|\sum_{p\in\mathcal E_P\cap I}w(p)\bigl(I_{q,d}(p)-1/16\bigr)\right|.
$$

原两步概率满足

$$
 \mathbb E_P K_P^2\ll_{\epsilon,J}\frac{(\log P)^6}{P},
 \qquad
 \mathbb P_P(K_P>P^{-1/4})
 \ll_{\epsilon,J}\frac{(\log P)^6}{\sqrt P}. \tag{4}
$$

第一个界是**未除以 $L_Y$ 的绝对误差**的二阶矩。它同时覆盖有很多、很少乃至零对孪生素数的区间。

为证 (4)，在单个辅助素数族内取 $\sigma\in\{+1,-1\}$，写

$$
 A_r^\sigma(p)=\frac{(1+\sigma\chi_r(p))(1+\sigma\chi_r(p+2))}{4},
 \qquad
 F_r^\sigma(p)=A_r^\sigma(p)-\tfrac14
 =\frac{\sigma\chi_r(p)+\sigma\chi_r(p+2)+\chi_r(p)\chi_r(p+2)}4.
$$

删除第二节的例外对后，若 $p\ne t$，$F_r^\sigma(p)F_r^\sigma(t)$ 展开成九项。各项涉及二至四个互异的素数端点，每项均可用 (2)，故乘积的期望绝对值 $\ll P^{-3}$。对角项用 $|F_r^\sigma|\le1$。取在该步抽样前已固定的实权 $|v_p|\le w(p)$，对任一连续块 $D$，在余下素数对上得到

$$
 \mathbb E_r\left|\sum_{p\in D}v_pF_r^\sigma(p)\right|^2
 \ll \sum_{p\in D}w(p)^2+P^{-3}\left(\sum_{p\in D}w(p)\right)^2. \tag{5}
$$

把这些素数对按 $p$ 递增排列，补零至长度为二的幂。任一前缀可写成至多 $O_J(\log P)$ 个二进块的并，每层至多取一块。用 Cauchy–Schwarz 后对各层的全部块求和，每层的块互不相交。再把区间写成两个前缀之差，恢复例外对，得到

$$
 \mathbb E_r\sup_I\left|\sum_{p\in\mathcal E_P\cap I}v_pF_r^\sigma(p)\right|^2
 \ll_J w_{\max}^2+(\log P)^2
 \left(\sum_{p\in\mathcal E_P}w(p)^2+
 P^{-3}\left(\sum_{p\in\mathcal E_P}w(p)\right)^2\right).
$$

空集和单元素情形直接由 (3) 控制。把求和放宽到全部整数即得 $\sum w^2\ll_J(\log P)^4/P$、$\sum w\ll_J(\log P)^3$。第一步取 $v_p=w(p)$、$\sigma=+1$。给定 $q$ 后，第二步取 $v_p=w(p)A_q^+(p)$、$\sigma=-1$。恒等式

$$
 I_{q,d}(p)-\tfrac1{16}
 =A_q^+(p)\bigl(A_d^-(p)-\tfrac14\bigr)
  +\tfrac14\bigl(A_q^+(p)-\tfrac14\bigr) \tag{6}
$$

给出 (4) 的第一式，第二式由 Markov 不等式推出。这里用的是固定 $q$ **以后**对 $d$ 的估计，完整保留了 $q\mid d+1$。一个辅助素数在两个端点同时取指定符号，给出系数 $1/4$，第二次条件平均再给一个 $1/4$。式 (4) 直接计入实际出现的孪生素数。

还可以逐项辨认剩下的算术量。展开 $16I_{q,d}$，令 $S,T\subset\{0,2\}$ 标记一个单项式的端点。除常数项外，下面的 $15$ 个和以概率 $1-O_{\epsilon,J}((\log P)^6/\sqrt P)$ **同时**满足

$$
 \max_{(S,T)\ne(\varnothing,\varnothing)}\
 \sup_{I\subset[P,H]\,\text{为区间}}
 \left|\sum_{p\in\mathcal E_P\cap I} w(p)
     \prod_{h\in S}\chi_q(p+h)\prod_{h\in T}\chi_d(p+h)\right|
 \le P^{-1/4}.
$$

若 $T$ 非空，固定 $q$，把 $q$ 的符号都放进 $v_p$，对 $d$ 重做 (5)。两个不同素数对给出非空的平方自由端点乘积，仍可用 (2) 和一次例外删除。若 $T$ 为空，则用 $q$ 这一族。对 $15$ 个失败概率求和即可。这样展开后，唯一需要新的下界的项是常数特征所乘的真实素数对和 $\sum w(p)$。

## 四　同一个原好窗口上的平均

因为 $\mathbb P_P(\mathcal G_P)\ge1/4$，(4) 的失败概率在原解析好对中也趋零。在这些好对里可以选择 $q,d$，使下式对**该对的每一个原好窗口**都成立，特别可先选对，再执行原来的上同调选窗步骤。

$$
 \left|J_{q,d}(Y)-\frac{M(Y)}{16}\right|
 \le \frac{P^{-1/4}}{L_Y}. \tag{7}
$$

若要多对 $q,d$ 在完全相同的 $Y$ 上取平均，还可保留原证书覆盖的整段 $[Y,Y^2]$。回到上同调选窗引理。对固定 ζ 数据，存在与 $P,q,d$ 无关的整数 $R\ge0$ 和常数 $0<\eta\le C_t$，此处 $C_t=401$。令 $B=C_t/\eta$。原证明可取 $J=2^{R+2}B^R$。给定外端与素数对后，有集合 $\mathcal X=\mathcal X(P,q,d)$，至多含 $R$ 个例外扭转指标。在基准范围内，只要每个 $h\in\mathcal X$ 都满足 $h\le x^{C_t}\Rightarrow h\le x^\eta$，证明中该点所需的塔估计就成立。其余上同调估计在这个条件下逐点拼合。这里的扭转指标与第二节解析特征的例外 $s_*$ 来自不同的定理。[例外集合及窗口常数 ThmBabOf](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/AutomorphicVarieties/LayeredIsomorphisms/ThmBabOf.lean) · [从任意合格窗口得到联合估计](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/AutomorphicVarieties/PrimeSetBoundDataOfPrimitiveChar.lean) · [还原到原始上同调窗口](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/AutomorphicVarieties/CMFieldExponentialBounds/AdmissibleOfPageN.lean)

在抽取任何辅助素数之前列出 $R+1$ 个起点

$$
 Y_i=P^{(2B)^i},\qquad 0\le i\le R. \tag{8}
$$

每个 $[Y_i,Y_i^2]$ 都在 $[P,H]$ 中。若某个 $h$ 在这段的某处违反上述蕴含，就有 $Y_i^\eta<h\le Y_i^{2C_t}$。由于 $Y_{i+1}^\eta=Y_i^{2C_t}$，这些左开右闭的指标区间两两不交。$R$ 个指标至多影响其中 $R$ 段。故对每个解析好对，至少有一个 (8) 的完整原始窗口合格。这里保留了平方右端，物理子窗 $[Y_i,Y_i^{3/2}]$ 保留在该证书范围内。

设 $\mathcal G_{P,i}$ 为同时满足 $\mathcal G_P$ 和第 $i$ 段完整上同调条件的候选对，$a_i=\mathbb P_P(\mathcal G_{P,i})$。同一辅助对可以属于多个 $\mathcal G_{P,i}$。由覆盖性，$\max_i a_i\ge1/[4(R+1)]$。对任意 $a_i>0$，由 (4)、$L_{Y_i}\gg\log P$ 及条件 Cauchy–Schwarz，

$$
 \mathbb E_P\!\left[
 \left|J_{q,d}(Y_i)-M(Y_i)/16\right|\ \middle|\ \mathcal G_{P,i}\right]
 \ll_{\epsilon,J}\frac{(\log P)^2}{\sqrt{P\,a_i}}. \tag{9}
$$

因此每个充分大的 $P$ 都有一个共同的 $Y_i$，在该窗的原好对上取平均，筛前质量的 $1/16$ 与筛后质量只差 $o(1)$。若另行指定某段，只要 $P a_i/(\log P)^4\to\infty$，同样可用 (9)。仅知道某个指定起点适用于一对辅助素数时，须另查其 $a_i$。

## 五　还差哪一个严格下界

为把窗口选择也写成可检验的公式，在 $\mathcal G_P$ 中对每个好对选 (8) 中第一个合格下标 $i(q,d)$，令 $\nu_i(P)=\mathbb P_P(i(q,d)=i\mid\mathcal G_P)$。它们是和为 $1$ 的权。由 (4) 的一致估计，直接在整个 $\mathcal G_P$ 上条件化即得

$$
 \mathbb E_P\!\left[J_{q,d}(Y_{i(q,d)})\mid\mathcal G_P\right]
 =\frac1{16}\sum_{i=0}^R\nu_i(P)\,M(Y_i)
   +O_{\epsilon,J}\!\left(\frac{(\log P)^2}{\sqrt P}\right). \tag{10}
$$

所以对这种平均，固定正下界准确地归结到

$$
 b=\limsup_{P\to\infty}\sum_{i=0}^R\nu_i(P)M(Y_i)>0. \tag{11}
$$

式 (11) 中的权允许依赖辅助对的选窗过程。若 $b>0$，沿无穷多个基准可选择好对，使所选窗口的 $\limsup J_{q,d}\ge b/16$。对所有选择都有 $0\le J_{q,d}(Y)\le M(Y)$。这里的估计能独立应用于任何预先给定的真孪生素数子集，包括空集。筛选前的正质量需要额外的算术下界。[尺度平均与共同好窗口](twin-prime-scale-averages.md)

## 六　从筛后质量到两个截点的二次统计

给出一个可独立验证的特征反演。置 $h=\exp(-\sqrt{\log Y}/2)$，$x=(1+h)X$，$Y\le X\le Z_Y$ 取整数，并记 $I_X=\{n:X<n<n+2\le x\}$。取任意二次幂 $Q$，使 $8\mid Q/2$ 且 $Q/2>2(x+2)^2$。以 $\mathbb E_Q$ 记导子恰为 $Q$ 的本原偶特征平均。对支撑在 $1\le m\le x+1$ 的实系数 $a_m$，令
$$
D_X(\chi)=\Re\sum_{\substack{n\in I_X\\ n\ {\rm odd}}}
             \overline{\chi((n+2)n^{-1})},\qquad
P_a(\chi)=\sum_{m\ {\rm odd}}a_m\chi(m).
$$

这里逆元模 $Q$ 取得。当 $X$ 充分大、$x+1<2X$ 时，有
$$
\mathbb E_Q D_X(\chi)|P_a(\chi)|^2
 =\sum_{\substack{n\in I_X\\n\ {\rm odd}}} a_n a_{n+2}. \tag{12}
$$

先对模 $Q$ 的所有偶特征取平均。正交性要求两个奇数之比模 $Q$ 等于 $\pm(n+2)/n$ 或其逆元。乘开分母，所有正整数乘积之和小于 $Q/2$，所以负号无解，正号是通常的有理数等式。由于 $(n,n+2)=1$，正整数解只能是 $kn,k(n+2)$。支撑强制 $k=1$。取实部所给的两个次序各带 $1/2$。同样的论证对模 $Q/2$ 成立。全体偶特征的半数从模 $Q/2$ 提升，余下半数正是本原偶特征，故 (12) 成立。

在 (12) 取 $a_m=\vartheta(m)V(m)\kappa_X(m)$，其中 $\vartheta(m)=\log m$ 当 $m$ 为素数，否则取零，$V(m)=\mathbf1_{(m/q)=+1,\,(m/d)=-1}$，$\kappa_X(m)=h^{-1}(x/m-1)_+$。也可以把 $m\le X$ 的系数换成同一斜坡在两个截点的差，或任意实权，(12) 仍然成立。用 $1/(L_YhX^2)$ 对左边累加，把所得实数记为 $\gamma(Y;q,d)$。则对每一对辅助素数一致有
$$
 \gamma(Y;q,d)=c_h J_{q,d}(Y)+o(1),\qquad
 c_h=h^{-3}\int_1^{1+h}\left((1+h)/v-1\right)^2dv
 =\tfrac13+O(h). \tag{13}
$$

事实上，交换 $X$ 与素数对起点 $p$ 的求和，内段 $2Y<p\le Z_Y$ 所得权为 $c_h/p+O(1/(hp^2))$，还需乘两端的对数及四符号指示。把第二个斜坡换成第一个，点态误差为 $O(1/(hp))$；在 $p/(1+h)\le X<p$ 内对 $\kappa_X(p)^2/X^2$ 用单调性将和换成积分，再以 $v=p/X$ 换元即得。内段误差除以 $L_Y$ 后趋于零，因为 $h^{-1}=Y^{o(1)}$。两端 $[Y,2Y]$ 与 $[Z_Y,2Z_Y]$ 对归一化质量均只贡献 $O(1/L_Y)$，此处用孪生素数上界筛 $\sum_{u<p\le2u,\;p+2\ {\rm prime}}w(p)\ll1$。[上界筛及其常数](https://arxiv.org/abs/2109.02851)

沿一列所选的 ζ 好窗口定义 $\Gamma_1=\limsup\gamma(Y;q,d)$。由 (13)，$\Gamma_1=\frac13\limsup J_{q,d}(Y)$。因此 (11) 的 $b>0$ 将给出可选的一列窗口满足 $\Gamma_1\ge b/48$。已证明的 $1/16$ 与 $1/3$ 分别来自两次条件平均和斜坡积分。$b>0$ 仍是开放的固定差二点估计。
