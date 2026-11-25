---
title: "从 theta 提升链的重数到分剩余类的素数计数"
date: "2026-11-25"
abstract: "Alpöge–Claude 循环层的标量只记录每一特征阶上的总和。本文计算其核，在原稿的 Bernoulli 选点中移动加重的陪集，同时约束一个指定的复特征和。保留原有坏 Frobenius 轨道预算后，原陪集误差、保留的均方估计、针对给定算术量选点，分别使近似误差在指数小于 3/5、3/4、1 时达到 EH 精度。真实提升链上的带符号主估计仍是开放步骤。"
---

# 从 theta 提升链的重数到分剩余类的素数计数

Alpöge–Claude 的证明在哪一步失去了区分剩余类的能力？答案可以定位到一个有限群上的傅里叶变换。证明给一个陪集比其他陪集更大的权重，随后把变换得到的系数替换为同一个实数。移动这个被加权的陪集，系数就会带上 EH 所需的特征相位。原稿输出的计数精度可以使近似误差在 $\theta<3/5$ 时得到控制。回到选择顶端点集的步骤，保留已证明的均方界，可以达到 $\theta<3/4$。在选择中再同时约束一个指定的算术量，可把**抽样误差**控制在每个固定的 $\theta<1$。平移后算术量本身的平均上界仍待证明。

## 1. 重数由哪一步产生

在斯捷潘诺夫的曲线方法中，辅助函数在有理点的高阶消失使每个点对零因子的次数作多次贡献。Alpöge 将数域中的对应机制描述为，在正定酉群之间连续作 $m$ 次 theta 提升。源码中的单步算子是有限加权核，链算子是这些算子的复合。其复数内积和伴随在这一层仍有定义。提升的整性和局部可除性给出格指数的估计，Rallis 内积公式参与范数计算，最后的除子不等式按层累加局部项、对数余体积和误差。这里的“重数”是算术除子估计中的累积增益，与单个 Dirichlet 零点的解析重数是两个量。[作者的说明，第 3–4 页](https://alpo.ge/remark.pdf) · [核、伴随与链复合](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Formalism/ThetaLiftIndependence/Statement.lean#L10-L45) · [提升塔与内积公式的衔接](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/ThetaCorrespondence/ThetaPairLifts/LocalTypeDictionary/ThetaTowerCompatibility.lean#L12-L16) · [除子不等式中的逐层求和](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Formalism/DivisorInequality/StatementObj.lean#L59-L75)

“链长约为 $\sqrt{x}$”有精确的版本。归档设 $\delta_1=x^{-\epsilon/8}$、$\tau=\epsilon/100$，要求
$$
 x^{-\tau}\sqrt{x\delta_1}\le m\le x^\tau\sqrt{x\delta_1}.
$$
同一组参数还要求 $k\ge m\sqrt{x}\log x$，并推出 $m^4\le n$。改变链长时需要同时重新核对这些秩和误差预算。即使把标量不等式的误差完全消除，其对特征的观测方式也仍由下一节的权重决定。[尺度定义](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Formalism/LocalRingCardinality/DataParams%E2%82%82.lean#L27-L45) · [`SizeM` 与 `SizeKRel`](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Formalism/SharpRestriction/ElevenLe.lean#L11-L43) · [链长四次方的界](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Formalism/LocalRingCardinality/FourthPowerRankLeNF.lean#L8-L19)

## 2. 原有循环层能区分哪些零点

对本原非主特征 $\chi$，原稿的统计量为
$$
 g_\chi(t)=\sum_{n\ge1}\Lambda(n)\chi(n)(t/n-1)_+
             +t\,\frac{L'}{L}(1,\chi).
$$
循环子群 $H$ 的阶记为 $k$。原稿的 `Level.ofSub` 把 `psK` 设为 `none`，故 `Level.Gm` 对所有非主特征使用 $1/(3k+1)$ 的权重。减去共同的主项 $G_\zeta(t)$，乘回 $3k+1$，所得为 $\sum_{\chi\in H,\chi\ne1}\Re g_{\chi^*}(t)$。符号 $\chi^*$ 表示诱导该特征的本原特征。$H$ 关于复共轭封闭，所以这个和也等于相应复数的和。原链在所选窗口提供对 `Gm` 的单侧上界。[解析量的定义](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Formalism/SharpRestriction/RamificationFrobenius/TotientGaloisClasses/Analytic.lean#L11-L26) · [`Gm` 的权重](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Formalism/SharpRestriction/RamificationFrobenius/TotientGaloisClasses/ExplicitFormulaErrorOfHostConditions.lean#L30-L41) · [循环层的选择](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Formalism/SharpRestriction/RamificationFrobenius/TotientGaloisClasses/ExplicitFormulaErrorOfHostConditions.lean#L105-L123) · [窗口上的不等式](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Formalism/SharpRestriction/RamificationFrobenius/TotientGaloisClasses/TopWin.lean#L23-L32)

**命题 1。** 对素数模数 $p$，令 $g=p-1$。在 $v_1=0$、$v_{\bar\chi}=\overline{v_\chi}$ 的实向量空间上，所有循环层的和作为线性观测，其秩为 $\tau(g)-1$，核维数为 $g-\tau(g)$。

**证明。** 每个 $d\mid g$ 对应一个阶 $d$ 的子群。该层记录所有阶整除 $d$ 的系数之和。对除数作 Möbius 反演恢复每个给定阶上的总和。不同阶的系数支撑不交，故 $d>1$ 时这些总和线性无关。起始空间的实维数为 $g-1$，得证。例如模 $47$ 时，$45$ 个实自由度中只能观测到 $3$ 个。

这一点在真实零点上有直接解释。对 $\Re s>1$，逐项积分得
$$
 \int_1^\infty g_\chi(t)t^{-s-1}\,dt
 =-\frac{(L'/L)(s,\chi)}{s(s-1)}
   +\frac{(L'/L)(1,\chi)}{s-1}.
$$
在非平凡零点 $\rho$ 处取亚纯延拓的留数，再乘 $-\rho(\rho-1)$，正好得到重数 $m_\chi(\rho)$。在临界线 $\rho=1/2+i\gamma$ 上，各零点对未平移循环层的留数贡献同号，权为 $1/(1/4+\gamma^2)$。如果精确知道整个层函数，就可以确定每个特征阶在该高度的总重数。若给出全部平移，则对 $U_q=(\mathbb Z/q\mathbb Z)^\times$、$g=\varphi(q)$ 定义
$$
 F_q(b,t)=\frac1g\sum_{\chi\bmod q,\,\chi\ne1}
                 \overline{\chi(b)}\,g_{\chi^*}(t),\qquad b\in U_q.
 \tag{1}
$$
它在实端点取实值。其 Mellin 留数乘 $-g\rho(\rho-1)$ 后等于 $\sum_{\chi\ne1}\overline{\chi(b)}m_{\chi^*}(\rho)$。对 $b$ 反演就逐特征恢复该高度的重数。这是精确数据的可辨识性结论，窗口上的单侧上界没有给出这些留数。

**一个保留四阶峰值的例子。** 上述核里甚至有系数全实、对取逆不变的数据。固定 $1/2<\theta<1$，$D\asymp x^\theta$，取素数 $D<p\le2D$、$p\equiv5\bmod6$，以原根的指数 $j\bmod g$ 表示剩余类。取二次幂 $D^2/x\le m<2D^2/x$，令 $I$ 从 $\lfloor g/8\rfloor$ 开始、含 $m$ 个连续指数。把归一化的 Rudin–Shapiro 符号 $c s_j/\sqrt m$ 放在 $I$ 上（$c=1/16$），在 $-I$ 上放其镜像，在 $3I,-3I$ 放相应的相反数。因为 $m=o(g)$ 且 $3$ 是模 $g$ 的单位，四份支撑两两不交。所得序列 $Z$ 满足
$$
 \|Z\|_\infty=c/\sqrt m,\qquad
 \frac1g\sum_j|Z(j)|^4=\frac{4c^4}{gm}.
$$
以 $B(k)$ 表示第一份的未归一化傅里叶系数，则总系数为 $v_k=2\Re(B(k)-B(3k))$。由 Rudin–Shapiro 多项式的单位圆上界，有 $|v_k|\le4\sqrt2 c<1$。乘 $3$ 在每个给定特征阶上都是置换，所以每一阶的总和、每一个循环层的观测都严格为零。[Rudin 的有界系数多项式构造](https://doi.org/10.1090/S0002-9939-1959-0116184-5)

施加前篇的导子和阶筛选仍得到同样的幂次。素数模数上每个非主特征都是本原的，归零低阶整层仍保留所有循环和为零。被删系数至多 $V\tau(g)=x^{o(1)}$ 个，故每个剩余类的值只改变 $x^{o(1)}/g=o(m^{-1/2})$。归一化四阶矩仍为 $\gg1/(gm)$。前篇的短关系四元组只消耗 $g^{-2}x^{o(1)}=o(1/(gm))$。由固定模 $6$ 的素数定理，在这些素数上剩余四阶和之和为 $\gg x/(D^2\log D)$。这是有限傅里叶数据的构造，说明只改善未平移标量的精度会遗漏什么。实际零点仍具有模型未施加的算术约束。[前篇的四阶归一化与筛选](common-windows-and-four-character-bounds.md)

## 3. 在提升链取标量之前移动陪集

在归档的 `Frame.Shares` 中，设 $C=A/A_0$、$|C|=k$。非空顶端点集 $TT$ 在 $C$ 上的经验分布是 $\nu(u)=\#\{t\in TT:t\bmod A_0=u\}/|TT|$。原稿的目标分布是 $\lambda_0(u)=(3+\mathbf 1_{u=0})/(3k+1)$，每个陪集的误差至多 $\delta=t^{-1+\epsilon}$，此处 $t$ 表示几何尺度。原稿 `tops_sum` 通过三角不等式将单个非平凡傅里叶系数的误差控制为 $k\delta$。其证明在这一步仍有完整的复数特征和。[`Shares`、`EAtH` 与原始误差](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Formalism/SharpRestriction/RamificationFrobenius/TotientGaloisClasses/ExplicitFormulaErrorOfHostConditions.lean#L135-L151) · [`tops_sum` 中取绝对值的位置](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/PrimitiveCosets/TopsSum.lean#L54-L77)

把额外的一份权重移到 $b_0\in C$，目标成为 $\lambda_{b_0}(u)=(3+\mathbf 1_{u=b_0})/(3k+1)$。用 $\widehat\nu(\chi)=\sum_{u\in C}\nu(u)\overline{\chi(u)}$。对每个非平凡特征，
$$
 (3k+1)\widehat\lambda_{b_0}(\chi)=\overline{\chi(b_0)}. \tag{2}
$$
这个等式只涉及有限群。把一个满足 `Shares` 的点集平移，确实给出相应的平移陪集计数。要将它代入原几何论证，还须使局部可除性、整格与范数比较同时成立。

**引理 2（对全部频率估计误差）。** 设 $e=\nu-\lambda_{b_0}$，$\nu$ 是概率分布，$(k^{-1}\sum|e|^2)^{1/2}\le\delta$。则
$$
 \sum_{\chi\ne1}|\widehat e(\chi)|^2=k\sum_{u\in C}|e(u)|^2\le k^2\delta^2.
 \tag{3}
$$
因而对 $|v_\chi|\le M$，用近似权重替代 (2) 后，未归一化的整层误差至多
$$
 (3k+1)\left|\sum_{\chi\ne1}\widehat e(\chi)v_\chi\right|
 \le (3k+1)k\sqrt{k-1}\,\delta M . \tag{4}
$$
这是 Parseval 恒等式及 Cauchy 不等式。给每一频率分别使用 `tops_sum` 再求和会多损失约 $\sqrt k$。

在仅知上述大小条件时，(4) 的 $k^{5/2}$ 幂次不能降低。取 $C=\mathbb Z/k\mathbb Z$，$k$ 为奇素数，在零陪集的目标分布上加 $e(u)=\delta(\frac uk)$，其中 $0<\delta\le3/(3k+1)$。这仍是概率分布。对全部 $k-1$ 个非零频率，二次 Gauss 和给出 $|\widehat e|=\delta\sqrt k$。也可直接平方，令 $u=av$ 并对 $v\ne0$ 求和，$a=1$ 的内和为 $k-1$，其余为 $-1$。选 $v_\chi=\overline{\widehat e(\chi)}/|\widehat e(\chi)|$，左侧恰为 $(3k+1)(k-1)\sqrt{k}\,\delta$。这个最优性只针对单个有限群上的误差引理。

## 4. 由平移循环层得到 EH 的判据

令 $G=\widehat{U_q}$，$g=|G|$。为兼容非循环单位群，对每个非平凡循环子群 $H\le G$ 取 $C_H=\widehat H$，其阶记为 $k$。剩余类 $b$ 在 $C_H$ 中的像是 $b_H:\chi\mapsto\chi(b)$。设在端点 $t$ 给出概率分布 $\nu_{H,b,t}$，它与目标 $\lambda_{b_H}$ 的均方根误差至多 $t^{-1+\epsilon}$。原来的逐陪集精度蕴含这一大小条件。端点、子群和剩余类之间允许选用不同的近似。定义
$$
 T_H(b,t)=(3k+1)\sum_{\chi\in H,\,\chi\ne1}
          \widehat\nu_{H,b,t}(\chi)\,g_{\chi^*}(t),\qquad T_{\{1\}}=0,
$$
$$
 \widetilde F_q(b,t)=\frac1g
 \sum_{\substack{H\le G\\H\ {\rm cyclic}}}\ 
       \sum_{K\le H}\mu(|H|/|K|)\,T_K(b,t). \tag{5}
$$
这些量均为实数，因为共轭特征成对出现。若各分布恰为目标分布，(5) 等于 (1)。事实上，内层 Möbius 反演只留下恰生成 $H$ 的特征，每个非主特征恰出现一次。对循环 $G$ 也可以直接用整个 $G$ 的平移层。

**定理 3（原精度的稳定范围）。** 写 $L=\log x$。固定 $A>0$、$0<\theta<3/5$，取
$$
 Q\le x^\theta,\quad h=L^{-A-4},\quad y_0=h x,\qquad
 0<\epsilon<\min\bigl(1/52,(3-5\theta)/2\bigr).
$$
设 (5) 的误差条件在 $y_0\le t\le(1+h)x$ 成立。用同一差分算子
$$
 \mathcal D_h f(y)=\frac{f((1+h)y)-(1+h)f(y)}{h}.
$$
记 $\psi(y;q,b)=\sum_{n\le y,\ n\equiv b\bmod q}\Lambda(n)$ 和 $E(x,Q)=\sum_{q\le Q}\max_{b\in U_q}\sup_{2\le y\le x}|\psi(y;q,b)-y/\varphi(q)|$。以 Dirichlet-GRH 为输入，有
$$
 \left|E(x,Q)-
  \sum_{q\le Q}\max_{b\in U_q}\sup_{y_0\le y\le x}
        |\mathcal D_h\widetilde F_q(b,y)|\right|
 \ll_{A,\theta,\epsilon} x L^{-A}. \tag{6}
$$
这给出了全部模数、最大剩余类、最大端点的判据。要推出 EH，仍需证明 (6) 中平移中心化量的总和本身为 $O(xL^{-A})$，对每个 $A>0$ 成立。

**证明的第一步，控制循环层合并时的损失。** 令 $M(t)=\max_{\chi\ne1}|g_{\chi^*}(t)|$（空集取零）。GRH 下 $M(t)\ll\sqrt t\log^2(2qt)$。具体地，对本原非主特征，$\Psi_\chi(v)=\sum_{n\le v}\Lambda(n)\chi(n)\ll\sqrt v\log^2(2qv)$。分部求和给出 $(L'/L)(1,\chi)=-\sum_{n\ge1}\Lambda(n)\chi(n)/n$，从而 $g_\chi(t)=-t\int_t^\infty\Psi_\chi(v)v^{-2}dv$。这也解释了原稿中线性项的作用。[Montgomery–Vaughan，等差数列中素数的 GRH 估计](https://public.websites.umich.edu/~hlm/math775/ch20.pdf)

用 (4) 代入 (5)。对 $d=|H|$，有 $\sum_{K\le H}|\mu(d/|K|)|\,|K|^{5/2}\le\zeta(5/2)d^{5/2}$。又
$$
 \sum_{H\ {\rm cyclic}}|H|^{5/2}
 =\sum_{\chi\in G}\frac{\operatorname{ord}(\chi)^{5/2}}
                         {\varphi(\operatorname{ord}(\chi))}
 \ll g^{5/2}\log(2g).
$$
这里按子群的生成元计数，并用 $d/\varphi(d)\le\zeta(2)\sum_{j\mid d}1/j\ll\log(2d)$，前一不等式来自逐素因子比较。于是逐模数、逐剩余类、逐端点有
$$
 |\widetilde F_q(b,t)-F_q(b,t)|
 \ll g^{3/2}\log(2g)\,t^{-1+\epsilon}M(t). \tag{7}
$$
在两个端点用 (7)，再对 $q\le Q$ 求和，全部扰动为
$$
 \ll_A Q^{5/2}x^{-1/2+\epsilon} L^{O_A(1)} . \tag{8}
$$
因 $t\ge hx$ 且 $h^{-1}$ 是对数幂，差分保留 (8) 的幂次。这一步不要求两个端点恰好用了同一批辅助对象。由 $\theta,\epsilon$ 的选择，(8) 比 $x$ 小一个固定幂。例如 $Q=x^{4/7}$ 时，其幂次为 $13/14+\epsilon$。

**第二步，由斜坡权回到素数计数。** 对单个 $(t/n-1)_+$ 作 $\mathcal D_h$，在 $n\le y$ 上权为 $1$，在 $y<n\le(1+h)y$ 上为 $((1+h)y/n-1)/h\in[0,1]$，在其余处为零。原稿的线性项被精确消去。对 $\chi^*$ 换回模 $q$ 的 $\chi$，只涉及 $p\mid q$ 的幂，每一素数在上述权下至多贡献 $\log(2x)$。经过 (1) 的归一化后，此误差每个模数为 $O(\omega(q)L)$。特征正交性把差分化为单个既约剩余类的计数减去既约类的实际均值。去掉短区间尾部的误差至多 $2(hx/q+1)\log(2x)$。用 $\zeta$ 的 RH 以及 $p\mid q$ 的素数幂界，把实际均值换为 $y/g$ 的代价是 $O((\sqrt x L^2+\omega(q)L)/g)$。

最后 $\sum_{q\le Q}1/\varphi(q)\ll\log(2Q)$，可由 $q/\varphi(q)=\sum_{d\mid q}\mu^2(d)/\varphi(d)$ 推出。对 $y<y_0$，锐截止的误差直接用 $(y_0/q+1)L+y_0/g$。合计这些算术误差为 $O(hxL^2+QL^2+\sqrt{x}L^3)$。它们对端点和剩余类一致，故可先取最大值。结合 (8) 得 (6)。

**可以怎样改进几何输入。** 如果在阶 $k$ 的层将陪集频率的均方根误差改进为 $t^{-1+\epsilon}k^{-\beta}$，$0\le\beta\le1/2$，(4) 中的 $k^{5/2}$ 就降为 $k^{5/2-\beta}$。在两次子群求和中作同样替换（$\zeta(5/2-\beta)$ 仍收敛），(8) 变成 $Q^{5/2-\beta}x^{-1/2+\epsilon}L^{O_A(1)}$。因此可容许的误差范围扩展至 $\theta<3/(5-2\beta)$，其中 $\epsilon$ 充分小。例如 $\beta=1/4$ 给出 $\theta<2/3$。这项计算只改善陪集近似的预算。

## 5. 直接改进原稿的 Bernoulli 选点

`select` 的证明先估计全部陪集计数误差的**平方和**，最后才输出各个陪集的上界。`topsGroup` 使用独立 Bernoulli 变量，其成功概率在 $[3\alpha/4,\alpha]$ 之间。对坏 Frobenius 轨道使用的 `not_good_le` 只要求这个区间条件。因此可以在选点前指定加重的陪集，也可以在同一次选择中加入一个复数线性统计量。[`select` 的二阶矩及见证的选择](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/ThickFrames/select.lean#L93-L158) · [`topsGroup` 使用轨道上各点概率的地方](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/ThickFrames/topsGroup.lean#L89-L123) · [`not_good_le` 的量词](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/ThickFrames/NotGoodLe.lean#L84-L90)

**定理 4（保留坏轨道预算的平移选点）。** 设一个 $N$ 元集合等分为 $k\ge2$ 个陪集，以有限阿贝尔群 $C$ 标记，给定 $b_0\in C$。在陪集 $b_0$ 上以概率 $\alpha$ 独立选点，其他陪集的概率为 $3\alpha/4$，其中 $0<\alpha\le1$。写 $D_0=\alpha N$，设 $D_0\ge128$。设 $B(T)$ 为所选点集的非负坏轨道损失，$\mathbb E B\le B_0$。给定任意复向量 $(v_\chi)_{\chi\in\widehat C,\chi\ne1}$，存在**同一个**点集 $T$ 满足 $D_0/2\le |T|\le2D_0$、$B(T)\le4B_0$，而且其经验分布 $\nu$ 同时满足
$$
 \left(\frac1k\sum_{u\in C}|\nu(u)-\lambda_{b_0}(u)|^2\right)^{1/2}
       \ll (kD_0)^{-1/2}, \tag{9}
$$
$$
 \left|(3k+1)\sum_{\chi\ne1}
   \big(\widehat\nu(\chi)-\widehat\lambda_{b_0}(\chi)\big)v_\chi\right|
       \ll (3k+1)\frac{(\sum_{\chi\ne1}|v_\chi|^2)^{1/2}}{\sqrt{D_0}}. \tag{10}
$$
常数均为绝对常数。(9) 在选定 $T$ 后对所有测试向量有用，(10) 针对事先固定的向量。后者也可同时对两个预先指定的向量成立，只需调整绝对常数。

**证明。** 用 $C_u$ 表示陪集计数，$\mu_u=\mathbb E C_u$，$\mu=\sum_u\mu_u$。有 $3D_0/4\le\mu\le D_0$ 且 $\mu_u/\mu=\lambda_{b_0}(u)$。写 $Y_u=C_u-\mu_u$。独立性给出 $\mathbb E\sum_u Y_u^2\le D_0$。再令
$$
 W(u)=\sum_{\chi\ne1}\overline{\chi(u)}v_\chi,\qquad
 w_0=\sum_u\lambda_{b_0}(u)W(u)=W(b_0)/(3k+1).
$$
若 $I_a$ 是一个点的入选指标，$\pi_a=\mathbb E I_a$，$c(a)$ 是其陪集，置
$$
 R=\sum_a(I_a-\pi_a)(W(c(a))-w_0).
$$
按独立性与 Parseval，利用 $W$ 的均值为零，得到
$$
 \mathbb E |R|^2\le D_0\!\left(\sum_{\chi\ne1}|v_\chi|^2+|w_0|^2\right)
                  \le 2D_0\sum_{\chi\ne1}|v_\chi|^2 .
$$
点集大小 $S=|T|$ 的方差也至多为 $D_0$。对 $\sum Y_u^2$、$(S-\mu)^2$、$|R|^2$ 和 $B$ 使用 Markov 不等式，分别取各自期望上界的八倍、八倍、四倍、四倍。四个坏事件的概率之和至多 $3/4$。故可同时取 $\sum Y_u^2\le8D_0$、$|S-\mu|^2\le8D_0$、$|R|^2\le8D_0\sum|v_\chi|^2$、$B(T)\le4B_0$。期望为零的非负项在概率一的集合上为零。有两个向量时，将两个 $|R|^2$ 各自的倍数改为八即可。

取这样的 $T$。有 $|S-\mu|\le\sqrt{8D_0}\le D_0/4$，故 $D_0/2\le S\le2D_0$。又 $\nu-\lambda=(Y-\lambda(S-\mu))/S$，而 $\|\lambda\|_2\le1$，遂得 (9)。最后 $R/S=\sum_u(\nu(u)-\lambda(u))W(u)$，得 (10)。原稿的轨道损失具体可取 $B(T)=\sum_p\omega_p\,\#\{\text{属于 }T\text{ 的坏轨道标记点}\}$。在 `topsGroup` 的轨道长度、大小与权重条件下，其估计给出 $B_0=C_b t^{-2}D_0\sum_p\omega_p$，$C_b$ 依赖轨道参数。选择任何 $b_0$ 都保持逐点概率区间，故同一坏轨道预算适用。这里改进的是点集与轨道的选择步骤。[原稿的轨道条件和损失定义](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/ThickFrames/topsGroup.lean#L64-L78) · [轨道概率的数值估计](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/ThickFrames/orbitTerms.lean#L97-L105)

**推论 5（对 EH 近似误差的影响）。** 原选择层的输入为 $\alpha=t^{-\epsilon/8}$、$N\ge t^{2+8\epsilon}/2$，故 $D_0\ge t^{2+63\epsilon/8}/2$。由 (9)，每个陪集的误差也至多 $O(D_0^{-1/2})$，对大 $t$ 即满足原稿输出的 $t^{-1+\epsilon}$ 精度。原证明在将计数归一化时使用了 $64C_o^4\le t^\epsilon,\ k\le C_o$。定理 4 将总点数的方差单独纳入选择，这个归一化步骤只需 $D_0\ge128$。轨道长度和其他几何输入仍须分别核对。若各循环子群和端点具备这一选点空间及上述轨道输入，(9) 将第 4 节的 $\beta$ 取到 $1/2$。因此对单纯的陪集近似，(6) 在 $\theta<3/4$ 仍成立。[密度与大小指数的定义](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/ThickFrames/MarkedLow.lean#L23-L29) · [原来的四次方条件](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/ThickFrames/topsGroup.lean#L12-L34)

进一步在每个 $H,b,t$ 选择 $v_\chi=g_{\chi^*}(t)$。由 (10)，每层误差至多 $O(k^{3/2}M(t)t^{-1})$。将定理 3 中 $5/2$ 次幂的两次子群求和都换成 $3/2$ 次幂，因 $\zeta(3/2)<\infty$，合并后单模数误差为 $O(g^{1/2}\log(2g)\,M(t)t^{-1})$。作两端点差分并对模数求和，(8) 改进为
$$
 O_A\!\left(Q^{3/2}x^{-1/2} L^{O_A(1)}\right). \tag{11}
$$
当 $Q=x^{1-\eta}$、$\eta>0$ 固定时，(11) 为 $x^{1-3\eta/2}$ 乘以对数幂。例如 $Q=x^{3/5}$ 时，模数求和后的这项误差仅为 $x^{2/5}$ 乘以对数幂。在这一选择下，(6) 的近似等式覆盖每个 $\theta<1$。每个剩余类和端点可以选自己的 $T$，误差常数一致。这个存在结论利用已指定的 $g_{\chi^*}(t)$ 作为测试向量。

**仍需从提升链证明的估计。** 原链在窗口上给含共同 $G_\zeta$ 的单侧标量上界。EH 需要对 (5) 中同一批点集的平移中心化量作带符号估计，使 (6) 中的差分总和为 $O_A(xL^{-A})$。相减两个单侧上界不能得到这种估计。可检验的对象是提升链内积中原测试向量与其平移的混合项，须证明它们给出 $\Re(\widehat\nu(\chi)g_{\chi^*})$，其中既有实部的乘积，也有虚部的乘积。另一方面，需要在 EH 端点将几何商群识别为相应的 $\widehat H$，并实现这些特征的导子。原辅助数域在大导子上仍有 $4r\le\lfloor P^{\epsilon/5}\rfloor$ 等限制。定理 4 在给定有限空间上选择点集，构造适用的大导子数域仍是独立问题。[链的加权内积和伴随](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Formalism/ThetaLiftIndependence/Statement.lean#L18-L39) · [辅助数域在大导子上的条件](common-windows-and-four-character-bounds.md)

## 6. 有限复核

模 $5$ 的两个四阶特征（Conrey 编号 $2,3$）的第一正零点纵坐标分别为 $6.183578195451$ 和 $4.132903705213$。前篇的临界带计数和变号区间表明，这两个高度各只属于一个非主特征。原循环标量给二者相同的权，在剩余类 $b=2$ 处，(1) 的归一化前留数系数分别为 $-i$ 与 $+i$。程序还在 $p=383$、$m=16$ 对命题 1 后的例子作了整数检验。每个循环层的和严格为零，而 $gmM_4/c^4=4$。

取端点 $240,270$，使用 120 位球算术从 Dirichlet $L$ 函数的一阶 Taylor 展开求 $(L'/L)(1,\chi^*)$。对模 $8,15,35,45$ 的每个既约类，先由全部循环层反演 (1)，再作差分并扣去诱导特征缺失的素数幂，与从素数幂直接算出的中心化加权计数核对。所有差值的球均含零，模的上界均小于 $9\times10^{-31}$。另给两个端点加入相反的有理数陪集误差，逐项验证了 (3) 及 (4) 求和后的界。在阶 $7,17,41$ 的群上还复核了产生 (4) 最优幂次的 Gauss 和。定理 4 的独立抽样示例取 $N=1024,k=4$，同时控制两个测试向量及相邻点对数这一示例损失。所得点集大小为 $437$，相邻点对数为 $178$，其期望为 $168$。

复算用 Python 3、`python-flint==0.9.0`，在本目录运行 `python -B verify-theta-translates.py`。另附 Lean 4.33.1 / Mathlib 4.33.1 对两条代数与概率引理的核验，分别是有限多个期望预算下的同时选点，以及 (10) 所用的中心化恒等式。两个定理仅使用 `propext`、`Classical.choice`、`Quot.sound`。从原稿的轨道输入推导 (9)–(11) 的证明见正文。[复算程序](verify-theta-translates.py) · [共用的零点与特征程序](verify.py) · [Lean 引理](ThetaSelection.lean)
