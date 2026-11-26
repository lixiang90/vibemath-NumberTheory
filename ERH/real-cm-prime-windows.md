---
title: "全实域和 CM 域中成对的二次 Hecke 素数窗口"
date: "2026-11-26"
abstract: "固定全实基域，在同一辅助素理想上同时控制全实二次特征及其 CM 扭转的斜坡和，保留各自真实的对数导数主项。对几乎所有辅助素理想对，六个特征在同一幂次窗口都有平方根级估计。CM 一侧给出平均 Faltings 高度的有限素理想和近似。"
---

# 全实域和 CM 域中成对的二次 Hecke 素数窗口

设 $F$ 是固定的全实数域，$d=[F:\mathbb Q]$。对 $\mathcal O_F$ 的非零整理想写 $N\mathfrak n$ 为绝对范数。对本原非平凡有限阶 Hecke 特征 $\theta$ 写 $\ell_\theta=L_F'(1,\theta)/L_F(1,\theta)$，在整除导子的素理想处将特征值延拓为零。定义

$$
 R_\theta(x)=\sum_{N\mathfrak n\le x}\Lambda_F(\mathfrak n)\theta(\mathfrak n)
                 \left(\frac{x}{N\mathfrak n}-1\right),\qquad
 G_\theta(x)=R_\theta(x)+x\ell_\theta . \tag{1}
$$

其中 $\Lambda_F(\mathfrak p^j)=\log N\mathfrak p$，在其余理想处为零。Alpöge–Claude 的 Lean 文件在有理数域的非主特征分支使用的就是这个符号和这个主项。这里允许理想和各实位参与特征，并把估计的量词落实到一族真正的数域扩张。[原证明的 Gchi 定义](https://github.com/leanprover-community/riemann-hypothesis/blob/9c9f665ed79ce0fa8da61bd803b7a4bf1ad62ed0/src/Development/AutomorphicArithmetic/DedekindLFunctions/CharFamily.lean#L59) · [原 Lean 证明的核验范围](https://math.mit.edu/~drew/grh/)

## 同一窗口中的全实扩张与 CM 扩张

预先指定 $F$ 的任意 CM 二次扩张 $E_0/F$，记其二次 Hecke 特征为 $\eta_0$。例如可取 $E_0=F(i)$。构造 Fisher–Friedberg 在 $F$ 上的二次 Hecke 族 $\xi_{\mathfrak q}$，在辅助集合中纳入 $\eta_0$ 的全部有限分歧处，并取辅助理想 $\mathfrak c$ 为其有限导子的倍数。素理想 $\mathfrak q\nmid\mathfrak c$ 所对应的 $\xi_{\mathfrak q}$ 在所有实位平凡，有限导子为 $\mathfrak q\mathfrak c_{\mathfrak q}$，其中 $\mathfrak c_{\mathfrak q}\mid\mathfrak c$。它对应一个全实二次扩张 $R_{\mathfrak q}/F$。乘积 $\eta_0\xi_{\mathfrak q}$ 在所有实位为奇，对应双二次扩张 $R_{\mathfrak q} E_0/F$ 中的第三个二次子域，它是 CM 域。两个新域都在 $\mathfrak q$ 分歧。以下特征乘积均取本原代表。[Goldmakher–Louvel，定义、定理 1.1 与第 2 节的导子计算](https://arxiv.org/pdf/1112.1642)

**窗口定理。** 固定 $0<a<1$、$0<\delta<1/2$ 和 $\varepsilon>0$。对充分大的 $Q$，在 $Q<N\mathfrak q\le2Q$ 中除去
$$
 O_{F,\mathfrak c,a,\delta,\varepsilon}(Q^{1-2a\delta+\varepsilon}) \tag{2}
$$
个素理想外，对 $\theta=\xi_{\mathfrak q}$ 和 $\theta=\eta_0\xi_{\mathfrak q}$ **同时**有
$$
 |G_\theta(x)|\le x^{1/2+\delta}
 \qquad (Q^a\le x\le Q,\ x\in\mathbb R). \tag{3}
$$
还可以选两个不同的素理想 $\mathfrak q,\mathfrak r$，两者的范数都在 $(Q,2Q]$。除去 $O_{F,\mathfrak c,a,\delta,\varepsilon}(Q^{2-2a\delta+\varepsilon})$ 个有序对外，(3) 在同一窗口对式 (4) 的六个特征同时成立。
$$
 \xi_{\mathfrak q},\quad \xi_{\mathfrak r},\quad \xi_{\mathfrak q}\xi_{\mathfrak r},
 \qquad \eta_0\xi_{\mathfrak q},\quad \eta_0\xi_{\mathfrak r},\quad
 \eta_0\xi_{\mathfrak q}\xi_{\mathfrak r}. \tag{4}
$$

例如 $Y=Q^{2/3}$ 给出 $[Y,Y^{3/2}]$。取 $\delta=1/10$，单个辅助素理想的例外数至多为 $O(Q^{13/15+\varepsilon})$。允许预先给 $\mathfrak q,\mathfrak r$ 分别指定一个**固定**的窄射线类。每类有 $\asymp Q/\log Q$ 个素理想，故取足够小的 $\varepsilon$ 时，上述两条结论仍覆盖密度一。这里用的是固定基域上的射线类素理想定理。

## 大筛与真实对数导数之间的尾项

我们证明一个对固定扭转一致的均方估计。设 $\omega$ 是 $F$ 的有限阶 Hecke 特征，$N\mathfrak f_\omega\le Q^A$，$A$ 固定。对每个 $\mathfrak q\nmid\mathfrak c\mathfrak f_\omega$ 取 $\theta_{\mathfrak q}=\omega\xi_{\mathfrak q}$。$\mathfrak q$ 的分歧保证 $\theta_{\mathfrak q}$ 非主，并有 $N\mathfrak f_{\theta_{\mathfrak q}}\ll_\mathfrak c Q^{A+1}$。对 $2\le V\le Q$，

$$
 \sum_{Q<N\mathfrak q\le2Q}^{\prime}
 \ \sup_{V\le x\le2V}|G_{\theta_{\mathfrak q}}(x)|^2
 \ll_{F,\mathfrak c,A,\varepsilon} Q^{1+\varepsilon}V . \tag{5}
$$

撇号表示与 $\mathfrak c\mathfrak f_\omega$ 互素的素理想。证明用到数域二次大筛，以及 Cureton–Thorner 的全高度无零半平面。

**第一步。** 对 $W\ge2$，在 $W<N\mathfrak p\le2W$ 中按范数排序，考虑任意前缀上的一次素理想和，系数为 $\omega(\mathfrak p)\log N\mathfrak p$，略去 $\mathfrak p\mid\mathfrak c$。乘上 $\xi_{\mathfrak q}(\mathfrak p)$ 后，对 $\mathfrak q$ 求平方和，前缀最大值的 $\ell^2$ 范数至多为
$$
 \ll_{F,\mathfrak c,\epsilon} (QW)^\epsilon\sqrt{(Q+W)W}. \tag{6}
$$
对固定的系数子集，这是 Goldmakher–Louvel 定理 1.1 的对偶形式。二次特征为实数，交换两边长度并取转置即可得到这里的指标顺序。外和可以由素理想扩大到平方自由理想。每个前缀可分成 $O_F(\log(2W))$ 个二进区块，同一层各区块的系数平方和相加至多为 $O_F(W\log^2(2W))$。按层作三角不等式即得 (6)。这项估计只用 $|\omega(\mathfrak p)|\le1$，对导子随 $Q$ 变化的 $\omega$ 同样有效。[数域二次大筛定理 1.1](https://arxiv.org/abs/1112.1642)

**第二步。** 写 $b_\theta(n)=\sum_{N\mathfrak n=n}\Lambda_F(\mathfrak n)\theta(\mathfrak n)$。在横坐标 $63/64$ 上，
$$
 \left|\frac{L'_F}{L_F}(63/64+it,\theta_{\mathfrak q})\right|
 \ll_{F,\mathfrak c,A}\log\big(2Q(2+|t|)^d\big).
$$
事实上 Cureton–Thorner 的推论 1.2 给出 $\beta\le31/32$，所以这条线与全部非平凡零点有至少 $1/64$ 的水平距离。将完备函数的 Hadamard 对数导数在 $\sigma+it$ 与 $2+it$ 相减，每个零点的差由 $O((1+|t-\gamma|^2)^{-1})$ 控制，对 $63/64\le\sigma\le2$ 一致。单位高度区间的零点数为 $O_F(\log(2Q^{A+1}(2+|t|)^d))$，Gamma 项由 Stirling 公式控制，$2+it$ 处用欧拉乘积。取 $Y_0=Q^{64}$，对指数核作 Mellin 反演，把积分从 $\Re w=1$ 移至 $-1/64$，只跨过 $\Gamma(w)$ 在 0 处的极点，遂得
$$
 -\ell_{\theta_{\mathfrak q}}
  =\sum_{n\ge1}\frac{b_{\theta_{\mathfrak q}}(n)}n e^{-n/Y_0}
      +O_{F,\mathfrak c,A}(Q^{-1}\log(2Q)). \tag{7}
$$
这一步将大筛可处理的平滑和接上了 (1) 中的真实主项。[Cureton–Thorner，定理 1.1 和推论 1.2](https://hexagonmath.org/2610.00200v2) · [该无零半平面的 Lean 形式化](https://github.com/AxiomMath/DedekindZetaZeroFree) · [Thorner–Zaman，Hecke 零点的局部计数，引理 2.5](https://msp.org/ant/2019/13-5/ant-v13-n5-p02-p.pdf)

**第三步。** 因为 $|b_\theta(n)|\le d\Lambda_{\mathbb Q}(n)$，在 $n\le x$ 时以 1 替代 $e^{-n/Y_0}$，对 (1) 的代价为 $O_F(x^2/Y_0)$。于是 $G_\theta(x)$ 等于式 (8) 的两项之和，余项为 $O_{F,\mathfrak c,A}(x^2/Y_0+xQ^{-1}\log(2Q))$。
$$
 -\sum_{n\le x}b_\theta(n),\qquad
 -x\sum_{n>x}\frac{b_\theta(n)}n e^{-n/Y_0}. \tag{8}
$$
在每个区块 $[W,2W]$，二次及以上素理想幂的总绝对权为 $O_F(\sqrt W\log^2(2W))$，可按有理素数幂计数。固定的 $\mathfrak p\mid\mathfrak c$ 在一次项中只贡献 $O_{F,\mathfrak c}(1)$。它们纳入 $\ell^2$ 范数后至多多出 $O_F(\sqrt Q)$ 倍，均在 (6) 的量级内。因而 (6) 也适用于含所有素理想幂的最大前缀。

对 (8) 第一项将 $W\le2V$ 的区块相加。对第二项在 $W=2^j V$ 上分部求和，单调权的大小与总变差均至多为 $1/W$，首个不完整区块可用两段前缀之差。超过 $Y_0^2$ 的部分由指数衰减直接略去。余下区块对 $\ell^2$ 范数的总贡献至多为
$$
 \ll Q^\epsilon V
 \left(\sum_{V\le W\le Q}^{\rm dyad}\sqrt{Q/W}
       +\sum_{Q<W\le Y_0^2}^{\rm dyad}1\right)
 \ll Q^{2\epsilon}\sqrt{QV}.
$$
这里可在 (6) 先取任意小的指数，因为 $W\le Q^{128}$。第一项和 (8) 的显式误差也有同一上界。重选 $\epsilon$ 得 (5)。

对 $Q^a\le x\le Q$ 以二进尺度分段，用 (5) 和 Markov 不等式将例外数界为 $O(Q^{1+\epsilon}\sum_{V\ge Q^a/2}^{\rm dyad}V^{-2\delta})$，即 (2)。取 $\omega=1,\eta_0$ 得成对的全实及 CM 扩张。对每个保留下来的 $\mathfrak q$，再以 $\mathfrak r$ 为变量，分别取 $\omega=1,\eta_0,\xi_{\mathfrak q},\eta_0\xi_{\mathfrak q}$。这些扭转的有限导子范数均为 $O_{F,\mathfrak c}(Q)$，故 (5) 的常数对第一步所选的 $\mathfrak q$ 一致。求和即得 (4) 的例外数。若预先指定有限多个 $E_0/F$，可把各自的有限分歧处一并纳入辅助集合，对例外集取有限并。相同幂次界遂对这些 CM 扭转同时成立，常数容许依赖所指定的有限族。

## CM 高度给出的主项

设 $E/F$ 是任意 CM 二次扩张，$\theta=\eta_{E/F}$，$C=|D_E|/D_F$。这里 $D_F>0$ 为绝对判别式。记 $\bar h_E$ 为 Andreatta–Goren–Howard–Madapusi Pera 定理 A 所用归一化下，$2^d$ 个 CM 型的稳定 Faltings 高度平均。该定理与 Yuan–Zhang 的独立证明给出
$$
 \bar h_E=-\frac12\frac{L'_F}{L_F}(0,\theta)
          -\frac14\log C-\frac d2\log(2\pi).
$$
完备因子为 $C^{s/2}\Gamma_{\mathbb R}(s+1)^d$，其中 $\Gamma_{\mathbb R}(s)=\pi^{-s/2}\Gamma(s/2)$。对函数方程取对数导数，利用 $\psi(1/2)+\psi(1)=-2\gamma-2\log2$，得到
$$
 \bar h_E=\frac12\ell_\theta+\frac14\log C
           -\frac d2\gamma-d\log(2\pi). \tag{9}
$$
式中 $\gamma$ 为欧拉常数。任意两个同一全实基域上的 CM 二次扩张相减，最后两项都精确消去。[Andreatta 等，定理 A 及其高度归一化](https://www.math.mcgill.ca/goren/PAPERSpublic/Colmez_v3.5.pdf) · [Yuan–Zhang 的平均 Colmez 定理](https://arxiv.org/abs/1507.06903)

**推论。** 对窗口定理保留下来的每个 CM 特征，以 $E$ 表示对应的 CM 域。在 (3) 的全部实数尺度上，
$$
 \left|\bar h_E+\frac{R_\theta(x)}{2x}
       -\frac14\log C+\frac d2\gamma+d\log(2\pi)\right|
 \le \frac12 x^{-1/2+\delta}. \tag{10}
$$
因此一个有限素理想和以明确的精度给出平均高度。比如在 $x=Q$，误差至多为 $Q^{-1/2+\delta}/2$。这是 (3) 与 (9) 的直接合取。

## 对全高度 ERH 还缺的输入

对全实 $F$ 及 CM 二次扩张 $E/F$，$\zeta_E=\zeta_F L_F(\eta_{E/F})$，两个因子均在开临界带全纯。因此 ERH$(E)$ 等价于 ERH$(F)$ 与相对 Hecke 函数的 RH 同时成立。特别地，若所有 CM 域满足 ERH，则给任意全实 $F$ 取 $F(i)$ 就得到 ERH$(F)$。全实域的部分与每个固定 CM 二次特征的部分都需要算术输入。

原 Lean 证明的周期命题在给定曲线理想、模型和嵌入指标后，对定义域的共轭嵌入求平均。式 (9) 对一个 CM 域的全部 CM 型求平均。要应用于高维提升链，还需识别链所用行列式线的度量、有限处正规化，以及它实际采用的 CM 型权重。伴篇在一个具体六次域算出了两类高度之间严格为正的差值，也给出了单位平衡的定量界。[原始周期等式的量词](https://github.com/leanprover-community/riemann-hypothesis/blob/9c9f665ed79ce0fa8da61bd803b7a4bf1ad62ed0/src/Development/AbelianVarietyTheory/BlockTheory/GlobalSignature/MeanPeriodLogarithms.lean#L49) · [六次域的两类高度与零点](cm229-heights-and-zeros.md)

(3) 和 (4) 允许导子随着窗口移动，并在固定射线类中保留足够多的辅助素理想。对预先指定的全实域或 CM 域，仍需把辅助特征的估计转成同一个固定因子在任意远窗口上的估计，并满足提升链另外指定的局部分歧、周期整除和上同调条件。这里的例外数也没有给出随 $Q$ 稀薄到任意程度的指定射线类中的素理想对。[原 Lean 的辅助宿主命题](https://github.com/leanprover-community/riemann-hypothesis/blob/9c9f665ed79ce0fa8da61bd803b7a4bf1ad62ed0/src/Development/AutomorphicArithmetic/CharacterBoxes/SupplyChi.lean#L46) · [Robin Zhang 对原论证各步作用的导读](https://www.ihes.fr/~rzhang/files/GRH%20note.pdf)
