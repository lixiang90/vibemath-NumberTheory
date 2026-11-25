---
title: "共同窗口能把四特征零点和压到多小？"
date: "2026-11-25"
abstract: "对预先固定的本原特征族，合并 Alpöge–Claude 证明中的例外尺度，得到共同窗口和辅助素数。同一导子、同一特征核只需选窗一次。对 Elliott–Halberstam 归约中的真实零点和，在每个指定端点证明剩余四阶项的模数块上界。即使四个特征都取素数模数上的最大阶，有限傅里叶模型也能在对数因子内达到该上界。"
---

# 共同窗口能把四特征零点和压到多小？

对预先固定的四个本原 Dirichlet 特征，原证明的几何窗口可以同时选取。可以共用辅助素数，甚至可以处理四个不同的层数归一化常数。若端点由 Elliott–Halberstam（EH）问题预先指定，还有一个定量差距。设模数 $D<q\le2D$，其中 $D\asymp x^{1/2+\delta}$。以归档的 Dirichlet-GRH 为输入，本文对同一端点的剩余四特征和证明了 $x^{-2\delta+o(1)}$ 上界。用四阶矩推出该范围的 EH，需要在 $x^{-4\delta}$ 的尺度上获得足够的节省。一个有限傅里叶构造给出了前一上界在现有范数信息下的下界例子。

下文的“窗口”特指原证明控制辅助数域上同调群的区间。“端点”是等差数列素数计数或相应显式公式的自变量。几何部分从固定版本的 Lean 归档中的中间引理出发，零点和估计另用 Dirichlet-GRH、Brun–Titchmarsh 不等式及特征正交性。[原始发布](https://alpo.ge/remark.pdf) · [归档核验记录](https://math.mit.edu/~drew/grh/) · [前篇的有限零点归约](low-zeros-and-character-relations.md)

## 1. 固定有限族的共同窗口

固定 $0<\epsilon<1/52$、$C_t=401$ 和非空有限本原特征族 $\mathcal F$。用 $X$ 表示几何构造的外层尺度，用 $u,v$ 表示原代码名为 `q,d` 的两个辅助素数。取一个共同的层数函数 $e(t)$，满足
$$
 t^{2+8\epsilon}\le 2^{e(t)}\le 2t^{2+8\epsilon}\qquad(t\ge1).
$$
原稿的 `VAt` 表示对由 $\psi,e(t),u,v$ 指定的每个数域 $E$，都有 `VVle` 给出的加权上同调长度界，权的截断点为 $\lfloor t^{401}\rfloor$，上界为 $\exp(t^\epsilon/2)$。[窗口、数域与 `VAt` 的定义](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/CharacterBoxes/ExistsPageN.lean#L12-L50)

**定理 1。** 存在 $J_{\mathcal F}\ge2$ 和门槛 $X_0$，使对每个 $X\ge X_0$，令 $P=X^{1/J_{\mathcal F}}$，则对原稿中每一对满足 $\mathrm{WinPair}(\epsilon,P)$ 的 $u,v$，存在同一个 $W$ 使
$$
 P\le W,\quad W^2\le X,\quad
 \mathrm{VAt}(\psi,\epsilon,401,u,v,e(t),t)
 \quad(\psi\in\mathcal F,\ W\le t\le W^2). \tag{1}
$$
$J_{\mathcal F}$ 和 $X_0$ 可以依赖该有限族。$W$ 在给定 $X,u,v$ 后选取。

**证明。** 记 $M_\psi$ 为原稿的阿贝尔域数据 `AbF.ofChar`，$C_K\ge1$、$R_0(M_\psi)$、$R_\ell(L_\psi(\ell))$ 取自 `thmBab_of` 所作的同一组选择。最后两项分别控制素数 $2$ 和有限个奇素数处的例外整数个数。记
$$
 R_*=\sum_{\psi\in\mathcal F}
 \left(R_0(M_\psi)+\sum_{\ell\in B_\psi}R_\ell(L_\psi(\ell))\right),
 \qquad \eta=\frac{\epsilon}{2C_K},\quad B=\frac{401}{\eta},
 \quad J_{\mathcal F}=2(2B)^{R_*}. \tag{2}
$$
即使不同特征产生相同的例外，也可以使用这个上界。[常数的选择与单特征的并集](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/AutomorphicVarieties/LayeredIsomorphisms/ThmBabOf.lean#L18-L38)

这里需要回到选窗之前。`partA` 与 `partB` 对**任意给定的** $J>0$ 成立。在两处同时代入 (2) 的 $J_{\mathcal F}$、同一 $X$、同一 $e(t)$ 和截断函数 $\lfloor t^{401}\rfloor$。后两个函数满足原稿 `admitted` 的条件，可取 $A_1=4$。对每个 $\psi$，$\mathrm{WinPair}$ 的大小与同余条件在 $X$ 充分大时保证 `PrimePair`，且
$$
 \log(uv)\le P^{\epsilon/5}(1+\log^2 P)\le P^{\epsilon/4}
       =X^{\epsilon/(4J_{\mathcal F})}.
$$
在 `partA` 中按 `TwoIn` 选择相应的底域，在 `partB` 中沿用原稿给各奇素数选出的底域及秩界。把有限多个门槛取最大值，可令 $P$ 超过它们。这样得到至多 $R_*$ 个例外整数的并集 $\mathcal E$。在 $P\le t\le X$ 上，只要 $(t^\eta,t^{401}]$ 不含 $\mathcal E$ 中的整数，每个 $\psi$ 的 `BaAt` 与 `BbAt` 就同时成立。例外集允许依赖 $X$ 和辅助素数。[`partA` 的量词](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/AutomorphicVarieties/LayeredIsomorphisms/PartA.lean#L24-L32) · [`partB`](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/AutomorphicVarieties/LayeredIsomorphisms/PartB.lean) · [例外条件 `Good`](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/AutomorphicVarieties/LayeredIsomorphisms/SecondBound.lean#L46-L48)

令 $U=\log X$。对每个例外整数 $n>1$ 置 $a=(\log n)/401$，须避开的 $\log t$ 构成半开区间 $[a,Ba)$。若没有起点落在 $(U/(2B),U]$，则 $[U/2,U]$ 已避开全部禁区。否则删去该段内的全部起点，并把当前上限改为 $U/(2B)$，删去的禁区在新上限右侧。至多重复 $R_*$ 次即得到 $[w,2w]$，其中 $w\ge(\log X)/[2(2B)^{R_*}]$。取 $W=\exp w$ 即得所需范围。这样也将原实变量引理分母中的 $4(2B)^{R_*}$ 改进为 $2(2B)^{R_*}$。这个因子在连续起点模型中是最优的。对可选的 $w$，每个禁区恰排除 $[a/2,Ba)$。令 $R_*$ 个这样的区间首尾相接，从 $(1+\xi)U/[2(2B)^{R_*}]$ 开始，它们覆盖至 $(1+\xi)U/2$ 的左极限。令 $\xi\downarrow0$ 即可。[原稿的实变量引理](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/AutomorphicVarieties/LayeredIsomorphisms/Window1.lean#L6-L34)

还能计算共同窗口的个数。预先指定整数 $N>R_*$，在上述构造中改用 $J_N=2(2B)^{N-1}$。令 $U_j=(\log X)/(2B)^j$，$0\le j<N$。一个起点 $a$ 能破坏候选区间 $[U_j/2,U_j]$ 当且仅当 $U_{j+1}<a\le U_j$。这些半开区间两两不交，故预先列出的 $N$ 个候选中至少有 $N-R_*$ 个对整族同时有效。取 $N=R_*+1$ 就恢复 (2)。增大 $N$ 的代价由 $J_N$ 明确给出。[有限族选窗引理的 Lean 证明](CommonWindow.lean)

还需核对从 `BaAt`、`BbAt` 到 `VAt` 的传递。原稿各特征的 `IsBWindow` 带有自己的 $J_\psi$，所以不能直接把新窗口代入该谓词。`thmHV` 的证明在单点 $t$ 使用 `PointHV`，所用量是 $\lfloor t^{401}\rfloor$、`BbAt` 的奇素数界以及辅助素数的大小。固定 $\mathcal F$ 后，$u,v$ 最终均超过 $2X^{401}+1$ 和各特征有限列表中的全部常数，故这一逐点推导在共同窗口成立。从 `BaAt` 取 $2$ 处的界，从 `PointHV` 取加部的界，再取原稿已对所有点给出的减部界。因 $\log(uv)\le P^{\epsilon/4}\le t^{\epsilon/4}$，原稿最后的导子估计和误差吸收逐点适用，得到 (1)。[`PointHV` 的实际调用](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/AutomorphicVarieties/HeightAndVanishing/ThmHV.lean#L22-L40) · [三部分合并为 `VVle`](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/AutomorphicVarieties/CMFieldExponentialBounds/AdmissibleOfPageN.lean#L95-L129)

**共用辅助素数与端点。** 取各导子的最小公倍数 $m$，用 `PairEx` 在 $[P,X]=[P,P^{J_{\mathcal F}}]$ 上选一对满足 `PairGood` 的素数。若 $r\mid m$，模 $r$ 的特征可提升到模 $m$，而与 $2muv$ 互素的整数也与 $2ruv$ 互素。由 `Induces` 的定义，模 $m$ 的扭转界随即涵盖模 $r$ 的扭转。零点盒的底部导子条件 $N\mid8r$ 推出 $N\mid8m$，允许高度随基准模数增大。因此同一对素数对各导子均满足 `PairGood`。[`PairGood` 的两组量词](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/CharacterBoxes/PairGood.lean#L20-L35) · [`Induces` 的定义](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Formalism/SharpRestriction/RamificationFrobenius/TotientGaloisClasses/Analytic.lean#L23-L29)

各特征在 `winPairs` 中允许有自己的常数 $c_i\ge1$。令 $k_i=c_i^{1/(2+8\epsilon)}$，$k=\max_i k_i$。当 $W$ 足够大，所有 $s\in[kW,(kW)^{3/2}]$ 都满足 $W\le s/k_i\le W^2$ 且 $s\le W^2\le X$。在 (1) 中分别用 $t=s/k_i$，取层数 $e_i=e(s/k_i)$。于是 $s^{2+8\epsilon}\le c_i2^{e_i}\le2s^{2+8\epsilon}$。又由 $s/k_i\ge k_i^{400}$ 得 $s^{400}\le(s/k_i)^{401}$。以 (1) 控制这个较短的截断，便对所有 $i$ 同时得到原稿的 $\mathrm{WinDatum}(\psi_i,\epsilon,J_{\mathcal F},c_i,s,P,u,v,e_i)$。这里 $s,P,u,v$ 真正相同，层数和数域可分别取值。[`WinDatum` 的定义](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/CharacterBoxes/TwistForm.lean#L20-L28) · [原有的缩放步骤](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/CharacterBoxes/WinPairs.lean#L55-L84)

**相同特征核只需一组例外。** 若两个本原特征的导子相同，在该导子的单位群上也有相同的核，则对固定的 $e,u,v$，它们的 `IsChiField` 完全等价。定义中涉及特征的两处正是导子的整除关系及 $\psi(a)=1$。非单位上的特征值为零，单位上的这个条件只取决于核。因此 `VAt` 的结论对这两个特征相同，`PairGood` 也只依赖导子。应用 (2) 时，对每一对不同的“导子、核”取一个代表即可。在固定素数模数 $p$ 上，单位群为循环群，给定阶 $d\mid p-1$ 的特征具有同一个核。全部 $p-2$ 个非主特征因此只需 $\tau(p-1)-1$ 个代表，全部阶恰为 $p-1$ 的特征只需一个。这个结论对每个固定 $p$ 成立，$X_0$ 仍允许依赖 $p$。[`IsChiField` 的核条件](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/CharacterBoxes/ExistsPageN.lean#L19-L26)

## 2. 为什么同一窗口尚不能覆盖 EH 的大导子

定理 1 先固定特征族，给出的 $J_{\mathcal F}$ 随例外总数呈指数增长。直接在素数计数的端点 $s$ 运行原辅助数域构造，还遇到一个与窗口间距无关的限制。源码 `supplyOfAdmOf` 在调用 `fieldChi` 时要求
$$
 4r\le b=\lfloor P^{\epsilon/5}\rfloor\le s^{\epsilon/5},
 \qquad \epsilon/5<1/260, \tag{3}
$$
其中 $r$ 是本原导子。前篇零点和保留 $r>\sqrt{x}/(\log x)^{A+4}$，而其所有端点 $s\le x$。充分大的 $x$ 下，两者的范围没有交集。[调用 `fieldChi` 的精确条件](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/ShortNorms/supplyOfAdmOf.lean#L83-L95)

这个调用可以局部改细。令 $o=\operatorname{ord}(\psi)$，保留 `ordAny` 返回的 $e_1$ 及 `fieldChi` 的辅助素数、分裂条件。在 `fieldChi` 的证明中可以改用
$$
 e\ge\max(e_1,2),\qquad 4r\le2^e,\qquad o^{64}\le2^e,\qquad
 4\ell\le b\quad(\ell\mid r,\ \ell\text{ 为奇素数}). \tag{4}
$$
理由是 $H.\mathrm{QS}\mid\operatorname{lcm}(2^{e+2},r)$，所以它的每个奇素因子都整除 $r$。(4) 同时给出这些素数的大小界和使用分裂条件所需的范围。原证明另用 $e\ge4r$ 和 $e\ge64o$ 来保证的两个幂不等式已直接列入 (4)。沿用原来的场和框架即可得到 `fieldChi` 的各项输出。特别地，$r$ 的奇素因子的**指数**不再进入对 $b$ 的要求。[原证明使用这三处界的位置](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/CharacterBoxes/FieldChi.lean#L10-L58)

但在该框架内 $r\mid H.\mathrm{QS}$，而 `FieldRows.small_le` 本身就要求 $4\ell\le b$ 对 $H.\mathrm{QS}$ 的每个奇素因子成立。因此对奇素数导子 $r$，(3) 的限制依然存在。另一条 `FieldRows.o_pow` 要求 $o^{64}\le2^e$。在同一端点处 $2^e\le2s^{2+8\epsilon}$，故可直接容纳的阶满足 $o\le2^{1/64}s^{7/208}$。这些是把该辅助数域构造接到 EH 端点时必须处理的两个具体条件。[`KerChi` 的整除条件](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/CharacterBoxes/TwistForm.lean#L9-L18) · [`FieldRows` 中的素因子与阶](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/CoprimeRanks/SecondSetNoRelation.lean#L13-L23)

## 3. 每个指定端点上的四阶上界

沿用前篇的截断和记号。固定 $A>0$，令 $L=\log x$，$h=L^{-A-4}$，$T=4h^{-1}\sqrt L$，$J=\lceil L^{A+4}\rceil$，$y=jx/J$，$1\le j\le J$。取
$$
 R_x=\sqrt{x}L^{-A-4},\qquad
 K=\lceil(A+6)\log L/\log2\rceil,\qquad V=x^{1/[4(K+2)]}.
$$
只计 $\omega(q)\le K$。对模 $q$ 的特征，若本原导子 $r_\chi>R_x$ 且 $\operatorname{ord}(\chi)>V$，在 Dirichlet-GRH 下定义
$$
 z_\chi(y)=\sum_{\substack{\rho=1/2+i\gamma\ \text{of }L(s,\chi^*)\\|\gamma|\le T}}
 \frac{y^{i\gamma}\exp(h^2\rho^2/2)}{\rho}.
$$
其余特征系数为零。写 $g=\varphi(q)$、$G$ 为模 $q$ 的特征群，并置
$$
 Z_q(a,y)=g^{-1}\sum_{\chi\in G}\overline{\chi(a)}z_\chi(y),
 \qquad M_j(q,y)=g^{-1}\sum_{a\in(\mathbb Z/q\mathbb Z)^\times}|Z_q(a,y)|^j.
$$
令 $H_q(y)$ 为下式中删去短关系四元组后剩余项之和的实部。
$$
 M_4(q,y)=g^{-4}\sum_{\chi,\eta,\nu\in G}
 z_\chi z_{\chi\eta\nu}\overline{z_{\chi\eta}z_{\chi\nu}}.
$$
这里的短关系指存在 $0<\max(|m|,|n|)\le\exp(\sqrt L)$ 使 $\eta^m\nu^n$ 的导子或阶不超过 $\exp(\sqrt L)$。前篇证明每个模数被删项的绝对值之和至多为 $g^{-2}x^{o(1)}$。[前篇定理 1、定理 2](low-zeros-and-character-relations.md)

**定理 2。** 设 $Q\le x^{1-\varepsilon}$，$\varepsilon>0$ 固定，置 $L_1=L\log(2L)$。则对每个实端点 $x/J\le y\le x$（包括上述全部网格点）及 $1\le D\le Q$ 一致有
$$
 \sum_{\substack{D<q\le\min(2D,Q)\\\omega(q)\le K}}|H_q(y)|
 \ \ll_{A,\varepsilon}\ \min\{L_1^4,L_1^2\,y/D^2\}. \tag{5}
$$
当 $D\asymp x^{3/5}$ 时，(5) 给出 $O_{A,\varepsilon}(x^{-1/5}L_1^2)$ 的绝对值上界。若仅超出平方根一个对数因子，取 $D\asymp\sqrt{x}L^b$，则上界为 $O_{A,\varepsilon}(L^{2-2b}\log^2(2L))$。这两个推论都适用于每个指定端点。

**证明。** 先保留全部非主特征，记相应的傅里叶和为 $Z^{\rm all}_q$。前篇的平滑显式公式及诱导特征修正，在单个既约剩余类中把 $-\sqrt y\,Z^{\rm all}_q$ 写成减去实际均值的平滑计数加 $O_A(L^2)$。平滑与锐截止之差为 $O((hy/q+1)L)$，对剩余类取均值也满足此界。因 $y\ge x/J$ 且 $q\le x^{1-\varepsilon}$，充分大的 $x$ 有 $y/q\ge x^{\varepsilon/2}$。

由 Brun–Titchmarsh 不等式，锐截止计数中素数的一部分至多 $2y\log y/(g\log(y/q))\ll_\varepsilon y/g$。对高次素数幂 $p^k$，同余式 $v^k\equiv a\bmod q$ 的解数至多 $k^{K+1}$，因为单位群至多有 $K+1$ 个生成元。每个解类中至多有 $\sqrt{y}/q+1$ 个可能的底数。对 $2\le k\le L/\log 2$ 求和并计入权重，得到 $\exp(O_A((\log L)^2))(1+\sqrt{y}/q)=o(y/q)$。因此对所有既约剩余类和它们的均值都有同一上界。显式公式遂给出 $\|Z^{\rm all}_q\|_\infty\ll_{A,\varepsilon}\sqrt y/g$。[Montgomery–Vaughan，定理 2](https://doi.org/10.1112/S0025579300004708)

按单位长度分组，$|\gamma|\le T$ 中每组含 $O_A(L)$ 个零点。因为 $T$ 是 $L$ 的固定次幂，且 $1/|\rho|\ll(1+|\gamma|)^{-1}$，每个未筛零点系数的模均为 $O_A(L_1)$。令 $C_R(q)=\sum_{r\mid q,\,r\le R}\varphi(r)$。导子不超过 $R_x$ 的特征至多有 $C_{R_x}(q)$ 个，阶不超过 $V$ 的至多有 $x^{1/4}$ 个。逐模数扣去这些系数，再用 Parseval 恒等式，得到
$$
 \|Z_q\|_\infty\ll_{A,\varepsilon}
 \min\!\left(L_1,\frac{\sqrt y+L_1(C_{R_x}(q)+x^{1/4})}{g}\right),
 \qquad M_2=g^{-2}\sum_\chi |z_\chi|^2\ll_A L_1^2/g .
$$
这里也用到了 $|\exp(h^2\rho^2/2)|\le\exp(h^2/8)$。[单位高度的 Dirichlet 零点计数](https://doi.org/10.1090/mcom/3599)

保留 $C_R(q)$ 的除数权可以把对数损失固定下来。写 $f(n)=(n/\varphi(n))^3$。有 $\sum_{n\le z}f(n)^j\ll_j z$（$j=1,2$），因为将 $f^j$ 写成常函数 $1$ 与积性函数 $u_j$ 的卷积时，$u_j$ 只支撑在无平方因子数上，且 $u_j(p)=(1-1/p)^{-3j}-1=O_j(1/p)$，故 $\sum_d u_j(d)/d<\infty$。特别地 $\sum_{n\le z,\,d\mid n} f(n)\ll zf(d)/d$。展开平方，记 $[r,s]$ 为最小公倍数，得到
$$
 \begin{aligned}
 \sum_{D<q\le2D} f(q)C_R(q)^2
 &\ll D\sum_{r,s\le R}\frac{\varphi(r)\varphi(s)f([r,s])}{[r,s]}\\
 &\ll D\sum_{r,s\le R}(r,s)f(r)f(s)
 \ll D R^2\log(2R).
 \end{aligned}
$$
末步展开 $(r,s)=\sum_{d\mid r,\,d\mid s}\varphi(d)$，并用 $\sum_{d\le R} f(d)^2/d\ll\log(2R)$。这里 $f(dm)\le f(d)f(m)$。

现在用 $M_4\le\|Z_q\|_\infty^2M_2$。在模数块上 $\sum 1/g\ll1$，因为 $q/g\le f(q)$。用逐模数上界的第二项和刚证的除数估计，得到
$$
 \sum_q M_4(q,y)\ll_{A,\varepsilon}\frac{L_1^2}{D^2}
 \left(y+L_1^2\bigl(R_x^2\log(2R_x)+\sqrt{x}\bigr)\right)
 \ll_{A,\varepsilon}\frac{L_1^2\,y}{D^2}.
$$
末步用 $R_x^2/y\ll L^{-A-4}$ 与 $\sqrt{x}/y\ll L^{A+4}/\sqrt{x}$。用上界的第一项则得 $\sum_q M_4\ll_A L_1^4$。最后，当 $2D\le R_x$ 时全部 $z_\chi$ 为零。其余模数块满足 $D>R_x/2$，被删短关系项的绝对值总和由前篇控制为 $D^{-1}\exp(O_A(\sqrt L\log L))$。这一量为 $o(1)$，与 $y/D^2$ 的比值亦为 $o(1)$，因为 $D/y\le Jx^{-\varepsilon}$。用三角不等式即得 (5)。

**距离 EH 的量级。** 由 Hölder 不等式，模数块的最大剩余类误差至多由 $\sqrt{x}D^{3/4}(2DJ\max_y\sum_q M_4(q,y))^{1/4}$ 控制，另有前篇已估计的截断误差。使该式比 $x$ 小一个固定幂，一种充分的目标是 $\sum_q M_4\ll x^{-\sigma}x^2/D^4$，其中 $\sigma>0$。当 $D\asymp x^{1/2+\delta}$ 时，(5) 与这个目标分别位于 $x^{-2\delta+o(1)}$ 和 $x^{-4\delta-\sigma}$。在 $\delta=1/10$ 的例子中，前篇对大误差的必要下界为 $x^{-2/5}$ 除以对数幂，与 (5) 相容。

用两个特征之间的相关量可以写出下一步所需的估计。令 $C_q(\eta,y)=g^{-1}\sum_\chi z_\chi(y)\overline{z_{\chi\eta}(y)}$，则由正交性
$$
 M_4(q,y)=g^{-2}\sum_{\eta\in G}|C_q(\eta,y)|^2.
$$
因此所需的额外节省可以在同一端点上，以这些扭转相关量的平方和来检验。

## 4. 上界所用信息的限度

**命题 3。** 固定 $1/2<\theta<1$，$D\asymp x^\theta$，考虑整个 $(D,2D]$。在每个素数模数上存在满足共轭对称的有限傅里叶系数，服从上一节的导子和阶筛选、$|z_\chi|\le1$、$\max_a |Z_q(a)|\ll\sqrt x/D$，使其剩余四阶和满足
$$
 \sum_{D<p\le2D} H_p\gg \frac{x}{D^2\log D}. \tag{6}
$$
这里的系数由下述多项式构造。这个命题检验逐特征界、剩余类上界与正交性所能推出的结论。真实 Dirichlet 零点的特征间依赖仍可用于改进 (5)。

**证明。** 取 $D^2/x\le m<2D^2/x$ 的二次幂，$c=1/8$。对素数 $p\in(D,2D]$，令 $g=p-1$，取原根 $b$，当 $x$ 足够大时 $m<g$。取 Rudin–Shapiro 符号 $s_0,\ldots,s_{m-1}\in\{-1,1\}$。从 $P_1=Q_1=1$ 开始递推 $P_{2n}=P_n+w^nQ_n$、$Q_{2n}=P_n-w^nQ_n$，则 $P_m=\sum_j s_jw^j$。在单位圆上 $|P_m|^2+|Q_m|^2=2m$。[Rudin 的原文](https://doi.org/10.1090/S0002-9939-1959-0116184-5)

在剩余类 $b^j$ 上放置实数 $F(b^j)=cs_j/\sqrt m$（$0\le j<m$），其余取零。对应的未归一化特征系数是 $cP_m(\chi(b))/\sqrt m$，模至多 $c\sqrt2$。素数模数上每个非主特征都是本原的。阶至多 $V$ 的特征至多 $V\tau(p-1)=x^{o(1)}$ 个，将这些系数归零后，每个 $F(a)$ 最多改变 $x^{o(1)}/g=o(m^{-1/2})$。所得 $Z_p$ 仍具有共轭对称的系数，$M_4\gg1/(gm)$。前篇短关系的计数只用系数的模上界，故此处也至多损失 $g^{-2}x^{o(1)}=o(1/(gm))$。于是每个 $H_p$ 均有正下界 $\gg1/(gm)$。由素数定理，该块有 $\gg D/\log D$ 个素数，得到 (6)。在 $y=x$ 时，(5) 与 (6) 只相差对数因子。

**命题 4（四个特征的核也可以相同）。** 在命题 3 的有限傅里叶模型中，可以进一步要求每个非零系数的特征都具有阶 $p-1$，仍有
$$
 \sum_{D<p\le2D} H_p
 \gg_\theta \frac{x}{D^2(\log D)^3[\log(2\log D)]^4}. \tag{7}
$$
每个素数模数上所有参与四元组的特征因而都有同一个核，依第 1 节会给出相同的 `IsChiField` 谓词。这个例子说明，即使对同一核，若只使用 (5) 证明中的系数界和傅里叶范数界，仍能出现 (7) 的大小。要改善**幂次**，需要从零点或素数计数得到额外约束。

**证明。** 写 $g=p-1$、$w=2^{\omega(g)}$、$\Delta=\varphi(g)/g$。由除数配对和 Brun–Titchmarsh 不等式，
$$
 \sum_{D<p\le2D} w
 \le \sum_{D<p\le2D}\tau(p-1)
 \le 2\sum_{d\le\sqrt{2D}}\pi(2D;d,1)
 \ll \frac{D}{\log D}\sum_{d\le\sqrt{2D}}\frac1{\varphi(d)}
 \ll D .
$$
最后一步可由第 3 节 $f$ 的均值分块求得。由素数定理，至少 $\gg D/\log D$ 个素数满足 $w\le C\log D$，其中 $C$ 是足够大的绝对常数。只在这些模数上构造系数，其余模数取零。将 $g$ 的不同素因子按大小编号，第 $i$ 个至少是 $i+1$，故 $\Delta\ge1/(1+\omega(g))\gg1/\log(2\log D)$。

取二次幂 $w^2D^2/x\le m<2w^2D^2/x$。因 $\theta<1$，有 $m<g$。在指数坐标 $j\bmod g$ 上使用命题 3 的 $F$，记 $I=\{0,\ldots,m-1\}$，故 $\sum_j|F_j|^2=c^2$。令 $\Pi$ 为保留与 $g$ 互素的频率的正交投影，恰对应阶 $g$ 的特征。由 Möbius 反演，这个投影是对 $d\mid g$、带系数 $\mu(d)$ 的投影之和，其中第 $d$ 项保留 $d\mid k$ 的频率 $k$。后者在指数坐标上等于对 $j,j+g/d,\ldots,j+(d-1)g/d$ 取平均，故 $\|\Pi\|_{\infty\to\infty}\le w$。

以 $F^{(t)}_j=e^{2\pi itj/g}F_j$ 平移频谱，$0\le t<g$。对 $t$ 取平均，每一频率留存的比例都是 $\Delta$，故某个 $t$ 满足 $\|\Pi F^{(t)}\|_2^2\ge\Delta c^2$，此处用计数测度。置 $H=\Pi F^{(t)}$。因为 $\langle H,F^{(t)}\rangle=\|H\|_2^2$ 且 $F^{(t)}$ 支撑在 $I$，Cauchy 不等式给出 $\sum_{j\in I}|H_j|^2\ge\Delta^2c^2$。取 $H$ 的实部与虚部中在 $I$ 上平方和较大的一个作为 $Z_p$。由互素频率集合关于取负对称，所取实序列的傅里叶系数仍只支撑在阶 $g$ 的特征上，系数模仍至多 $c\sqrt2$。又有 $\|Z_p\|_\infty\le wc/\sqrt m\le c\sqrt x/D$，并且
$$
 M_4(p)\ge \frac{1}{gm}\left(\sum_{j\in I}|Z_p(j)|^2\right)^2
 \ge \frac{\Delta^4c^4}{4gm}.
$$
短关系项的绝对值总和仍至多 $g^{-2}x^{o(1)}$。它是上述下界的 $o(1)$ 倍，因为 $m/g\ll (\log D)^2D/x$ 而 $\theta<1$ 固定。对选定的素数求和，再代入 $w,\Delta$ 的界，得到 (7)。这个模型甚至保留了单个剩余类计数的非负性约束。由于 $g\le2D$，在 $y=x$ 时有 $x/g-\sqrt{x}Z_p(a)\ge(1-2c)x/g>0$。

从原数域的阶限制还能得到一个具体判断。置 $\alpha=7/208$。每个素数模数 $p\in(D,2D]$ 上，阶至多 $2^{1/64}x^\alpha$ 的特征最多 $x^{\alpha+o(1)}$ 个，因为每个给定阶至多有该阶那么多个特征，且阶必整除 $p-1$。乘积相等的四元组一旦指定其中一个位置和另外两个因子，第四个就已确定。因此含有至少一个这种低阶因子的四元组，对整个素数模数块上真实四阶和的总绝对贡献至多为 $x^{\alpha+o(1)}/D$。若 $D\le x^\theta$、$\theta<(2-\alpha)/3=409/624$，该量比 $x^2/D^4$ 小一个固定幂。特别在 $\theta=3/5$ 时，若只对素数模数求和也出现前篇意义下的大误差，相应四阶超额必须有一部分来自**四个因子都超出原数域阶上限**的组合。(5) 的模数范围则包括合数。

## 5. 有限核对

对模 $47$ 的全部 $44$ 个非实特征，在 $|\gamma|\le12$ 内按特征共计 $590$ 个零点。程序对每对共轭特征中的一个，用球算术沿上下水平边从 $\Re s=2$ 连续跟踪 $L$ 函数的辐角。每个子段的值域位于一个开半平面内，因此可累加端点比值的主辐角。由函数方程及 Gamma 因子得到整个临界带的计数，临界线上的变号区间数与之相等。另一个特征的纵坐标由共轭反射得到。

取数值宽度 $h=0.30$，只在此有限检验中把短关系阈值取为 $S=2$，保留阶大于 $2$ 的特征。下表将直接对剩余类求得的四阶矩与特征侧展开比较。两种算法的最大绝对差小于 $3\times10^{-18}$。

| 端点 $y$ | 完整四阶矩 $M_4$ | 去掉短关系后的 $H_{47}$ | $(\max_a\lvert Z\rvert)^2 M_2$ |
| --- | ---: | ---: | ---: |
| $256$ | $0.001433800$ | $+0.000242316$ | $0.002414748$ |
| $512$ | $0.001083679$ | $+0.000138192$ | $0.002291185$ |
| $1024$ | $0.001646381$ | $+0.000397608$ | $0.002695170$ |

另以有理数检验 $24\,576$ 组区间配置，含禁区的两个边界。命题 3 的构造取 $p=12289,m=64,c=1/8$，删去阶不超过 $4$ 的系数后，$gmM_4/c^4=0.99869804$。命题 4 的投影取 $p=1009,m=64$，保留 $288$ 个阶为 $1008$ 的频率，得到 $gmM_4/c^4=0.00767790$。回到特征侧，所有被排除频率的系数模均小于 $2.1\times10^{-16}$。这些有限计算检验了求和的符号和归一化。(5)–(7) 的渐近结论由正文证明给出。

复算使用 Python 3、`python-flint==0.9.0` 和同目录的 `verify.py`。运行 `python -B verify-common-windows.py` 可再现上表与两种傅里叶构造。`CommonWindow.lean` 在 Lean 4.33.1 与 Mathlib 4.33.1 下通过编译，两个公开定理的公理依赖均为 `propext`、`Classical.choice`、`Quot.sound`。它验证实变量共同选窗及窗口个数。从上同调输入到 `VAt` 的传递和四阶估计由正文给出证明。[数值复算](verify-common-windows.py) · [共用的零点隔离程序](verify.py) · [Lean 文件](CommonWindow.lean)
