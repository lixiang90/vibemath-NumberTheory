---
title: "同一端点的大导子辅助数域与 theta 测试向量的平移内积"
date: "2026-11-25"
abstract: "在预先给定的端点构造大导子圆分辅助域，在同一分裂素数处用局部 Weil 测试向量精确实现算术混合内积，并给出与矩阵秩及剩余次数无关的范数。原 theta 链在固定水平型上的子群平移只产生一个特征相位。对短区间中几乎每个共同平移，经典大筛无条件给出超过二分之一的模数指数。沿提升链推进完整 EH，仍缺全局匹配和高导子的相关估计。"
---

# 同一端点的大导子辅助数域与 theta 测试向量的平移内积

取 $Q=x^{3/5}$。在任意长为 $H\ge x^{1/5+\delta}$ 的整数区间里，对几乎每个 $b$，把同一个 $b$ 用于所有与它互素的 $q\le Q$，等差数列中素数的误差绝对值之和是 $O_A(x/\log^A x)$。每个模数可以另选最坏计数端点。这是经典均值定理与乘法大筛的无条件推论。要允许每个模数另选最坏**剩余类**，还需要控制误差在单位群上的集中。

Alpöge–Claude 的证明用 theta 提升的次数实现斯捷潘诺夫方法中重数所起的作用。有限处平移能读取哪些特征，还取决于向量在水平子群上的变换律。在原稿指定的水平型上，沿原紧开子群平移只给一个相位。本文在固定的计数端点构造圆分辅助域，再为同一模数的所有特征找一个共同分裂素数。在此处的局部 Weil 表示中，一个仅用矩阵左上角元素的测试函数能精确给出等差数列计数所需的全部平移内积。把这些向量配成全局 theta 链，并估计其提升后的权，是尚未解决的两步。[Alpöge 的说明，第 3–4 页](https://alpo.ge/remark.pdf) · [形式化证明的独立复核记录](https://math.mit.edu/~drew/grh/) · [前篇关于斜坡权与平移的研究](theta-chains-and-translated-cyclic-tests.md)

## 1. 在给定端点选择数域

固定 $C\ge0,\ 0<\eta<1,\ 0<\epsilon<1/52,\ 0<\kappa\le1$，先选
$$
x/(\log(2x))^C\le t\le x,\qquad P=t^\kappa,\qquad B=P^{\epsilon/5},\qquad Q\le x^{1-\eta}.
$$
以下断言对充分大的 $x$ 一致。$q$ 表示等差数列模数，$r\mid q$ 表示本原非主特征 $\chi$ 的导子，$o$ 表示其阶。为区别于 $q$，用 $u,v$ 表示源文件 `WinPair` 中的两个辅助素数。`IsChiField` 指定的是圆分扩张里的逐点固定子。`WinDatum` 还要求辅助素数的解析估计和一个算术上同调界。[固定子、素数窗口及上同调窗口的原定义](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/CharacterBoxes/ExistsPageN.lean#L19-L48) · [`WinDatum` 的各项要求](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/CharacterBoxes/TwistForm.lean#L9-L26)

**命题 1（共享辅助素数）。** 使用 Dirichlet GRH 的算术级数素数定理，可在 $t$ 已固定以后选出一个 $u>Q$ 及至多 $O_{\epsilon,\kappa}(\log(2x))$ 个 $v>Q$。每对 $(u,v)$ 都满足 `WinPair` 在 $P$ 的同余和大小条件。对每个 $q\le Q$，表中有一个 $v=v_q$ 使
$$
 \left(\frac{-v}{\ell}\right)=1\qquad(\ell\mid q,\ \ell\text{ 为奇素数}).
 \tag{1}
$$
同一个 $v_q$ 适用于这个 $q$ 的全部特征。

**证明。** 在 $[e^B/2,e^B]$ 中取素数 $u\equiv3\bmod8$，固定模数的素数定理保证其存在。置
$$
 M=16u\prod_{3\le\ell\le B,\ \ell\ {\rm prime}}\ell,\qquad Y=\exp(B\log^2 P).
$$
在 $[Y/2,Y]$ 中取满足 $v\equiv7\bmod16$ 且 $v\equiv-1\bmod(M/16)$ 的素数。它们均满足 `WinPair`，并对所有奇素数 $\ell\le B$ 满足 (1)。给定 $q$，把其大于 $B$ 的不同素因子之积记为 $D$，个数记为 $s$。有 $s\le\log Q/\log B=O_{\epsilon,\kappa}(1)$。因 $u>Q$，对每个 $\ell\mid D$ 可独立指定 $-v$ 为模 $\ell$ 的非零平方。中国剩余定理给出固定模 $M$ 类的 $\varphi(D)$ 个既约细分类，其中恰有 $2^{-s}\varphi(D)$ 个可用。

对 $m=MD$，GRH 给出每个既约类在 $[Y/2,Y]$ 中的素数个数为
$$
 \frac{\operatorname{Li}(Y)-\operatorname{Li}(Y/2)}{\varphi(m)}
       +O\!\left(\sqrt Y\log^2(Ym)\right).
$$
由于 $\log m=O(B+\log x)=o(\log Y)$，误差与主项之比一致趋于零。故在原模 $M$ 类中的素数上均匀抽样，每个 $q$ 成功的概率至少为同一个常数 $c_{\epsilon,\kappa}>0$。独立抽取 $k$ 次，某个 $q$ 未被覆盖的概率至多 $Q(1-c_{\epsilon,\kappa})^k$。选 $k=O_{\epsilon,\kappa}(\log(2x))$ 使它小于 $1$，去除重复项即得。这里使用 GRH 的地方是对增长模数 $MD$ 的素数计数。[E. S. Lee 的 GRH 一致显式公式及等差数列素数定理](https://arxiv.org/abs/2301.13457)

**命题 2（显式的圆分固定域）。** 写 $F_\chi=\mathbb Q(\zeta_r)^{\ker\chi}$，$K_e^+=\mathbb Q(\zeta_{2^{e+2}})^+$，$2^\nu\parallel r$，并约定 $\nu=0$ 时相应圆分域为 $\mathbb Q$。定义
$$
 d_\chi=[F_\chi\cap\mathbb Q(\zeta_{2^\nu})^+:\mathbb Q],\qquad
 c_\chi=o/d_\chi,\qquad
 e_\chi=\left\lceil\log_2(t^{2+8\epsilon}/c_\chi)\right\rceil .
$$
取 $e=e_\chi$ 和命题 1 的 $u,v_q$。下列数域满足 `IsChiField` 的固定子条件，
$$
 E_\chi=K_e^+F_\chi\,\mathbb Q(\sqrt{-u},\sqrt{-v_q}),\qquad
 [E_\chi:\mathbb Q]=4c_\chi 2^e,\qquad
 t^{2+8\epsilon}\le c_\chi 2^e<2t^{2+8\epsilon}. \tag{2}
$$
它是阿贝尔 CM 域。本原特征具有相同核时，这一构造给同一个域。对同一 $q$ 还可把所有 $E_\chi$ 同时嵌入一个域
$$
 K^+_{\lceil\log_2 t^{2+8\epsilon}\rceil}\mathbb Q(\zeta_q,\sqrt{-u},\sqrt{-v_q}).
 \tag{3}
$$

**证明。** 因 $c_\chi\le r\le x^{1-\eta}$，有 $2^e\ge t^{2+8\epsilon}/r\ge r$。特别地 $e+2\ge\nu$。圆分域的交公式给出 $\mathbb Q(\zeta_r)\cap K_e^+=\mathbb Q(\zeta_{2^\nu})^+$，故 $[K_e^+F_\chi:\mathbb Q]=2^e o/d_\chi$。两个负二次判别式分别是 $-u,-v_q$，二者互素且与 $\operatorname{lcm}(2^{e+2},r)$ 互素。圆分伽罗瓦群的中国剩余分解使复合域次数再乘 $4$。复共轭在这个阿贝尔扩张中非平凡，其不动域全实且指数为 $2$。在层 $N=\operatorname{lcm}(2^{e+2},r)uv_q$ 上，单位 $a$ 逐点固定复合域当且仅当
$$
 a\equiv\pm1\bmod 2^{e+2},\qquad \chi(a)=1,\qquad
 (a/u)=(a/v_q)=1.
$$
这正是所需固定子。取整给出 (2) 的窗口，$e_\chi\le\lceil\log_2 t^{2+8\epsilon}\rceil$ 给出 (3)。奇导子时 $d_\chi=1$。偶导子需要保留交次数，例如 $r=16$ 时能出现 $d_\chi>1$。原稿 `ordAny` 证明框架的阶最终具有 $c2^e$ 的形状；这里的圆分次数计算给出固定子构造中 $c_\chi$ 的取值及足够的稳定范围。[原稿关于阶的计算](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/CharacterBoxes/OrdAny.lean#L84-L105)

**进入原提升链还需要什么。** `KerChi` 要求 $r\mid\mathrm{QS}$。`FieldRows.small_le` 要求对每个奇素数 $\ell\mid\mathrm{QS}$ 有 $4\ell\le b$，该处使用 $b=\lfloor P^{\epsilon/5}\rfloor$。所以 $r$ 含大于 $B/4$ 的奇素因子时，(1) 即使在所有这些素数处成立，也无法满足这条大小约束。`FieldRows.o_pow` 另要求 $o^{64}\le2^e$；对 (2) 的奇导子参数，这要求 $o^{65}<2t^{2+8\epsilon}$。例如 $t=x$、素导子 $r\asymp x^{3/5}$、$\chi$ 忠实，(2) 仍给次数在 $[4t^{2+8\epsilon},8t^{2+8\epsilon})$ 的域，两条链上要求都失效。原构造还逐特征选上同调好窗口。上述代数构造解决了固定端点上的固定子与次数问题，完整的 `Host` 及共同上同调界仍待构造。[`FieldRows` 的原始假设](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/CoprimeRanks/SecondSetNoRelation.lean#L13-L23) · [`fieldChi` 额外使用的阈值](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/CharacterBoxes/FieldChi.lean#L20-L36) · [在好窗口上调用它的位置](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/ShortNorms/supplyOfAdmOf.lean#L83-L96)

## 2. 辅助数域中可计算的混合内积

记 $U_r=(\mathbb Z/r\mathbb Z)^\times$，$g_r=\varphi(r)$。在圆分域的各个复嵌入上取内积 $\langle f,h\rangle=g_r^{-1}\sum_{a\in U_r}\overline{f(a)}h(a)$，定义 $U_b f(a)=f(ba)$。令 $S$ 为 $\langle\chi\rangle$ 的全部生成元，它们都具有导子 $r$。对 $d\mid o$，以 $T_d$ 表示到 $F_\chi$ 中唯一次数 $d$ 子域的归一化相对迹，定义在 $\mathbb Q(\zeta_r)$ 上。置
$$
 P_S=\sum_{d\mid o}\mu(o/d)T_d,\qquad z_S=P_S\zeta_r,\qquad
 v_S=\sqrt{g_r/r}\,z_S .
$$
$z_S$ 属于 $F_\chi$，$g_rz_S$ 是代数整数。$v_S$ 的实数归一化在复嵌入空间中进行。

**命题 3（混合内积与其大小）。** 给 $b,n\in U_r$，令 $h$ 为根单位 $\chi(nb^{-1})$ 的阶。则
$$
 \langle U_bv_S,U_nv_S\rangle
 =\frac1{g_r}\sum_{\xi\in S}\overline{\xi(b)}\xi(n)
 =\frac{|S|}{g_r}\frac{\mu(h)}{\varphi(h)}. \tag{4}
$$
特别地，$h$ 含平方因子时内积精确为零。该内积用 $F_\chi$、$E_\chi$ 或 (3) 的全部嵌入计算都相等，均按次数归一化。只改变辅助素数 $u,v_q$，(4) 的值保持相同。

**证明。** $T_d$ 是对 $\ker(\chi^{o/d})$ 取平均。它保留 $\langle\chi\rangle$ 中阶整除 $d$ 的特征。莫比乌斯反演表明 $P_S$ 是到 $S$ 张成空间的正交投影。每个迹分母整除 $g_r$，给出整性。向量 $\zeta_r$ 在频率 $\xi\in S$ 上的系数是 $\tau(\overline\xi)/g_r$。本原 Gauss 和的模是 $\sqrt r$，故 $v_S$ 在每个频率的系数模方都是 $1/g_r$。$U_n$ 的本征值为 $\xi(n)$，遂得 (4) 的第一式。第二式是 Ramanujan 和公式，也可由对阶整除 $d$ 的所有特征求和再反演得到。子域的每个复嵌入向扩域延拓的次数相同，归一化内积因此与扩域无关。[Montgomery–Vaughan，第 19 章的本原 Gauss 和及乘法大筛](https://personal.science.psu.edu/rcv4/571s25/montgomery-vaughanII.pdf)

这些内积可以逐项还原素数计数。写 $\Lambda$ 为 von Mangoldt 函数，$\psi(y;q,b)=\sum_{n\le y,\ n\equiv b\bmod q}\Lambda(n)$。以 $\mathcal S(r)$ 表示导子 $r$ 的全部本原特征按上述生成元集合构成的分拆，令 $W_{r,S}^{[q]}(y)=\sum_{n\le y,\ (n,q)=1}\Lambda(n)U_n v_S$。对每个 $b\in U_q$，
$$
 \psi(y;q,b)-\frac{\sum_{n\le y,\ (n,q)=1}\Lambda(n)}{\varphi(q)}
 =\frac1{\varphi(q)}
   \sum_{\substack{r\mid q\\r>1}}g_r\sum_{S\in\mathcal S(r)}
       \langle U_b v_S,W_{r,S}^{[q]}(y)\rangle . \tag{5}
$$
这是 (4) 和模 $q$ 特征正交性的直接结果。限制 $(n,q)=1$ 精确处理了诱导到 $q$ 后删去的素数幂。在预先指定的 $y=t$，(3) 是这一等式中全部算术测试向量的共同嵌入空间。

## 3. 链上已有的平移公式能读出多少相位

原稿实际证明了一个单层混合 Petersson 公式。对满足该层度量假设的纯张量 $f$，它的形式为
$$
 {\rm Pet}(R(b_\infty k)f,f)
       =\lambda_{\rm lev}(k)\,{\rm ip}(u,u)
         \prod_\tau\langle v_\tau,\pi_\tau(b_\infty)v_\tau\rangle ,
       \qquad k\in K_{\rm lev}. \tag{6}
$$
这里 $b_\infty$ 是无穷处平移，$K_{\rm lev}$ 是有限处紧开子群，$\lambda_{\rm lev}$ 是该水平的特征。源码的 Petersson 约定对左变量线性，(4) 的约定对右变量线性，交换两变量只使特征共轭。源码还给出了 theta 积分对右平移的相容性，以及指定标准数据的 Petersson 自内积公式。[原链逐层的平移矩阵系数](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AbelianVarietyTheory/Proofs/DecoratedLatticeProofs/ThetaChainStagePeterssonValueRealLevelTranslate.lean#L13-L37) · [theta 积分与右平移](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/ThetaCorrespondence/ThetaPairLifts/LocalTypeDictionary/ThetaIntegralCompRightTranslate.lean#L13-L18) · [标准数据的 Rallis 范数公式](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/ThetaCorrespondence/ThetaPairLifts/LocalTypeDictionary/ThetaLiftPeterssonNormEqRallisIntegral.lean#L38-L47)

**命题 4（水平子群平移的精度下限）。** 固定 (6) 的无穷处分量。若将算术单位群经群同态 $U_r\to K_{\rm lev}$ 实现为有限平移，则该层在 $U_r$ 上只给出一个特征的标量倍数。若线性组合的各项共使用 $m$ 个不同的水平特征拉回，令所得函数为 $F$，并以 $K_S(a)=\langle v_S,U_av_S\rangle$ 为目标。则
$$
 \left(\frac1{g_r}\sum_{a\in U_r}|F(a)-K_S(a)|^2\right)^{1/2}
 \ge \frac{\sqrt{\max(|S|-m,0)}}{g_r}. \tag{7}
$$
这对固定形状的测试向量如何归一化都是必要条件。

**证明。** (6) 中只有 $\lambda_{\rm lev}(k)$ 依赖有限平移。每个拉回的水平特征贡献一个频率，而 (4) 在 $S$ 上的每个傅里叶系数均为 $1/g_r$。在 $F$ 遗漏的频率上应用 Parseval 即得 (7)。允许选择 $S$ 中任意 $m\le |S|$ 个频率并自由定标时能取等号。例如模 $17$ 的一个忠实特征有 $8$ 个生成元，$K_S(1)=1/2$、$K_S(-1)=-1/2$，其余 $14$ 个单位上的值全为零。单水平特征的最小均方误差为 $7/256$。

**分裂处的局部构造。** 可以在 (2) 的域上另选一个有限处，在真正的局部 Weil 表示中同时产生 (4) 的全部频率。在命题 1 的 GRH 输入下，对给定 $q$，这个有理素数 $w=w_q$ 可对所有 $r\mid q$ 共用。令 $m=\exp U_q$（单位群的指数）、$L=\log(2x)$，置
$$
 Z=(Q+1)^2(B+L)^4 L^4.
$$
在 $[Z,2Z]$ 中可选 $w\equiv1\bmod m$ 且 $(-u/w)=1$。为证此点，在 $n\equiv1\bmod m$ 上以 $(1+\chi_{-u}(n))/2$ 加权求和 $\Lambda(n)$，其中 $\chi_{-u}$ 为判别式 $-u$ 的二次特征。因 $u>Q\ge m$，模 $m$ 的每个特征与 $\chi_{-u}$ 的乘积都非主。特征正交性和 GRH 使区间主项为 $Z/(2\varphi(m))$，误差为 $O(\sqrt Z\log^2(2Zmu))$。删去高次素数幂的代价为 $O(\sqrt Z\log^2(2Z))$。有 $\sqrt Z\ge m(B+L)^2 L^2$ 及 $\log(2Zmu)\ll B+L$，主项占优。又有 $2Z<u<v_q$，故 $w$ 与 (2) 中各个导子互素。若 $F=E_\chi^+$，则 $E_\chi=F(\sqrt{-u})$。$-u$ 已在 $\mathbb Q_w$ 中成为平方，故 $F$ 在 $w$ 上的**每个**有限处都在 $E_\chi/F$ 中分裂。此估计控制有理素数 $w$，其在 $F$ 中的剩余次数可以较大。

取这样一个有限处 $\mathfrak p$，以 $k$ 表示其剩余域。在分裂酉对的局部 Schrödinger 模型中，测试函数属于矩阵空间上的 Schwartz–Bruhat 空间，目标群的对角单位通过矩阵乘法作用。取两个任意正秩，把整数矩阵格的 Haar 体积归一为 $1$。由 $o\mid m\mid w-1$，可将 $\chi(U_r)$ 识别为 $\mathbb F_w^\times$ 中的 $o$ 阶子群 $C_o$，并以 Teichmüller 提升将它放入目标群的第一对角坐标。若选择分裂时带有酉特征因子，下述不交支撑及内积结论依然成立。[Moss–Trias 的 II 型局部 Weil 矩阵模型](https://arxiv.org/abs/2312.12031) · [原稿分裂处模型所用的额外匹配条件](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AbelianVarietyTheory/Proofs/DecoratedLatticeProofs/IsThetaFunSplitPlaceSchroedingerModelStatement.lean#L10-L17)

从 $\mathbb F_w^\times/C_o$ 的每个陪集选一个代表，组成 $T$。令 $\Phi$ 为下列集合的 $L^2$ 归一化示性函数，
$$
 \left\{X\text{ 为整数矩阵}:
       \operatorname{Tr}_{k/\mathbb F_w}(\overline X_{11})\in T\right\}.
$$
剩余域迹满射，各纤维等大，故该集合的体积为 $(w-1)/(wo)$。对 $h\in C_o$ 的平移将迹条件变为落在 $hT$ 中（相反的作用约定给 $h^{-1}T$）。这些集合两两不交。因此直接得到
$$
 \langle\omega(h)\Phi,\omega(h')\Phi\rangle={\bf1}_{h=h'},\qquad
 \|\Phi\|_\infty=\sqrt{\frac{wo}{w-1}}. \tag{7a}
$$
这两个量都与矩阵秩及剩余次数 $[k:\mathbb F_w]$ 无关。任何由 $o$ 个格内不交平移生成的单位向量都有 $\|\Phi\|_\infty\ge\sqrt o$，因为每个支撑的体积至多 $1/o$。

记 $P^{\rm loc}$ 为到 $C_o$ 的忠实特征之和的正交投影，$\Phi_S=\sqrt{o/g_r}\,P^{\rm loc}\Phi$，记 $\iota(a)$ 为 $\chi(a)$ 所确定的对角单位。(7a) 表明每个特征在 $\Phi$ 中的投影范数平方恰为 $1/o$；这也可从投影算子 $o^{-1}\sum_{h\in C_o}\overline{\alpha(h)}\,\omega(h)$ 直接算出。于是
$$
 \langle\omega(\iota(b))\Phi_S,\omega(\iota(n))\Phi_S\rangle
   =\frac1{g_r}\sum_{\xi\in S}\overline{\xi(b)}\xi(n). \tag{7b}
$$
这给出了从 (4) 的平移张成空间到局部 Weil 空间的等距映射，因为两族向量的 Gram 矩阵相同。(5) 中的混合内积可以逐项换成这些局部内积。原稿还定义了按水平提升选取代表的有限商作用 `actG0`，有代表及指定变换型时证明了该作用的公式。这给出了比较上述局部向量与全局测试数据的具体位置。[水平变换律、有限商平移与作用公式](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Formalism/LocalToGlobal/IndefiniteReal/SubgroupCharExtensionOneApply.lean#L65-L115) · [存在水平提升时的作用公式](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Formalism/LocalToGlobal/IndefiniteReal/SubgroupCharExtensionOneApply.lean#L223-L235)

**局部向量经过全局提升后还缺的估计。** 要把 (7b) 用于全局 theta 链，须先证明这处指定的局部向量能与其余位置的数据、整性及水平条件同时匹配。设匹配并提升后从 $S$ 的频率空间到实际链的复合映射为 $J$。若它与 $U_r$ 的酉作用相容，则 $J^*J$ 在各特征线上作用为 $\lambda_\xi=\|J\xi\|^2$（特征取单位范数）。此时提升后的 (7b) 为 $g_r^{-1}\sum_{\xi\in S}\lambda_\xi\overline{\xi(b)}\xi(n)$。对 $n$ 作 (5) 的加权和，每个频率都乘上 $\lambda_\xi$，其中一些也可能为零。全局的非消失、各权的大小以及加权和在模数平均中的抵消都需要新的输入。

若只考察一个生成元集合 $S$，可精确度量上述权的误差。写 $Z_\xi=\sum_{n\le y,\ (n,q)=1}\Lambda(n)\xi(n)$，以 $\Delta_b$ 表示 (5) 中这个内积在提升后与提升前的差（不含外侧的 $g_r/\varphi(q)$）。在上述相容映射存在的假设下，再用一次 Parseval 得
$$
 \frac1{g_r}\sum_{b\in U_r}|\Delta_b|^2
 =\frac1{g_r^2}\sum_{\xi\in S}|\lambda_\xi-1|^2|Z_\xi|^2 .
$$
在局部构造中也可以只投影到一个特征 $\alpha$，再乘 $\sqrt{o/g_r}$。令 $\xi=\alpha\circ\iota$（把对角群认作 $C_o$），它的混合内积为 $\overline{\xi(b)}\xi(n)/g_r$。由此可单独指定四特征相关式中的每一个频率。

## 4. 对所有模数和几乎每个共同平移的估计

以下直接对算术内积求平均，证明中的常数独立于几何匹配问题。

定义
$$
 \mathcal E_b(x,Q)=\sum_{\substack{q\le Q\\(q,b)=1}}
    \sup_{0\le y\le x}\left|\psi(y;q,b)-\frac{y}{\varphi(q)}\right|,
 \qquad L=\log(2x).
$$
$y$ 可以随 $q$、$b$ 变化。取上确界也涵盖计数函数跳跃处的左极限。

**定理 5（短区间中的共同平移）。** 固定 $A>0,\ 0<\eta<1$，$1\le Q\le x^{1-\eta}$。对含 $H\ge1$ 个连续整数的任意区间 $I$，无条件地有
$$
 \frac{\#\{b\in I:\mathcal E_b(x,Q)>xL^{-A}\}}{H}
 \ll_{A,\eta} L^{C_A}\left(\frac Qx+\frac{Q^2}{xH}\right). \tag{8}
$$
$C_A$ 只依赖 $A$。取 $Q=x^\theta$，任意 $1/2<\theta<1$ 及 $H\ge x^{2\theta-1+\delta}$、$\delta>0$，右端趋于零，且对区间位置一致。这个推论使用 Vaughan 型特征均值估计和乘法大筛。Maynard 关于固定剩余类与具有适当大小因子的模数的结果采用另一组量词。[Drappeau–Fiorilli，引理 2.1](https://doi.org/10.1112/tlm3.12030) · [Montgomery–Vaughan，定理 19.16](https://personal.science.psu.edu/rcv4/571s25/montgomery-vaughanII.pdf) · [Maynard 关于固定剩余类的定理](https://arxiv.org/abs/2006.06572)

**证明。** 写 $\Psi_\xi(y)=\sum_{n\le y}\Lambda(n)\xi(n)$，先用本原特征和替代模 $q$ 的诱导特征和。对每个特征，这一替换的误差至多 $\omega(q)L$，因为每个 $p\mid q$ 的所有幂合计至多 $L$。在傅里叶反演的 $1/\varphi(q)$ 归一化下，对模数求和的误差是 $O(QL^2)$。主特征均值与 $y/\varphi(q)$ 之间的总误差为 $O_N(xL^{-N}+QL^2)$，由带任意对数节省的素数定理、$\sum_{q\le Q}1/\varphi(q)\ll L$ 给出。这两项对所有 $b$、所有端点一致。

取 $D=\sqrt x/L^B$，其中 $B=B(A)$ 充分大。唯一导子 $r=1$ 的项已经计入主项。其余 $r\le D$ 的本原特征满足
$$
 \sum_{1<r\le D}\frac1{\varphi(r)}
       \sum_{\xi\bmod r}^{*}\sup_{y\le x}|\Psi_\xi(y)|
       \ll_A xL^{-A-3}. \tag{9}
$$
事实上，对 $r\le L^K$ 用 Siegel–Walfisz，并取足够多的对数节省。该估计对上确界同样成立，$y\le\sqrt x$ 时直接用平凡界。对 $L^K<r\le D$，Drappeau–Fiorilli 引理 2.1 给出 $L^{O(1)}(x/L^K+D\sqrt x+x^{5/6})$ 的上界，取足够大的 $K,B$ 即可。又因 $\varphi(rm)\ge\varphi(r)\varphi(m)$，有
$$
 \sum_{\substack{q\le Q\\r\mid q}}\frac1{\varphi(q)}\ll\frac L{\varphi(r)}.
$$
故 (9) 的特征在全部 $q\le Q$ 中诱导出的总贡献一致为 $O_A(xL^{-A-2})$。

保留 $r>D$。对任何整数区间 $J\subseteq[1,\lfloor x\rfloor]$ 定义
$$
 T_q(a,J)=\frac1{\varphi(q)}
 \sum_{\substack{r\mid q\\r>D}}\sum_{\xi\bmod r}^{*}\overline{\xi(a)}
                       \sum_{n\in J}\Lambda(n)\xi(n),\qquad
 M_q(a)=\max_{m\le x,\ m\in\mathbb N}|T_q(a,[1,m])|.
$$
对每个二进区间 $R<q\le2R$，有
$$
 \sum_{R<q\le2R}\sum_{a\in U_q}M_q(a)^2\ll_A x R L^{C'_A}. \tag{10}
$$
为证此式，将 $[1,\lfloor x\rfloor]$ 补零到二次幂长度，每个前缀是 $O(L)$ 个标准二进子区间的不交并。对端点取最大后用 Cauchy，$M_q(a)^2\ll L\sum_J|T_q(a,J)|^2$，其中 $J$ 遍历全部标准区间。固定 $J$，模 $q$ 的字符正交性把对 $a$ 的平方和化为 $1/\varphi(q)$ 乘以对应本原特征和的模方之和。再对 $q$ 求和并固定导子段 $S<r\le2S$（可从 $S=D$ 起逐次加倍）。用上式的诱导权上界及 $r>S$，这段至多
$$
 \frac{L}{S}\sum_{S<r\le2S}\frac r{\varphi(r)}
   \sum_{\xi\bmod r}^{*}\left|\sum_{n\in J}\Lambda(n)\xi(n)\right|^2
 \ll L\left(\frac{x}{S}+S\right)\sum_{n\in J}\Lambda(n)^2 .
$$
最后一步是乘法大筛。这里 $S\ge D$，因而 $x/S\le S L^{2B}$。对 $S<2R$ 加倍求和只产生 $O(R)$ 的几何和。每个给定长度的二进区间彼此不交，$\sum_{n\le x}\Lambda(n)^2\le xL^2$。对 $O(L)$ 个长度求和即得 (10)。当 $2R\le D$ 时左边为零。

每个模 $q$ 的类在 $I$ 中至多出现 $H/q+1$ 次。应用 Cauchy，再按 $q$ 加倍分段并用 (10)，得到
$$
 \frac1H\sum_{b\in I}\left(\sum_{\substack{q\le Q\\(q,b)=1}}M_q(b)\right)^2
 \le Q\sum_{q\le Q}(q^{-1}+H^{-1})\sum_{a\in U_q}M_q(a)^2
 \ll_A x Q L^{C''_A}(1+Q/H).
$$
对充分大的 $x$，先前剥离的总误差小于 $xL^{-A}/2$。对余项使用 Chebyshev 不等式得到 (8)。证毕。

## 5. 哪一种抵消还可以提高完整的分布指数

**二阶矩的集中例。** 固定 $1/2<\theta<1$，取 $R=x^\theta$、$1\le H\le R^2/x$、$I=\{1,\ldots,H\}$。对每个素数 $R<q\le2R$，在 $U_q$ 上先令 $f_q(a)=x\varepsilon_a/(4R)$（$a\in I$），其余为零，然后减去 $U_q$ 上的均值；$\varepsilon_a\in\{-1,1\}$。恒有 $\sum_{a\in U_q}|f_q(a)|^2\le H(x/(4R))^2\le x/16$。取共同的一列独立随机符号，对复特征和的实、虚部应用有界独立变量的指数尾估计，再对至多 $O(R^2)$ 个模数与特征的组合取并，存在选取使所有非主频率都满足
$$
 \left|\sum_{a\in U_q} f_q(a)\chi(a)\right|
 \ll (x/R)\sqrt{H\log(2R)}\ll\sqrt{x\log(2R)} .
$$
另一方面，$H/R\to0$，故每个 $b\in I$ 和每个上述素数 $q$ 都有 $|f_q(b)|\ge x/(8R)$。仅对这些素数求和已有 $\sum_q|f_q(b)|\gg x/\log R$。甚至可在每个 $f_q$ 上加回 $x/(q-1)$ 得到非负数组。这个构造说明，即使附上逐特征平方根界，单靠上述二阶矩信息也允许整个长度 $R^2/x$ 的区间同时集中。这些数组只用作估计方法的反例。

**完整 EH 的一个足够条件。** 沿用定理 5 的高导子余项。设在 $D/2<R\le Q$ 的所有二进模数段中能证明
$$
 \sum_{R<q\le2R}\sum_{a\in U_q}M_q(a)^4
       \ll_A x^{2+\sigma} L^{C_A},\qquad \sigma\ge0. \tag{11}
$$
每个 $q$ 对剩余类取最大，Hölder 仍给出该段总和至多 $O_A(R^{3/4}x^{1/2+\sigma/4}L^{C_A/4})$。结合 (9)，对 $Q=x^\theta$，只要 $3\theta+\sigma<2$，任意对数节省都成立。特别地，在 $\theta=3/5$ 处，$\sigma<1/5$ 已足够。

可以在 (11) 的左端用 $M_q(a)^4\ll L^3\sum_J |T_q(a,J)|^4$，把端点上确界化为共同的二进区间。若 $A_\xi(J)$ 是模 $q$ 的特征由导子 $>D$ 的本原特征给出的区间和，其余特征的系数取零，则精确地有
$$
 \sum_{a\in U_q}|T_q(a,J)|^4
 =\frac1{\varphi(q)^3}\!
 \sum_{\xi_1\xi_2=\xi_3\xi_4}
 A_{\xi_1}(J)A_{\xi_2}(J)\overline{A_{\xi_3}(J)A_{\xi_4}(J)} . \tag{12}
$$
(12) 保留了有乘积关系的不同特征。在显式公式中，它们对应同一区间端点处的四个零点和。要用 theta 链改善完整 EH，需要实现能读取这些频率的向量，或直接给这一相关和的模数平均建立足够的界。本文得到的算术迹恒等式、单水平特征的下限和共同平移定理都未给出 (11)。[前篇关于共同端点处四特征零点和的分析](common-windows-and-four-character-bounds.md)

若要在共同尺度上直接观测零点，可将 (12) 中的 $\mathbf1_J$ 换成 $f(n/t)$，其中 $f\in C_c^\infty(0,\infty)$。设 $\widehat f(s)=\int_0^\infty f(v)v^{s-1}\,dv$。对本原非主特征 $\xi$、$t\ge1$，Mellin 反演给出
$$
 \sum_n\Lambda(n)\xi(n)f(n/t)
 =-\sum_{\rho_\xi}t^{\rho_\xi}\widehat f(\rho_\xi)
   -{\bf1}_{\xi(-1)=1}\widehat f(0)
   +O_f\!\left(t^{-1/2}\log(2r)\right). \tag{13}
$$
这里 $\rho_\xi$ 取遍非平凡零点并计重数。将 $-L'/L$ 的积分线从 $\Re s=2$ 移到 $-1/2$，恰有偶特征在零点 $s=0$ 的额外留数；新积分线上的对数导数由函数方程及 $\Re(1-s)=3/2$ 上的绝对收敛控制。$\widehat f$ 的快速衰减保证零点和绝对收敛。把 (13) 代入平滑版 (12)，每个零点和中的 $t$ 完全相同。把这种相关和从平滑权传回 (11) 的各个区间还需要一致的端点控制。

## 6. 可复算的有限检验

程序从根单位出发，按子域的相对迹构造 $z_S$，直接取嵌入内积。它用 Ramanujan 和与特征和分别核对 (4)，并对模 $5,8,12,15,16,17,25,35,45,56,63,65$ 的每个既约类，在端点 $217,512$ 将 (5) 与直接枚举的素数幂计数核对。共 $440$ 个差值的复数球全部含零，模的上端均小于 $1.2\times10^{-53}$。程序还枚举圆分伽罗瓦群的像检验 (2) 中的交次数修正，并在模 $17$ 上直接核对 (7) 的 $7/256$。

同余覆盖的有限算例取 $u=131$ 及小素数 $3,5,7$。$v=24318839$ 和 $16836119$ 两个素数已覆盖所有 $q\le105$ 的条件 (1)。在 $r=17,w=113$ 的二次剩余域扩张算例里，$16$ 个平移的支撑各占 $7/113$，投影后的每个混合内积都与 (4) 吻合；另核验了阶 $4,12,20$ 的局部向量。程序还在 $q=15,28,35,45$、$D=5$、$J=(73,217]$ 上分别从剩余类和与四特征乘积关系计算 (12)。这些有限算例检验同余、归一化和恒等式，渐近结论由正文证明。用 Python 3 和 `python-flint==0.9.0` 运行 `python -B verify-cyclotomic-mixed-periods.py` 即可复核。[复算程序](verify-cyclotomic-mixed-periods.py)
