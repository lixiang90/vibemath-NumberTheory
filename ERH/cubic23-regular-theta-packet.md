---
title: "判别式 −23 的三条 CM theta 特征线与提升链的局部输入"
date: "2026-11-25"
abstract: "构造同级同权的三条正则 CM theta 特征线。它们的 Frobenius 幂迹差的平方和，在剩余次数与 3 互素时精确检出三循环素数。给出三个混合 Rankin–Selberg 因子的消去公式，并证明固定 CM 域的分裂资料无法替代这三条特征线。"
---

# 判别式 −23 的三条 CM theta 特征线与提升链的局部输入

设 $K=\mathbb Q(\alpha)$，$\alpha^3-\alpha-1=0$，$M$ 为其 $S_3$ 分裂域，$k=\mathbb Q(\sqrt{-23})$。$M/k$ 是 Hilbert 类域，$\mathrm{Cl}(k)=C_3$。固定一个非平凡理想类特征 $\psi$，记 $\chi=\chi_{-23}$，$L=L_k(s,\psi)=\zeta_K(s)/\zeta(s)$。记 $C$ 为 Frobenius 是三循环的有理素数集，$F$ 为完全分裂素数集，以下均略去 23。$\chi(p)=1$ 时，两种情况分别是 $p$ 上方的素理想属于非主类、主类。[域、二次型与 Artin 因子](cubic23-zeros.md)

原均方将问题归结为 $C$ 的素数权之和。这里构造同一权空间中的三个 Hecke 本征形式。即使先把 Frobenius 平方，它们的局部迹仍给出 $18\,\mathbf1_C(p)$。单纯在原来的惰性位置乘一个三次单位根，会丢失这条信息。[原均方的精确对角项](cubic23-class-character-mean-square.md)

## 在原提升链的哪一步加入类特征

Alpöge 的说明把长度约 $\sqrt X$ 的一串 theta 提升与 Stepanov 方法的重数联系起来。源码中可定位以下不同的命题。这一区分决定了可以替换哪项输入。[作者的说明，第 1 页](https://alpo.ge/remark.pdf) · [上游结论的独立核验](https://math.mit.edu/~drew/grh/)

| 源码中的接口 | 对 −23 推广的意义 |
| --- | --- |
| `TLiftK` 在 CM 域 $E$ 与其最大全实子域上规定分裂特征的限制、无穷型和权的不等式，然后给出提升存在性。 | 在 $k/\mathbb Q$，$\psi$ 在来自 $\mathbb Q$ 的 idèle 上平凡，无穷型为零。把已有的合适特征乘 $\psi^j$ 保留这两项条件。完整提升还须满足命题中的表示与权假设。[提升的定义](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/Foundations/GaloisRepresentations/ThetaParseval/HalfIndexLattices/TLiftK.lean#L12) |
| `chainPar_twist_chT` 将下一步的参数写成最初参数的扭转加上若干段参数。`StepExactNS` 在同一步上计算两个提升后向量的内积，还要求一个分裂位置满足 `NSAt`、指定的标准 $L$ 值及局部积分。 | 跟踪参数可保留特征下标。若原向量正交，其内积乘同一个标量仍为零。不同本征线的两个 Hecke 系数之差，须先经共同的线性泛函读出，再作平方。[参数递推](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AbelianVarietyTheory/Proofs/DecoratedLatticeProofs/ThetaLiftChainHeckeParameterTwistNextTargetCharStatement.lean#L10) · [内积公式的完整条件](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AbelianVarietyTheory/Definitions/LatticeFrames/StepExactNS.lean#L89) · [导出伴随标量的条件](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AbelianVarietyTheory/Proofs/DecoratedLatticeProofs/TwistedThetaChainStepAdjointScalarStatement.lean#L20) |
| 一项深度估计要求奇素数、未分歧以及 `InertInK`。其中 `InertInK` 指 $-H.d$ 在 $\mathbb F_p$ 中非平方。另一项在 $k_0\le mp$ 下控制深度的振幅。 | 这是对指定辅助虚二次域的惰性条件。区间宽度和振幅界可用于比较局部深度，应用时还须知道所比较的整数在目标类上有多少整除性。[惰性条件](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AbelianVarietyTheory/BlockTheory/DivisorBase/TypedPartSplitPrimes.lean#L40) · [因子深度区间](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AbelianVarietyTheory/Proofs/DivisorProofs/FactorDepthBoundedOddInertPrimeCmTypeStatement.lean#L10) · [振幅界](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AbelianVarietyTheory/Proofs/DivisorProofs/DepthOscillationLargeInertPrimeStatement.lean#L10) |

链的满射性另有有限维性、特征向量张成和参数独立性等假设。最终组装还用到 theta 积分与局部匹配、余体积估计及窗口输入，所输出的 `EveryCharacter` 量化的是有理数域上的 Dirichlet 特征。下文给出新的局部参数和矩阵元；将它们实现于整条链，并满足混合矩阵元的度量条件，仍需证明。[满射命题](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AbelianVarietyTheory/Proofs/DecoratedLatticeProofs/ThetaLiftChainSurjectiveOfAdjointNonvanishingStatement.lean#L10) · [最终组装](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/ThetaCorrespondence/ThetaPairLifts/ThirdPlaceSeries/EveryCharacterOfLiteratureOfVWin.lean#L810) · [EveryCharacter](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/CharacterBoxes/EveryCharacter.lean#L13)

## 惰性位置的限制有多强

若反分圆的奇数阶特征 $\phi$ 满足 $\phi^c=\phi^{-1}$，则对每个未分歧且由复共轭 $c$ 固定的素理想 $\mathfrak P$ 有 $\phi(\mathfrak P)^2=1$，故 $\phi(\mathfrak P)=1$。特别地，若取 $k$ 本身为辅助虚二次域，惰性 $p$ 的素理想就是主理想 $(p)$。在这样的地方，三阶扭转给三条特征线完全相同的值。即便在其他地方，用单位根乘单个非零代数整数也保持其所有有限赋值。

还有一项关于固定辅助域的定量限制。令 $E/\mathbb Q$ 是有限个 CM 域和全实域的正规闭包的合成域，$G=\operatorname{Gal}(E/\mathbb Q)$。设 $a(p)$ 是 $S_3$ 标准表示之迹，在完全分裂、换位、三循环上分别为 $2,0,-1$。则对任意复值类函数 $v$ 有

$$
 \lim_{x\to\infty}\frac{1}{\pi(x)}
 \sum_{\substack{p\le x\\ p\ {\rm unramified\ in}\ EM}}
 |a(p)-v(\operatorname{Frob}_{p,E})|^2
 =1+\frac1{|G|}\sum_{h\in G}|v(h)|^2. \tag{1}
$$

**证明。** CM 域每个嵌入都把自身的复共轭送到像上的复共轭。于是复共轭在其正规闭包中居中，在全实正规闭包上为恒等，故在 $G$ 中也居中。$E\cap M$ 在 $\mathbb Q$ 上正规，只可能是 $\mathbb Q,k,M$。最后一种会使复共轭在商群 $S_3$ 中居中，但它在那里是一个换位。若交为 $\mathbb Q$，在每个 $h$ 上方标准迹均值为零。若交为 $k$，偶陪集上的三个迹为 $2,-1,-1$，奇陪集上的三个迹均为零。把 Chebotarev 密度定理应用到纤维积 $\operatorname{Gal}(EM/\mathbb Q)$，交叉项为零，$a$ 的全群平方均值为 $(4+2)/6=1$，即得 (1)。这只限制**固定**的辅助分裂资料。$k$ 的 Hilbert 类域 $M$ 本身就有非中心的复共轭，理想类特征能读取其中的额外信息。

## 三条正则 theta 特征线

写 $\omega=(1+\sqrt{-23})/2$。对 $k$ 的非零分式理想 $\mathfrak a$，选取 $\mathfrak a^3=(\delta)$，生成元只差符号。取 $m\ge1$，$d=6m$，$j=0,1,2$，定义

$$
 h_{j,m}(\mathfrak a)=\psi(\mathfrak a)^j\delta^{2m},\qquad
 f_{j,m}(z)=\sum_{\mathfrak a\subseteq\mathcal O_k,\ \mathfrak a\ne0}
 h_{j,m}(\mathfrak a)e^{2\pi i zN\mathfrak a}. \tag{2}
$$

偶次幂消除了符号。$h_{j,m}$ 乘法，在 $(x)$ 上等于 $x^d$，因此有限导子为 1。Hecke–Shimura 定理给出级 23、权 $d+1$、特征 $\chi$ 的归一化 CM 尖点新形式。三个诱导参数互异，因为底层特征相差 $\psi^j$，而共轭把无穷型 $(d,0)$ 换为 $(0,d)$。令 $\nu(\mathfrak a)=\delta/\bar\delta$，幺正化后的特征为 $\xi_j=\nu^m\psi^j$。[Hecke 特征与新形式的定理 1.4](https://arxiv.org/pdf/math/0511228)

令 $t_{j,r}(p)$ 是 $f_{j,m}$ 的两个幺正 Satake 根的 $r$ 次幂之和。以下恒等式对每个 $r\ge1$ 成立，

$$
 \boxed{\displaystyle
 \sum_{0\le i<j\le2}\big(t_{i,r}(p)-t_{j,r}(p)\big)^2
   =18\,\mathbf1_C(p)\,\mathbf1_{3\nmid r}.} \tag{3}
$$

**证明。** 若 $p$ 在 $k$ 分裂，选 $\mathfrak p\mid p$，设 $w=\psi(\mathfrak p)$、$z=\nu(\mathfrak p)^m$。则 $t_{j,r}=z^rw^{jr}+\bar z^rw^{-jr}$。若 $w^r=1$，三个实数相同；否则其和为零，平方和为 $6$，由 $\sum_{i<j}(t_i-t_j)^2=3\sum t_i^2-(\sum t_i)^2$ 得 $18$。若 $p$ 惰性，第一幂迹为零、幺正行列式为 $-1$，所有幂迹与 $j$ 无关。证毕。

特别地，向一个独立的辅助虚二次域限制，惰性位置的剩余次数为 2，Frobenius 取平方后仍能区分 $C$。例如 $13$ 和 $59$ 都在 $\mathbb Q(\sqrt{-19})$ 惰性，却分别属于 $C$ 和 $F$，式 (3) 在两点分别为 $18$ 与 $0$。若改用剩余次数 3，三个通道的差全部消失。对任意个正则权作总质量为 1 的加权平均，(3) 的右端也完全相同。因而这种权平均自身无法改善前稿的平方对角主项。

写 $a_j(p)$ 为 (2) 的普通系数，则 $p^{d/2}t_{j,1}=a_j(p)$。$p^d t_{j,2}=a_j(p)^2-2\chi(p)p^d$，而 Fourier 系数 $a_j(p^2)=a_j(p)^2-\chi(p)p^d$。三个形式都满足 $a_1(f_{j,m})=1$，这里 $a_1$ 指取第一 Fourier 系数的线性泛函。故 $a_1(T_p f_{i,m}-T_p f_{j,m})$ 就是可作平方的标量 $a_i(p)-a_j(p)$。这一步给出共同读数；沿提升后度量控制这个泛函，需要额外的估计。

## 用素数 2 区分三条线所需的分母

固定 $2$ 上方由 $\delta=1+\omega$ 生成其三次幂的素理想，调整 $\psi$ 的选向一次。令 $s=\sqrt{69}$。权 7 和权 13 的 $T_2$ 特征值三元组依次为

$$
 x=(-7,(7-3s)/2,(7+3s)/2),\qquad
 y=(-79,(79+21s)/2,(79-21s)/2).
$$

事实上 $\delta^2=-5+3\omega$，$\delta^4=-29-21\omega$，乘三次单位根并与共轭相加即得。以 $\mathbf1,x,x^2$ 为列的行列式绝对值为 $135s$，以 $\mathbf1,x,y$ 为列时为 $576s$，其中 $x^2$ 逐项平方。由于 $\gcd(135,576)=9$，这些列在 $\mathcal O_{\mathbb Q(\sqrt{69})}[1/(3\cdot23)]$ 上生成整个三维自由模。证明是分别将两个矩阵乘伴随矩阵，再对 $135,576$ 用 Bézout 恒等式；$s$ 在该环中可逆。这给三个已由 $\psi$ 标记的坐标统一的分母控制。两种权仍在各自的模形式空间中。[整数复算](verify_theta23_stepanov.py)

## 一条在整数格上闭合的 Hecke 算子链

固定权 7，令 $S=\langle f_{0,1},f_{1,1},f_{2,1}\rangle$，$w=(f_{0,1}+f_{1,1}+f_{2,1})/3$，$\ell(f)=a_1(f)$。平均后 $w$ 等于主类的调和 theta 级数 $\Theta_6$。因 $\ell(w)=1$，$P(x)=\ell(x)w$ 是秩一投影。令 $b_p=\ell(T_p w)=a_p(\Theta_6)$。对未分歧 $p$，有如下**算子恒等式**，
$$
 P T_p P=b_pP,\qquad
 P\big(T_p+p^n(1+T)\operatorname{Id}\big)P
   =\big(b_p+p^n(1+T)\big)P . \tag{4}
$$
因此在一串不同的素数处交替作用 Hecke 算子和同一个 $P$，乘积恰为各处右侧标量之积乘 $P$。这从实际模形式空间中给出了后篇的整数多项式。若略去步间的 $P$，系数 $a_{p_1p_2}(w)$ 会混入不同理想类的配对。

这条链可以在每一个有限处取整，包括 3。令 $R=\mathcal O_{\mathbb Q(\sqrt{69})}$。三个形式的 Fourier 系数都是 $R$ 中的代数整数，因为这些系数是 $k(\sqrt{-3})$ 中的代数整数且由复共轭固定。在它们的系数域线性张成空间中，取所有 Fourier 系数都在 $R$ 的格 $\mathcal L$。三个形式在格中，而有限选取三个 Fourier 系数就能区分三维空间中的向量，因此 $\mathcal L$ 嵌入 $R^3$，是满秩有限生成格。$w$ 的系数是有理整数，$\ell(\mathcal L)\subset R$，故 $P$ 保持 $\mathcal L$ 且 $Rw$ 为直和项。好素数处 Hecke 公式 $a_t(T_p f)=a_{pt}(f)+\chi(p)p^6a_{t/p}(f)$ 也保持 $\mathcal L$，约定非整数下标的系数为零。这里在格中没有逐次投影造成的分母增长，尽管 $P$ 在三个本征向量的基下含有 $1/3$。

在阿基米德处，令 $N_i=\|f_{i,1}\|_{\rm Pet}^2$。由于不同的归一化新形式 Petersson 正交，投影的算子范数为
$$
 \|P\|_{\rm Pet}^2=\|w\|^2\|\ell\|^2
 =\frac19\left(\sum_iN_i\right)\left(\sum_iN_i^{-1}\right). \tag{5}
$$
这是与素数、深度和步数无关的一个固定常数。含 $m$ 个 Hecke 因子的乘积逐次估计时，$m+1$ 次投影额外贡献的对数至多 $(m+1)\log\|P\|$。后篇给出这个开销进入高斯窗口后的大小。这条链在经典 CM 模形式空间内闭合。原酉群提升链中相应整截面的构造及上游余体积估计如何作用于它，仍然是另一项算术问题。[分组后的 Hasse 导数与整体高度](cubic23-cyclotomic-jets.md)

## 哪一个混合矩包含目标 L 函数

记 $\pi_j=\operatorname{Ind}_{k}^{\mathbb Q}\xi_j$ 为幺正参数，上标 $(23)$ 表示略去 23 的 Euler 因子。在 $\Re s>1$ 有

$$
 \frac{\prod_{j=0}^2 L^{(23)}(s,\pi_j\times\pi_j)}
      {\prod_{0\le i<j\le2} L^{(23)}(s,\pi_i\times\pi_j)}
   =\left(\frac{\zeta_k^{(23)}(s)}{L^{(23)}(s)}\right)^3 . \tag{6}
$$

对指标按模 3 运算。诱导表示的张量积分解为 $\operatorname{Ind}(\nu^{2m}\psi^{i+j})\oplus\operatorname{Ind}(\psi^{i-j})$，因为 $\xi_j^c=\xi_j^{-1}$。第一项在三条对角和三个异指标对中各取遍一次，Euler 因子逐素数约去。第二项在对角给 $\zeta_k$，异指标给 $L$，因为 $\psi$ 与 $\psi^{-1}$ 的 Hecke $L$ 函数相同。这个消去在每个权都成立。因此只分别控制三个对角范数会遗下恰好含有 $L$ 的混合项。

后篇从 (2) 作主类投影，在 $p^2$ 系数上构造两个权的整数，求出其全部公因子，并检验相应 Stepanov 重数能否给出所需的高度界。[两个 theta 权的整除次数与高度](cubic23-stepanov-content.md)
