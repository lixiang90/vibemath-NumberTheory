---
title: "全体二面体四次域的零点中心带与素数分裂区间"
date: "2026-11-24"
abstract: "假设有理数上的普通 RH 和二次 Dirichlet-RH。在按判别式排序的八阶二面体四次域族中，证明密度一的零点中心带和整段素数分裂估计，并给出固定二次基域时的零点密度界。"
---

# 全体二面体四次域的零点中心带与素数分裂区间

设 $\mathcal F(X)$ 为所有满足 $X<|D_K|\le 2X$ 的四次数域 $K$ 的 $\mathbb Q$-同构类，其中 $D_K$ 为域判别式，伽罗瓦闭包 $M$ 的群为八阶二面体群 $D_4$。计数时允许二次子域、分歧素数和实嵌入数全部变动。本文对普通 ζ 函数及有理数上的本原二次 Dirichlet 特征假设 RH。对新增的二维因子，零点密度估计与迹的素数和均不需要此假设。

**中心带结论。** 给定 $0<\eta<1/2$ 和 $0<h<2\eta/3$，当 $X\to\infty$ 时，$\mathcal F(X)$ 中比例趋于 1 的 $K$ 满足

$$
\left.
\begin{gathered}
  \rho\text{ 是 }\zeta_K\text{ 或 }\zeta_M\text{ 的非平凡零点},\\
  |\Im\rho|\le X^h
\end{gathered}
\right\}
\quad\Longrightarrow\quad
\left|\Re\rho-\tfrac12\right|<\eta . \tag{1}
$$

例如，可同时取中心带 $2/5<\Re s<3/5$ 和高度 $X^{1/20}$。例外域数为 $o_{\eta,h}(X)$。二面体四次域总数在这个判别式区间为 $cX+o(X)$，且 $c>0$。这是按四次域计数的结论，同一闭包可以来自不同的四次域。[Cohen、Diaz y Diaz、Olivier，推论 1.4 和 2.2](https://doi.org/10.1023/A:1016310902973)

**素数分裂结论。** 固定 $0<\alpha<1$ 与 $\delta>0$。对比例趋于 1 的 $K$，同时对每个 $X^\alpha\le x\le X$ 以及 $G=\mathrm{Gal}(M/\mathbb Q)$ 的每个共轭类 $C$，有

$$
\sum_{\substack{p\le x,\ p\nmid D_K\\ \mathrm{Frob}_p\in C}}\log p
=\frac{|C|}{8}x+O_{\alpha,\delta}\big(x^{1/2+\delta}\big). \tag{2}
$$

取 $C=\{1\}$ 就得到完全分裂素数的主项 $x/8$。可让同一个域同时满足 (1) 和 (2)。此处先按判别式取域族，再要求同一成员同时满足整个尺度区间上的估计。

## 一　把全体四次域的新增因子写成二次特征

每个 $K$ 有唯一二次子域 $k$。写 $d=|D_k|$，以 $\psi$ 表示 $K/k$ 的本原二次 Hecke 特征，$\mathfrak f_\psi$ 为其有限导子。导子判别式公式给出 $|D_K|=d^2N\mathfrak f_\psi$，其中 $N$ 是理想的绝对范数。令 $F=M^{[G,G]}$，则 $F/\mathbb Q$ 是双二次域。四根的置换表示分解为 $k$ 的两维置换表示与 $G$ 唯一的二维不可约表示 $\tau$。该不可约表示也是从 $k$ 诱导的 $\psi$。四个一维表示在 $F$ 上因子化，故

$$
L(s,\tau)=L_k(s,\psi),\qquad
\zeta_K(s)=\zeta_k(s)L_k(s,\psi),\qquad
\zeta_M(s)=\zeta_F(s)L_k(s,\psi)^2 . \tag{3}
$$

其中 $\zeta_k,\zeta_F$ 的全部非平凡零点由二次 Dirichlet-RH 控制，普通 $\zeta$ 因子也计入。$L_k(s,\psi)$ 整且对应一个二面体 $\mathrm{GL}_2/\mathbb Q$ 尖点表示。下文对这个新增因子的估计不依赖任何 RH。只有在回到完整的 $\zeta_M,\zeta_K$ 以及分裂素数主项时才使用有理数上的输入。[Chen An，§2 与定理 5.1 的证明](https://www.numdam.org/item/10.5802/jtnb.1109.pdf)

## 二　固定二次基域时的大筛

固定 $k$，置 $R=X/d^2$。取 $R$ 充分大，令 $\mathbb E$ 表示对 $R<N\mathfrak f_\psi\le2R$ 的全部本原二次 Hecke 特征求和，再除以 $R$。此测度的总质量为 $O_k(1)$。若 $c(\mathfrak n)$ 与 $\psi$ 无关且对每个 $\epsilon>0$ 有 $|c(\mathfrak n)|\ll_{k,\epsilon}(N\mathfrak n)^\epsilon$，则

$$
\mathbb E\max_{D\le b\le2D}
\left|\sum_{D<N\mathfrak n\le b}
c(\mathfrak n)\psi(\mathfrak n)(N\mathfrak n)^{-u-iv}\right|^2
\ll_{k,\epsilon}(2RD)^\epsilon
\left(D^{1-2u}+R^{-1}D^{2-2u}\right). \tag{4}
$$

这里 $u$ 在固定紧区间，$v$ 为任意实数，$D\ge1$。范数 1 的单位理想可单列。

固定有限素理想集 $S$，包括 2 上方的素理想，并使 $S$-整数环的理想类群平凡。Goldmakher–Louvel 第 2 节的 Fisher–Friedberg 构造为每个与 $S$ 互素的平方自由理想 $\mathfrak a$ 赋予二次特征 $\chi_{\mathfrak a}$，其导子在 $S$ 外恰为 $\mathfrak a$。给定任意实际的 $\psi$，取其导子在 $S$ 外的部分为 $\mathfrak a$。奇剩余特征处的分歧二次特征在惯性群上有唯一的非平凡限制，所以 $\psi/\chi_{\mathfrak a}$ 在 $S$ 外处处不分歧。这种商只有 $O_k(1)$ 种。事实上由 Kummer 理论，其平方类在 $S$ 外赋值全偶，只有有限的 $S$-单位模平方与理想类群的 2-挠需要计入。在实位允许所有符号。按商分组以后，对与 $S$ 互素的理想，把这个固定商吸入系数即可。每组 $\mathfrak a$ 的范数至多 $2R$。原文定理 1.1 对平方自由的两个指标给出 $(R+B)(2RB)^\epsilon\sum|c_{\mathfrak m}|^2$ 的未归一化上界。对矩阵取伴随即可对调两个指标。[Goldmakher–Louvel，定理 1.1、§2 的导子公式](https://arxiv.org/pdf/1112.1642)

分解 $\mathfrak n=\mathfrak h\mathfrak b^2\mathfrak m$，其中 $\mathfrak h$ 支持在 $S$，$\mathfrak b,\mathfrak m$ 与 $S$ 互素，$\mathfrak m$ 平方自由。对固定的 $\mathfrak h,\mathfrak b$，特征值 $\psi(\mathfrak h)\psi(\mathfrak b)^2$ 的模至多为 1，可在每个 $\psi$ 的内和外提出。令 $B=D/(N\mathfrak h\,(N\mathfrak b)^2)$。非空区段有 $B\ge1/2$，其中的理想数为 $O_k(B)$。上段大筛给出均方根上界

$$
\ll_{k,\epsilon}(2RD)^\epsilon D^{-u}\big(\sqrt B+B/\sqrt R\big).
$$

对最大前缀，将按整数范数排列的前缀分成二进制区段，同层不交，Cauchy 不等式与层数只引入对数损失。最后在均方根范数下相加。用 $\sum_{N\mathfrak b\le 2D}(N\mathfrak b)^{-1}\ll_k\log(2D)$、$\sum_{\mathfrak b}(N\mathfrak b)^{-2}<\infty$ 及 $\sum_{\mathfrak h}(N\mathfrak h)^{-1/2}<\infty$，即得 (4)。每个 $\mathfrak a$ 只对应有界多个 $\psi$，也给出前述总质量界。

## 三　从特征平均到实际零点

设 $0<\xi<\eta/8$，$\sigma_0=1/2+\eta-\xi$，$r=\min(\xi/2,1/64)$。以下取 $\sigma_0\le\sigma=\Re s\le1+r$、$t=\Im s$、$|t|\le H+1$、$1\le H\le R^{1/2}$。令 $Y=\lfloor\sqrt R\rfloor$ 及 $Q=\sqrt R(2+|t|)$。用理想 Möbius 函数 $\mu_k$ 定义

$$
P_\psi(s)=\sum_{N\mathfrak n\le Y}\psi(\mathfrak n)(N\mathfrak n)^{-s},
\qquad
J_\psi(s)=\sum_{N\mathfrak n\le Y}\mu_k(\mathfrak n)\psi(\mathfrak n)(N\mathfrak n)^{-s}.
$$

完备因子为 $q_\psi^{s/2}\gamma(s)$，其中 $q_\psi=dN\mathfrak f_\psi$，$\gamma(s)=\prod_{j=1}^2\Gamma_{\mathbb R}(s+\kappa_j)$，$\kappa_j\in\{0,1\}$。复位给出 $\{\kappa_1,\kappa_2\}=\{0,1\}$。二次扩张的两个 Dedekind 函数的函数方程给出根数 $+1$。取偶核 $e^{w^2}$，在完备函数的竖线积分中越过 $w=0$ 并对左边使用函数方程，就得到精确的两项近似函数方程。若 $n=N\mathfrak n$，主和的权与对偶和的整个权分别为

$$
\begin{aligned}
V_\psi(n)&=\frac1{2\pi i}\int_{(2)}
 \frac{e^{w^2}}{w}\frac{\gamma(s+w)}{\gamma(s)}
 \left(\frac{n}{\sqrt{q_\psi}}\right)^{-w}dw,\\
U_\psi(n)&=\frac{q_\psi^{(1-2s)/2}}{2\pi i}\int_{(2)}
 \frac{e^{w^2}}{w}\frac{\gamma(1-s+w)}{\gamma(s)}
 \left(\frac{n}{\sqrt{q_\psi}}\right)^{-w}dw .
\end{aligned}
$$

故 $L_k(s,\psi)=\sum_{\mathfrak n}\psi(\mathfrak n)[n^{-s}V_\psi(n)+n^{s-1}U_\psi(n)]$。以下权估计也对 $\log n$ 的一阶导数成立，隐含常数允许依赖 $k,\eta,\xi$。在 $n\le Q$，主和积分线移至 $-1/2$，得 $V_\psi(n)-1=O((n/Q)^{1/2})$。取 $\kappa=\eta/4$，对偶和留在正线 $\Re w=\kappa$，得 $|U_\psi(n)|\ll Q^{1-2\sigma}(Q/n)^\kappa$。这里 $1-\sigma+\kappa>0$，覆盖了实偶 Gamma 因子在低高度的情形。对 $n>Q$，移到任意固定的大正数 $A$，相应上界为 $O_A((Q/n)^A)$ 与 $O_A(Q^{1-2\sigma}(Q/n)^A)$。这些界由 Stirling 公式给出，高斯核控制积分变量的大虚部。

逐个 $\psi$ 对权分部求和，再用 (4) 控制最大部分和。对 $n\le Y$ 的主和换权误差，在 $D\asymp n$ 的段上，其均方根至多为

$$
\ll_{k,\epsilon}(2RQ)^\epsilon (D/Q)^{1/2}
      \big(D^{1/2-\sigma_0}+R^{-1/2}D^{1-\sigma_0}\big).
$$

由于 $\sigma_0<1$，两项都从 $D\asymp Y$ 给出最大贡献。单位理想的误差为 $O_k(Q^{-1/2})$。对 $n>Y$ 的主和，第一项从最短段给出 $Y^{1/2-\sigma_0}$，第二项从 $D\asymp Q$ 给出 $R^{-1/2}Q^{1-\sigma_0}$。对偶和在小端的两项，与 $D$ 有关的幂是 $\sigma-1/2-\kappa$ 和 $\sigma-\kappa$，均为正，故相加后分别为 $O((2RQ)^\epsilon Q^{1/2-\sigma})$ 和 $O((2RQ)^\epsilon R^{-1/2}Q^{1-\sigma})$。大端用上述快衰减。由 $Y\asymp\sqrt R$、$\sqrt R\le Q\ll R$ 得到

$$
\mathbb E|L_k(s,\psi)-P_\psi(s)|^2
 \ll_{k,\eta,\xi,\epsilon} R^\epsilon
 \left(Y^{1-2\sigma_0}+R^{-1}Q^{2-2\sigma_0}\right). \tag{5}
$$

在 $P_\psi J_\psi-1$ 中，范数至多 $Y$ 的项由理想卷积恒等式完全消去。余下的系数由理想除数函数控制且范数至多 $Y^2$。$J_\psi^2$ 的系数也有这个界。逐段应用 (4)，得到

$$
\begin{aligned}
\mathbb E|P_\psi J_\psi-1|^2
 &\ll R^\epsilon\big(Y^{1-2\sigma_0}+R^{-1}(Y^2)^{2-2\sigma_0}\big),\\
\mathbb E|J_\psi|^4&\ll R^\epsilon .
\end{aligned} \tag{6}
$$

第二式使用 $D\le Y^2\le R$，于是在 (4) 右侧第二项不超过第一项。式 (5)–(6) 的两个平方误差为 $O(R^{-\eta+O(\xi)+\epsilon})$。在 $L_kJ_\psi-1=(L_k-P_\psi)J_\psi+(P_\psi J_\psi-1)$ 中用 Hölder 不等式，第一项的 $4/3$ 次矩可由两个矩的 $2/3$ 与 $1/3$ 次幂之积控制。先指定误差幂再选充分小的 $\xi$ 和大筛损失，可得

$$
\mathbb E|L_k(s,\psi)J_\psi(s)-1|^{4/3}
 \ll_{k,\eta,\epsilon} R^{-2\eta/3+\epsilon}. \tag{7}
$$

记 $N_\psi(\eta,H)$ 为右侧 $\Re\rho\ge1/2+\eta,\ |\Im\rho|\le H$ 内按重数计的非平凡零点数。以每个这样的零点为中心作半径 $r$ 的圆盘，圆盘在上述矩形内。次调和均值不等式给出 $|L_kJ_\psi-1|^{4/3}$ 在每个盘上的面积积分至少为 $\pi r^2$。含重数的盘重叠次数为 $O_k(\log(2R(H+2)))$。事实上在 $2+iv$ 欧拉积有一致正下界，在同心半径 3 的圆上函数方程和凸性界给出 $(2R(2+|v|))^{O_k(1)}$ 的上界，Jensen 公式即可数出半径 2 盘内的零点。对 (7) 作面积积分得到

$$
\sum_{R<N\mathfrak f_\psi\le 2R} N_\psi(\eta,H)
  \ll_{k,\eta,\epsilon}(1+H)R^{1-2\eta/3+\epsilon}.
  \tag{8}
$$

对固定 $k$ 和 $H=X^h=d^{2h}R^h$，取 $\epsilon<2\eta/3-h$，有右侧偏线零点的域数为 $o_k(X)$。二次特征自对偶，函数方程给出左侧同一结论。由 (3) 得到两种完整的 Dedekind 函数的断言。

## 四　为何可以放开二次子域

式 (8) 的常数依赖 $k$。下面直接控制大二次子域的域数，使这个依赖不影响密度结论。若 $A_k(V)$ 为范数不超过 $V$ 的非零整理想数，对 $V\ge1$ 有一致估计

$$
A_k(V)=\sum_{m\le V}\chi_{D_k}(m)\lfloor V/m\rfloor
 \ll V\log(2d).
$$

将取整替换为 $V/m$ 的总误差至多 $V$。特征 $\chi_{D_k}$ 在整个周期的和为零，对 $m>d$ 分部求和、对 $m\le d$ 取绝对值，就得所述上界。Minkowski 界使每个理想类有范数 $O(\sqrt d)$ 的代表，所以类数 $h_k\ll\sqrt d\log(2d)$。[Milne，定理 4.3–4.4](https://www.jmilne.org/math/CourseNotes/ANT.pdf)

给定 $K=k(\sqrt a)$，把 $a$ 的有限素理想赋值按奇偶性记为一个平方自由理想。每个奇赋值处必分歧，故此理想的范数不超过相对判别式的范数。每一个奇偶赋值向量至多有 $|\mathcal O_k^\times/\mathcal O_k^{\times 2}|\,|\mathrm{Cl}(k)[2]|$ 个平方类原像。这由 $(a)=\mathfrak b^2$ 时的映射 $a\mapsto[\mathfrak b]$ 得到，其核是单位模平方。二次数域的单位模平方数至多 4。于是相对判别式范数至多 $V$ 的二次扩张数

$$
\ll h_k A_k(V)\ll V\sqrt d\,\log^2(2d).
$$

以 $V=2X/d^2$ 代入，若 $V<1$ 则无贡献。每个正整数至多是两个二次域判别式的绝对值。因此对任意 $B\ge2$，

$$
\#\{K\in\mathcal F(X):|D_k|>B\}
\ll X\sum_{n>B}n^{-3/2}\log^2(2n)
\ll X B^{-1/2}\log^2(2B). \tag{9}
$$

先固定 $B$，对有限个 $d\le B$ 使用 (8)，再对 $B$ 取极限，得 (1)。这样得到全族的 $o(X)$ 例外界。固定 $k$ 时 (8) 另给出幂次节省。

## 五　用同一大筛直接检验整段素数窗口

对未分歧有理素数 $p$，

$$
a_\psi(p)=\operatorname{tr}\tau(\mathrm{Frob}_p)
=\sum_{\substack{\mathfrak p\mid p\\N\mathfrak p=p}}\psi(\mathfrak p).
$$

在 (4) 取 $u=v=0$ 并只留绝对剩余次数为 1 的素理想，系数取 $\log N\mathfrak p$。固定 $k$ 时，$V\le2X$ 意味着 $V\ll_k R$。对直到 $V$ 的区段在均方根范数中相加，得到前缀最大值平方的 $\mathbb E$ 平均 $\ll_{k,\epsilon}X^\epsilon V$。在 $X^\alpha\ll V\le2X$ 的二进制尺度，把阈值取为 $V^{1/2+\delta}$，Markov 不等式给出 $o_k(X)$ 个例外。排除这些域后，整个窗口内的带迹加权素数和是 $O_\delta(x^{1/2+\delta})$。删去分歧有理素数的代价至多为 $2\log |D_K|$，因为 $K$ 与 $M$ 的分歧有理素数集合相同。

把每个共轭类的指标按 $D_4$ 的不可约特征标展开，除 $a_\psi(p)$ 之外只剩四个一维特征。其中平凡特征给出主项，其余三个是 $F$ 的二次 Dirichlet 特征。它们的导子都是基本判别式的绝对值，奇素因子均分歧于 $K$，故导子至多为 $8|D_K|$。二次 Dirichlet-RH 给这些和以 $O(\sqrt x\log^2(2Xx))$ 的统一误差。在 $x\ge X^\alpha$ 时可吸入 (2) 的误差。对有限个 $k$ 使用以上估计，再用 (9)，即得 (2)。[Dirichlet-GRH 下素数和的显式估计](https://arxiv.org/abs/2003.01925)

式 (2) 也有一个短区间推论。固定 $0<\delta<1/4$，若 $X^\alpha\le x\le X/2$ 且 $x^{1/2+2\delta}\le y\le x$，在 $x+y$ 与 $x$ 两处相减可得

$$
\sum_{\substack{x<p\le x+y,\ p\nmid D_K\\\mathrm{Frob}_p\in C}}\log p
=\frac{|C|}{8}y+O_{\alpha,\delta}(x^{1/2+\delta}).
$$

因此，对上述密度一的域，每一个这样的区间都有在 $M$ 中完全分裂的素数。误差与主项的比值是 $O_{\alpha,\delta}(x^{-\delta})$。勒让德猜想所问的相邻平方区间长恰为平方根量级，此处保留了正幂余量。

## 六　由全族结论恢复二次 Dirichlet-RH

(1) 还有一个逆向蕴含。如果对所有固定的 $\eta$ 都有 (1)，就可恢复普通 RH 和所有二次 Dirichlet 特征的 RH。理由是每个固定二次域 $k$ 都在正比例的 $\mathcal F(X)$ 中出现。Cohen–Diaz y Diaz–Olivier 推论 1.2 给出相对判别式范数在 $(R,2R]$ 内的二次扩张数渐近于 $b_kR$，其中 $b_k=2^{-r_2(k)}\operatorname{Res}_{s=1}\zeta_k(s)/\zeta_k(2)>0$。在这些扩张中，若四次域对 $\mathbb Q$ 亦伽罗瓦，相对导子在 $k$ 的自同构下不变。在 $2d$ 外，奇导子平方自由，位于分裂素数上的两素理想须成对出现，惰性素理想本身的范数为平方。因此这部分导子形如 $n\mathcal O_k$，其中 $n\ll_k\sqrt R$。有限个坏处的导子指数有界，每个这样的导子对应 $O_k(1)$ 个二次特征，因而共只排除 $O_k(\sqrt R)$。其余每个 $D_4$ 四次域恰由两个 $k$-共轭的扩张给出。于是包含 $k$ 的四次域数渐近于 $b_k X/(2d^2)$。$\zeta_k$ 的任一固定偏线零点都由 (3) 出现在这些 $\zeta_K$ 中，取比其离线距离更小的 $\eta$ 即与 (1) 矛盾。[相对计数与同构重数，推论 1.2、2.2](https://doi.org/10.1023/A:1016310902973)

固定 $\eta$ 与 $h$ 时，(1) 允许某个指定的域属于例外集。固定非阿贝尔因子的全部高度仍须单独控制。Chen An 研究固定双二次子域的族，在靠近 $1$ 的区域得到带例外估计的无零矩形。本文采用按四次域判别式排序的全族，证明使用固定二次基域的大筛与 (9) 的尾界。[Chen An 的无零区域及其域族](https://www.numdam.org/item/10.5802/jtnb.1109.pdf) · [一个固定三次域因子的有限高度证书](artin-108-zero-certificate.md)
