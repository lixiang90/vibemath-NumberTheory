---
title: "CM 型轮换的 Frobenius 次数与判别式 229 的约化斜率"
date: "2026-11-26"
abstract: "计算全实三次域与虚二次域复合后八个 CM 型在六种未分歧 Frobenius 类型下的 Newton 斜率。即使斜率全部为二分之一，固定整环作用切换 r 个切空间符号仍须付出至少 p 的 r 次方的同源次数，特殊纤维上的部分 Frobenius 达到此界。另给出多二次塔平均 CM 高度的精确加法及变权轮换误差。"
---

# CM 型轮换的 Frobenius 次数与判别式 229 的约化斜率

设 $F=\mathbb Q(\alpha)$，$\alpha^3-4\alpha+1=0$，$N$ 为其全实 $S_3$ 正规闭包。取虚二次域 $K$ 及 $E=FK$，则 $NK$ 的群为 $S_3\times C_2$。八个 CM 型由三个符号 $(s_1,s_2,s_3)$ 表示，两个同号型是诱导型，六个混合型是一个 Galois 轨道。原周期公式固定几何对象。要在提升链中轮换这些类型，需要同时了解稳定高度、各素数处的 $p$-可除群和施加于其上的整环作用。[前篇的八个 CM 型和高度差](cm229-heights-and-zeros.md) · [原周期等式的量词](https://github.com/leanprover-community/riemann-hypothesis/blob/9c9f665ed79ce0fa8da61bd803b7a4bf1ad62ed0/src/Development/AbelianVarietyTheory/BlockTheory/GlobalSignature/MeanPeriodLogarithms.lean#L49)

## 一张完整的 Newton 斜率表

令 $p$ 在 $NK$ 不分歧。在固定一个上方素处之后，记 Frobenius 为 $(g,t)$，其中 $g\in S_3$，$t=+1$ 表示 $p$ 在 $K$ 分裂，$t=-1$ 表示惰性。在拥有完整 $\mathcal O_E$ 作用的 CM 阿贝尔簇的良约化上工作，必要时先作使复乘有定义且获得良约化的有限扩张。对六个嵌入 $(j,s)$，Frobenius 的作用是 $(j,s)\mapsto(g(j),ts)$。对每条轨道 $C$，Newton 斜率是 $|\Phi\cap C|/|C|$，重数为 $|C|$。采用对偶约定会将斜率取补，以下多重集保持相同。[Blake 的 CM 型与分解群公式，定理 1.2](https://arxiv.org/abs/1311.5454)

下表每一行将八个 CM 型按完整斜率多重集分类。指数表示重数，$O=(0^3,1^3)$、$S=((1/2)^6)$，另记 $P=(0^2,(1/2)^2,1^2)$、$T=(0,(1/2)^4,1)$、$C=((1/3)^3,(2/3)^3)$。

| $g$ 的循环型 | $t$ | $O$ | $S$ | $P$ | $T$ | $C$ |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 单位元 | $+1$ | 8 | 0 | 0 | 0 | 0 |
| 单位元 | $-1$ | 0 | 8 | 0 | 0 | 0 |
| 换位 | $+1$ | 4 | 0 | 0 | 4 | 0 |
| 换位 | $-1$ | 0 | 4 | 4 | 0 | 0 |
| 三循环 | $+1$ | 2 | 0 | 0 | 0 | 6 |
| 三循环 | $-1$ | 0 | 8 | 0 | 0 | 0 |

例如 $g=(12),t=-1$ 时，固定指标 3 上方的两个符号构成一条二元轨道，各型均在其中选一次。其余两条轨道是 $\{(1,+),(2,-)\}$ 和 $\{(1,-),(2,+)\}$。若 $s_1=s_2$，它们各选一次。若 $s_1\ne s_2$，一条选两次，另一条选零次。由此得到换位与惰性相遇时的两种斜率。$g$ 为三循环且 $t=-1$ 时只有一条六元轨道，必选三次。其余行由相同轨道计算得到。配套 Lean 文件枚举六个置换、两种惰性或分裂取值与八个 CM 型，以普通内核归约验证了全部 96 种组合。[Lean 有限轨道证书](CM229Frobenius.lean)

取 $K=\mathbb Q(i)$、$p=7$ 即有一个具体反例。模 7 的分解为 $(X-4)(X^2+4X+5)$，二次因子的判别式 $3\bmod7$ 非平方。$7$ 在 $\mathbb Q(i)$ 中惰性，且在 $F$ 有一个一次素理想。该一次局部分量确实总有斜率 $1/2$，但八个 CM 型中有四个在整株阿贝尔簇上还出现 0 与 1。就所有未分歧有理素数而言，使八种类型同时超奇异的密度为 $3/12=1/4$，对应 $t=-1$ 且 $g$ 为单位元或三循环。后者在 $F$ 没有一次素理想。若还要求存在一次素理想，则只剩密度 $1/12$ 的单位元一行。[三次域数据](https://www.lmfdb.org/NumberField/3.3.229.1) · [Chebotarev 密度定理](https://www.jmilne.org/math/CourseNotes/CFT.pdf)

## 相同斜率下切换 CM 型的最小局部次数

**命题。** 设一般全实基域 $F$ 的次数为 $d$，$E=FK$，$p$ 在 $F$ 完全分裂且在虚二次域 $K$ 惰性，在两域均不分歧。固定 $\overline{\mathbb F}_p$ 上各局部整环因子的标号。若两个带 $\mathcal O_E$ 作用的良约化具有相差 $r$ 个符号的 Lie 型，则对每个保持这一作用的同源 $\phi$，都有

$$
 v_p(\deg\phi)\ \ge\ r. \tag{1}
$$

在特殊纤维上，给定起点和要翻转的 $r$ 个位置，存在达到等号的局部 $p$-可除群同源，其靶有指定的 Lie 型。起点若来自阿贝尔簇，也可取其有限子群概形的商，从而在特殊纤维上实现这一操作。

**证明。** 此时 $\mathcal O_E\otimes\mathbb Z_p\simeq\prod_{j=1}^d\mathbb Z_{p^2}$，其中 $\mathbb Z_{p^2}=W(\mathbb F_{p^2})$。各幂等元切出高度二、维数一且斜率均为 $1/2$ 的连通因子 $G_j$。在非分歧处将底环扩至严格 Hensel 环，嵌入对应的幂等元在好模型上分裂，Lie 模各特征空间的秩从特征零保持到特殊纤维。固定标号后，$\mathbb F_{p^2}$ 在 $\operatorname{Lie}G_j$ 上按两个嵌入 $\lambda$、$\lambda^p$ 之一作用。若源与靶选取不同嵌入，一个交换此作用的切映射 $D_j$ 必须为零，因为可取 $a\in\mathbb F_{p^2}$ 使 $\lambda(a)\ne\lambda(a)^p$。若该因子上的同源高度为零，它就是同构，切映射必可逆。故每个翻转位置的同源高度至少为一，总次数的 $p$ 赋值是各因子高度之和，得到 (1)。

在需要翻转的分量用相对 Frobenius $G_j\to G_j^{(p)}$，其余分量用恒等。相对 Frobenius 的自然性保证在靶上取扭转后的作用时，映射交换 $\mathbb Z_{p^2}$ 的作用。切空间嵌入变为 $\lambda^p$。一维形式参数上的相对 Frobenius 将参数送到 $p$ 次方，核的长度为 $p$，且包含在 $G_j[p]$ 中。所选各核的乘积于是为阶 $p^r$ 的 $\mathcal O_E$ 稳定子群，取商得到所称的靶。这证明了局部最小次数的可达性。[Moonen 关于 $p$-可除群、Frobenius 与 Verschiebung 的讲义](https://www.math.ru.nl/~bmoonen/BookAV/TateBT.pdf)

在三次例子里，$p=347$ 满足这些条件，$f\bmod347$ 的三个根是 $9,92,246$，且 $347\equiv3\bmod4$。若逐次走遍模整体共轭的四个 CM 型再回到初始符号，每条边至少翻转一个符号，四次局部同源的次数乘积必被 $p^4$ 整除。取

$$
 (+++)\longrightarrow(++-)\longrightarrow(+--)
 \longrightarrow(+-+)\longrightarrow(+++)
$$

可使每次只翻一处，特殊纤维上的部分 Frobenius 正好实现次数乘积 $p^4$。末端恢复的是标号后的 Lie 型。若只要求稳定高度按 $1:3$ 配平，可以重复同一混合型，四个不同类型的遍历条件随之改变。要把特殊纤维的同源用于特征零提升链，仍需处理与指定极化、周期格及 Taylor 系数相容的整模型。特别是相对 Frobenius 在改动位置的切映射为零，不能当作保持一阶系数的可逆换坐标。[原 HypShape 使用的整系数、范数关系及形式群条件](https://github.com/leanprover-community/riemann-hypothesis/blob/9c9f665ed79ce0fa8da61bd803b7a4bf1ad62ed0/src/Formalism/LocalToGlobal/TorsionIntegrality/UniformGaussModels/UniformGaussModelCondition.lean#L85)

伴篇把固定 $E=F(i)$ 的 ERH 归结到同时在两个移动虚二次域惰性、在 $F$ 完全分裂的有理素数。对这些移动宿主，所有 CM 型已经位于表中第二行。因而在这个判据里，可集中处理 (1) 所展示的整结构，以及相应的几何主项。[双惰性素数的固定域转移](winpair-inert-transfer.md)

## 多二次塔中的平均 Colmez 高度

这个恒等式适用于任意全实基域 $F$。取全实的多二次伽罗瓦扩张 $R/F$，次数 $m=2^k$，以及固定 CM 二次扩张 $E_0/F$，特征为 $\eta_0$。令 $M=R E_0$，$X=\widehat{\operatorname{Gal}(R/F)}$，$E_\chi/F$ 为奇特征 $\eta_0\chi$ 对应的 CM 二次扩张。以 $\bar h$ 表示前篇归一化下的平均稳定 Faltings 高度。则有精确加法

$$
 \boxed{\displaystyle\bar h_M=\sum_{\chi\in X}\bar h_{E_\chi}} . \tag{2}
$$

事实上，$M$ 的最大实子域为 $R$，绝对次数为 $2m[F:\mathbb Q]$。由 Artin 诱导及阿贝尔正则表示分解，

$$
 L_R(s,\eta_{M/R})=\prod_{\chi\in X}L_F(s,\eta_0\chi),\qquad
 \frac{|D_M|}{D_R}=\prod_{\chi\in X}\frac{|D_{E_\chi}|}{D_F}.
$$

第二式来自导子判别式公式，取本原局部因子，容许各扩张在相同有限素数分歧。应用平均 Colmez 定理，$s=1$ 处的对数导数相加，判别式比值的对数相加，由次数决定的常数项也因 $[R:\mathbb Q]=m[F:\mathbb Q]$ 而相加。因此 (2) 包括全部常数。特别地，当 $R/F$ 为二次时，复合域在 $F$ 上的次数为 4，且有 $\bar h_{RE_0}=\bar h_{E_0}+\bar h_{E_1}$，$E_1$ 为另一个 CM 子域。所有右端的相对特征在无穷位仍为奇。固定偶部由伴篇的素数计数恢复。[平均 Colmez 定理 A](https://www.math.mcgill.ca/goren/PAPERSpublic/Colmez_v3.5.pdf) · [平均高度在 $s=1$ 的完整常数](real-cm-prime-windows.md) · [Milne 的导子判别式公式](https://www.jmilne.org/math/CourseNotes/CFT.pdf)

## 不等权的四步与高秩 Gram 矩阵

对 $F(i)$，设 $\Delta=h_M-h_I=1.1323555935686836618893\ldots$，其中本节的 $h_M$ 表示混合型高度。若第 $j$ 步的权是 $a_j$，取 $t_j=1$ 表示选诱导型，其余取 0。逐型高度主项与全型平均主项的差为 $-\Delta\sum_j a_j(t_j-1/4)$。将诱导步安排在每四步中的第二步，前缀差 $D_j=\sum_{k\le j}(t_k-1/4)$ 的四个循环值为 $-1/4,1/2,1/4,0$。分部求和给出

$$
 \left|\sum_{j=1}^{m}a_j(t_j-1/4)\right|
 \le |a_m|\,|D_m|+\frac12\sum_{j<m}|a_{j+1}-a_j|. \tag{3}
$$

若 $4\mid m$，端点项消失。任何从第一步起计的安排，在第二步的诱导步数均与 $2/4$ 相差至少 $1/2$，故前缀差常数已最优。等权完整块确实精确消去高度差，对变动的微分阶数则由 (3) 支付总变差。

行列式单位平衡控制平均斜率。列之间的夹角另有自由度，甚至在整系数和行列式一的情形仍然如此。取 $B_t=\left(\begin{smallmatrix}1&t\\0&1\end{smallmatrix}\right)$、$t\in\mathbb N$，在任意秩 $k\ge2$ 时直加单位矩阵。在每个无穷位 $\det(B_t^*B_t)=1$，而两列各自归一化之后的相关系数为 $t/\sqrt{1+t^2}$，趋于 1。给整格乘一标量单位保持该系数。这个例子说明还必须利用实际选取的列映射及其轨道和。若确实另证各局部算子满足 $\|K_v-I\|\le\epsilon_v$ 且全局算子分解为其张量积，则逐项展开张量积可得 $\|\bigotimes_v K_v-I\|\le\prod_v(1+\epsilon_v)-1$。

在原 Lean 中，高秩的数值估计压低了一个已经识别出的退化轨道误差。把它送入算子结论的抽象定理仍使用绝对收敛、行正交、主轨道及退化轨道的估计。这里已经给出了特殊纤维的局部同源及其精确次数，将这些同源的 Taylor 整除收益和固定列 Gram 界接入数域上的几何不等式仍是独立的数学问题。[原 PairBoundAbstract 的输入](https://github.com/leanprover-community/riemann-hypothesis/blob/9c9f665ed79ce0fa8da61bd803b7a4bf1ad62ed0/src/Formalism/HermitianCurvature/StatementPB.lean#L192) · [退化项的数值估计 t_NumNG](https://github.com/leanprover-community/riemann-hypothesis/blob/9c9f665ed79ce0fa8da61bd803b7a4bf1ad62ed0/src/Formalism/HermitianCurvature/HermitianCurvatureProofs/GrowthOfBasic.lean#L157)

在含有证书的目录运行下列命令可核验整张有限轨道表以及模 7、37、347 的三个例子。证书只依赖 Lean 标准库。约化、同源、解析转移及高度公式的证明分别在两篇正文中。[Lean 源文件](CM229Frobenius.lean)

~~~sh
lean +leanprover/lean4:v4.35.0 CM229Frobenius.lean
~~~
