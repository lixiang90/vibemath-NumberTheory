---
title: "一个辅助导子的收缩零点窗与固定差二"
date: "2026-11-25"
abstract: "沿用 Alpöge–Claude 的辅助素数尺度和 Dirichlet-GRH 输入。一个辅助素数导子的全部奇特征在同一段收缩高度窗内的零点，即可用二次式恢复整个物理窗口的孪生素数计数。在加性频率一侧，分母达到平方根除以对数幂的主弧承载渐近一半的平方质量，其固定差二系数给出 Hardy–Littlewood 奇异级数。由此将原共同好窗口的严格正性归结为补弧上一个带符号的二阶系数。另给出同一导子族在两个相近高度的协方差极限。"
---

# 一个辅助导子的收缩零点窗与固定差二

比较最近三篇，最值得沿用的是**先保留特征标签和端点相位，再作二点平均**。第 28 篇在四个零点窗中读出素数和互补因子，第 29 篇用辅助素数换号测出陈氏半素数的质量，第 30 篇用多区间平方和标定这种波动。后两种办法能解释正号出现的频率，但扣除半素数后仍要证明一个严格超额。本篇保留第 28 篇的标签，改变取导子的步骤。原来的辅助素数远大于物理端点，可以固定其中一个，对整个物理区间使用同一组零点。[四窗的因子公式](https://lixiang90.github.io/vibemath-NumberTheory/zeros-exact-distribution/articles/gap-two-separated-zero-windows.html) · [换号后的零点矩](https://lixiang90.github.io/vibemath-NumberTheory/zeros-exact-distribution/articles/gap-two-chen-auxiliary-switch.html) · [多区间的波动与门槛](https://lixiang90.github.io/vibemath-NumberTheory/zeros-exact-distribution/articles/gap-two-window-quadratic-variation.html)

这样得到两个不同层次的结论。给定一批带特征标签的零点，式 (5) 能检出物理区间中单独一对孪生素数。为了从分布定理推出孪生素数的存在，还需要估计这个二次式。式 (10) 算出了其中承载渐近一半平方质量的主弧，式 (16) 将原窗口所缺的严格下界定位到补弧。下文采用归档的 Dirichlet-GRH 作为输入，其余推导是纸面证明。[归档定理的陈述与验证方式](https://github.com/leanprover-community/riemann-hypothesis#theorem-statements)

## 1. 固定导子，在同一高度窗同时读出所有小整数

在原 ζ 构造中固定 $0<\epsilon<1/52$ 和幂指数 $J$。令 $P\to\infty$，$H=P^J$，$z=P^{\epsilon/5}$，$z'=z\log^2P$。原构造先取 $e^{z/2}<Q\le e^z$、模 $8$ 余 $3$ 的素数，再取 $e^{z'/2}<D\le e^{z'}$、模 $16$ 余 $7$ 的素数。后一素数在 $Q$ 和每个不超过 $z$ 的奇素数模下余 $-1$。原几何条件可在选出这一对以后确定一个完整的 $[Y,Y^2]\subset[P,H]$。以下估计对每一对和所有这些 $Y$ 同时成立。原来的好对条件控制指定的二次扭转。这里扩展到模 $Q$ 或 $D$ 的全部奇特征时，用的是归档最终的 Dirichlet-GRH 定理。[辅助素数的定义和尺度](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/MeanSquareBounds/PrimeSetTwoBounds.lean#L37) · [选取辅助对和完整窗口的次序](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/CharacterBoxes/WinPairs.lean#L32)

取 $v=Q$ 或 $D$，$L=\log v$，任取实窗心 $h$，其中 $|h|\le\exp(P^C)$、$C>0$ 固定。给定误差指数 $a>0$，设

$$
D_h=2+\log(v(3+|h|)),\quad
k=\lceil(a+3)\log(2P)+3\log D_h\rceil,\quad
B=L/4,\quad b=L/(16k),\quad R=32k/L .
\tag{1}
$$

取函数 $g=\mathbf1_{[-B,B]}*((2b)^{-1}\mathbf1_{[-b,b]})^{*2k}$。它在 $|u|\le L/8$ 上为 $1$，在 $|u|>3L/8$ 上为零。用 Fourier 约定 $\widehat g(t)=\int g(u)e^{itu}\,du$，得到

$$
\widehat g(t)=\frac{2\sin(Bt)}t\left(\frac{\sin(bt)}{bt}\right)^{2k},\qquad
S_\chi(h)=\sum_{|\gamma_\chi-h|\le R}\widehat g(\gamma_\chi-h),\qquad
Z_{v,h}(c)=-\frac{2}{v-1}\sum_{\chi(-1)=-1}\overline{\chi(c)}S_\chi(h).
\tag{2}
$$

在可去奇点取连续值。对模素数 $v$ 的每个奇特征，$\tfrac12+i\gamma_\chi$ 遍历其 $L$ 函数在临界带内的零点，按重数计。复特征自身的正负高度都须保留。半宽 $R=O_{a,J,C,\epsilon}(\log P/\log v)=o(1)$。当窗心固定时，这约为 $O(\log P)$ 个平均零点间距。全部 $(v-1)/2$ 个特征都参与投影，这个族随 $P$ 增长很快。

**命题 1。** 对所有整数 $2\le c\le H+2$，一致有

$$
Z_{v,h}(c)=\frac{\Lambda(c)}{\sqrt c}\,c^{-ih}+O(P^{-a}).
\tag{3}
$$

证明先将 $S_\chi$ 扩展到所有高度，记为 $S_\chi^\infty$。平移后的显式公式给出

$$
S_\chi^\infty(h)=A_{v,h}
-\sum_{m\ge2}\frac{\Lambda(m)}{\sqrt m}g(\log m)
\big(\chi(m)m^{-ih}+\overline{\chi(m)}m^{ih}\big).
\tag{4}
$$

同导子、同奇偶性的导子项和 Gamma 项合为共同实数 $A_{v,h}$。平移的符号可直接从围道积分看出。令 $H_h(s)=\int g(u)e^{(s-1/2-ih)u}\,du$，在 $\Re s>1$ 展开对数导数，Fourier 反演给出 $g(\log m)m^{-ih}$。左侧直线经函数方程反射后得到共轭项。$H_h$ 在每条固定竖直带内按 $O((1+|\Im s|)^{-2k-1})$ 衰减，换线可用避开零点的一列高度进行。未平移的归一化可对照 Hughes–Rudnick 式 (2.1)。[Hughes–Rudnick，第 2 节](https://arxiv.org/pdf/math/0208230#page=4)

奇特征正交性将 $\chi(m)$ 投影成 $\mathbf1_{m\equiv c}-\mathbf1_{m\equiv-c}$，将 $\bar\chi(m)$ 投影成 $\mathbf1_{mc\equiv1}-\mathbf1_{mc\equiv-1}$，同余均模 $v$。对所取 $c$ 共同项消失。支撑给出 $m\le v^{3/8}$，而 $H+2<v^{1/8}$ 对一切充分大的 $P$ 成立。于是第一项只剩实际整数 $m=c$，同时 $1<mc<v-1$ 排除全部对偶项。此时 $g(\log c)=1$。

最后估计删掉的高度。原始 Dirichlet $L$ 函数在任意实高度 $t$ 的单位区间中有 $O(\log(v(3+|t|)))$ 个零点。此界可由对称零点计数公式作差得到，在原点附近直接用固定高度计数。在 $2^jR<|\gamma-h|\le2^{j+1}R$ 上，$|\widehat g|\ll L\,2^{-2k(j+1)}$，零点数为 $O((1+2^jR)(D_h+\log(2+2^jR)))$。因为 $R\le1$ 最终成立，总尾和为 $O(LD_h4^{-k})=O(P^{-a})$。这对每个特征都成立，投影的权重绝对值之和为 $1$。因此窗心和 $c$ 都能在事后选定。[Bennett、Martin、O’Bryant、Rechnitzer，定理 1.1](https://arxiv.org/pdf/2005.02989#page=2)

## 2. 一个零点二次式恢复固定差二

记 $\mu$ 为 Möbius 函数。取 $v,v'\in\{Q,D\}$，窗心可相同或不同。对任意 $\mathcal I\subset[2,H]\cap\mathbb Z$ 和实权 $|w_n|\le1$，定义

$$
\begin{aligned}
A_{v,h}(n)&=\frac1{\log n}\sum_{\substack{r^j=n\\r\ge2,\ j\ge1}}
       \mu(j)\sqrt r\,r^{ih} Z_{v,h}(r),\\
\mathcal T(w)&=\Re\sum_{n\in\mathcal I} w_n A_{v,h}(n)\overline{A_{v',h'}(n+2)} .
\end{aligned}
\tag{5}
$$

给定 $A>0$，在 (1) 取 $a>A+4J+4$。则对于全部上述选择，包括依赖所取好对的区间和权重，

$$
\mathcal T(w)=\sum_{\substack{n\in\mathcal I\\n,n+2\ {\rm prime}}}w_n+O(P^{-A}).
\tag{6}
$$

事实上 $\sum_{r^j=n}\mu(j)\Lambda(r)=\log n\,\mathbf1_{\mathbb P}(n)$。当 $n=p^m$ 时左边是 $\log p\sum_{j\mid m}\mu(j)$，含两种素因子时每项均为零。可选根的数目不超过 $\log n/\log2$，因此 $A_{v,h}(n)=\mathbf1_{\mathbb P}(n)+O(\sqrt H P^{-a})$，对 (5) 至多 $H$ 行求和即得 (6)。构造 (5) 只需列举整数和它们的整数次根。对于 $n+2=rs$、$r\ne s$ 都是素数的陈氏行，(4) 在 $n+2$ 标签上的算术项精确为零，高次素数幂则由 Möbius 和逐项消去。

展开 (5)，一对奇特征 $\chi\bmod v,\psi\bmod v'$ 的系数除去共同归一化 $4/((v-1)(v'-1))$ 后为

$$
\sum_{n\in\mathcal I}\frac{w_n}{\log n\log(n+2)}
 \sum_{r^j=n,\ s^\ell=n+2}\mu(j)\mu(\ell)\sqrt{rs}\,r^{ih}s^{-ih'}
 \overline{\chi(r)}\psi(s).
\tag{7}
$$

它乘上 $S_\chi(h)S_\psi(h')$，再对两特征求和取实部。取 $v'=v,h'=h$ 即为同一批零点的二次式。两个谱标签通过整数差 $2$ 相连，所有区间共用同一向量 $(S_\chi(h))_\chi$。

- 取 $w_n=1$，充分大时四舍五入直接得到计数。取 $P=2^j$，$\mathcal I=(P,2P]\cap\mathbb Z$，并让第二端 $n+2$ 可以越过 $2P$。孪生素数无限当且仅当无穷多个这样的区间满足 $\mathcal T(1)>1/2$。这里的 $1/2$ 是整数计数的舍入界。
- 在原窗 $I_Y=(2Y,Y^{3/2}/2]$，取 $w_n=\log n\log(n+2)I_{Q,D}(n)/n$。其中 $I_{Q,D}(n)$ 要求 $n,n+2$ 在 $Q$ 处都是二次剩余、在 $D$ 处都是二次非剩余。除以 $L_Y=\sum_{Y\le m\le Y^{3/2}/2}1/m$（$m$ 取整数）后，得到前文的 $J_{Q,D}(Y)+O(P^{-A}/L_Y)$。原好窗可以在取定 $Q,D$ 以后选。该序列上的 $\Gamma_1$ 等于 $J_{Q,D}$ 上极限的三分之一。[共同好窗口的定义和归一化](https://github.com/lixiang90/vibemath-NumberTheory/blob/main/Goldbach_TwinPrimes/twin-prime-auxiliary-signs.md)

## 3. 同导子双高度的协方差能看到什么

在模 $v$ 的奇特征上取均匀平均，协方差按总体平均定义。两窗取共同的 $k$，用两高度中较大的 $D_h$ 定它。对两个容许的窗心 $h,h'$，

$$
\operatorname{Cov}_\chi(S_\chi(h),S_\chi(h'))=
2\sum_{m\le v^{3/8}}\frac{\Lambda(m)^2}{m}\,g(\log m)^2
\cos((h'-h)\log m)+O(LP^{-a}+P^{-2a}).
\tag{8}
$$

在 (4) 中减去特征均值后，只有素数幂和剩下。乘开两份和，$mn\le v^{3/4}<v-1$ 排除了 $\chi(mn)$ 的平凡和负号同余。含 $\chi(m)\bar\chi(n)$ 的项只留下整数对角 $m=n$。这就证明了全高度的精确恒等式。对角和为 $O(L^2)$，故逐特征删尾的 $L^2$ 误差经 Cauchy–Schwarz 至多造成 (8) 的误差。

若 $h'=h+s/L$、$s$ 固定，则

$$
\frac{\operatorname{Cov}_\chi(S_\chi(h),S_\chi(h+s/L))}{L^2}
\longrightarrow 2\int_0^{1/4}u\cos(su)\,du .
\tag{9}
$$

由素数定理及分部求和，$\sum_{m\le x}\Lambda(m)^2/m\sim(\log x)^2/2$。核 $g(Lu)$ 是把区间 $[-1/4,1/4]$ 的指示函数与一个方差 $1/(384k)$ 的概率分布卷积所得。它在离开两个端点的闭集上一致趋于该指示函数。先删去端点的固定小邻域，用上述部分和控制其质量，再令邻域收缩，即得 (9)。方差常数是 $1/16$，取 $s=4\pi$ 的相关系数趋于 $-4/\pi^2$。这是限制 Fourier 支撑下的平滑、带标签的同族零点统计，其显式公式背景见 Hughes–Rudnick。式 (8) 的算术项只含 $\Lambda(m)^2$，固定差二的交叉相位位于 (7)。这解释了为什么仅加强这种平均仍难给出孪生下界。[Hughes–Rudnick，线性统计量的矩](https://arxiv.org/abs/math/0208230)

## 4. 固定差二的主弧占多少平方质量

取 $P\le U<2U\le H$，并记 $e(t)=\exp(2\pi it)$。用同一导子、同一零点窗定义

$$
F_v(\alpha)=\sum_{U<n\le2U}\log n\,A_{v,h}(n)e(n\alpha).
$$

令 $\lambda=\log U$，$W=\lambda^2$，$K=\lfloor\sqrt U/\lambda^8\rfloor$。在圆 $\mathbb R/\mathbb Z$ 上，以所有既约 $a/q$、$q\le K$ 为中心，取半宽 $W/U$ 的弧，其并集记为 $\mathfrak M(U)$，补集记为 $\mathfrak m(U)$。$2K^2W<U$ 保证诸弧不交。则

$$
\begin{aligned}
\frac1U\int_{\mathfrak M(U)} |F_v(\alpha)|^2\cos(4\pi\alpha)\,d\alpha
 &=\mathfrak S(2)+O(\lambda^{-2}),\\
\int_{\mathfrak M(U)} |F_v(\alpha)|^2\,d\alpha
 &=U\big(\tfrac12\log U-8\log\log U+O(1)\big),\\
\mathfrak S(2)&=2\prod_{p>2}\left(1-\frac1{(p-1)^2}\right).
\end{aligned}
\tag{10}
$$

全圆的平方质量是 $U\log U+O(U)$。因此渐近一半的平方质量落在 (10) 这些弧上。Hardy–Littlewood 的奇异级数以及平方根分母下的一半现象有经典先例。Goldston 在 Farey 弧上研究了素数指数和的主弧近似。这里需要的是**由同一窄高度窗恢复的二次式在指定弧集上的第二个 Fourier 系数**，以及对原选窗可用的误差。[Hardy–Littlewood，Partitio Numerorum III](https://archive.ymsc.tsinghua.edu.cn/pacm_download/117/5327-11511_2006_Article_BF02403921.pdf) · [Goldston，特别是引理 2](https://bibliotekanauki.pl/articles/1207182.pdf)

**证明。** 由 (3)–(6)，可把 $F_v$ 换为 $F(\alpha)=\sum_{U<p\le2U}\log p\,e(p\alpha)$。逐点误差为 $O(H^{3/2}\log H\,P^{-a})$，在任意可测子集上对平方积分的影响都为 $O(P^{-A})$。在 $a/q+\beta$ 附近，近似函数为 $M_q(\beta)=\mu(q)V(\beta)/\varphi(q)$，其中 $V(\beta)=\int_U^{2U} e(t\beta)\,dt$。对每个 $q\le K$，

$$
\sum_{a\bmod q}^{*}|F(a/q+\beta)-M_q(\beta)|^2
\ll qU\lambda^4(1+U|\beta|)^2.
\tag{11}
$$

为说明 $q$ 的幂次，在每个可逆剩余类 $b$ 上，令 $E_b$ 为该类带相位的素数和减 $V(\beta)/\varphi(q)$。Dirichlet-GRH 给出各特征在 $(U,x]$ 上去掉主项后的 $\vartheta$ 和为 $O(\sqrt U\log^2(qU))$，对 $U\le x\le2U$ 一致。分部求和加上相位后多出 $1+U|\beta|$。对乘法特征使用 Parseval，得到 $\sum_b|E_b|^2\ll U\lambda^4(1+U|\beta|)^2$。这些 $p>q$ 都在可逆类上。再用加法 Parseval，全部分子 $a\bmod q$ 的误差平方和为 $q\sum_b|E_b|^2$。因 $\sum_{b\bmod q}^{*}e(ab/q)=\mu(q)$ 对既约 $a$ 成立，限制到这些分子即得 (11)。[Lee，Dirichlet-GRH 下的素数定理](https://arxiv.org/abs/2301.13457)

积分并对 $q$ 求和，得到

$$
\int_{\mathfrak M}|F-M|^2\ll K^2 W^3\lambda^4\ll U\lambda^{-6}.
\tag{12}
$$

这里 $M(a/q+\beta)=M_q(\beta)$。另一方面 $\int_{\mathfrak M}|M|^2\ll U\log U$，故在有界相位下将 $|F|^2$ 换为 $|M|^2$ 的积分误差是 $O(U\lambda^{-5/2})$。

记 $c_q(d)=\sum_{a\bmod q}^{*}e(ad/q)$。带相位的模型积分分解为 $\sum_{q\le K}\mu(q)^2c_q(2)/\varphi(q)^2$ 乘上 $\Re\int_{-W/U}^{W/U}|V(\beta)|^2e(2\beta)\,d\beta$。后一因子是 $U-2+O(U/W)$，因为全实线积分给出区间与其平移 $2$ 的重叠长度，且 $|V(\beta)|\ll\min(U,1/|\beta|)$。对于平方自由的 $q$，$|c_q(2)|=1$。由 $\sum \mu(q)^2 q^{1/2}/\varphi(q)^2<\infty$，延长求和的代价是 $O(K^{-1/2})$。其 Euler 乘积恰为 (10) 的 $\mathfrak S(2)$。

去掉相位，则算术因子变为 $\sum_{q\le K}\mu(q)^2/\varphi(q)=\log K+O(1)$，全实线积分变为 $U$。为核对对数系数，该级数的 Dirichlet 生成函数除以 $\zeta(1+s)$，余下局部因子为 $1+(p^{-s}-p^{-2s})/(p(p-1))$。乘积在 $\Re s>-1/2$ 绝对收敛，在 $s=0$ 值为 $1$。与调和和卷积就给出所述估计。结合 (12) 和 Parseval 即证 (10)。

## 5. 原共同窗口还缺哪一个单边估计

若只考虑二倍段，在 $U$ 的无穷子列上证明

$$
\frac1U\int_{\mathfrak m(U)}|F_v(\alpha)|^2\cos(4\pi\alpha)\,d\alpha
\ \ge -\mathfrak S(2)+\delta,\qquad \delta>0\ \text{固定},
\tag{13}
$$

就有 $\gg_\delta U/\log^2 U$ 对孪生素数两端同在该段。这比单对检测要求更高。$\cos(4\pi\alpha)$ 可正可负。补弧上剩余的无符号质量仍有约 $\tfrac12 U\log U$，所以 (10) 中的一半本身不能给 (13) 所需的符号。

为保持原来的归一化，设 $Z_Y=Y^{3/2}/2$。从 $u=2Y$ 起，将 $(2Y,Z_Y]$ 依次切成 $B_j=(u_j,b_j]$，其中 $b_j=\min(2u_j,Z_Y)$。每段在本身的 $\mathfrak m(u_j)$ 上积分，统一使用先前所取的 $v,h$。定义

$$
G_j(\alpha)=\sum_{n\in B_j}\frac{\log n}{\sqrt n}A_{v,h}(n)e(n\alpha),\qquad
\mathcal E_v(Y)=\frac1{L_Y}\sum_j\int_{\mathfrak m(u_j)}
       |G_j(\alpha)|^2\cos(4\pi\alpha)\,d\alpha .
\tag{14}
$$

若 $M(Y)=L_Y^{-1}\sum_{p\in I_Y,\ p+2\ {\rm prime}}\log p\log(p+2)/p$，其中 $p$ 取素数，则一致有

$$
M(Y)=\mathfrak S(2)\frac{\log(Z_Y/(2Y))}{L_Y}
       +\mathcal E_v(Y)+O((\log Y)^{-2}).
\tag{15}
$$

证明时将 (11) 的分部求和权乘以 $t^{-1/2}$，右端除以 $u_j$，对任意上端 $b_j\le2u_j$ 都成立。模型积分使用 $V_j(\beta)=\int_{u_j}^{b_j}t^{-1/2}e(t\beta)\,dt$，有 $|V_j|\ll\min(\sqrt{u_j},1/(\sqrt{u_j}|\beta|))$。其全实线带相位积分是 $\int_{u_j}^{b_j-2}dt/\sqrt{t(t+2)}=\log(b_j/u_j)+O(1/u_j)$，空重叠时积分取零。由 (10) 的证明，每段主弧误差为 $O((\log u_j)^{-2})$。原窗中跨段或跨右端点的行，每个边界至多两行，总权为 $O(\log^2Y/Y)$。把 $\sqrt{p(p+2)}$ 换回原分母 $p$ 的总代价也为 $O(\log^2Y/Y)$。主项对数相加即为 $\log(Z_Y/(2Y))$，共有 $O(\log Y)$ 段，而 $L_Y\asymp\log Y$，证得 (15)。

前文在基准 $P$ 给出有限个候选 $Y_i$，以及原解析好对选择第一个合格完整窗口的条件概率 $\nu_i(P)$。因为 (15) 对每一选择都一致，$\sum_i\nu_i=1$ 且 $\log(Z_Y/(2Y))/L_Y=1+O(1/\log Y)$，原第 18 篇式 (11) 严格等价于

$$
\limsup_{P\to\infty}\sum_i\nu_i(P)\mathcal E_v(Y_i)>-\mathfrak S(2).
\tag{16}
$$

这里每个 $P$ 甚至可以只用预先选定的一个可用导子 $v$，所有 $\mathcal E_v(Y_i)$ 共用它的零点。如果左端比右端大 $\delta>0$，前文的好对平均给出可选见证上的 $\Gamma_1\ge\delta/48$。本文算出了 (15) 的整个主弧项。补弧的严格超额、$\Gamma_1>0$ 及孪生素数无限性仍待证明。[原辅助对平均与共同好窗口](https://github.com/lixiang90/vibemath-NumberTheory/blob/main/Goldbach_TwinPrimes/twin-prime-auxiliary-signs.md)

## 6. 实际零点与主弧的独立复算

在小导子 $43,71$ 分别使用全部 $21,35$ 个奇特征。先由 Hurwitz ζ 值及离散特征变换取得实 Hardy 函数，在临界线上找变号区间，再用 Arb 的逐特征 Hardy 函数求根。另在 $\Im s=\pm 51.37$ 的水平边上以两种网格延拓 $\arg L$，由函数方程、Gamma 因子和辐角原理算得两族在这段高度内各有 $1680,3080$ 个零点，逐特征与找到的根数相符。这里的网格及求根是浮点核对。

用小尺度核 $B=1.55,\ b=0.16,\ k=2$，半宽 $48$，窗心分别为 $1.25,-1.75$。其算术支撑只含 $m\le8$。在这两个小导子上，部分标签位于平台边缘，程序逐标签除以已知的 $g(\log r)>0$。对 $43$ 的标签 $r\le5$ 和对 $71$ 的标签 $r\le7$，$8r<v-1$，故 (3) 的同余排除仍然成立。原好对的全套选窗条件和渐近收缩率由前面的证明处理。下表由两个高度窗中的实际零点直接得到。

| 左端 $n$ | 右端 $n+2$ | 式 (5) 的单行实部 |
| ---: | ---: | ---: |
| 2 | 4 | −0.00000423 |
| 3 | 5 | 0.99999937 |
| 4 | 6 | 0.00000000 |
| 5 | 7 | 1.00037814 |

总和 $2.00037328$ 舍入为 $2$。若两端都只用导子 $71$、窗心 $-1.75$ 的同一批零点，得到 $2.00036973$，舍入同样给出 $2$。两个导子分别使用 (4) 的算术展开复核，中心化零点和的最大绝对误差为 $2.32\cdot10^{-5},2.67\cdot10^{-5}$。在导子 $71$ 的 $-1.75,0.85$ 两个高度，零点给出协方差 $-1.088418148$，(8) 右侧对角和给出 $-1.088414114$。

下表单独检验 (10) 的频率分解。在所列二倍区间上取 $W=4$ 和给定有限分母上限 $K$。一种算法用区间内真实素数的指数和逐弧积分，另一种用全部整数差的自相关、Ramanujan 和以及每段弧的精确 Fourier 积分。两者在各行吻合。表中的有限参数用于复算分解，渐近定理中的 $K,W$ 由第 4 节指定。

| $U$ | $K$ | 主弧占全圆平方质量 | 主弧差二系数 $/U$ | 补弧差二系数 $/U$ | 合计 $/U$ |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 8192 | 16 | 0.454365 | 1.265736 | −0.051957 | 1.213779 |
| 32768 | 32 | 0.456571 | 1.218649 | 0.037922 | 1.256571 |
| 131072 | 64 | 0.463781 | 1.285850 | 0.015131 | 1.300981 |

这三个区间的补弧系数有正有负。式 (13) 所要求的是它在无穷子列上高于同一个负阈值。程序可用 Python、NumPy、SciPy 和 python-flint 运行。[零点定位、特征投影及双重积分复算程序](https://github.com/lixiang90/vibemath-NumberTheory/blob/main/zeros-exact-distribution/gap-two-auxiliary-conductor-bands/verify.py)

