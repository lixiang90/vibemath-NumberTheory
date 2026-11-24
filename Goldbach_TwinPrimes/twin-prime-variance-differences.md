---
title: "固定差 2 的方差差分与双截断误差"
date: "2026-11-24"
abstract: "证明截断除数和的两个一阶误差在三个相关和的差分中精确抵消。在 Dirichlet-GRH 下，给出孪生素数与哥德巴赫相关和的幂次误差分解及仍需控制的带符号余项。"
---

# 固定差 2 的方差差分与双截断误差

假设普通 RH 和 Dirichlet-GRH。Goldston–Yıldırım 的截断除数和可用两个独立截点 $R,S$ 测试固定差 $2$。把两个单截断项与一个双截断项一起展开，一阶截断误差恰好抵消。下面在 $R=S=X^{1/4}$ 给出 $O_\varepsilon(X^{3/4+\varepsilon})$ 的余项分解，再将两端未评价的除数推至 $X^{1/2-o(1)}$。Gou–Wang–Xu 由比值猜想推出固定移位渐近。此处只评价能从显式算术和中分离的部分。仍需估计同一对整数的两个大除数余项的乘积。[Goldston–Yıldırım](https://arxiv.org/abs/math/0111212) · [Gou–Wang–Xu](https://arxiv.org/abs/2611.18825)

## 一　方差差分里已经可以算到什么精度

取整数 $X\ge 3$，$I_X=(X,2X]\cap\mathbb N$，$\Lambda$ 为 von Mangoldt 函数，$\mu$ 为 Möbius 函数。设

$$
A_2(X)=\sum_{n\in I_X}\Lambda(n)\Lambda(n+2),\qquad
J_X(H)=\int_X^{2X}\big|\psi(x+H)-\psi(x)-H\big|^2\,dx .
$$

沿用 Goldston–Yıldırım 的截断除数和，定义

$$
L_R(n)=\sum_{\substack{d\mid n\\d\le R}}\mu(d)\log(R/d),\qquad
E_R(n)=\Lambda(n)-L_R(n),\qquad
\mathcal E_{R,S}(X)=\sum_{n\in I_X}E_R(n)E_S(n+2).
$$

对 $n>1$，Möbius 反演给出 $E_R(n)=-\sum_{d\mid n,\ d>R}\mu(d)\log(d/R)$。特别地，当 $1<n\le R$ 时，$E_R(n)=0$。允许两端截点 $R,S$ 独立，是下面消去整段除数贡献所需要的自由度。[Goldston–Yıldırım 的截断除数和与二点相关，第 1、5 节](https://arxiv.org/html/math/0111212)

**命题 1。** 若 $2\le R,S\le X$，则对每个 $0<\varepsilon<1/2$，

$$
\begin{aligned}
 A_2(X)&=\mathfrak S(2)X+\mathcal E_{R,S}(X)\\
 &\quad+O_\varepsilon\!\left(
 X(RS)^{-1/2+\varepsilon}
 +\sqrt X(R+S)\log^2(2X)+RS\right),\\
 \mathfrak S(2)&=2\prod_{p>2}\left(1-\frac1{(p-1)^2}\right)=1.3203236\ldots .
\end{aligned} \tag{1}
$$

取 $R=S=X^{1/4}$，令 $b_n=\Lambda(n)-1$。展开三个窗口的平方，对角项与差 $1$ 的项在内点逐个消去。差分减去 $2\sum_{I_X} b_nb_{n+2}$ 后，精确留下的边界项为
$b_{2X+2}^2-b_{X+2}^2+2(b_{2X+1}b_{2X+2}-b_{X+1}b_{X+2})=O(\log^2 X)$。对中心化产生的两个一次和用 RH 的素数定理，便得到

$$
\boxed{\quad
\frac{J_X(3)-2J_X(2)+J_X(1)}2
=(\mathfrak S(2)-1)X+\mathcal E_{R,R}(X)+O_\varepsilon(X^{3/4+\varepsilon}).
\quad} \tag{2}
$$

从 $A_2$ 中删去任一端为高次素数幂的项只花 $O(\sqrt X\log^2 X)$。因此 $\mathcal E_{R,R}(X)=o(X)$ 将给出二进区间上的孪生素数渐近式。只求无穷多对，则 $\limsup_{X\to\infty}\mathcal E_{R,R}(X)/X>-\mathfrak S(2)$ 就已足够。这个严格不等式还会给出无穷多个区间各有 $\gg X/\log^2 X$ 对。现有孪生素数上界筛给出 $A_2(X)\ll X$，故 (1) 已给出双侧量级 $\mathcal E_{R,R}(X)=O(X)$。非负性只给出 $\mathcal E_{R,R}(X)\ge-\mathfrak S(2)X-O_\varepsilon(X^{3/4+\varepsilon})$，这里缺的是与该端点的正距离。[Lichtman 的孪生素数上界](https://arxiv.org/abs/2109.02851)

## 二　为什么要把三个可计算项一起估

逐项恒等式
$\Lambda\Lambda'-EE'=\Lambda L'+L\Lambda'-LL'$
中，撇号表示在另一个整数取值。前两项各只需一个算术级数中的素数定理，第三项只需数同余类。本节一并处理 $n+2$ 与 $k-n$。前者取 $k=2$，后者将用于第三节的 $k=N$。定义

$$
\begin{aligned}
m_k(R)&=\sum_{\substack{d\le R\\(d,k)=1}}
 \frac{\mu(d)}{\varphi(d)}\log(R/d),\\
m_{2,k}(R,S)&=\sum_{\substack{d\le R,\ e\le S\\(d,e)\mid k}}
 \frac{\mu(d)\mu(e)}{[d,e]}\log(R/d)\log(S/e).
\end{aligned}
$$

这里 $[d,e]$ 表示最小公倍数。对互素类，Dirichlet-GRH 给 $\psi(y;d,a)=y/\varphi(d)+O(\sqrt y\log^2(2yd))$。对非互素类，只有底数整除 $d$ 的素数幂可出现，它们在 $y\le 2X+2$ 的总权重为 $O(\omega(d)\log(2X))$。利用 $\sum_{d\le R}\log(R/d)\ll R$，两个混合项的总误差为 $O(\sqrt X(R+S)\log^2(2X))$。两个除数条件由中国剩余定理合并，恰在 $(d,e)\mid k$ 时相容，密度为 $1/[d,e]$，计整误差总共 $O(RS)$。三项合计的系数为 $m_k(R)+m_k(S)-m_{2,k}(R,S)$。同样的误差适用于端点在 $[2,2X+2]$ 内的子区间，故分部求和也允许任意固定光滑权。[Ernvall-Hytönen–Palojärvi 的 Dirichlet-GRH 误差估计](https://arxiv.org/abs/2003.01925)

**两次截断的一阶误差恰好抵消。** 对每个素数 $p$ 写 $x=p^{-1-s},\,y=p^{-1-t}$。定义

$$
H_k(s,t)=\prod_p\left(1+
 \frac{(p\mathbf1_{p\mid k}-1)xy}{(1-x)(1-y)}\right)
 =\sum_{u,v\ge1}\frac{h_k(u,v)}{u^s v^t}. \tag{3}
$$

$H_k(s,t)/(\zeta(1+s)\zeta(1+t))$ 的局部因子是 $1-x-y+p\mathbf1_{p\mid k}xy$，正好给出 $m_{2,k}$ 中不带对数权的系数。另一方面，令 $t=0$，再乘以 $1-x$，在 $p\mid k$ 时局部因子为 $1$，在 $p\nmid k$ 时为 $1-p^{-s}/(p-1)$。所以 $H_k(s,0)/\zeta(1+s)$ 正好给出 $m_k$ 对应的系数。设

$$
D(y)=\sum_{n\le y}\frac{\mu(n)}n\log(y/n)\quad(y\ge1),\qquad D(y)=0\quad(0<y<1).
$$

RH 下 $\sum_{n\le y}\mu(n)\ll_\varepsilon y^{1/2+\varepsilon}$。在 $s=1$ 用 $1/\zeta(s)$ 的值 $0$ 与导数 $1$ 作分部求和，得到 $D(y)=1+O_\varepsilon(y^{-1/2+\varepsilon})$。对 $0<y<1$ 这个误差界同样成立。[Soundararajan 对 Möbius 部分和的估计](https://arxiv.org/abs/0705.0723)

分别卷积两个变量，$m_{2,k}(R,S)=\sum h_k(u,v)D(R/u)D(S/v)$。两个单变量的式子各只保留一个 $D$。三式相加减，得出精确恒等式

$$
m_k(R)+m_k(S)-m_{2,k}(R,S)
=H_k(0,0)
-\sum_{u,v}h_k(u,v)\big(1-D(R/u)\big)\big(1-D(S/v)\big). \tag{4}
$$

最后一项要求**两个**截断误差同时出现。对 $0<\eta<1/2$，令 $\sigma=1/2-\eta$。级数 $\sum |h_k(u,v)|(uv)^\sigma$ 的局部因子在 $p\nmid k$ 时是 $1+O_\eta(p^{-1-2\eta})$，在 $p\mid k$ 时是 $1+O_\eta(p^{-2\eta})$。这是把 (3) 分母中的两个几何级数逐个取绝对值得到的。任给 $\delta>0$，第二类局部因子对充分大的 $p$ 至多为 $p^\delta$，故整个乘积 $\ll_{\eta,\delta}k^\delta$。于是 (4) 的尾项为 $O_{\eta,\delta}(k^\delta (RS)^{-1/2+\eta})$。在 $s=t=0$ 代入 (3)，偶数 $k$ 的值为

$$
H_k(0,0)=\mathfrak S(k)=
 \mathfrak S(2)\prod_{\substack{p\mid k\\p>2}}\frac{p-1}{p-2}.
$$

取 $k=2$ 即证 (1)。若将三个系数分开取绝对值，只得到 $R^{-1/2+\eta}+S^{-1/2+\eta}$ 级的截断界。取相同截点时，保留 (4) 的符号把这种估法可得的总误差 $O_\varepsilon(X^{5/6+\varepsilon})$ 改到 (2) 的 $O_\varepsilon(X^{3/4+\varepsilon})$。这里的幂次节省属于 $\Lambda L'+L\Lambda'-LL'$ 这个组合。

### 把两端的大除数都推到平方根附近

令 $r=X^{1/4}$，$B=\sqrt X\exp(-\sqrt{\log X})$。对充分大的 $X$，$r<B$。由 (1)，一致地对 $r\le R,S\le B$ 和每个固定 $A>0$ 有

$$
\mathcal E_{R,S}(X)=A_2(X)-\mathfrak S(2)X+O_A\big(X/\log^A X\big). \tag{5}
$$

因此设 $F=L_B-L_r=E_r-E_B$，则 $\sum F(n)E_B(n+2)$、$\sum E_B(n)F(n+2)$、$\sum F(n)F(n+2)$ **分别**是 $O_A(X/\log^A X)$。例如第一项是 $\mathcal E_{r,B}-\mathcal E_{B,B}$，第三项是四个不同截点组合的二阶差分。以下只须处理 $\mathcal E_{B,B}$。这一结论把 (2) 的幂次误差换成任意固定对数次节省，换来的是两个原始除数都大于 $B=X^{1/2-o(1)}$。

## 三　固定和的对应式与一项自动消失的边界贡献

令 $N$ 为偶数，$G(N)=\sum_{m+n=N}\Lambda(m)\Lambda(n)$ 按次序计数。第二节证明直接用于 $2\le n\le N-2$ 与 $N-n$，联合同余类的条件仍为 $(d,e)\mid N$。取 $R=N^{1/4}$，在 (4) 中给 $N^\delta$ 留出任意小的指数，得到对每个偶数一致的

$$
G(N)=N\mathfrak S(N)+\sum_{2\le n\le N-2}E_R(n)E_R(N-n)
 +O_\varepsilon(N^{3/4+\varepsilon}). \tag{6}
$$

也可把 $R$ 换成 $B_N=\sqrt N\exp(-\sqrt{\log N})$，误差换成 $O_A(N/\log^A N)$。区间长度与 $N$ 的有界差由误差吸收。要得逐偶数表示，只需对每个充分大偶数证明这个余项 $\ge(-\mathfrak S(N)+\delta)N$，其中 $\delta>0$ 固定。若余项为 $o(N)$，就给出 Hardy–Littlewood 主项。Gou–Wang–Xu 在比值侧使用带反射符号 $\chi(-1)$ 的 Dirichlet 特征族平均。(6) 只评价这一个固定和所选出的算术组合。[Gou–Wang–Xu，第 3 节](https://arxiv.org/html/2611.18825v1#S3)

固定和有一个额外简化。令 $C=N/2$，$U_N=L_C-L_{B_N}$。每对 $n,N-n$ 至少有一个数不超过 $C$，所以 $E_C(n)E_C(N-n)=0$ **逐项成立**。对 $C<n\le N-2$，又有 $E_C(n)=-\mu(n)\log(n/C)$。由对称性，(6) 在截点 $B_N$ 的整个余项精确地成为

$$
\sum_{n=2}^{N-2}U_N(n)U_N(N-n)
 -2\sum_{C<n\le N-2}\mu(n)\log(n/C)\,U_N(N-n). \tag{7}
$$

这两项仍需逐偶数控制。孪生素数一侧允许选取无穷多个好窗口，哥德巴赫一侧须覆盖每个充分大的偶数，之后还须覆盖有效界以下的范围。

固定差 2 也可在 ζ 的零点对上用带符号的三频带核直接选出。前者隔离算术端的大除数，后者给出零点侧的等式，两者留下的正质量问题相同。[从零点对分离差 2](../zeros-exact-distribution/zeta-zero-pairs-gap-two.md)
