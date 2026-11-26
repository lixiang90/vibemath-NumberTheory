---
title: "三阶 Hecke 互素对的整窗平滑与有限高度判据"
date: "2026-11-27"
abstract: "对 p 大于 T 的三分之二次幂的原系数双线性和，证明以 T 的 8/15 次幂为半径作局部平均，整窗绝对误差受 W 乘以任意小次幂控制。结合有限条带内的 Perron 估计，证明两个 Hecke 因子在 Y 的 7/10 次幂的 100 倍高度内满足相应无零条件，即可控制该窗口的原和。"
---

# 三阶 Hecke 互素对的整窗平滑与有限高度判据

继续取 $E=\mathbb Q(\alpha,i)$，$\alpha^3-4\alpha+1=0$。研究 $M/B$ 的三阶特征所给的两支次数二因子 $L_\nu$，$\nu=+,-$。在每个块 $T=2^jY$ 上，将原互素对的算术和写成实变量函数

$$
 Q_T(x)=\frac12\sum_{\nu=+,-}
 \sum_{\substack{h>z,\ p>T^{2/3}\ {\rm prime}\\2x/r<hp\le2x}}
 c_{z,\nu}(h)a_\nu(p)\log p . \tag{1}
$$

这里 $a_\nu,u_\nu$ 是 $L_\nu,1/L_\nu$ 的系数，$c_{z,\nu}=\mathbf1_{h=1}-u_{\nu,\le z}*a_\nu$，$z=\max(1,\lfloor T^{1/6}/\log^3(2T)\rfloor)$，$r=2^{1/\lceil\log^2(2Y)\rceil}$。记 $X=Y^{3/2}$，$W=\lceil Y^3\rceil$，$k_n=\lfloor W/n^{3/2}\rfloor$。对 $Y\le n<X$ 取含 $n$ 的块 $T\le n<2T$，并写 $\|G\|_Y=W^{-1}\sum_n k_n|G_T(n)|$，$P(G;Y)=W^{-1}\sum_n k_n[G_T(n)]_+$。一个块内平移 $x$ 时，始终用该块自己的参数。[原和及其幂指数](cm229-type-ii-large-prime.md) · [余因子的精确消失与低零点](cm229-coprime-cofactors.md)

**整窗平滑定理。** 定义归一化后的对数局部平均

$$
 Q_{T,K}^{\rm av}(x)=\frac12\int_{-1}^1 e^{-w/(2K)}Q_T(xe^{w/K})\,dw. \tag{2}
$$

取 $K=T^{7/15+\xi}$，$0\le\xi\le1/60$。对每个 $\varepsilon>0$，所有充分大的实数 $Y$ 都满足

$$
 \|Q-Q^{\rm av}\|_Y\ll_\varepsilon Y^{-\xi/2+\varepsilon}. \tag{3}
$$

$\xi=0$ 时，平均半径约为 $T^{8/15}$，未归一化的绝对误差在 $WY^{o(1)}$ 预算内。$\xi=1/60$ 时，半径约为 $T^{31/60}$，误差为 $O_\varepsilon(WY^{-1/120+\varepsilon})$。因 $|[v]_+-[w]_+|\le|v-w|$，两种尺度均可在**先平均、后取正部**的顺序下研究原问题。式 (3) 使用 Hu 预印本的 $7/8$ 无零半平面，不要求本域 ERH。[Hu，定理 1.1 与 1.2](https://hexagonmath.org/2610.00216v2)

**单窗的有限高度判据。** 固定 $0<\lambda<1/20$。若两支在下列矩形内均无零点

$$
 1/2+\lambda<\Re\rho<1,\qquad |\Im\rho|\le 100Y^{7/10}, \tag{4}
$$

则该窗口的**原和**满足

$$
 P(Q;Y)\le\|Q\|_Y\ll_{\lambda,\varepsilon}
                    Y^{15\lambda/2+\varepsilon}. \tag{5}
$$

因此，对给定的 $\eta>0$，只要本窗满足 (4)，其中 $\lambda$ 固定且充分小于 $\eta$，就有 $\sum k_n[Q_T(n)]_+\ll_\eta WY^\eta$。高度来自 $X^{7/15}=Y^{7/10}$。此前高度 $64$ 的证书已遍历两支的低零点。(4) 所要求的高度随窗口增长，目前只得到这个条件判据。[低零点证书](cm229-heights-and-zeros.md) · [全双线性和的有限高度显式公式](cm229-type-ii-zeros.md)

前篇给出 $P(Q;Y)$ 的完整窗口幂指数 $\delta=\tfrac32(\beta_*-1/2)$，其中 $\beta_*$ 是 $L_+L_-$ 的非平凡零点实部的上确界。由 (3)，下式也是无条件的

$$
 \lim_{Y\to\infty}\frac{\log(1+P(Q^{\rm av};Y))}{\log Y}=\delta .
$$

若 $\delta>0$，在 (3) 中取 $\varepsilon<\delta/2$ 即可传递下界。若 $\delta=0$，上界与非负性便给出极限。于是，局部平均以后的和仍须解决临界线右侧的零点问题。由 $\zeta_E=\zeta L(s,\chi_{-4})L_+L_-$ 及有理 Dirichlet GRH，这恰是本例 $E$ 的 ERH。[原和的等价判据](cm229-type-ii-large-prime.md) · [Alpöge–Claude 有理 Dirichlet 定理的核验范围](https://math.mit.edu/~drew/grh/)

就未经平滑的**全部正部**而言，代入 Hu 的界得到 $O_\varepsilon(WY^{9/16+\varepsilon})$。采用 Liu 的 $\beta_*\le0.8733326374459051924$ 则有 $\delta<0.56$，从而在所有充分远的完整窗口上有 $O(WY^{0.56})$。式 (3) 在高频误差上达到了任意小正幂，全和仍有上述指数差距。[Liu，定理 7.33](https://hexagonmath.org/2611.00741v1)

## 在对数变量上估计高频

固定 $T$。将 $h,p$ 各自按二进区间 $[D,2D),[P,2P)$ 切开，并与 (1) 原有的两个截断取交。保留乘积范围与 $T/2<x<4T$ 中某一短区间相交的矩形，再将每个矩形作为有限多项式使用。这样得到的函数在所需的 $x$ 及其平移处等于原和，同时在整个对数轴上有紧支撑。每个非空矩形都满足

$$
 DP\asymp T,\qquad T^{2/3}\ll P\ll T^{5/6}\log^3(2T).
$$

在临界线 $s=1/2+it$ 上，记这两个多项式为 $C_D(t)=\sum_h c_{z,\nu}(h)h^{-s}$ 和 $B_P(t)=\sum_p a_\nu(p)(\log p)p^{-s}$。该矩形对 $e^{-v/2}Q_T(e^v)$ 的 Fourier 变换，除去半和的系数，恰为

$$
 m_r(t)C_D(t)B_P(t),\qquad m_r(t)=2^s(1-r^{-s})/s,\qquad
 |m_r(t)|\ll(1+|t|)^{-1}. \tag{6}
$$

它可由对每个 $q=hp$ 积分 $\mathbf1_{\{q/2\le e^v<rq/2\}}e^{-sv}$ 直接得到。对 $|t|\le 2T^2$，置 $b=7/8$，有

$$
 |B_P(t)|\ll_\varepsilon P^{b-1/2}T^\varepsilon. \tag{7}
$$

对 $-L_\nu'/L_\nu$ 的截断级数用前篇式 (3)，先在 $1/2+\varepsilon'$ 上估计两个前缀，再在 $[P,2P]$ 分部求和移回 $1/2$。在临界线删去 $p^j$、$j\ge2$ 的代价至多 $O(\log^2(2T))$，因为 $|b_\nu(p^j)|\le2\log p$，从 $j=2$ 起可作几何求和。取 $\varepsilon'$ 充分小，Hu 的界 $\beta_*\le b$ 给出 (7)。上下端点可以在矩形内任意截断。[截断对数导数的竖线估计](cm229-type-ii-large-prime.md)

对 $k=2,3$ 和 $1\le U\le 2T^2$，Dirichlet 多项式均值定理给出

$$
 \int_{-U}^U |C_D(t)|^{2k}\,dt\ll_\varepsilon (U+D^k)T^\varepsilon,\qquad
 \int_{-U}^U |B_P(t)|^2\,dt\ll_\varepsilon (U+P)T^\varepsilon .
$$

第一式将 $C_D^k$ 看成长为 $O_k(D^k)$ 的多项式。由 $|c_{z,\nu}|\le \mathbf1_{h=1}+d_4(h)$，其卷积系数除去 $n^{-1/2}$ 后为 $O_{k,\epsilon}(n^\epsilon)$，支撑于 $[D^k,(2D)^k]$，从而系数平方和为 $O_\varepsilon(T^\varepsilon)$。第二式同理。对 $I(U)=\int_{-U}^U |C_D(t)B_P(t)|^2dt$ 用 Hölder，将 $B_P$ 的 $2k/(k-1)$ 次矩以二次矩和 (7) 界住，得到

$$
 I(U)\ll_\varepsilon T^\varepsilon(U+D^k)^{1/k}
      P^{(2b-1)/k}(U+P)^{1-1/k}. \tag{8}
$$

这里始终使用原余因子和同一支的素数迹。[Montgomery–Vaughan，Dirichlet 多项式均值推论 3](https://personal.science.psu.edu/rcv4/personal/Publications/s2-8-1-73.pdf)

置 $K_0=T^{7/15}$。在 $P\le T^{4/5}$ 用 $k=2$，在 $P>T^{4/5}$ 用 $k=3$。具体地，令 $v_P=\log P/\log T$、$\kappa=7/15$、$\theta=1-b=1/8$。因 $K_0<P$ 且 $D\asymp T/P$，(8) 在 $U=K_0$ 给出 $I(K_0)\ll T^\varepsilon K_0^2$ 所需的两个指数不等式为

$$
 \frac{v_P(k-2\theta)}{2k-1}\le\kappa,\qquad
 \frac12-\frac{v_P\theta}{k}\le\kappa. \tag{9}
$$

当 $2/3\le v_P\le4/5$，左端分别为 $7v_P/12$ 与 $1/2-v_P/16$。当 $4/5\le v_P\le5/6$，分别为 $11v_P/20$ 与 $1/2-v_P/24$。两段均满足 (9)，端点的常数及对数幂由 $T^\varepsilon$ 吸收。对 $K_0\le U\le 2T^2$，(8) 的两个括号相比 $U=K_0$ 各至多增长 $U/K_0$ 倍，其指数之和为 $1$。于是 $I(U)\ll_\varepsilon T^\varepsilon U K_0$。对更高频率，直接将 $C_DB_P$ 看成长度 $O(T)$ 的多项式，其归一化系数平方和也是 $O_\varepsilon(T^\varepsilon)$，给出 $I(U)\ll_\varepsilon (U+T)T^\varepsilon$。

对 $K_0\le K\le T^{7/15+1/60}$，用 (6) 在 $K\le |t|\le T^2$ 按二进区间求和，再加上 $|t|>T^2$ 的尾部，有

$$
 \int_{|t|\ge K}|m_r(t)C_D(t)B_P(t)|^2dt
 \ll_\varepsilon T^\varepsilon K_0/K . \tag{10}
$$

(2) 在变量 $v=\log x$ 中是对 $e^{-v/2}Q_T(e^v)$ 作宽 $2/K$ 的均匀平均。差的乘子为 $1-\sin(t/K)/(t/K)$，当 $|t|\le K$ 时，它与 $m_r(t)$ 的乘积模长至多 $C/K$。这段频率的平方积分至多 $C I(K)/K^2$，也满足 (10) 的右端。由 Plancherel，卷积差在对数轴上的二范数至多 $T^\varepsilon(K_0/K)^{1/2}$。对多对数个矩形相加，再在长度有界的对数块上用 Cauchy–Schwarz，即得相同的块积分上界。

最后从积分回到整数。移动 $x$ 不超过两个单位时，(1) 两端各跨过 $O(1)$ 个乘积下标，每个下标的系数为 $O_\varepsilon(T^\varepsilon)$。对 (2) 取同一平移作平均，这个变动界仍成立。在每个整数 $n$ 旁取 $[n,n+1]$，用 $k_n/W\le n^{-3/2}$，从而块和至多是常数乘相应的 $\int |Q_T-Q^{\rm av}_{T,K}|x^{-3/2}dx$ 加 $O_\varepsilon(T^{-1/2+\varepsilon})$。积分用上述对数块估计即可，最后一个不满块同样适用。对 $T$ 求和、缩小中间使用的 $\varepsilon$，得到 (3)。若 $X$ 为整数并把 $n=X$ 也计入，单个端点的归一化费用至多 $O_\varepsilon(X^{-1/2+\varepsilon})$。

同一计算可以直接代入其他无零界。设 $4/5\le b\le7/8$ 且已知 $\beta_*\le b$，则将矩的分界改取 $v_P=3/(2(1+b))$，可取 $\kappa=b/(1+b)$。在该分界，(9) 中 $k=2$ 的第一项与 $k=3$ 的第二项相等，其余两个外端点的限制均为 $(1+2b)/6\le b/(1+b)$。于是相应的平均半径指数为 $1/(1+b)$，下节的单窗高度指数为 $3b/(2(1+b))$。例如采用 Liu 预印本给出的上界 $b=0.8733326374459051924$，高度指数可取 $0.699288$。[Liu，定理 7.33](https://hexagonmath.org/2611.00741v1)

在原系数上，剩余问题也可以完全用有限和表述。保留上述矩形，令 $F_T(t)=\tfrac12\sum_{\nu,D,P}C_{D,\nu}(t)B_{P,\nu}(t)$。由 (6)、(10) 与 $|m_r(t)|\ll\min(\log r,(1+|t|)^{-1})$，若能对每个 $\epsilon>0$ 证明

$$
 \int_{|t|<T^{7/15}}\min\big((\log r)^2,(1+|t|)^{-2}\big)
                   |F_T(t)|^2dt\ll_\epsilon T^\epsilon, \tag{10a}
$$

则原和的整窗绝对值与正部都达到 $WY^{o(1)}$ 的未归一化预算。(10a) 允许两支特征、不同余因子和不同大素数在取模方前抵消。直接证明 (10a) 仍是开放的问题。有限高度无零信息给出下节的定量界。

## 低频与有限高度判据的证明

固定 (4) 中的 $\lambda$。对一块记 $K_0=T^{7/15}$。在 $\Re s\ge1/2+2\lambda$、$|\Im s|\le50K_0$ 的固定宽度条带内，局部对数导数公式和单位高度的 $O(\log(2K_0))$ 零点数给出 $L_\nu'/L_\nu\ll_\lambda\log(2K_0)$。在该矩形内，水平积分至绝对收敛线 $\Re s=2$ 给出一致的解析分支 $\log L_\nu(s)\ll_\lambda\log(2K_0)$。对 $\log L_\nu(s)e^{(s-it_0)^2}$ 在左右边 $1/2+2\lambda,2$ 之间插值。右边为 $O(1)$，两条水平边在 $|t_0|\le20K_0$ 时有指数衰减，故在 $\Re s\ge1/2+3\lambda$ 得 $|\log L_\nu(s)|\ll_\lambda (\log(2K_0))^{1-c_\lambda}$，其中 $c_\lambda>0$。这也可直接在矩形上除以 $(\log(2K_0))^{(2-s)/(2-(1/2+2\lambda))}$ 后用最大模原理证明。因而 $L_\nu^{\pm1}\ll_{\lambda,\varepsilon}T^\varepsilon$ 在缩小的矩形内成立。右侧绝对收敛区用欧拉乘积即可。[本例的局部零点计数和对数导数](cm229-critical-logarithm.md)

以下 Perron 截断的高度取 $K_0$。对系数 $f=u_\nu,b_\nu,c_{z,\nu}$，其绝对值皆有与 $z$ 无关的除数函数及对数幂上界。端点取 $D'=\lfloor D\rfloor+1/2$，$D\ge1$。在 $1/2\le\sigma\le1$，初始线 $\Re w=1-\sigma+1/\log(2T)$ 上，截断 Perron 的误差至多

$$
 O_\varepsilon\!\left(T^\varepsilon
          \left(D^{1-\sigma}/K_0+D^{-\sigma}\right)\right). \tag{11}
$$

事实上，对 $D'/2<j<2D'$ 按 $|j-D'|$ 求和，用 $\min(1,CD'/(K_0|j-D'|))$。对其余项，绝对值级数在 $1+1/\log(2T)$ 仅花对数幂。本文的端点至多 $O(T^{5/6}\log^3(2T))$，其平方根为 $o(K_0)$。移线时，两条水平边的 $D'^{\Re w}/K_0$ 因子也有相同余量。

先对 $1/L_\nu(s+w)$ 用 (11)，得到

$$
 |U_{z,\nu}(s)|\ll_{\lambda,\varepsilon}T^\varepsilon
 \quad(\Re s\ge1/2+4\lambda,\ |\Im s|\le3K_0).
$$

对 $\Re s<1$ 移至 $\Re w=1/\log(2T)$，竖线积分的 $1/|w|$ 只付对数，经过的点均在上述局部矩形内。$\Re s\ge1$ 时可直接用绝对系数截到 $z$。由此 $1-U_{z,\nu}L_\nu$ 与 $-L_\nu'/L_\nu$ 在 $\Re s\ge1/2+4\lambda$、$|\Im s|\le 2K_0$ 都为 $O_{\lambda,\varepsilon}(T^\varepsilon)$。回到 $s=1/2+it$、$|t|\le K_0$，对这两个级数分别用 (11) 并将积分移至 $\Re w=5\lambda$，全程在 $w=0$ 右方。相减两个前缀，并从第二个级数扣去素数幂，便得到

$$
 |C_D(t)|\ll_{\lambda,\varepsilon}D^{5\lambda}T^\varepsilon,\qquad
 |B_P(t)|\ll_{\lambda,\varepsilon}P^{5\lambda}T^\varepsilon.
$$

由于 $DP\asymp T$，低频部分在 (6) 中的二范数至多为 $T^{5\lambda+\varepsilon}$。对高频使用无条件的 (10)，取 $K=K_0$。按上节相同的 Plancherel、离散化及整窗求和，得到 (5)。

## 在同一个完整窗口复算平滑损失

取 $Y=4096$，用 $160$ 位实球算术遍历 $4096\le n<262144$，$K=T^{7/15}$。表中各数是严格下端点，上端点为所列数加 $10^{-6}$。原截断在六个块上均为 $1$。第二行专用 $z=\lfloor T^{1/6}\rfloor$ 检验 $z>1$ 时的卷积，六块的截断依次为 $4,4,5,5,6,7$。

| 截断 | $P(Q;Y)$ | $P(Q^{\rm av};Y)$ | $\lVert Q-Q^{\rm av}\rVert_Y$ |
| --- | ---: | ---: | ---: |
| 原式 | 0.532182 | 0.339071 | 0.655849 |
| $\lfloor T^{1/6}\rfloor$ | 0.390309 | 0.261359 | 0.480923 |

计算利用 $S(v)=\sum_{q\le v}A_q$ 与 $D(v)=\sum_{q\le v}A_q/\sqrt q$，其中 $A_q$ 是 (1) 的乘积系数。在同一块上，$Q_T(x)x^{-3/2}$ 的一个连续原函数为

$$
 J_T(x)=2\sqrt2\big(D(2x)-r^{-1/2}D(2x/r)\big)-2x^{-1/2}Q_T(x).
$$

每个乘积 $q$ 在 $q/2$ 与 $rq/2$ 的跳跃分别抵消，区间内部求导即可验证。于是 (2) 等于 $\tfrac12 K\sqrt n\,[J_T(ne^{1/K})-J_T(ne^{-1/K})]$。在每块的首、中、末三点，程序另外逐个积分各乘积的指示区间，两种算法的实球均相交。未平滑的两行同时复现了前篇用另一方向的卷积截断算出的正部。数值检验的是权重、端点与局部平均。关于任意远窗口的结论来自 (3) 与 (5) 的证明。[复算程序](verify_cm229_coprime_smoothing.py)
