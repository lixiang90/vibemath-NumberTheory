---
title: "判别式 7168 的四次 CM 域与判别式 229 三次域的 CM 扩张"
date: "2026-11-27"
abstract: "对一个非正规四次 CM 域、它的反射域和正规闭包，算出共同的二面体 Hecke 因子，并认证正高度 256 以内的全部零点，计重分别为 804、855、1822。与此前六次例在同一截断高度和素数范围上比较，得到全尺度显式公式界、相同的反射域相对高度，以及主项为 x/32 的整窗正部判据。"
---

# 判别式 7168 的四次 CM 域与判别式 229 三次域的 CM 扩张

取 $F=\mathbb Q(\sqrt2)$，$K=F(\sqrt{-(3+\sqrt2)})$。这个四次域有一个与自身不同的反射域 $K^r$。两域的全高度 ERH 完全等价，准 ERH 的最右零点边界也相同，因为两者中唯一超出有理 Dirichlet 定理的因子是同一个 $L$ 函数。与此前 $F_3=\mathbb Q[u]/(u^3-4u+1)$、$E=F_3(i)$ 的比较统一取高度 $256$，以及有理素数上限 $10^8$。[六次域的分解与零点](cm229-qrh-boundary.md) · [有理 Dirichlet 定理的独立核验](https://math.mit.edu/~drew/grh/)

| 量 | 新四次域 $K$ | 此前的六次域 $E$ |
| --- | ---: | ---: |
| 实子域判别式 | $8$ | $229$ |
| CM 域绝对判别式 | $7168$ | $3356224$ |
| 正规闭包及其次数 | $D_4$，$8$ | $S_3\times C_2$，$12$ |
| 非阿贝尔因子的导子 | $896$ | $229,\ 3664$ |
| 实子域在 $0<\Im\rho\le256$ 的零点数 | $306$ | $552$ |
| $0<\Im\rho\le256$ 的非平凡零点 | $804$，均临界且单 | $1275$，均临界且单 |
| 扣除低零点后 $X^{b-1}$ 项的上界系数 | $0.020155$ | $0.040372$ |
| 双惰性平均的主项 | $x/32$ | $x/48$ |

新例的实子域是二次域，其全高度 ERH 已由有理 Dirichlet 定理给出。旧例的实子域本身含有导子 $229$ 的非阿贝尔因子。表中 $b=\theta_G=0.8733326374459050661\ldots$ 取自 Liu 关于一般 Hecke 函数的预印本。两个 CM 域在所有高度都可应用 $\Re s>b$ 的无零半平面。本文的计算把新例在 $|\Im s|\le256$ 的边界推到 $1/2$。在更高处，把固定常数降到 $0.87333$ 还缺一个明确的素理想和估计。[Liu 定理 7.33](https://hexagonmath.org/2611.00741v1) · [临界常数区间](hecke-qrh-contact.md) · [两例的七次幂自由理想和判据](cm-powerfree-boundary.md)

## 整环、反射域与共同的非阿贝尔因子

在复数中取 $\alpha=i\sqrt{3+\sqrt2}$、$\beta=i\sqrt{3-\sqrt2}$，故 $\alpha\beta=-\sqrt7$。$\alpha$ 的极小多项式为 $P(T)=T^4+6T^2+7$，两个实嵌入下 $-(3+\sqrt2)$ 均为负。在 $F$ 上，$3+\sqrt2$ 的范数为 $7$，其生成的素理想处相对方程为 Eisenstein。在 $\pi=\sqrt2$ 处，$\alpha-1$ 满足 $T^2+2T+4+\sqrt2=0$，同样为 Eisenstein。其余素理想处相对幂基的判别式 $-4(3+\sqrt2)$ 是单位。于是 $\mathcal O_K=\mathcal O_F[\alpha]=\mathbb Z[\alpha]$，相对判别式范数为 $16\cdot7=112$，$D_K=8^2\cdot112=7168$。

正规闭包为 $M=K(\sqrt7)$。因 $K\cap\mathbb R=F$ 且 $\sqrt7\notin F$，其次数为八。对 $\alpha,\beta$ 独立变号并交换，给出整个 $D_4$。同时变号是中心元 $c$，也是复共轭。令 $r$ 将 $(\alpha,\beta)$ 送到 $(\beta,-\alpha)$。则

$$
 M^{\langle c\rangle}=B=\mathbb Q(\sqrt2,\sqrt7),\qquad
 M^{\langle r\rangle}=k=\mathbb Q(\sqrt{14}).
$$

四种 CM 型各从 $\{\alpha,-\alpha\}$、$\{\beta,-\beta\}$ 中选一根，构成一个伽罗瓦轨道。选 $\alpha,\beta$ 的型由交换两根所稳定。型迹 $\alpha+\beta$ 的平方是 $-6-2\sqrt7$，它生成该稳定子的四次不动域。因此 $K^r=\mathbb Q(\sqrt{-(6+2\sqrt7)})$，实子域为 $\mathbb Q(\sqrt7)$。$K$ 仅有实二次子域 $F$，故所有型均为本原型。[Streng，四次 CM 型与反射域，第 7 节](https://www.sagemath.org/files/thesis/streng-thesis-2010.pdf)

设 $\lambda$ 为 $M/k$ 的一个忠实四阶 Hecke 特征，$\eta$ 为 $K/F$ 的二次特征，$\rho$ 为 $D_4$ 在两根所张成平面上的表示。诱导给出 $\rho=\operatorname{Ind}_{\langle r\rangle}^{D_4}\lambda$。$K$ 的四点置换表示在 $\pm$ 根的两对上取和，得到 $1\oplus\chi_8$，在每对上取差，得到 $\rho$。反射域同理给 $1\oplus\chi_{28}\oplus\rho$。用 $L_D=L(s,\chi_D)$ 表示基本判别式 $D$ 的有理因子，则包括所有分歧处都有

$$
 \begin{aligned}
 H(s)&=L_k(s,\lambda)=L_k(s,\bar\lambda)=L_F(s,\eta)=L(s,\rho),\\
 \zeta_K&=\zeta L_8 H,\qquad \zeta_{K^r}=\zeta L_{28}H,\\
 \zeta_M&=\zeta L_8L_{28}L_{56}H^2.
 \end{aligned} \tag{1}
$$

式中 $\bar\lambda$ 与 $\lambda$ 在 $k/\mathbb Q$ 的自同构下互换。$H$ 整且自对偶，有限导子由判别式给出 $q=7168/8=896$。它的完备因子为 $896^{s/2}\Gamma_{\mathbb R}(s+1)^2$，函数方程符号为 $+1$，可对 $K,F$ 的完备函数作商。导子判别式公式还给 $D_{K^r}=28\cdot896=25088$、$D_M=8\cdot28\cdot56\cdot896^2=10070523904$。对每个 $1/2\le d<1$，(1) 中三个数域的 $\zeta$ 在 $\Re s>d$ 无零都等价于 $H$ 在此半平面无零。陈述无零性时均扣除 $s=1$ 的极点。[Artin 形式与类域互反](https://www.jmilne.org/math/CourseNotes/CFT.pdf)

## 高度 256 的完整计数

对 $p\nmid14$，令 $a=(2/p)$、$g=(7/p)$。写 $P_p(X)=\det(1-X\rho(\mathrm{Frob}_p))$。若有一个符号为负，迹为零。若 $a=g=1$，取 $t^2=2\bmod p$，迹为 $2(-(3+t)/p)$。两种 $t$ 给同一符号，因为 $(3+t)(3-t)=7\bmod p$。其中迹 $+2$ 对应完全分裂，$-2$ 对应 $c$。行列式总为 $ag$。在分歧处 $P_2(X)=1$、$P_7(X)=1-X$，由 $P\bmod p$ 的不同不可约因子与实二次域的 Euler 因子作商可得。程序对 $p\le65536$ 逐一比较两种局部计算。

| 因子 | $\zeta$ | $L_8$ | $L_{28}$ | $L_{56}$ | $H$ |
| --- | ---: | ---: | ---: | ---: | ---: |
| 正高度至 $256$ 的零点数 | $111$ | $195$ | $246$ | $274$ | $498$ |

$H$ 的首末正零点高度分别落在 $(5039604,5039605)/2^{22}$ 与 $(1073091971,1073091972)/2^{22}$。五支全部区间宽 $2^{-22}$，跨支或同支的相邻区间端点至少相隔 $1111/2^{22}$。因此 $\zeta_K,\zeta_{K^r}$ 分别有 $804,855$ 个单零点。两者恰好共享 $\zeta$ 和 $H$ 的 $609$ 个正高度零点。$\zeta_M$ 在同一范围有 $1324$ 个不同零点，计重 $1822$，其中 $498$ 个恰为二重。完整计数还排除了实区间 $(0,1)$ 内的零点。[全部隔离区间](cm7168-zeros256.csv) · [核验程序](verify_cm7168.py)

有理 $\zeta$ 的零点按序号由 FLINT 隔离，第 $111,112$ 个分居高度 $256$ 两侧。三个二次特征用直接的复区间 $L$ 值计算。对 $H$，在 $\Re s=176$ 用前 $12288$ 个 Dirichlet 系数，$|a(n)|\le d_2(n)$ 控制余项。取傅里叶周期 $1024$、步长 $1/8$、采样半径 $1792$，保留 $2000$ 个正频率。用 $2688$ 位 Arb 算术在顶边以 $1/128$ 的步长、$32$ 次插值和 Cauchy 余项逐段排除原点，顶边最小坐标余量大于 $2.78\cdot10^{-173}$。在 $\Re s=2$ 以绝对收敛 Euler 乘积固定辐角，由函数方程及共轭对称补齐矩形。辐角原理给 $H$ 的总重数 $996$，恰被正负高度的临界线异号区间用尽。这一步与表中的有理因子计数均独立于全高度的准 ERH 输入。[傅里叶余项与围道推导](tetra3136-prime-bias.md) · [FLINT 的有序零点算法](https://flintlib.org/doc/algorithms/special-functions/zeta-zeros.html)

## 全尺度误差与 CM 高度

令 $\mathcal Z_J(X)=\sum_\rho X^{\rho-1}/(\rho(\rho-1))$，遍历 $\zeta_J$ 的全部非平凡零点并计重。对所有 $X\ge2$ 有

$$
 \begin{aligned}
 |\mathcal Z_K(X)|&\le 2.492753 X^{-1/2}+0.020155 X^{b-1},\\
 |\mathcal Z_{K^r}(X)|&\le 2.830066 X^{-1/2}+0.020155 X^{b-1},\\
 |\mathcal Z_E(X)|&\le 4.70862 X^{-1/2}+0.040372 X^{b-1}.
 \end{aligned} \tag{2}
$$

最后一行取自此前的六次例。新四次域在 $1024\le X\le\exp(1350000)$ 已有 $|\mathcal Z_K(X)|<0.897 X^{-0.12667}$。具体地，置 $d=0.87333$，用有向舍入计算 $2.492753\cdot1024^{1/2-d}+0.020155\exp(1350000(b-d))<0.897$。这个更长的区间是素理想截断尺度 $X$ 的范围。全高度边界的改进要求让 $X$ 无界。[六次例相同正规化的界](cm229-qrh-boundary.md)

为得 (2)，记 $\ell=H'(1)/H(1)$，$u=\sum_{0<\gamma\le256}2/(1/4+\gamma^2)$ 仅对 $H$ 求和。Hadamard 乘积及函数方程给出 $H$ 的总质量

$$
 M_H=\sum_\rho\Re\frac{1}{\rho(1-\rho)}
 =2\ell+\log896-2\log\pi+2\psi(1)
$$

其中 $\psi=\Gamma'/\Gamma$。由 $1/8<1-b\le\Re\rho\le b<7/8$，其余零点的绝对权之和至多为 $\sqrt{1+(3/(4\cdot256))^2}(M_H-u)<0.020155$。有理因子的全部零点都按 $X^{-1/2}$ 加权。其总质量为 $M_\zeta=2+\gamma-\log(4\pi)$ 及 $M_D=2L_D'(1)/L_D(1)+\log(D/\pi)+\psi(1/2)$，其中 $D=8,28$。程序算得 $M_\zeta+M_8+u<2.492753$、$M_\zeta+M_{28}+u<2.830066$。两域相减时 $H$ 全部消去，故对任意 $X\ge2$ 还有 $|\mathcal Z_K(X)-\mathcal Z_{K^r}(X)|<0.808178 X^{-1/2}$。[常数与三处素数截断的独立交叉检验](verify_cm7168.py)

取 Andreatta–Goren–Howard–Madapusi Pera 平均 Colmez 定理的高度正规化。每个域内的四种型均共轭，故各型高度均等于本域的平均值。两域的相对函数、导子和次数也都相同，所以每个型有同一个精确高度

$$
 \bar h_K=\bar h_{K^r}
 =\tfrac14\log896-\gamma-2\log(2\pi)+\tfrac12\ell
 =-2.8341764746916777079\ldots .
$$

这里 $\ell=-0.56138356043037299592\ldots$。由相对显式公式，两域的有限素理想相对高度逼近 $\widehat h_{\rm rel}(X)$ 也对每个实数 $X>1$ 逐点相等，因为按有理范数汇集的系数均由 $-H'/H$ 给出。此前六次域则有诱导型和混合型两个轨道，其高度相差 $1.13235559\ldots$。这给出了使用 CM 型比较两例时一个全尺度的区别。[平均 Colmez 定理 A](https://arxiv.org/abs/1508.00178) · [有限素理想高度公式](cm-qrh-height-formula.md) · [六次域两种型的高度](cm229-heights-and-zeros.md)

## 完全分裂素数与整窗正部

令 $\Theta_+(x),\Theta_-(x)$ 分别为 $p\le x$、$\mathrm{Frob}_p=1,c$ 的未分歧素数的 $\log p$ 之和。两者主项均为 $x/8$，其和来自在 $B$ 完全分裂的素数，有 $x/4+O(\sqrt x\log^2(2x))$。令 $\Psi_M$ 为数域 von Mangoldt 函数的累积和。对每个未分歧的 $p$，$p^j$ 的系数为 $8\log p$ 当且仅当 $\mathrm{Frob}_p^j=1$。满足 $g^2=1$ 的元素恰为中心的两个元素和四个反射，即 $G/\langle c\rangle$ 中四类的三类。因这个商为双二次群，有理 Dirichlet 定理把对应的平方项计至精度 $O(x^{1/4}\log^2(2x))$。分歧项为 $O(\log(2x))$，三次以上素数幂给 $O(x^{1/3}\log^2(2x))$，故

$$
 \Theta_+(x)-x/8
 =\frac{\Psi_M(x)-x}{8}-\frac34\sqrt x+O(x^{1/3}\log^2(2x)). \tag{3}
$$

在先前的双惰性素数平均 $A_Y$ 中，把固定素数集合取为此处的完全分裂素数。取相同的 $0<e<1/52$、$z=Y^{e/5}$。原转移给 $4A_Y(x)=\Theta_+(x)-\Theta_+(z)+O_N(Y^{-N})$，一致于整个 $[Y,Y^{3/2}]$。因此主项是 $x/32$，平方项为 $-3\sqrt x/16$。六次例的相应两项为 $x/48$ 与 $-\sqrt x/6$。[双惰性平均的定义与转移定理](winpair-inert-transfer.md) · [六次例的平方项](cm229-one-sided-primes.md)

设 $\beta_K$ 为 $H$ 的非平凡零点实部上确界，置 $P_K(Y)=\int_Y^{Y^{3/2}}[A_Y(x)-x/32]_+\,dx/x^{3/2}$。则

$$
 \lim_{Y\to\infty}\frac{\log(1+P_K(Y))}{\log Y}
 =\frac32(\beta_K-\tfrac12). \tag{4}
$$

具体地，写 $T(x)=\sum_{n\le x}b(n)$，其中 $\sum b(n)n^{-s}=-H'/H(s)$。有理素数处的迹给 $T=2(\Theta_+-\Theta_-)+O(\sqrt x\log^2(2x))$。用 $\Theta_++\Theta_-$ 的估计以及转移，有 $T(x)=16(A_Y(x)-x/32)+O(\sqrt x\log^2(2x))$。在 $\Re s>\beta_K+\epsilon$ 对对数导数作截断 Perron 移线，给 $T(x)=O_\epsilon(x^{\beta_K+\epsilon})$，从而得到 (4) 的上界。

若 $1/2<\theta<\beta_K$，则 $T(e^v)$ 的 Laplace 变换 $-H'(s)/(sH(s))$ 在 $\Re s>\theta$ 有极点，且在正实轴没有极点，后者由完整零点计数及 $H(1)\ne0$ 得到。固定 $1<a<3/2$，原窗口引理给出对每个充分大的 $Y$ 都有 $\int_{Y^a}^{Y^{3/2}}[T(x)]_+\,dx/x^{1+\theta}>1$，否则无界反例序列将排除这个极点。移到 $A_Y-x/32$ 的积分误差趋于零。乘 $Y^{a(\theta-1/2)}$ 并分别令 $a\uparrow3/2,\theta\uparrow\beta_K$，就得 (4) 的下界。当 $\beta_K=1/2$，下界由正性直接得到。这也证明新域的 ERH 等价于对每个 $\eta>0$ 在任意远的完整窗口上取得 $P_K(Y)\le Y^\eta$。[原窗口引理及六次例的末段论证](cm229-binomial-slopes.md) · [式 (4) 对应的单变量理想和](cm-powerfree-boundary.md)

最后用同一区间 $2\le x\le10^8$ 比较 $J(x)=(\Theta_+(x)-x/n)/\sqrt x$，其中 $n$ 为相应正规闭包次数。

| 素数统计 | 新例 $n=8$ | 六次例 $n=12$ |
| --- | ---: | ---: |
| 完全分裂素数个数 | $719724$ | $479605$ |
| 首个和最后一个跳后高于主项的素数 | $47,\ 28806233$ | $37,\ 7207729$ |
| 最后一次回落至主项的实数位置 | $28806382.557387\ldots$ | $7207871.987665\ldots$ |
| $\max J$（取到处） | $0.37236742\ldots$（$9791$） | $0.46195081\ldots$（$873157$） |
| $\inf J$（从左趋于） | $-2.14203573\ldots$（$91009$） | $-2.15418485\ldots$（$55511413$） |
| $J(10^8)$ | $-0.89838774\ldots$ | $-0.93830924\ldots$ |

在新例的 $1439448$ 个 $B$ 中完全分裂的候选素数上，独立比较二次剩余判据与 $T^p=T\bmod(p,P(T))$。两次跳跃之间 $J$ 严格递减，用 $128$ 位 Arb 检查每次跳跃及端点即可确定全范围的上下确界。表中最后一次回落都仅指核验范围内的值。全高度尚须证明 (4) 的右端为零，或为准 ERH 目标独立压低该指数。[新例素数核验](verify_cm7168_primes.py) · [六次例素数核验](verify_cm229_split_primes.py)

在 ERH 目录下用 Python 与 python-flint 0.9.0 复算。第二条从欧拉系数重建零点表。

~~~sh
python verify_cm7168.py
python verify_cm7168.py --regenerate
python verify_cm7168_primes.py
~~~
