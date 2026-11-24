---
title: "固定偶数的哥德巴赫表示与低高度 Dirichlet 零点对"
date: "2026-11-24"
abstract: "在 Dirichlet-GRH 下，对模 N 的奇特征平均，精确求出同编号零点的贡献，给出固定偶数的带符号零点对判据、统一截高和常数的严格数值界。"
---

# 固定偶数的哥德巴赫表示与低高度 Dirichlet 零点对

本文假设 Dirichlet $L$ 函数的 GRH。下列推导使用中心线上的实际纵坐标，允许重零点、中心零点和非本原特征。文中 $\log$ 为自然对数。算术起点来自李嘉旻、刘建亚的奇特征正交公式，零点侧从 Mellin 显式公式重新推导。[李–刘原文](https://mrc.sdu.edu.cn/english/info/1018/1514.htm)

## 1. 一个可以直接检验的充分条件

取两个固定光滑窗，分别围住 $2N/5$ 和 $3N/5$。对每个偶数 $N>2$，定义 $P(N)$ 为这两个窗内**真正素数**的加权表示数，精确定义见下节。定理给出

$$
\boxed{\displaystyle \frac{P(N)}N=B_f-V_{T_N}(N)+O_A((\log N)^{-A}).} \tag{1}
$$

这里 $A>0$ 可任意固定，$B_f$ 是只依赖所选窗口的绝对常数，满足 $2\sqrt6/5<B_f<1.061$。$V_T(N)$ 是模 $N$ 全体奇特征中，同一个本原 $L$ 函数的**两个不同编号的零点**形成的实数和，公式见 (7)。所需的纵坐标高度为

$$
|\gamma|\le T_N=C_A s_N\bigl(\log(2+s_N)\bigr)^2,\qquad s_N=\log\log N, \tag{2}
$$

其中 $N$ 充分大，常数 $C_A$ 由固定窗口与精度决定。若能证明 **$V_{T_N}(N)\le0.9862$ 对每个充分大的偶数都成立**，每个这样的 $N$ 就有 $\gg N/\log^2 N$ 种平衡素数表示。这里用了第 7 节的严格数值界。更一般地，一个与 $N$ 无关的正余量 $B_f-V_{T_N}(N)\ge\eta$ 就足够。本文证明 (1)、(2) 以及这个蕴含，第 8 节给出两个具体模数的有限高度例子。逐偶数的正余量仍是待证命题。[在短区间平均后可证的正余量](goldbach-short-intervals.md)

## 2. 窗口和精确的算术恒等式

记 $c_1=2/5,\ c_2=3/5,\ b=1/20$。取一个非负、偶、积分为 1 且支撑于 $[-b,b]$ 的光滑密度 $f$，要求 $f(0)>0$。第 5 节给出满足 (2) 的具体选法。定义

$$
w_j(u)=u^{-1/2}f(\log(u/c_j))\quad(u>0),\qquad
P(N)=\sum_{p+q=N}\log p\log q\,w_1(p/N)w_2(q/N),
$$

其中 $p,q$ 为素数。两窗支撑分别在 $(1/3,1/2)$、$(1/2,2/3)$ 内，端点均离这些区间边界一段固定距离。因此每一种被计入的无序表示恰有一个次序。

设 $\Lambda$ 为 von Mangoldt 函数，$X_N^-=\{\chi\bmod N:\chi(-1)=-1\}$，$|X_N^-|=\varphi(N)/2$。令

$$
A_j(\chi,N)=\sum_{m\ge1}\Lambda(m)\chi(m)w_j(m/N),\qquad
\langle F\rangle_N=\frac2{\varphi(N)}\sum_{\chi\in X_N^-}F(\chi).
$$

展开奇特征的正交关系。在两个窗口中 $m=n$ 不可能，$m\equiv -n\bmod N$ 又等价于 $m+n=N$。故有精确等式

$$
-\Re\langle A_1\overline{A_2}\rangle_N
=\sum_{\substack{m+n=N\\(m,N)=1}}\Lambda(m)\Lambda(n)w_1(m/N)w_2(n/N). \tag{3}
$$

两端为不同素数时，$p\mid N=p+q$ 会推出 $p=q$，从而所有目标素数对自动互素。$N$ 以内指数至少为 2 的素数幂至多 $O(\sqrt N)$ 个，每个只有一个搭档，总权至多 $O_f(\sqrt N\log^2 N)$。这是把 (3) 换成 $P(N)$ 的全部代价。

## 3. 从奇特征到纵坐标

设 $\chi^*$ 为诱导 $\chi$ 的本原奇特征，导子 $q\mid N$。写其非平凡零点为 $1/2+i\gamma_r$，$r$ 遍历含重数的编号，两种符号的高度都取。对复特征，单个 $L(s,\chi^*)$ 的零点不必关于实轴成对。取 Fourier 约定 $\widehat f(t)=\int f(v)e^{itv}\,dv$，Mellin 约定 $\widetilde w(s)=\int_0^\infty w(u)u^{s-1}\,du$。令

$$
U_j(\chi,N)=\sum_r e^{i\gamma_r\log(Nc_j)}\widehat f(\gamma_r).
$$

对任一支撑在 $(1,\infty)$ 内的光滑紧支函数 $W$，奇本原特征的显式公式为

$$
\sum_m\Lambda(m)\chi^*(m)W(m)
=-\sum_r\widetilde W(1/2+i\gamma_r)-\sum_{k\ge0}\widetilde W(-1-2k). \tag{4}
$$

可从 $\Re s>1$ 上的 $-L'/L$ 作 Mellin 反演，向左移线。奇特征的平凡零点恰为 $-1,-3,\ldots$，每个留数带负号。支撑下端大于 1，使左移后的积分趋零。具体可取避开零点的水平边、在负偶整数上取左边，用函数方程及分部积分控制各边。各零点和绝对收敛。全篇用到的局部计数界为 $\#\{r:t\le\gamma_r\le t+1\}\ll\log(2q(2+|t|))$，常数对导子一致。[Koukoulopoulos，第 11 章](https://dms.umontreal.ca/~koukoulo/documents/publications/primes.pdf) · [Bennett–Martin–O’Bryant–Rechnitzer 的零点计数](https://arxiv.org/abs/2005.02989)

在 (4) 中取 $W(m)=w_j(m/N)$。非平凡零点项为 $-\sqrt N\,U_j$，平凡零点总和为 $O_f(N^{-1})$。由本原特征换回 $\chi$ 时，只涉及 $p\mid N,\ p\nmid q$。在固定倍比窗口内，每个 $p$ 只有 $O_f(1)$ 个幂次，故删去的项总计 $O_f(\sum_{p\mid N}\log p)=O_f(\log N)$。从而，逐特征一致有

$$
\frac{A_j(\chi,N)}{\sqrt N}=-U_j(\chi,N)+O_f(N^{-1/2}\log N),
\qquad |U_j(\chi,N)|\ll_f\log(2N). \tag{5}
$$

将 (5) 代回 (3)，再扣除高次幂，得到

$$
\frac{P(N)}{N}=-\Re\langle U_1\overline{U_2}\rangle_N
+O_f(N^{-1/2}\log^2 N). \tag{6}
$$

## 4. 精确求出同一零点的贡献

令 $L=\log(c_2/c_1)=\log(3/2)$，$G=f*f$。$G$ 是 $[-2b,2b]$ 上的偶概率密度。**把 (4) 再用一次**，这次令 $W(x)=x^{-1/2}G(\log x-L)$。它的支撑落在 $[e^{L-2b},e^{L+2b}]\subset(1,2)$。左侧恰为零。非平凡零点的 Mellin 值为 $e^{i\gamma_r L}\widehat f(\gamma_r)^2$，平凡零点可求和成一个几何级数。因此对每一个本原奇特征都有

$$
\sum_r e^{i\gamma_r L}\widehat f(\gamma_r)^2=-B_f,\qquad
B_f=\int_{-2b}^{2b}G(z)K(L+z)\,dz,\qquad
K(v)=\frac{e^{-3v/2}}{1-e^{-2v}}.
$$

这个等式甚至确定了左边的虚部为零。$K(v)=\sum_{k\ge0}e^{-(2k+3/2)v}$ 在 $v>0$ 严格凸。Jensen 不等式与两端点的弦分别给出

$$
\frac{2\sqrt6}{5}=K(L)<B_f
\le\frac{K(L-2b)+K(L+2b)}2<1.061.
$$

现在展开 (6)。其中 $r=s$ 的全部项已有精确值。按下式定义其余项的截断

$$
V_T(N)=\frac2{\varphi(N)}\sum_{\chi\in X_N^-}
\sum_{\substack{r\ne s\\|\gamma_r|,\,|\gamma_s|\le T}}
\widehat f(\gamma_r)\widehat f(\gamma_s)
\cos\!\left(\gamma_r\log(2N/5)-\gamma_s\log(3N/5)\right). \tag{7}
$$

每个内和的两个编号都属于同一个 $L(s,\chi^*)$，不同编号处即使纵坐标相等也予保留。因 $f$ 为实偶函数，变换为实数。带上 (6) 的负号，同一编号贡献的是正数 $B_f$，所以 $P(N)/N=B_f-V_\infty(N)+O_f(N^{-1/2}\log^2 N)$。一般的零间距分布只观测 $\gamma_r-\gamma_s$。这里还出现 $\log N$ 处的振荡，以及两个不同窗口中心的位置。

## 5. 降低所需的零点高度

令 $a_j=\kappa/[j(\log(2j))^2]$，以 $\sum_{j\ge1}a_j=b$ 确定 $\kappa>0$。令 $f$ 为独立的 $[-a_j,a_j]$ 均匀随机变量之和的密度。级数绝对收敛，特征函数为

$$
\widehat f(t)=\prod_{j\ge1}\frac{\sin(a_jt)}{a_jt}.
$$

当 $|t|$ 充分大时，前 $\lfloor c|t|/\log^2(2+|t|)\rfloor$ 个因子各自满足 $a_j|t|\ge2$（取 $c>0$ 充分小），故各自的绝对值至多 $1/2$。其余因子的绝对值至多 1，遂有

$$
|\widehat f(t)|\ll \exp[-c'|t|/\log^2(2+|t|)].
$$

Fourier 反演给出 $C^\infty$ 密度。又 $\sum_{j\ge2}1/[j\log^2(2j)]\le 1/\log 2<1/\log^2 2$，第一个均匀区间比所有剩余半宽之和还宽。因而 $f$ 在零点附近恒为 $1/(2a_1)$。这也保证短区间平均所用的积分主项严格为正。

将高度按单位区间分组，用第 3 节的零点计数，再将衰减指数中的常数略微减小，得一致上界

$$
\sum_{|\gamma_r|>T}|\widehat f(\gamma_r)|
\ll_f \log(2N(2+T))\,e^{-c''T/\log^2(2+T)}\quad(T\ge2).
$$

例如可先按 $T,2T,4T,\ldots$ 分块，每块的长度因子由该指数吸收。整个单和的绝对权重为 $O_f(\log(2N))$。至少一个高度超出 $T$ 的双和因此有同样的指数衰减，只多一个对数因子。取 (2) 中的 $C_A$ 充分大，就得到 (1)。这个截高对每个 $\chi$ 一致。模 $N$ 的奇特征总数仍为 $\varphi(N)/2$，纵坐标还要精确到能分辨 $\gamma\log N$ 的相位。

若 (7) 有固定正余量，则由 (1) 得 $P(N)\gg N$，每对素数贡献至多 $\|w_1\|_\infty\|w_2\|_\infty\log^2 N$，即得第 1 节的表示数下界。只存在一对素数时，它的归一化贡献至多为 $O_f(\log^2 N/N)$，小于 (1) 当前的误差尺度。故固定余量是推出大量表示的充分条件。要覆盖全部偶数，还须让有效起点接上有限核验。

## 6. 用完备函数的对数导数反演

对导子 $q$ 的本原奇特征写 $\Xi_\chi(s)=(q/\pi)^{(s+1)/2}\Gamma((s+1)/2)L(s,\chi)$。在 GRH 下，Hadamard 分解及函数方程给出，$a>0$ 时

$$
K_\chi(a,t)=\Re\frac{\Xi'_\chi}{\Xi_\chi}(1/2+a+it)
=\sum_r\frac{a}{a^2+(t-\gamma_r)^2}.
$$

对上述哥德巴赫窗口有精确反演

$$
U_j(\chi,N)=\frac{N^a}{\pi}\int_{\mathbb R}K_{\chi^*}(a,t)
N^{it}\widetilde w_j(1/2+a+it)\,dt. \tag{8}
$$

逐个零点验证即可。频率 $y=\log(Nu)$ 在窗口上为正，Poisson 核在该频率的变换为 $e^{i\gamma_r y}e^{-ay}$。把 $N^au^a$ 乘回后恰好抵消 $e^{-ay}$。先对单个核计算，再用快速衰减及零点计数求和，保证收敛。取 $a=1/\log N$ 时，前因子 $N^a=e$。

于是同一个 $V_T$ 的研究可以从零点表入手，也可以从同一奇特征的两个窗口积分入手。两个积分先相乘（第二个取共轭），再对 $\chi$ 平均。要由这条路线取得逐偶数进展，具体目标是控制 (7) 的带符号和。对各项分别取绝对值会丢掉所需的相位。[Rubinstein 对原方法的解析说明](https://www.math.uwaterloo.ca/~mrubinst/alpoge_claude_grh/zeta-paper_v2.pdf)

## 7. 所选窗口的常数可以严格包围

对第 5 节选定的 $f$，有 **$0.98621344<B_f<0.98621390$**。因此第 1 节的 $0.9862$ 确实留有固定余量。每个偶数原本就满足 $V_{T_N}(N)\le B_f+O_A((\log N)^{-A})$，这是 $P(N)\ge0$ 的结果。要得到大量素数表示，需要把上界压到 $B_f$ 以下。

计算只需要正项。写 $\lambda_k=2k+3/2$，则

$$
B_f=\sum_{k\ge0}e^{-\lambda_k L}
\prod_{j\ge1}\left(\frac{\sinh(\lambda_k a_j)}{\lambda_k a_j}\right)^2 .
$$

用积分判别法将 $\sum_{j>J}1/[j\log^2(2j)]$ 夹在 $1/\log(2(J+1))$ 和 $1/\log(2J)$ 之间。$j>J$ 的平方和满足 $\sum a_j^2\le\kappa^2/[J\log^4(2J)]$。由幂级数逐项比较，$1\le\sinh(x)/x\le e^{x^2/6}$。故可以同时界住每个无穷乘积的尾部。$k>K$ 的尾部再用 $G$ 的支撑及 $\sum_{k>K}e^{-\lambda_k(L-2b)}$ 作上界。以下代码用 112 比特 Arb 球算术执行了这些界，两个断言端点均向外取整。

```python
from flint import arb, ctx

ctx.prec = 112
J, K = 256, 54
b = arb(1) / 20
L = (arb(3) / 2).log()
u = [1 / (j * arb(2*j).log()**2) for j in range(1, J+1)]
mass = sum(u, arb(0))
klo = b / (mass + 1 / arb(2*J).log())
khi = b / (mass + 1 / arb(2*(J+1)).log())
tail2 = khi**2 / (J * arb(2*J).log()**4)
low = high = arb(0)
for k in range(K+1):
    t = 2*k + arb(3)/2
    ml = mh = arb(1)
    for uj in u:
        x, y = t*klo*uj, t*khi*uj
        ml *= x.sinh()/x
        mh *= y.sinh()/y
    discount = (-t*L).exp()
    low += discount*ml**2
    high += discount*mh**2 * (t**2*tail2 / 3).exp()
v = L - 2*b
high += (-(2*(K+1)+arb(3)/2)*v).exp() / (1-(-2*v).exp())
assert arb("0.98621344") < low and high < arb("0.98621390")
print(low, high)
```

## 8. 两个取不同符号的实际模数

沿用同一个 $f$ 和式 (7)，可以严格给出 **$V_{20000}(64)>B_f$，$V_{20000}(4096)<-1.83$**。这里通过算术侧及零点总数上界来包围有限零点和，无须逐个近似纵坐标。这两个例子展示了同一个测试在不同模数上如何改变符号。

令 $N=2^k\ge16$。每个奇特征的本原导子都含因子 2，诱导时没有额外 Euler 因子。记

$$
D_j=\frac1N\int_{-b}^b f(v)\,
\frac{(c_j e^v)^{-3/2}}{1-(Nc_j e^v)^{-2}}\,dv,\qquad
S_\Lambda(N)=\frac1N\sum_{\substack{m+n=N\\m\ {\rm odd}}}
\Lambda(m)\Lambda(n)w_1(m/N)w_2(n/N).
$$

式 (4) 精确给出 $\sqrt N U_j=-A_j-D_j$。对 $A_j$ 作奇特征平均，只可能留下 $m\equiv\pm1\bmod N$，这两个剩余类均在两窗之外。故 $\langle A_j\rangle_N=0$，展开后有

$$
V_\infty(N)=B_f-S_\Lambda(N)+D_1D_2/N,\qquad
(3/2)^{3/2}/N<D_j<6/N.
$$

对 $64$，下窗只含一个奇数 $25$，其搭档 $39=3\cdot13$ 的 $\Lambda$ 值为零。对 $4096$，$1637+2459=4096$ 是两素数之和。第一项距 $2N/5$ 至多 $N/300$，两窗都落在密度的平坦部分，故仅这一对便给出 $S_\Lambda(4096)\ge 200\log(1637)\log(2459)/4096>2.82$。[平坦部分的下界](goldbach-short-intervals.md)

最后控制在 $T=20000$ 处作硬截断的误差。对 $q\le4096,\ x\ge1$，Bennett 等人的推论 1.2 给出宽松的显式上界 $N_q(x)\le2x\log(8192x)$，这里计入 $|\gamma|\le x$ 的全部零点。取前 $m$ 个 sinc 因子，令 $H_m=\prod_{j=1}^m a_j^{-1}$。由 Stieltjes 分部积分，

$$
R_m(x):=2mH_m x^{1-m}
\left(\frac{\log(8192x)}{m-1}+\frac1{(m-1)^2}\right)
\ \ge\sum_{|\gamma|>x}|\widehat f(\gamma)|\qquad (m>1).
$$

用 $m=2,x=180$ 得每个特征的总绝对权重小于 7300，再用 $m=23,x=20000$ 得 $|V_{20000}-V_\infty|<0.000005$。它小于模 64 时的正修正 $D_1D_2/64$。以下球算术及素性检查接续第 7 节代码，可复算两个结论。[Bennett 等人，推论 1.2](https://arxiv.org/pdf/2005.02989#page=2)

```python
from math import isqrt

def zero_tail(m, x):
    coeff = arb(1)
    for uj in u[:m]:
        coeff /= klo * uj
    z = arb(x)
    integral = (8192*z).log()/(m-1) + arb(1)/(m-1)**2
    return 2*m * coeff * z**(1-m) * integral

whole = 360 * arb(8192*180).log() + zero_tail(2, 180)
truncation = 2 * whole * zero_tail(23, 20000)
assert whole < 7300 and truncation < arb("0.000005")
assert (arb(3)/2)**3/64**3 > truncation

center = arb(64)*2/5
assert 24 < center*(-b).exp() < 25 and 26 < center*b.exp() < 27
p, q = 1637, 2459
assert p+q == 4096 and abs(arb(p)/4096-arb(2)/5) < arb(1)/300
assert all(r % k for r in (p, q) for k in range(2, isqrt(r)+1))
pair_weight = 200 * arb(p).log() * arb(q).log() / 4096
assert high + 36/arb(4096)**3 - pair_weight + truncation < arb("-1.83")
```
