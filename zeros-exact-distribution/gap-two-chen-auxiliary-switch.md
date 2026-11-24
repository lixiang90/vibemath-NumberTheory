---
title: "辅助素数换号后的零点矩与陈氏半素数"
date: "2026-11-25"
abstract: "在 Alpöge–Claude 两阶段抽取辅助素数时翻转一个端点的二次符号。两族 Dirichlet 零点给出的二次和，其二阶矩以系数 1/128 测出陈氏定理留下的半素数。若半素数占正比例，四个因子符号类的波动联合趋向独立正态。由此得到真孪生素数占正比例的等价概率判据，并算出只在部分共同好窗口上检验时所需的均方下降量。实际的严格下降仍待证明。"
---

# 辅助素数换号后的零点矩与陈氏半素数

前三种思路中，从 ζ 零点逐个恢复素数能给出精确位置，用 ζ 零点纵坐标差能检测非常稀疏的孪生素数，固定剩余类的 Dirichlet 显式公式则容许把 $p+2$ 分解并在算术级数内求平均。以改进下界为目标，我选择第三条。它已有陈氏定理保证的候选数量，只需衡量其中多少是半素数。下文给出原选窗中这部分质量对应的波动尺度和概率门槛。[ζ 零点对的固定差判据](https://github.com/lixiang90/vibemath-NumberTheory/blob/main/zeros-exact-distribution/zeta-zero-pairs-gap-two.md) · [固定剩余类与互补因子](https://github.com/lixiang90/vibemath-NumberTheory/blob/main/zeros-exact-distribution/gap-two-complementary-divisors.md)

设 $p$ 是素数，$x=p+2$。陈氏定理保证在每个大二倍区间内，有充分多的 $p$ 使 $x$ 为素数或两个大素数之积。前篇把剩余半素数写为 $x=rs$，其较小素因子满足 $x^a\le r\le x^b$，其中可固定 $3/11<a<b<1/2$ 而保住下界的数量级。直接估计这些半素数时，同一个 $r$ 可出现许多次。本文在 Alpöge–Claude 原来的辅助素数 $Q,D$ 上保留 $r$ 的符号，另在较小端点 $p$ 处翻转 $D$ 的符号。这样作差后，每行都有一个只属于该行的素数 $p$。[前篇，第 5 节](https://lixiang90.github.io/vibemath-NumberTheory/zeros-exact-distribution/articles/gap-two-separated-zero-windows.html) · [Green–Tao，定理 6.1 与脚注 3](https://www.numdam.org/article/JTNB_2006__18_1_147_0.pdf#page=20)

写 $\mu=(\log U)^2/U$，$T(U)$ 为 $U<p\le2U$ 中真孪生素数对数，$B(U)$ 为上述内部半素数的个数。我们构造一个实的零点二次和 $\Delta_\sigma$，其中 $\sigma\in\{\pm1\}^2$ 指定 $r$ 在 $Q,D$ 两处的符号。对每个固定 $\sigma$，

$$
\mathbb E_{Q,D}\Delta_\sigma^2
=\frac{\mu^2}{128}B(U)+O\!\left(U^{-1-a/2}(\log U)^5\right).
\tag{1}
$$

若沿一列尺度有 $\mu B(U)\ge c>0$，则 $\Delta_\sigma/(\mu\sqrt{B(U)/128})$ 按原来的辅助素数抽样收敛到实标准正态分布。因而，若反复有正比例的辅助对使 $|\Delta_\sigma|=o(\sqrt\mu)$（具体量词见第 5 节），则这些尺度上 $T(U)\gg U/\log^2U$。这个判据给出应当在零点上检验的量。本文证明其分布和推论，尚未证明所需的小值事件。

## 1. 原两阶段抽样的算术输入

固定 ζ 目标、$0<\epsilon<1/52$ 和原选窗所用的幂指数 $J$。令基准 $P\to\infty$，外端 $H=P^J$，$z=P^{\epsilon/5}$，$z'=z(\log P)^2$。先在 $e^{z/2}<Q\le e^z$ 且 $Q\equiv3\pmod8$ 的素数中均匀抽取。给定 $Q$ 后，在 $e^{z'/2}<D\le e^{z'}$ 且 $D\equiv7\pmod {16}$、$D\equiv-1$ 模 $Q$ 及所有奇素数 $\ell\le z$ 的素数中均匀抽取。记 $\lambda_V(n)=(n/V)$。下文对确定的 $P\le U,\ 2U\le H$ 作估计，各常数可依赖固定的 $\epsilon,J,a,b$。原第二阶段的可容许性验证使用参数 $6\epsilon/5$，其候选集对每个第一阶段候选 $Q$ 均非空。[原候选定义](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/MeanSquareBounds/PrimeSetTwoBounds.lean#L37) · [第二阶段及模数验证](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/MeanSquareBounds/CharOf.lean#L56)

归档的定理 $\mathrm{QH.twistK}$ 有两个任意的固定指数 $A,C>0$。在第一族，或在先固定 $Q$ 的第二族，记当前素数为 $V$、相应基模为 $m$。除至多一个指标 $e_*$ 外，对奇的平方自由数 $1<u\le P^A$、$(u,m)=1$，有

$$
|\mathbb E_V\lambda_V(u)|\le P^{-C}.
\tag{2}
$$

这是对原定理使用二次互反律，因各族中 $V\bmod4$ 固定。第二族的 $e_*$ 可以依赖 $Q$。当 $u$ 含有第二基模中的小素因子时，这些因子的符号由 $D$ 的固定剩余类决定，先取出即可。本文出现的物理端点、因子和固定阶数的乘积都小于辅助素数，或在使用 (2) 时小于 $P$ 的某个固定次幂。每一阶矩可各选一次 $A,C$。[原定理 twistK 的量词](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/MeanSquareBounds/twistK.lean#L145)

有一个对重复标签有用的推论。将一列复系数按平方自由标签 $u$ 合并，设所有标签与基模互素、两两乘积处在 (2) 的范围内。写 $\operatorname{sf}(k)$ 为 $k$ 中以奇数次出现的素数之积。由于 $v=\operatorname{sf}(ue_*)$ 对 $u$ 单值且为对合（例外能出现时 $e_*$ 为平方自由数），有

$$
\mathbb E_V\left|\sum_u c_u\lambda_V(u)\right|^2
\ \le\ 2\sum_u|c_u|^2+ P^{-C}\left(\sum_u|c_u|\right)^2 .
\tag{3}
$$

证明只需展开平方。$\operatorname{sf}(uv)=1$ 给出对角项。例外项的绝对和由 Cauchy 不等式至多再给一份对角项。其余项用 (2)。所以遇到一个例外时，也可以直接保留全部行。后面要求主项常数时，还要仔细数哪些行可能相连。

## 2. 两个零点窗和只翻转一个端点的差

只求和素数 $p\equiv5\pmod6$，充分大时这包含 $T(U)$ 和 $B(U)$ 的全部候选。记 $\ell=\log x$，$a_d\equiv2d^{-1}\pmod p$，$\mathcal D(p)=\{d:x^a\le d\le x^b,\ (d,6)=1\}$。因 $b<1/2$，所有这样的 $d$ 均小于 $\sqrt x$，取 $b_x(d)=\ell^2/(2\log d\log(x/d))$。

简述前篇的投影。令 $\eta$ 为模 $3$ 的非主特征，$\omega$ 取 $1$ 或 $\chi_8=(8/\cdot)$。对每个模 $p$ 的奇特征 $\chi$，使用本原 $L(s,\chi\eta\omega)$。以 $N$ 表示窗口权的阶数，$J_0=\int_{\mathbb R}e^{-\cosh t}\,dt$。令 $V$ 为 $2N$ 个独立的 $[-4,4]$ 均匀变量之和，再加上独立的 $2j$，其中 $|j|<N$、概率为 $(N-|j|)/N^2$。对指定窗心 $h$ 写

$$
\begin{aligned}
H_N(s)&=\frac{2\pi}{32J_0}\,\mathbb E\exp[-\cos((s+iV)/32)],\\
D_h&=2+\log(24p(3+|h|)),\quad
R_h=10N-2+32\log\!\left(64\,[\log(24p)+\log D_h]\right),\\
F_\omega(c,y,h)&=-\frac2{p-1}
\sum_{\chi(-1)=-1}\overline{\chi(c)}
\sum_{\substack{L(\rho,\chi\eta\omega)=0,\ 0<\Re\rho<1\\|\Im\rho-h|\le R_h}}
y^\rho H_N(\rho-ih).
\end{aligned}
\tag{4}
$$

内层按重数求和。下式左因子一律使用 $F_1(\cdot,\cdot,h_1)$，右因子一律使用 $F_{\chi_8}(\cdot,\cdot,h_2)$，简记为 $F_1,F_8$。两窗可预先分离，也可随 $p,Q,D$ 选取。定义

$$
\begin{aligned}
A(n)&=(1+\lambda_Q(n))/2,\qquad C(n)=(1-\lambda_D(n))/2,\\
I(p)&=A(p)A(x)C(x)C(p),\qquad W(p)=A(p)A(x)C(x),\\
\phi_\sigma(d)&=(1+\sigma_1\lambda_Q(d))(1+\sigma_2\lambda_D(d))/4,\\
E_p&=\ell^{-2}\Re\!\left[\chi_8(x)F_1(2,x)\overline{F_8(2,x)}\right],\\
K_{p,\sigma}&=\ell^{-2}\Re\!\sum_{d\in\mathcal D(p)}
  b_x(d)\phi_\sigma(d)\chi_8(x)\chi_8(d)
  F_1(d,d)\overline{F_8(a_d,x/d)} .
\end{aligned}
\tag{5}
$$

在 $K_{p,\sigma}$ 里遍历整个整数区间，不用输入 $x$ 的因子表。前篇式 (6) 至 (9) 的远处代表估计先按乘积合并，再用 $\sum_{d\mid m}\Lambda(d)\Lambda(m/d)\le(\log m)^2$。因此可插入任意 $|\phi_\sigma(d)|\le1$ 的权。置 $e_p=4^{-N}+p^{-2}$，分别对两项得到

$$
E_p=\frac{\Lambda(x)^2}{\ell^2}+O(e_p),\qquad
K_{p,\sigma}=\sum_{\substack{d\in\mathcal D(p)\\ d\mid x}}
\frac{b_x(d)}{\ell^2}\phi_\sigma(d)\Lambda(d)\Lambda(x/d)+O(e_p).
\tag{6}
$$

两个误差对全部外部符号和两个窗心一致。素数 $x$ 为 $(E,K)$ 贡献算术值 $(1,0)$。对于内部半素数 $x=rs$，算术值为 $(0,\phi_\sigma(r)/2)$。导子分别是 $3p,24p$，每个二次项来自两个不同且非共轭的函数。涉及幂次精度时以下取 $N=\lceil2\log p/\log4\rceil$，于是 $e_p\ll p^{-2}$。对 $|h_j|\le\exp(p^c)$、$c>0$ 固定，所需半窗宽为 $O_c(\log p)$。[前篇的投影、逐项误差与统一截高](https://lixiang90.github.io/vibemath-NumberTheory/zeros-exact-distribution/articles/gap-two-separated-zero-windows.html)

在 $(U,2U]$ 中求和，用常权 $\mu=(\log U)^2/U$。定义

$$
\begin{aligned}
L_\sigma&=\mu\sum_p I(p)(E_p+K_{p,\sigma}),\\
R_\sigma&=\mu\sum_p W(p)(1-C(p))K_{p,\sigma},\\
\Delta_\sigma&=\mu\sum_p W(p)(2C(p)-1)K_{p,\sigma}
               =-\mu\sum_p W(p)\lambda_D(p)K_{p,\sigma}.
\end{aligned}
\tag{7}
$$

由定义，$L_\sigma-R_\sigma=\mu\sum_p I(p)E_p+\Delta_\sigma$。在 $R_\sigma$ 中只翻转了 $D$ 在 $p$ 处的符号，且只计因子项。两个和在完全相同的尺度和零点上求值。$\Delta_\sigma$ 的二阶矩是四次零点统计量。

## 3. 重复因子、例外扭转和系数 1/128

令 $\mathcal B(U)$ 是内部半素数所对应的 $p$ 集。对 $p\in\mathcal B(U)$，写 $x=rs$，$r<s$。在同一个大二倍区间内，三类素数 $p,r,s$ 的大小范围两两分离。对某个固定的 $r$，需对 $s,rs-2$ 同时筛素数。在长 $O(U/r)$ 的区间中，对每个小奇素数去掉 $0,2r^{-1}$ 两个类，维数二上界筛给出 $O(U/(r\log^2 U))$ 行。固定 $s$ 时同样至多有 $O(U/(s\log^2 U))$ 行，因为 $U/s\gg U^a$。因此一个素数标签关联的半素数行数最多为

$$
d_{\max}\ll \frac{U^{1-a}}{\log^2 U}+1 .
\tag{8}
$$

可用 Montgomery–Vaughan 的算术大筛，在变量区间长 $V$ 时只筛到 $V^{1/24}$、水平取 $V^{1/3}$。禁类数恒为二，分母 $\gg\log^2 V$。这里固定因子总大于筛限，故常数不依赖因子。[Montgomery–Vaughan，推论 1](https://personal.science.psu.edu/rcv4/personal/Publications/large_sieve.pdf#page=2) · [前篇第 5 节的逐因子估计](https://lixiang90.github.io/vibemath-NumberTheory/zeros-exact-distribution/articles/gap-two-separated-zero-windows.html)

除 $\mathcal B(U)$ 外，(6) 中有因子贡献的 $x$ 都被某个至少为 $U^a$ 的高次素数幂整除。由高次幂计数 $O(\sqrt v\log(2v))$ 和部分求和，这样的行至多 $O(U^{1-a/2}\log U)$ 个。每行的算术因子总权 $O(1)$。另有 $O(\sqrt U\log U)$ 个高次幂 $x$ 需要在端点项 $E_p$ 中处理。连同 $T(U)+B(U)\ll U/\log^2U$，全部非零算术行有 $O(\mu^{-1})$ 个。

先只对 $\mathcal B(U)$ 计算 (7) 的最后一式。由于 $\lambda_D(x)=\lambda_D(r)\lambda_D(s)$，其算术贡献等于

$$
-\frac{\mu}{2}\sum_{p\in\mathcal B(U)}
 G_p(Q)\,\lambda_D(p)\,J_p(D),
\quad
\begin{aligned}
G_p(Q)&=\mathbf1_{(\lambda_Q(p),\lambda_Q(r),\lambda_Q(s))=(1,\sigma_1,\sigma_1)},\\
J_p(D)&=\mathbf1_{(\lambda_D(r),\lambda_D(s))=(\sigma_2,-\sigma_2)} .
\end{aligned}
\tag{9}
$$

**对角项。** 固定 $Q$，在 $D$ 上展开 $J_p$。常数项是 $1/4$。其余三项使用 (2)，只有 $r,s,rs$ 中等于第二族例外指标的项可能较大。对固定 $Q$，这些例外至多涉及 $O(d_{\max})$ 行。再对 $Q$ 展开 $G_p$，常数项为 $1/8$，七个非空标签中若有一个等于第一族例外，它也至多涉及 $O(d_{\max})$ 行。三个素数均大于 $z$。故对角主项是 $\mu^2 B(U)\,(1/2)^2(1/4)(1/8)=\mu^2 B(U)/128$。

**交叉项。** $p\ne p'$ 时，固定 $Q$ 后，每个由 (9) 得到的 $D$ 标签的平方自由核都含有 $p,p'$，并且只有这两个素因子大于 $U$。故对于该 $Q$ 唯一可能的例外指标，最多只有一对无序的 $\{p,p'\}$ 会命中。所有普通交叉项共为 $O(P^{-C})$，因为每行展开系数的绝对值之和有绝对上界，归一化总质量有界。例外交叉项为 $O(\mu^2)$。若先抹去符号 $\lambda_D(p)$，仅凭 $r=r'$ 就可以有许多相关行。

**高次幂及截断。** 把有因子贡献的其余行放回。它们的对角项总计至多 $O(\mu^2 U^{1-a/2}\log U)$。对每个实际因子 $d\mid x$，$d,x/d\ll U^{1-a}<U$。纵使二者为同底素数幂，交叉项的大素因子仍恰好是 $p,p'$。基模内的小素因子取出后这一点不变。对同一 $p$ 合并全部展开的绝对权仍为 $O(1)$。端点 $E_p$ 没有进入 (9)。从算术量返回 $\Delta_\sigma$ 的逐块绝对误差为 $O(U^{-2}\log U)$，因为这里只需 $\#\{p\}\ll U/\log U$。以上各项和 $\mu^2 d_{\max}$ 全部落在 (1) 的余项内，取 $C>4J+10$ 即可。类似地展开一次和，最后阶段的例外只涉及至多一个大素数 $p$，故 $\mathbb E\Delta_\sigma=O(\mu)$。在 (1) 的精度内也可把二阶矩写成方差。

还可以单独核对 $L_\sigma$ 的均值。真孪生素数要满足四个端点符号，每对的常数项为 $1/16$。半素数还指定 $r$ 在两处的符号，各族分别固定三个素数符号，常数项为 $1/64$，因子权再乘 $1/2$。先固定 $Q$ 对 $D$ 使用 (2)，再平均 $Q$。每个例外至多涉及 (8) 的行数，端点为素数的情形至多涉及常数行。高次幂按行数估计，即得

$$
\begin{aligned}
\mathbb E L_\sigma&=\mu\left(\frac{T(U)}{16}+\frac{B(U)}{128}\right)
           +O\big(U^{-a/2}(\log U)^3\big),\\
\mathbb E L_\sigma-\mu^{-1}\mathbb E\Delta_\sigma^2
&=\frac{\mu T(U)}{16}+O\big(U^{-a/2}(\log U)^3\big).
\end{aligned}
\tag{10}
$$

这里的期望均使用先选 $Q$、再依该 $Q$ 选 $D$ 的原概率。对四种 $\sigma$，首行均有固定正下界。若陈氏下界全由半素数承担，第二行仍可以趋于零。

## 4. 原辅助素数抽样中的正态极限

**定理。** 设 $P\to\infty$，预先给定 $P\le U$, $2U\le P^J$，并沿这列尺度有 $\mu B(U)\ge c>0$。固定 $\sigma$，则 (7) 中用零点定义的实数满足

$$
\frac{\Delta_\sigma}{\mu\sqrt{B(U)/128}}\ \Longrightarrow\ \mathcal N(0,1).
\tag{11}
$$

随机性来自抽取辅助素数。零点属于由物理端点 $p$ 决定的两族 Dirichlet 函数，窗心可以预先指定；(6) 的一致性也容许窗心依赖抽到的辅助对。证明使用固定阶矩，每一阶都须保留原来的两阶段次序。

**独立符号下的极限。** 暂给每个素数两枚独立均匀的符号，分别供第一、第二阶段使用，且不同素数彼此独立。简记 $B=B(U)$。给定所有第一阶段符号及 $r,s$ 的第二阶段符号，令 $M=\sum_{p\in\mathcal B(U)}G_pJ_p$。余下的 $M$ 个端点 $p$ 各带一枚独立均匀的第二阶段符号。每行入选的概率是 $1/32$，两行只有共享因子才可能在入选事件上相关。由 (8)，

$$
\mathbb E M=B/32,\qquad \operatorname{Var} M\ll B(1+d_{\max}),\qquad
d_{\max}/B\ll_c U^{-a}.
\tag{12}
$$

故 $32M/B\to1$ 依概率成立。条件于上述符号，把 (9) 除以 (11) 的分母，其特征函数恰为 $\cos(t\sqrt{32/B})^M$，趋向 $\exp(-t^2/2)$。也可直接得到全部固定阶矩。$M$ 个独立均匀符号之和的奇数阶矩为零，第 $2j$ 阶矩为 $(2j-1)!!M^j+O_j((1+M)^{j-1})$，因首项来自 $j$ 对互异的下标。用 $0\le M/B\le1$ 和 (12) 即得标准正态的矩。

**保留一个例外时的高阶矩界。** 固定 $Q$，取任意由 $Q$ 确定的 $m$ 行子集 $\mathcal A\subset\mathcal B(U)$。设 $S_\mathcal A=\sum_{\mathcal A}G_p(Q)\lambda_D(p)J_p(D)$。对每个固定整数 $j\ge1$ 和任意固定 $C>0$，有

$$
\mathbb E_{D\mid Q}|S_\mathcal A|^{2j}
\ll_j m^j+P^{-C}m^{2j}.
\tag{13}
$$

将 $J_p$ 展成四个符号单项式，其系数绝对值之和为 $1$。取 $2j$ 行的任一项，所得 $D$ 标签在 $(U,\infty)$ 内的素因子，恰为出现奇数次的那些 $p$。平凡标签要求这些 $p$ 的集合为空，有 $O_j(m^j)$ 个有序组。例外标签若能命中，这个集合也已由该例外指定。如果含 $2v$ 个素数，先给每个指定的素数放入一次，再把剩下的 $2j-2v$ 个位置配对，可能性至多 $O_j(m^{j-v})$。若例外的所有素因子均小于 $U$，仍用空集的界即可。余下标签用 (2)，乘积至多为 $P^{O(jJ)}$。在独立的第二阶段同样有 (13)，且可省略第二项。

**从真实抽样转到独立符号。** 对目标阶数 $k$，取 $j\ge k$，为 (2) 选择足够大的固定 $A,C$，例如控制 $2j$ 阶时取 $A,C>(4j+4)(J+1)$。先固定真实 $Q$。如果第二阶段的非平凡例外可能出现在至多 $k$ 行的乘积中，选它的一个素因子，暂时删去含该素数的全部行。至多删 $d_{\max}$ 行，选择只依赖 $Q$。余下行的乘积没有这个例外，(2) 给出实际 $D$ 与独立符号 $D$ 的第 $k$ 阶矩之差至多 $O(P^{-C}B^k)$。若例外不能出现，直接比较全部行。

在两种分布下，将删去的行加回所引起的误差都可控制。由 (13)，全部行的 $L^k$ 范数为 $O_k(\sqrt B)$，删去部分的 $L^k$ 范数为 $O_k(\sqrt{d_{\max}})$，常数中已吸收 $P^{-C}$ 项。在 $|x^k-y^k|\le k|x-y|(|x|+|y|)^{k-1}$ 中用 Hölder 不等式，除以 $B^{k/2}$ 后，误差为 $O_k(\sqrt{d_{\max}/B})=o(1)$。这一步对每个 $Q$ 一致。

此时 $D$ 已独立。再对第一阶段的唯一可能例外删去含其某个素因子的行。固定独立 $D$，将第一阶段各 $G_p$ 的八项展开，(2) 以同样的精度给出独立 $Q$ 的矩。在真实和独立的 $Q$ 下，加回删去行时都先对独立 $D$ 用 (13)。由此，(9) 标准化后的每个固定阶矩均趋于上述正态矩。正态分布由其矩唯一确定。

**返回零点。** 有高次幂因子的其余行，按第 3 节的唯一大端点论证，其换号和的二阶矩为 $O(\mu^2 U^{1-a/2}\log U+\mu^2+P^{-C})$。除以 $\mu^2 B$ 后趋零。零点截断的确定性误差 $O(U^{-2}\log U)$ 也比 (11) 的标准差小，由此得到 (11)。

原证明的解析好对事件 $\mathcal G_P$ 满足 $\mathbb P_P(\mathcal G_P)=1-o(1)$。事实上固定所需的 $\epsilon,\tau,J,N_0$ 后，原定理中的常数在任意固定舍弃比例 $\eta>0$ 之前选定，故通过比例的下极限至少为 $(1-\eta)^2$，再令 $\eta\downarrow0$。因此 (11) 也适用于条件分布 $\mathbb P_P(\cdot\mid\mathcal G_P)$。对另加某个指定上同调窗口条件的子集，其概率可只有固定正比例，须使用下一节给出的量词。[原两阶段的保留比例](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/MeanSquareBounds/PairShare.lean#L22) · [转成原解析窗口条件](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/MeanSquareBounds/pairEx.lean#L146)

**四个符号类的联合分布。** (11) 中的四个坐标联合趋于四个独立的实标准正态。事实上，在独立符号模型中给定前述符号后，每个入选的 $p$ 只能分给 $\sigma=(\lambda_Q(r),\lambda_D(r))$ 这一格。各格要相加的第二阶段端点符号集合互不相交。各自的 $32M_\sigma/B$ 由 (12) 同时依概率趋于 $1$。条件特征函数为四个余弦幂的乘积，极限为 $\exp(-\sum_\sigma t_\sigma^2/2)$。对任意固定的实线性组合 $\sum_\sigma t_\sigma\Delta_\sigma$，每行的符号展开绝对权至多乘上 $\sum_\sigma|t_\sigma|$，其中的第二阶段大素数仍恰为 $p$。故 (13)、删一素因子以及矩的转移逐项适用。由 Cramér–Wold 判别法即得联合结论。

## 5. 一个固定概率门槛

取同一零点族可计算的 $m_\sigma=\mathbb E L_\sigma$。由陈氏下界和 (10)，充分大的尺度上 $m_\sigma\asymp1$。引入

$$
V_\sigma=\frac{\Delta_\sigma}{\sqrt{\mu m_\sigma}},\qquad
v_U=\frac{B(U)}{8T(U)+B(U)}.
\tag{14}
$$

$V_\sigma$ 的定义只需要 (4)、(5)、(7) 的零点和、已知的素数导子及辅助抽样。引入右端的算术比值是为了说明它测到的量。记 $\Phi$ 为实标准正态分布函数。对任意固定 $t>0$，有

$$
\mathbb P_P\bigl(|V_\sigma|\le t\bigr)
=2\Phi(t/\sqrt{v_U})-1+o(1),
\tag{15}
$$

当 $v_U=0$ 时右端主项取 $1$。证明可对任意子列抽取 $v_U$ 收敛的子列。其极限若为正，因 $\mu(8T+B)\asymp1$，(11) 可用，且换算到 (14) 的方差趋于这个极限。若极限为零，(1) 和 (10) 直接给 $V_\sigma\to0$ 依概率成立。两种情形的 $[-t,t]$ 边界皆无极限原子。这个论证也给出平滑版 $\mathbb E\cos(t V_\sigma)=\exp(-t^2 v_U/2)+o(1)$。

特别地，沿任意预先指定的尺度序列，

$$
\limsup \mu T(U)>0
\quad\Longleftrightarrow\quad
\limsup \mathbb P_P\bigl(|V_\sigma|\le1\bigr)>2\Phi(1)-1
=0.68268949\ldots .
\tag{16}
$$

这是因为 $1-v_U=8T/(8T+B)$，且 $\mu(8T+B)$ 有绝对上下界。右侧要求超过常数一个固定幅度。半素数和真孪生素数可以各占正比例，(16) 都适用。也可取单一频率，用 $\limsup\mathbb E\cos V_\sigma>e^{-1/2}$ 作等价判据。将概率和期望限制到 $\mathcal G_P$ 的条件分布仍成立，归一化的 $m_\sigma$ 沿用完整抽样。

若四格同时检验，(15) 的主项为 $[2\Phi(t/\sqrt{v_U})-1]^4$，左端事件取 $\max_\sigma|V_\sigma|\le t$。因而在 (16) 中，四格同时落入 $[-1,1]$ 的门槛降为 $(2\Phi(1)-1)^4=0.21721653\ldots$。这是同一批辅助对上四个联合坐标的概率。

还有一个仅需对部分实际辅助对检验的充分条件。若对某列尺度存在 $\tau_k\downarrow0$、固定 $\eta>0$ 和事件 $G_k$，使

$$
\mathbb P_{P_k}\bigl(G_k\cap\{|\Delta_\sigma|\le\tau_k\sqrt{\mu_k}\}\bigr)\ge\eta ,
\tag{17}
$$

则 $\mu_k B(U_k)\to0$ 且 $T(U_k)\gg U_k/\log^2 U_k$。否则取 $\mu_k B\ge c>0$ 的子列，(11) 的连续极限使其中小值的总概率趋零。若 $\mathbb P(G_k)\ge\eta$ 且 $\mathbb E(\Delta_\sigma^2\mid G_k)=o(\mu_k)$，Markov 不等式也给出 (17)，把右边常数缩小一半即可。这里容许 $G_k$ 来自原证明在两阶段之后选出的共同窗口。仅持有某个正比例事件上的条件概率时，不能把 (16) 的 $0.68268949\ldots$ 原样当作该事件内部的门槛。

**在共同窗口上只需固定幅度的均方改进。** 对质量至少为 $\eta\in(0,1)$ 的任何事件 $G_k$，设 $S=\sum_\sigma V_\sigma^2$。令 $F_d$ 为 $d$ 自由度的卡方分布函数，定义

$$
b_4(\eta)=\eta^{-1}\int_0^\eta F_4^{-1}(u)\,du
         =\frac{4F_6(F_4^{-1}(\eta))}{\eta}.
\tag{17a}
$$

沿 $v_{U_k}\to v$ 的每个子列，

$$
\liminf_{k\to\infty}\mathbb E(S\mid G_k)\ \ge\ v\,b_4(\eta).
\tag{17b}
$$

由联合极限，$v>0$ 时 $S$ 依分布趋于 $v\chi^2_4$。写 $F_S$ 为实际分布函数；对任意 $y\ge0$，$\mathbb P(S>y\mid G_k)\ge(1-F_S(y)/\eta)_+$。积分并用 Fatou 引理得到 (17b)。$v=0$ 时非负性足够。对卡方密度有 $xF'_4(x)=4F'_6(x)$，于是 (17a) 的两个表达相等。在连续极限中选取最低的 $\eta$ 分位恰好取到这个常数，因而只知道事件比例时门槛不能再提高。

例如保留四分之一辅助对，$b_4(1/4)=1.17309886\ldots$。若在质量至少 $\eta$ 的共同窗口事件上证明 $\mathbb E(S\mid G_k)\le b_4(\eta)-\delta$，某固定 $\delta>0$ 对无穷多个预先指定的二倍区间成立，由 (17b) 得这些区间 $\limsup v_{U_k}\le1-\delta/b_4(\eta)$，故 $\liminf\mu_kT(U_k)>0$。这里仅需一个固定的均方缺口，也容许半素数有正密度。$\eta$ 取原选窗实际能保留的比例时可直接使用同一公式。$\mathcal G_P$ 的通过比例及例外扭转指标的个数本身均未给出这项均方上界。

为何原均方引理尚差一步，可以从它的量词看清。原式对一个所有辅助素数共用的系数列 $|\kappa(n)|\le1$，在一族辅助素数上给 $Cy\log y$ 的均方上界。(7) 还同时使用两端 $p,p+2$、因子 $d$ 的符号及两族零点；第二阶段的系数依赖第一阶段的选择。即使用 (2) 逐个比较所有固定阶矩，算出的方差中仍留下 $B/(8T+B)$。为使 (16) 或 (17b) 严格超过现有水平，须另外估计保留这些端点和因子条件后的有符号和，或其特征函数。原引理没有给出可直接与 (17a) 比较的首项常数。[原均方引理的定义](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/MeanSquareBounds/PrimeSetTwoBounds.lean#L24)

## 6. 同一个 ζ 好窗口上的四个符号类

用原来的物理窗口 $I_Y=(2Y,Y^{3/2}/2]$，其中 $P\le Y,\ Y^2\le H$。令 $L_Y=\sum_{Y\le n\le Y^{3/2}/2}1/n\asymp\log Y$，$w(p)=\log p\log(p+2)/p$。将 (7) 中的 $\mu\sum_p$ 改为 $L_Y^{-1}\sum_{p\in I_Y}w(p)$，记作 $L_\sigma(Y),R_\sigma(Y),\Delta_\sigma(Y)$。用相同的加权法对真实孪生、内部半素数分别计数，记作 $M_T(Y),M_B(Y)$。在真孪生素数项再乘上四个端点符号的指示函数 $I(p)$，所得为 $J_{Q,D}(Y)$。陈氏下界逐个用到窗内的完整二倍区间，上界筛用于覆盖全窗，得到 $M_T(Y)+M_B(Y)\asymp1$。

**窗口估计。** 除去原抽样中概率 $O_J(1/\log P)$ 的辅助对，以下三式对全部上述 $Y$ 和四个 $\sigma$ 同时成立，因而允许在给定 $Q,D$ 之后再选 $Y$。

$$
\begin{aligned}
L_\sigma(Y)&=\tfrac1{16}M_T(Y)+\tfrac1{128}M_B(Y)
             +O_J\big(Y^{-a/2}(\log Y)^2\big),\\
R_\sigma(Y)&=\tfrac1{128}M_B(Y)+O_J\big(Y^{-a/2}(\log Y)^2\big),\\
L_\sigma(Y)-R_\sigma(Y)
 &=J_{Q,D}(Y)+O_J\big(Y^{-1/2}(\log Y)^2\big)\\
 &=\tfrac1{16}M_T(Y)+O_J\big(Y^{-1/2}(\log Y)^2\big).
\end{aligned}
\tag{18}
$$

**一致性的证明。** 取 $V=2^jP$，在 $(V,2V]$ 考察任意子区间的未归一化加权和，最后一段可在 $H$ 处截断。对真实半素数把 $(p,r,s)$ 在两族中的六个符号全部展开。常数项在每个指定符号格的比例为 $1/64$。任一其余单项式，若第二阶段非空，就先固定 $Q$，把其第一阶段符号当作模长至多 $1$ 的系数；若只有第一阶段非空，就用第一阶段。对当前非空的标签子集 $S\subset\{p,r,s\}$，按标签乘积合并系数，用 (3)。

当 $S=\{r\}$ 时，(8) 的逐因子估计给每个纤维的总绝对权至多 $C/r$。其平方之和至多 $C\sum_{r\ge V^a,\ r\ \mathrm{prime}}r^{-2}\ll V^{-a}/\log V$。当 $S=\{s\}$ 时同理，使用 $s\gg V^{1-b}$。其余非空子集均唯一确定一行，因为 $p$ 唯一而 $rs=p+2$，三类素数的大小范围互不相交。这时平方权之和为 $O((\log V)^2/V)$。真孪生素数的四个符号同样展开，常数比例为 $1/16$，每个非空子集都唯一确定一行；两个端点分属模 $6$ 的不同类。

为估计 $\Delta_\sigma$，先固定 $Q$。对 (6) 中所有非零因子行展开第二阶段。每项都有一个且只有一个大于 $V$ 的素数标签 $p$，即使小因子是素数幂或有部分符号由基模固定也是如此。每行全部展开项的绝对权不超过 $Cw(p)$，合并标签后的平方权之和因此至多 $C\sum w(p)^2\ll (\log V)^2/V$，和遍及 $O(V/\log^2V)$ 个非零端点。普通的非空字符项与例外标签都已计入 (3)。

为取子区间的上确界，在整数端点上作二进剖分。每个前缀由每层至多一块组成。先用 Cauchy 不等式支付一个层数，再对各层所有块取二阶矩的和。由于使用的是绝对权，同层各子块的平方纤维和合计至多为整段的平方纤维和，(3) 中绝对权总和的平方也有相同性状。两个前缀之差给出任意子区间。总计损失 $O((\log V)^2)$ 因子。故半素数格的中心化加权和，其子区间上确界的二阶矩为 $O(V^{-a}\log V)$。真孪生素数格的中心化加权和以及全部因子行的换号和，对应的界均为 $O(V^{-1}(\log V)^4)$。用 (2) 时取 $C>4J+10$ 即可。

分别以 $V^{-a/2}(\log P)^3$ 和 $V^{-1/2}(\log P)^3$ 作阈值。Chebyshev 不等式对 $O_J(\log P)$ 个块及有限种符号一起应用，总失败概率为 $O_J(1/\log P)$。原始 $L_\sigma,R_\sigma$ 中高次幂因子的确定性贡献每块至多 $O(V^{-a/2}(\log V)^3)$，端点高次幂贡献至多 $O(V^{-1/2}(\log V)^3)$。零点截断每块误差为 $O(V^{-2}\log V)$。从 $V\asymp Y$ 起按几何级数求和，再除以 $L_Y$，即得前两行。第三行直接在 $L_\sigma-R_\sigma=L_Y^{-1}\sum w I E+\Delta_\sigma(Y)$ 使用较强的换号界，保住了平方根量级。

原几何选窗保留有限个候选起点 $Y_i=P^{(2C_t/\eta)^i}$，$0\le i\le R_0$，取外端 $J$ 足够大。每个解析好对至多有 $R_0$ 个例外扭转指标，而各段的可能冲突区间 $(Y_i^\eta,Y_i^{2C_t}]$ 两两不交。因此每个解析好对至少有一个完整 $[Y_i,Y_i^2]$ 合格，某个固定 $i=i(P)$ 至少被 $1/[4(R_0+1)]$ 比例的辅助对共享。与 (18) 的事件相交后仍保留正比例。由此在同一原 ζ 好窗口内，四个 $L_\sigma(Y)$ 同时具有固定正下界。[原几何选窗的有限例外](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/AutomorphicVarieties/LayeredIsomorphisms/ThmBabOf.lean#L76) · [同一辅助对保留全窗口](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/CharacterBoxes/WinPairs.lean#L32) · [原共同好窗口的质量与归一化](https://github.com/lixiang90/vibemath-NumberTheory/blob/main/Goldbach_TwinPrimes/twin-prime-auxiliary-signs.md)

在这批窗口上，前文的 $\Gamma_1$ 等于 $\frac13\limsup (L_\sigma(Y)-R_\sigma(Y))=\frac1{48}\limsup M_T(Y)$，上极限在所选序列上取。单个二倍区间满足 (16) 的严格不等式只提供一段的质量。若沿一列共同好窗口，每个窗口有 $\gg\log Y$ 个不相交的完整二倍区间，在其上 (15) 的概率以同一个固定正量超过 $2\Phi(1)-1$（取 $t=1$），则每段有 $\gg1$ 的未归一化真孪生权重，遂得 $\Gamma_1>0$。现有输入证明了 (18) 与候选的正下界，尚未给出这个严格超额。

若只研究 $\Gamma_1$ 的极限值，在 $|h_j|\le\exp(p^c)$ 的同样范围，可取 $N=\lceil(A+2)\log\log p/\log4\rceil$，窗宽降到 $O_{A,c}(\log\log p)$。此时 (18) 另加 $O_A((\log Y)^{-A})$ 的误差。第 4、5 节检验单个二倍区间的微小波动时使用第 2 节的较高精度。

## 7. 数值复核

取 $53+2=5\cdot11$ 和 $83+2=5\cdot17$，共享因子 $5$。给两阶段的五个素数标签逐一赋予符号，共 $2^{10}=1024$ 种。再加入 $89+2=7\cdot13$，共 $2^{16}=65536$ 种，已有两格同时非空的情形。逐项核实四格均值、对角二阶矩系数 $1/128$、六个非对角二阶矩为零，以及给定筛选结果后的整符号和四阶矩 $3M^2-2M$。程序还检查了 (3) 的例外标签对合。

另取 $48$ 个实际的 $Q\equiv3\bmod8$，每个配 $24$ 个 $D\equiv7\bmod16$、$D\equiv-1\bmod 105Q$ 的实际素数，依原次序平均。因子区间在这项有限实验中取 $[x^{3/10},x^{9/20}]$。样本用于核对符号、归一化和高次幂，所取参数没有验证原构造的渐近尺寸或上同调条件。表中为便于与整数计数比较，把每行的 $\mu$ 先除掉。$M_1$ 为四格 $L_\sigma/\mu$ 的样本均值，$M_2$ 为四格 $(\Delta_\sigma/\mu)^2$ 的样本均值。末两列用固定的 $(-,+)$ 一格，算术行保留全部高次幂，正态比较值取 (15) 的 $t=1/4$。

| $U$ | 真孪生 $T$ | 内部半素数 $B$ | $128M_2$ | $16(M_1-M_2)$ | 小值比例 | 正态比较值 |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1,000 | 26 | 10 | 8.658 | 23.762 | 0.796 | 0.757 |
| 10,000 | 137 | 86 | 95.992 | 136.299 | 0.651 | 0.646 |
| 100,000 | 936 | 547 | 563.681 | 937.898 | 0.638 | 0.662 |

对最后一行，仅保留真实内部半素数并按 (11) 标准化，固定 $(-,+)$ 一格的二、四阶矩分别为 $1.008,3.387$，样本在零点的原子质量为 $0.100$。第一行相应的原子质量高达 $0.796$，说明小样本下的离散性很明显。表中没有用于推断无穷多尺度上的严格超额。

最后取 $Q=1000003$，$D=47880143639$，在高度 $320,960$ 上分别用导子 $3p,24p$ 的围道积分计算 (5)。$p=23$ 时 $x=25$ 的端点量为 $1/4$、因子量为零；$p=29$ 时 $x=31$ 的端点量为 $1$；$p=53,83$ 时端点量均为零，同一 $(-,+)$ 因子格均为 $1/2$，而 $D$ 在 $53,83$ 处的符号相反。四个例子与直接算术值的最大绝对偏差依次为 $4.26\cdot10^{-12},1.33\cdot10^{-11},6.61\cdot10^{-16},4.29\cdot10^{-15}$。计算使用两层积分网格，并在各窗心用 Arb 独立核对对数导数。这是围道数值复核，未逐个分离非平凡零点。[复算程序及运行参数](https://github.com/lixiang90/vibemath-NumberTheory/blob/main/zeros-exact-distribution/gap-two-chen-auxiliary-switch/verify.py)
