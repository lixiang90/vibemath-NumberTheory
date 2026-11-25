---
title: "分端点的三段零点与陈氏候选的显式下界"
date: "2026-11-25"
abstract: "在固定差二的两个端点分别用二次特征分组，从三个共用的收缩零点窗口读出两个方向各自的半素数质量。对每个方向，粗筛指数为 3/11 时，真孪生与半素数邻项的归一化质量之和下极限大于 0.026434。将不同因子比例的下界叠加，得到单方向的三阶零点读数与 0.4138875 倍孪生质量之和的显式下界 0.0107210。文中给出已证不等式容许孪生质量为零的例子。"
---

# 分端点的三段零点与陈氏候选的显式下界

最近三篇中，最有希望改变孪生素数判据的是按半素数的因子大小做切换筛。第 35 篇选出了需要控制的因子比例，第 36 篇开始估计它的密度，第 37 篇用可分解线性筛给出了更小的上界，也给出了把零点读数反演到整数精度的方法。再提高反演精度，只会更准确地测出窗口里已有多少对。本篇转而给这些窗口里的候选数建立可复算的下界，并把需要压低的零点读数拆成两个方向。这样只要其中一个方向出现严格缺口，就能推出真孪生。[因子比例](https://lixiang90.github.io/vibemath-NumberTheory/zeros-exact-distribution/articles/gap-two-chen-factor-ratios.html) · [切换筛](https://lixiang90.github.io/vibemath-NumberTheory/zeros-exact-distribution/articles/gap-two-switched-cubic-bound.html) · [可分解筛与零点反演](https://lixiang90.github.io/vibemath-NumberTheory/zeros-exact-distribution/articles/gap-two-linear-sieve-cubic-inversion.html)

较早讨论的三种显式公式思路也说明了这个选择。黎曼零点的精确重构能找出每个已给窗口的孪生对，零点差的等价判据能识别任意稀疏的孪生序列。固定平移的 Dirichlet 相关级数直接记录差二，但它在 $\Re s=1$ 上所需的边界行为本身已含很强的素数对信息。下文采用第三种思路中的固定平移，转到可用算术分布定理估计的一族二次 Dirichlet 特征，再问这些估计能强迫哪些零点读数出现。[Korevaar 对固定素数对的 Tauber 型判据](https://doi.org/10.1016/S0019-3577(05)80013-8)

## 1. 同一窗口上的两个方向

固定原构造允许的 $\epsilon,J$，取 $H=P^J$，$[Y,Y^2]\subset[P,H]$，令 $P\to\infty$。原 $\zeta$ 构造先取辅助素数 $Q$，再在同余候选中取 $d_0$，最后选几何好窗口。具体地，$z=P^{\epsilon/5},\ z'=z(\log P)^2$，$e^{z/2}<Q\le e^z,\ Q\equiv3\bmod8$，$e^{z'/2}<d_0\le e^{z'},\ d_0\equiv7\bmod16$，$d_0$ 模 $Q$ 及不超过 $z$ 的各个奇素数皆余 $-1$。令 $G_n$ 表示 $n,n+2$ 在 $Q$ 处均为二次剩余、在 $d_0$ 处均为二次非剩余。以下算术不等式在去掉 $G_n$ 时无条件成立。转移到原候选与零点族时，采用归档的 Dirichlet-GRH 及第 33 至 36 篇的扭转估计。[原构造中的候选符号](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/MeanSquareBounds/PrimeSetTwoBounds.lean#L41) · [先选辅助对再选窗](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/CharacterBoxes/WinPairs.lean#L5)

写 $a=3/11,\ b=1/2$。在 $\mathcal I_Y=(2Y,Y^{3/2}/2]$ 取 $n(n+2)$ 平方自由、所有素因子都大于 $n^a$ 的行。行质量为

$$
 \mu_Y(n)=\frac{G_n\log n\log(n+2)}{nL_Y},\qquad
 L_Y=\sum_{Y\le k\le Y^{3/2}/2}\frac1k\sim\frac{\log Y}{2}.
$$

$J_Y$ 为两端均为素数的质量。正方向把 $(v,m)=(n,n+2)$，负方向把 $(v,m)=(n+2,n)$。当 $v$ 为素数、$m=st$ 且 $s<t$ 是素数时，把该行质量计入 $B_{\pm,Y}$ 在 $u=\max(a,\log s/\log m)$ 处的点质量。故 $J_Y+B_{+,Y}([r,b])$ 计数 $p$ 与 $p+2$，$J_Y+B_{-,Y}([r,b])$ 计数 $p$ 与 $p-2$，两个式子都保留真孪生。令

$$
 f(u)=2u(1-u),\qquad C_{\pm,Y}=\int_a^b f(u)\,dB_{\pm,Y}(u).
 \tag{1}
$$

**主要结论。** 以下界对两个方向分别成立。在原两阶段均匀抽样概率为 $1-o(1)$ 的见证集上，对允许的物理窗口一致成立，因而适用于其中由原几何条件选出的共同好窗口。

| 较小因子的比例下限 $r$ | 每个方向 $\liminf\bigl(J_Y+B_{\pm,Y}([r,b])\bigr)$ 的严格下界 |
| --- | ---: |
| $3/11$ | $0.026434$ |
| $11/40$ | $0.023652$ |
| $7/25$ | $0.017144$ |
| $2/7$ | $0.009396$ |
| $36/125$ | $0.006323$ |
| $29/100$ | $0.003645$ |
| $117/400$ | $0.000314$ |

第 2 节用三个共用零点窗口构造 $\widehat J,\widehat C_+,\widehat C_-$，分别为 $J_Y,C_{+,Y},C_{-,Y}$ 加 $o(1)$。在同一见证集上有

$$
 \liminf\bigl(\widehat C_{\pm}+0.4138875\,\widehat J\bigr)>0.0107210 .
 \tag{2}
$$

例如，若在无穷多个共同好窗口上，某一方向有 $\widehat C_{\pm}\le0.0107210-\delta$，其中 $\delta>0$ 固定，则在这些窗口上 $J_Y\ge\delta/0.4138875-o(1)$。即使每窗取两个方向中较小的读数，这个判据也成立。第 35 篇与原 $\Gamma_1$ 的换算给 $\Gamma_1\ge\delta/(3\cdot0.4138875)$。实际是否存在这样的缺口仍是待证的不等式。[四个辅助符号与原归一化](https://github.com/lixiang90/vibemath-NumberTheory/blob/main/Goldbach_TwinPrimes/twin-prime-auxiliary-signs.md)

## 2. 在两个端点分别给导子分组

固定物理见证后，让 $d$ 遍历同一 $Q$ 下的第二阶段全部候选，每个导子只观测 $L(s,(\cdot/d))$。沿用第 33 篇的显式公式。具体取 $B=z/8,\ k=\lceil(\log P)^2\rceil,\ b_0=B/(4k),\ R=2/b_0$。以 $2k$ 个质量为一的 $[-b_0,b_0]$ 方盒和 $\mathbf1_{[-B,B]}$ 卷积得 $g$，Fourier 相位取 $e^{-itu}$。令 $T_d(h)$ 为

$$
 -\frac12\left(\sum_{|\gamma_d-h|\le R}\widehat g(\gamma_d-h)-\log d\right)
$$

减去其导子族均值，$\tfrac12+i\gamma_d$ 的零点按重数计。三个中心取 $h_0=0,\ h_1=4R,\ h_2=8R$，窗口两两分离，且 $R\log H\to0$。将所有观测导子按单个平方自由整数 $m$ 的符号 $(m/d)=+1,-1$ 分组。令 $\Delta_1(T;m)$ 为两个条件均值正组减负组的一半，$\Delta_2(T,T';m)$ 为相应两个条件协方差之差的一半。在两组内分别中心化后计算协方差。[条件矩与独立符号比较](https://lixiang90.github.io/vibemath-NumberTheory/zeros-exact-distribution/articles/gap-two-quadratic-conductor-covariance.html)

写 $T_i=T_d(h_i)$，定义

$$
 P_i(m)=\frac{\sqrt m}{\log m}\Delta_1(T_i;m),\qquad
 V(m)=\frac{\sqrt m}{(\log m)^2}\Delta_2(T_1,T_2;m).
 \tag{3}
$$

对所有保留行的两个端点，置 $e_Y=(R\log H)^2+Y^{-a}$。任意固定 $A>0$ 下，一致有

$$
 \begin{aligned}
 P_i(m)&=\mathbf1_{m\ \mathrm{prime}}\bigl(1+O(e_Y)\bigr)+O_A(P^{-A}),\\
 V(m)&=\begin{cases}
 f(\log s/\log m)+O(e_Y)+O_A(P^{-A}),&m=st,\ s<t\ \mathrm{primes},\\
 O_A(P^{-A}),&\text{其余保留端点}.
 \end{cases}
 \end{aligned}
 \tag{4}
$$

证明只用前篇显式公式的一、二阶项。素数标签 $\ell$ 在 $T_i$ 中的系数是 $(\log\ell)/\sqrt\ell$ 乘 $1+O((R\log H)^2+1/\ell)$。一次幂的余弦给第一项误差，更高奇数次幂给第二项。在独立符号模型中，单素数分组固定了这个标签，其他符号的条件协方差在两组相同。双素数分组时，一阶均值之差为零，协方差之差的一半留下 $s,t$ 的两个排列，恰给 $2\log s\log t/(\log m)^2$。三个及以上不同素因子的符号乘积与二次以下的每个单项式正交。第 33 篇的矩比较可在乘 $\sqrt m$ 之前预取任意大的负幂精度，故得到 (4)。把 $\log s/\log m$ 在 $a$ 处截断只多出 $O(1/(Y\log Y))$。

以 $\sum_n$ 表示对保留行求和。三个零点统计为

$$
 \widehat J=\sum_n\mu_Y(n)P_0(n)P_1(n+2),\qquad
 \widehat C_{\pm}=\sum_n\mu_Y(n)P_0(v)V(m).
 \tag{5}
$$

这里先在每个端点自己的条件族内取平均，再把两个结果相乘。两端均为素数时 $V(m)$ 消失。只有 $v$ 是素数、$m$ 是半素数的指定方向留下 (1) 的权。把 (4) 代入 (5) 即得三个 $o(1)$ 近似，粗界 $\sum_n\mu_Y(n)\ll(\log Y)^2$ 已够吸收误差。前篇用整行乘积的符号分组，其匹配三阶均值满足 $\overline C=\widehat C_++\widehat C_-+o(1)$。[匹配均值的定义](https://lixiang90.github.io/vibemath-NumberTheory/zeros-exact-distribution/articles/gap-two-coherent-zero-heights.html)

这三个窗口也能读出表中的尾部质量。取 $P_0(v)>1/2$，再取 $V(m)>1/4$，对所有充分大的 $P$，这逐行精确选出指定方向的素数与半素数，因为 $f(a)>0.39$。另取固定 $a<r<1/2$ 时，把第二个截点改为 $f(r)$。由 (4) 与第 37 篇的区间密度上界，可能分错的半素数落在 $r$ 附近宽度趋零的区间内，归一化质量是 $o(1)$。因此按行截断求和可以只用三个高度。第 37 篇按线性观测先求和再统一取整，适用于另一种计算顺序。

## 3. 粗筛序列与带一个指定素数的模数

暂去物理符号 $G_n$。在 $X<v\le(1+h)X,\ h=1/\log X$ 上固定一个方向，令 $v=p$ 为素数，$m=p\pm2$。以下计数以 $4C_2hX/\log^2 X$ 为单位，其中 $C_2=\prod_{\ell>2}(1-(\ell-1)^{-2})$。设 $d_*=2/7$，$a\le c<d_*$。若只要求 $m$ 的最小素因子至少为 $X^c$，得到下界系数

$$
 H(c)=\frac74\log 3-2\log 2-\log\frac{3c}{2-5c}.
 \tag{6}
$$

在 $z_s=X^{1/7}$ 处对素数序列作下界线性筛。奇素数的局部密度为 $1/(\ell-1)$，筛积等价于 $2C_2e^{-\gamma}/\log z_s$。Bombieri–Friedlander–Iwaniec 对固定剩余类的可分解权给严格小于 $4/7$ 的水平。下界筛函数在 $2\le s\le4$ 为 $2e^\gamma\log(s-1)/s$，令水平上趋 $4/7$ 即得初项 $7\log 3/4$。素数计数可在小格两端作差，所需的误差节省可预先取为任意对数幂。[BFI，定理 10](https://doi.org/10.1007/BF02399204) · [Iwaniec 的上下界线性筛](https://doi.org/10.4064/aa-37-1-307-320) · [可分解上界筛，定理 2.10](https://arxiv.org/pdf/2109.02851#page=8)

Buchstab 恒等式从这里扣掉最小素因子 $q$ 在 $[X^{1/7},X^c)$ 的行。写 $q\asymp X^u$，在 $q\mid m$ 后上筛其余小素数，可以使用总模数水平的严格上端

$$
 \theta(u)=\min\bigl((1+u)/2,\ 1-3u/2\bigr).
 \tag{7}
$$

下面直接给出加上 $q$ 后的模数分配。Pascadi 命题 4.4 的短边为 $X^\eta$。把残余的可分解权分成两段，三段模数的指数记为 $t_1,t_2,t_3$。忽略严格余量时需 $t_1\le\eta,\ 2\eta+t_2+2t_3\le1,\ 2\eta+5t_2+2t_3\le2$。取略小于 $\theta(u)$ 的 $\theta$ 且 $\theta>2u$，按下面两种情况安排。

| 短边 | $(t_1,t_2,t_3)$ | 指定素数去向 |
| --- | --- | --- |
| $\eta\le\theta-u$ | $(\eta,u,\theta-u-\eta)$ | $q$ 单独放在第二段 |
| $\eta>\theta-u$ | $(\eta,\theta-\eta,0)$ | $q$ 和指数为 $\eta-u$ 的残余段合并 |

第一行的两个左边为 $2\theta-u$ 和 $2\theta+3u$。第二行有 $\eta\ge u$，两个左边是 $\theta+\eta$ 和 $5\theta-3\eta$。前者小于一，因为 $\theta\le5/8,\ \eta\le1/3+O(\varepsilon_0)$，后者小于 $2\theta+3u$。严格处理时将第一段实际目标从 $\eta$ 降为 $\eta-\varepsilon_0$，也用这个新目标选行，两个不等式仅损失 $O(\varepsilon_0)$。先留出水平余量，再按命题 4.4 取充分小的 $\varepsilon_0$。在 $u\le c<2/7$ 的紧区间上 $\theta(u)-2u$ 有固定正余量。[Pascadi，命题 4.4、式 (4.15) 与第 5 节](https://arxiv.org/html/2505.00653v2#S4) · [固定余类素数分解，Maynard](https://arxiv.org/pdf/2006.07088#page=7)

为说明这张表能作用在实际余项上，把 $u$ 分成固定细格 $[u_0,u_1]$，在格内用 $X^{u_0}$ 上筛、固定残余权水平为 $X^{\theta-u_1}$，总水平 $\theta$ 严格取在格内允许范围内。每个可分解权按选定大小作两项卷积，第二行把一项与指定素数的指示卷积。可以把残余若干卷积因子逐项限于素因子小于 $X^{u_0}$ 的整数，因为在给定乘积处这等于对原筛权作同一限制。于是合并系数中 $q$ 是唯一一个至少为 $X^{u_0}$ 的素因子，有界性保留。模数均为奇数。素数端用三阶 Heath–Brown 恒等式和相同的小素数预筛。分解项一旦有变量在 $[X^{2\varepsilon_0},X^{1/3+O(\varepsilon_0)}]$，上表和命题 4.4 就适用。短边的 Siegel–Walfisz 条件、对数权分部求和及长边的除数界正是原素数分解中的条件。其余项只有一个或两个大光滑变量，其他变量的积 $X^{O(\varepsilon_0)}$。一个大变量时用原 Type I 估计，两个时用除数函数的 $2/3-O(\varepsilon_0)$ 水平，覆盖 $5/8$。两种估计都容许上述有界模权与区间端点。分箱和二倍分段的对数幂损失可由原余项吸收。[Pascadi，定理 1.3(i) 证明中的三类项](https://arxiv.org/html/2505.00653v2#S4) · [Maynard，除数函数的区间估计](https://arxiv.org/pdf/2006.07088#page=7)

格的主项有 $\sum_q1/(q-1)$，极限为 $\int du/u$。残余筛参数 $s=(\theta(u)-u)/u$ 在 $[1,3]$，上界函数 $F(s)=2e^\gamma/s$，故每格扣除量趋于 $du/[u(\theta(u)-u)]$。取细格和严格余量的极限即可。从 $1/7$ 积到 $1/4$ 得 $2\log2$，从 $1/4$ 积到 $c$ 得 $\log(3c/(2-5c))$，这就证明 (6)。这也是 Ratliff 修订稿数值程序中 $H$ 项的一半，差别来自归一化。[Ratliff v4 的基项和积分定义](https://zenodo.org/records/24441524/files/certify_integrals.py?download=1)

## 4. 扣掉三个素因子，再扣指定区间的半素数

因为 $c>1/4$，通过粗筛的 $m$ 至多含三个素因子。含平方的行总数至多

$$
 \sum_{X^c\le j\le\sqrt{2X}}\bigl(O(hX/j^2)+1\bigr)
 \ll hX^{1-c}+\sqrt X.
$$

对三个不同素因子，按大小排列其对数比例为 $u\le v'\le w,\ u+v'+w=1$。写

$$
 \lambda(t)=\max\!\left(\frac35,\min\!\left(1-t,\frac{2+3t}{5}\right)\right),
 \qquad l_c(w)=\max(c,1-2w).
$$

筛相邻素数端所得的三因子上界系数为

$$
 \begin{aligned}
 T(c)&=\iint_{\substack{c\le u\le v'\le w\\ w=1-u-v'}}
           \frac{du\,dv'}{u v' w\,\lambda(w)}\\
 &=\int_{1/3}^{1-2c}
     \frac{\log((1-w-l_c(w))/l_c(w))}
          {w(1-w)\lambda(w)}\,dw .
 \end{aligned}
 \tag{8}
$$

固定中、短素数的相对宽度为 $h/(\log X)^3$，将第三个因子的区间向外补齐，使整数积覆盖原物理小格。由素数定理，体积是 $hX/\log X$ 乘 $\iint du\,dv'/(uv'w)+o(1)$，其中 $o(1)$ 放在归一化体积内。排序只数一次。因子相等边界或格子与边界相交的贡献，用外扩见方的比例格再缩小到边界处理。将三个素数放进矩形后，可按需要把任意两个的指示卷积成一条绝对值至多二的序列。积的短区间由矩形覆盖实现。

对相邻整数用可分解上界线性筛，筛限为 $D^{1/2}$，筛积乘 $F(2)$ 给 $4C_2/\log D$。最小素数边 $X^u$ 满足 $u\le1/3$，Pascadi 命题 5.2 对上界线性筛的分箱权给任何严格小于 $3/5$ 的水平。其三组分配配合命题 4.4，短边就是这个实际素数，满足粗糙性和 Siegel–Walfisz 条件，长边取另两个的卷积。也可选最大素数 $X^w$ 作短边，因为 $w\le1-2c<1/2$。将任一个可分解权拆成两段，指数趋于 $w,\lambda-w$，第三段取一，命题 4.4 给 $\lambda<\min(1-w,(2+3w)/5)$。取两者较好值即 (8)。主特征项中的互素限制可移除，若有一个矩形因子 $\ell\mid d$，则该格中积为 $\ell$ 的倍数的个数至多 $X^{o(1)}(hX/\ell+1)$，且 $\ell\ge X^{c-o(1)}$。对模数的除数权除以 $\varphi(d)$ 求和仅损失对数幂。实际使用时先在固定小比例格上给指数留余量，再令格宽趋零。改变筛限 $X^c$ 为 $n^c$ 或 $m^c$ 可先把 $c$ 略增再取极限，(6)、(8) 连续。[Pascadi，命题 5.2 及其对上筛权的应用](https://arxiv.org/html/2505.00653v2#S5) · [固定差二切换时的矩形与归一化](https://lixiang90.github.io/vibemath-NumberTheory/zeros-exact-distribution/articles/gap-two-linear-sieve-cubic-inversion.html)

给定 $a\le c<2/7$ 和 $c\le r\le1/3$，再扣掉小因子在 $[c,r)$ 的双素数乘积。第 37 篇单方向的密度上界在这段用 $\lambda=3/5$。因此去符号单方向下界系数为

$$
 K(c,r)=H(c)-T(c)-\frac53\log\frac{r(1-c)}{c(1-r)} .
 \tag{9}
$$

将区间行数乘回 $\log^2 X/X$，在 $\mathcal I_Y$ 的完整对数小格求和，再除 $L_Y$。两端不足一格的部分可以舍去，完整小格的 $\sum h/L_Y\to1$。对两个方向、每个固定截点都能使用第 35、36 篇的同时均方估计。在每格中，各行的 $\{n,n+2,n(n+2)\}$ 平方类空间两两交为零，四个符号均通过的概率为 $1/16+O_A(P^{-A})$，不同行的指示协方差为 $O_A(P^{-A})$。在所有对数小格对偏差取并集，再拼成允许的物理窗口。这个步骤对按算术性质预先指定的行集同样适用，正、负方向可以在同一见证上使用。因此

$$
 \liminf \bigl(J_Y+B_{\pm,Y}([r,b])\bigr)\ \ge\
 \max\!\left(0,\frac{C_2}{4}K(c,r)\right).
 \tag{10}
$$

这里的 $C_2/4$ 来自 $4C_2/16$，真孪生在每个方向只计一次。表中固定 $c=\min(r,561/2000)$，始终严格小于 $2/7$。因此即使 $r>2/7$，也无需在 $u=2/7$ 处使用临界模数估计。Ratliff v4 在其选定截点使用了与 (9) 对应的三个筛法积分 $H,I,J$，另作了有符号的 Buchstab 改进。本篇用基础积分对一组截点和两个方向同时取界，再把它们作用到原几何窗口的三个零点观测。[v4 的五个积分区与系数](https://zenodo.org/records/24441524/files/certify_integrals.py?download=1)

## 5. 三阶读数的数值下界

取 $e=117/400$，$r_j=a+(e-a)j/128,\ 0\le j\le128$，$k(r)$ 表示 (10) 在表中选定 $c$ 后的右边。令 $\alpha_0=f(a),\ \alpha_j=f(r_j)-f(r_{j-1})$。对每个 $u$，$\sum_{r_j\le u}\alpha_j\le f(u)$。所以在每一方向，

$$
 C_{\pm,Y}+f(e)J_Y
 \ \ge\ \sum_{j=0}^{128}\alpha_j\bigl(B_{\pm,Y}([r_j,b])+J_Y\bigr)
 \ \ge\ \sum_{j=0}^{128}\alpha_j k(r_j)-o(1).
 \tag{11}
$$

用复球区间积分在 $w=(1-c)/2,3/8,2/5$ 处分开计算 (8)。例如 $H(a)\approx0.284963,\ T(a)\approx0.124792$，扣除三素因子行后还剩约 $0.160171$，再乘 $C_2/4$ 就得到首行下界。$C_2$ 的部分乘积取到 $10^6$，尾乘积在 $1-1/999999$ 与 $1$ 之间。两个方向的 (11) 右端之和在 $(0.0214420,0.0214422)$ 内。由于 $f(e)=0.4138875$，(2) 随即成立。与前篇的匹配均值合并，亦有

$$
 \liminf(\overline C+0.827775 J_Y)>0.0214420.
 \tag{12}
$$

作为比较，单用 $r=a$ 的精确下界再乘 $f(a)$，两方向合起来只给约 $0.020973$ 的三阶下界（在 $J_Y\to0$ 时）。(11) 增加的部分来自较平衡的半素数必须承担一部分质量。

## 6. 怎样使用外部累积计数

还有一项独立的转移可用。设一方向在 $t$ 以下的候选累积数为 $A_r(t)$，已知对每个充分大的 $t$，

$$
 A_r(t)\ge(2C_2b_r-o(1))\,t/\log^2 t,\qquad
 A_r(t)\ll_r t/\log^2t .
$$

这里的候选是素数 $p$ 和 $m=p\pm2$，$m$ 至多含两个素因子且 $P^-(m)>m^r$。上界由平方根以下固定正幂的二维上界筛给出，或由孪生上界加第 37 篇的半素数上界给出。若原计数允许素数平方，其个数为 $O(\sqrt t)$，可以去掉。则在 (10) 相同的共同好窗口上，

$$
 \liminf\bigl(J_Y+B_{\pm,Y}([r,b])\bigr)\ge C_2 b_r/8 .
 \tag{13}
$$

这是由物理窗口的幂次宽度带来的。取 $U=2Y,\ V=Y^{3/2}/2,\ W(t)=\log^2t/t$。分部求和给 $W(V)A_r(V)-W(U)A_r(U)+\int_U^V(-W'(t))A_r(t)\,dt$。边界为 $O_r(1)$，正积分至少为 $(2C_2b_r-o(1))\log(V/U)-O_r(1)$。移两个整数、换回行权的误差除以 $L_Y$ 后都趋零，且 $\log(V/U)/L_Y\to1$。最后对这个固定候选集作同样的 $1/16$ 符号转移。整个论证只要求累积计数的下界，可以用于还没有逐个对数小格下界的定理。

Ratliff 于 2026 年 11 月 22 日发布的 v4 修订稿对 $m=p-2$ 给出下列参数。作者已撤回 v3 关于 $1/3$ 与 $0.3844197\ldots$ 的原主张。本文把 v4 的累积计数作为此段额外输入，(2) 和第一张表只用第 3、4 节的筛法。本文复算的是 (8) 至 (13) 的常数，未独立重跑 v4 的大规模筛权证书。[修订稿及作者更正](https://zenodo.org/records/24441524) · [v4 四组下界的证书汇总程序](https://zenodo.org/records/24441524/files/certify_ledger.py?download=1)

| $r$ | v4 的 $b_r$ | (13) 给负方向的严格下界 |
| ---: | ---: | ---: |
| $0.296$ | $0.0766$ | $0.0063210$ |
| $0.3000$ | $0.0129$ | $0.0010645$ |
| $0.3008$ | $0.00025$ | $0.0000206300$ |
| $0.300816$ | $0.000002$ | $0.00000016504$ |

给 (11) 的网格加上这四个截点。对负方向的每个 $r_j$，同时取所有 $r\ge r_j$ 的 (13) 与 (10) 的最大值。正方向仍用 (10)。按 (11) 的共同网格计算，得到 $\liminf(\overline C+0.841302936576 J_Y)>0.0214754$。较大的系数来自用 $2f(0.300816)$ 同时付两个方向的孪生权。

## 7. 数值核对与一个仍然容许的分布

有限样本用四符号见证 $(2965867,12456641399)$，直接分解端点，得到下表。数值用于核对方向、因子比例和物理归一化。

| $Y$ | $r$ | $J_Y$ | $J_Y+B_{+,Y}([r,b])$ | $J_Y+B_{-,Y}([r,b])$ |
| ---: | --- | ---: | ---: | ---: |
| $8192$ | $3/11$ | $0.0642431$ | $0.1207345$ | $0.1230795$ |
| $8192$ | $117/400$ | $0.0642431$ | $0.1162696$ | $0.1192159$ |
| $32768$ | $3/11$ | $0.0644709$ | $0.1275123$ | $0.1262471$ |
| $32768$ | $117/400$ | $0.0644709$ | $0.1220582$ | $0.1209062$ |

另取 $800<d\le12000,\ d\equiv23\bmod48$ 的全部 82 个素导子，沿用第 33 篇的二次 theta 积分与辐角计数，共得 4762 个数值零点。取中心 $0,0.09,0.18$，对 $n=65,71,77,83,89$ 的两个方向直接按各端点对模 $d$ 的 Euler 判别值分组。二十个乘积读数中，两端一阶读数的乘积与有限族逐素数显式公式最多相差 $0.00004480$，一阶乘二阶的读数最多相差 $0.00006222$。小导子实验使用相互重叠的宽截尾窗口检验两侧实现。渐近的分离性与统一误差由 (3)、(4) 给出。[完整复算程序](https://github.com/lixiang90/vibemath-NumberTheory/blob/main/zeros-exact-distribution/gap-two-explicit-chen-zero-lower-bound/verify.py)

可以具体检验目前离孪生无穷性还有哪一步。令抽象极限质量 $J=0$，两个方向各在 $[1/3,2/5]$ 上放常密度 $0.9$，每个方向得到质量 $0.06$。第 37 篇的单向允许密度为 $C_2/[4\lambda(u)u(1-u)]>1$。对本文的 $a\le c<2/7,\ c\le r\le1/3$，(10) 的右边至多 $C_2H(a)/4<0.048$，(13) 的右边也小于 $0.06$。故这个分布同时满足所有这些上下界。它给出的是不等式的相容性检验。要证明孪生无穷，仍须在真实的导子族上压低至少一个方向的半素数读数，或者给出另一个能排除这种分布的算术约束。式 (2) 把前一条路所需的量和固定阈值都明确了。
