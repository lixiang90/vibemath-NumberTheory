---
title: "陈氏半素数的因子比例与固定差二的三阶零点阈值"
date: "2026-11-25"
abstract: "把固定差二的半素数按较小因子的对数比例分开，求出其加权分布的显式上界。在原共同好窗口中，这把判别孪生素数所用的三阶零点阈值从直线提高为一条有显式系数的曲线。作为纯筛法推论，每个大二倍区间保留正比例的陈氏下界时，可以同时剔除最偏斜和最接近平方的两段因子比例。还可在高度稀疏的预设二次符号族上保留好窗口，并改变中心零点和的四阶累积量。"
---

# 陈氏半素数的因子比例与固定差二的三阶零点阈值

最近三篇中，第 32 篇把差二信息集中到一个补弧频率，但该频率上的主项还需下界。第 33 篇用二阶、三阶条件零点矩分别识别真孪生与半素数，二者之和已有陈氏定理给出的下界。第 34 篇找到了三阶统计中不会随高度平均消失的一项，并算出了每个半素数在其中的权重。就改善孪生素数判据而言，我选择后两篇相接的这条路线，因为可以用独立的算术估计约束这个权重。[固定频率](https://lixiang90.github.io/vibemath-NumberTheory/zeros-exact-distribution/articles/gap-two-translated-quadratic-characters.html) · [条件零点矩](https://lixiang90.github.io/vibemath-NumberTheory/zeros-exact-distribution/articles/gap-two-quadratic-conductor-covariance.html) · [匹配相位](https://lixiang90.github.io/vibemath-NumberTheory/zeros-exact-distribution/articles/gap-two-coherent-zero-heights.html)

具体地，若一端是素数、另一端是 $st$，$s<t$ 为素数，写 $\theta=\log s/\log(st)$。三阶统计的权重是 $f(\theta)=2\theta(1-\theta)$。在允许范围 $[3/11,1/2]$，其最小值为 $48/121$。上篇把这个最小值用于每一行。下面证明，若半素数有固定正质量，它们在最小值附近放不下。这给出一个严格高于旧阈值的零点判据，还给出一个不涉及零点的陈氏定理推论。关于原导子与零点的结论使用归档的 Dirichlet-GRH。[归档定理](https://github.com/leanprover-community/riemann-hypothesis#theorem-statements)

## 1. 固定差二的三个量

沿用 $\zeta$ 的原构造。固定 $0<\epsilon<1/52$ 及原证书的 $J$，令 $H=P^J,\ z=P^{\epsilon/5},\ z'=z(\log P)^2$。先取素数 $e^{z/2}<Q\le e^z,\ Q\equiv3\bmod8$，再取素数 $e^{z'/2}<d_0\le e^{z'}$，满足 $d_0\equiv7\bmod16$，模 $Q$ 以及模每个不超过 $z$ 的奇素数都余 $-1$。原几何步骤为这对素数选择 $[Y,Y^2]\subset[P,H]$。写 $\chi_v(n)=(n/v)$，$G_n$ 指示 $n,n+2$ 在 $Q$ 处均为二次剩余、在 $d_0$ 处均为二次非剩余。[原候选及扭转](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/MeanSquareBounds/PrimeSetTwoBounds.lean#L37) · [取对后选窗](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/CharacterBoxes/WinPairs.lean#L32)

令 $a=3/11$，$I_Y=(2Y,Y^{3/2}/2]$，$L_Y=\sum_{Y\le j\le Y^{3/2}/2}1/j\sim(\log Y)/2$，分母对整数求和。保留 $n\in I_Y$ 中 $n(n+2)$ 平方自由且每个素因子都大于 $n^a$ 的行。$\mathcal T$ 为两端均为素数的行，$\mathcal B$ 为一端是素数、另一端是不同素数之积 $st$ 的行，两个方向均计。令 $w(n)=\log n\log(n+2)/n,\ v_n=G_n w(n)/L_Y$。定义

$$
 J_Y=\sum_{\mathcal T}v_n,\qquad b_Y=\sum_{\mathcal B}v_n,\qquad
 A_Y=J_Y+b_Y,\qquad M_Y^0=\sum_{\mathcal B}f(\theta_n)v_n .
 \tag{1}
$$

用上标 $\circ$ 表示把 $G_n$ 换成 $1$。若 $st=n+2$，则 $\theta_n$ 可能比 $a$ 小 $O(1/(n\log n))$。在定义有限测度 $\nu_Y=\sum_{\mathcal B}v_n\delta_{\max(a,\theta_n)}$ 时将它截到 $a$。这一处理给 $M_Y^0$ 带来的误差是 $O(1/(Y\log Y))$。

## 2. 半素数能占据多窄的因子比例

令 $C_2=\prod_{p>2}(1-(p-1)^{-2})$，乘积对素数取。对任何区间 $I\subset[a,1/2]$，一致有

$$
 \nu_Y^\circ(I)\le 32 C_2\int_I\frac{du}{u(1-u)^2}
               +O\!\left(\frac{\log\log Y}{\log Y}\right).
 \tag{2}
$$

这个系数包括两个方向。以下推导保留常数，因为两阶段符号筛后可得到足够小的数值界。

取一个短段 $U<n\le(1+h)U$，其中 $h=1/\log Y$，$U$ 在物理窗口的范围内。固定较小素因子 $s$，则较大因子 $t$ 和 $st+2$ 或 $st-2$ 必须同时为素数。把奇数 $t=2j+1$ 代入，$j$ 的区间长度为 $W=hU/(2s)+O(1)$。对奇素数 $\ell\ne s$，两线性式正好排除两个不同剩余类。筛水平取 $Z=\sqrt W/\log W$，所数的两个素数都大于 $Z$。Montgomery–Vaughan 大筛推论 1 给每个方向的计数上界

$$
 \frac{W+O(Z^2+1)}{L_s(Z)},\qquad
 L_s(Z)=\sum_{\substack{d\le Z\\(d,2s)=1}}
          \mu^2(d)\prod_{\ell\mid d}\frac{2}{\ell-2}.
$$

仅排除 $2$ 的分母 $L_*(Z)=(\log Z)^2/(8C_2)+O(\log Z)$。事实上其系数的 Dirichlet 级数是 $\zeta(1+v)^2F(v)$，$F(0)=1/(4C_2)$。在奇素数处，$F$ 的因子是 $(1+2\ell^{-v}/(\ell-2))(1-\ell^{-1-v})^2$，在 $2$ 处是 $(1-2^{-1-v})^2$。乘积在零点左边的一个小邻域绝对收敛。与 $\sum_{m\le x}\tau(m)/m=(\log x)^2/2+O(\log(2x))$ 作卷积即得所述误差。再排除 $s$ 至多减去 $2L_*(Z)/(s-2)$，而 $s\gg U^a$。所以每个方向的上界是

$$
 \left(16C_2+O\!\left(\frac{\log\log U}{\log U}\right)\right)
 \frac{hU}{s(\log(U/s))^2}.
 \tag{3}
$$

此处从 $32C_2$ 到 $16C_2$ 的因子 $1/2$ 来自 $t$ 为奇数。$\log(U/s)\ge(\tfrac12-o(1))\log U$ 保证误差对 $s$ 一致。[Montgomery–Vaughan，推论 1](https://personal.science.psu.edu/rcv4/personal/Publications/large_sieve.pdf#page=2)

在短段上 $w(n)\le(1+O(h+1/U))\log^2 U/U$。若 $I=[u,v]$，则 $\log s/\log U$ 落在 $[u,v]$ 向两端各延长 $O(1/\log U)$ 的区间内，这也覆盖在 $a$ 截断的行。由素数倒数的 Mertens 公式作分部求和，对这些 $s$ 求 $s^{-1}(1-\log s/\log U)^{-2}$ 的和，上界为 $\int_u^v dt/[t(1-t)^2]+O(1/\log U)$，也适用于退化区间。把 (3) 的两个方向相加，每个短段的上界是 $h$ 乘以 (2) 的右边。用公比 $1+h$ 的短段覆盖 $I_Y$，末端可补齐，所有 $h$ 之和为 $(1+O(1/\log Y))L_Y$。这证明 (2)。对 $U<n\le2U$ 用同一短段证明并去掉 $w$，得到单方向的区间计数 $\ll (|I|+\log\log U/\log U)U/\log^2U$。

**一个纯算术推论。** 存在固定的 $\eta>0,c_0>0$，每个充分大的 $(U,2U]$ 中至少有 $c_0U/\log^2U$ 个素数 $p$，满足 $p+2$ 是素数，或 $p+2=st$ 且

$$
 (p+2)^{3/11+\eta}<s<(p+2)^{1/2-\eta},\qquad s<t\ \text{均为素数}.
 \tag{4}
$$

Green–Tao 所引的定量陈氏定理先给下界 $c_1U/\log^2U$。平方端至多有 $O(\sqrt U)$ 个。用刚证明的单方向上界删去 $[a,a+\eta]$ 与 $[1/2-\eta,1/2]$，损失至多 $(C\eta+o(1))U/\log^2 U$。取 $0<\eta<\min((1/2-a)/3,c_1/(2C))$ 即可。这个推论同时保留了真孪生和因子比例落在内段的陈氏行。$\eta$ 未作数值优化。[Green–Tao，定理 6.1 及脚注 3](https://arxiv.org/pdf/math/0405581#page=19)

## 3. 原来的好窗口与预设二次符号

先核对选窗的一个量词。在最终的 Dirichlet-GRH 下，$\zeta$ 构造的**每一对**上述同余候选最终都满足原来的解析条件 $\mathrm{PairGood}(1,\epsilon,P,H,Q,d_0)$。原条件涉及本原扭转，至少含 $Q,d_0$ 之一。相应新素数处的非平凡分量保证特征非主，导子 $N$ 整除 $8Qd_0$。所需无零矩形的左边界为 $1-1/400$。令 $\psi_\chi(u)=\sum_{n\le u}\Lambda(n)\chi(n)$。即使 $u<N$，也有 $|\psi_\chi(u)|\ll\sqrt u\log^2(3Nu)$。在截断显式公式中取高度 $u\ge2$，GRH 给 $|u^\rho|=\sqrt u,\ |\rho|\ge1/2$，在单位高度段中用 $O(\log(N(|t|+2)))$ 的零点计数累加即可。对由本原特征诱导的模 $m$ 特征，漏掉的 Euler 因子对 $\psi$ 贡献 $O(\log(2m)\log(2u))$，同样可吸收。原中心化斜坡满足

$$
 G_\chi(y)=-y\int_y^\infty\psi_\chi(u)\frac{du}{u^2},\qquad
 |G_\chi(y)|\ll\sqrt y\,\log^2(3Ny)\le y^{1/2+\epsilon/2}
 \quad(P\le y\le H).
 \tag{5}
$$

最后一步用 $(z')^2=o(P^{\epsilon/2})$。又把 $-L'/L(1,\chi)$ 的级数截在 $P$，分部求和给尾部 $O(P^{-1/2}\log^2(3NP))=o(1)$，前段绝对值至多 $\sum_{n\le P}\Lambda(n)/n=(1+o(1))\log P$。这给原来的 $4\log y$ 界。几何步骤随后仍对所选的每一对使用它自己的有限例外集合。这里把最终 GRH 用于研究原构造的自由度。[PairGood 的三个条件](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/CharacterBoxes/PairGood.lean#L15) · [原讲义的零点计数、截断显式公式 Z2–Z3 及其后的斜坡积分 III.23](https://www.math.uwaterloo.ca/~mrubinst/alpoge_claude_grh/zeta-paper_v2.pdf#page=67) · [几何例外与选窗](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/AutomorphicVarieties/LayeredIsomorphisms/ThmBabOf.lean#L72)

还可以预先给 $d_0$ 开出一整列二次符号。固定 $\gamma>0,\ \gamma\log2<1/2$，令 $T=\gamma z'\log z'$。在所有 $z<p\le T$ 的素数上任意指定 $\chi_{d_0}(p)=\sigma_p$，记其个数为 $r$。当 $P$ 充分大时 $T<P^a$。每个 $Q$ 下，这个子族在原候选中的比例为

$$
 2^{-r}+O(X^{-1/2+o(1)}),\qquad X=e^{z'},\qquad r\sim\gamma\log X .
 \tag{6}
$$

证明时，模 $M=16Q\prod_{3\le p\le z}p$ 的原候选是一个可逆类，乘积只取素数。对与 $M$ 互素的非平凡奇平方类 $u$，二次互反及模 $M$ 的特征展开使素数和中的每个扭转均为非主。GRH 给归一化偏差 $\ll\varphi(M)\log X\,\log^2(3XMu)/\sqrt X$。可用的 $u$ 包括任意一组指定素数之积再乘固定次的物理或观测标签，再取平方自由核。因为 $\log M=O(z),\ \log u=O(T+z+\log H)$，该偏差是 $X^{-1/2+o(1)}$。展开 $\prod_{z<p\le T}(1+\sigma_p\chi_{d_0}(p))/2$，所有系数的绝对值之和为 $1$，即得 (6)。对完全位于 $T$ 以上、与 $M$ 互素的非平凡相关标签，给定这些符号后的偏差至多 $O(2^rX^{-1/2+o(1)})$，仍小于任意固定 $P^{-D}$。第一阶段模 $8$ 的对应估计对物理标签同样成立。[两阶段特征和的推导](https://lixiang90.github.io/vibemath-NumberTheory/zeros-exact-distribution/articles/gap-two-quadratic-conductor-covariance.html)

**符号筛后的因子比例。** 先均匀取 $Q$，再对每个 $Q$ 均匀取满足上述指定的 $d_0$。也可以取全部原候选。在每一种预设下，以概率 $1-o(1)$，对所有 $P\le Y,\ Y^2\le H$ 及所有 $I\subset[a,1/2]$ 同时有

$$
 A_Y\ge c>0,\qquad
 \nu_Y(I)\le 2C_2\int_I\frac{du}{u(1-u)^2}+O(e_P)
       \le 11|I|+O(e_P),\qquad e_P=\frac{\log\log P}{\log P}.
 \tag{7}
$$

常数与预设的符号值无关，$\gamma,\epsilon,J$ 固定。无任何见证概率限制时，由 $0\le G_n\le1$ 及 (2) 可把最后的 $11$ 换成 $175$。确实，$1/[u(1-u)^2]$ 在所述区间的最大值为 $8$，仅用 $3,5,7$ 就给 $C_2\le175/256$。两阶段筛保留的 $1/16$ 将 $256C_2$ 变为 $16C_2<11$。

证明 (7) 不能简单地将无条件事件与 (6) 的稀疏事件相交。应在指定符号之后重新算行的均方。每行的非平凡平方类是 $n,n+2,n(n+2)$，均由大于 $T$ 的素数生成。在同一个短段中，不同行的端点不相同，因为保留行均有 $n\equiv5\bmod6$。两端乘积又比该段任一端点大，且随 $n$ 严格递增。故任意两行的二维平方类空间只有零向量相交。在每个阶段的独立素数符号模型中两行的端点符号独立。用 (6) 后的特征估计依次比较两阶段，得对同一短段内的任意子集和权 $w_n$，

$$
 \mathbb E\left|\sum_n w_n(G_n-1/16)\right|^2
 =\frac{15}{256}\sum_n |w_n|^2+
 O_D\!\left(P^{-D}\big(\sum_n |w_n|\big)^2\right).
 \tag{8}
$$

为允许筛后选 $Y$，把 $[P,H+2]$ 用公比 $1+1/\log P$ 的网格覆盖，角度另用步长至多 $1/\log P$、点数 $O(\log P)$ 的网格。在每一物理格、每一角度前缀对半素数用 (8)，并对孪生行单独用一次。每格中未归一化的 $\sum w(n)^2\ll_J(\log P)^4/P$。取误差阈值 $(\log P)^{-3}$，Chebyshev 不等式及 $O_J((\log P)^3)$ 次取并给总失败概率 $O_J((\log P)^{13}/P)$，常数也可依赖固定的 $\epsilon,\gamma$。完整网格的误差可直接相加。窗口两端的不完整格中，半素数的总权由 (3) 求和给 $O(1/\log P)$，孪生行对 $n,n+2$ 用同一个二维上界筛也给 $O(1/\log P)$。将角度两端向外取到网格后用 (2)，就得到 (7) 的测度界。同时有 $A_Y=A_Y^\circ/16+O(e_P)$。每个 $I_Y$ 含 $\gg\log Y$ 个完整二倍区间，每段上定量陈氏定理给 $\gg1$ 的未归一化权，故有一致的 $c>0$。原几何证明的有限窗口覆盖还允许在每种预设的族内选出一段，供固定正比例的见证共同使用。[原共同窗口的有限覆盖](https://github.com/lixiang90/vibemath-NumberTheory/blob/main/Goldbach_TwinPrimes/twin-prime-auxiliary-signs.md)

## 4. 从因子比例到更高的零点阈值

记 $\beta=f(a)=48/121$。对密度上界 $D$ 定义严格递增的函数

$$
 F_D(x)=\int_0^x f\!\left(\min\left(\tfrac12,a+\tfrac{t}{D}\right)\right)dt .
$$

在 $0\le x\le D(1/2-a)$ 上，两端各自给出的极值函数为

$$
 F_D(x)=\beta x+\frac{1-2a}{D}x^2-\frac{2}{3D^2}x^3,\qquad
 U_D(x)=\frac{x}{2}-\frac{2}{3D^2}x^3 .
 \tag{9}
$$

在第 3 节的共同好窗口取 $D=11$。写 $b=b_Y$，则

$$
 F_D(b)-O(e_P)\le M_Y^0\le U_D(b)+O(e_P).
 \tag{10}
$$

式 (7) 给 $b\le D(1/2-a)+O(e_P)$，端点外的 $O(e_P)$ 差别可直接吸收。在任意原好见证上，(10) 仍成立，可取 $D=175$。证明可按质量从左到右排序。设总质量的第 $t$ 个单位处的比例为 $\vartheta(t)$，$0<t<b$。若所有区间的质量不超过 $D|I|+E$，则

$$
 a+(t-E)_+/D\le\vartheta(t)\le
 \tfrac12-(b-t-E)_+/D .
$$

这里 $x_+=\max(x,0)$，广义分位点按几乎处处取值。对递增的 $f$ 从 $0$ 积到 $b$，$E=O(e_P)$ 造成 $O(e_P)$ 的误差。两边分别把密度 $D$ 紧靠左端和右端排列，积分就是 (9)。因此在只知道常数密度上界时，这两条曲线均可取到。保留 (7) 中随 $u$ 变化的密度还可用同一分位数公式进一步收紧。

把它写回零点上。取 $B=z/8,\ K_*=\lceil(\log P)^2\rceil,\ b_*=B/(4K_*),\ R=2/b_*$，以及 $g=\mathbf1_{[-B,B]}*((2b_*)^{-1}\mathbf1_{[-b_*,b_*]})^{*2K_*}$。对观测导子 $d$ 定义 $S_d(h)=\sum_{|\gamma_d-h|\le R}\widehat g(\gamma_d-h)$，零点为 $1/2+i\gamma_d$，按重数计，Fourier 相位为 $e^{-itu}$。令 $T_d(h)$ 为 $-(S_d(h)-\log d)/2$ 在观测族内的中心化。固定 $N=n(n+2)$，以 $\Delta_3(h_1,h_2,h_3;N)$ 记条件 $\chi_d(N)=+1$ 与 $-1$ 时三个 $T_d$ 的联合三阶累积量之差的一半，正组在前。物理见证 $Q,d_0,Y$ 先固定，观测 $d$ 遍历相同 $Q$ 下的第二阶段候选。设 $h_k=2\pi k/(8\log Y),\ \xi=3R$。第 34 篇所用的匹配三阶量为

$$
 C_k=\frac1{L_Y}\sum_{n}
 \frac{2G_n\cos(h_k\log n)}{\log n}\,
 \Delta_3(0,\xi,h_k;n(n+2)),
 \tag{11}
$$

求和遍历第 1 节的全部平方自由粗糙行。三个零点窗口两两分离。对固定 $0<\sigma<a$，任选 $m$ 个不同的 $1\le k\le Y^\sigma$，记其平均为 $\overline C$，则前篇式 (6)、(7) 给

$$
 |\overline C-M_Y^0|\ll_\sigma
 m^{-1/2}+(R\log H)^2+Y^{-a}.
 \tag{12}
$$

这连同 (10) 给出新判据。存在常数 $K_0$，令 $E_Y=K_0(e_P+m^{-1/2}+(R\log H)^2+Y^{-a})$，在上述共同好窗口上有

$$
 \boxed{\displaystyle
 J_Y\ge\left[A_Y-F_{11}^{-1}\!\left((\overline C+E_Y)_+\right)\right]_+.}
 \tag{13}
$$

任意固定 $x>0$ 都有 $F_{11}(x)>\beta x$。若 $A_Y\ge c$，$m\to\infty$，且在无穷多个共同好窗口上 $\overline C\le F_{11}(c)-\delta$，其中 $\delta>0$ 固定，则 $J_Y\ge2\delta-o(1)$，从而原归一化的 $\Gamma_1\ge2\delta/3$。因 $F_{11}'\le1/2$，这个最后的下界直接来自 (13)。原来需要把匹配均值压到 $\beta c$ 以下。如果只求无穷性，$F_{11}(A_Y)-\overline C>E_Y$ 在无穷多个窗口成立就足以保证每窗至少一对，盈余可以趋零。证明至此仍缺这些窗口上的**匹配三阶均值上界**。补弧的失配小值、陈氏下界或零点主方差单独都不给这条估计。[共用高度的平方和界](https://lixiang90.github.io/vibemath-NumberTheory/zeros-exact-distribution/articles/gap-two-coherent-zero-heights.html) · [四符号质量与 $\Gamma_1$](https://github.com/lixiang90/vibemath-NumberTheory/blob/main/Goldbach_TwinPrimes/twin-prime-auxiliary-signs.md)

## 5. 改变哪些零点数据可以保留同一判据

式 (6) 还允许把**观测**导子 $d$ 全部限制到指定子族，物理见证与观测平均仍分开取。记

$$
 a_h(p)=\frac{\log p}{\sqrt p}\sum_{j\ge0}p^{-j}
    g((2j+1)\log p)\cos((2j+1)h\log p).
$$

显式公式把全零点和化为导子项、族内相同的无穷处与偶次幂项，以及 $-2\sum_p a_h(p)\chi_d(p)$。截到半宽 $R$ 的误差对 $|h|\le H+1$ 为 $O_A(P^{-A})$。对 $p\le z$，原同余已固定 $\chi_d(p)=+1$。若再将 $z<p\le T$ 全部指定为 $+$ 或全部为 $-$，就得到两族零点均值的位移，以及任一指定族中心窗的四阶累积量

$$
 \begin{aligned}
 \mathbb E_+[S_d(0)-\log d]-\mathbb E_-[S_d(0)-\log d]
   &=-4\sum_{z<p\le T}a_0(p)+O_A(P^{-A})
      \sim-8(\sqrt T-\sqrt z),\\
 \kappa_4(T_d(0))&\sim-2\,\frac{(\log T)^3}{T}.
 \end{aligned}
 \tag{14}
$$

第二式未知符号的部分从 $p>T$ 开始。固定阶矩可与独立正负号比较，因为每个非主标签的条件偏差为 $O(2^rX^{-1/2+o(1)})$，且 $\sum|a_h(p)|\le\exp(3B/4+O(\log B))=X^{o(1)}$。使用截尾误差时可用 Hölder 不等式，线性和的每个固定阶范数是 $O_j(B+1)$，这由其偶次矩和 $\sum_p a_h(p)^2\ll B^2$ 得到。独立符号的四阶累积量为 $-2$，故主项为 $-2\sum_{p>T}a_0(p)^4$。在核的平台上 $a_0(p)=(\log p/\sqrt p)(1+O(1/p))$，平台外的四次和指数小，素数定理给 (14)。其绝对值比原族以 $z$ 为下限时缩小了约 $T/z\asymp(\log P)^3$ 倍。[显式公式与原族的累积量](https://lixiang90.github.io/vibemath-NumberTheory/zeros-exact-distribution/articles/gap-two-quadratic-conductor-covariance.html) · [Hughes–Rudnick，第 2 节](https://arxiv.org/pdf/math/0208230#page=4)

每个物理行的全部素因子都大于 $T$。再按 $\chi_d(N)=\pm1$ 分组时，两组在指定族中的比例仍各为 $1/2+O_A(P^{-A})$。在独立模型中，中心化后只保留 $p>T$ 的线性和，三阶累积量半差仍只在 $N$ 恰有三个素因子时留下三个算术系数的排列和。两因子情形可同时翻转所有自由符号消去三阶项，更多因子的标记则与次数三的多项式正交。矩比较误差可取到小于任意 $P^{-A}$，对所用高度一致。因此第 34 篇在排列和上作的相位估计仍给 (12)，(13) 也成立。这样可以显著改变 (14) 的零点统计，但这个操作对所需的匹配三阶统计仍给同一个算术中心 $M_Y^0$，至消失误差。它指出继续选扭转时需要控制的是涉及大素因子本身的相关。

## 6. 有限复算

先用有理数在六个等宽比例区间中取密度 $0,11/2,11$，逐一积分全部 $729$ 个密度函数，(9) 的上下界均成立，并各有七个取等样本。再取 $Q=2965867,\ d_0=12456641399$，逐行分解整数、应用四个实际二次符号，计算 (1) 的量。下表将**算术中心值** $M_Y^0$ 代入旧、新反函数，不计渐近误差，便于比较判据本身。

| $Y$ | 筛后真孪生 / 半素数行 | 实际 $J_Y$ | $M_Y^0$ | $[A_Y-M_Y^0/\beta]_+$ | $[A_Y-F_{11}^{-1}(M_Y^0)]_+$ |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1024 | 8 / 31 | 0.0267754 | 0.0630466 | 0.0039710 | 0.0064670 |
| 8192 | 192 / 367 | 0.0642431 | 0.0535546 | 0.0445688 | 0.0463838 |
| 32768 | 1218 / 2266 | 0.0644709 | 0.0579875 | 0.0431116 | 0.0452318 |

最后单独检验符号选择对零点的作用。取 $800<d\le24000,\ d\equiv23\bmod48$ 的素导子。以 $B=1.55,\ b_*=0.30$ 及四个均匀核卷积，算术支撑止于 $e^{2.75}<17$。除已固定的 $2,3$ 外，支撑内恰好只有 $5,7,11,13$ 四个素数。四处全为 $+$ 与全为 $-$ 的导子分别有 $10,9$ 个。用二次 theta 积分求零点并与辐角计数核对，合计使用 $|\gamma|\le20.76748$ 内的 $1158$ 个零点。两族 $S_d(0)-\log d$ 的均值差实测为 $-1.65451925$，由四个素数独立算出的全零点公式预测为 $-1.65448394$，两个族内各自的最大极差均小于 $0.000009$。这是浮点截尾实验。实际几何窗口由第 3 节保证，渐近阈值需要第 4 节的误差项。[可复算的程序](https://github.com/lixiang90/vibemath-NumberTheory/blob/main/zeros-exact-distribution/gap-two-chen-factor-ratios/verify.py)
