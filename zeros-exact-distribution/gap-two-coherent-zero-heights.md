---
title: "共同好窗口内三段分离零点的相关与半素数的端点相位"
date: "2026-11-25"
abstract: "对固定差二，在每个原幂次窗口取一组所有整数对共用的零点高度。对任意固定 σ 小于 3/11，在 Y 的 σ 次幂个高度上，匹配端点相位后的二阶和三阶统计分别集中在真孪生质量和加权半素数质量附近，总平方偏差有界。三段零点窗两两不交。这给出一个可在观察数据后挑选低值高度的孪生素数判据。另给出原四符号筛选的精确二阶主项，以及相位偏离匹配值时的响应公式。"
---

# 共同好窗口内三段分离零点的相关与半素数的端点相位

第 33 篇把孪生素数与半素数邻项分别写成二阶、三阶条件零点统计。陈氏定理给两者的总质量一个固定正下界。要进一步估计三阶部分，我们考察第三个零点窗移到许多不同高度时的行为。若统计经常很小，能否断定半素数贡献很小？

答案取决于计算时保留的相位。若在第 $k$ 个高度乘回整数端点 $n$ 的相位，就有一个在几乎所有这些高度上都能测到的半素数项。把相位乘上任意固定的不同倍数，三阶统计本身就在绝大多数高度上很小。下文还算出两种行为的过渡，并给出在观察数据后挑选低值高度仍然有效的孪生素数判据。零点部分使用归档的 Dirichlet-GRH。原四符号筛选的均方上界还有一个直接使用原扭转引理的证明。[前篇的条件协方差与三阶累积量](https://lixiang90.github.io/vibemath-NumberTheory/zeros-exact-distribution/articles/gap-two-quadratic-conductor-covariance.html) · [归档的 Dirichlet-GRH 陈述](https://github.com/leanprover-community/riemann-hypothesis#theorem-statements)

## 1. 一个幂次窗口只取一组零点高度

沿用 $H=P^J,\ z=P^{\epsilon/5},\ z'=z(\log P)^2$，$0<\epsilon<1/52$ 和 $J$ 固定。第一素数 $Q\equiv3\bmod8$ 且 $z/2<\log Q\le z$。给 $Q$ 后，第二素数 $d$ 遍历 $z'/2<\log d\le z'$ 中满足 $d\equiv7\bmod16$、模 $Q$ 及每个不超过 $z$ 的奇素数均余 $-1$ 的候选。固定原好见证 $Q,d_0,[Y,Y^2]$，其中 $[Y,Y^2]\subset[P,H]$。原证明在取过 $Q,d_0$ 后得到至多固定个数的例外扭转指标，用这些指标避开的区间选 $Y$。这个步骤限制 $Y$ 与指标的相对大小。[原候选的定义](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/MeanSquareBounds/PrimeSetTwoBounds.lean#L37) · [先取例外集合再取好窗口](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/AutomorphicVarieties/LayeredIsomorphisms/ThmBabOf.lean#L72)

记 $I_Y=(2Y,Y^{3/2}/2]$，$L=\log Y$，$L_Y=\sum_{Y\le j\le Y^{3/2}/2}1/j$（对整数求和），$a=3/11$，$w(n)=\log n\log(n+2)/n$。取 $I_Y$ 中 $N_n=n(n+2)$ 平方自由且每个素因子均大于 $n^a$ 的行，记作 $\mathcal R_Y$。令 $G_n$ 指示两端在 $Q$ 处均为二次剩余、在 $d_0$ 处均为非剩余，$v_n=G_n w(n)/L_Y$。真孪生行组成 $\mathcal T$。恰有一端为素数 $p$、另一端为两个不同素数之积 $st$ 的保留行组成 $\mathcal B$，取 $s<t$，左右两个方向均计。写

$$
J_Y=\sum_{\mathcal T}v_n,\qquad A_Y=J_Y+\sum_{\mathcal B}v_n,\qquad x_n=\frac{\log n}{L}.
\tag{1}
$$

在每个充分大的 $P$，有固定正比例的原好见证，其每个好窗口均有 $A_Y\ge c>0$，也可取正比例见证共同合格的一段。这里用前篇定理 3，所用的算术输入是每个大二倍区间的定量陈氏定理。[Green–Tao，定理 6.1 及脚注 3](https://arxiv.org/pdf/math/0405581#page=19)

对新的观测导子 $d$，只取 $L(s,\chi_d)$，其中 $\chi_d(m)=(m/d)$。令 $B=z/8,\ K_*=\lceil(\log P)^2\rceil,\ b=B/(4K_*),\ R=2/b$。取 $g=\mathbf1_{[-B,B]}*((2b)^{-1}\mathbf1_{[-b,b]})^{*2K_*}$，$\widehat g(t)=\int g(u)e^{-itu}\,du$。对实窗心 $h$ 记

$$
S_d(h)=\sum_{|\gamma_d-h|\le R}\widehat g(\gamma_d-h),\qquad
T_d(h)=-\tfrac12\bigl(S_d(h)-\log d-\mathbb E_d[S_d(h)-\log d]\bigr).
$$

零点按重数计。$\mathbb E_d$ 对给定 $Q$ 的全体第二阶段候选取均匀平均。把这批导子按 $\chi_d(N)=+1,-1$ 分组，记 $D(h_1,\ldots,h_j;N)$ 为两组 $T_d(h_i)$ 联合累积量之差的一半，正组在前。二阶是条件协方差，三阶是组内中心化三变量的乘积均值。这些观测可以逐个用于已经固定的原几何见证。

置 $\xi=3R$。给定 $0<\sigma<a$ 和整数 $1\le K\le Y^\sigma$，所有物理行共用

$$
h_k=\frac{2\pi k}{8L},\qquad 1\le k\le K .
$$

以 $0,\xi,h_1,\ldots,h_K$ 为中心的半宽 $R$ 的各窗在充分大时两两不交，因为 $RL\to0$。在三阶统计中也只配来自三个不同区间的零点。若希望全部观测都趋近实轴，可取 $K=\lfloor\sqrt L\rfloor$。对 $0\le\lambda\le3/2$ 定义

$$
\begin{aligned}
V_k(\lambda)&=\frac1{L_Y}\sum_{n\in\mathcal R_Y}
G_n\cos(\lambda h_k\log n)\, D(0,h_k;N_n),\\
C_k(\lambda)&=\frac1{L_Y}\sum_{n\in\mathcal R_Y}
\frac{2G_n\cos(\lambda h_k\log n)}{\log n}\,D(0,\xi,h_k;N_n).
\end{aligned}
\tag{2}
$$

在相同的零点数据上改变 $\lambda$，只需改变由 $n$ 直接算出的外部权。

## 2. 原四道符号能改变多少物理质量

先固定确定性的 $Y$。对 $\mathcal R_Y$ 上任意预先给定的复系数 $u_n$，有

$$
\mathbb E_{Q,d_0}\left|\sum_n u_n(G_n-1/16)\right|^2
=\frac{15}{256}\sum_n|u_n|^2+
O_D\!\left(P^{-D}\Big(\sum_n|u_n|\Big)^2\right)
\tag{3}
$$

对每个固定 $D>0$ 成立。概率为先均匀取 $Q$，再按 $Q$ 均匀取 $d_0$。即使两行的半素数共享了一个素因子，(3) 仍给同一个系数。这把前篇按共享因子数得到的均方上界精细到单行方差。

**证明。** 在这段上，每行的三个非平凡平方类都以平方自由整数 $n,n+2,N_n$ 表示。一行的两端互素。所有保留的 $n$ 模 $6$ 余 $5$，故不同的行没有相同端点。所有 $N_n$ 单调递增，最小的 $N_n$ 又大于整个物理窗口的最大端点。因此不同两行生成的两个二维平方类空间交集为零。对独立素数符号，这使两行的四个端点符号在每一抽样阶段相互独立。单行两阶段全过的概率为 $1/16$，两行同时全过的概率为 $1/256$。真实的两阶段平均先固定 $Q$ 使用一次素数特征和界，再对 $Q$ 使用同一界。端点平方类的非平凡乘积均与相应基模互素，至多是 $P$ 的固定次幂。前篇式 (2) 的 GRH 界给每个期望任意负幂误差，即得 (3)。不能在任意跨度内略掉标签大小条件，例如跨段的两个左端 $5,35$ 满足 $5\cdot7=35$。[前篇的两阶段特征和界](https://lixiang90.github.io/vibemath-NumberTheory/zeros-exact-distribution/articles/gap-two-quadratic-conductor-covariance.html) · [Lee 的 GRH 素数特征和界](https://arxiv.org/abs/2301.13457)

只调用原归档中容许一个例外的扭转引理，(3) 的等号右边仍可改为 $O(\sum|u_n|^2)+O_D(P^{-D}(\sum|u_n|)^2)$ 的上界。事实上，对一个阶段指定同一符号 $\delta=\pm1$ 的两端指示函数减 $1/4$，得到 $(\delta\chi(n)+\delta\chi(n+2)+\chi(N_n))/4$。把所有标签合在一个和中。若例外平方类为 $e$，内积中能互相匹配的标签至多是自身及自身乘 $e$ 的平方自由核。后者是一个对合，其贡献由 Cauchy–Schwarz 控制在一份系数平方和内。第二阶段可先把已知的第一阶段指示函数乘入 $u_n$，再用一次第一阶段的同样估计。这说明原引理本身足以控制筛选造成的均方损失。[原扭转引理的固定指数与例外量词](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/MeanSquareBounds/twistK.lean#L145)

以 $A_Y^\circ$ 表示在 (1) 中把 $G_n$ 换为 $1$ 的总质量，$J_Y^\circ$ 同理。普通二维上界筛给二倍段 $U<n\le2U$ 中的保留行数 $\ll U/\log^2 U$。于是 $\sum_{\mathcal R_Y}(w(n)/L_Y)^2\ll 1/Y$。对 $A_Y$，在第一个完整内段用陈氏下界，这个平方和限制到 $\mathcal T\cup\mathcal B$ 上仍 $\gg1/Y$。因此 $A_Y-A_Y^\circ/16$ 的均方根恰为 $Y^{-1/2}$ 量级。在概率为 $\alpha$ 的共同好见证集上取条件平均，偏差至多 $O((\alpha Y)^{-1/2})$。

若 $Y$ 也依赖见证，则

$$
\left\|\sup_{\substack{Y\ge P\\Y^2\le H}}
 \left|A_Y-A_Y^\circ/16\right|\right\|_2
\ll \frac{\log P}{\sqrt P}.
\tag{4}
$$

对 $J_Y$ 以及下节的加权半素数质量也有这个上界。对 $V\le Y\le2V$，把所有可能的 $n$ 排在 $(2V,(2V)^{3/2}/2]$ 内。在这段中仍无端点与另一行的乘积重合。任意前缀拆成至多 $O(\log V)$ 个二进块，每层块互不相交。对前缀和用 Cauchy–Schwarz，再在每层用 (3)，得到最大部分和的二阶界 $\ll(\log V)^2\sum |u_n|^2$，加任意负幂误差。此处未除以 $L_Y$ 的系数以 $O(w(n))$ 为界且不依赖 $Y$，平方和为 $O((\log V)^2/V)$。取两个前缀的差，除以 $L_Y\asymp\log V$，最后对 $V=2^jP$ 的最大误差平方求和，即得 (4)。它允许在筛选后选窗，窗口内原本有多少半素数仍由算术决定。

## 3. 在多少个高度上能看见同一个量

记

$$
a_h(v)=\frac{\log v}{\sqrt v}\sum_{j\ge0}v^{-j}
 g((2j+1)\log v)\cos((2j+1)h\log v).
$$

在 $n\in\mathcal B$ 的三个素因子中，若 $v$ 为其中一个、$u_1,u_2$ 为另两个，定义

$$
b_{n,v}=\frac{G_n\log v}{L_Y\sqrt v\log n}
\big(a_0(u_1)a_\xi(u_2)+a_\xi(u_1)a_0(u_2)\big).
\tag{5}
$$

此处 $p$ 始终指整端点素数。令 $M_Y=\sum_{\mathcal B} b_{n,p}$，$\theta_n=\log s/\log(st)$。因为 $|\xi|\log H=o(1)$，三个 $b_{n,v}$ 均非负且 $\ll v_n$。有

$$
b_{n,p}=v_n\left(2\theta_n(1-\theta_n)
  +O((RL)^2+Y^{-a})\right),\qquad
\frac{48}{121}(A_Y-J_Y)-o(1)\le M_Y\le\tfrac12(A_Y-J_Y)+o(1).
\tag{6}
$$

这用到了 $\theta_n\ge3/11-O(1/(n\log n))$ 和 $\sqrt{pst}=\sqrt{n(n+2)}$。在对多项式个高度的总平方偏差作估计时，(5) 保留两个固定窗中的真实系数。

**定理 1。** 设 $\mu_2=\sum_{\mathcal T}v_n\delta_{x_n}$，$\mu_3=\sum_{\mathcal B}b_{n,p}\delta_{x_n}$，总质量分别为 $J_Y,M_Y$。以 $Z_{2,k}=V_k,\ Z_{3,k}=C_k$ 表示 (2)。对 $j=2,3$，所有原好见证、$0\le\lambda\le3/2$ 和 $1\le K\le Y^\sigma$ 一致有

$$
\sum_{k=1}^K\left|Z_{j,k}(\lambda)-
 \int\cos\!\left(\frac{2\pi k(\lambda-1)x}{8}\right)d\mu_j(x)\right|^2
\ll_{\sigma,\epsilon,J}1 .
\tag{7}
$$

特别地，$\lambda=1$ 时，二阶统计与 $J_Y$、三阶统计与 $M_Y$ 的总平方偏差有界。对于 $\lambda\ne1$，则有

$$
\sum_{k=1}^K|Z_{j,k}(\lambda)|^2
\ll_{\sigma,\epsilon,J}1+\min(K,|\lambda-1|^{-1}).
\tag{8}
$$

故允许零点高度的数目随 $Y$ 作幂次增长，对任意固定容差，越过容差的高度个数仍可有与 $Y$ 无关的上界。在 (7) 中换了端点的相位倍数时，须同时换掉被减去的中心值。

**先展开零点统计。** 第 33 篇的条件累积量公式在 $|h|\le Y^\sigma+1$ 上仍具有任意负幂精度。显式公式截尾的单位高度零点数只需把 $\log(d(3+|t|))$ 换成 $\log(d(3+|h|+|t-h|))$，而 $\log(d(3+Y^\sigma))=O(z')$。在算术一侧，余弦系数的绝对值界与 $h$ 无关。固定次数的矩比较和中心化可照用。因而在 $\mathcal R_Y$ 上，二阶统计只留 $\omega(N_n)=2$ 的行，等于两个 $a_h$ 的排列和。三阶统计只留 $\omega(N_n)=3$ 的三个 $a_h$ 的排列和。$\omega=2$ 时，三阶项在条件符号模型同时翻转全部符号下为奇函数，主项为零。[前篇定理 1 及截尾证明](https://lixiang90.github.io/vibemath-NumberTheory/zeros-exact-distribution/articles/gap-two-quadratic-conductor-covariance.html) · [Dirichlet 零点计数](https://arxiv.org/abs/2005.02989)

对随 $k$ 移动的一个算术系数用
$a_{h_k}(v)=(\log v/\sqrt v)(\cos(h_k\log v)+O(1/v))$，其余两个固定系数保留在 (5) 中。二阶项的两个素数都在 $n$ 或 $n+2$，可以把它们的相位换成 $\log n$。利用 $2\cos x\cos y=\cos(x-y)+\cos(x+y)$，可得如下主项，两式各自的总误差同为 $O(Y^{-a}+K/(YL))+O_D(P^{-D})$，

$$
\begin{aligned}
V_k(\lambda)&=\sum_{\mathcal T}v_n\!
 \left[\cos\big(h_k(\lambda-1)\log n\big)+\cos\big(h_k(\lambda+1)\log n\big)\right]
 +\mathrm{err},\\
C_k(\lambda)&=\sum_{\mathcal B}\sum_{v\in\{p,s,t\}} b_{n,v}
 \left[\cos\big(h_k(\lambda\log n-\log v)\big)+
       \cos\big(h_k(\lambda\log n+\log v)\big)\right]
 +\mathrm{err}.
\end{aligned}
\tag{9}
$$

三阶式中 $v=p$ 时也可换 $\log p$ 为 $\log n$，代价已计入上述误差。求和使用 $\sum_{\mathcal R_Y}w(n)/L_Y\ll1$ 的上界筛。由于 $K\le Y^\sigma,\ \sigma<a$，所有误差的平方和为 $o(1)$。在第二式中留下 $v=p$ 的差相位，在第一式中留下两端共用的差相位，便是 (7) 的两个积分。下面处理其余相位。

## 4. 半素数因子的相位在细区间内有多少质量

把一个余弦写成 $\cos(2\pi k u)$。对于 $v=s,t$，其角度为 $u=(\lambda\log n\pm\log v)/(8L)$。在这些角度上放质量 $b_{n,v}$。同样，在两种统计的整素数端点和相位 $u=(\lambda+1)x_n/8$ 上分别放质量 $v_n$ 或 $b_{n,p}$。每种所得正测度 $\nu$ 的总质量为 $O(1)$。本节证明，对圆周上的任意区间 $I$ 均有

$$
\nu(I)\ll_\sigma |I|+Y^{-\sigma}.
\tag{10}
$$

特别是即使 $\lambda$ 落在允许的半素数因子比例之内，因子也无法在单个角度形成固定质量的原子。沿趋于无穷的 $Y$ 所取的任何弱极限都有有界密度。这里的幂次分辨率使 (7) 能同时覆盖多于 $\log Y$ 个高度。

**短区间筛。** 只须证明长 $\delta_0=Y^{-\sigma}$ 的区间所获质量为 $O_\sigma(\delta_0)$，再作覆盖。各角度在实线上有界，折回圆周只需考虑常数个逆像。令 $\Delta=\delta_0 L$。在一个物理二倍段 $U<n\le2U$ 中，将 $n$ 分到长度 $\asymp\Delta U$ 的小段，端部可向外补齐。一个小段内 $\log n$ 的变化为 $O(\Delta)$。固定目标角度区间后，$\log v$ 被限制在长 $O(\Delta)$ 的区间。因此 $v$ 落在某个 $V_0\asymp v$ 附近，候选区间长 $O(\Delta V_0)$。只要该区间有解，就有 $U^a\ll V_0\ll U^{1-a}$。

区间内素数 $v$ 的个数 $\ll_\sigma\Delta V_0/\log U$。对固定的 $v$，另一个因子 $t'$ 可在长 $O(\Delta U/V_0)$ 的区间中取值，且 $t',vt'\pm2$ 都须为素数。这给 $\ll_\sigma \Delta U/(V_0\log^2 U)$ 个选择，两个方向只多常数倍。两个上界都是**平移短区间**的上界筛。两个筛的长度均可向外取到至少常数乘 $U^{a-\sigma}$。在长度 $W$ 的区间内取大筛水平 $W^{1/3}$，只用奇素数 $\ell\le W^{1/24}$。二维时每个这样的素数恰禁两个类，因为 $v$ 和 $t'$ 均远大于筛限，$\ell$ 不整除 $2v$。分母可如下估计。积 $\prod_{\ell}(1-2/\ell)^{-1}\gg\log^2 W$，将展开项归一后，$\log d$ 的期望是 $\sum_\ell 2\log\ell/\ell=(1/12+o(1))\log W$，故保留 $d\le W^{1/3}$ 仍有固定比例。这给二维分母 $\gg\log^2 W$。一维时每个奇素数禁一个类，分母为 $\gg\log W$。这些估计的一切常数可统一选取。[Montgomery–Vaughan 的算术大筛，推论 1](https://personal.science.psu.edu/rcv4/personal/Publications/large_sieve.pdf#page=2)

两界相乘，再对 $O(1/\Delta)$ 个小段求和，整个二倍段的行数 $\ll_\sigma\Delta U/\log^3 U$。每行质量 $\ll\log^2 U/(UL)$。在 $I_Y$ 的 $O(L)$ 个二倍段上相加，因 $\log U\asymp L$，得到 $O_\sigma(\delta_0)$。这对 $s,t$ 两个因子均适用，无须假设它们独立分布。

整素数端的和相位更简单。$1\le\lambda+1\le5/2$，所以长 $\delta_0$ 的角度区间只对应一个对数宽度 $O(\Delta)$ 的 $n$ 段。对其中的 $n,n+2$ 用上面的二维筛，其小素因子已经由 $\mathcal R_Y$ 排除。这段至多横跨常数个相邻物理二倍段，所获质量也是 $O(\delta_0)$。同一论证直接用于 $x_n$，给 $\mu_2,\mu_3$ 的区间质量上界 $O(|I|+\delta_0)$。以上界都适用于真实四符号条件的任意保留子集，因为相应权非负。

**从区间质量到总平方偏差。** 若正测度 $\nu$ 满足 (10)，则

$$
\sum_{|k|\le K}|\widehat\nu(k)|^2\ll (1+K\delta_0)\nu(\mathbb T).
\tag{11}
$$

用阶数 $2K+1$ 的 Fejér 核，其 Fourier 系数在 $|k|\le K$ 时至少为 $1/2$。该核在圆周上至多为 $CK/(1+K\|u\|)^2$。把与 $\nu$ 的卷积分为中心的 $1/K$ 区间及左右二倍环带，代入 (10) 求和，得点态上界 $O(1+K\delta_0)$。再积分一次即得 (11)。对 (9) 中除整端差相位外的全部测度应用 (11)，得到 (7)。当 $\lambda\ne1$，把 $\mu_j$ 按 $x\mapsto(\lambda-1)x/8$ 推前，区间上界变为 $O(|I|/|\lambda-1|+\delta_0)$。同样的核证明与平凡的 $O(K)$ 界合用，就得 (8)。

## 5. 原好窗口上的正下界与失配时的响应

设所在窗口满足第 1 节的 $A_Y\ge c$。由 (6)、(7)，对每个固定 $\eta>0$、每个充分大的 $Y$，除至多 $O_{\sigma,\eta}(1)$ 个 $k$ 外，都有

$$
V_k(1)+\frac{121}{48}C_k(1)\ge c-\eta .
\tag{12}
$$

常数还允许依赖固定的原证书参数。这个正下界同时发生在同一幂次窗口的许多组实际零点上。若其中真孪生质量趋零，则除有界多个例外高度外，三阶项自己就须满足 $C_k(1)\ge (48/121)c-\eta$。相反，对固定 $\lambda\ne1$，即使真孪生质量未知，(8) 也保证 $|C_k(\lambda)|\ge\eta$ 的高度只有有界多个。

更定量地，若有 $m\ge1$ 个不同的 $k\le K$ 满足 $C_k(1)\le t$，则由 (7) 的平方和界和 (6) 得

$$
J_Y\ge A_Y-\frac{121}{48}
 \left(t+\frac{C_\sigma}{\sqrt m}\right)
 -O\big((RL)^2+Y^{-a}\big).
\tag{13}
$$

可以在看完零点和几何证书后选这 $m$ 个高度。若在无穷多个上述原窗口上有 $m\to\infty$ 且 $t\le(48/121)(c-\delta)$，其中 $\delta>0$ 固定，则 $\limsup J_Y\ge\delta$，进而 $\Gamma_1\ge\delta/3$。式 (8) 实际证明的众多小值出现在失配相位。所需的是式 (13) 中 $\lambda=1$ 的小值。目前没有得到它们在无穷多个窗口上出现的下界。[原来的 $\Gamma_1$ 与四符号孪生质量](https://github.com/lixiang90/vibemath-NumberTheory/blob/main/Goldbach_TwinPrimes/twin-prime-auxiliary-signs.md)

相位偏离多少才会影响观测也可计算。若 $K\to\infty,\ K\le Y^\sigma$ 且 $\lambda=1+\tau/K$，则对每个固定 $T$，在 $|\tau|\le T$ 上一致有

$$
\frac1K\sum_{k=1}^K Z_{j,k}(1+\tau/K)=
\int \operatorname{sinc}\!\left(\frac{2\pi\tau x}{8}\right)d\mu_j(x)
+O_T\left(\frac1K+Y^{-a}+\frac{K}{YL}\right),\qquad
\operatorname{sinc}(u)=\frac{\sin u}{u}.
\tag{14}
$$

在零处取连续值 $1$。当 $\lambda$ 如此接近 $1$ 时，背景的全部角度都与整数保持固定距离。对于 $s,t$ 的差相位，这来自 $\log v/\log n\le1-a+o(1)$。各和相位也有同样的间隙，且都小于 $1/2$。因此每个背景余弦在 $1\le k\le K$ 上的平均为 $O(1/K)$，由等比级数直接可见。整端差相位的平均是 $\int_0^1\cos(2\pi\tau x u/8)\,du$ 的 Riemann 和，误差 $O_T(1/K)$。代回 (9) 得证。在任何有物理质量弱极限的子列上，(14) 给出该极限的积分变换。随 $\tau$ 变化的曲线还确定物理极限测度，因为它在零点的偶次导数给出所有偶矩，而 $x$ 限在 $[1,3/2]$，$x\mapsto x^2$ 在这里单射。

### 从任选的一块连续高度检测稀疏孪生行

若高度是连续的一块，还可以得到相对误差。选 $m<k_0,\ k_0+m\le K$，记 $\overline V,\overline C$ 为 $V_{k_0+i}(1),C_{k_0+i}(1)$ 在 $-m\le i\le m$ 上以 $\binom{2m}{m+i}/4^m$ 取的加权平均。令 $\rho=\cos^2(\pi/40)<1$。一致有

$$
\begin{aligned}
|\overline V-J_Y|&\ll\left(2^{-m}+Y^{-1}+
        \left(\frac{K}{YL}\right)^2\right)J_Y+O_D(P^{-D}),\\
|\overline C-M_Y|&\ll\left(\rho^m+Y^{-a}+
        \left(\frac{K}{YL}\right)^2\right)(A_Y-J_Y)+O_D(P^{-D}).
\end{aligned}
\tag{15}
$$

在角度 $u$ 上，二项权的响应恰为 $\cos(2\pi k_0u)\cos^{2m}(\pi u)$。孪生行的和相位在 $[1/4,3/8+o(1)]$ 内，故至多剩 $2^{-m}$。半素数行除整端差相位以外的五个相位到整数的距离最终均至少为 $1/40$，因而贡献至多为原行权的 $O(\rho^m)$。整端差相位至多为 $O(1/(YL))$，用 $|1-\cos t|\le t^2/2$。二阶展开中略去的奇幂相对该孪生行为 $O(Y^{-1})$，三阶展开中略去的移动项奇幂相对该半素数行为 $O(Y^{-a})$。其他类型的行只留下条件累积量比较的任意负幂误差。这样在主要质量极小时也保留了 (15) 的精度。

例如取 $m\to\infty$，只要这批高度容得下相应的块，在一个给定充分大的窗口中，$\overline V>Y^{-2}$ 当且仅当四符号筛后至少有一对真实孪生素数。单对的质量 $\gg L/Y^{3/2}$ 已高于这个阈值。判定覆盖的只是给定窗口及见证，是否有无穷多个原好窗口通过，仍须对真实的零点统计作额外估计。

## 6. 数值复核

先穷举共享素因子 $5$ 的两行 $(53,55)$、$(83,85)$ 的两套素数符号，两行同时过筛的概率恰为 $1/256$。把第二行换成跨尺度例子中的 $(35,37)$，与 $(5,7)$ 同时过筛的概率为零。这复核了 (3) 为什么要核对乘积标签。实际抽样另取 $Y=128$ 的 45 个陈氏或孪生行，均匀遍历 $1027\le Q\le4096$ 中 $Q\equiv3\bmod8$ 的 101 个素数，再在 $\sqrt{60000000}<d_0\le60000000$ 中遍历 $d_0\equiv7\bmod16,\ 3Q\mid d_0+1$ 的素数。共有 10,070 对。先在各 $Q$ 下平均再平均，过筛行数的均值、方差为 $2.802192,2.579650$。独立平方类模型给 $2.812500,2.636719$。

**实际零点。** 对 $800<d\le4200$、$d\equiv23\bmod48$ 的全部 27 个素导子，以二次 theta 级数定位 1,416 个零点，逐导子用两种水平边网格的辐角计数核对。有限核参数取 $B=4.4,b=1.05,K_*=2$，三个窗心为 $0,7.01,14.13747391\ldots$，半宽均为 $3.45$，两两不交。用 $\lambda=1$ 比较 (2) 的单行贡献，暂省共同的窗归一化与见证筛选。下表的算术列在同一有限导子族上直接算条件矩，包含该族实际的符号相关。

| $n$ | 从零点算二阶 | 从素数和算二阶 | 从零点算三阶 | 从素数和算三阶 |
| ---: | ---: | ---: | ---: | ---: |
| 17 | 0.0410673 | 0.0410331 | 0.0849248 | 0.0848519 |
| 53 | −0.2283665 | −0.2284353 | −0.1412175 | −0.1410849 |

**检验整组高度。** 另取素数 $Q=2965867,d_0=12456641399$，其满足 $Q\equiv3\bmod8,\ d_0\equiv7\bmod16,\ 3\cdot5\cdot7Q\mid d_0+1$。在三个物理窗口穷举两端及其因子，再在独立符号模型中用排列和及最高九次的奇素数幂计算 96 个频率。取第二个固定算术高度为 $L^{-3}$。表中“高段”对 $49\le k\le96$ 平均，$M^{(0)}=\sum_{\mathcal B}2\theta_n(1-\theta_n)v_n$ 供比较。

| $Y$ | 筛后孪生 / 半素数行 | $J_Y$ | 高段 $V(1)$ | $M^{(0)}$ | 高段 $C(1)$ | 高段 $C(0.4)$ |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1024 | 8 / 31 | 0.0267754 | 0.0265976 | 0.0630466 | 0.0665107 | 0.0006741 |
| 8192 | 192 / 367 | 0.0642431 | 0.0641855 | 0.0535546 | 0.0543926 | 0.0041950 |
| 32768 | 1218 / 2266 | 0.0644709 | 0.0644854 | 0.0579875 | 0.0582874 | 0.0001658 |

在 $Y=32768$ 的 96 个频率上，三阶值低于 $M^{(0)}/2$ 的次数，在 $\lambda=1$ 时为 5，在 $\lambda=0.4$ 时为 86。取 $k_0=64,m=16$，(15) 的两个块平均为 $0.0644709,0.0592994$。在同一窗口取 $\tau=2$ 时，两种统计在 96 个频率上的均值为 $0.0289945,0.0254660$，以一次素数幂计算的 sinc 预测为 $0.0298164,0.0270795$。这组较大样本核对相位、累积量和因子算术，前一组才逐个求了零点。两组有限实验均未验证原上同调选窗条件，式 (12) 的共同好窗口由解析证明保证。[复算程序](https://github.com/lixiang90/vibemath-NumberTheory/blob/main/zeros-exact-distribution/gap-two-coherent-zero-heights/verify.py)
