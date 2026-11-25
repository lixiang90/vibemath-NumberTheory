---
title: "切换筛与固定差二三阶零点均值的显式上界"
date: "2026-11-25"
abstract: "在 Alpöge–Claude 的共同好窗口上，先固定两个素因子的乘积，再筛相差二的另一端。切换后的算术密度与三阶零点统计的因子权重相消，给出按因子比例分段的显式上界。结合 Fouvry–Radziwiłł 的分布定理，匹配三阶均值的极限上界小于 0.29997。文中还给出同一组分离零点窗口对各因子比例区间的读数，以及据此下估真孪生质量的公式。"
---

# 切换筛与固定差二三阶零点均值的显式上界

陈氏定理保证有很多素数 $p$ 使 $p+2$ 为素数或两个素数之积。对于后一种可能，令 $p+2=st$，$s<t$。第 35 篇先固定 $s$，对 $t$ 与 $st-2$ 同时用上界筛。这次让 $s,t$ 都遍历素数，只对 $st-2$ 用一次筛。这样利用了素数乘积在算术级数中的分布。另一方向 $st+2$ 也可以一起计算。[第 35 篇的因子比例与阈值](https://lixiang90.github.io/vibemath-NumberTheory/zeros-exact-distribution/articles/gap-two-chen-factor-ratios.html) · [陈氏切换方法的背景](https://arxiv.org/abs/2405.19063)

这一改变给出第 34 篇所需的匹配三阶均值的一个实际数值上界。还可以在同一批零点上单独读取任何一段因子比例，将半素数在不同区间的质量分别扣除。本文使用归档的 Dirichlet-GRH 做特征和估计，并把 Fouvry–Radziwiłł 的定理用于一个较窄的因子范围。最后的孪生素数严格下界仍有定量差距。[归档定理](https://github.com/leanprover-community/riemann-hypothesis#theorem-statements) · [第 34 篇的匹配统计](https://lixiang90.github.io/vibemath-NumberTheory/zeros-exact-distribution/articles/gap-two-coherent-zero-heights.html)

## 1. 所估计的零点均值

取原 $\zeta$ 构造的 $0<\epsilon<1/52,\ H=P^J,\ z=P^{\epsilon/5},\ z'=z(\log P)^2$，$J$ 固定。先均匀取满足 $e^{z/2}<Q\le e^z,\ Q\equiv3\bmod8$ 的素数，再均匀取 $e^{z'/2}<d_0\le e^{z'}$ 的素数，要求 $d_0\equiv7\bmod16$，模 $Q$ 和模各个不超过 $z$ 的奇素数都余 $-1$。几何步骤随后选择 $[Y,Y^2]\subset[P,H]$。$G_n$ 指示 $n,n+2$ 在 $Q$ 处均为二次剩余，在 $d_0$ 处均为二次非剩余。[原两阶段候选](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/MeanSquareBounds/PrimeSetTwoBounds.lean#L37) · [取对后选窗](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/CharacterBoxes/WinPairs.lean#L32)

写 $a=3/11,\ b=1/2,\ \mathcal I_Y=(2Y,Y^{3/2}/2]$，$L_Y=\sum_{Y\le k\le Y^{3/2}/2}1/k\sim(\log Y)/2$，分母对整数求和。保留 $n\in\mathcal I_Y$ 中 $N_n=n(n+2)$ 平方自由且全部素因子大于 $n^a$ 的行。记 $w(n)=\log n\log(n+2)/n,\ v_n=G_nw(n)/L_Y$。真孪生行的总质量为 $J_Y$。$\mathcal B$ 是一端为素数、另一端 $m=st$ 为不同素数之积的行，两种方向都计。设

$$
 \theta_n=\max\!\left(a,\frac{\log s}{\log m}\right),\quad f(u)=2u(1-u),\quad
 \nu_Y=\sum_{\mathcal B}v_n\delta_{\theta_n},\quad
 d\tau_Y=f\,d\nu_Y,\quad A_Y=J_Y+\nu_Y([a,b]).
 \tag{1}
$$

在 $a$ 处截断至多改动原比例 $O(1/(n\log n))$。上标 $\circ$ 表示去掉 $G_n$。对典型候选对，所有这些物理窗口均有 $A_Y\ge c>0$，并有固定正比例见证可以共用一段几何好窗口。这用到定量陈氏定理和原选择步骤中平方类的两阶段均方界。第 28 篇已有将因子比例收缩到固定内段的陈氏推论，本文精细化的是每一段的质量和它在零点中的读数。[陈氏定理的定量形式，定理 6.1](https://arxiv.org/pdf/math/0405581#page=19) · [原共同窗口与筛选比例](https://lixiang90.github.io/vibemath-NumberTheory/zeros-exact-distribution/articles/gap-two-chen-factor-ratios.html) · [第 28 篇第 5 节](https://lixiang90.github.io/vibemath-NumberTheory/zeros-exact-distribution/articles/gap-two-separated-zero-windows.html)

固定物理见证后，另取相同 $Q$ 下的观测导子 $d$。以 $\Delta_3(h_1,h_2,h_3;N)$ 表示三个零点和在 $(N/d)=+1,-1$ 两组中的三阶联合累积量之差的一半，正组在前。第 4 节给出完整定义。写 $\xi=3R$，$R$ 为零点窗半宽，$h_k=2\pi k/(8\log Y)$。匹配的三阶统计为

$$
 C_k=\frac1{L_Y}\sum_{n}\frac{2G_n\cos(h_k\log n)}{\log n}
          \Delta_3(0,\xi,h_k;N_n).
 \tag{2}
$$

和遍历全部保留行。对固定 $0<\sigma<a$，从 $1\le k\le Y^\sigma$ 任选 $m_0\to\infty$ 个不同高度作平均。第 34 篇的平方和界给 $\overline C=\tau_Y([a,b])+o(1)$，允许看完数据再挑选高度。

令 $C_2=\prod_{\ell>2}(1-(\ell-1)^{-2})$，乘积对素数取。前篇的密度界给 $\tau_Y(I)\le4C_2\int_I du/(1-u)+o(1)$。本篇证明，在两阶段概率趋于一的见证集上，对所有物理窗口及比例区间 $I\subset[a,b]$ 同时有

$$
 \tau_Y(I)\le 2C_2|I|+O\!\left(\frac{\log\log P}{\log P}\right),\qquad
 \nu_Y(I)\le C_2\int_I\frac{du}{u(1-u)}
                   +O\!\left(\frac{\log\log P}{\log P}\right).
 \tag{3}
$$

因此全段的匹配均值至多为 $5C_2/11+o(1)$。利用大小接近的两因子在更大模数上的分布还能证明

$$
 \limsup\overline C
 \le C_2\left(\frac5{11}-\frac1{56}+\log\frac{57}{56}\right)<0.29997 .
 \tag{4}
$$

这里取 $P\to\infty$，见证来自上述高概率集合，$Y$ 可依见证而定。几何上可取它们共同合格的窗口。式 (3) 两个上界都是对非负测度的控制，全段的 $\nu_Y$ 上界是 $C_2\log(8/3)+o(1)$。下文给出位置与因子比例的联合版本。

## 2. 交换筛选次序并保留常数

先对未加四符号的行计数。在物理小段 $U<n\le(1+h)U$ 上，先取 $h=1/\log U$。为后面容许更细的定位，也同时处理 $h=U^{-\kappa}$，其中 $0<\kappa<1/4$ 固定。对数短段以下记为 $\kappa=0$。误差常数可依赖 $\kappa,\epsilon,J$。写 $e_U=\log\log U/\log U$。

**将素数乘积分成矩形。** 对每个方向，限定 $m=st=n$ 或 $n+2$，要求另一端 $m\pm2$ 为素数。设 $\Delta=h/(\log U)^3$。按公比 $1+\Delta$ 将较小素因子 $s$ 分到互不相交的 $[S,(1+\Delta)S)$。取与目标比例相容的所有格，可向外补齐。每格的较大因子 $t$ 取在 $[(U-2)/((1+\Delta)S),((1+h)U+2)/S]$，得到原集合的矩形覆盖。矩形的个数是 $O(h^{-1}(\log U)^4)$。若目标比例为 $I$，所用的 $\log s/\log U$ 落在 $I$ 两端各扩 $O(1/\log U)$ 的范围内，$s\gg U^a,\ s\ll\sqrt U$。设 $W$ 为这些矩形内有序素数对的总数，尚未筛 $st\pm2$。则

$$
 W\le\frac{hU}{\log U}
       \left(\int_I\frac{du}{u(1-u)}+O_\kappa(1/\log U)\right).
 \tag{5}
$$

事实上，固定一个素数 $s$ 时，同格的 $t$ 区间长至多 $(h+O(\Delta+1/U))U/S$。对该区间用 GRH 素数定理，误差 $O(\sqrt{U/S}\log^2 U)$，与主计数尺度相比至多 $O(U^{\kappa-1/4}(\log U)^4)$，对 $\kappa=0$ 也趋零快于对数负幂。把 $S$ 换回 $s$，再对素数倒数作分部求和，便得到 $1/[u(1-u)]$。这里角度端点的误差也是乘在 $hU/\log U$ 上的。

**一个剩余类的误差。** 对奇平方自由 $q$，令 $A_q$ 计数上述有序对中 $st\equiv\mp2\bmod q$ 的对，$A_1=W$。在模 $q$ 的特征展开中，一个非主特征的矩形和是两个素数特征和之积。GRH 在每个因子上给平方根误差，故在 $q\le U$ 上

$$
 A_q=\frac{W}{\varphi(q)}+
 O\!\left(h^{-1}\sqrt U(\log U)^8+
      \frac{\omega(q)(hU^{1-a}+1)}{\varphi(q)}\right).
 \tag{6}
$$

第一项已经包含所有矩形，因为在每个矩形上 $\sqrt S\sqrt{U/S}=\sqrt U$。第二项来自主特征对 $(st,q)=1$ 的要求。固定 $s\mid q$ 时只涉及一个 $s$ 格，故损失 $O(hU/s+1)$ 对。固定 $t\mid q$ 时，所有这些矩形的乘积都在长 $O(hU)$ 的同一段内，$s$ 的格互不相交，故至多损失 $O(hU/t+1)$ 对。这使 (6) 对幂次短段仍成立。所用的素数特征和界由截断显式公式取与端点同高的截断并按单位段加总零点得到，支持特征模数大于短因子的情形。[讲义的零点计数与截断显式公式 Z2、Z3](https://www.math.uwaterloo.ca/~mrubinst/alpoge_claude_grh/zeta-paper_v2.pdf#page=67)

**筛相差二的端点。** 取 Selberg 上界筛半径

$$
 R_s=U^{1/4-\kappa}(\log U)^{-100},\qquad
 G(R_s)=\sum_{\substack{d\le R_s\\2\nmid d}}\mu^2(d)
                  \prod_{\ell\mid d}\frac1{\ell-2}
       =\frac{\log R_s}{2C_2}+O(1).
 \tag{7}
$$

局部密度 $g(\ell)=1/(\ell-1)$ 来自 (6) 的一个可逆剩余类。最优平方权 $\lambda_d$ 的 $\lambda_1=1,\ |\lambda_d|\le1$，且 $\sum_{d,e}\lambda_d\lambda_e g([d,e])=1/G(R_s)$，其中 $g$ 在平方自由数上按乘法延拓。为算 (7)，把系数的 Dirichlet 级数写成 $\zeta(1+v)F(v)$。在 $v=0$，$2$ 处因子为 $1/2$，各奇素数处为 $(1+1/(\ell-2))(1-1/\ell)$，故 $F(0)=1/(2C_2)$。余乘积在零点左侧的一个邻域绝对收敛，与调和和作卷积给 $O(1)$ 余项。模 $2$ 的条件已经由 $st\pm2$ 为奇数计入。[Selberg 平方权及其最小值](https://www.math.tau.ac.il/~rudnick/courses/sieves2015/SelbergSieve_lecture_2.pdf)

对矩形内每一对，将另一端为素数的指示用 $(\sum_{d\mid st\pm2}\lambda_d)^2$ 控制，因该端大于 $R_s$。展开后，(6) 的第一个误差最多乘 $R_s^2$，给 $O(hU/(\log U)^{20})$。对于第二个，给定 $q=[d,e]$ 最多有 $3^{\omega(q)}$ 种分配，且 $\sum_{q\le R_s^2}\mu^2(q)3^{\omega(q)}(1+\omega(q))/\varphi(q)\ll(\log U)^5$，从 Euler 乘积即得。第二项也可吸收。因此一个方向的计数至多为 $W/G(R_s)+O(hU/(\log U)^{20})$。

将两个方向相加，乘上 $w(n)\le (1+O(h+1/U))\log^2 U/U$，得到

$$
 \sum_{\substack{U<n\le(1+h)U\\n\in\mathcal B,\ \theta_n\in I}} w(n)
 \le h\left(\frac{16C_2}{1-4\kappa}\int_I\frac{du}{u(1-u)}
             +O_\kappa(e_U)\right).
 \tag{8}
$$

这里 $\mathcal B$ 在单段中指同样的行类型，不要求该段已包含于 $\mathcal I_Y$。$f$ 在 $[a,b]$ 上递增。将 $f$ 写为区间左端的值加上 $f'$ 的积分，对截短的子区间用 (8)，即可在左边插入 $f(\theta_n)$。右边的积分随之变为 $2|I|$。这就是三阶权重消去算术分母的地方。

## 3. 允许筛后选择位置和因子比例

从 $U_0=P$ 开始递推 $U_{j+1}=(1+h_j)U_j$ 至超过 $H+2$，$h_j$ 按上一节的约定取值。对接触 $\mathcal I_Y$ 的每格任选一个比例区间 $I_j\subset[a,b]$，也可选空集。以 $\tau_Y((I_j))$ 记只计入各格中 $\theta_n\in I_j$ 的 (1) 式三阶质量。对每个固定 $0\le\kappa<1/4$，以概率 $1-o(1)$，所有窗口及所有这类选择同时满足

$$
 \tau_Y((I_j))\le \frac{2C_2}{1-4\kappa}\,
           \frac{\sum_j h_j|I_j|}{L_Y}+O_\kappa(e_P).
 \tag{9}
$$

不同格的选择可以依赖已看到的数据。若每个格只追踪一个宽度趋零的因子比例区间，即使物理格长取 $U^{1-\kappa}$，它们合计至多携带 $o(1)$ 的三阶质量。

证明用两阶段均方界。在同一格内，每一行的三个非主平方类为 $n,n+2,n(n+2)$。不同保留行的端点互异，因为 $n\equiv5\bmod6$，且乘积超出本格端点范围。原算术筛的单行通过率为 $1/16+O_D(P^{-D})$，两行通过率为 $1/256+O_D(P^{-D})$。在指定小素数符号的稀疏族内同样成立，因为所有物理因子都大于指定范围。于是对一个格中任意预先给定的比例前缀，$\sum w(n)(G_n-1/16)$ 的二阶矩 $\ll h_j(\log U_j)^4/U_j+O_D(P^{-D})$。对 $O(\log P)$ 个等距角度网格的前缀用 Chebyshev 不等式，阈值取 $h_j/\log P$。再对物理格取并，只需累加 $\sum_j (h_jU_j)^{-1}$，这可用 $\int_P^{H+2}du/(u h(u))^2$ 及首格项控制。总失败概率为 $O_{\kappa,J,\epsilon}((\log P)^{10}P^{-1+2\kappa})$，任意负幂误差可预先取足。用 (8) 将角度端点向外移动至网格，每格产生 $O_\kappa(h_j e_P)$ 误差。对权 $wf$ 可直接使用同样的二阶矩，或用分部求和。最后 $\sum_jh_j=(1+O(e_P))L_Y$，零头格可补齐，得到 (9)。未筛选的对应上界由 (8) 对每个格直接求和，主系数为 (9) 的 $16$ 倍。[原均方与稀疏子族](https://lixiang90.github.io/vibemath-NumberTheory/zeros-exact-distribution/articles/gap-two-chen-factor-ratios.html)

取 $\kappa=0$ 且每格使用同一 $I$，即得 (3)，其中未乘 $f$ 的式子用 (8) 直接转移。还可令 $x_n=\log n/\log Y$。对于物理位置区间 $E\subset[1,3/2]$，只在接触 $E$ 的格中使用 $I$，其余格取空，便有

$$
 \sum_{\mathcal B:\,x_n\in E,\,\theta_n\in I}v_nf(\theta_n)
 \le 4C_2\,|E|\,|I|+O(e_P).
 \tag{10}
$$

因此这两个坐标的任意共同弱极限均有面积密度上界。仅有两个边际分布的密度界时，质量仍可能集中在一条曲线上。式 (9) 也控制了在许多短格中分别挑选的窄带。

## 4. 用相同的分离零点窗口读取各段质量

令 $B=z/8,\ L=\log Y,\ K_*=\lceil(\log P)^2\rceil,\ b_*=B/(4K_*),\ R=2/b_*$。取 $g=\mathbf1_{[-B,B]}*((2b_*)^{-1}\mathbf1_{[-b_*,b_*]})^{*2K_*}$，Fourier 变换用相位 $e^{-itu}$。给定第二阶段素导子 $d$，定义

$$
 S_d(h)=\sum_{|\gamma_d-h|\le R}\widehat g(\gamma_d-h),\qquad
 T_d(h)=-\tfrac12\big(S_d(h)-\log d-\mathbb E_d[S_d(h)-\log d]\big).
$$

这里 $1/2+i\gamma_d$ 为本原二次 $L(s,\chi_d)$ 的零点，按重数计。$\Delta_3$ 由三个 $T_d$ 分组后各自中心化，计算三变量乘积均值，再取两组之差的一半。物理筛始终使用 $d_0$。按第 35 篇也可以在预先给小素数指定符号的观测族中计算。[条件累积量与显式公式](https://lixiang90.github.io/vibemath-NumberTheory/zeros-exact-distribution/articles/gap-two-quadratic-conductor-covariance.html) · [Hughes–Rudnick 的 Dirichlet 显式公式](https://arxiv.org/pdf/math/0208230#page=4)

取偶数 $2\le K\le Y^{D_0}$，$D_0>0$ 固定。所有行共用窗心 $0,\xi=3R,\ h_0^*=6R,\ h_k^*=2\pi k/(8L)\ (1\le k\le K)$。任意两窗的间距最终都大于 $2R$。例如取趋于无穷的 $K=o(\log Y)$，则这些窗心全部趋向零。以下系数只用已知的整数 $n$、$Y$ 和待检验的比例区间，$\theta_n$ 无须输入。

对行 $n$ 置 $x=\log n/(8L)$。在 $\mathbb R/\mathbb Z$ 上取偶函数

$$
 B_{x,I}(t)=\tfrac12\sum_{J'=I,\,1-I}
       \big(\mathbf1_{xJ'}(t)+\mathbf1_{-xJ'}(t)\big),
 \qquad 1-I=\{1-u:u\in I\}.
$$

对空区间取零。将其与单位积分的 Jackson 核 $\mathcal J_k(t)$ 卷积，此核正比于 $(\sin(\pi kt)/\sin(\pi t))^4$，连续延拓至整数。用 $k=K/2+1$，写 $B_{x,I}*\mathcal J_k=\sum_{j=0}^K c_j(x,I)\cos(2\pi jt)$。设 $I(n)$ 为该行所在物理格选定的比例区间。定义

$$
 Z_Y((I_j))=\frac1{L_Y}\sum_{n}\frac{G_n}{\log n}
       \sum_{j=0}^K c_j(x,I(n))\,\Delta_3(0,\xi,h_j^*;N_n).
 \tag{11}
$$

此和只要求对 $n(n+2)$ 作第 1 节的粗筛、平方自由性检查及二次符号运算，输入中不含哪一端为素数、哪一端为半素数的标记。对每个原见证、所有 $Y$ 和所有逐格选择，一致有

$$
 |Z_Y((I_j))-\tau_Y((I_j))|\ll_{\kappa,D_0,\epsilon,J}
 K^{-1}+e_P+\big((R\log H)^2+Y^{-a}\big)\log(2K).
 \tag{12}
$$

**证明。** 独立二次符号模型中，三阶差在 $\omega(N_n)=3$ 时等于三个素因子的算术系数在三个高度上的排列和。若 $\omega(N_n)=2$，同时翻转全部自由符号消去三阶项。若 $\omega(N_n)>3$，次数三的项无法含全部标签。由原扭转的 GRH 特征和界比较固定阶矩，以上规则用于真实零点的误差为 $O_A(P^{-A})$。对 $|h_j^*|\le Y^{D_0}+1$，单位高度段零点数是 $O(\log(d(3+|h_j^*|+|t-h_j^*|)))$，同样的截尾界仍一致。算术系数中的一次素数幂是 $(\log v/\sqrt v)\cos(h_j^*\log v)$，其余奇次幂相对基本尺度为 $O(1/v)$。因此在三因子行上，$\Delta_3/\log n$ 除以 $w(n)$，主项为 $f(\theta_n)$ 乘以 $\sum_{v\in\{p,s,t\}}\cos(h_j^*\log v)$，$p$ 表示整素数端。逐项误差是 $O((R\log H)^2+Y^{-a})$，另有可任取负幂的绝对误差。对 $j=0$ 把余弦换成 $1$ 也在此误差内。有限个区间指示的 Fourier 系数为 $O(1/(1+j))$，故累积代价至多再乘 $\log(2K)$。[排列和及移动高度下的误差](https://lixiang90.github.io/vibemath-NumberTheory/zeros-exact-distribution/articles/gap-two-coherent-zero-heights.html) · [Dirichlet 零点计数](https://arxiv.org/abs/2005.02989)

在素因子 $v$ 处，上述余弦恰好读取平滑函数在 $\log v/(8L)$ 处的值。$p$ 的位置与 $x$ 相差 $O(1/(YL))$，到这四段的距离一致大于正常数。$s,t$ 分别在 $x\theta_n,\ x(1-\theta_n)$ 的 $O(1/(YL))$ 范围，且 $x$ 有一致的正上下界。因此两个因子落在选定区间时各给 $1/2$。区间相接处用平滑误差控制。具体地，$\mathcal J_k(t)\ll K(1+K\|t\|)^{-4}$，这是由正弦商及其四次方积分 $\asymp k^3$ 得到的。离边界距离 $r$ 处的尾积分 $\ll(1+Kr)^{-3}$。在每个物理格中，边界坐标可随格改变，而每个格最多有固定个数的边界。用 (8) 对这些边界的 $O(1/(YL))+O(2^m/K)$ 邻域逐层求和，误差是 $O_\kappa(h_j(K^{-1}+e_P))$，此处 $h_j$ 表示物理格的相对宽度。最后除以 $L_Y$ 求和即得 (12)。对所有选定区间同时成立的原因是 (8) 已对每格的所有区间一致。

这样，(3)、(9)、(10) 均可转成实际三阶零点统计的上界。后两式还允许按物理位置追踪不同的因子比例，零点数据只取一次。

## 5. 因子接近平方根时的额外节省

令 $\delta_0=1/112$。当 $u=\log s/\log(st)=1/2-\delta$ 且 $0<\delta<\delta_0$ 时，Fouvry–Radziwiłł 的定理 1.1 给固定剩余类中两个因子乘积超过平方根水平的分布估计。写

$$
 v(u)=
 \begin{cases}
  1,&a\le u\le b-\delta_0,\\
  1/(2(1-u)),&b-\delta_0<u\le b .
 \end{cases}
$$

在 (3) 的同一高概率事件上，对任意固定 $I\subset[a,b]$，有

$$
 \limsup_{P\to\infty}\tau_Y(I)\le2C_2\int_I v(u)\,du .
 \tag{13}
$$

(10) 的共同弱极限也由 $4C_2 v(u)\,dx\,du$ 控制。下面只在 $h=1/\log U$ 的格上使用额外输入，这样矩形总数是对数的固定次幂。

在一个矩形中，记长因子尺度为 $N$、短因子尺度为 $M$、$X=MN\asymp U$。系数为限制到各自子区间的素数指示。长侧满足定理要求的 Siegel–Walfisz 条件，对端点作差即可保留一致性。另加与整数 $r$ 互素的要求，因系数只取素数而至多删去 $\omega(r)$ 项。对小模使用 Siegel–Walfisz 或 GRH，对于大于任意指定高对数幂的模数可直接用区间的整数计数。于是参考定理适用于固定的两个剩余类 $+2,-2$。对于 $q\asymp Q'$，其不带权的绝对误差和为 $O_A(X/(\log X)^A)$，适用范围为

$$
 N^{56/23}X^{-17/23+\eta}\le Q'\le NX^{-\eta}
 \tag{14}
$$

其中 $\eta>0$ 固定且可任意小。若 $\log N/\log X=1/2+\delta+O(1/\log U)$，左端的幂指数为 $11/23+(56/23)\delta+\eta+O(1/\log U)$。在 $0<\delta<1/112$ 的任意紧子段中，可取统一余量使左端低于 $U^{1/2-\eta'}$。因而 (6) 对较小模数的估计能接上 (14)。[Fouvry–Radziwiłł，定理 1.1 及其参数范围](https://arxiv.org/pdf/1812.00562#page=3)

还要核对平方权引入的模数重数。对一个矩形令 $r_q$ 为其剩余类计数与主特征项之差。矩形内每个整数至多有两种有序的两素数分解，所以对所用模数 $|r_q|\ll U/\varphi(q)$。在奇平方自由模上，$\sum_{q\le U^{0.52}}9^{\omega(q)}|r_q|\ll U(\log U)^9$。由 Cauchy–Schwarz，(14) 给 $\sum 3^{\omega(q)}|r_q|\ll_B U/(\log U)^B$，其中 $B$ 可任取，只需将输入的 $A$ 取足。所有矩形和模数二倍段的损失均为对数幂。主特征与 $W/\varphi(q)$ 之间仍用 (6) 的第二项。

固定一个很短的比例段 $I_0$，其闭包在 $(b-\delta_0,b)$ 内。现在可取 $R_s^2=U^\lambda$，其中 $1/2<\lambda<1-\sup I_0$。于是 (7) 的筛主系数由 $8C_2/\log U$ 降为 $4C_2/(\lambda\log U)$。先固定紧子段及有限分割，令 $P\to\infty$，再让分割变细、$\lambda$ 趋近各段上限。紧子段两端未覆盖的质量由 (3) 控制，长度趋零时消失。于是得到 (13)。这一步给极限不等式，(3) 的有限尺度误差单独保留。积分 $2C_2v(u)$ 即得 (4)。在 $a,b$ 附近的点质量也由 (3) 排除。

**去掉辅助符号时的无条件结论。** 对 $h=1/\log U$，(8) 的算术上界以及 $\limsup\tau_Y^\circ(I)\le32C_2\int_I v(u)\,du$ 也可由经典分布定理无条件得到。此时两边素因子区间的相对宽度都是对数负幂，(5) 可用带任意对数节省的素数定理。对每个矩形，在 $q\le \sqrt U/(\log U)^{B_1}$ 以内改用 Bombieri–Friedlander–Iwaniec 定理 0(b) 的双线性大筛估计，误差绝对值之和为 $O_A(U/(\log U)^A)$。两因子尺度均在固定幂次范围，长侧素数指示在相对宽度 $\asymp1/\log U$ 的区间上满足其条件 $(\mathrm A_2)$。事实上，该指示的二范数乘 $\sqrt N$ 只比 $N$ 小一个对数因子，可在 Siegel–Walfisz 的任意对数节省中吸收。矩形个数和模数重数分别用 (5) 前的计数和上面的 Cauchy–Schwarz 估计控制。回到 (7) 时，把 $100$ 换为充分大的固定数使 $R_s^2$ 落在有效范围，首项及 $O(e_U)$ 不变。在接近平方根的内段再接 (14)，便得所述极限。这一论证只涉及去掉 $G_n$ 的整数计数。带四个符号的选择与零点表示仍使用第 1、3、4 节的输入。[Bombieri–Friedlander–Iwaniec，定理 0(b) 与式 (1.5)](https://doi.org/10.1007/BF02399204) · [Fouvry–Radziwiłł，定理 1.1 后的拼接说明](https://arxiv.org/pdf/1812.00562#page=3)

## 6. 可从总量和分段读数推出的孪生下界

给定三阶总质量，$\nu_Y([a,b])=\int d\tau_Y/f$。由于 $1/f$ 递减，在 (3) 的密度上界下，要产生尽可能多的半素数质量，最有利的分布是将三阶质量以密度 $2C_2$ 从左端开始填充。这个极值有初等闭式。对 $I=[l,r]$ 以及实数 $m$，记

$$
 t_I(m)=\min\left(r,l+\frac{\max(m,0)}{2C_2}\right),\qquad
 B_I(m)= C_2\log\frac{t_I(m)(1-l)}{l(1-t_I(m))}.
 \tag{15}
$$

取固定划分 $I_1,\ldots,I_k$ 覆盖 $[a,b]$，相邻格按半开规则处理。令 $Z_i$ 为 (11) 对全物理窗使用同一 $I_i$ 的读数，$E$ 为 (12) 右侧与 $e_P$ 之和。则

$$
 J_Y\ \ge\ \left[A_Y-\sum_{i=1}^k B_{I_i}(Z_i)-O_\kappa(kE)\right]_+ .
 \tag{16}
$$

证明可直接从前缀开始。在 $I=[l,r]$ 内，任意前缀 $[l,t]$ 的三阶质量至多是 $\min(\tau_Y(I),2C_2(t-l))+O(e_P)$。将 $1/f$ 作部分积分，得到将质量从 $l$ 往右填满所得的上界，正好是 (15) 的对数。该函数对 $m$ 的 Lipschitz 常数至多 $1/f(a)=121/48$，故能代入各段的实际零点读数。对有限原子，边界误差也在 $O(e_P)$ 内。式 (16) 允许只在重要的几段取更好的零点上界，其余段保留已证的筛上界。合并相邻两段只会放宽允许的质量分布，所给孪生下界相应变弱或相等。

即使只知道 (3)，也有 $J_Y\ge[A_Y-C_2\log(8/3)-O(e_P)]_+$。有了三阶总读数，若半素数质量是 $x\le C_2\log(8/3)$，所需匹配三阶均值的最小值为

$$
 F(x)=2C_2\left(\frac{a e^{x/C_2}}{1-a+a e^{x/C_2}}-a\right)
 \quad\text{至消失误差}.
 \tag{17}
$$

这是 (3) 允许的最小值。式 (13) 在靠近 $1/2$ 处还能收紧。因 $C_2/[u(1-u)]<11$ 在整个 $[a,b]$ 成立，(17) 对正的可容许 $x$ 严格高于第 35 篇采用常数 $11$ 的下阈值。若已知 $A_Y\ge c$，$0<c<C_2\log(8/3)$，且在无穷多个共同好窗口上 $\overline C\le F(c)-\delta$，$\delta>0$ 固定，则 $J_Y\ge2\delta-o(1)$。这是因为 $F'(x)=f(a e^{x/C_2}/(1-a+a e^{x/C_2}))\le1/2$。若 (16) 的分段盈余大于误差，哪怕盈余趋零，也能在该窗保证至少一对真孪生。

目前陈氏定理给出的 $c>0$ 尚不能与 (4) 合用推出这类盈余。可以直接检查已有不等式留下的自由度。例如在 $[0.32,0.40]$ 上给 $\tau$ 一个充分小的正常密度，取 $\nu=\tau/f,\ A=\nu([a,b]),\ J=0$，算术质量关系及上述密度上界同时成立。这是对不等式的检验。要证明 $\Gamma_1>0$，仍须在真实的共同好窗口上压低 (16) 的半素数总贡献，或取得足以超过它的 $A_Y$ 下界。[原 $\Gamma_1$ 与共同好窗口](https://github.com/lixiang90/vibemath-NumberTheory/blob/main/Goldbach_TwinPrimes/twin-prime-auxiliary-signs.md)

## 7. 数值核对

以有理数计算筛半径 $R_s=31$ 的最优平方权，所得二次型恰为 $5124735/18004813=1/G(31)$。在 $37\le s<65,\ 149\le t<239$ 的素数矩形中，$st-2$、$st+2$ 为素数的对数分别为 $19,13$，有限平方权和分别为 $34.001931,31.758294$。另用 100 比特 Arb 算到 $10^6$ 的 Euler 乘积，省略部分用 $\prod(1-u_i)\ge1-\sum u_i$ 及 $\sum_{n>N}(n-1)^{-2}\le1/(N-1)$ 夹住。这给 (4) 右端严格的区间 $(0.29996,0.29997)$。

取有限四符号见证 $(Q,d_0)=(2965867,12456641399)$，按端点 $a,0.35,0.43,0.50$ 分三段。下表用实际整数分解算出各段三阶质量，再代入 (15)。所列反推值比较的是同一数据在整段与分段时的数值，未扣渐近误差。

| $Y$ | 三段 $\tau_Y$ 的算术值 | 整段反推 | 三段反推 | 实际 $J_Y$ |
| ---: | --- | ---: | ---: | ---: |
| 8192 | 0.0180205，0.0218060，0.0137281 | 0.0501517 | 0.0594709 | 0.0642431 |
| 32768 | 0.0194090，0.0208172，0.0177613 | 0.0495965 | 0.0597486 | 0.0644709 |

独立用排列和、每个素数截至九次的奇幂及 Jackson 余弦系数计算共同高度读数。$Y=32768,\ K=2048$ 时，得到 $0.0193703,0.0209551,0.0178331$。这部分在算术一侧检验角度与归一化。另对 $800<d\le12000,\ d\equiv23\bmod48$ 的全部 $82$ 个素导子，用二次 theta 积分及辐角计数找到 $|\gamma|\le20.76748$ 内的 $4762$ 个零点。对 $n=17,53,83$ 各取三个比例段，用这些零点和同一有限族的逐素数显式公式分别计算，九个读数的最大差小于 $0.00000353$。该小导子实验用截尾、相互重叠的窗口来比较两侧浮点计算。渐近定理的窗口分离、筛常数和误差由前面各节的证明给出。[复算程序](https://github.com/lixiang90/vibemath-NumberTheory/blob/main/zeros-exact-distribution/gap-two-switched-cubic-bound/verify.py)
