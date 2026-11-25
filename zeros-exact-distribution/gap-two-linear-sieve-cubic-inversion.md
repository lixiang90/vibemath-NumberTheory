---
title: "可分解线性筛与固定差二零点计数的反演"
date: "2026-11-25"
abstract: "对原辅助导子族的固定差二三阶零点统计，用 Iwaniec 线性筛和 Pascadi 的双线性分布估计，把匹配均值的极限上界降到 0.24858 以下。最小因子比例区间的密度上界下降六分之一。另给出公共零点窗口的反演式。在有界高度内用对数个分离窗口校正奇数次素数幂，使二阶和三阶条件累积量分别以任意固定负幂精度计数孪生与半素数邻项。"
---

# 可分解线性筛与固定差二零点计数的反演

固定差二的陈氏候选有两类。两端均为素数时，二阶条件零点统计出现信号。一端为素数、另一端为两个素数之积时，三阶统计出现信号。后一类能占陈氏下界的多少，决定了这种方法对孪生素数能说什么。第 36 篇得到匹配三阶均值小于 $0.29997$ 的极限上界。本文把它降至 $0.24858$ 以下，并给出从同一批零点直接读取两类整数对数目的公式。[第 33 篇的条件矩](https://lixiang90.github.io/vibemath-NumberTheory/zeros-exact-distribution/articles/gap-two-quadratic-conductor-covariance.html) · [第 36 篇的切换筛](https://lixiang90.github.io/vibemath-NumberTheory/zeros-exact-distribution/articles/gap-two-switched-cubic-bound.html)

算术改进来自两个互相配合的选择。先遍历两个真素因子，再选能分解成卷积的筛权，以便在筛相邻端点时使用超过平方根的模数。零点反演则用到原候选素数上的真实 Jacobi 符号。它们允许在同一物理窗口比较许多不同的整数行。这里的零点属于一族二次 Dirichlet $L$ 函数。零点部分以原归档的 Dirichlet-GRH 和前篇的显式公式为输入。新的计数和密度不等式是本文的推导。[原定理陈述](https://github.com/leanprover-community/riemann-hypothesis#theorem-statements) · [原稿中的辅助素数与符号条件](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/MeanSquareBounds/PrimeSetTwoBounds.lean#L41)

## 1. 窗口与两个上界

固定 $0<\epsilon<1/52$ 和原几何步骤给出的 $J$。设 $H=P^J,\ z=P^{\epsilon/5},\ z'=z(\log P)^2$，取 $[Y,Y^2]\subset[P,H]$，记 $L=\log Y$。原 $\zeta$ 候选先在 $e^{z/2}<Q\le e^z$ 中取 $Q\equiv3\bmod8$ 的素数，再在 $e^{z'/2}<d_0\le e^{z'}$ 中取素数 $d_0\equiv7\bmod16$，且模 $Q$ 和各个不超过 $z$ 的奇素数都余 $-1$。记 $G_n$ 为 $n,n+2$ 在 $Q$ 处均为二次剩余、在 $d_0$ 处均为二次非剩余的指示。两阶段概率分别在当阶段的候选中均匀取。物理见证取定以后，零点观测遍历同一 $Q$ 下第二阶段的素导子 $d$。前篇的共窗结论允许这些算术估计在几何步骤选出的好窗口上使用。[原选窗量词](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/CharacterBoxes/WinPairs.lean#L5) · [共窗的均方转移](https://lixiang90.github.io/vibemath-NumberTheory/zeros-exact-distribution/articles/gap-two-chen-factor-ratios.html)

取 $a=3/11,\ b=1/2$。在 $\mathcal I_Y=(2Y,Y^{3/2}/2]$ 内保留 $N_n=n(n+2)$ 平方自由且所有素因子均大于 $n^a$ 的行，记为 $\mathcal R_Y$。写 $L_Y=\sum_{Y\le k\le Y^{3/2}/2}1/k\sim L/2$，$k$ 遍历整数，$w(n)=\log n\log(n+2)/n$。以 $G_nw(n)/L_Y$ 为行质量。$J_Y$ 为真孪生质量，$\mathcal B_Y$ 为一端是素数 $p$、另一端是 $m=st$ 的行，其中 $s<t$ 是素数，两种方向都计。置

$$
 u_n=\max(a,\log s/\log m),\quad f(u)=2u(1-u),\quad
 \nu_Y=\sum_{n\in\mathcal B_Y}\frac{G_nw(n)}{L_Y}\delta_{u_n},
 \quad d\tau_Y=f\,d\nu_Y,\quad A_Y=J_Y+\nu_Y([a,b]).
 \tag{1}
$$

第 34 篇的 $\overline C$ 是三个分离零点窗口的条件累积量乘回整数端点相位后，对不同高度的平均。在该篇条件下 $\overline C=\tau_Y([a,b])+o(1)$。$A_Y$ 有陈氏定理和辅助符号转移给出的固定正下界。写 $C_2=\prod_{\ell>2}(1-1/(\ell-1)^2)$，乘积遍历素数。以下函数中的指数是针对本题特定筛权的可用分布水平。

$$
 \lambda(u)=\max\!\left(\frac35,\ \min\!\left(1-u,\frac{2+3u}{5}\right)\right).
 \tag{2}
$$

**密度上界。** 在两阶段概率趋于一的见证集上，所有允许的物理窗口以及 $I\subset[a,b]$ 中的所有区间同时满足

$$
 \tau_Y(I)\le C_2\int_I\frac{du}{\lambda(u)}+o(1),\qquad
 \nu_Y(I)\le \frac{C_2}{2}\int_I\frac{du}{\lambda(u)u(1-u)}+o(1).
 \tag{3}
$$

这里 $P\to\infty$，包括在上述事件内依见证挑选 $Y,I$。设 $C_\tau,C_\nu$ 为 (3) 两个全段积分。按 $1/3,\ 3/8,\ 2/5$ 分段有

$$
 \limsup\overline C\le C_\tau
 =C_2\left(\frac{53}{198}+\frac83\log\frac{25}{24}\right)<0.24858,
 \qquad \limsup\nu_Y([a,b])\le C_\nu<0.53640 .
 \tag{4}
$$

第 36 篇对应的上界分别约为 $0.299970,\ 0.647506$，后者取自该篇式 (3)。本篇在 $[3/11,1/3]$ 上给三阶密度 $5C_2/3$，前篇在同一区间给 $2C_2$。两种统计的分段差别对反推孪生质量有影响。[匹配零点均值的定义与平方和估计](https://lixiang90.github.io/vibemath-NumberTheory/zeros-exact-distribution/articles/gap-two-coherent-zero-heights.html) · [定量陈氏定理，定理 6.1](https://arxiv.org/pdf/math/0405581#page=19)

## 2. 小因子大小决定的线性筛水平

先去掉 $G_n$。在 $X<n\le(1+h)X,\ h=1/\log X$ 上，把每个方向的 $s,t$ 放进素数矩形。$s$ 格的相对宽度取 $h/(\log X)^3$，$t$ 区间向外补齐使所有目标乘积均被覆盖。给定比例段 $I$，矩形个数 $O((\log X)^5)$，矩形内素数对总数 $W$ 满足

$$
 W\le\frac{hX}{\log X}
       \left(\int_I\frac{du}{u(1-u)}+O(1/\log X)\right).
 \tag{5}
$$

这就是第 36 篇式 (5)。两边的相对宽度均为对数负幂，经典素数定理也能给出 (5)。每对的相邻端点是 $st\pm2$。对奇素数 $\ell$，所筛的同余条件 $st\equiv\mp2\bmod\ell$ 占可逆剩余类的 $1/(\ell-1)$。设 $g(d)=1/\varphi(d)$ 对奇平方自由数作乘法延拓。采用 Iwaniec 的可分解上界线性筛，水平记为 $D=X^\lambda$，筛限取 $Z_s=D^{1/2}$。所有真正的素数端点均通过该筛。筛函数满足 $F(2)=e^\gamma$，而

$$
 V(Z_s)=\prod_{2<\ell<Z_s}\left(1-\frac1{\ell-1}\right)
 \sim\frac{2C_2e^{-\gamma}}{\log Z_s},\qquad
 V(Z_s)F(2)\sim\frac{4C_2}{\lambda\log X}.
 \tag{6}
$$

筛法给主项 $WV(Z_s)(F(2)+O(\eta_s))$，余项为对剩余类误差加权求和，$\eta_s>0$ 可任意小且先固定。权分解为仅依赖精度的有限个良好可分解序列。对任意 $D_1D_2=D$，每个序列可写成两条有界序列的卷积，分别支撑于 $d_i\le D_i$。限制模数的素因子为小于 $Z_s$ 的奇素数可乘到各卷积因子上，因为筛权已支撑在平方自由数上。这次所需的正是这种余项结构。[Lichtman，第 2.10 定理及式 (2-4)、(2-5)](https://msp.org/ant/2025/19-1/ant-v19-n1-p01-s.pdf#page=9)

在一个矩形内设 $N\asymp X^u$ 为较小素因子的尺度，$M\asymp X/N$。Pascadi 命题 4.4 对固定剩余类 $a_0=\pm2$ 估计一个模数为 $q_1q_2q_3$ 的双线性和。把该类的指示减去 $\mathbf1_{(st,q)=1}/\varphi(q)$，在定理的参数范围内，总误差为 $O_A(X(\log X)^{-A})$。三个模数因子都允许带有界系数。较小素因子一侧的素数指示满足 Siegel–Walfisz 条件，限制到子区间后对两端作差仍然成立。它也满足该命题的最小素因子条件。定理所需的三个带余量的不等式为

$$
 Q_1\le N X^{-\delta},\qquad
 N^2Q_2Q_3^2\le X^{1-15\delta},\qquad
 N^2Q_2^5Q_3^2\le X^{2-40\delta}.
 \tag{7}
$$

令 $q_3=1$，取 $Q_1$ 略小于 $N$，$Q_2=D/Q_1$。令 $t=\log Q_2/\log X$。令严格余量趋于零时，(7) 化为 $2u+t\le1,\ 2u+5t\le2$。这给出水平指数的开区间上端

$$
 \lambda_2(u)=u+\min(1-2u,(2-2u)/5)
              =\min(1-u,(2+3u)/5).
 \tag{8}
$$

具体地，在 $u<1/2$ 的闭比例段上先选 $\lambda<\inf\lambda_2$ 且大于段内的 $u$，再把 $Q_1$ 指数从 $u$ 减去足够小的正常数，在 $Q_2$ 中补回。分细比例段可以满足这一要求。与段端点以及模数二倍分割有关的固定常数都由余量吸收。在 $u=3/8$，(8) 达到 $5/8$。单用这个两因子分解，在全段积分已给三阶上界 $C_2[(5/3)\log(275/248)+\log(5/4)]\approx0.261015$，最后一小段用第 36 篇的密度界后取极限。第 3 节只需在 $\lambda_2<3/5$ 的范围继续改善它。[Pascadi，命题 4.4 与式 (4.15)](https://arxiv.org/pdf/2505.00653v2#page=20)

## 3. 两个素数的卷积在全段的基准水平

还需证明对 $a\le u\le b$ 可以统一使用每个严格小于 $3/5$ 的水平。这里两个变量都带素数指示。Pascadi 关于单个素数的定理 1.3(ii) 要在其证明中的双线性步骤上使用，才适用于本题。Lichtman 命题 5.5 给出了将此类论证推广至短区间素数乘积的做法。以下说明替换分布定理后的参数和矩形端点。[Lichtman，命题 5.5](https://msp.org/ant/2025/19-1/ant-v19-n1-p01-s.pdf#page=27) · [Pascadi，命题 5.2 及其应用](https://arxiv.org/pdf/2505.00653v2#page=24)

**分解两个素数权。** 用端点平滑的大区间覆盖每个矩形。两边各外扩相对宽度 $(\log X)^{-12}$，取等于一的光滑上包络，支撑在外扩后区间内。扩边相对两边原宽度至多 $O((\log X)^{-8})$，素数定理给 (5) 同样的误差。以 $\Lambda(v)/\log v$ 计每个因子，并加 $P^-(v)\ge z_0=X^{1/(\log\log X)^3}$，其中 $P^-$ 表示最小素因子。所有原素数保留。含非一次素数幂的额外对至多 $X^{1-a/2+o(1)}$，因为在 $v\asymp V$ 中高次幂仅 $O(\sqrt V\log V)$ 个。对其同余项求模数和，每对至多多付 $\tau(st\pm2)^{O(1)}=X^{o(1)}$，主特征项的代价为对数幂。

分别对两个 $\Lambda$ 使用三阶 Heath–Brown 恒等式。每项至多有十二个变量，每个带 Möbius 权的变量不超过其原端点尺度的立方根乘一个常数，故至多 $X^{(1-a)/3+o(1)}$。其余变量为常数权或对数权，所有变量均保留 $P^-\ge z_0$ 条件，约定 $P^-(1)=\infty$。把变量放入二倍区间。对两个原端点的光滑包络作 Mellin 反演，可将它们分配给各变量。包络在对数坐标的各阶导数为固定对数幂，分部积分允许以另一固定对数幂高度截断，误差小于任意指定对数负幂。两次频率积分的总损失也是对数幂。见下述截断估计对绝对系数的控制。

若某一变量尺度 $N_0$ 位于 $[X^{2\delta},X^{1/3+\delta}]$，把它放在短边，其余变量并成长边。Pascadi 命题 5.2 对 $D=X^{3/5-50\delta}$ 在上界线性筛支撑内分配模数的素因子，三组分别满足该文式 (5.3)。使用该文第 5 节的分箱筛权，把同一箱的素因子按此分配，即得到命题 4.4 所需的三组系数。这里按原证明先在水平上留出余量，再取分箱和小素数预筛参数。预筛的小因子可以并入一组，消耗任意小的固定幂余量。命题 4.4 的参数可比这里的 $\delta$ 小任意固定倍。因此任给 $\lambda<3/5$，都可先选足够小的 $\delta$ 并保留全部严格余量。一般的可分解序列在这一步需要上界线性筛的具体支撑。[筛权的分箱与三组分配](https://arxiv.org/html/2505.00653v2#S5) · [Maynard 的上界筛权构造](https://arxiv.org/pdf/2006.07088#page=23)

这里有两项技术核对。第一，短边的粗糙数常数权与粗糙数 Möbius 权均满足所需的 Siegel–Walfisz 估计。可分别由 $1$ 或 $\mu$ 卷积一个支撑于 $z_0$-光滑数、绝对值至多一的序列得到。截去卷积除数 $e>N_0^{1/2}$ 的倒数和有任意对数节省，Rankin 界给出因子 $\exp(-\log N_0/(2\log z_0))$，其余 Euler 乘积为对数幂。余下部分用小模数上的 Siegel–Walfisz 估计。另加与 $d$ 互素的要求时，可在相应 Euler 因子中一并去掉 $p\mid d$，代价至多 $\tau(d)^{O(1)}$。大模数用整数计数。分部求和容许对数权与前述对数高度的 Mellin 相位。

第二，合并而成的长边系数被固定次除数函数控制。可将超过 $(\log X)^C$ 的系数按绝对值删去，代价 $\ll_A X(\log X)^{-A}$。例如剩余类一侧，先对每个整数积 $v$ 的模数约数求和，用 $\tau(v-a_0)^{O(1)}$ 控制。Cauchy–Schwarz、固定次除数矩上界以及长边大系数的高阶矩界给出任意所需的对数节省。主特征一侧用 $\sum_{q\le D}\mu^2(q)\tau(q)^{O(1)}/\varphi(q)\ll(\log X)^{O(1)}$。这也控制了平滑端点的 Mellin 截尾。除去一个对数幂后可以使用命题 4.4 的有界系数形式。

若没有上述中等变量，每个 Möbius 变量都小于 $X^{2\delta}$。每个原素数端点的积至少为 $X^{a-o(1)}$，因而各自至少有一个尺度超过 $X^{1/3+\delta}$ 的变量，否则六个小变量不足以达到该端点。两个端点合计至多有两个这样的变量，故恰好每端一个，都是常数权或对数权。其余变量的积 $M_0\ll X^{20\delta}$。Maynard 的除数函数分布估计允许模数到 $(X/M_0)^{2/3-3\delta'}$，对充分小的 $\delta,\delta'$ 覆盖 $X^{3/5}$。该证明对两个大变量各自的任意子区间成立，分部求和吸收两个对数权和 Mellin 相位。其模数误差是绝对值和，可容纳筛法中有限个有界权；有固定次除数权时也可按上一段截断。这穷尽了所有分解项。[Maynard，命题 5.2 及对子区间的证明](https://arxiv.org/pdf/2006.07088#page=14)

**回到物理窗。** 上述分布估计的主项在 $(st,q)=1$ 的对上求和。要换成 $W/\varphi(q)$，固定 $s\mid q$ 或 $t\mid q$ 时，因乘积总在同一长 $O(hX)$ 的区间内，每个素除子最多损失 $O(hX^{1-a}+1)$ 对。平滑扩边带来的重叠有界。再对 $q$ 的固定次除数权除以 $\varphi(q)$ 求和，损失比 $hX/(\log X)^2$ 小任意固定对数幂。对所有矩形加总分布误差时，原定理中的对数节省可预先取足。

对一个比例小段用 $\lambda<\inf\lambda(u)$，将 (5)、(6) 相乘，再乘两个方向和 $\log^2 X/X$。这样得到未筛符号的半素数质量上界，每个物理小格的主项为 $8hC_2\int_I du/[\lambda u(1-u)]$。在积分内再乘 $f$ 即得三阶质量。原两阶段均方界将两者在高概率事件上乘以 $1/16$。具体地，第 36 篇第 3 节已对每个物理格中的 $O(\log P)$ 个角度前缀同时控制了符号波动。可用旧上界补齐端点，再用本篇估计替换这些前缀差的算术上界。先固定有限比例分割、水平余量和筛精度，令 $P\to\infty$，再让分割变细、余量和筛误差趋零。$\lambda$ 连续且有正下界，端点由旧密度界控制，从而 (3) 的 $o(1)$ 对区间一致。去掉 $G_n$ 后 (3) 两个主系数均乘 $16$，这部分是经典筛法与分布定理的无条件推论。

同理，对相对宽度 $h_j=1/\log X_j$ 的物理格各自选择区间 $I_j$，三阶上界为 $(C_2/L_Y)\sum_jh_j\int_{I_j}du/\lambda(u)+o(1)$。物理对数坐标 $\log n/L$ 与 $u$ 的共同弱极限密度因而至多 $2C_2/\lambda(u)$。这里的求和可用于看完物理窗后挑选的区间。对幂次短物理格，第 36 篇另给了有显式尺度损失的上界。

## 4. 小因子比例给出的孪生阈值

当已知的半素数质量为 $x$，在 (3) 下，能产生它的最小三阶质量来自先填满最左的比例。这由 $f$ 递增及对前缀分部求和得到。定义

$$
 D(t)=\frac{C_2}{2}\int_a^t\frac{du}{\lambda(u)u(1-u)},\qquad
 T(t)=C_2\int_a^t\frac{du}{\lambda(u)},\qquad
 F_*(x)=T(D^{-1}(x))\quad(0\le x\le C_\nu).
 \tag{9}
$$

特别地，$0\le x\le x_0=(5C_2/6)\log(4/3)$ 时，$D^{-1}(x)\le1/3$，所以

$$
 F_*(x)=\frac{5C_2}{3}
 \left(\frac{a e^{6x/(5C_2)}}{1-a+a e^{6x/(5C_2)}}-a\right).
 \tag{10}
$$

和第 36 篇式 (17) 相比，小 $x$ 展开的一次项仍为 $2a(1-a)x$，二次项系数从 $a(1-a)(1-2a)/C_2$ 提高为其 $6/5$ 倍。例如 $x=0.10$ 时，阈值从约 $0.041003$ 升至 $0.041260$。若 $A_Y\ge c$，且在无穷多个共同好窗口上 $\overline C\le F_*(c)-\delta$，$0<c<C_\nu,\ \delta>0$ 固定，则 $J_Y\ge 2\delta-o(1)$。这里 $F_*'(x)=f(D^{-1}(x))\le1/2$。对分段读数也可从每段左端开始填满，得到与第 36 篇式 (16) 同样形式的改进。

## 5. 同一批低处零点直接读取半素数总质量

记 $B=z/8,\ k_*=\lceil(\log P)^2\rceil,\ b_*=B/(4k_*),\ R=2/b_*$。令 $g=\mathbf1_{[-B,B]}*((2b_*)^{-1}\mathbf1_{[-b_*,b_*]})^{*2k_*}$，$\widehat g$ 用相位 $e^{-itu}$。对观测导子 $d$，用 $1/2+i\gamma_d$ 表示它的本原二次 $L$ 函数的零点。以导子族均值中心化，得到

$$
 S_d(h)=\sum_{|\gamma_d-h|\le R}\widehat g(\gamma_d-h),\quad
 T_d(h)=-\tfrac12\left(S_d(h)-\log d-\mathbb E_d[S_d(h)-\log d]\right).
$$

按 $(N_n/d)=+1,-1$ 分成两族。$\Delta_k$ 为相应联合累积量之差的一半，正组在前。二阶是条件协方差，三阶是三个变量在组内中心化后的乘积均值。由前篇的显式公式，标签 $v$ 是素数时，在线性零点观测 $T_d(h)$ 中的系数为

$$
 a_h(v)=\frac{\log v}{\sqrt v}
 \sum_{\substack{j\ge1\\j\ \mathrm{odd}}}
 v^{-(j-1)/2}g(j\log v)\cos(jh\log v).
 \tag{11}
$$

当 $r=\omega(N_n)=k$ 时，$\Delta_k$ 的主项是这些标签与观测变量之间的 $k!$ 个排列乘积之和。$r>k$ 时主项消失。$r=2,\ k=3$ 时同时翻转独立模型的所有符号也使主项消失。在归档 Dirichlet-GRH 上对原候选族估计全部所需扭转，便知这些规则在真实导子族上的绝对误差可小于 $P$ 的任意固定负幂，对有界高度一致。[条件矩的排列公式及原族误差](https://lixiang90.github.io/vibemath-NumberTheory/zeros-exact-distribution/articles/gap-two-quadratic-conductor-covariance.html) · [Dirichlet 显式公式](https://arxiv.org/pdf/math/0208230#page=4)

先给一个窗心可全部趋零的简单公式。写 $\omega=22\pi/39,\ \rho=\cos(2\pi/13)$，$\xi=3R,\ h_0=6R,\ h_j=j\omega/L\ (1\le j\le K)$。置 $D_j(n)=\Delta_3(T_d(0),T_d(\xi),T_d(h_j);N_n)$。对物理素因子 $v$ 都有 $a\le\log v/L\le3/2$。因此 $t=\omega\log v/L\in[2\pi/13,11\pi/13]$，$|\cos t|\le\rho<1$。用

$$
 G(y)=\frac1{\arccos y}=\frac{2/\pi}{1-(2/\pi)\arcsin y}
      =\sum_{j\ge0}g_j y^j,\qquad
 \sum_{j=0}^K g_j\cos^j t=\sum_{j=0}^K b_{j,K}\cos(jt).
$$

两组系数均非负，因 $\arcsin y$ 的 Taylor 系数非负且 $\cos^j t$ 的余弦展开非负。取 $r=K/(K+1)$，由 $\arccos r\ge\sqrt{2(1-r)}$，得 $\sum b_{j,K}=\sum_{j\le K}g_j\ll\sqrt{K+1}$。对 $|y|\le\rho$，Taylor 尾项绝对值至多 $O_\rho(\sqrt{K+1}\rho^K)$，其最大值恰在 $y=\rho$。设

$$
 B_{K,Y}=\frac1{L_Y}\sum_{n\in\mathcal R_Y}\frac{G_n}{\log n}
 \left(\frac{\omega\log n}{2L}\sum_{j=0}^K b_{j,K}D_j(n)-\frac{D_0(n)}6\right).
 \tag{12}
$$

对所有原物理见证，$1\le K\ll L$ 时都有

$$
 B_{K,Y}=\nu_Y([a,b])
   +O\!\left(\sqrt{K+1}\,[\rho^K+(R\log H)^2+Y^{-a}]\right).
 \tag{13}
$$

为看清 $1/6$ 的来源，对三因子行展开 (11)。除以 $\log n$ 和 $w(n)$ 后，$D_j$ 的主项为 $f(u_n)\sum_{v=p,s,t}\cos(h_j\log v)$。将 $h_0$ 的余弦换成 $1$、将第一对近零高度换为零，误差 $O((R\log H)^2)$。略去其余奇数次素数幂的误差 $O(Y^{-a})$。而在三个标签处，

$$
 \sum_{v=p,s,t}\left(\frac{\log n}{2\log v}-\frac16\right)
 =\frac1{2u_n(1-u_n)}+O(1/(n\log n)).
$$

这是因为 $\log p/\log n$ 和 $\log(st)/\log n$ 都是 $1+O(1/(n\log n))$。整素数标签给出的 $1/2$ 与三份 $1/6$ 相消。所有三因子行的未加符号归一化总质量有界，由 (5)、(6) 或前篇上界即可。乘回权重求和得到 (13)。取 $K$ 为充分大常数乘 $\log\log P$，误差便小于 $(\log P)^{-A}$，其中 $A>0$ 可预先指定。此时所有窗心都趋向零。用该式估计总质量无需剖分因子比例。

## 6. 校正奇数次素数幂，以整数精度计数

可把 (13) 中与窗口半宽和高次幂有关的误差一起消去。下面只用 $O_{\sigma}(\log Y)$ 个彼此分离、物理行共用的零点窗口，其最高高度受仅依精度的常数控制。字母 $\sigma>0$ 表示所需的幂次精度。输入含 $n,N_n$、粗筛、平方自由检查以及原 Jacobi 符号。观测权只依赖窗口与预选精度。

**归一化计数公式。** 可以构造三组实线性零点观测 $U_{i,d},\ i=1,2,3$，使所有 $n\in\mathcal R_Y$ 同时满足

$$
 E_{k,Y}(n)=\frac{\sqrt{N_n}}{k! L^k}
        \Delta_k(U_{1,d},\ldots,U_{k,d};N_n)
        =\mathbf1_{\omega(N_n)=k}+O_\sigma(Y^{-\sigma}),\qquad k=2,3 .
 \tag{14}
$$

$U_{i,d}$ 对该窗的所有行相同，按各自的 $N_n$ 分组。每组用的零点窗口与其他组分离。此前单行的正主项已经能识别因子数。这里使每个合格行的主项都等于一，因此可以先对任意一批行求和，再作一次取整。

**奇数次幂的反演。** 取固定正奇数 $m$ 使 $a(m+1)/2>\sigma+3$。令 $b_m=3m/2,\ \omega_m=\pi/(a+b_m),\ \rho_m=\cos(\omega_m a)$，置

$$
 V_L(x)=\frac1x\sum_{\substack{j\le m\\j\ \mathrm{odd}}}
              \frac{\mu(j)}{j}e^{-(j-1)Lx/2}.
 \tag{15}
$$

将 $V_L(dx)$ 乘以 $e^{-(d-1)Lx/2}$，在奇数 $d\le m$ 上求和。对 $k=dj\le m$，$e^{-(k-1)Lx/2}/(kx)$ 的系数为 $\sum_{j\mid k}\mu(j)$。于是 $a\le x\le3/2$ 时

$$
 \sum_{\substack{d\le m\\d\ \mathrm{odd}}}e^{-(d-1)Lx/2}V_L(dx)
 =\frac1x+O_m(Y^{-a(m+1)/2}).
 \tag{16}
$$

**公共高度和系数的大小。** 取位移 $\ell_i=0,4R,8R$。对 $i=1,2,3$，在 $y=0$ 展开

$$
 H_i(y)=\frac{V_L(\arccos y/\omega_m)}
               {\cos(\ell_i L\,\arccos y/\omega_m)} .
 \tag{17}
$$

用主支。在 $|y|\le r<1$ 上写 $\arccos y=t+iv$，则 $0<t<\pi$ 且 $\cos^2t+\sinh^2v=|y|^2$。于是 $t\ge\arccos r$，同时 $|v|\le\operatorname{arsinh}(1)$。因 $R L\to0$，(17) 的分母在单位圆盘内一致远离零。指数项模长至多一，故 $|H_i(y)|\ll_m(1-r)^{-1/2}$。Cauchy 估计取 $r=j/(j+1)$，给第 $j$ 个系数 $O_m(\sqrt{j+1})$。记 $P_{i,K}$ 为前 $K$ 次 Taylor 多项式。写

$$
 P_{i,K}(\cos t)=\sum_{k=0}^K c_{i,k}\cos(kt),\qquad
 \sum_{k=0}^K|c_{i,k}|\ll_m(K+1)^{3/2}.
$$

最后一个估计利用每个单项式余弦展开系数的绝对值和为一。$|y|\le\rho_m$ 上的截断误差 $\ll_m\sqrt{K+1}\rho_m^K$。取 $K=\lceil(\sigma+3)L/|\log\rho_m|\rceil$。按

$$
 U_{i,d}=\sum_{k=0}^K c_{i,k}\,
 \frac{T_d(k\omega_m/L+\ell_i)+T_d(k\omega_m/L-\ell_i)}2
 \tag{18}
$$

取零点。由于 $T_d$ 是偶函数，每个高度按绝对值使用。重合中心合并以后，三个集合共有 $5K+3$ 个中心。同一基本高度附近最小间距 $4R$，相邻基本高度的窗心相距至少 $\omega_m/L-16R>2R$ 对充分大的 $P$ 成立。最大中心不超过 $(\sigma+3)\omega_m/|\log\rho_m|+o(1)$。

**验证每个标签的系数。** 在 (11) 中对两个位移取平均，得到 $\cos(\ell_i j\log v)P_{i,K}(\cos(\omega_m jx))$，其中 $x=\log v/L$。若 $j\le m$，则 $jx\in[a,b_m]$，$|\cos(\omega_m jx)|\le\rho_m$，且 $g(j\log v)=1$，因为 $m\log v\ll_m L\ll B$。由 (17) 乘回位移余弦，就得到 $V_L(jx)$ 加上上述截断误差。若 $j>m$，余弦多项式在整个实轴上至多 $O_m((K+1)^{3/2})$，几何尾至多再乘 $O(Y^{-a(m+1)/2})$。由 (16)，标签 $v$ 在每个 $U_{i,d}$ 中的系数都为

$$
 \frac{\log v}{\sqrt v}\big(1/x+O_\sigma(Y^{-\sigma-1})\big)
 =\frac{L}{\sqrt v}\big(1+O_\sigma(Y^{-\sigma-1})\big).
$$

最后在原条件矩的排列公式中替换系数。$r=k$ 时每个排列的主项为 $L^k/\sqrt{N_n}$，其余行由第 5 节的消失规则处理。累积量具有多线性，有限线性组合至多使原族比较的误差再乘 $O_m((K+1)^{3k/2})$。$Y\le P^J$ 允许在原公式中预取足够大的负幂，乘 $\sqrt{N_n}/L^k$ 后仍小于 $Y^{-\sigma}$。这就证明了 (14)。

**一次取整。** 取 $\sigma>3/2$。对 $\mathcal R_Y$ 的任意子集 $\mathcal S$，包括按窗口数据再选的集合，在所有充分大的 $Y$ 上有

$$
 \#\{n\in\mathcal S:\omega(N_n)=k\}
 =\left\lfloor\ \sum_{n\in\mathcal S}E_{k,Y}(n)+\frac12\ \right\rfloor ,
 \qquad k=2,3 .
 \tag{19}
$$

因为行数至多 $Y^{3/2}$，总误差 $O_\sigma(Y^{3/2-\sigma})=o(1)$。保留行两端互素，各大于一，故两种因子数恰好分别计数孪生行与半素数邻项行。取 $\mathcal S=\mathcal R_Y$ 时，所有足够大的真孪生都通过粗筛。取 $\mathcal S=\{n\in\mathcal R_Y:G_n=1\}$ 则保留原四符号。对 $G_nw(n)/L_Y$ 求和也可用 (14)，即使只采用 $\sum_{n\in\mathcal I_Y}w(n)/L_Y\ll L^2$ 的粗界，总误差亦可取任意固定负幂。

## 7. 数值与剩下的正性问题

使用 Arb 计算 $10^6$ 以内的 $C_2$ 部分乘积。尾部由 $1-\sum_{n>10^6}(n-1)^{-2}\ge1-1/(10^6-1)$ 与 $1$ 夹住。结合 (2) 的四段初等积分，得到 $0.24857<C_\tau<0.24858,\ 0.53639<C_\nu<0.53640$。两个区间来自带误差包络的计算。式 (12) 的简单核，在 $K=32,64,96$ 的最大绝对逼近误差分别约为 $0.0099543,\ 0.00015076,\ 0.000002557$。

取有限见证 $(Q,d_0)=(2965867,12456641399)$。先用整数分解列出物理行，再用算术侧的排列和、(15) 至 (18) 的系数计算。下表取 $m=5,\ K=1100,\ R=0.0001$，三组窗口在该高度网格上也分离。算术实验使用平台值 $g=1$，奇数次幂计至 25 次。整数目标由直接计数给出。这个实验检验反演系数和先求和再取整的归一化。

| $Y$ | 两端均为素数的实际对数 | 二阶求和 | 半素数邻项实际对数 | 三阶求和 |
| ---: | ---: | ---: | ---: | ---: |
| 8192 | 192 | 192.0000000 | 367 | 367.0004334 |
| 32768 | 1218 | 1218.0000000 | 2266 | 2266.0003086 |

另取 $800<d\le12000,\ d\equiv23\bmod48$ 的全部 $82$ 个素导子，以二次 theta 积分和辐角计数求出 $|\gamma|\le20.76748$ 内的 $4762$ 个数值零点。取 $Y=32$ 和 $n=65,71,77,83,89$，对每个整数直接用 $(n(n+2)/d)$ 分组。式 (12) 的十五个有限族读数与同一族逐素数显式公式的最大差约为 $0.000039835$。校正奇数次幂的二阶、三阶共二十个读数，最大差约 $0.000348632$。小导子实验以相互重叠的宽截尾窗口比较零点侧与算术侧。渐近的窗口分离和精度来自第 6 节。[完整复算程序](https://github.com/lixiang90/vibemath-NumberTheory/blob/main/zeros-exact-distribution/gap-two-linear-sieve-cubic-inversion/verify.py)

对于无穷性，(19) 把单个窗口含有孪生素数写成归一化协方差和越过 $1/2$。对去掉物理筛 $G_n$ 的所有大尺度窗口，任一充分大的孪生对都落在某个这样的窗口内，例如取 $P=Y=p/3,\ J\ge2$，物理窗便含有较小端点 $p$。观测导子和分组符号照常使用。因而在无界的一列这些窗口上越过 $1/2$ 与孪生素数无穷性等价。若坚持原共同好窗和全部四个物理符号，在无界的一列这些窗口上越过 $1/2$ 会是一个更强的充分条件。

目前缺少的仍是这些读数在无穷多个窗口上的正下界。即使取了本篇更小的 $C_\nu$，所有已证密度上界依然容许一个小的陈氏下界完全来自半素数。例如令 $\nu$ 在 $[a,t]$ 上取 (3) 的允许密度，取 $A=\nu([a,t]),\ J=0$，随 $t\downarrow a$ 能得到任意小的正 $A$。在真实窗口中，只要能使 $A_Y-B_{K,Y}$ 超过 (13) 的误差就能检出真孪生。若这份盈余在无穷多个共同好窗口上有固定正下界，便推出原归一化的 $\Gamma_1>0$。式 (3) 改善了这一比较的算术上界，式 (14) 则给出了当盈余很小时仍可用的测量精度。[原 $\Gamma_1$ 与辅助符号](https://github.com/lixiang90/vibemath-NumberTheory/blob/main/Goldbach_TwinPrimes/twin-prime-auxiliary-signs.md)
