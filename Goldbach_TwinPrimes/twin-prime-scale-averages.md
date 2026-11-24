---
title: "固定差 2 的基准尺度平均与共同好窗口的正质量条件"
date: "2026-11-24"
abstract: "精确计算基准平均中每对孪生素数的权重。若其双对数上平均为正，有限 Szemerédi 定理允许一次选定网格，在共同好窗口上传递正质量；另给出等价的移位 Dirichlet 级数判据。"
---

# 固定差 2 的基准尺度平均与共同好窗口的正质量条件

在基准尺度上取测度 $dP/(P\log P)$，一对真实孪生素数 $p,p+2$ 给长区间平均主项贡献 $\log(p+2)/p$。本文证明换序误差对区间位置与长度一致有界，并给出一条能让有质量的尺度通过例外扭转指标筛选的条件定理。

具体地，若下文 (1) 定义的双对数上平均 $\mathcal D>0$，则对每个 $0<c<\mathcal D$，存在一组固定的选窗常数和无穷多个共同好窗口，其中四符号筛后的质量至少为 $c/16-o(1)$；相应的斜坡质量上极限 $\Gamma_1\ge c/48$。这里 $\mathcal D>0$ 本身是一个尚未证明的固定差 2 下界。第 2 节区分预先给定网格和容许一次调整常数这两种量词。[共同好窗口、质量及斜坡统计的定义](twin-prime-auxiliary-signs.md)

下文“共同好窗口”取自 Alpöge–Claude 的 ζ 目标。输入是每个解析好对至多有 $R$ 个例外扭转指标，允许的粗估计指数可增大且 $R$ 固定，以及辅助对在共同好窗上的四符号转移。式 (2)、(3)、(6)、(8) 自身只涉及普通孪生素数。

## 一　每一对真素数在尺度平均中留下多少权

设 $a=3/2$，$\alpha=\log a$，$Y(t)=\exp(\exp t)$，$m(t)=M(Y(t))$。当 $P=Y(t)$ 时，$dt=dP/(P\log P)$，幂变换 $P\mapsto P^A$ 成为平移 $t\mapsto t+\log A$。设 $w(p)=\log p\log(p+2)/p$。$M(Y)$ 在 $2Y<p\le Y^a/2$ 的孪生素数对上累加 $w(p)=\log p\log(p+2)/p$，再除以 $L_Y=\sum_{Y\le k\le Y^a/2}1/k$，其中 $k$ 取整数。定义正测度

$$
 v_p=\log\log p,\qquad
 \mu=\sum_{\substack{p,\ p+2\ {\rm prime}}}
       \frac{\log(p+2)}{p}\,\delta_{v_p},\qquad
 \mathcal D=\limsup_{\substack{T\to\infty\\ L\to\infty}}
       \frac{\mu((T,T+L])}{L}. \tag{1}
$$

这里的联合上极限允许任选 $T,L$ 同时趋于无穷。原变量中的长度为 $\log\log z-\log\log x$，其中 $x=Y(T)$、$z=Y(T+L)$。

孪生素数上界筛给 $\sum_{u<p\le2u,\ p+2\ {\rm prime}}w(p)\ll1$，下文凡以 $p$ 求和均指素数。于是 $\mu([t,t+1])\ll1$。事实上，在 $Y(t)<p\le Y(t)^e$ 中，每个二进段贡献 $\ll1/\log u$，把这些段从左到右相加即得常数。把 $M$ 的两个端点改为 $Y,Y^a$ 只增减常数个二进段。再用 $L_Y=(a-1)\log Y-\log 2+O(1/Y)$ 以及 $\sum_{Y<p\le Y^a}w(p)\ll\log Y$，一致得到

$$
 m(t)=\frac1{a-1}\int_{(t,t+\alpha]} e^{v-t}\,d\mu(v)+O(e^{-t}),
 \qquad 0\le m(t)\le C\quad(t\text{ 充分大}). \tag{2}
$$

常数 $C$ 可取固定的上界筛常数。核的积分是 $\int_0^\alpha e^s\,ds/(a-1)=1$。对 $T\le t\le T+L$ 交换积分次序，内点 $T+\alpha\le v\le T+L$ 得到全部单位权。两条边界带长 $\alpha$，总质量有界。因此对 $L\ge1$、$T$ 充分大，

$$
 \boxed{\displaystyle
 \int_T^{T+L} m(t)\,dt
   =\sum_{\substack{Y(T)<p\le Y(T+L)\\p+2\ {\rm prime}}}
          \frac{\log(p+2)}p+O(1).} \tag{3}
$$

误差与 $T,L$ 无关。也可直接核对一对素数的原始归一化。若积分区间覆盖其全部可见尺度，令 $u=\log Y$，其精确贡献是

$$
 w(p)\int_{(\log p+\log2)/a}^{\log p-\log2}
                 \frac{du}{uL_{e^u}}
 =\frac{\log(p+2)}p\left(1+O\!\left(\frac1{\log p}\right)\right).
$$

这个展开用 $L_{e^u}=(a-1)u+O(1)$ 即得。由 (3)，$\mathcal D$ 等于无筛选窗口质量 $m$ 在长尺度区间上的上平均。换序后差仍为固定的 $2$，每一项仍要求 $p+2$ 是素数。[Lichtman，孪生素数的上界筛](https://arxiv.org/abs/2109.02851)

## 二　固定网格上例外指标的损失

记 $C_t=401$，$\eta$ 为上同调证明中小扭转可用粗上界的指数，$B=C_t/\eta$，$h=\log(2B)$。例外指标数至多为 $R$，基准常数为 $J=4e^{Rh}$。辅助素数构造的 $R+1$ 个候选窗口在双对数坐标的起点为 $t,t+h,\ldots,t+Rh$。令 $\nu_i(P)$ 为在原解析好对中选到第一个合格下标 $i$ 的条件概率。按下面权重选窗所得的筛前质量为

$$
 b_h(t)=\sum_{i=0}^R\nu_i(Y(t))m(t+ih),\qquad \sum_i\nu_i(Y(t))=1.
$$

在等距基准 $t_k=t_0+kh$ 上可作一次有限而精确的换序。取 $0\le k<N$，记 $\nu_{k,i}=\nu_i(Y(t_k))$，则

$$
 \sum_{k=0}^{N-1}b_h(t_k)=\sum_{j=0}^{N-1+R}c_jm(t_j),
 \quad c_j=\sum_{\substack{0\le i\le R\\0\le j-i<N}}\nu_{j-i,i},
 \quad 0\le c_j\le R+1,\quad \sum_j c_j=N. \tag{4}
$$

总权守恒，但某个 $c_j$ 可以为零。一个无需假设独立性的固定网格结论来自逐点不等式 $b_h(t)\ge\min_i m(t+ih)\ge\sum_i m(t+ih)-RC$。在实现 (1) 的长区间上积分这个可测下界，固定平移只损失 $O_{R,h,C}(1)$ 的边界积分，故

$$
 \limsup_{t\to\infty} b_h(t)\ge
 \max\{0,(R+1)\mathcal D-RC\}. \tag{5}
$$

原证明给的例外集依赖外端和已选的 $q,d$。其中二进塔的底域可以随辅助对改变。甚至只固定 $q$，允许的基准也须满足 $(\log q)^{5/\epsilon}\le P<(2\log q)^{5/\epsilon}$。这个范围在 $\log\log P$ 轴上的宽度趋于零。因此随 $P$ 积分时，不能把一份例外表用于所有基准。[例外集的实际选择顺序](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/AutomorphicVarieties/LayeredIsomorphisms/ThmBabOf.lean) · [辅助素数的原始大小条件](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/CharacterBoxes/ExistsPageN.lean)

**一个检验量词的例子。** 固定上述网格且 $R\ge1$。一个正整数扭转 $k>1$ 恰在 $Y(s)^\eta<k\le Y(s)^{2C_t}$ 时会破坏整段 $[Y(s),Y(s)^2]$ 所用的充分条件。因为此时可取 $x=\max(Y(s),k^{1/C_t})$，有 $x^\eta<k\le x^{C_t}$。这里原参数满足 $\eta<C_t$。对固定 $k$，这些起点组成长度为 $h$ 的半开区间

$$
 [\,\log\log k-\log(2C_t),\ \log\log k-\log\eta\,).
$$

在双对数轴上取一个周期为 $Q$ 的测试测度，每周期只在一段长度 $\theta$ 上有密度 $1$，其余为零。取 $0<\theta<h-\alpha$，$Q>\log(J/2)+\theta+\alpha$。它的长区间平均为 $\theta/Q>0$。用 (2) 的核平滑后，正质量只出现在每周期长度 $\theta+\alpha<h$ 的一段上。对任意基准 $t$，所有可能的完整原始窗口起点都在 $[t,t+\log(J/2)]$。这段至多碰到一段正质量区。若交集的闭包为 $[u,v]$，则 $v-u\le\theta+\alpha<h$。对充分晚的基准，可选整数

$$
 Y(v)^\eta<k\le Y(u)^{2C_t}.
$$

两端的对数之比至少为 $e^{h-\theta-\alpha}>1$，因而确实容得下整数。且 $k\le(P^J)^{C_t}$。这一项已能遮住该基准的整个正质量区。这个例子只使用了挑窗充分条件和例外基数，测试测度与指标都没有被赋予真实素数或数域的性质。它说明，固定 $h$ 时，仅从正的普通平均及“至多 $R$ 个”无法得到对任意大小正平均都有效的传递结论。$R=0$ 时只有一个窗口，式 (3) 已经给出所需传递。

## 三　允许一次性选取原证明的常数时，可以消除选窗损失

**条件命题。** 若 $\mathcal D>0$，则对每个固定的 $0<c<\mathcal D$，有一组原证明允许的固定常数和无界的一批基准 $P$，使辅助素数构造的全部 $R+1$ 个候选起点都满足 $M(Y_i)>c$。这批基准在 $dP/(P\log P)$ 下的上密度为正，这里的密度采用 (1) 的移动长区间。每个这样的基准都有一个正比例原辅助对共用的完整好窗口，因而选中窗口的筛前质量上极限至少为 $c$，且可选择辅助对使斜坡质量 $\Gamma_1\ge c/48$。

**先找同时有质量的尺度。** 令 $F_c=\{t:m(t)>c\}$。由 (2)、(3)，有固定 $\delta>0$，且有任意晚、长度趋于无穷的区间 $I$ 满足 $|I\cap F_c|\ge\delta|I|$。事实上可取平均大于 $(c+\mathcal D)/2$ 的区间，再用 $m\le C$。设 $R\ge1$，用有限 Szemerédi 定理取 $N=N(R+1,\delta/2)$，使 $\{0,\ldots,N-1\}$ 中密度至少为 $\delta/2$ 的任何集合都含长 $R+1$、公差为正的等差数列。[Szemerédi，正密度整数集中的等差数列](https://doi.org/10.4064/aa-27-1-199-245)

把 $I$ 切为长度 $Nh$ 的整段，末尾余段长小于 $Nh$。在每个整段左端 $u$，以 $0\le r<h$ 遍历 $u+r+jh$，$0\le j<N$。对足够长的 $I$，舍去余段后平均密度仍至少 $3\delta/4$。所以在所有“整段、相位 $r$”的组合中，至少 $\delta/4$ 的比例有 $\delta N/2$ 个点属于 $F_c$。每个这样的组合贡献至少一个数列，其公差为 $nh$，$1\le n\le N$。令

$$
 E_{n,c}=\bigcap_{i=0}^{R}(F_c-inh).
$$

对固定公差，每个数列起点在一个整段和一个相位下只对应一个格点。对相位积分，再对整段与公差求和，得 $\sum_{n=1}^N|I\cap E_{n,c}|\gg\delta |I|/N$。在一列越来越晚且越来越长的 $I$ 上，至少有**同一个** $n\le N$ 满足 $|I\cap E_{n,c}|/|I|\gg\delta/N^2$。从这里起把 $n$ 固定。用这些交集可以等价地表述普通尺度条件。

$$
 \mathcal D>0
 \quad\Longleftrightarrow\quad
 \exists n\ge1,\ c>0,\quad
 \limsup_{T,L\to\infty}\frac{|E_{n,c}\cap[T,T+L]|}{L}>0. \tag{6}
$$

反向只需 $E_{n,c}\subset F_c$ 及 (3)。$R=0$ 可直接取 $n=1$。

**再把这个公差放回原上同调估计。** 在归档的符号中，设 $C_K\ge1$ 是有限层粗上界的有效指数。原证明取

$$
 \eta=\epsilon/(2C_K),\qquad h=\log(4C_K C_t/\epsilon),
 \qquad J=4e^{Rh}.
$$

引理 kmono 表明，$C_K$ 增大后同一粗上界仍成立。在构造窗口常数时，例外个数来自固定的二进塔秩界及有限个奇素数塔秩界之和，这两组秩界在选粗指数以前已经给定。因此沿用同一 $R$，取

$$
 C_K'=C_K e^{(n-1)h},\qquad
 \eta'=\epsilon/(2C_K'),\qquad h'=nh,\qquad J'=4e^{Rnh}. \tag{7}
$$

这是一次性选择，$C_K',J'$ 随后对所有 $P,q,d$ 相同。原证明对任意固定外端指数选择辅助素数，取 $P^{J'}$ 为外端。新网格的起点就是 $Y(t+inh)$。其指标区间 $(Y_i^{\eta'},Y_i^{2C_t}]$ 两两不交，所有 $Y_i^2\le P^{J'}$。由原构造，对每个解析好对，至少一个这样的完整窗口有证书。各窗口的联合概率相加至少为解析好对概率 $1/4$，故某个窗口占全部候选辅助对的比例至少 $1/[4(R+1)]$。[粗上界的单调性 kmono](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/AutomorphicVarieties/LayeredIsomorphisms/Kmono.lean) · [秩界与粗指数的分别选取](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/AutomorphicVarieties/LayeredIsomorphisms/ThmBabOf.lean) · [将这些常数用于完整上同调窗口](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/AutomorphicVarieties/HeightAndVanishing/JointWinOfCounts.lean) · [固定外端指数的两阶段选素数](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/MeanSquareBounds/PairShare.lean)

对 $t\in E_{n,c}$，这 $R+1$ 个窗口的筛前质量都大于 $c$。因此即使用原来的“选第一个合格窗口”，各权 $\nu_i$ 加起来也给 $b_{h'}(t)>c$。在上述共同窗上用四符号转移的条件均方界，可选到筛后质量至少为 $c/16-o(1)$ 的辅助对。斜坡积分再给系数 $1/3$，证得所述 $\Gamma_1\ge c/48$。这里的结论属于 (7) 所给的固定网格。对先前已经指定的特定 $h$，本论证并未把 Szemerédi 定理给的整数 $n$ 强制成 $1$。[四符号转移式 (9) 及斜坡归一化式 (13)](twin-prime-auxiliary-signs.md)

反过来，对任何一组固定常数都有 $b_{h'}(t)\le\sum_{i=0}^R m(t+ih')$。故通过长基准区间上的**正平均**证正性时，(3) 也要求 $\mathcal D>0$。上述共同好窗的存在性只要求**正上极限**，稀疏的一列好窗口也可能达到它。

## 四　把缺少的下界写成一个固定移位 Dirichlet 级数估计

对实数 $s>0$ 定义

$$
 F(s)=\sum_{p,\ p+2\ {\rm prime}}\frac{\log(p+2)}{p^{1+s}},
 \qquad G(s)=\sum_{m\ge2}\frac{\Lambda(m)\Lambda(m+2)}{m^{1+s}} .
$$

式 (1) 的测度满足 $F(e^{-T})=\int e^{-e^{v-T}}\,d\mu(v)$。它与 $\mu((-\infty,T])$ 之差为 $O(1)$。在 $T$ 左边第 $j$ 个单位段，核与 $1$ 相差 $O(e^{-j})$，右边第 $j$ 个单位段的核至多为 $\exp(-e^j)$。由 $\mu$ 的单位段上界可求和。又有 $-F'(s)=G(s)+O(1)$，误差对 $0<s\le1$ 一致，因为含高次素数幂端点的贡献至多为常数倍的 $\sum_{\ell\ {\rm prime}}\sum_{k\ge2}(k\log\ell)^2/\ell^k<\infty$。两式结合给出

$$
 \boxed{\displaystyle
 \mathcal D=
 \limsup_{\substack{0<\tau<\sigma\to0\\
                     \log(\sigma/\tau)\to\infty}}
 \frac{1}{\log(\sigma/\tau)}
 \int_\tau^\sigma sG(s)\,\frac{ds}{s}.} \tag{8}
$$

这里只在绝对收敛的实半轴取值。上界筛给 $sG(s)\ll1$，系数非负给 $sG(s)\ge0$。第三节真正需要的是 (8) 的严格正性。若有 Hardy–Littlewood 的差 $2$ 渐近，会有 $sG(s)\to\mathfrak S(2)$。式 (8) 允许 $sG(s)$ 在不同尺度上波动，也只需要一列长区间。Gou–Wang–Xu 在整数情形由比值猜想推导的渐近比这里所需的输入更强。[Gou–Wang–Xu，固定移位素数对的条件渐近](https://arxiv.org/abs/2611.18825)

## 五　目前能在每一个共同好窗口上证明哪一种正性

取 Chen 素数集合 $\mathcal C$，含 $p+2$ 为素数的素数 $p$，以及 $p+2=r_1r_2$、$r_1,r_2$ 均为大于 $p^{3/11}$ 的素数时的 $p$。Green–Tao 第 6 节使用的 Chen 型定理给每个充分大的 $(u,2u]$ 至少 $\gg u/(\log u)^2$ 个这样的 $p$。每个完整二进段对 $\sum w(p)$ 的贡献遂 $\gg1$。在 $(2Y,Y^{3/2}/2]$ 内有 $\gg\log Y$ 个完整二进段，除以 $L_Y$ 得到对**所有**充分大的 $Y$ 成立的下界

$$
 C_{\mathcal C}(Y)=\frac1{L_Y}
       \sum_{\substack{2Y<p\le Y^{3/2}/2\\p\in\mathcal C}}w(p)
 \ \ge c_{\mathcal C}>0. \tag{9}
$$

特别地，原来的固定网格、任何第一个合格下标和任何 $P$ 平均都保留这个下界。[Green–Tao，定理 6.1 及 Chen 素数定义](https://jtnb.centre-mersenne.org/item/10.5802/jtnb.538.pdf#page=19)

令 $\lambda(n)=(-1)^{\Omega(n)}$，其中 $\Omega$ 按重数计素因子。以 $E_{\mathcal C}(Y)$ 记 (9) 中再乘 $\lambda(p+2)$ 的带符号和，使用相同的归一化。在这个确切的支撑上有

$$
 2M(Y)=C_{\mathcal C}(Y)-E_{\mathcal C}(Y). \tag{10}
$$

在真孪生素数上符号为 $-1$，在两个大素数的乘积上为 $+1$。若直接对 $P$ 取平均，可把尚需的算术输入写得更具体。令 $v_p=\log\log p$ 如前，在第 (1) 式的长区间 $I=(T,T+L]$ 上写

$$
 D_{\mathcal C}(I)=\sum_{\substack{p\in\mathcal C\\v_p\in I}}\frac{\log(p+2)}p,
 \qquad
 D_{\mathcal C,\lambda}(I)=\sum_{\substack{p\in\mathcal C\\v_p\in I}}
                  \lambda(p+2)\frac{\log(p+2)}p.
$$

对原变量中的每个完整二进段，Chen 型下界现在给 $\gg1/\log u$，求和得 $D_{\mathcal C}(I)\gg L$，对充分大的 $T,L$ 一致。于是只要在一列 $T,L\to\infty$ 上证得 $D_{\mathcal C,\lambda}(I)\le(1-\beta)D_{\mathcal C}(I)$，其中 $\beta>0$ 固定，恒等式 $2\mu(I)=D_{\mathcal C}(I)-D_{\mathcal C,\lambda}(I)$ 就给 $\mathcal D>0$。这只需一次单侧的奇偶估计。这个带符号和须保留 $\mathcal C$ 的因数限制。例外扭转指标的个数界留下了对这个受限和的估计。

## 六　有限核对

对实际孪生素数用两种方式计算 (2) 主项在 $Y$ 从 $50$ 到 $5000$ 时的积分。按每对素数的截断核积分，得 $0.966991191626858$。按全部端点把 $t$ 轴分段、逐段积分，得 $0.966991191626859$。直接在 $50<p\le5000$ 累计 $\log(p+2)/p$ 得 $0.936203734382581$。前两数核对了换序，第三数保留了 (3) 中的边界差。

下表枚举实际整数窗口，保留原来的 $2Y,Y^{3/2}/2,L_Y$。半素数列只计两个因子均超过 $p^{3/11}$ 的情形，以整数比较 $r^{11}>p^3$ 判断最小因子 $r$。

| $Y$ | 孪生素数对的数量 | 大因子半素数邻项数 | $M(Y)$ | $C_{\mathcal C}(Y)$ | $E_{\mathcal C}(Y)$ |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 100 | 9 | 14 | 0.584517 | 1.451081 | 0.282047 |
| 1000 | 221 | 208 | 0.981613 | 1.858684 | −0.104543 |
| 10000 | 4223 | 3910 | 1.075810 | 2.063120 | −0.088500 |

这些数值核对了 (10) 的符号和 $1/2$ 系数。所取三个 $Y$ 只用于算术核对。第三节是纸面条件推导，$\mathcal D>0$ 与 $\Gamma_1>0$ 的实际证明仍然开放。
