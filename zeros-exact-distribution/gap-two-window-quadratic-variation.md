---
title: "同一辅助素数对上的多区间零点波动与孪生素数"
date: "2026-11-25"
abstract: "在 Alpöge–Claude 的同一次辅助素数抽样中研究多个导子区间。由陈氏定理，在每个大二倍区间都能构造一个具有非退化正态极限的零点统计量。换用不同因子符号类，四个分离高度窗的零点乘积之和以真孪生素数权重的八分之一为中心值。其误差在任意预先指定的原窗口上可由同批零点标定并趋于正态。若用双对数半宽的窗口，所得符号判据已能检测趋零的加权质量。正号概率的严格超额仍待证明。"
---

# 同一辅助素数对上的多区间零点波动与孪生素数

前篇在一个二倍区间内翻转端点 $p$ 的二次符号，用 Dirichlet 零点和的二阶矩测出 $p+2=rs$ 的半素数质量。本篇固定同一对辅助素数，考察沿 $p$ 变化的多个区间。这样有两个实际问题。较早的素数端点可能成为较晚一行的因子，幂次窗口的起点又由辅助素数决定。若直接把不同区间当作独立样本，就会在孪生素数下界的严格正性处产生选择偏差。[前篇的单区间零点矩](https://lixiang90.github.io/vibemath-NumberTheory/zeros-exact-distribution/articles/gap-two-chen-auxiliary-switch.html)

将每个二倍区间等分成 $K$ 格，在每格对零点和取平方，再沿原物理窗口求和。本文得到三项结果。

- 对任意固定次数，多个格子的算术混合矩可与独立、方差已标定的正态变量比较，绝对误差节省 $P$ 的一个固定幂。转回零点时，端点模 $6$ 余 $5$ 这一标记仍然保留。取两个半宽 $O(\log\log P)$ 的高度窗，并叠加公共端点项，在每个大二倍区间都得到一个标准化后趋于正态的零点统计量。当两次观测使用四个分离高度时，同一符号类的相关系数趋于 $1$，不同符号类的相关系数测量孪生权重。
- 对一次真实的辅助素数抽样，从候选质量中减去上述平方和后，得到原窗内真孪生质量的 $1/16$，误差在 $L^4$ 中为 $O((K\log P)^{-1/2})$。该界先对所有可用窗口取上确界，因此可以随后再用上同调条件选窗。
- 也可直接把不同符号类在两对分离高度窗上的零点和相乘。这个四次和以真孪生质量的 $1/8$ 为中心值，第 5 节式 (18c) 在任意预先指定的窗口上给出由零点本身标定的正态误差。它为正的概率以 $1/2$ 为门槛，能检测低于陈氏主尺度的孪生质量。若允许在 $d$ 个不相交的候选窗口间看完数据再选择，零质量假设下相应的单边通过概率上限为 $1-2^{-d}+o(1)$。本文没有证明通过概率实际超过门槛。

## 1. 原窗口、谱权和算术质量

沿用 ζ 目标。固定原证书的幂指数 $J$ 和 $0<\epsilon<1/52$，令 $P\to\infty$，$H=P^J$，$z=P^{\epsilon/5}$，$z'=z\log^2P$。先均匀选取 $e^{z/2}<Q\le e^z$、模 $8$ 余 $3$ 的素数。给定 $Q$，再均匀选取 $e^{z'/2}<D\le e^{z'}$、模 $16$ 余 $7$，且在 $Q$ 和每个不超过 $z$ 的奇素数模下余 $-1$ 的素数。所有格子始终使用这同一次抽样。记 $\lambda_v(n)=(n/v)$。这里的分布以辅助素数为随机变量，高度窗中的零点按重数计。[两族的定义](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/MeanSquareBounds/PrimeSetTwoBounds.lean#L37) · [第二阶段对每个第一阶段素数可取](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/MeanSquareBounds/CharOf.lean#L56)

只用 $p\equiv5\bmod6$ 的素数端点，记 $x=p+2$。取前文已证存在的一组固定常数 $3/11<a<b<1/2$。$\mathcal T$ 表示 $x$ 也是素数的端点，$\mathcal B$ 表示 $x=rs$、$r,s$ 均为素数且 $x^a\le r\le x^b$ 的端点。陈氏定理及端部因子上界筛给出

$$
\#((\mathcal T\cup\mathcal B)\cap(U,2U])\asymp U/\log^2 U,\qquad
w(p)=\frac{\log p\log(p+2)}p .
\tag{1}
$$

这里的 $a,b$ 由下界选择。小规模复算另行给出具体指数。[Green–Tao，定理 6.1 与脚注 3](https://www.numdam.org/article/JTNB_2006__18_1_147_0.pdf#page=20) · [内部因子区间的推论与证明](https://lixiang90.github.io/vibemath-NumberTheory/zeros-exact-distribution/articles/gap-two-separated-zero-windows.html)

在网格 $(U,2U]$、$U=2^jP$ 上等分 $K$ 个左开右闭的区间 $C$，只用 $[P,H]$ 内所需的格。设 $n=\log P$，$K$ 为正整数。以下的四阶一致估计适用于 $K\le P^{a/100}$，正态极限使用 $\log K=o(\log P)$。物理窗口及其归一化为

$$
I_Y=(2Y,Y^{3/2}/2],\qquad L_Y=\sum_{Y\le m\le Y^{3/2}/2}\frac1m,\qquad
T_Y=\sum_{p\in\mathcal T\cap I_Y} w(p),\qquad M_T(Y)=T_Y/L_Y ,
\tag{2}
$$

其中 $P\le Y$, $Y^2\le H$，$m$ 取整数。记 $\mathcal I(Y)=\{C:C\subset I_Y\}$，因此首尾最多漏掉两个不完整格。记 $t_C=\sum_{\mathcal T\cap C}w(p)$，$b_C=\frac1{128}\sum_{\mathcal B\cap C}w(p)$。

为便于复算，给出零点量的定义。取窗心 $h_1,h_2$，可随 $p$ 和第一阶段 $Q$ 变化，在抽 $D$ 之前确定。令 $\eta$ 是模 $3$ 的非主特征，$\chi_8=(8/\cdot)$。对 $\omega=1,\chi_8$ 用相应的 $h=h_1,h_2$，定义

$$
F_\omega(c,y)=-\frac2{p-1}\sum_{\chi(-1)=-1}\overline{\chi(c)}
    \sum_{\substack{L(\rho,\chi\eta\omega)=0,\ 0<\Re\rho<1\\|\Im\rho-h|\le R_h}}
          y^\rho H_N(\rho-ih).
\tag{3}
$$

外层取模 $p$ 的奇特征。这里 $H_N(s)=2\pi(32J_0)^{-1}\mathbb E\exp[-\cos((s+iV)/32)]$，$J_0=\int_\mathbb R e^{-\cosh u}\,du$。$V$ 是 $2N$ 个独立 $[-4,4]$ 均匀变量之和，再加独立变量 $2j$，$\Pr(j=k)=(N-|k|)/N^2$，$|k|<N$。设 $D_h=2+\log(24p(3+|h|))$，$R_h=10N-2+32\log(64[\log(24p)+\log D_h])$。前后的两个函数分别有导子 $3p,24p$。简记 $F_{\chi_8}=F_8$，$\ell=\log x$，$a_d=2d^{-1}\bmod p$，$\mathcal D_p=\{d:x^a\le d\le x^b,\ (d,6)=1\}$。定义实系数

$$
\begin{aligned}
E_p&=\ell^{-2}\Re\big[\chi_8(x)F_1(2,x)\overline{F_8(2,x)}\big],\\
k_{p,d}&=\frac{\chi_8(x)\chi_8(d)}{2\log d\,\log(x/d)}
       \Re\big[F_1(d,d)\overline{F_8(a_d,x/d)}\big],\\
K_{p,\sigma}&=\sum_{d\in\mathcal D_p} k_{p,d}\,
     \frac{(1+\sigma_1\lambda_Q(d))(1+\sigma_2\lambda_D(d))}{4},
      \qquad \sigma\in\{-1,1\}^2.
\end{aligned}
\tag{4}
$$

记 $I_p$ 为 $(\lambda_Q(p),\lambda_Q(x),\lambda_D(x),\lambda_D(p))=(1,1,-1,-1)$ 的指示函数，$W_p$ 只要求前三个符号。每个物理格中构造

$$
\begin{aligned}
Z_{C,\sigma}&=-\sum_{p\in C}\sqrt{w(p)}\,W_p\lambda_D(p)K_{p,\sigma},&
H_C&=\tfrac14\sum_\sigma Z_{C,\sigma}^2,\\
A_C&=\sum_{p\in C} w(p)I_p\big(E_p+\tfrac14\sum_\sigma K_{p,\sigma}\big),&
\mathcal F_K(Y)&=\frac1{L_Y}\sum_{C\in\mathcal I(Y)}(A_C-H_C).
\end{aligned}
\tag{5}
$$

平方前取 $\sqrt{w(p)}$，使对角项恢复原窗口所用的 $w(p)$。每个 $Z$ 都在一个格内先求和，平方保留不同导子之间的交叉项。也可把 $m$ 个格依次编号为 $0,\ldots,m-1$，对各格的 $Z$ 乘 $e^{2\pi i j\nu/m}$ 后相加。按 $\nu$ 平均其模平方，由有限 Fourier 正交性恰得 $\sum_C Z_{C,\sigma}^2$。构造 (5) 只输入导子端点、特征和指定零点窗，无需输入 $p+2$ 的素性或因子表。

## 2. 在双对数宽的零点窗中控制截断误差

以 $k^0_{p,d}$ 表示 (4) 的直接算术项。若 $d\mid x$，它等于 $\Lambda(d)\Lambda(x/d)/(2\log d\,\log(x/d))$，否则为零。前篇逐项显式公式给出

$$
\left|E_p-\frac{\Lambda(x)^2}{\ell^2}\right|
+\sum_{d\in\mathcal D_p}|k_{p,d}-k^0_{p,d}|
\ll e,\qquad
N=\left\lceil\frac{\log(2K)+4\log\log(3P)}{\log4}\right\rceil,\qquad
e=4^{-N}+P^{-2}.
\tag{6}
$$

第二项中的绝对值可以放在 $d$ 求和之内。前篇证明在替换 $F_1(d,d)$、处理对偶项和按 $d m=j p\pm2$ 合并远端项时，分别使用绝对上界，其中合并远端项用到 $\sum_{d\mid k}\Lambda(d)\Lambda(k/d)\le\log^2 k$。这也说明误差系数可在固定 $Q$ 后先于 $D$ 的抽样给定。若 $|h_i|\le\exp(p^c)$、$c>0$ 固定，(3) 的半宽为 $O_{J,c}(\log(2K)+\log\log P)$。[显式公式及远端项的逐项证明](https://lixiang90.github.io/vibemath-NumberTheory/zeros-exact-distribution/articles/gap-two-separated-zero-windows.html)

**端点标记。** 在一个格 $C\subset(U,2U]$ 内，固定 $Q$。展开 $W_p\lambda_D(p)K_{p,\sigma}$ 的第二阶段符号后，每一项的平方自由标签都有一个且只有一个大于 $U$、模 $6$ 余 $5$ 的素因子，即 $p$。事实上 $d<\sqrt{x}<U$。若 $x$ 是合数，其素因子均小于 $U$，若 $x$ 是素数，则 $x\equiv1\bmod6$。从第二阶段基模中取出已固定符号的小素因子不影响这个标记。这个观察对每个谱误差 $k_{p,d}-k^0_{p,d}$ 同样成立，包括 $d\nmid x$ 的项。

这里用归档引理 $\mathrm{QH.twistK}$。在第一阶段以及给定 $Q$ 的第二阶段，各自的奇平方自由标签 $1<u\le P^A$，与相应基模互素时，除至多一个指标外满足 $|\mathbb E\lambda(u)|\le P^{-B}$，其中 $A,B$ 是任意预先固定的正数。由固定的辅助剩余类可把原定理中的 Jacobi 符号换成这里的方向。例外指标在第二阶段可以依赖 $Q$。[原引理及全部量词](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/MeanSquareBounds/twistK.lean#L145)

若每行关于 $D$ 的展开系数绝对和至多 $c_p$ 且都有上述标记，对固定整数 $m\ge1$，

$$
\mathbb E_{D\mid Q}\left|\sum_{p\in C}\mathrm{row}_p(D)\right|^{2m}
\ll_m \left(\sum_{p\in C}c_p^2\right)^m+
    P^{-B}\left(\sum_{p\in C}c_p\right)^{2m}.
\tag{7}
$$

平凡标签要求奇重数端点的集合为空，一个例外标签也只允许一个指定集合。给该集合的每个端点预留一次出现，再将其余位置配对。指定端点的每个权重至多为 $(\sum c_p^2)^{1/2}$，每对自由位置的权重和至多为 $\sum c_p^2$，排列的常数只依赖 $m$。所有普通非平凡标签用原扭转估计。所涉及的标签乘积均在 $P^{O(mJ)}$ 内，可以把 (7) 的末项在下文中压到任意固定负幂。

每格长 $U/K\gg U^{1-a/100}$。一维和二维上界筛分别给 $\sum_{p\in C}w(p)\ll n/K$ 及 $t_C+128 b_C\ll1/K$，前者对全部素数，后者筛去 $p,p+2$ 的小素因子。可在格长 $V$ 内对奇素数筛至 $V^{1/24}$，所用筛限小于内部因子下界。区间的平移不改变上界。对 $K=P^{o(1)}$ 同样适用。[Montgomery–Vaughan 的算术大筛，推论 1](https://personal.science.psu.edu/rcv4/personal/Publications/large_sieve.pdf#page=2)

令 $\widetilde Z_{C,\sigma}$ 只保留属于 $\mathcal B$ 的真实算术行。由 (7) 对 (6) 使用 $c_p\ll e\sqrt{w(p)}$，谱误差的 $L^m$ 范数至多 $O_m(e\sqrt{n/K}+P^{-10J})$。还须删去有高次素数幂的因子行。整段中这种 $p$ 至多有 $O(U^{1-a/2}\log U)$ 个，因为 $x$ 被某个至少为 $U^a$ 的高次幂整除。对这些高次幂逐个计倍数并作部分求和便得此界。每行算术总权有界，(7) 取 $c_p\ll\sqrt{w(p)}$。故若 $\kappa=a/8$，对每个固定 $m\ge1$ 一致有

$$
\|Z_{C,\sigma}-\widetilde Z_{C,\sigma}\|_m
 \ll_m P^{-\kappa}+e\sqrt{n/K},\qquad
\|\widetilde Z_{C,\sigma}\|_m\ll_m K^{-1/2}.
\tag{8}
$$

这一步使用了零点窗心在抽 $D$ 前确定。原来的物理起点 $Y$ 可以在抽完 $D$ 后决定，两种选择发生在不同处。

## 3. 不同导子区间的混合矩

给每个格和每个符号类引入独立实正态变量 $G_{C,\sigma}$，均值零，方差 $b_C$，允许方差为零。对任意固定 $k$ 和任意 $k$ 个可以重复的格与符号指标，有

$$
\left|\mathbb E\prod_{i=1}^k\widetilde Z_{C_i,\sigma_i}
       -\mathbb E\prod_{i=1}^k G_{C_i,\sigma_i}\right|
\ll_{k,J} P^{-\kappa}.
\tag{9}
$$

各格可位于 $[P,H]$ 中相距很远的二倍区间。(9) 外层的真实期望始终按原两阶段次序取。证明分为两步。

**先更换辅助符号。** 每个真实半素数行的贡献为

$$
-\tfrac12\sqrt{w(p)}\,\lambda_D(p)
 \mathbf1_{(\lambda_Q(p),\lambda_Q(r),\lambda_Q(s))=(1,\sigma_1,\sigma_1)}
 \mathbf1_{(\lambda_D(r),\lambda_D(s))=(\sigma_2,-\sigma_2)} .
\tag{10}
$$

在每个母区间内，一个固定素数可属于至多 $O(U^{1-a}/\log^2U+1)$ 个这样的三元组。固定因子 $r$ 或 $s$ 后，对另一因子 $v$ 和 $rv-2$ 或 $sv-2$ 作二维上界筛，变量区间的长至少为常数乘 $U^a$。逐因子的上界分别为 $O(U/(r\log^2 U))$ 和 $O(U/(s\log^2 U))$，固定端点只有一行。故在每个格内，所有含某个指定素数的行，其 $w$ 总质量为 $O(U^{-a})$。[逐因子筛和例外标签的处理](https://lixiang90.github.io/vibemath-NumberTheory/zeros-exact-distribution/articles/gap-two-chen-auxiliary-switch.html)

固定实际的 $Q$。若第二阶段允许的非平凡例外指标能出现在这 $k$ 个因子的展开中，选它的一个素因子，从全部相关格中删去含该素数的行。若不能出现，直接保留各行。由 (7)，每格删去部分的 $L^k$ 范数为 $O_k(P^{-a/2})$，这对独立均匀的第二阶段素数符号也成立。逐因子展开乘积之差，用 Hölder 不等式和 (8) 的单格矩界，删除和加回的代价为 $O_k(P^{-a/2})$。删行后各非平凡项都受扭转引理控制，整个展开的绝对系数和至多 $P^{O(kJ)}$，可把比较误差取为 $O(P^{-10J})$。这样就把第二阶段换成了对各素数独立抛符号，且这些符号独立于 $Q$。

对第一阶段的例外指标作同样删除。先固定独立的第二阶段符号，比较真实 $Q$ 与独立第一阶段符号。删除和加回时先平均独立的第二阶段，由 (7) 仍有同一代价。各真实半素数的 $p,r,s$ 均大于 $z$，在两阶段都无需处理预先固定符号的基模小素数。至此两套素数符号都已独立，而同一个素数在不同格出现时仍取同一符号值。

**从最大的格向下取条件期望。** 在独立模型中，固定全部第一阶段符号和不超过当前格左端的全部第二阶段素数符号。较早格中的所有变量已经确定。当前格的所有 $r,s$ 小于母区间左端 $U$，故 (10) 中两个指示函数也已确定，剩余的 $\lambda_D(p)$ 在新端点上独立且均匀。每行只进一个 $\sigma$，所以当前四格符号和条件独立。写 $v_{C,\sigma}$ 为其条件方差，它是入选行的 $w(p)/4$ 之和。每行入选概率为 $1/32$，故 $\mathbb E v_{C,\sigma}=b_C$。不共享素数的三元组，其入选事件独立。按上面的逐素数质量界逐行计共享标签，得到

$$
0\le v_{C,\sigma}\ll1/K,\qquad
\operatorname{Var}v_{C,\sigma}\ll P^{-a}.
\tag{11}
$$

给定这些符号，当前格内奇数阶矩为零。偶数阶 $2l$ 的条件矩与 $(2l-1)!!\,v_{C,\sigma}^l$ 之差为 $O_l(\max_{p\in C} w(p))$。将均匀正负号的偶数矩按位置配对，只有同一端点占至少四个位置时与正态配对有差异，其总差由 $\sum c_p^4\le(\max c_p^2)\sum c_p^2$ 控制，其中 $|c_p|\le\sqrt{w(p)}/2$。$l=1$ 时精确相等。条件于当前格左端后，四个符号类可以相乘。再用 (11) 和方差的有界性，将所有 $v_{C,\sigma}$ 换成 $b_C$ 的 $L^2$ 代价为 $O_k(P^{-a/2})$。它乘上较早格的因子时用 Cauchy 不等式，所需高阶矩由 (7) 保证。从最后一格向前迭代便证得 (9)。这样允许晚格的 $r$ 或 $s$ 正好等于早格的 $p$。

**零点分布推论。** 设 $K=P^{o(1)}$，在固定有限个不同的指标 $(C,\sigma)$ 上 $K b_C\to\beta_C$。则实际零点变量 $\sqrt K\,Z_{C,\sigma}$ 联合趋于相互独立的 $\mathcal N(0,\beta_C)$。由 (9)，任一固定阶矩的比较误差乘上 $K^{k/2}$ 仍趋零，(8) 乘上 $\sqrt K$ 后也趋零。用正态分布的矩唯一性即得。特别地，前篇单个二倍段的正态波动，对于预先指定的两个高度窗，半宽可从 $O(\log p)$ 收窄至 $O(\log\log p)$。

**每个大二倍区间上的非退化极限。** 加入一个四个符号类共用的端点量，用上标 $\boldsymbol h=(h_1,h_2)$ 记录两个谱窗，定义

$$
D_C^{\boldsymbol h}=-\sum_{p\in C}\sqrt{w(p)}\,W_p\lambda_D(p)E_p,\qquad
X_{C,\sigma}^{\boldsymbol h}=D_C^{\boldsymbol h}+Z_{C,\sigma}^{\boldsymbol h},
\qquad d_C=t_C/8 .
\tag{11a}
$$

这里 $D_C^{\boldsymbol h}$ 是观测值，$D$ 仍表示第二个辅助素数。真实端点部分只在孪生行上取算术值 $-\sqrt w\,W_p\lambda_D(p)$，每行入选 $W_p$ 的概率为 $1/8$。端点高次幂在母区间只有 $O(\sqrt U\log U)$ 行，其 $L^m$ 贡献及谱误差按 (7) 与 (8) 处理。将真实孪生端点作为第五个坐标加入 (9) 时，相应的比较变量独立于原四个，方差为 $d_C$。具体地，删除例外标签的任何一个素因子时，在一个格最多删一条孪生行，因为 $p,x$ 分属不同的模 $6$ 类。按格向下取条件期望时，再预先给出全部模 $6$ 余 $1$ 的素数的第二阶段符号。这样当前的孪生大端 $x$ 已知，而每个当前小端 $p$ 的第二阶段符号仍然独立。五类入选行的端点互不相交，孪生坐标的条件方差均值为 $d_C$，方差至多 $O(\max_{p\in C}w(p))$，因为不同孪生行的两组素数标签互不重合。余下混合矩证明与 (9) 相同。

取 $K=1$、$C=(U,2U]$。由 (1)，$d_C+b_C\asymp1$ 在每个充分大的二倍区间都成立。因此对每个固定 $\sigma$，沿任意允许的 $U\to\infty$，

$$
\frac{X_{C,\sigma}^{\boldsymbol h}}{\sqrt{d_C+b_C}}\ \Longrightarrow\ \mathcal N(0,1).
\tag{11b}
$$

证明对任何子列再抽取 $d_C,b_C$ 分别收敛的子列即可。由二阶矩比较，分母也可替为完全由零点给出的 $(\mathbb E|X_{C,\sigma}^{\boldsymbol h}|^2)^{1/2}$。这个正态极限无须预设孪生或半素数中的任何一类单独有正密度。

用另两个窗心 $\boldsymbol h'$ 观测同一区间。允许四个高度窗两两分离，例如中心取 $0,e^p,2e^p,3e^p$。两次观测的真实算术部分完全相同。由 (6)、(7)，每个符号类有 $\|X_{C,\sigma}^{\boldsymbol h}-X_{C,\sigma}^{\boldsymbol h'}\|_2\ll e\sqrt n+P^{-10J}$。标准化后两个中心化变量的 $L^2$ 距离的平方等于 $2(1-\operatorname{Corr})$，故

$$
\begin{aligned}
1-\operatorname{Corr}(X_{C,\sigma}^{\boldsymbol h},X_{C,\sigma}^{\boldsymbol h'})
  &\ll e^2 n+P^{-10J},\\
\operatorname{Corr}(X_{C,\sigma}^{\boldsymbol h},X_{C,\tau}^{\boldsymbol h'})
  &= r_C+o(1),\qquad
 r_C=\frac{d_C}{d_C+b_C}=\frac{16t_C}{16t_C+\sum_{p\in\mathcal B\cap C}w(p)}
 \quad(\sigma\ne\tau).
\end{aligned}
\tag{11c}
$$

第二行来自唯一的共有坐标 $D_C$。结合二维正态极限，在不同符号类上还有

$$
\Pr\!\left(X_{C,\sigma}^{\boldsymbol h}X_{C,\tau}^{\boldsymbol h'}>0\right)
   =\tfrac12+\frac{\arcsin r_C}{\pi}+o(1)
   \qquad(\sigma\ne\tau).
\tag{11d}
$$

对二维独立标准正态向量作两次相关投影，其方向均匀，异号两扇区占的角度比例为 $\arccos(r_C)/\pi$。取收敛子列处理变化的 $r_C$，极限边缘无零原子，亦包括 $r_C\to1$。由此得到一个跨四个高度窗的单段检验。同号概率总是至少 $1/2-o(1)$。它在无穷多段上固定超过 $1/2$，等价于这些段中可取一列 $t_C\gg1$，即每段 $\gg U/\log^2 U$ 对真孪生素数。此处的相关性来自对辅助符号加权的导子族观测，没有使用单个 ζ 函数的零点间距猜想。

## 4. 事后选出的幂次窗口

**一致估计。** 在原两阶段抽样下，

$$
\left\|\sup_{P\le Y,\,Y^2\le H}
  \left|\mathcal F_K(Y)-\frac{M_T(Y)}{16}\right|\right\|_4
       \ll_J \frac1{\sqrt{K\log P}}
       \qquad(1\le K\le P^{a/100}).
\tag{12}
$$

常数也可依赖固定的 $\epsilon,a,b$。若 $K=P^{o(1)}$，左侧范数可换为任意固定有限阶 $L^m$，常数还可依赖 $m$。

**平方和的误差。** 记 $\widetilde H_C=\tfrac14\sum_\sigma\widetilde Z_{C,\sigma}^2$。由 (9)，在任意固定阶的展开中，它可以与相互独立的 $H_C^\circ=(b_C/4)\chi^2_{4,C}$ 比较，每个单项式误差为 $O_k(P^{-\kappa})$。由此

$$
\mathbb EH_C^\circ=b_C,\qquad
\operatorname{Var}H_C^\circ=b_C^2/2,\qquad
\mathbb E(H_C^\circ)^2=3b_C^2/2 .
\tag{13}
$$

全范围至多有 $N_*\ll_J Kn$ 个格。对任意连续 $h$ 格，独立模型中的中心化和的四阶矩至多 $C h^2/K^4$，因为每个 $b_C\ll1/K$，独立展开仅保留四次同格或两个成对的格。转回 $\widetilde H$ 的误差为 $O(h^4P^{-\kappa})$。对 $K\le P^{a/100}$，有 $K^6 n^2P^{-a/8}=o(1)$，故同一个四阶上界仍然成立。对于 $K=P^{o(1)}$，每个固定 $2m$ 阶的对应界为 $O_m(h^m/K^{2m})$，因为含单次下标的项消失，误差 $O_m(h^{2m}P^{-\kappa})$ 可吸收。

对格子序列作二进剖分，必要时在末端补零。第 $l$ 层至多 $2^l$ 个段，每段长 $O(N_*2^{-l})$。该层最大绝对和的 $L^4$ 范数至多 $O(\sqrt{N_*}\,2^{-l/4}/K)$，对层求和即控制任意前缀。任意连续格可由两个前缀作差。因此，全体 $\mathcal I(Y)$ 上未归一化中心平方和的上确界，其 $L^4$ 范数是 $O_J(\sqrt{n/K})$。相同证明在 $2m>2$ 时用衰减因子 $2^{-l(1/2-1/(2m))}$，较低阶由高阶支配。

由 (8)，各格的两种 $Z$ 的固定阶范数均为 $O_m(K^{-1/2})$。将 $\widetilde H$ 全部换回真实零点的 $H$，用 Hölder 不等式，再对所有格的绝对差求和，代价至多

$$
O_{m,J}\!\left(n\sqrt K\,P^{-\kappa}+ e n^{3/2}\right)
=o\!\left(\sqrt{n/K}\right).
\tag{14}
$$

**事后挑选谱窗。** 若比较两份谱窗规定，算术行及高次幂完全相同，(8) 的差只余 $O_m(e\sqrt{n/K}+P^{-B})$，其中 $B$ 可任意固定。候选和的总谱误差为 $O(e n^2)$。对 (5) 中未归一化的 $R_K=L_Y\mathcal F_K$，由此得到

$$
\left\|\sup_Y |R_K^{\boldsymbol h}(Y)-R_K^{\boldsymbol h'}(Y)|\right\|_m
 \ll_{m,J} e n^2+P^{-10J}.
\tag{14a}
$$

甚至可以在给定 $Q$ 后预先列出至多 $(Kn)^B$ 份完整的逐 $p$ 谱窗规定，其中 $B$ 固定，再看 $D$ 来选一份。将 (14a) 对各份与同一参考规定比较，以固定的 $l>\max(m,4B)$ 阶取最大值，代价至多 $(Kn)^{B/l}$。由 (6)，总差仍为 $o_{L^m}(\sqrt{n/K})$，低于第 5 节的误差尺度。这里的数目限制针对完整规定。事后对每个 $p$ 独立任择窗心会产生更多完整规定。

**候选和与边界。** 在 $A_C$ 中先换回真实算术权。端点是素数时常数符号项给 $t_C/16$。半素数四格平均后常数项给 $b_C$。对其余符号单项式，在含有非空标签的最后阶段按平方自由核合并。若各系数合并为 $c_u$，包括唯一扭转例外时的二阶上界为 $2\sum|c_u|^2+P^{-B}(\sum|c_u|)^2$，因为例外只把每个核 $u$ 连到一个核。每个非空纤维的绝对权至多 $O(U^{-a})$，这是逐素数质量界和互补因子唯一分解的推论。孪生行的两个端点分属模 $6$ 的不同类，每个非空纤维只有一行。各格总权有界，故符号误差的 $L^2$ 范数为 $O(P^{-a/2})$。算术高次幂误差至多 $O(U^{-a/2}\log^3U)$，也包括端点为高次幂的情形。由确定性有界性，算术 $A_C-t_C/16-b_C$ 的 $L^4$ 范数为 $O(P^{-\kappa})$，任意固定阶也节省一个依赖阶数的正幂。

总共有 $O_J(Kn)$ 格。候选和中的确定性谱截断误差，在所有格上合计至多 $O(e n^2)$，这里用 $\sum_{P<p\le H}w(p)\ll_J n^2$。这些误差均为 $o(\sqrt{n/K})$，四阶时上述固定幂允许 $K\le P^{a/100}$。舍去的两端不完整格的真孪生权重至多 $O(1/K)$，这点在 $K$ 大于 $\log P$ 时也须保留。最后 $L_Y\asymp_J n$，得到 (12)。

**原上同调选择。** 原构造先取辅助素数，再在所许范围内取完整 $[Y,Y^2]$。固定参数后，有至多 $R_0$ 个几何例外扭转指标。可事先列出 $d=R_0+1$ 个起点 $Y_i=P^{(2C_t/\eta)^i}$，各自可能冲突的 $(Y_i^\eta,Y_i^{2C_t}]$ 两两不交。每个解析好对至少有一段可用。解析好对的概率趋于一，且某个 $i=i(P)$ 被正比例的辅助对共用。在任意概率至少为 $\alpha>0$ 的原好对事件上，即使 $Y$ 随辅助对变化，(12) 在条件期望中给出的绝对误差至多 $O_J(\alpha^{-1/4}/\sqrt{Kn})$。因此若无穷多个基准上，这样的共同好窗口有正比例的辅助对满足 $\mathcal F_K(Y)\ge c>0$，则可选到 $M_T(Y)/16\ge c-o(1)$。前篇原四个端点符号的同时估计，再给 $\Gamma_1\ge c/3$。目前 (12) 约束的是选窗偏差，$\mathcal F_K(Y)$ 的固定正下界仍未证得。[原有限例外及区间选择](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/AutomorphicVarieties/LayeredIsomorphisms/ThmBabOf.lean#L85) · [先取素数对再取原窗](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/CharacterBoxes/WinPairs.lean#L32) · [端点质量与原归一化](https://github.com/lixiang90/vibemath-NumberTheory/blob/main/Goldbach_TwinPrimes/twin-prime-auxiliary-signs.md)

## 5. 低于陈氏主尺度的孪生质量

以下取 $K=P^{o(1)}$，并在抽辅助素数前指定 $Y$。将真实误差的标准差参数及其零点估计分别写作

$$
s_Y^2=\tfrac12\sum_{C\in\mathcal I(Y)} b_C^2,\qquad
\widehat s_Y^2=\tfrac13\sum_{C\in\mathcal I(Y)} H_C^2,\qquad
R_K(Y)=L_Y\mathcal F_K(Y).
\tag{15}
$$

若 $T_Y=o(\log Y)$，则沿这列尺度有

$$
s_Y^2\asymp_J (\log Y)/K,\qquad
\widehat s_Y^2/s_Y^2\ \longrightarrow\ 1\ \text{依概率},\qquad
\frac{R_K(Y)-T_Y/16}{\widehat s_Y}\ \Longrightarrow\ \mathcal N(0,1).
\tag{16}
$$

若分母偶尔为零，商在该事件上可取零。由第二式，该事件概率趋于零。

**证明。** 对 $I_Y$ 内的完整二倍区间用 (1)，总候选权重 $\asymp\log Y$。扣去真孪生权和至多两个不完整格后，$\sum b_C\asymp\log Y$。上界 $b_C\ll1/K$ 和 Cauchy 不等式下界（格数 $O(K\log Y)$）给出 (16) 首式。特别地，$\max_C b_C/s_Y\ll1/\sqrt{K\log Y}\to0$，即使 $s_Y$ 本身趋于零也成立。

先对独立的 $H_C^\circ$ 作计算。中心化和的方差正是 $s_Y^2$。标准化后任意固定阶数的矩趋于实标准正态矩。事实上，令 $U_C=H_C^\circ-b_C$，则 $|\mathbb E U_C^l|\ll_l b_C^l$。展开时各下标出现一次的项为零，某个下标出现三次或更多时，由 $\sum(b_C/s_Y)^2=2$ 可将其总和界为 $O_l(\max b_C/s_Y)$。其余完全配对，给出正态矩。由 (9)，将每格换回 $\widetilde H_C$ 时第 $l$ 阶原始展开的误差为 $O_l((Kn)^l P^{-\kappa})$，再除以 $s_Y^l$ 依然趋零，因为 $K=P^{o(1)}$。用 (14)、候选和估计及端格的 $O(1/K)$，所有被删去的项合计为 $o_{L^2}(s_Y)$。这证明以 $s_Y$ 作分母的正态极限。

$H_C^\circ/b_C$ 在 $b_C>0$ 时服从形状 $2$、尺度 $1/2$ 的 Gamma 分布，故 $\mathbb E(H_C^\circ)^2/3=b_C^2/2$ 且 $\operatorname{Var}((H_C^\circ)^2/3)=7b_C^4/12$。由独立性及 (9)，用 $\widetilde H_C$ 估计 $s_Y^2$ 的均值误差为 $O(KnP^{-\kappa})$，方差为 $O(\sum b_C^4+(Kn)^2P^{-\kappa})=o(s_Y^4)$。将 $\widetilde H_C^2$ 换回 $H_C^2$，由 (8) 及 Hölder 不等式，单格 $L^1$ 误差为 $O(K^{-3/2}(P^{-\kappa}+e\sqrt{n/K}))$，合计为 $o(s_Y^2)$。所以数据给出的分母也可用。

**单窗的符号判据。** 对任意预先指定的窗口序列，不必预知 $T_Y$ 是大是小，

$$
\limsup_{P\to\infty}\frac{T_Y}{\sqrt{(\log Y)/K}}>0
\quad\Longleftrightarrow\quad
\limsup_{P\to\infty}\Pr_{Q,D}\{R_K(Y)>0\}>\tfrac12 .
\tag{17}
$$

若左侧的比值趋零，(16) 给右侧概率趋于 $1/2$。若沿子列比值趋向一个有限正数，则 $T_Y=o(\log Y)$，且 $T_Y/(16s_Y)$ 有固定正下界。取子列令此移位收敛，正态极限给出严格超额。若比值趋于无穷，则 (12) 直接给概率趋于 $1$。这些子列覆盖全部情形。原解析好对事件的概率趋于一，故在该事件上条件化保留 (17)。[固定误差比例之前已选定原解析常数](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/MeanSquareBounds/PairShare.lean#L22)

例如取 $K=\lceil(\log P)^3\rceil$。此时两个零点窗均只需双对数半宽。右边超过 $1/2$ 一个固定量，等价于无穷多窗上 $T_Y\gg1/\log Y$。用以证明 $\Gamma_1>0$ 的质量尺度是 $T_Y\gg\log Y$。因此 (17) 的严格超额已足以推出孪生素数无限性，并容许 $M_T(Y)\to0$。单个孤立素数对在这段中至多贡献 $O((\log Y)^2/Y)$，所以任意稀疏的孪生素数序列未必达到 (17) 的门槛。

**事后挑选的概率代价。** 原有限列表中的物理窗两两不交。更一般地，设有固定 $d$ 个预先确定、互不相交的 $I_{Y_i}$。若各自 $T_{Y_i}=o(\sqrt{(\log Y_i)/K})$，则

$$
\left(\frac{R_K(Y_i)}{\widehat s_{Y_i}}\right)_{i=1}^d
\ \Longrightarrow\ \mathcal N(0,I_d).
\tag{18}
$$

证明对任意固定的坐标线性组合重复 (16) 的矩比较。各窗使用互不相交的格，比较模型中各窗独立，(9) 已处理真实抽样中跨窗复用素数的影响。因而允许自由挑选最大值时，$\Pr(\max_i R_K(Y_i)/\widehat s_{Y_i}>x)\to1-\Phi(x)^d$。对任何必须同时通过原上同调条件的选窗规则，这给出无条件通过概率的上界。若仅在概率 $1-o(1)$ 的解析好对事件上运行，界仍成立。在 $x=0$ 处就是 $1-2^{-d}$。若只在通过指定几何窗口条件的正比例辅助对子集内求条件概率，还须计入该子集的概率。若在原有限列表上以固定幅度突破无条件上界，则至少一个候选窗沿子列达到 (17) 左端的质量，已蕴含孪生素数无限性。现有的选窗引理给出覆盖率，未提供这一符号偏向。

**直接用四个分离高度窗计孪生质量。** 将因子符号 $(-,-),(-,+),(+,-),(+,+)$ 依次编号 $1,2,3,4$。在格 $C$ 上，前两个 $X$ 采用谱窗对 $\boldsymbol h$，后两个采用谱窗对 $\boldsymbol h'$，所得实数简记为 $x_{C,i}$。例如四个高度中心在全部导子上统一取 $0,P,2P,3P$，充分大时四个谱窗两两分离。令

$$
\begin{aligned}
\mathcal P_K(Y)&=\sum_{C\in\mathcal I(Y)}x_{C,1}x_{C,4},\\
\widehat V_Y&=\sum_{C\in\mathcal I(Y)}
\left(\frac16\sum_{i<j}x_{C,i}^2x_{C,j}^2
-\frac13 x_{C,1}x_{C,2}x_{C,3}x_{C,4}\right).
\end{aligned}
\tag{18a}
$$

此式的分子在每个指定高度窗各使用一个零点，保留同格内不同导子的所有交叉项。方差估计的每个括号对任意实数取值都非负，因为按两对互补下标用算术几何平均，第一项至少为 $|x_{C,1}x_{C,2}x_{C,3}x_{C,4}|$。

令 $y_C=d_C+b_C$、$V_Y=\sum_{C\in\mathcal I(Y)}(y_C^2+d_C^2)$。在 (11a) 的比较模型中，各坐标为 $X_i^*=D^*+G_i$，公共正态项方差 $d_C$，四个独立标记项各有方差 $b_C$。正态配对给出

$$
\begin{aligned}
\mathbb E X_1^*X_4^*&=d_C,\qquad
\operatorname{Var}(X_1^*X_4^*)=y_C^2+d_C^2,\\
\mathbb E (X_i^*X_j^*)^2&=y_C^2+2d_C^2\quad(i\ne j),\qquad
\mathbb E\prod_{i=1}^4 X_i^*=3d_C^2 .
\end{aligned}
\tag{18b}
$$

现在对于每一列预先指定的 $Y$ 和每个 $K=P^{o(1)}$，均有

$$
V_Y\asymp_J(\log Y)/K,\qquad
\frac{\widehat V_Y}{V_Y}\ \longrightarrow\ 1\ \text{依概率},\qquad
\boxed{\quad
\frac{\mathcal P_K(Y)-T_Y/8}{\sqrt{\widehat V_Y}}
\ \Longrightarrow\ \mathcal N(0,1).\quad}
\tag{18c}
$$

**证明。** 由 (1)，不管孪生和半素数的相对数量如何，$\sum y_C\asymp\log Y$，且 $y_C\ll1/K$。Cauchy 上下界给出 $V_Y$ 的数量级。各格中心化乘积的每个固定 $l$ 阶绝对矩至多 $O_l(y_C^l)$，而 $\max y_C/\sqrt{V_Y}\to0$。在独立格中可按 (16) 的配对论证得到标准正态极限。真实算术行用第 3 节含第五个坐标的 (9) 逐项比较，$O_l((Kn)^lP^{-\kappa})$ 除以 $V_Y^{l/2}$ 仍趋零。由 (8) 及端点的同一估计，将真实算术行换回两个谱窗时，每个乘积的 $L^2$ 误差为 $O(K^{-1/2}(P^{-\kappa}+e\sqrt{n/K}))$，合计为 $o(\sqrt{V_Y})$。遗漏的边格至多造成 $O(1/K)$ 的 $T_Y$ 误差。

由 (18b)，比较模型中 (18a) 方差估计的每个括号均值为 $y_C^2+d_C^2$，方差至多 $O(y_C^4)$，因为它是四次多项式。不同格独立且 $\sum y_C^4=o(V_Y^2)$。真实算术模型只多出固定阶矩的 $P^{-\kappa}$ 型误差，转回零点的单格 $L^1$ 误差至多 $O(K^{-3/2}(P^{-\kappa}+e\sqrt{n/K}))$。总误差为 $o(V_Y)$，证毕。偶遇 $\widehat V_Y=0$ 时仍将商记为零。

对 $x_{C,1}x_{C,4}-d_C$ 作第 4 节的连续格最大值估计，还得到

$$
\left\|\sup_Y\left|\frac{\mathcal P_K(Y)}{L_Y}-\frac{M_T(Y)}8\right|\right\|_4
\ll_J (K\log P)^{-1/2}\qquad(1\le K\le P^{a/100}).
\tag{18d}
$$

这里比较模型的每格方差 $O(K^{-2})$，四阶矩 $O(K^{-4})$，混合矩展开和谱误差的界与 (12) 相同。于是 (17) 中的 $R_K$ 可直接替换成 $\mathcal P_K$，(18) 对不相交窗口也有完全相同的选择上界。特别地，若在预定窗口上证明 $\Pr(\mathcal P_K>0)\ge1/2+\delta$ 沿无穷多个基准成立，其中 $\delta>0$ 固定，即可推出孪生素数无限性。若在有正比例辅助对共享的原好窗口上证明 $\mathcal P_K/L_Y\ge c>0$ 的联合事件有固定正概率，(18d) 和原端点符号平均则给 $\Gamma_1\ge c/6$。二者都是待证的单边估计。

## 6. 复算与可观测的选择偏差

取两个真实半素数邻项 $53+2=5\cdot11$、$1217+2=23\cdot53$。较早的端点 $53$ 后来再次作为因子出现。穷举五个素数标签的两套独立符号，共 $1024$ 种。指定一格 $\sigma$，两个线性和的交叉矩精确为零。两个平方的混合矩与各自均值之积相比，当 $\sigma_1=-1$ 时比值为 $0$，当 $\sigma_1=1$ 时比值为 $2$。这直接展示了有限尺度独立性会失败的情形，也说明第 3 节需要先估计条件方差再转移矩。

再取真孪生 $(59,61)$ 和半素数邻项 $53+2=5\cdot11$，分别给行权 $4,9$。穷举另一组 $1024$ 种赋值，得到四个 $X$ 的对角方差均为 $4/8+9/128$，非对角协方差均为 $4/8$。在下表样本的最大二倍段 $(65536,131072]$ 上，换用 $(-,-)$ 和 $(+,+)$ 两格，真实辅助素数给出的样本相关系数为 $0.960589$，(11c) 的算术比较值为 $0.963878$。样本同号比例为 $0.905382$，(11d) 给出的比较值为 $0.914184$。

令 $Y=4096$，则 $I_Y=(8192,131072]$ 恰含四个完整二倍区间。以 $[x^{3/10},x^{9/20}]$ 作有限样本的因子范围，数得真孪生 $1349$ 对、内部半素数 $826$ 个。对同一批 $48$ 个真实 $Q$ 及每个 $Q$ 下的 $24$ 个真实 $D$ 作计算，保留所有端点和因子的素数幂。这里 $Q\equiv3\bmod8$、$D\equiv7\bmod16$、$D\equiv-1\bmod105Q$，数值例子仅采用这些同余。表中各项未除以 $L_Y$。真值 $T_Y/16=0.223982$、$\sum b_C=0.017073$，$H$ 表示 $\sum_C H_C$。

| 每段格数 $K$ | 样本 $\mathbb E H$ | 样本标准差 | 由 $b_C$ 算的 $s_Y$ | 样本 $\sqrt{\mathbb E\widehat s_Y^2}$ | 样本 $\mathbb E R_K(Y)$ |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 0.018339 | 0.007179 | 0.006039 | 0.006737 | 0.222702 |
| 4 | 0.018226 | 0.003703 | 0.003068 | 0.003498 | 0.222815 |
| 16 | 0.018246 | 0.002467 | 0.001605 | 0.002047 | 0.222794 |
| 64 | 0.018225 | 0.001950 | 0.000896 | 0.001429 | 0.222816 |

对 (18a) 在同一批样本上另算如下，算术中心值为 $T_Y/8=0.447963123$。

| 每段格数 $K$ | 样本 $\mathbb E\mathcal P_K$ | 样本标准差 | 由算术质量算的 $\sqrt{V_Y}$ | 样本 $\sqrt{\mathbb E\widehat V_Y}$ |
| ---: | ---: | ---: | ---: | ---: |
| 1 | 0.442277 | 0.330187 | 0.323057 | 0.320051 |
| 4 | 0.442134 | 0.161355 | 0.165746 | 0.163700 |
| 16 | 0.444005 | 0.094359 | 0.084379 | 0.088509 |
| 64 | 0.441479 | 0.054144 | 0.044819 | 0.051667 |

细分后样本标准差下降，在 $K=64$ 时仍超过极限模型数值的两倍。这组端点有限、辅助素数取样也只覆盖有限同余条件，没有验证原渐近上同调条件。独立平方模型另外可以精确复算选择概率。具体取常数 $b_C>0$、每窗 $m$ 格和三窗，单窗正号概率为形状 $2m$、尺度 $1$ 的 Gamma 分布在 $2m$ 处的分布函数。当 $m=4,64,1024$ 时，三个窗中取最大值的正号概率分别为 $0.907064,0.883610,0.877191$，趋向 $7/8$。这里使用分布函数直接求值。

最后在物理格 $(50,70],(70,90]$ 的全部五个素数端点 $p\equiv5\bmod6$ 上，取 $Q=1000003,D=47880143639$、谱阶数 $N=2$，分别在 $320,960$ 和 $1920,2880$ 两组窗心用围道积分算 (4)。程序逐个导子核实四个高度窗两两分离。两次观测中，逐 $p$ 的 $\sum_d|k_{p,d}-k^0_{p,d}|$ 最大值分别为 $1.887\cdot10^{-15}$ 和 $1.783\cdot10^{-15}$。两格的 $\sum H_C$ 均为 $0.033544782$，换高度后四种 $Z$ 的最大改变量为 $8.882\cdot10^{-16}$，带公共端点的 $X$ 的最大改变量为 $1.567\cdot10^{-10}$。

采用第一对高度的 $(-,-)$ 格和第二对高度的 $(+,+)$ 格相乘，(18a) 的算术复算值为 $0.257589508775$。程序独立比较四窗围道值，同时核对两层积分网格、算术展开、各窗心处 Arb 的 $L'/L$ 和两个物理格上的 Fourier 平方恒等式。上述围道计算未逐个定位零点。笔记给出纸面推导，数值用于复核常数和有限实现。[复算程序及运行方法](https://github.com/lixiang90/vibemath-NumberTheory/blob/main/zeros-exact-distribution/gap-two-window-quadratic-variation/verify.py)
