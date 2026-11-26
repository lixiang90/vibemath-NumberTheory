---
title: "原互素对在对数频带上检出完整窗口的精确幂指数"
date: "2026-11-27"
abstract: "对判别式 229 的实三次域及其 CM 扩张，证明式（10a）中任意固定常数倍 log T 的频带已检出整窗累计正部的精确幂指数。余因子可同时限于 z 与常数倍 z log² T 之间，每块只需在该频带内取整数采样。给出素数平方引起的显式修正，并用区间算术检验两支特征交叉项的符号。"
---

# 原互素对在对数频带上检出完整窗口的精确幂指数

令 $F=\mathbb Q(\alpha)$、$\alpha^3-4\alpha+1=0$、$E=F(i)$。$M/B=F^{\rm gal}(i)/\mathbb Q(\sqrt{229},i)$ 的三阶特征给出 $L_B(s,\psi)=L_+(s)L_-(s)$。写 $a_\nu,u_\nu,b_\nu$ 为 $L_\nu,1/L_\nu,-L_\nu'/L_\nu$ 的系数。令 $\beta_*$ 为这两支的非平凡零点实部的上确界，$\alpha_*=\beta_*-1/2$，$\delta=3\alpha_*/2$。下文证明，**对每个固定 $c>0$，取 $|t|<c\log(2T)$ 就能检出 $\delta$，还可同时只用 $z<h\ll z\log^2(2T)$ 的余因子**。所有等式保留两支特征和两类余因子的原符号。[原频带与式（10a）](cm229-coprime-smoothing.md) · [两类余因子的零点系数](cm229-signed-zero-modes.md)

## 完整窗口上的频带定理

取 $X=Y^{3/2}$，$T=2^jY<X$，$j\ge0$，$m_Y=\lceil\log^2(2Y)\rceil$，$r=2^{1/m_Y}$，$v=\log r$。取 $Z=T^{1/6}/\log^3(2T)$，以及原截断 $z=\max(1,\lfloor Z\rfloor)$。原和为

$$
 Q_T(x)=\tfrac12\sum_{\nu=+,-}
 \sum_{\substack{h>z,\ p>T^{2/3}\ {\rm prime}\\2x/r<hp\le2x}}
 c_{z,\nu}(h)a_\nu(p)\log p,\qquad
 c_{z,\nu}=\mathbf1_{h=1}-u_{\nu,\le z}*a_\nu .
$$

将 $h,p$ 分成二进区间 $[D,2D),[P,2P)$，与各自的截断取交。按原式保留乘积范围与 $T/r<hp<8T$ 有交的矩形，保留后用整个矩形。这样在 $T/2<x<4T$ 还原原和，其余处也定义了一个对数变量上紧支撑的延拓。写 $s=1/2+it$，令

$$
 F_T(t)=\tfrac12\sum_{\nu,D,P}
 \left(\sum_h c_{z,\nu}(h)h^{-s}\right)
 \left(\sum_p a_\nu(p)(\log p)p^{-s}\right),\qquad
 {\cal E}_T(A;K)=\int_{-K}^K\min\{v^2,(1+|t|)^{-2}\}|A(t)|^2\,dt .
$$

再令 $H=D_0 Z\log^2(2T)$，$D_0$ 为下文确定的充分大绝对常数。$F_T^\Phi$ 在上述 $h$ 和中乘 $\Phi_T(h/H)$，$F_T^{\le U}$ 则保留 $h\le U$。这里 $\Phi_T(y)=\Pr(e^S\ge y)$，$S$ 是 $12\lceil\log(2T)\rceil+20$ 个独立均匀变量之和，每个变量的区间取 $[0,\log2/(12\lceil\log(2T)\rceil+20)]$。因此 $F_T^\Phi$ 是对 $H\le U\le2H$ 的硬截断作正权平均。充分大的 $T$ 有 $2H<z^2$，所留非零项只有素数余因子与全部素因子 $\le z$ 的余因子。[原系数的精确消失条件](cm229-coprime-cofactors.md)

**定理。** 固定 $c>0$，令 $K_c=K_c(T)=c\log(2T)$。$N_T$ 可取以下三个量中的任何一个。

$$
 {\cal E}_T(F_T;K_c)^{1/2},\qquad
 {\cal E}_T(F_T^\Phi;K_c)^{1/2},\qquad
 \max_{H\le U\le2H}{\cal E}_T(F_T^{\le U};K_c)^{1/2}.
$$

对每个选择，沿全部充分大的实数 $Y$ 都有精确极限

$$
 \lim_{Y\to\infty}\frac{\log(1+\max_{T=2^jY<X}N_T)}{\log Y}=\delta .
 \tag{1}
$$

前篇已证明 $\log(1+W^{-1}\sum_{Y\le n<X}k_n[Q_T(n)]_+)/\log Y\to\delta$，其中 $W=\lceil Y^3\rceil$、$k_n=\lfloor W/n^{3/2}\rfloor$。所以以上任一量的整窗最大值在任意远窗口上对每个 $\eta>0$ 满足 $O_\eta(Y^\eta)$ 的上界，当且仅当原正部满足该预算，也当且仅当 $E$ 满足 ERH。最后一步用 $\zeta_E=\zeta L(s,\chi_{-4})L_+L_-$、函数方程与有理 Dirichlet GRH。若 $\delta>0$，每个充分远的窗口都检出增长。这里 $t$ 是有限 Dirichlet 多项式的实频率。定理给的是随 $Y$ 的精确指数。对单一窗口仍须区分滤波后的量与原整数正部。[原正部的指数定理](cm229-type-ii-large-prime.md) · [Alpöge–Claude 的有理 Dirichlet 定理](https://math.mit.edu/~drew/grh/)

## 正权平均后的长余因子消去

对实变量函数 $U(x)$ 写 $f_U(w)=e^{-w/2}U(e^w)$。记 $G_a f(w)=\mathbb E f(w+\sqrt a\,N)$，$N$ 为标准正态变量，取 $a=(1024\log(2Y))^{-1}$。冻结每个块自己的参数。对 $R=r,2$ 置 $J_R(x)=\tfrac12\sum_{\nu}\sum_{2x/R<j\le2x}b_\nu(j)$。$Q_T^\Phi$ 表示在 $Q_T$ 中使用上述短余因子权。对 $T\le x\le2T$，一致有

$$
 \begin{aligned}
 G_af_{Q_T}(\log x)-G_af_{Q_T^\Phi}(\log x)&=O(T^{-8}),\\
 G_af_{Q_T}(\log x)-G_af_{J_r}(\log x)&=O_\epsilon(T^\epsilon z^{\alpha_*}).
 \end{aligned} \tag{2}
$$

第一式是带符号和的点态消去。它把这一平均中 $h$ 的上限明确给到 $2H$。第二式的误差在整窗累加只花 $O_\epsilon(Y^{\delta/6+\epsilon})$。高斯在原变量中的典型宽度约为 $T/\sqrt{\log(2Y)}$。对这种宽平均，以下证明保留正部的增长指数。前篇对较窄平均另有直接控制原整数绝对误差的定理。[原整数和的绝对平滑误差](cm229-coprime-smoothing.md)

为证 (2)，起初让 $h,p$ 取其整个无限范围。归一化 Mellin 反演核为

$$
 K_x(s)=x^{s-1/2}\frac{2^s(1-r^{-s})}{s}
                \exp\big(a(s-1/2)^2/2\big). \tag{3}
$$

在 $\Re s=2$ 各级数绝对收敛。写 $C_{z,\nu}=1-U_{z,\nu}L_\nu$、$U_{z,\nu}=\sum_{d\le z}u_\nu(d)d^{-s}$。对 $|t|\le K=256\log(2T)$，将 $C_{z,\nu}(2+it)$ 换成短余因子多项式只花 $O(T^{-20})$，$D_0$ 取为适合此精度和高度的固定常数。具体地，在前篇的光滑 Mellin 反演中移至 $\Re w=-A$，$A=2\lceil\log(2T)\rceil$。在 $\Re s=2$ 时从偶整数 $A=4$ 起使用函数方程，反射侧绝对收敛，两份伽马商仍由 $[C(1+A+|t+\Im w|)]^{2A+O(1)}$ 控制。$d$ 的求和以 $\sum_{d\le z}d_2(d)d^{A-2}\ll z^{A-1}\log(2z)$ 控制。余下的几何因子为 $(Cz(A+|t|)^2/H)^A$，增大 $D_0$ 即得所需精度。超过 $K$ 的 Mellin 积分有 $\exp(-aK^2/4)\le(2T)^{-16}$，其余高斯因子可积，足以吸收初始线的 $O(T^{3/2})$。换回矩形延拓仅影响 $|\log x'-\log x|\ge\log2$ 的高斯尾。原系数与矩形系数均逐个受除数函数乘对数控制，故归一化函数在 $x'\ge1/2$ 为 $O_\epsilon((2x')^{1/2+\epsilon})$。高斯尾的费用小于 $O(T^{-8})$。这证明第一式。[短余因子逼近的核与移线常数](cm229-coprime-cofactors.md)

对第二式逐支写 $B_\nu=-L_\nu'/L_\nu$，$B_{\nu,\mathrm{pp}}$ 为全部 $p^k,\ k\ge2$ 的级数，$P_\nu$ 为 $p\le T^{2/3}$ 的有限素数多项式。无限双和的级数恰为 $C_{z,\nu}(B_\nu-B_{\nu,\mathrm{pp}}-P_\nu)$。由于

$$
 C_{z,\nu}B_\nu=B_\nu+U_{z,\nu}L_\nu',\qquad
 C_{z,\nu}P_\nu=P_\nu-U_{z,\nu}L_\nu P_\nu ,
$$

$U_{z,\nu}L_\nu'$ 及 $C_{z,\nu}P_\nu$ 在 (3) 中可移至任意固定负线 $-A_0$。$L_\nu,L_\nu'$ 在该线上至多按 $|t|$ 的固定幂增长，新积分至多为 $O_{A_0}(T^{-A_0-1/2}(zT^{2/3})^{A_0+1}\log^{O_{A_0}(1)}(2T))$。因 $zT^{2/3}\le T^{5/6}$，可使它为 $O(T^{-8})$。最后 $B_{\nu,\mathrm{pp}}$ 在 $\Re s>1/2$ 绝对收敛。在 $\sigma=1/2+\epsilon'$ 处，用前篇的逆系数界 $|U_{z,\nu}(\sigma+it)|\ll_{\epsilon'}(2z)^{(\beta_*-\sigma)_++\epsilon'}(2+|t|)^{\epsilon'}$、$L_\nu$ 的多项式增长与高斯衰减，就得到 (2)。[截断逆系数的 Perron 估计](cm229-type-ii-large-prime.md)

## 保留正部并缩小频带

先说明 $J_r$ 的平均仍检出 $\delta$。在固定比例 $R=2$ 和固定方差 $1$ 下，全对数轴上的 $e^{w/2}G_1 f_{J_2}(w)$ 的 Laplace 变换为

$$
 e^{(s-1/2)^2/2}\frac{2^s-1}{2s}
 \left(-\frac{(L_+L_-)'(s)}{(L_+L_-)(s)}\right).
$$

该函数在 $w<0$ 的积分给出整函数，因为 $J_2(e^w)$ 在充分负的 $w$ 为零，高斯尾二次指数衰减。$w\ge0$ 时原实函数为 $O(e^w)$，这由 $|b_\nu(j)|\le2\Lambda(j)$ 得到。高度 $64$ 的完整计数已排除实非平凡零点，上式在正实轴无极点。在开临界带内，高斯和 $2^s-1$ 均不消去零点极点。于是可将 Alpöge–Claude 的完整窗口 Landau 引理用于这个固定函数。[窗口引理的 Lean 定理](https://github.com/leanprover-community/riemann-hypothesis/blob/9c9f665ed79ce0fa8da61bd803b7a4bf1ad62ed0/src/Development/AutomorphicArithmetic/DedekindLFunctions/LandauTheorems/landauWindow.lean#L125) · [两支实零点的排除](cm229-heights-and-zeros.md)

若 $\alpha_*>0$，固定 $1/2<\theta<\beta_*$、$1<d<e<3/2$，令 $L_Y=\log Y$。由该引理的逆否命题，所有充分大的 $Y$ 都有

$$
 \int_{dL_Y}^{eL_Y}[G_1f_{J_2}(w)]_+\,dw>Y^{d(\theta-1/2)}. \tag{4}
$$

事实上，若无界序列上左端不超过右端，乘上 $e^{-(\theta-1/2)w}$ 得带权积分至多 $1$。以 $dL_Y$ 为窗口左端、$e/d>1$ 为固定比，便满足窗口引理的任意小指数预算，排除了 $\Re\rho>\theta$，与 $\beta_*>\theta$ 矛盾。

精确分割给出 $f_{J_2}(w)=\sum_{j=0}^{m_Y-1}r^{-j/2} f_{J_r}(w-jv)$。写 $G_1=G_{1-a}G_a$，这些算子都保序、质量均为 $1$。在 (4) 的积分内，用正部次可加性并作换元。点 $w-jv$ 与整窗两端相距 $\gg_{d,e} L_Y$。对 $G_{1-a}$ 的卷积，窗外部分总计 $O_{d,e}(\exp(-c_{d,e}L_Y^2+O(L_Y)))$，因为 $|G_af_{J_r}(u)|\ll e^{u/2}$。窗内部分至多 $m_Y\int_{L_Y}^{3L_Y/2}[G_a f_{J_r}]_+\,du$。令 $d\uparrow3/2$、$\theta\uparrow\beta_*$ 得到下指数 $\delta$。反向由 $\sum_{j\le x}b_\nu(j)\ll_\epsilon x^{\beta_*+\epsilon}$ 得到，后式是对 $B_\nu$ 的截断 Perron 估计。$\alpha_*=0$ 时利用非负性即足。故 $\log(1+\int [G_af_{J_r}]_+)/\log Y\to\delta$，积分遍及整个窗口。再用 (2)，因 $z\le T^{1/6}$，分块积分 $[G_a f_{Q_T}]_+$ 或 $[G_af_{Q_T^\Phi}]_+$，两者均有同一极限。此处先作平均再取正部，每个可能不满的末块只积分到 $X$。

回到有限多项式。每个矩形有 $DP\asymp T$。对 $|t|\le T^2$、任意 $\epsilon>0$，有

$$
 |C_{D,\nu}(t)|\ll_\epsilon D^{\alpha_*}T^\epsilon,\qquad
 |B_{P,\nu}(t)|\ll_\epsilon P^{\alpha_*}T^\epsilon,\qquad
 |F_T(t)|\ll_\epsilon T^{\alpha_*+\epsilon}. \tag{5}
$$

具体地，在 $\Re s\ge\beta_*+\epsilon'$，已知 $L_\nu^{\pm1}$、$B_\nu$ 在多项式高度内为 $O(T^\kappa)$，可取任意 $\kappa>0$。前者由在无零半平面内对 $\log L_\nu$ 作三线插值得到，$U_{z,\nu}$ 由同一半平面内的截断 Perron 公式估计。对 $C_{z,\nu}$ 或 $B_\nu$ 从临界线的加权前缀出发，初始 $\Re w=2$、高度 $T^{10}$，把积分移至 $\Re w=\alpha_*+\epsilon'$。取半整数端点，系数分别由 $\mathbf1_{n=1}+d_4(n)$ 与 $2\Lambda(n)$ 控制，截断与水平边的误差可吸收。新竖线付 $D^{\alpha_*+\epsilon'}$ 或 $P^{\alpha_*+\epsilon'}$ 及任意小幂。后者去掉素数幂的临界线费用为 $O(\log^2(2T))$。这些结论对区间内所有截断端点一致，单调权 $\Phi_T$ 由分部求和处理。取 $\epsilon',\kappa$ 充分小并合计 $O(\log^2(2T))$ 个矩形得到 (5)，且对 $F_T^\Phi,F_T^{\le U}$ 同样成立。$|t|>T^2$ 时三个有限和均有平凡界 $O_\epsilon(T^{1/2+\epsilon})$。[竖线与前缀界的详细推导](cm229-type-ii-large-prime.md)

矩形延拓的 $f_{Q_T}$ 的 Fourier 变换为 $m_r(t)F_T(t)$，其中 $m_r(t)=2^s(1-r^{-s})/s$，$s=1/2+it$。有 $|m_r(t)|\ll\min(v,(1+|t|)^{-1})$。短余因子亦同。高斯平均再乘 $e^{-at^2/2}$。在 $c\log(2T)\le|t|\le T^2$ 的 Fourier 逆变换中，由 (5) 抽出 $e^{-a c^2\log^2(2T)/4}\le(2T)^{-c^2/4096}$，剩下的积分只付对数。$|t|>T^2$ 的尾比任意负幂更小。因此每块删去这部分所得的点态误差为 $O_{\epsilon,c}(T^{\alpha_*-c^2/4096+\epsilon})$。较低频率由 Plancherel 的二范数上界 $O({\cal E}_T^{1/2})$ 控制，在长度至多 $\log2$ 的块上也控制一范数。

当 $\delta>0$，删去部分的整窗上指数至多为 $\max(0,\tfrac32(\alpha_*-c^2/4096))<\delta$，故上述正部下界强制 $\max_T{\cal E}_T^{1/2}$ 的下指数至少为 $\delta$。上指数直接来自 (5) 和 $T\le X$。$\delta=0$ 时用同一上界即可。最后，$F_T^\Phi$ 是硬截断多项式的正权平均，加权二范数的三角不等式给出 $\|F_T^\Phi\|\le\max_{H\le U\le2H}\|F_T^{\le U}\|$。这证明 (1) 的全部三式。第三式允许端点随块变化。证明在 $\delta>0$ 时用的是严格小于 $\delta$ 的误差指数，尚不能给其余频率一个独立的任意小幂上界。

## 只在整数频率求值

对 $A=F_T$ 或 $F_T^\Phi$，在 (1) 中还可以取

$$
 N_T=\max_{0\le j\le c\log(2T),\ j\in\mathbb Z}|A(j)|.
 \tag{5a}
$$

硬截断也成立，左侧需另对 $H\le U\le2H$ 取最大。尤其对指定的短余因子权 $\Phi_T$，每个完整窗口只有 $O_c((\log Y)^2)$ 个待估的显式有限和。固定任一 $c>0$，本例 $E$ 的 ERH 等价于

$$
 \forall\eta>0,\quad\forall Y_0\ge2,\quad\exists Y\ge Y_0:\quad
 \max_{\substack{T=2^kY<X,\ k\ge0\\0\le j\le c\log(2T),\ j\in\mathbb Z}}
 |F_T^\Phi(j)|\le Y^\eta .
$$

由 (1) 和 (5a)，若这个条件成立，同样的预算也适用于所有充分远的完整窗口。

还可以逐支应用 (1) 和 (5a)。令 $\beta_\nu=\sup_{L_\nu(\rho)=0,\ 0<\Re\rho<1}\Re\rho$，在有限和里只保留这一支，仍取原权 $1/2$。此时极限为 $\delta_\nu=\tfrac32(\beta_\nu-1/2)$。证明中 (2) 的逆系数和 (5) 的前缀界均可使用本支的 $\beta_\nu$。窗口引理所用的固定函数改取本支的 $J_2$，其非实零点仍留下极点，实零点已由同一证书排除。其余估计与采样核逐支适用。因 $\zeta_F=\zeta L_+$，正支的判据恰对应全实三次域 $F$ 的 ERH。两支合并后，短余因子和的**精确指数等于 $\max(\delta_+,\delta_-)$**，因此即便在有限块上发生交叉抵消，也无法把完整窗口的指数压到较大单支指数以下。

事实上，非空矩形中任何两个乘积的比值小于 $4$。由其选取规则，所有 $q=hp$ 都在 $T/(4r)<q<32T$ 内。令 $g(t)=(2T)^{it}A(t)$，则组成 $g$ 的频率 $\xi=\log(q/(2T))$ 全满足 $|\xi|\le\Omega=\log16<\pi$。这给出了以间隔 1 采样所需的余量。以下给出从无限采样到有限采样的误差。

取 $b=1/4$，$\Omega+b<\pi$，令 $\mathrm{sinc}(u)=\sin(u)/u$，原点取值为 $1$。对任意整数 $m\ge2$ 均有恒等式

$$
 g(t)=\sum_{j\in\mathbb Z}g(j)\,\mathrm{sinc}(\pi(t-j))
               \mathrm{sinc}\big(b(t-j)/m\big)^m . \tag{5b}
$$

只须验证 $g(t)=e^{-i\xi t}$。函数 $\omega\mapsto e^{-it\omega}$ 在 $[-\pi,\pi]$ 的 Fourier 系数为 $\mathrm{sinc}(\pi(t-j))$，采用基底 $e^{-ij\omega}$。其对称 Fourier 级数在 $|\omega|\le\Omega+b$ 一致收敛于原函数，可由分部求和或 Dirichlet 核得到。令 $S$ 为 $m$ 个独立 $[-b/m,b/m]$ 均匀变量之和。在 $\omega=\xi+S$ 处代入级数，乘 $e^{itS}$ 后取期望即得 (5b)，所得核绝对可和。

对 $|t|\le K_c/2$，令 $R=K_c/2$、$m=\lfloor bR/2\rfloor$，充分大时 $m\ge2$。$|j|\le K_c$ 的核绝对值之和为 $O_c(\log\log(3T))$，由第一个 $\mathrm{sinc}$ 因子即得。$|j|>K_c$ 时令 $u=|t-j|\ge R$，核至多为 $C(1+u)^{-1}(R/(2u))^m$，总和为 $O(2^{-m})$。在 $|j|\le T^2$ 用 (5)，在 $|j|>T^2$ 用平凡界，后一段还含 $(O(R/T^2))^{m-1}$。因 $|A(-j)|=|A(j)|$，这给出

$$
 \max_{|t|\le K_c/2}|A(t)|\ll_c \log\log(3T)\,N_T
          +O_{\epsilon,c}(T^{\alpha_*-c(\log2)/32+\epsilon}).
$$

将 (1) 用于 $c/2$，频带能量的平方根至多为左端乘 $O(\sqrt{\log(2T)})$。若 $\delta>0$，误差的整窗幂指数严格较小，故 (5a) 保留下指数 $\delta$。上界由 (5) 得到，亦覆盖 $\delta=0$。对硬截断的界一致，取最大即可。这个采样结论估计的是带原符号的 $A(j)$，其中 $j$ 是实频率上的整数点，与 Hecke 零点的纵坐标是两个变量。

## 单窗的高斯平均与素数平方

使用同一归一化，可以把 (2) 中的 $O_\epsilon(T^\epsilon z^{\alpha_*})$ 再分解到 $O_\epsilon(T^\epsilon)$ 精度。令 $\mathcal Z$ 为两支非平凡零点的并多重集，$\mathcal Z_+$ 只取正支，所有和均按重数计。置 $D_z(s)=C_{z,+}(s)+C_{z,-}(s)$，$K=256\log(2T)$。则在 $T\le x\le2T$，

$$
 \begin{aligned}
 G_af_{Q_T}(\log x)={}&-\tfrac12
 \sum_{\substack{\rho\in\mathcal Z,\ \Re\rho>1/2\\ |\Im\rho|\le K}} K_x(\rho)
 -\tfrac14 D_z(1/2)K_x(1/2)\\
 &+\tfrac14\sum_{\substack{\rho\in\mathcal Z_+,\ \Re\rho>17/20\\|\Im\rho|\le 2K}}
 D_z(\rho/2)K_x(\rho/2)+O_\epsilon(T^\epsilon).
 \end{aligned}\tag{6}
$$

左端换成 $G_a f_{Q_T^\Phi}$ 也成立。第一和的重数使用两支的和。第二项来自素数平方中的平凡表示，其符号由实数 $D_z(1/2)=2-\sum_\nu U_{z,\nu}(1/2)L_\nu(1/2)$ 决定。第三项只用 $L_+$，因为二次扭转的平方在未分歧素数上为 $1$。

证明从 (2) 的级数恒等式开始。对 $A=\zeta,\chi_{229},+$，分别令 $B_A$ 表示 $\zeta,L(s,\chi_{229}),L_+$ 的负对数导数。用 $\mathrm{tr}(\rho_{\rm std}(g^2))=1-\mathrm{sgn}(g)+\mathrm{tr}(\rho_{\rm std}(g))$ 可得

$$
 B_{\nu,\mathrm{pp}}(s)=B_\zeta(2s)-B_{\chi_{229}}(2s)+B_+(2s)+E_\nu(s).
$$

其中 $E_\nu$ 在 $\Re s>1/3$ 绝对收敛，因已匹配所有未分歧素数的平方项，有限分歧素数可用各自的收敛局部级数。对 $\sigma<1/2$，对前述逆系数前缀作 Abel 分部求和，仍有 $U_{z,\nu}(\sigma+it)\ll_{\sigma,\epsilon}(2z)^{\beta_*-\sigma+\epsilon}(2+|t|)^\epsilon$。因此 $C_{z,\nu}E_\nu$ 移至 $\sigma=1/3+\epsilon'$ 后，在 (3) 中的费用为 $O_\epsilon(T^{-11/144+\epsilon})$，这里用了 $\beta_*\le7/8$。[平方迹与素数幂分解](cm229-prime-carry-drift.md) · [Hu 的无零半平面](https://hexagonmath.org/2610.00216v2)

将 $B_\nu(s)$ 和 $D_z(s)B_A(2s)$ 的各个 Mellin 积分移至 $\Re s=-1/4$。单位高度内的零点数为 $O(\log(2+|t|))$，可沿与零点纵坐标相距 $\gg1/\log(2+|t|)$ 的一列水平边移线。高斯使其积分趋零。固定负线上函数方程给出多项式增长，$U_{z,\nu}$ 的上述界使新线费用为负幂。原点的平凡项同样为负幂。$\zeta$ 在 $2s=1$ 的留数因缩放给 $+1/2$，各非平凡零点对负对数导数在 $s=\rho/2$ 给 $-m(\rho)/2$。乘原式的 $-D_z/2$ 得 (6) 的两个平方项。$B_\nu(s)$ 自身给第一和。对任一非平凡零点 $\rho=\beta+i\gamma$，在 $|\gamma|\le2K$ 有

$$
 |D_z(\rho/2)K_x(\rho/2)|
 \ll_\epsilon T^{\beta/2-1/2+(\beta_*-\beta/2)/6+\epsilon}.
$$

当 $\beta\le17/20$，主指数至多为 $0$，因为 $(5/6)(17/40)+(7/8)/6=1/2$。所有有理 $\zeta$ 与 $\chi_{229}$ 的零点均在此范围。$B_\nu(s)$ 中 $\Re\rho\le1/2$ 的项也只花 $O_\epsilon(T^\epsilon)$。两个保留和之外的高度尾用 (3) 的高斯及单位零点数估计，费用为负幂。这证明 (6)。[有理 Dirichlet 临界线定理](https://math.mit.edu/~drew/grh/) · [Hecke 显式公式与零点计数](https://aareyanmanzoor.github.io/assets/articles/lagarias-odlyzko.pdf)

例如，若两支在 $|\Im\rho|\le512\log(2X)$ 内的零点均有 $\Re\rho\le1/2+\lambda$，其中 $0<\lambda<1/20$，则第三和为空。第一和在各块为 $O_{\lambda,\epsilon}(T^{\lambda+\epsilon})$。又 $K_x(1/2)=2\sqrt2(1-r^{-1/2})>0$。于是整窗内**高斯平均后的**归一化正部积分至多为

$$
 O_{\lambda,\epsilon}\left(Y^{3\lambda/2+\epsilon}
                  +v\sum_{T=2^jY<X}[-D_z(1/2)]_+\right).
$$

该单窗估计显式保留由原截断决定的中心项。关于原整数和，(1) 给的是沿完整窗口的极限与等价判据，从一个固定窗口的这份短零点表还推不出其未平滑正部上界。

## 两支及两类余因子的有限能量

取 $Y=T=4096$ 的首块。实际 $z=1$，$K_0=T^{7/15}=48.5029\ldots$，$\log(2T)=9.0109\ldots$。记 $E_+,E_-$ 为各自含 $1/2$ 权的单支能量，$E$ 为先加两支再平方的原能量。两表各项都落在所列小数与其加 $10^{-9}$ 之间。

| $K$ | $E_+$ | $E_-$ | $E$ |
| --- | ---: | ---: | ---: |
| $\log(2T)$ | 0.005179261 | 0.003736268 | 0.008022517 |
| $T^{7/15}$ | 0.051100428 | 0.025969693 | 0.089883597 |

同一组系数在对数频带的交叉项 $E-E_+-E_-$ 约为 $-0.000893012$，在原频带则约为 $+0.012813475$。再用同一块与同一素数下限，专取 $z=4$、$4<h\le16$ 以使两类余因子同时出现。把素数、光滑两类与两支组成四个部分。$D$ 为四个单独能量之和，$E_h$ 先在每类余因子内加两支再计能量，$E_\nu$ 先在每支内加两类余因子再计能量，$E$ 对四部分全加后计算。

| $K$ | $D$ | $E_h$ | $E_\nu$ | $E$ |
| --- | ---: | ---: | ---: | ---: |
| $\log(2T)$ | 0.005084165 | 0.003938653 | 0.003327013 | 0.003375333 |
| $T^{7/15}$ | 0.028500458 | 0.025665250 | 0.030764813 | 0.028801352 |

第二表的对数频带内 $E<D$，扩至原频带后 $E>D$。这两个有限检验排除了相应平方展开中交叉项恒为非正的简单界。原截断的长期结论由 (1) 给出。

在 (5a) 取 $c=1$ 时，首块只用 $j=0,\ldots,9$。第一表的原多项式在 $j=9$ 取到最大模 $5.966869266\ldots$，第二表的有限截断在 $j=0$ 取到最大模 $2.791384673\ldots$。对 $t=3/7$ 独立检验 (5b)，取 $m=14$ 和 $|j|\le160$。用其余整数上 $|g(j)|\le M$ 及核的幂衰减，插值尾上界分别小于 $6.5\times10^{-6}$ 和 $2.1\times10^{-6}$，两次直接求和与插值的差均落在该界内。

数值用 160 位球算术。此块内权恒为 $v^2$。在 $[0,K]$ 分段，每段长 $\Delta\le2$，用 20 点 Gauss–Legendre 规则，再以偶性计及负频率。若 $M$ 为各部分系数绝对质量之和的上界，$\Omega=\log(q_{\max}/q_{\min})$，则能量是频率不超过 $\Omega$ 的有限指数和的积分。对指数在每段中点展开至 39 次，用 40 阶导数界及积分和正权求积各自的余项，全部截断误差至多 $4Kv^2M^2(\Omega\Delta/2)^{40}/40!$。四组计算中该界均小于 $8\times10^{-25}$。程序还在 $t=0,3/7$ 直接按乘积模四投影与分支乘法作独立对照。[区间复算程序](verify_cm229_low_band.py) · [零点处两个余因子分量的独立复算](cm229-signed-zero-modes.md)
