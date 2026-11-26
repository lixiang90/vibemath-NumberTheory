---
title: "在原 WinPair 同余条件下用双惰性素数恢复固定数域的素数和"
date: "2026-11-26"
abstract: "在原 Lean 证明规定的两种指数尺度上选择辅助素数。双惰性正权平均以任意负幂精度恢复固定素数集合的有限和，几乎每个单独宿主对也有平方根级转移误差。结合原窗口引理，得到任意固定数域按一次素理想计数的单侧稀疏窗口判据。在判别式 229 的六次 CM 域，只需对一种素数给出主项为 x/48 的上界。"
---

# 在原 WinPair 同余条件下用双惰性素数恢复固定数域的素数和

固定一个有理素数集合 $S$，写 $\Theta_S(x)=\sum_{p\le x,\,p\in S}\log p$，约定 $2\notin S$。两个虚二次宿主 $\mathbb Q(\sqrt{-q}),\mathbb Q(\sqrt{-d})$ 给出非负的计数

$$
 I_{q,d}(x)=
 \sum_{\substack{p\le x,\ p\in S\\(-q/p)=(-d/p)=-1}}\log p. \tag{1}
$$

下面证明，在原证明指定的同余和大小条件下，可以从 (1) 恢复固定的 $\Theta_S(x)$。平均恢复时，须加回一个完全明确的小素数项。选单独一对宿主时，几乎每一对都能在整个幂次窗口上同时达到平方根级的精度。这个结论对 $S$ 的分布没有假设。固定域的 ERH 最终归结为几何一侧对 (1) 本身的单侧估计。

## 选择满足全部 WinPair 条件的两个素数

取 $0<e<1/52$、固定 $A>1$，令 $Y$ 趋向无穷，置

$$
 z=Y^{e/5},\quad P=\prod_{3\le p\le z}p,\quad
 U=\exp(3z/4),\qquad D=\exp\!\big(3z(\log Y)^2/4\big).
$$

这里及下文乘积中的 $p$ 都是素数。先从 $(U,2U]$ 中按权 $\log q$ 选取 $q\equiv3\bmod 8$ 的素数。给定 $q$，再从 $(D,2D]$ 中按权 $\log d$ 选取满足 $d\equiv7\bmod16$、$d\equiv-1\bmod qP$ 的素数。每一步将权除以本步的总权。以 $\mathbb E$ 和 $\mathbb P$ 表示这种先选 $q$、再选 $d$ 的期望与概率。两族非空性见下文证明。对充分大的 $Y$，所选素数都大于 $Y^A$。

这些素数对逐条满足原 Lean 定义的 $\operatorname{WinPair}(e,Y,q,d)$，包括 $q\mid d+1$ 和每个奇素数 $p\le z$ 都整除 $d+1$。两个素数的对数分别在 $(z/2,z]$ 和 $(z(\log Y)^2/2,z(\log Y)^2]$ 内。若硬把两个奇素数放在同一整数二进区间 $(Q,2Q]$，则 $q\mid d+1$ 将迫使 $2q\le d+1\le 2Q+1$，已经矛盾。这里采用两个尺度正是为了保留原有同余条件。[原始 WinPair 与 VWinAt 的定义](https://github.com/leanprover-community/riemann-hypothesis/blob/9c9f665ed79ce0fa8da61bd803b7a4bf1ad62ed0/src/Development/AutomorphicArithmetic/CharacterBoxes/ExistsPageN.lean#L70) · [此前固定全实基域上的理想窗口](real-cm-prime-windows.md)

**转移定理。** 对所有 $Y\le x\le Y^A$ 和所有预先指定的 $S$，有以下两项结论。

1. 对每个 $B>0$，一致地有
   $$
   \Theta_S(x)=\Theta_S(z)+4\mathbb E I_{q,d}(x)+O_{e,A,B}(Y^{-B}). \tag{2}
   $$
2. 给定 $\delta>0$，除去总概率 $O_{e,A,\delta}(Y^{-2\delta}\log^4(2Y))$ 的素数对，其余每对都同时满足
   $$
    |4 I_{q,d}(x)-\Theta_S(x)+\Theta_S(z)|\le x^{1/2+\delta}
       \qquad(Y\le x\le Y^A). \tag{3}
   $$

各常数对 $S$ 一致，(3) 的例外集可以依赖于 $S$。任意预先指定的有限组 $S$ 可对例外集取并。式 (2) 也适用于随 $Y$ 给定、与所抽取的 $q,d$ 无关的集合。两式用到已核验的有理 Dirichlet GRH。[原 Dirichlet GRH 的独立核验](https://math.mit.edu/~drew/grh/)

**证明。** 对 $z<p\le Y^A$ 写 $u_p=(-q/p)$、$v_p=(-d/p)$。对这一区间内由一或两个互异素数组成的非空集合 $T$，令 $m_T=\prod_{p\in T}p$。则

$$
 \left|\mathbb E_q\prod_{p\in T}u_p\right|
 +\sup_q\left|\mathbb E_{d\mid q}\prod_{p\in T}v_p\right|
 \ll_{e,A} e^{-c z} \tag{4}
$$

对某个绝对常数 $c>0$ 和充分大的 $Y$ 成立。事实上，除一个固定符号外，两个乘积都是变量 $q$ 或 $d$ 的模 $m_T$ 非主二次特征。对 $q$ 的同余类模 8 展开特征，对 $d$ 的同余类模 $h_q=16qP$ 展开特征。由于 $m_T$ 与这些模数互素，乘上来自 $m_T$ 的特征后每一项仍非主。Dirichlet GRH 的显式公式给出模 $m$ 的非主特征在 $(t,2t]$ 上的素数对数权和为 $O(\sqrt t\log^2(2tm))$。这条标准界由截断显式公式、单位高度零点数 $O(\log(m(2+|t'|)))$ 以及去除素数高次幂得到，也适用于由较小导子诱导的特征。正交展开分母为 $\varphi(h_q)$，误差没有再乘 $\varphi(h_q)$。特别地，两个抽样分母分别为

$$
 U/4+O(\sqrt U\log^2(2U)),\qquad
 D/\varphi(h_q)+O(\sqrt D\log^2(2Dh_q)).
$$

因为 $\log h_q=O(z)$、$m_T\le Y^{2A}$，分子误差除以分母后，对 $q$ 至多为 $O(U^{-1/2}\log^2(2UY^{2A}))$，对 $d$ 至多为 $O(\varphi(h_q)D^{-1/2}\log^2(2Dh_qY^{2A}))$。前者与后者均满足 (4)，也验证了分母为正。[Montgomery–Vaughan 关于各模特征的 GRH 估计，第 17 章](https://personal.science.psu.edu/rcv4/personal/Publications/MNTI/17.0_pp_419_451_Conditional_estimates.pdf) · [Ernvall-Hytönen–Palojärvi 的显式界](https://arxiv.org/abs/2003.01925)

取 $X_p=(1-u_p)(1-v_p)=4\mathbf1_{\{u_p=v_p=-1\}}$。先对 $d$ 再对 $q$ 求期望，(4) 给出 $\mathbb E X_p=1+O(e^{-cz})$。对互异的 $p,t$ 展开乘积，同样给出 $\mathbb E[(X_p-1)(X_t-1)]=O(e^{-cz})$。这是因为关于任一变量的每个非恒定项都是 (4) 中一阶或二阶的乘积。对 $p=t$ 直接有 $|X_p-1|\le3$。另一方面，$3\le p\le z$ 时 $-d\equiv1\bmod p$，所以 $I_{q,d}$ 确实遗漏了此范围内的全部素数。用 $\sum_{p\le x}\log p\ll x$ 相加，误差 $O(xe^{-cz})$ 对 $x\le Y^A$ 小于每个 $Y$ 的固定负幂，得到 (2)。

为证 (3)，对 $z<p\le 2V$（同时截于 $Y^A$）的任意系数 $b_p$，二阶估计给出

$$
 \mathbb E\left|\sum_p b_p(X_p-1)\right|^2
 \le 9\sum_p |b_p|^2 + O\!\left(e^{-cz}\Big(\sum_p|b_p|\Big)^2\right).
$$

按素数大小排序，用二进区块分解所有前缀。每个前缀在每层最多用一个区块，Cauchy 不等式给出一个对数因子。每层各区块互不相交，再对层数求和给出第二个对数因子。取 $b_p=\mathbf1_S(p)\log p$，由上式可得前缀最大绝对值的均方 $\ll V\log^4(2V)$，指数小的误差可一致吸收。对 $V=2^jY\le Y^A$ 用 Markov 不等式并加总几何级数，便得 (3)。

## 任意有限正规数域的稀疏窗口判据

设 $M/\mathbb Q$ 是固定的有限伽罗瓦扩张，$n=[M:\mathbb Q]$。在 (1) 取 $S$ 为在 $M$ 完全分裂的未分歧奇素数。固定 $A=3/2$ 及一个 $0<e<1/52$。则 $\operatorname{ERH}(M)$ **当且仅当**对每个 $\delta>0$，都有一个常数 $C_\delta$ 及任意大的 $Y$，使

$$
 \mathbb E I_{q,d}(x)\ \ge\ \frac{x}{4n}-C_\delta x^{1/2+\delta}
     \qquad(Y\le x\le Y^{3/2}). \tag{5}
$$

也可以删去期望，要求至少十分之一总权的宿主对各自在整窗满足 (5) 的下界，用同一个 $C_\delta$。在 ERH 下，这种宿主对的总权实际上趋于 1。

为证明逆向，记 $\Psi_M(x)=\sum_{N\mathfrak a\le x}\Lambda_M(\mathfrak a)$。由于扩张正规，未分歧的一次素理想只来自完全分裂素数，每个贡献 $n$ 次。其余素理想幂以及分歧素数总共给出

$$
 \Psi_M(x)=n\Theta_S(x)+O_M\!\big(\sqrt x\log^2(2x)\big). \tag{6}
$$

式 (2) 将 (5) 转成固定的 $\Theta_S$ 的同精度下界。若已知十分之一总权上的下界，则 (3) 保证其中至少有一对也满足转移估计。反向在 ERH 下，由素理想定理 $\Psi_M(x)=x+O_M(\sqrt{x}\log^2(2x))$ 得到两种宿主表述。此处 $\Theta_S(z)\ll z=o(\sqrt Y)$。[Grenié–Molteni 的素理想定理](https://arxiv.org/abs/1312.4463)

最后说明为何可以任意稀疏地选择 $Y$。取可测实函数 $g(u)=e^u-\Psi_M(e^u)$。它在 $u\ge0$ 上是 $O_M(e^{2u})$。对 $\Re s>1$ 的 Laplace 变换是

$$
 \frac1{s-1}+\frac{\zeta_M'(s)}{s\zeta_M(s)}.
$$

$s=1$ 处的极点相消。右半平面其余各极点恰是 $\zeta_M$ 的零点 $\rho$，若重数为 $m$，留数为 $m/\rho$，实极点的留数严格为正。这逐项验证了原抽象窗口引理的极点条件。由 (5)、(6)，对任意小的 $\delta>0$，在任意远的 $[V,3V/2]$ 上有 $g(u)\le C'_\delta e^{(1/2+\delta)u}$，其中 $V=\log Y$。取此右端为非负上界，乘 $e^{-u/2}$ 后的积分至多为 $O_\delta(Ve^{3\delta V/2})$。给定 $\eta>0$，先取 $3\delta/2<\eta$，再取足够大的 $V$。这正是原引理在 $\theta=1/2$ 的窗口条件。因而每个零点 $\Re\rho\le1/2$，函数方程给出 ERH。[原 Lean 极点与窗口条件](https://github.com/leanprover-community/riemann-hypothesis/blob/9c9f665ed79ce0fa8da61bd803b7a4bf1ad62ed0/src/Development/AutomorphicArithmetic/DedekindLFunctions/LandauZeroFreeCriterion.lean#L57) · [原始 landau_window 定理](https://github.com/leanprover-community/riemann-hypothesis/blob/9c9f665ed79ce0fa8da61bd803b7a4bf1ad62ed0/src/Development/AutomorphicArithmetic/DedekindLFunctions/LandauTheorems/landauWindow.lean#L125)

**任意固定数域的推论。** 设 $K/\mathbb Q$ 是数域，$r_K(p)$ 是 $K$ 中范数等于 $p$ 的素理想数。在奇的未分歧素数上定义

$$
 J_{K;q,d}(x)=\sum_{\substack{p\le x,\ p\nmid 2D_K\\(-q/p)=(-d/p)=-1}}
                r_K(p)\log p .
$$

则 $\operatorname{ERH}(K)$ 等价于按照 (5) 的窗口与误差量词成立 $\mathbb E J_{K;q,d}(x)\ge x/4-C_{K,\delta}x^{1/2+\delta}$。十分之一总权的版本同样成立。这包括每个指定的非正规全实域与非正规 CM 域。证明只需对至多 $[K:\mathbb Q]$ 个固定集合 $S_j=\{p\nmid 2D_K:r_K(p)\ge j\}$ 应用 (2)、(3) 并相加。由于剩余次数或幂指数大于一的素理想幂只贡献 $O_K(\sqrt{x}\log^2(2x))$，这些固定素数和之和即在此误差内给出 $\Psi_K(x)$。对 $\zeta_K$ 重复刚才的 Laplace 变换证明即可。在非正规域中，一个一次素理想即可有贡献，$p$ 在 $K$ 完全分裂时才有 $r_K(p)=[K:\mathbb Q]$。

## 判别式 229 的六次 CM 域只需一个正权上界

取 $F=\mathbb Q(\alpha)$，$f(\alpha)=\alpha^3-4\alpha+1=0$，$E=F(i)$，$N$ 为 $F$ 的正规闭包，$M=N(i)$。此时 $n=12$，$S$ 的成员恰是满足 $p\nmid 2\cdot229$、$p\equiv1\bmod4$ 且 $f$ 在 $\mathbb F_p$ 上有三个根的素数。对这个 $S$，$\operatorname{ERH}(E)$ 与 (5) 等价，而且也等价于把 (5) 整条换成**上界**

$$
 \mathbb E I_{q,d}(x)\ \le\ \frac{x}{48}+C_\delta x^{1/2+\delta}
 \qquad(Y\le x\le Y^{3/2}). \tag{7}
$$

上界同样可以只要求十分之一总权的宿主对在整窗满足。这里的目标是 $E$ 本身，其中同时包含全实因子 $L_+$ 及其固定奇扭转 $L_-$。由前篇的 $\zeta_E=\zeta L(\chi_{-4})L_+ L_-$、$\zeta_M=\zeta_{\mathbb Q(\sqrt{229},i)}(L_+ L_-)^2$，以及 Dirichlet GRH，两个数域的 ERH 等价。[两域分解及两个因子的实零点排除](cm229-heights-and-zeros.md)

为得到上界方向，令 $L=L_+L_-$，$-L'/L=\sum b(m)m^{-s}$，$T(x)=\sum_{m\le x}b(m)$，并记 $S_0$ 为在 $\mathbb Q(\sqrt{229},i)$ 完全分裂的未分歧奇素数。$S_3$ 的标准迹在单位元、换位、三循环上分别是 $2,0,-1$，等于 $3\mathbf1_{g=1}-\mathbf1_{\operatorname{sgn}g=1}$。再乘 $1+\chi_{-4}(p)$，得到一次项 $b(p)=(6\mathbf1_S(p)-2\mathbf1_{S_0}(p))\log p$。有限个分歧处及全部高次素数幂的代价为 $O(\sqrt x\log^2(2x))$。已知阿贝尔域的 ERH 给出 $\Theta_{S_0}(x)=x/4+O(\sqrt x\log^2(2x))$。故

$$
 T(x)=6\big(\Theta_S(x)-x/12\big)+O(\sqrt x\log^2(2x)).
$$

取 $g(u)=T(e^u)$，其 Laplace 变换为 $-L'(s)/(sL(s))$。实部为正的零点 $\rho$ 给出留数 $-m/\rho$。前篇的完整矩形计数已排除两个因子的全部实非平凡零点，$\Re s\ge1$ 由 Hecke 非消失性。因此这里没有需要正留数的实极点，原窗口引理同样适用于这个符号。应用 (2) 或 (3) 就得到 (7) 的逆向。这解释了为何可以把此前两个模 4 素数和的双侧判据压缩成一个正权稀疏窗口上界。

目前已有的全高度无零区可以产生一个实际的误差指数。对任意 $\lambda>0$，在每个充分远的上述窗口上有

$$
 \mathbb E I_{q,d}(x)=\frac{x}{48}+O_\lambda(x^{31/32+\lambda}).
$$

事实上，Cureton–Thorner 给出固定正规域 $M$ 的 $\zeta_M$ 在 $\Re s>31/32$ 无零。在 $31/32+\lambda/2\le\Re s\le2$ 远离极点 1 时，对数导数由局部零点计数及 Hadamard 乘积界为 $O_{M,\lambda}(\log(2+|\Im s|))$，此处只需考虑小的 $\lambda>0$。在半整数 $y$ 对 $\Psi_M(y)$ 用高度 $y^2$ 的截断 Perron 公式，移至上述左边界，跨过 $s=1$ 的留数 $y$，余项为 $O_{M,\lambda}(y^{31/32+\lambda/2}\log^2(2y))$。系数在 $n$ 处的绝对值至多为 $12\Lambda(n)$，取半整数消除了截断点的跳跃问题，再由单调性还原所有实 $x$。结合 (2)、(6) 即得所示估计。由 (3)，同一误差指数也适用于总权趋于 1 的单独宿主对。[Cureton–Thorner 的全高度无零半平面](https://hexagonmath.org/2610.00200v2)

每个计入 (7) 的 $p$ 在 $F$ 完全分裂，在两个虚二次宿主都惰性。于是对于 $F(\sqrt{-q})$ 和 $F(\sqrt{-d})$，各自的八个 CM 型在这样的 $p$ 处都具有纯 $1/2$ 的 Newton 斜率。相应的局部整模型还需处理切空间作用和 Taylor 系数。伴篇给出这一步的具体障碍。[CM 型、Newton 斜率与整系数标架](cm229-frobenius-and-heights.md)

## 惰性一次素理想进入原斜坡和的校正项

设 $K$ 是固定数域，$\eta$ 是 $K$ 的本原二次特征，$t=N\mathfrak p$。写 $R_\eta,R_1$ 为前篇定义的斜坡和，$J_\eta(x)=2\sum_{\eta(\mathfrak p)=-1}\log t\,(x/t-1)_+$。在分歧处令 $\eta(\mathfrak p)=0$。直接从全素理想幂过渡到一次惰性素理想时，有如下统一公式

$$
 \begin{aligned}
 R_1(x)-R_\eta(x)&=J_\eta(x)+xB_\eta-D_\eta(x),\\
 B_\eta&=2\sum_{\eta(\mathfrak p)=-1}\frac{\log t}{t^3(1-t^{-2})}
       +\sum_{\mathfrak p\mid\mathfrak f_\eta}\frac{\log t}{t-1},\\
 0\le D_\eta(x)&\ll_K x^{1/3}\log^2(2x)
                  +\log(2N\mathfrak f_\eta)\log(2x).
 \end{aligned} \tag{8}
$$

在惰性处，$1-\eta(\mathfrak p)^j$ 对奇数 $j$ 为 2，对偶数 $j$ 为 0。在分歧处每个 $j\ge1$ 都给出 1。将 (8) 第一行未计入 $J_\eta$ 的每个斜坡写作 $x/t^j-\min(1,x/t^j)$，线性部分就是绝对收敛的 $xB_\eta$。其余正项之和就是 $D_\eta$。高奇次幂 $t^j\le y,\ j\ge3$ 的总对数权至多为 $O_K(y^{1/3}\log^2(2y))$，可将 $\mathfrak p$ 按下方有理素数分组求和。若这部分的累积权为 $H(y)$，则它对 $D_\eta(x)$ 的贡献精确等于 $x\int_x^\infty H(y)y^{-2}dy$。每个分歧素理想的误差由几何级数界为 $O(\log(2x)+\log t)$，得到 (8)。

对**固定**的 $\eta$ 还能确定首个低阶项。二次扩张的素理想定理给出 $2\sum_{t\le y^{1/3},\,\eta(\mathfrak p)=-1}\log t\sim y^{1/3}$。$j\ge5$ 的总权为 $O_K(y^{1/5}\log^2(2y))$。代入上面的积分可得 $D_\eta(x)\sim \tfrac32 x^{1/3}$。对随 $Y$ 移动的宿主使用 (8) 的统一上界。特别地，将本文的有理二次特征拉回固定的 $K$ 时，$\log(2N\mathfrak f_\eta)\ll_K z(\log Y)^2$，故在整窗 $D_\eta(x)=o(\sqrt x)$。

若记 $\gamma_K$ 为 $\zeta_K'/\zeta_K$ 在 1 处的 Laurent 常数项，$\ell_\eta=L_K'(1,\eta)/L_K(1,\eta)$，则规范化后的差为

$$
 G_K(x)-G_\eta(x)=J_\eta(x)-x\log x+
          (1+\gamma_K-\ell_\eta+B_\eta)x-D_\eta(x).
$$

这里 $G_K=R_1-x\log x+(1+\gamma_K)x$。这给出了将惰性一次项送入几何斜率不等式时必须匹配的线性系数。若再乘固定特征 $\omega$，需要保留原始乘积与本原化的局部差。例如在 $\omega=\eta$ 时，分歧处逐点相乘给出零，本原化后的平凡特征在同处给出一。[前篇的斜坡归一化](real-cm-prime-windows.md) · [Chebotarev 密度与素理想定理](https://www.jmilne.org/math/CourseNotes/CFT.pdf)

## 仍需几何给出的估计

有限多个宿主的惰性指示函数无法逐点覆盖所有素数。完全分裂于它们合成域的素数使所有这些指示函数同时为零。例如 $q=11,d=263,p=37$，有 $q\mid d+1$，且 $10^2\equiv-11,\ 12^2\equiv-263\bmod37$。三次多项式在 37 的根为 $8,13,16$，所以这个 $p$ 恰在固定目标 $S$ 中。正确的四象限指示函数是 $(1\pm u_p)(1\pm v_p)/4$。若固定一对与 $2D_M$ 互素的不同辅助素数，每个非平凡二次子域都在 $q$ 或 $d$ 分歧，故宿主的双二次合成域与 $M$ 线性无交。由 Chebotarev 定理，在 $M$ 完全分裂的四个象限各占全体素数的 $1/(4n)$。式 (2) 和 (3) 进一步处理了宿主随窗口移动时的误差。

原代码的 $\operatorname{PairGood}$ 还要求对一族扭转的对数导数和零点盒作统一控制，$\operatorname{VWinAt}$ 对相应数域与子窗口另有量词。固定同余类的素数定理本身不能满足这些条件，$Y^{3/2}/\exp(Y^e)\to0$ 也表明仅有 $\exp(Y^e)$ 量级的例外数上界无法保证多项式尺度上的宿主存在。本文已经处理了真实 $\operatorname{WinPair}$ 中增长的同余模。要推出全高度 ERH，仍须从几何侧给出 (5) 或 (7) 中含有正确主项的单侧界。式 (2) 和 (3) 对任何 $S$ 都成立，本身没有为这个固定 $S$ 提供相对 $x/48$ 的抵消。[原 PairGood 的完整条件](https://github.com/leanprover-community/riemann-hypothesis/blob/9c9f665ed79ce0fa8da61bd803b7a4bf1ad62ed0/src/Development/AutomorphicArithmetic/CharacterBoxes/PairGood.lean#L55)

对 (7) 的后续研究已证明短区间上的严格上界，给出一亿以内固定域计数的逐点证书，并将原幂次窗口上的充分条件降到超额正部的一个加权积分。[主项为 x/48 的素数上界能在哪些区间证明？](cm229-one-sided-primes.md)
