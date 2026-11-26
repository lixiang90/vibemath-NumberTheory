---
title: "判别式 229 的最优格映射在临界线上的对数渐近式"
date: "2026-11-27"
abstract: "在固定六次 CM 域的 ERH 下，证明沿所有完整幂次窗口，最小奇异值代价除以 W log Y 收敛到同一个正常数。以零点留数的有限三角和逼近该常数，在已隔离的最低零点上施加 80 次 Fejér 核，得到大于 0.61070 的下界。结合临界线外零点的窗口定理，该对数规模的下界无条件成立。"
---

# 判别式 229 的最优格映射在临界线上的对数渐近式

设 $F=\mathbb Q(\alpha)$，$\alpha^3-4\alpha+1=0$，$E=F(i)$，$M=F^{\rm gal}(i)$，$B=\mathbb Q(\sqrt{229},i)$。写 $\zeta_M=\zeta_B L^2$，其中 $L=L_+L_-$ 是前篇给出的两个二维自对偶 Artin 因子之积。取 $W=\lceil Y^3\rceil$，$C_*(Y)$ 为从 $[Y,Y^{3/2}]$ 的整值多项式格到标准格的最小正奇异值代价，使用未除以 $[B:\mathbb Q]$ 的度数。伴篇已从混合 CM 型构造了达到此值的对角整同构。[整同构、度量与代价的定义](cm229-optimal-splitting-map.md) · [两个因子及其低零点](cm229-heights-and-zeros.md)

**定理。** 若 $E$ 满足 ERH，则存在常数 $c_{229}$，使沿全部实数 $Y\to\infty$ 有

$$
 \frac{C_*(Y)}{W\log Y}\longrightarrow c_{229},
 \qquad c_{229}>0.61070. \tag{1}
$$

若 $E$ 有临界线外的非平凡零点，则此比值趋于 $+\infty$。因而在两种情形下，每个充分远的完整窗口均满足 $C_*(Y)>0.61070\,W\log Y$。本文使用已核验的有理 Dirichlet GRH，关于 $E$ 的 ERH 只用于推出 (1) 的有限极限。[有理 Dirichlet GRH 的核验](https://math.mit.edu/~drew/grh/) · [有临界线外零点时的完整窗口指数](cm229-binomial-slopes.md)

## 从奇异值到临界线上的平均

设 $-L'/L(s)=\sum_{m\ge1} b(m)m^{-s}$，$D(x)=\sum_{m\le x}\sum_{d\mid m}b(d)$。记 $H(x)=D(2x)-2D(x)$，$G(u)=e^{-u/2}H(e^u)$，并写 $U=\log Y$。$H(n)=2\log r_n$，其中 $r_n$ 是伴篇的有理生成元。由完整窗口的离散化和 $k_n=\lfloor W/n^{3/2}\rfloor$ 的取整界，

$$
 \frac{C_*(Y)}{2W}
     =\int_U^{3U/2}[G(u)]_+\,du+O(1). \tag{2}
$$

具体地，在 $[n,n+1)$ 的后半段 $H$ 比 $H(n)$ 多出 $\sum_{d\mid2n+1}b(d)=O(\log(2n))$，故以离散和替换连续正部积分花费 $O(Y^{-1/2}\log(2Y))$。整数权的取整误差至多 $W^{-1}\sum_{n<Y^{3/2}}|H(n)|=O(1)$，这里 $H(n)=O(n)$。因此 (1) 中的比值与 $[G]_+$ 在 $[U,3U/2]$ 上的平均之差为 $O(1/U)$。[前篇公式 (2)–(5)](cm229-binomial-slopes.md)

对 $\Re s>1$ 有

$$
 \Phi(s)=\int_1^\infty H(x)x^{-s-1}\,dx
       =\frac{(2^s-2)\zeta(s)}s\left(-\frac{L'}L(s)\right). \tag{3}
$$

乘子在 $s=1$ 可去。开始假设 $\operatorname{ERH}(E)$。此时 $L$ 的非平凡零点都在临界线上，前篇的 Perron 估计给出 $H(x)\ll_\varepsilon x^{1/2+\varepsilon}$。故 (3) 的积分对 $\Re s>1/2$ 绝对收敛，解析延拓给出同一个 $\Phi$。每个不同的零点 $\rho=1/2+i\gamma$ 定义一个系数

$$
 A_\gamma=-m_\rho\frac{(2^\rho-2)\zeta(\rho)}{\rho}, \tag{4}
$$

其中 $m_\rho$ 为它在乘积 $L$ 中的总重数。若 $\zeta(\rho)=0$，置 $A_\gamma=0$ 正好对应 (3) 在该处的可去奇点。低零点证书排除了 $L(1/2)=0$，因此没有频率 0，且 $A_{-\gamma}=\overline{A_\gamma}$。[零点隔离及实轴的论证](cm229-heights-and-zeros.md)

## 临界线附近的均方尾界

令 $\sigma=1/2+\tau$，$0<\tau\le1/4$。对每个整数 $m$，有

$$
 \int_m^{m+1}\left|\frac{L'}{L}(\sigma+it)\right|^2dt
       \ll \tau^{-1}\log^2(2+|m|). \tag{5}
$$

事实上，在 $\sigma+it$ 和 $2+it$ 相减完备函数的 Hadamard 对数导数，结合 Gamma 因子的 Stirling 估计及右侧欧拉乘积，得到 $L'/L(\sigma+it)=\sum_{|t-\gamma|\le1}m_\rho/(\sigma+it-\rho)+O(\log(2+|t|))$。远处零点的差由 $(1+|t-\gamma|^2)^{-1}$ 控制，单位区间的零点总重数为 $O(\log(2+|t|))$。对近处和用 Cauchy–Schwarz，再对每个重数项积分 $(\tau^2+(t-\gamma)^2)^{-1}$，上界为 $\pi/\tau$，从而得到 (5)。这步允许多个零点具有相同或很接近的纵坐标。

在 $1/2\le\sigma\le3/4$、$|t|\ge2$，Riemann zeta 的经典凸性界给出 $|\zeta(\sigma+it)|\ll(1+|t|)^{1/4+1/16}$。结合 (3)、(5)，单位区间对 $\tau\int|\Phi|^2$ 的贡献至多为 $O((1+|m|)^{-11/8}\log^2(2+|m|))$。因此对于 $T\ge3$，

$$
 \tau\int_{|t|\ge T}|\Phi(1/2+\tau+it)|^2dt\ll T^{-1/4}. \tag{6}
$$

这里的常数只依赖固定的两个因子。指数 $1/4$ 足以使高频尾在取长窗口极限后趋于零。[Titchmarsh，Riemann zeta 的临界带增长](https://sites.math.rutgers.edu/~zeilberg/EM18/TitchmarshZeta.pdf)

## 有限留数和给出同一个极限

取 $T$ 避开零点的纵坐标，按不同频率组成实三角多项式

$$
 P_T(u)=\sum_{|\gamma|<T} A_\gamma e^{i\gamma u}.
$$

在 $u\ge0$ 对 $e^{-\tau u}(G-P_T)$ 使用 Plancherel 等式。由 (3)，其平方积分等于

$$
 \frac1{2\pi}\int_{\mathbb R}
 \left|\Phi(1/2+\tau+it)-
     \sum_{|\gamma|<T}\frac{A_\gamma}{\tau+i(t-\gamma)}\right|^2dt . \tag{7}
$$

对固定 $T$，在 $|t|\le T$ 减去极部后余项一致有界。在 $|t|>T$，有限个有理项的平方积分随 $\tau\downarrow0$ 保持有界，因为截断点与每个保留的极点有正距离。将 (7) 乘 $\tau$，用 (6) 控制余尾，令 $\tau=1/U$。在 $[U,3U/2]$ 上有 $e^{-2\tau u}\ge e^{-3}$，故

$$
 \limsup_{U\to\infty}\frac2U\int_U^{3U/2}|G(u)-P_T(u)|^2du
       \ll T^{-1/4}. \tag{8}
$$

对任何有限 $T$，$[P_T]_+$ 在这些区间上的均值有极限 $\mu_T$。确实，它是有限多个相位 $e^{i\gamma u}$ 的连续函数，在相应环面上用三角多项式一致逼近后，每个非零频率的区间平均均趋于 0。此论证允许频率间有有理关系。由于正部是 1-Lipschitz，由 (8) 及 Cauchy–Schwarz，$[G]_+$ 的区间平均与 $\mu_T$ 之差的绝对值的上极限为 $O(T^{-1/8})$。令 $T\to\infty$ 知 $\mu_T$ 为 Cauchy 族，且

$$
 c_{229}=\lim_{T\to\infty}\mu_T
        =\lim_{U\to\infty}\frac2U\int_U^{3U/2}[G(u)]_+du . \tag{9}
$$

结合 (2) 即得 (1) 中的极限。该推导使用 (3) 的具体乘子，在临界线上保留了 $\zeta(\rho)$ 可能造成的抵消。一般素数误差项的均方几乎周期方法可参见 Akbary、Ng 和 Shahabi。[素数误差的极限分布与均方逼近](https://arxiv.org/abs/1306.1657)

## 一个可核验的正常数下界

每个 $P_T$ 的平均为 0。由 (8)，$G$ 的平均也趋于 0，因此 $|G|$ 的平均趋于 $2c_{229}$。固定一个零点频率 $\gamma_0$。一旦 $T>|\gamma_0|$，$P_T(u)e^{-i\gamma_0u}$ 的平均趋于 $A_{\gamma_0}$。再次使用 (8) 并令 $T\to\infty$，有

$$
 2c_{229}\ge
 \left|\lim_{U\to\infty}\frac2U\int_U^{3U/2}G(u)e^{-i\gamma_0u}du\right|
 =|A_{\gamma_0}|.
$$

取两个因子合在一起的第一个正高度，它属于 $L_-$，是 $L$ 的单零点，且

$$
 \frac{13347347}{16777216}<\gamma_0<\frac{13347348}{16777216}.
$$

对整段在 (4) 中作复球求值，得到 $0.48564<|A_{\gamma_0}|/2<0.48565$。还可以利用该零点的整数倍。对 $2\le j\le80$，$j\gamma_0$ 仍在高度 64 以下，与全部 192 个正零点区间均不相交。将 $\gamma_0$ 的区间乘 $j$ 再逐一比较端点，最小间隔为 $13542/2^{24}$。由 (8)，$G$ 在这些频率上的平均傅里叶系数全部为零。[零点 CSV](cm229-zeros64.csv) · [绕数及单性核验](verify_cm229.py)

用次数为 80、积分归一为 1 的 Fejér 核卷积 $\operatorname{sgn}(\cos v)$，得到实三角多项式 $Q(v)$。正核保证 $|Q|\le1$，一次余弦的系数为 $(4/\pi)(80/81)$。取 $\varphi=\arg A_{\gamma_0}$，在 $G(u)Q(\gamma_0 u+\varphi)$ 的区间平均中，其余谐波的贡献均为零。因为 $|G|\ge GQ$，遂有

$$
 c_{229}\ \ge\ \frac{2}{\pi}\frac{80}{81}|A_{\gamma_0}|
        \ >\ 0.61070 .
$$

程序用整数检验全部谐波的间隔，并将上式中间的量严格包在 $(0.61070,0.61072)$ 内。这个区间给出 $c_{229}$ 的下界，$c_{229}$ 本身的大小还取决于其余零点。[留数和 Fejér 下界的区间程序](verify_cm229_optimal_map.py)

最后，若 $\operatorname{ERH}(E)$ 不成立，由函数方程和有理 Dirichlet GRH，$\beta_*=\sup_{L(\rho)=0,\,0<\Re\rho<1}\Re\rho>1/2$。完整窗口定理给出 $\log(1+C_*(Y)/W)/\log Y\to\tfrac32(\beta_*-1/2)>0$，故 $C_*(Y)/(W\log Y)\to+\infty$。于是 $WY^{o(1)}$ 在任意远窗口上可达，当且仅当 (1) 在全部充分远窗口上成立。现在仍须独立控制伴篇公式 (5) 中分裂类与惰性类的大素数带符号和，才能判定本例属于哪种情形。[完整窗口的零点判据](cm229-binomial-slopes.md) · [只含平方根以上素数的格映射](cm229-optimal-splitting-map.md)
