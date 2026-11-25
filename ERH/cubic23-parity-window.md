---
title: "判别式 −23 的平方自由理想奇偶与趋零斜率的剩余区间"
date: "2026-11-25"
abstract: "在所有正的小斜率上证明全整数系数均方和不分奇偶的理想类估计。完全分裂素数的单侧下界可进一步归约到具有偶数个素理想因子的平方自由理想，其范数指数落在宽 O(√λ) 的区间内，余因子的范数至多为 Y 的 O(√λ) 次幂。"
---

# 判别式 −23 的平方自由理想奇偶与趋零斜率的剩余区间

设 $k=\mathbb Q(\sqrt{-23})$，$\psi$ 是其理想类群的非平凡三阶特征，$L=L_k(s,\psi)=\sum a(m)m^{-s}$。共轭作用将 $\psi$ 变成 $\psi^2=\bar\psi$，且保持理想范数，故 $L_k(s,\psi^2)=L(s)$。记 $Y=e^u$，$v=\lambda u$，$0<\lambda\le1/16$，$\sigma_\lambda=1/2+\lambda/8$，并沿用归一化高斯 $k_v(t)=(2\pi v)^{-1/2}e^{-t^2/(2v)}$。以下结果在每个固定斜率上对任意大的 $Y$ 成立。[域与权一形式的识别](cubic23-zeros.md) · [二次扭转均方中的素数对角线](cubic23-class-character-mean-square.md)

## 两种均可达到平方根精度的平均

对全部正整数，平方系数满足

$$
 D(s)=\sum_{m\ge1}\frac{a(m)^2}{m^s}
 =\frac{\zeta_k(s)L(s)}{\zeta(2s)(1+23^{-s})},\qquad \Re s>1. \tag{1}
$$

对未分歧 $p$，以 $\alpha,\beta$ 记标准表示 $\rho$ 的局部特征值，$\alpha\beta=\chi_{-23}(p)$。平方系数的生成函数等于 $\rho\otimes\rho$ 的局部因子乘 $1-p^{-2s}$。这可将 $a(p^j)=\sum_{h=0}^j\alpha^h\beta^{j-h}$ 平方求和验证，包括 $\alpha=\beta$ 的情形。由 $C_3$ 特征表或 $S_3$ 的标准表示，$\rho\otimes\rho=1\oplus\mathrm{sgn}\oplus\rho$。在 23 处 $a(23^j)=1$，上式额外的 $1+23^{-s}$ 恰使局部因子为 $(1-23^{-s})^{-1}$。特别在未分歧素数上有精确关系 $a(p)^2=1+\chi_{-23}(p)+a(p)$。[局部因子及诱导表示](cubic23-zeros.md)

令 $R=\operatorname{Res}_{s=1}D(s)>0$。对每个 $\varepsilon>0$、$u$ 充分大，

$$
 \sum_{m\ge1} a(m)^2 k_v(u-\log m)
 =R e^{u+v/2}+O_{\lambda,\varepsilon}(e^{(\sigma_\lambda+\varepsilon)u}).
 \tag{2}
$$

取 $c=1/2+\delta<1$，$\delta>0$ 足够小，使 $c+\lambda c^2/2<\sigma_\lambda+\varepsilon$。由 (1)，将 $D(s)e^{su+vs^2/2}$ 的反演线移到 $c$，只越过极点 1。$1/\zeta(2s)$ 在该线上绝对收敛，$1+23^{-s}$ 一致远离零，分子在固定竖带内有多项式增长。$v\ge1$ 后高斯给出所需积分界。这证明了所有整数上的系数均方估计。

在理想一侧，记 $\mu_k$ 为理想 Möbius 函数。由 Euler 因子逐项相乘，

$$
 \sum_{\mathfrak d}\frac{\mu_k^2(\mathfrak d)\psi(\mathfrak d)}
                              {(N\mathfrak d)^s}=\frac{L(s)}{L(2s)},\qquad
 \sum_{\mathfrak d}\frac{\mu_k(\mathfrak d)\psi(\mathfrak d)}
                              {(N\mathfrak d)^s}=\frac1{L(s)}. \tag{3}
$$

这两式起初在 $\Re s>1$ 成立。给每个有序对 $(\mathfrak d,\mathfrak m)$ 赋权 $\log N\mathfrak m\, k_v(u-\log N(\mathfrak d\mathfrak m))$。记 $U_v(u)$ 为在平方自由 $\mathfrak d$ 和所有整理想 $\mathfrak m$ 上，以 $\Re\psi(\mathfrak d\mathfrak m)$ 乘此权所得的和。其变换是 $-L'(s)L(s)/L(2s)$。$2c>1$，所以 $1/L(2s)$ 在上一段的竖线绝对收敛。移线时无极点，于是

$$
 |U_v(u)|\ll_{\lambda,\varepsilon} e^{(\sigma_\lambda+\varepsilon)u}. \tag{4}
$$

式 (2)、(4) 同时适用于任何给定的趋零正斜率序列，常数可随斜率变化。两式的证明只需 Hecke 函数的经典解析延拓与函数方程，以及绝对收敛半平面内的 Euler 积。[此函数的 Mellin 积分与函数方程](cubic23-zeros.md)

## 偶数个素理想因子留下的范围

设
$$
 h_\lambda=\sqrt{\lambda+3\lambda^2/4},\qquad
 \eta_-=1+\lambda-h_\lambda,\qquad \eta_+=1+\lambda+h_\lambda.
$$

给定一个理想类 $C$，令 $W_C(u)$ 为下述有序对的**非负**加权总量，

$$
 \begin{gathered}
 \mathfrak d\ \text{平方自由},\qquad
 \omega_k(\mathfrak d)\equiv0\pmod2,\qquad
 Y^{\eta_-}\le N\mathfrak d\le Y^{\eta_+},\\
 1\le N\mathfrak m\le Y^{2h_\lambda},\qquad
 [\mathfrak d\mathfrak m]=C,\qquad
 \text{权}=(\log N\mathfrak m) k_v(u-\log N(\mathfrak d\mathfrak m)).
 \end{gathered} \tag{5}
$$

$\omega_k$ 数的是不同素理想。两个因子允许有公共素理想，所有分歧处均计入。记 $C_0$ 为主类，$C_1,C_2$ 为互逆的两个非主类。共轭给出 $W_{C_1}=W_{C_2}$。

**定理。** 令 $P_v$ 为六次闭包中完全分裂的未分歧有理素数的 $\log p$ 高斯和。在已核验的 Dirichlet GRH 下，对每个固定 $\lambda,\varepsilon>0$（$\lambda\le1/16$）有

$$
 6P_v(u)-e^{u+v/2}
 =4\big(W_{C_0}(u)-W_{C_1}(u)\big)
   +O_{\lambda,\varepsilon}(e^{(\sigma_\lambda+\varepsilon)u}). \tag{6}
$$

去掉 (5) 的偶数限制、保留平方自由和两个范数限制时，相应两类之差本身已经是 $O_{\lambda,\varepsilon}(e^{(\sigma_\lambda+\varepsilon)u})$。因此，所需的单侧下界只须控制选择偶数个素理想因子后留下的类差。例如取 $\lambda=2^{-30}$，(5) 中的余因子已经满足 $N\mathfrak m\le Y^{0.00006104}$。

**证明。** 定义实值函数

$$
 B_v(w)=\sum_{\mathfrak m}\psi(\mathfrak m)\log N\mathfrak m\,
                              k_v(w-\log N\mathfrak m).
$$

共轭理想使虚部相消。它的双边 Laplace 变换是 $-L'(s)e^{vs^2/2}$，在任何固定竖线 $\Re s=c$ 均可反演。由整函数 $L'$ 在固定竖带的多项式增长，对 $v\ge1$ 有

$$
 |B_v(w)|\le C_c e^{cw+vc^2/2},\qquad w\in\mathbb R. \tag{7}
$$

记 $S_L^v$ 为 $-L'/L$ 的系数按同一高斯加权的和。在 $\Re s>1$ 用 (3) 乘以 $-L'(s)$ 再反演，得到

$$
 S_L^v(u)=\sum_{\mathfrak d}\mu_k(\mathfrak d)\Re\psi(\mathfrak d)
                        B_v(u-\log N\mathfrak d).
$$

若偶数、奇数个素理想因子的平方自由 $\mathfrak d$ 各自给出无 Möbius 符号的贡献 $V_e,V_o$，便有 $S_L^v=V_e-V_o=2V_e-U_v$。同时用六次闭包的分解、$\zeta_k$ 的 Dirichlet GRH 及高次素数幂界，可得 $6P_v-e^{u+v/2}=2S_L^v+O(e^{\sigma_\lambda u})$。[素数与素数幂之间的误差](cubic23-gaussian-threshold.md)

以下把 $V_e$ 截到 (5)。固定域 $k$ 的范数不超过 $x$ 的整理想数为 $O_k(x)$，在这里也可由三个正定二元二次型的格点计数直接得到。[三个二次型](cubic23-zeros.md) 对 $N\mathfrak d<Y^{\eta_-}$，在 (7) 取 $c_-=(\eta_--1)/\lambda<1$。即使对 $\mathfrak d$ 的系数逐项取绝对值，删掉的贡献也至多为常数乘

$$
 \exp\big((c_-+\lambda c_-^2/2+\eta_-(1-c_-))u\big)
   =\exp\big((\eta_--(\eta_--1)^2/(2\lambda))u\big)
   = e^{\sigma_\lambda u}.
$$

对于 $N\mathfrak d>Y^{\eta_+}$ 取 $c_+=(\eta_+-1)/\lambda>1$。此时 $\sum_{N\mathfrak d>D}(N\mathfrak d)^{-c_+}\ll_\lambda D^{1-c_+}$，同一计算在右端也恰给出 $e^{\sigma_\lambda u}$。这两个端点是 (7) 导出的指数与所需平方根指数的交点。

余下若 $N\mathfrak m>Y^{2h_\lambda}$，则 $N(\mathfrak d\mathfrak m)>Y^{\eta_+}$。范数乘积不超过 $x$ 的有序对按 $\log N\mathfrak m$ 加权共为 $O_k(x\log^2(2x))$，由上面的理想计数求和即得。再对高斯用 $k_v(t)\le(2\pi v)^{-1/2}e^{c_+t+vc_+^2/2}$ 并分部求和，这一截尾至多 $O_\lambda((1+u)^2e^{\sigma_\lambda u})$。全部截断对奇、偶或两者的并集都有效，因为对外层理想只用模长界。截后 $V_e$ 正是 $\sum_C\Re\psi(C)W_C=W_{C_0}-W_{C_1}$。结合 (4) 并吸收多项式因子即证 (6)。同样截断 (4) 即得不分奇偶的断言。

在原整数窗口判据中，容差指数可以任意小。因此 (6) 可以直接用 $W_{C_0}-W_{C_1}$ 替换素数误差，只须在应用余项时先取更小的容差。等价地，在伴篇的斜率平均密度公式中作同样替换，再令容差趋零。需要证明的算术超额集中在 (5) 的偶数层。[整数窗口判据](cubic23-gaussian-threshold.md) · [需要严格超过的一半密度](cubic23-class-character-mean-square.md)

复算程序以 $\eta(z)\eta(23z)$ 的五边形展开取得 $a(m)$，独立用素数在三次域中的分解构造带符号和不带符号的平方自由理想系数。它逐项核验 (1)、(3) 到范数 32768，并将 $\log$ 分别换成 $2,5,23,59$ 的加性赋值，精确核对 $2V_e-U$ 的卷积系数，包括分歧素数 23。[整数复算程序](verify_class23_moments.py)
