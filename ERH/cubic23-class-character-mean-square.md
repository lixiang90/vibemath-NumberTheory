---
title: "三阶理想类特征的二次扭转均方与二分之一密度判据"
date: "2026-11-25"
abstract: "在判别式 −23 的三次域中，两个共轭类特征之差的二次扭转均方精确检出完全分裂素数的短缺。对几乎每个正斜率，即使趋于零，平方根级单侧下界也在至少一半的整数对数尺度上成立。严格超过一半的部分给出最小平滑斜率的反演公式。"
---

# 三阶理想类特征的二次扭转均方与二分之一密度判据

令 $K=\mathbb Q(\alpha)$，$\alpha^3-\alpha-1=0$，$k=\mathbb Q(\sqrt{-23})$，$M$ 为六次分裂域。$M/k$ 是 Hilbert 类域，$\mathrm{Cl}(k)\simeq C_3$。取非平凡类特征 $\psi$，记 $L=L_k(s,\psi)=\zeta_K/\zeta$。下文使用已核验的 Dirichlet GRH。对 $L$ 的零点只用函数方程、显式公式及已有的低零点证书。[域和类特征](cubic23-zeros.md) · [上游定理核验](https://math.mit.edu/~drew/grh/) · [斜率与零点的公式](cubic23-gaussian-threshold.md)

本篇证明，原均方中的非对角项可以在平方根精度以内全部扣除。三阶类特征决定了余下的对角项。对几乎每个 $t\in[1,2]$，沿斜率 $t2^{-j}$ 可以先取得至少一半尺度的单侧下界。把这个密度严格推过一半，就足以排除相应的偏线零点。

## 类特征给出的对角项

设 $u=\log Y$，$v=\lambda u$，$0<\lambda\le1/16$，$k_v(t)=e^{-t^2/(2v)}/\sqrt{2\pi v}$，写 $w_p=(\log p)k_v(u-\log p)$。令 $Q_v(u)$ 为在 $k$ 中分裂的有理素数的 $w_p$ 之和，$P_v(u)$ 为在 $M$ 完全分裂者的和，均除去 23。对每个在 $k$ 分裂的 $p$，选其上一素理想 $\mathfrak p$，令

$$
 d_p=\psi(\mathfrak p)-\overline{\psi(\mathfrak p)}.
$$

主类时 $|d_p|^2=0$，其余两类时 $|d_p|^2=3$。交换 $\mathfrak p,\overline{\mathfrak p}$ 只改变符号。置 $X=Y^2$，$U=Y^{10}$，取辅助素数族 $\mathcal Q=\{q:\sqrt U<q\le U,\ q\equiv3\bmod8\}$。定义

$$
 M_\lambda(u)=\frac1{|\mathcal Q|}\sum_{q\in\mathcal Q}
 \left|\sum_{\substack{p\le X\\\chi_{-23}(p)=1}}
        d_p\sqrt{w_p}\left(\frac pq\right)\right|^2 . \tag{1}
$$

**定理 1。** 对 $1\le v\le u/16$、$u$ 充分大，一致有

$$
 M_\lambda(u)=3\big(Q_v(u)-P_v(u)\big)+O(e^{-u}(1+u)^4). \tag{2}
$$

常数与素理想的选向无关。因此 $M_\lambda-2Q_v$ 正好测量 $Q_v-3P_v$，误差趋零。

**证明。** 展开平方。$p=r$ 时 $(p/q)^2=1$。若 $p\ne r$，将 $q\equiv3\bmod8$ 的指标展开为模 8 的四个特征。由二次互反律，$q\mapsto(pr/q)$ 在奇数上由导子整除 $8pr$ 的二次特征给出。$pr$ 至少有一个奇素因子，其上有非平凡分歧，乘模 8 特征后仍非主。Dirichlet GRH 的素数显式公式和分部求和，连同 $|\mathcal Q|\asymp U/\log U$，给出一致界

$$
 \left|\frac1{|\mathcal Q|}\sum_{q\in\mathcal Q}(pr/q)\right|
 \ll U^{-1/2}\log^3(8UX^2).
$$

因 $\sum_{p\le X}|d_p|\sqrt{w_p}\ll X\sqrt{\log X}$，非对角项共为 $O(Y^{-1}(1+u)^4)$。截尾用 $k_v(t)\le (2\pi v)^{-1/2}e^{8t+32v}$ 及 $\sum_{p>X}\log(p)p^{-8}\ll X^{-7}$，得到 $O(Y^{-4})$。剩余的平方对角线给出 (2)。[Dirichlet 特征素数和的 GRH 界](https://arxiv.org/abs/2003.01925)

原代码的 `QH.MeanSq` 对任意有界系数取绝对上界，`meanSq_abs` 在展开后估计成对出现的 Jacobi 符号，`meanSqK` 则用于指定的大导子区间和斜坡权。这里沿用其展开顺序，将测试系数改为 $d_p\sqrt{w_p}$，保留完全相同的素数这一对角线。公式 (2) 可以在多项式大小的辅助导子上得到，使用的是上游最终的 Dirichlet GRH。(2) 给出针对三阶类特征与该测试权的主项和加性误差。[原均方的定义](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/MeanSquareBounds/PrimeSetTwoBounds.lean) · [原平方展开与绝对估计](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/MeanSquareBounds/MeanSqAbs.lean) · [原窗口与辅助导子定理](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/MeanSquareBounds/meanSqK.lean)

仅用 $|d_p|^2\le3$，(2) 只能给出至多 $3Q_v$ 的主尺度。所需上界中心是 $2Q_v$。例如在只施加局部三次单位根条件的模型中，若每个分裂素数都赋非平凡值，对角就等于 $3Q_v$。要取得 $2Q_v$ 加平方根误差，须用到真实理想类在素数间的分布。下一节给出可以容许多少取样尺度失效。

## 对一列趋零斜率已有多少个好尺度

置 $\sigma_\lambda=1/2+\lambda/8$，$E_\lambda(n)=e^{(\sigma_\lambda+\varepsilon)n}$，其中 $\varepsilon>0$。对 $1\le n\le N$ 的整数，定义下列两种好尺度的比例，$v=\lambda n$，

$$
 \begin{aligned}
 G_{\lambda,\varepsilon}(N)&=\frac1N\#\left\{n:
      P_v(n)\ge e^{n+v/2}/6-E_\lambda(n)\right\},\\
 H_{\lambda,\varepsilon}(N)&=\frac1N\#\left\{n:
      M_\lambda(n)\le2Q_v(n)+E_\lambda(n)\right\}.
 \end{aligned} \tag{3}
$$

若 (1) 中辅助素数族为空，约定该点的 $M_\lambda(n)=0$。每个固定正斜率下，有限开头对密度无影响。对 $L$ 的非平凡零点 $\rho=\beta+i\gamma$，设

$$
 r(\varepsilon)=\sup\!\left(\{0\}\cup
 \left\{\frac{2\beta-1-2\varepsilon}{\gamma^2-\beta^2+1/4}:
                      \beta>1/2+\varepsilon\right\}\right).
 \tag{4}
$$

分母为正，因为全部偏线零点的绝对高度都大于已认证的高度。$r(\varepsilon)\uparrow\lambda_{23}$ 当 $\varepsilon\downarrow0$。[原阈值的定义与显式公式](cubic23-gaussian-threshold.md)

**定理 2。** 固定 $\varepsilon>0$。对几乎每个 $0<\lambda<1/16$，(3) 的两个自然密度都存在，且均为

$$
 \begin{cases}
 1/2,& 0<\lambda<r(\varepsilon),\\
 1,& \lambda>r(\varepsilon).
 \end{cases} \tag{5}
$$

对几乎每个 $t\in[1,2]$，沿同一列 $\lambda_j=t2^{-j}$（$j\ge6$），**每个斜率、每个 $\varepsilon>0$ 的好尺度下密度都至少为 $1/2$**。此断言对是否存在偏线零点一律成立。

**证明。** 写 $\Phi_\lambda(\rho)=\rho+\lambda\rho^2/2$。由前稿的移线式，$6P_v(n)-e^{n+v/2}$ 等于 $-2\sum_\rho e^{n\Phi_\lambda(\rho)}+O_\lambda(e^{\sigma_\lambda n})$。另外 $2Q_v(n)=e^{n+v/2}+O_\lambda(e^{\sigma_\lambda n})$，这由 $\zeta_k$ 的 Dirichlet GRH、扣除高次素数幂及 23 得到。与 (2) 合用可知，均方超出 $2Q_v$ 的量等于素数短缺的三倍加 $O_\lambda(e^{\sigma_\lambda n})$。[显式公式及高次幂的统一控制](cubic23-gaussian-threshold.md)

某个 $\Re\Phi_\lambda(\rho)>\sigma_\lambda+\varepsilon$ 恰好对应 $\lambda<r(\varepsilon)$。高处零点的 $\Re\Phi_\lambda(\rho)$ 一致趋于负无穷，故只要超过阈值，最大值 $B$ 由有限组零点达到。实轴上完备函数严格为正。将共轭对视为一组，不同组的 $\beta+\lambda(\beta^2-\gamma^2)/2$ 是不同的仿射函数，至多交于一点。去掉可数个交点，再去掉使 $(1+\lambda\beta)\gamma/(2\pi)$ 为有理数的可数个斜率。在其余斜率上，最大项来自唯一一对共轭零点。其实部与其余项有严格间隙。按重数 $m$ 计，素数短缺满足

$$
 e^{-Bn}\big(e^{n+v/2}/6-P_v(n)\big)
   = (2m/3)\cos\!\big(n(1+\lambda\beta)\gamma\big)+o(1).
$$

余项由 $\sum_\rho e^{u_0\Re\Phi_\lambda(\rho)}<\infty$（任取 $u_0>0$）及指数间隙控制。对均方短缺 $M_\lambda(n)-2Q_v(n)$，余弦系数为 $2m$。当 $B>\sigma_\lambda+\varepsilon$，归一化后的容许误差趋零。无理旋转在圆周上等分布，余弦正负各占一半。先排除 $|\cos|\le\eta$ 的弧再令 $\eta\downarrow0$，得到 (5) 第一行。

若 $\lambda>r(\varepsilon)$，每项的实部严格小于 $\sigma_\lambda+\varepsilon$。任何可能达到等号的零点必须有 $\beta>1/2+\varepsilon$，其式 (4) 的商才可能等于正的 $\lambda$。由高处衰减，最大实部亦严格小于阈值，故两个误差最终都在允许范围内。对每个有理 $\varepsilon>0$ 同时剔除零测集，再对各 $j$ 拉回到变量 $t$。任意实 $\varepsilon>0$ 的好尺度包含每一个更小有理容差的好尺度，这证明了趋零列的断言。

## 用超过一半的密度恢复阈值

对任意 $0<a<b\le1/16$，(3) 中任取其一记作 $J_{\lambda,\varepsilon}(N)$。由有界收敛及 (5)，

$$
 \lim_{\varepsilon\downarrow0}\lim_{N\to\infty}
 \frac1{b-a}\int_a^b J_{\lambda,\varepsilon}(N)\,d\lambda
 =1-\frac{|(a,b)\cap(0,\lambda_{23})|}{2(b-a)}. \tag{6}
$$

因此，只须对一列 $b_j=2a_j\downarrow0$ 证明 (6) 左侧严格大于 $1/2$，便有 $\lambda_{23}=0$。所需超额允许随 $j$ 任意衰减。用类特征表达时，(6) 左侧是有限均方 (1) 落在 $2Q_v$ 加容许误差以内的频率，再对斜率取平均和极限。公式 (2) 已将辅助素数平均造成的误差压低到 $O(Y^{-1}\log^4Y)$，严格超额仍取决于平方对角线的算术分布。[把对角线进一步化为平方自由理想的奇偶问题](cubic23-parity-window.md)

对每个正斜率还可给出全尺度的带指数误差。取 $T=8192$，若 $\Re\rho\le b\le1$ 是已知零点上界，则对 $1\le\lambda u$、$0<\lambda\le1/16$，

$$
 |M_\lambda(u)-2Q_v(u)|
 \ll e^{\sigma_\lambda u}+(4+T^2)
          \exp\big((b+\lambda(b^2-T^2)/2)u\big). \tag{7}
$$

由 (2) 和素数显式公式，左侧由 $L$ 的零点和及平方根级余项控制。证书内的零点均在中心线上。对更高零点用 $\sum_\rho(4+\gamma^2)^{-1}<5$ 以及 $(4+t^2)e^{-vt^2/2}$ 在 $t\ge T$ 递减，即得第二项。在 $\lambda=2^{-26}$ 上取已形式化的 $b=31/32$，第二项的指数已经低于 $\sigma_\lambda$。再向零取斜率时，(7) 仍给出在每个充分大尺度成立的估计，(5) 则给出达到平方根精度的尺度比例。[高度 8192 与斜率的数值界](cubic23-height-8192.md) · [零点总量界](cubic23-gaussian-bias.md)

作为有限核对，取 $p\le101$、$\chi_{-23}(p)=1$，权全为 1，对非主类统一取 $d_p=i\sqrt3$。11 个分裂素数中有 9 个非主类，对角线精确等于 27。对 $q\equiv3\bmod8$ 在 $(10^3,10^4]$、$(10^4,10^5]$、$(10^5,10^6]$ 上分别取平均，均方为 $2491/89,\ 28971/1049,\ 12919/479$。程序独立累加平方与非对角和，并用有理数核对两者之差。自然密度 (5) 由上述证明得到。[精确复算程序](verify_class23_moments.py) · [目前已保证全尺度下界的斜率区间](cubic23-height-8192.md)

后续研究在级 23 的三维 CM 模形式空间中构造了保持整数格的 Hecke 算子链。用两种权可求出三循环素数的全部指定公因子；固定权 7 后，首个非零 Hasse 导数与这条均方对角线之间的额外高度已降至平方根量级。[theta 特征线与整数格](cubic23-regular-theta-packet.md) · [局部公因子](cubic23-stepanov-content.md) · [整体导数高度与剩余上界](cubic23-cyclotomic-jets.md)
