---
title: "判别式 −23 的三次域在高度 1024 内的全部非平凡零点位置"
date: "2026-11-25"
abstract: "将权一 η 乘积的 Mellin 积分旋转到固定复射线，再用带解析误差界的梯形求积和一次卷积验证符号。新增二维因子在正高度 1024 以内恰有 1845 个零点，全部是中心线上的单零点。"
---

# 判别式 −23 的三次域在高度 1024 内的全部非平凡零点位置

设 $K=\mathbb Q(\alpha)$，$\alpha^3-\alpha-1=0$，$L(s)=\zeta_K(s)/\zeta(s)$。**在 $0<\Re s<1$、$0<\Im s\le1024$ 内，$L$ 恰有 1845 个零点，均为中心线上的单零点。** 其中 1786 个位于此前证书的高度 64 之上。每个正纵坐标落在端点分母为 64 的一个开区间内。最后一个在 $(65514/64,65515/64)$。实轴没有该因子的非平凡零点，负高度由复共轭得到。[全部区间 CSV](zeros23-1024.csv) · [可独立生成区间并验证围道的程序](verify_zeros23_1024.py)

| 高度上限 | 64 | 128 | 256 | 512 | 1024 |
| --- | ---: | ---: | ---: | ---: | ---: |
| 二维因子的正高度零点数 | 59 | 146 | 348 | 809 | 1845 |

两个计数独立进行。围道按重数得到 1845，中心线网格得到 1845 次严格变号。所有较低的 59 个新区间也各自包含前稿给出的宽 $2^{-24}$ 的窄区间。另在 $t=0,8,32,64$、$\sigma=1/2,3/2$ 的八个点，两种独立求值所得复区间逐一相交。结合 $\zeta_K=\zeta L$ 和六次闭包的分解 $\zeta_M=\zeta_{\mathbb Q(\sqrt{-23})}L^2$，已核验的 Dirichlet GRH 便给出 $K$ 与 $M$ 的所有非平凡零点在 $|\Im s|\le1024$ 内的位置。单重性针对 $L$。在 $\zeta_M$ 中同一因子出现两次。[高度 64 的独立算法及数据](cubic23-zeros.md) · [上游 Dirichlet 定理的核验](https://math.mit.edu/~drew/grh/)

## 把指数衰减从数值积分中提出

令 $C=\sqrt{23}/(2\pi)$，$H(z)=\eta(iz/\sqrt{23})\eta(i\sqrt{23}\,z)$。在 $\Re z>0$ 内有 $H(z)=z^{-1}H(1/z)$，$H(\bar z)=\overline{H(z)}$，并有展开 $H(z)=\sum_{n\ge1}a(n)e^{-nz/C}$，$|a(n)|\le d(n)\le n$。写 $\Lambda(s)=C^s\Gamma(s)L(s)$。对折前的整函数积分为 $\Lambda(s)=\int_{\mathbb R}H(e^u)e^{su}du$，在实轴严格为正，满足 $\Lambda(s)=\Lambda(1-s)$。[Zagier 关于此 η 乘积的恒等式](https://www.its.caltech.edu/~matilde/Zagier123ModularForms.pdf)

在 $s=\sigma+it$、$t>0$ 时，$\Gamma(s)$ 的指数因子约为 $e^{-\pi t/2}$。取 $\alpha=\pi/2-1/16$，向上移动 $u$ 的积分线，得到

$$
 F(s)=e^{-i\alpha(s-1/2)}\Lambda(s)
     = e^{i\alpha/2}\int_{\mathbb R} H(e^{u+i\alpha})e^{su}\,du . \tag{1}
$$

被积函数在所扫过的水平带内解析，两端由 $H$ 的展开和变换式双指数衰减。式 (1) 在中心线上恰为 $e^{\alpha t}\Lambda(1/2+it)$，保留了实值和符号。计算时只剩 $e^{-t/16}$ 的指数衰减。Booker 对一般 Artin $L$ 函数给出了通过旋转 Mellin 积分并用 Fourier 方法批量求值的框架。下文针对导子 23 指定射线、所有截断参数和可复核余量。[Booker，Artin 函数的数值检验](https://arxiv.org/abs/math/0507502)

## 梯形求积的两个误差

令 $d=1/32$、$h=1/1024$、$J=9/h=9216$、$T=1024$。对 $|\beta-\alpha|\le d$ 有 $\cos\beta/C>d$。在 $u\ge0$，系数界给出

$$
 |H(e^{u+i\beta})|\le (1-e^{-d})^{-2}e^{-d e^u}.
$$

在负半轴用 $H(z)=z^{-1}H(1/z)$。于是 $-1/2\le\sigma\le3$ 时有一致积分界

$$
 \int_{\mathbb R}|H(e^{u+i\beta})|e^{\sigma u}du
 \le\frac{4}{(1-e^{-d})^2d^3}< 2^{29}=:B. \tag{2}
$$

这里用两边的最大指数 $\max(\sigma,1-\sigma)\le3$，并以 $\int_0^\infty y^2e^{-dy}dy=2/d^3$ 估计。对 $H(e^{u+i\alpha})e^{\sigma u}$ 的 Fourier 积分向上或向下移 $id$，其在频率 $\xi$ 处的模至多为 $Be^{-d|\xi|}$。由 Poisson 求和，步长 $h$ 时非零频移贡献的误差对 $|t|\le T$ 至多为

$$
 E_1=\frac{2B\exp(-d(2\pi/h-T))}{1-\exp(-2\pi d/h)}
     <4.06\cdot10^{-65}. \tag{3}
$$

在 $1/2\le\sigma\le2$ 时舍去 $|j|>J$ 的网格。因 $e^{2u}e^{-d e^u}$ 在 $u\ge9$ 递减，以从 9 起的积分控制离散尾和。用 $(1-e^{-d})^{-2}<4/d^2$ 得

$$
 E_2=8(de^9+1)e^{-de^9}/d^4<2.28\cdot10^{-101}. \tag{4}
$$

每次复数求值对实部、虚部都向外加上 $E_1+E_2$，中心线仅需实部。全部计算采用 288 比特球算术。[解析函数的梯形求积](https://doi.org/10.1137/130932132) · [FLINT 的模函数实现](https://flintlib.org/doc/acb_modular/index.html)

## 一次卷积计算整条网格的符号

记 $H_j=H(e^{jh+i\alpha})$，置 $P(z)=H_0/2+\sum_{j=1}^JH_j z^j$。由倒数变换合并正负采样点，(1) 的有限梯形和为

$$
 h\left\{e^{i\alpha/2}P(e^{hs})
       +\overline{e^{i\alpha/2}P(e^{h(1-\bar s)})}\right\}. \tag{5}
$$

在中心线上两项互为共轭。取 $t=m/64$，$0\le m\le65536$。设 $\lambda=h/64$、$v_j=e^{i\lambda j^2/2}$，则 $e^{i\lambda jm}=v_j v_m\overline{v_{m-j}}$。将 (5) 的系数先乘 $e^{jh/2}v_j$，再与 $\overline{v_{m-j}}$ 作一次多项式卷积，即得所有 $t$ 的值。每个加入误差后的实区间均有严格符号。此步给出 1845 对相邻、符号相反的端点，区间表由程序计算生成。[批量求值和符号检验源码](verify_zeros23_1024.py) · [有向舍入的球算术](https://python-flint.readthedocs.io/en/latest/)

## 把线上变号与全带计数相配

在矩形 $[-1,2]+i[0,T]$，底边的 $\Lambda$ 为正。右边 $\log L(2+it)$ 从欧拉积取连续分支，其虚部绝对值至多 $2\log\zeta(2)<\pi/2$。上边先从 $1/2+iT$ 向右追踪 $F$ 的辐角增量 $\Delta_F$。还原 $\Lambda$ 须加回 $3\alpha/2$，因为该段实部变化了 $3/2$。由函数方程及共轭对称性，矩形绕数为

$$
 \frac{T\log C+\Im\log\Gamma(2+iT)+\arg L(2+iT)-3\alpha/2-\Delta_F}{\pi}. \tag{6}
$$

为验证整条上边，将 $[1/2,2]$ 分为 384 段，间距 $1/256$。在上边固定高度，先把 $e^{ijhT}$ 乘入系数，随后只在**实数** $e^{h\sigma}$ 和 $e^{h(1-\sigma)}$ 处计算多项式，避免长复数递推导致区间膨胀。式 (1)、(2) 给出上边每个单位圆盘内 $|F|\le B$。对每段某一端点的较大实或虚坐标 $g$ 作 40 阶前向或反向 Newton 插值。段内与该端点的差至多为

$$
 \sum_{j=1}^{40}\frac{|\Delta^j g|}{j} + \frac{B\,40!}{256^{41}} .
$$

这来自 $|\binom vj|\le1/j$（$0\le v\le1$）、Cauchy 导数界与实插值余项。每段的起始坐标绝对值扣去该界，最小值仍大于 $1.676\cdot10^{-28}$。于是每段的像都在同一开半平面内，端点商的主辐角可直接相加。式 (6) 的区间包含在 $1845\pm1.4\cdot10^{-35}$ 内。其唯一整数等于线上变号数，故全部零点均已计入且重数为一。这个论证同时核验了上边无零点。

在仓库根目录运行 `python -B ERH/verify_zeros23_1024.py` 可从 η 函数起重新计算绕数与所有区间。附加 `--csv` 可输出完整端点。程序使用 `python-flint==0.9.0`。把这个有限高度换成实际素数估计，需给更高零点的总贡献一个独立的界。[素数平方次项与零点尾界](cubic23-gaussian-bias.md)
