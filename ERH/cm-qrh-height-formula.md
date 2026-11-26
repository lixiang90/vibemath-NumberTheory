---
title: "CM 平均高度的有限素理想公式与准黎曼半平面的等价判据"
date: "2026-11-27"
abstract: "对每个 CM 域构造相对与全域两种有限素理想逼近。同一平均 Faltings 高度的两个误差分别等于相对 Hecke 因子与完整 Dedekind ζ 函数的绝对收敛零点和。给出带显式判别式常数的全高度误差界，并刻画改进无零边界所需的几何估计。"
---

# CM 平均高度的有限素理想公式与准黎曼半平面的等价判据

令 $E$ 为任意 CM 域，$F=E^+$，$d=[F:\mathbb Q]$，$D_E,D_F$ 均表示绝对判别式的绝对值，$C=D_E/D_F$。设 $\eta=\eta_{E/F}$，$L=L_F(s,\eta)$，$\zeta_E=\zeta_F L$。记 $\bar h_E$ 为对 $2^d$ 个 CM 型的稳定 Faltings 高度平均，采用 Andreatta–Goren–Howard–Madapusi Pera 定理 A 的正规化。对任意固定 $1/2\le b<1$，本文给出等价关系

$$
 \zeta_E(s)\ne0\quad(\Re s>b,\ s\ne1)
 \quad\Longleftrightarrow\quad
 |\widehat h_{\rm all}(X)-\bar h_E|
     \ll_{E,\epsilon} X^{b-1+\epsilon}\quad(\forall\epsilon>0,\ X\ge2). \tag{1}
$$

下文的 $\widehat h_{\rm all}$ 只用 $E$ 的有限素理想、判别式和基域的 Euler–Kronecker 常数。若使用只含相对分裂信息的 $\widehat h_{\rm rel}$，同一判据检测 $L$ 的零点。Ihara 已对数域素理想与 Euler–Kronecker 常数证明显式公式。这里采用其权乘 $1-1/X$，与 CM 平均高度相配，并同时写出两个误差的归一化与逆推。[Ihara，§1.1–1.3](https://www.math.uci.edu/~mfried/conflist-rims/yihara10-14-07.pdf) · [平均 Colmez 定理 A](https://arxiv.org/abs/1508.00178) · [Yuan–Zhang 的独立证明](https://arxiv.org/abs/1507.06903)

## 两个可以计算的逼近

记 $\gamma$ 为欧拉常数，$\gamma_K$ 为 $\zeta_K$ 在 $1$ 处 Laurent 展开的常数项除以留数。令 $B_\eta=-L'/L$ 以及

$$
 h_0=\tfrac14\log C-\tfrac d2\gamma-d\log(2\pi),\qquad
 c_E=\log D_E-2d(\gamma+\log(2\pi)).
$$

相对特征在全部实位为奇。完备因子为 $C^{s/2}\Gamma_{\mathbb R}(s+1)^d$。将平均 Colmez 公式用函数方程移至 $1$，得到

$$
 \bar h_E=h_0-\tfrac12 B_\eta(1)
         =h_0+\tfrac12(\gamma_E-\gamma_F),\qquad
 B_\eta(0)=2\bar h_E+\tfrac12\log C+d\log(2\pi). \tag{2}
$$

最后一个等式直接来自 $0$ 处的平均 Colmez 公式。$1$ 处的 Euler–Kronecker 常数之差由 $\zeta_E=\zeta_F L$ 得到。[高度正规化及函数方程的推导](real-cm-prime-windows.md)

对 $X>1$，置 $q=X^{-1}$。只保留范数幂不超过 $X$ 的项，定义

$$
 S_\eta(X)=\sum_{N\mathfrak p^k\le X}
  \eta(\mathfrak p)^k\log N\mathfrak p
       \left(\frac1{N\mathfrak p^k}-q\right),\qquad
 S_E(X)=\sum_{N\mathfrak P^k\le X}
    \log N\mathfrak P\left(\frac1{N\mathfrak P^k}-q\right).
$$

第一和取 $F$ 的素理想，分歧处 $\eta(\mathfrak p)=0$。第二和取 $E$ 的全部素理想。两式的 $k$ 都从 $1$ 起，端点权均为零。定义两个实函数

$$
 A(q)=q\,\operatorname{arctanh}q+\tfrac12\log(1-q^2),\qquad
 D(q)=q+(1-q)\log(1-q).
$$

则所需的有限量为

$$
 \begin{aligned}
 \widehat h_{\rm rel}(X)
   &=\frac{h_0-S_\eta(X)/2-q(\tfrac14\log C+\tfrac d2\log(2\pi))
                   -dA(q)/2}{1+q},\\
 \widehat\gamma_E(X)
   &=\frac{\log X-1-S_E(X)+q((d-1)(1+\log X)-c_E)-dD(q)}{1+q},\\
 \widehat h_{\rm all}(X)&=h_0+\tfrac12(\widehat\gamma_E(X)-\gamma_F).
 \end{aligned} \tag{3}
$$

**带显式常数的结论。** 对 $X\ge2$ 及每个 CM 域均有

$$
 |\widehat h_{\rm all}(X)-\bar h_E|
 \le \frac{30}{7}(\log D_E+3)\frac{X^{-1/8}}{1+X^{-1}}. \tag{4}
$$

更一般地，在任一已证无零半平面 $\Re s>b$ 下，可将式 (4) 的 $30/7$ 和 $-1/8$ 分别换为 $(1+b)/(4b(1-b))$ 和 $b-1$。这允许使用 Liu 的一般数域界并得到 $X^{-0.1266673625540949\ldots}$ 量级，或以任何明确的有理上端点代入。这里的判别式常数对次数、分歧和 CM 型都没有隐含依赖。[Hu，定理 1.2](https://hexagonmath.org/2610.00216v2) · [全域边界与有理端点](hecke-qrh-contact.md)

## 零点等式与判据的证明

两种逼近的共同形式为

$$
 2(1+X^{-1})(\widehat h_{\mathcal L}(X)-\bar h_E)
   =\sum_{\mathcal L(\rho)=0}^{*}
             \frac{X^{\rho-1}}{\rho(\rho-1)},\qquad
 (\widehat h_{\mathcal L},\mathcal L)
       =(\widehat h_{\rm rel},L),\ (\widehat h_{\rm all},\zeta_E). \tag{5}
$$

星号仅取非平凡零点，计重数。对每个固定的域，两个级数均绝对收敛。两种误差相减恰好留下 $\zeta_F$ 的零点，包括它与相对因子有共同零点时属于基域的重数。

为证 (5)，对 $B=-\mathcal L'/\mathcal L$ 积分

$$
 \frac1{2\pi i}\int_{(2)} B(s)\frac{X^{s-1}}{s(s-1)}\,ds .
$$

逐项反演就是 (3) 所用的 $S_\eta$ 或 $S_E$。当 $\mathcal L=L$ 时，$1,0$ 的留数依次为 $B_\eta(1),-qB_\eta(0)$。相对特征非主，$L$ 在 $1$ 及 $0$ 非零，平凡零点全部在 $-1,-3,\ldots$，每处重数 $d$。一切非平凡零点的留数总和为式 (5) 右端的负数。平凡零点给 $-d A(q)$，这是因为 $A(q)=\sum_{k\ge0}q^{2k+2}/((2k+1)(2k+2))$。代入 (2) 并收集 $\bar h_E$ 给出相对式。

当 $\mathcal L=\zeta_E$ 时，$B_E(s)=1/(s-1)-\gamma_E+O(s-1)$，故 $1$ 处留数为 $\log X-1-\gamma_E$。$E$ 的无穷因子是 $\Gamma_{\mathbb C}(s)^d$，其中 $\Gamma_{\mathbb C}(s)=2(2\pi)^{-s}\Gamma(s)$。其函数方程和 $\psi(s)=-1/s-\gamma+O(s)$ 给出

$$
 B_E(s)=-(d-1)/s+\gamma_E+c_E+O_E(s)\qquad(s\to0).
$$

因此 $0$ 处留数为 $q((d-1)(1+\log X)-\gamma_E-c_E)$。在 $-1,-2,\ldots$ 的平凡零点每处重数 $d$，对应总和为 $-dD(q)$，因为 $D(q)=\sum_{k\ge1}q^{k+1}/(k(k+1))$。代入 (3) 就得到全域式。

上述移线可先移至 $\Re s=-2N-1/2$，它与两类平凡零点均有固定距离。用单位高度 $O_E(\log(2+|t|))$ 的零点数避开零点选取水平边，再令高度趋于无穷。对数导数的局部展开给这些水平边上 $O_{E,N}(\log^2(2+|t|))$ 的界。负线上由函数方程、伽马对数导数的反射公式以及右侧的绝对收敛，增长至多为 $O(\log D_E+d\log(2+N+|t|))$。固定 $X>1$ 再令 $N\to\infty$，负线积分趋于零。$|\rho(\rho-1)|^{-1}$ 对大高度是 $O((1+|\Im\rho|)^{-2})$，这也保证 (5) 的绝对收敛。[Hecke 显式公式与零点计数](https://aareyanmanzoor.github.io/assets/articles/lagarias-odlyzko.pdf) · [Ihara 的数域显式公式](https://www.math.uci.edu/~mfried/conflist-rims/yihara10-14-07.pdf)

若右半平面无零，函数方程给所有非平凡零点 $1-b\le\beta=\Re\rho\le b$。写 $\gamma_\rho=\Im\rho$，直接平方可见

$$
 |\rho(\rho-1)|\ge\beta(1-\beta)+\gamma_\rho^2,\qquad
 \frac1{|\rho(\rho-1)|}\le
       \frac{1+b}{b(1-b)}\,\Re\frac1{2-\rho}.
$$

将极点消去后的 $\xi_E(s)=s(s-1)D_E^{s/2}\Gamma_{\mathbb C}(s)^d\zeta_E(s)$ 取 Hadamard 乘积。关于 $1/2$ 对称配对，线性指数项由函数方程消失，故 $\sum_\rho\Re(2-\rho)^{-1}=\xi_E'(2)/\xi_E(2)$。显式求导后，$\zeta_E'(2)/\zeta_E(2)\le0$，并且 $\psi(2)-\log(2\pi)<0$，得到

$$
 \sum_\rho\Re\frac1{2-\rho}\le\tfrac12\log D_E+\tfrac32 .
$$

代回 (5) 给出 (4) 及一般 $b$ 的结论。对于 $L$ 用无极点的完备函数取相同乘积，右端可用 $\tfrac12\log C+2d$ 代替，因为 $|L'/L(2)|\le-d\zeta'/\zeta(2)<2d$，且 $\Gamma_{\mathbb R}'(3)/\Gamma_{\mathbb R}(3)<0$。这给每个固定相对二次特征同样的幂次。

反过来，对 (5) 的右端在 $X\ge1$ 作 Mellin 变换，在 $\Re w>0$ 得到

$$
 \sum_{\rho}^{*}\frac1{\rho(\rho-1)(w+1-\rho)} .
$$

它在极点集以外局部正规收敛，在 $w=\rho-1$ 处留数为 $m(\rho)/(\rho(\rho-1))\ne0$。若 (1) 的右端对每个 $\epsilon>0$ 成立，则变换由积分在 $\Re w>b-1$ 全纯，排除了每个 $\Re\rho>b$ 的零点。这同时证明全域与相对两个判据。进一步，若想避开预先求 $\gamma_F$ 或 $\bar h_E$，可以测试纯有限量 $\widehat\gamma_E(2X)-\widehat\gamma_E(X)$。对固定 $b<1$，同样的 $O_{E,\epsilon}(X^{b-1+\epsilon})$ 对所有 $X\ge2$ 成立，也与 (1) 等价。事实上已有 (4) 给 $\widehat\gamma_E(X)\to\gamma_E$，取 $0<\epsilon<1-b$ 后在 $X,2X,4X,\ldots$ 将差求和为收敛几何级数。

因此，若目标是把一般 CM 域边界推进至 $0.87333$，一个几何证明须独立给每个固定 $E$ 的全域误差 $O_{E,\epsilon}(X^{-0.12667+\epsilon})$，或给出上述成倍截断差的同阶界。全型平均本身是 (2) 的精确常数。单位平衡能控制无穷位的形状。从这样的控制推出所需衰减率仍缺一个随 $X$ 定量的估计。

对判别式 $229$ 的例子，可逐项扣除低零点，再把剩余谱质量转换成 (5) 的误差常数。由此得到对所有 $X>1$ 的两项界，且在 $1024\le X\le\exp(10^6)$ 两个域的零点和都满足 $0.92\,X^{-0.12667}$ 的上界。[两个具体域的全尺度估计](cm229-qrh-boundary.md)

## 两个域上的有限检验

下表给出 $\widehat h-\bar h_E$ 的严格区间下端，上端为该数加 $10^{-9}$。六次域取 $F=\mathbb Q(\alpha)$、$\alpha^3-4\alpha+1=0$、$E=F(i)$，有 $D_F=229$、$D_E=3356224$。两域的平均高度分别为 $-1.65710651619148\ldots$ 和 $-4.12205285339793\ldots$，采用 (2) 的正规化。

| $E$ | $X$ | 相对误差 | 全域误差 |
| --- | ---: | ---: | ---: |
| $\mathbb Q(i)$ | 256 | 0.000440208 | 0.000681058 |
| $\mathbb Q(i)$ | 4096 | -0.000296446 | -0.000297794 |
| $\mathbb Q(i)$ | 65536 | 0.000015160 | -0.000009353 |
| 六次域 | 256 | 0.016262349 | 0.017202168 |
| 六次域 | 4096 | -0.020238363 | -0.017016608 |
| 六次域 | 65536 | 0.003245038 | 0.003836463 |

程序从剩余次数与二次特征分别构造素数幂系数，包括 $2,229$ 的分歧项。对六次域的偶、奇 Artin 因子分别在 $1$ 处求对数导数，以复现基域常数与平均高度。还独立核对了两种误差之差所对应的基域零点和。使用最前面 48 个有理 ζ 正高度零点，六次域另加既有证书中 82 个偶因子零点。其余零点由总量 $\sum_\rho\Re(1/(\rho(1-\rho)))=2\xi_F'(1)/\xi_F(1)$ 扣掉已用的正权来估计。在 $1/8\le\Re\rho\le7/8$，复权模长至多为其实部的 $4/\sqrt7$ 倍，由 $|\gamma_\rho|/(\beta(1-\beta)+\gamma_\rho^2)\le1/(2\sqrt{\beta(1-\beta)})$ 得到。整个检验只用一般的 $7/8$ 半平面。[区间复算程序](verify_cm_qrh_heights.py) · [六次域的局部因子与隔离证书](cm229-heights-and-zeros.md)
