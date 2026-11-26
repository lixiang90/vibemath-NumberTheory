---
title: "三阶双线性和在完整窗口上的有限高度零点公式"
date: "2026-11-27"
abstract: "对原式（9）的带符号双线性和，证明临界线上及其左侧所有零点的归一化整窗绝对值贡献为 O(log log Y)。给出依赖窗口的有限高度公式及现有无零半平面下的对数节省。在本域 ERH 下得到绝对值上界和带符号总和的衰减，并用最低零点给出累计正部的下界。"
---

# 三阶双线性和在完整窗口上的有限高度零点公式

沿用 $F=\mathbb Q(\alpha)$、$\alpha^3-4\alpha+1=0$、$E=F(i)$ 及 $M/B=F^{\rm gal}(i)/\mathbb Q(\sqrt{229},i)$。令 $L=L_B(s,\psi)=L_+L_-$，其中 $\psi$ 为非平凡三阶特征。这里研究的是原式（9）的 $R_Y(n)$，其短区间宽度为 $v=\log r=(\log2)/\lceil\log^2(2Y)\rceil$。记 $X=Y^{3/2}$，$A=\log Y$，$W=\lceil Y^3\rceil$，$k_n=\lfloor W/n^{3/2}\rfloor$，以及

$$
 P_R(Y)=\frac1W\sum_{Y\le n<X}k_n[R_Y(n)]_+,\qquad
 N_R(Y)=\frac1W\sum_{Y\le n<X}k_n|R_Y(n)|.
$$

伴篇已对一整段同时满足 $h,\ell>z$ 的项证明多对数上界，并给出了两因子都在四分之一次幂以上的替代卷积。本篇估计全部因子的总贡献。[原始双线性和](cm229-cubic-hecke-carries.md) · [可删除的因子区间](cm229-type-ii-ranges.md)

**定理 1。** 对 $K=\sqrt X$、$H=K\sqrt{8\log(2X)}$ 和 $0\le\lambda\le1/2$，用重数计 $L$ 的非平凡零点 $\rho=\beta+i\gamma$，置

$$
 Z_{Y,\lambda}(t)=
 -\frac12\sum_{\substack{\beta>1/2+\lambda\\64<|\gamma|\le H}}
 2^\rho\,\frac{1-e^{-v\rho}}{\rho}
 \exp\!\left(\frac{\rho^2}{2K^2}\right)e^{(\rho-1/2)t}.
 \tag{1}
$$

共轭零点使此和实值。对 $f(x)=[x]_+$ 或 $f(x)=|x|$ 都有

$$
 \left|\frac1W\sum_n k_n f(R_Y(n))-
 \int_A^{3A/2}f(Z_{Y,\lambda}(t))\,dt\right|
 \ll 1+
 \frac{\log(2+A)}{A}\,X^\lambda
 \min(1+A,1+1/\lambda). \tag{2}
$$

取 $\lambda=0$ 时，最小值取 $1+A$。界对所有充分大的实数 $Y$ 一致。特别地，对某个窗口，只要高度 $H$ 以内没有 $\beta>1/2+1/\log X$ 的零点，该窗口便有 $N_R(Y)\ll\log(2+\log Y)$。若允许 $Y^{o(1)}$ 的预算，则还可取 $\lambda=1/\sqrt{\log X}$，把比临界线向右伸出这么宽的一条带也计入误差。已认证的高度 $64$ 以内没有临界线外零点。[低零点的全计数与隔离](cm229-heights-and-zeros.md)

对一般 $Y$，采用 Liu 的无零半平面，取其上界 $b=0.8733326374459051924$，将 $\lambda=b-1/2$ 代入 (2)，有限和为空。得到全双线性和的**绝对值**上界

$$
 N_R(Y)\ll 1+
 Y^{\frac32(b-1/2)}\frac{\log(2+\log Y)}{\log Y}. \tag{3}
$$

若 $E$ 满足 ERH，则 $N_R(Y)\ll\log(2+\log Y)$。若 $L$ 存在临界线右侧零点，$P_R(Y)$ 沿所有充分远窗口的精确幂指数仍为 $\tfrac32(\beta_*-1/2)>0$，其中 $\beta_*=\sup_\rho\Re\rho$。这是前篇式（8）、（10）与正部的 1-Lipschitz 性的直接推论。因 $\zeta_E=\zeta L(s,\chi_{-4})L_+L_-$ 及有理 Dirichlet GRH，这两个情形涵盖本域 ERH 的两种可能。[Liu 的无零界，定理 7.33](https://hexagonmath.org/pdf/2611.00741v1) · [完整窗口指数](cm229-cubic-hecke-carries.md) · [有理 Dirichlet GRH 的核验](https://math.mit.edu/~drew/grh/)

## 从算术和到平滑显式公式

设 $-\!L'/L=\sum b(q)q^{-s}$，$T_L(x)=\sum_{q\le x}b(q)$，$J(x)=(T_L(2x)-T_L(2x/r))/2$。伴篇的卷积估计 (4) 用到原截断 $D=z^2\le T^{1/3}$，每个二进块的归一化误差为 $O(T^{-1/18}\log^C(2T))$。从 $Y$ 起求和得到

$$
 W^{-1}\sum_n k_n|R_Y(n)-J(n)|\ll Y^{-1/18}\log^C(2Y). \tag{4}
$$

$|b(q)|\le4\Lambda(q)$，故 $|J(x)|\ll x$，在一个单位段内的变动为 $O(\log(2x))$。对上述两种 $f$，以 $\int_Y^X f(J(x))x^{-3/2}dx$ 替换 $\sum_n n^{-3/2}f(J(n))$ 只花 $O(Y^{-1/2}\log(2Y))$，包含两个端点。用 $k_n/W$ 替换 $n^{-3/2}$ 至多再花 $O(W^{-1}X^2)=O(1)$。

取标准实高斯变量 $G$。对任意 $K\ge2$ 置 $J_K(x)=\mathbb E[J(xe^{G/K})]$。若平移的对数为 $w$，对单个端点 $ax$（$a=2$ 或 $2/r$），交换 $q$ 和 $x$ 的求和积分。跨过 $q$ 的 $x$ 对权 $x^{-3/2}dx$ 的费用至多为 $Cq^{-1/2}|w|e^{|w|/2}$，且 $q\le2Xe^{|w|}$。用 $\sum_{q\le t}|b(q)|q^{-1/2}\ll\sqrt t$，即得

$$
 \int_Y^X |J_K(x)-J(x)|x^{-3/2}dx\ll \sqrt X/K. \tag{5}
$$

$J_K$ 的 Mellin 核是 $\tfrac12(2x)^s(1-e^{-vs})e^{s^2/(2K^2)}/s$，它在 $s=0$ 可去。在 $\Re s>1$ 与 $-L'/L$ 相乘后作逆变换，移线至 $-1/2$。非平凡零点给出 (1) 中的负号和重数。在 $s=0$ 可能遇到平凡零点，其贡献为 $O(1)$。新线上由函数方程有 $L'/L\ll\log(2+|\Im s|)$，积分为 $O(x^{-1/2}\log^2(2K))$。移线可在避开零点的高度取极限，单位高度内的零点数是 $O(\log(2+|\gamma|))$，高斯因子保证零点和绝对收敛。因而 $e^{-t/2}J_K(e^t)$ 等于 (1) 中遍历全部非平凡零点的和，加上 $O(e^{-t/2}+e^{-t}\log^2(2K))$。[Lagarias–Odlyzko，Hecke 零点计数与显式公式，第 5、7 节](https://aareyanmanzoor.github.io/assets/articles/lagarias-odlyzko.pdf) · [本域的对数导数与单位高度零点界](cm229-critical-logarithm.md)

现取 $K=\sqrt X$。对 $|\gamma|>H$，把高斯衰减一分为二，一半至多为 $e^{-H^2/(4K^2)}=(2X)^{-2}$。另一半按单位高度将 $\log(2+|\gamma|)/(1+|\gamma|)$ 的相应积分截断在 $O(\log^2(2K))$ 内。即使以 $\sqrt X$ 上界所有的 $e^{(\beta-1/2)t}$，整个高尾的绝对值积分也趋于零。接下来只须处理 (1) 中删掉的左侧零点。

## 短区间对零点的均方权重

每个零点的系数除去 $e^{(\rho-1/2)t}$ 后有绝对值上界

$$
 C\min\{v,(1+|\gamma|)^{-1}\}, \tag{6}
$$

因为 $(1-e^{-v\rho})/\rho=\int_0^v e^{-w\rho}dw$。对 $\beta\le1/2+\lambda$，取任一单位段 $a\le t\le a+1$、$A\le a\le3A/2$。用一个固定的光滑紧支撑函数大于等于该段平移到 $[0,1]$ 后的指示函数。对零点和展开平方积分。两次分部积分将一对高度的核界为 $C(1+|\gamma-\gamma'|)^{-2}$，因平移后实指数在 $[-1/2,1/2]$ 内，常数一致。按整数高度 $j$ 分组，以该可和核作卷积，得到这一单位段上的平方积分至多为

$$
 C e^{2\lambda a}
 \sum_{j\in\mathbb Z}\log^2(2+|j|)
        \min\{v,(1+|j|)^{-1}\}^2
 \ll e^{2\lambda a}v\log^2(2/v). \tag{7}
$$

这一步只用零点的局部总数，允许重零点和彼此任意接近的高度。取平方根后，对单位段求和，$\sum_a e^{\lambda a}\ll X^\lambda\min(1+A,1+1/\lambda)$。由 $\sqrt v\asymp1/A$ 及 (4)、(5) 得 (2)。余下的 $64$ 以下零点均已包含在 (7)。更具体地，这 192 对共轭零点各由 (6) 给出 $O(v)$，合计整窗只花 $O(Av)=O(1/\log Y)$。[低零点证书](verify_cm229.py)

## 不取正部的总和与取正部的下界

**定理 2。** 在本域 ERH 下有

$$
 \left|\frac1W\sum_n k_n R_Y(n)\right|
       \ll \frac{\log^2(2+\log Y)}{\log^2Y}. \tag{8}
$$

对实际的正部，则在两种零点情形下都有

$$
 \liminf_{Y\to\infty}(\log Y)P_R(Y)
    \ \ge\ \frac{80\log2}{81\pi\sqrt2}
    =0.1540868434954109\ldots . \tag{9}
$$

这是窄区间双线性和的下界。平方根以上全部进位素数的正部还有另一个、规模为 $W\log Y$ 的下界。此前得到的等价关系比较的是两者的幂指数。[原始进位素数的对数下界](cm229-prime-carry-drift.md)

**证明。** 先假设 ERH。为在更小的尺度上传递 (4)，在 (5) 中取 $K=X^2$。由定理 1 的绝对值上界、(4) 及 $k_n/W\ge(1-Y^{-3/4})n^{-3/2}$，取整误差现在至多 $O(Y^{-3/4}\log(2+A))$。连续端点和新积分线的误差仍为负幂。于是无论对恒等函数还是对正部，都可用 $e^{-t/2}J_K(e^t)$ 在 $[A,3A/2]$ 的相应积分替换归一化窗口和，误差为 $O(Y^{-1/20})$。

写临界线上的系数为

$$
 a_\gamma(v,K)=-2^{\rho-1}\frac{1-e^{-v\rho}}{\rho}
                       e^{\rho^2/(2K^2)}
$$

并按重数求和。已认证 $L(1/2)\ne0$ 且最低正高度 $\gamma_0>0$。对每个频率，$\left|\int_A^{3A/2}e^{i\gamma t}dt\right|\le2/|\gamma|$。由 (6) 及单位高度计数，将零点在 $|\gamma|=1/v$ 处分开，所得绝对上界为 $O(v\log^2(2/v))$。此即 (8)。

为证 (9)，取最低零点 $\rho_0=1/2+i\gamma_0$。它是单零点，且 $|a_{\gamma_0}(v,K)|/v\to1/\sqrt2$。取次数 80 的归一化 Fejér 核与 $\mathbf1_{\{\cos\theta>0\}}$ 作圆周卷积，记为 $Q$。则 $0\le Q\le1$，其一次复 Fourier 系数为 $\widehat Q(1)=80/(81\pi)$。证书逐一排除了 $j\gamma_0$（$2\le j\le80$）与任何零点高度重合，最小隔离间距为 $13542/2^{24}$。

固定 $0\le w\le1$ 光滑且紧支撑于 $(1,3/2)$。用 $w(t/A)Q(\gamma_0t+\arg a_{\gamma_0})$ 对实值零点和加权，该积分是累计正部的下界。在 $Q$ 的有限 Fourier 展开中只有 $\pm\gamma_0$ 与零点共振，贡献为 $2 A(\int w)\widehat Q(1)|a_{\gamma_0}|$。其余频率使用 $|\widehat w(A(\gamma-j\gamma_0))|\ll_w (1+A|\gamma-j\gamma_0|)^{-2}$。由 (6)、有限间隔及单位高度计数，全部非共振贡献为 $O_w(v/A)$。令 $\int w$ 从下方趋于 $1/2$，用 $A^2v\to\log2$ 即得 (9)。若 ERH 不成立，则前篇的正幂下界已经给出 $(\log Y)P_R(Y)\to\infty$。[零点隔离区间](cm229-zeros64.csv) · [80 个谐波与常数的复算](verify_cm229_type_ii.py)
