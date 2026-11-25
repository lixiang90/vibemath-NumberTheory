---
title: "判别式 −23 的全尺度素数估计与最小高斯平滑量"
date: "2026-11-25"
abstract: "给出适用于一切尺度的完全分裂素数平方根估计。所需的最小方差斜率等于偏线零点的一个显式上确界，目前上界小于 9.537×10⁻⁷。对每个固定正斜率，只须在任意远处、长度任意慢地趋于无穷的整数区段检查下侧误差。"
---

# 判别式 −23 的全尺度素数估计与最小高斯平滑量

设 $K=\mathbb Q(\alpha)$，$\alpha^3-\alpha-1=0$，$M$ 是其六次分裂域。给定 $u=\log Y$ 和方差 $v>0$，记

$$
 k_v(t)=\frac{e^{-t^2/(2v)}}{\sqrt{2\pi v}},\qquad
 P_v(u)=\sum_{p\ \text{在 }M\text{ 完全分裂}}\log p\,k_v(u-\log p).
$$

素数均未分歧，每个只计一次。此处允许高斯在 $\log p$ 方向上的宽度 $\sqrt v$ 随尺度增长。这样既能把估计推到任意大的 $Y$，也能明确写出须进一步消去的平滑量。以下采用已核验的 Dirichlet GRH。[原 Dirichlet 定理](https://math.mit.edu/~drew/grh/) · [此域的素数分解与二维因子](cubic23-zeros.md)

## 每个实数尺度上的显式界

**定理 1。** 取 $T=1024$。对每个实数 $u\ge30$ 以及每个满足下式的 $v$，

$$
 v\ \ge\ v_*(u):=\max\left(1,\frac{u+99}{T^2-1}\right),
$$

都有

$$
 \left|P_v(u)-\frac{e^{u+v/2}}6\right|<22\,e^{u/2+v/8}. \tag{1}
$$

对 $u\le1048476$ 可以取固定方差 $v=1$。再往右走，高斯在对数坐标上的标准差 $\sqrt{v_*(u)}$ 渐近等于 $\sqrt{\log Y}/\sqrt{1048575}$。主项和误差的归一化包含高斯的指数矩，分别对应 $e^u$ 和 $e^{u/2}$ 的高斯平均。[高度 1024 的全带证书](cubic23-height-1024.md) · [固定方差下更精细的素数平方偏差](cubic23-gaussian-bias.md)

**证明。** 对 $F=L=\zeta_K/\zeta$ 或 $F=\zeta_k$，其中 $k=\mathbb Q(\sqrt{-23})$，记 $S_F^v$ 为 $-F'/F$ 的系数按同一高斯加权的和。伴篇的显式公式和零点总量估计给出

$$
 \begin{aligned}
 S_F^v(u)&=\mathbf1_{F=\zeta_k}e^{u+v/2}-\mathbf1_{F=L}
          -\sum_\rho e^{\rho u+v\rho^2/2}+R_F^v(u),\\
 |R_F^v(u)|&\le
 \frac{50\sqrt{2\pi/v}+8/v}{2\pi}\,e^{-u/2+v/8},\qquad
 \sum_\rho(4+\gamma^2)^{-1}<
 \begin{cases}5&F=L,\\3&F=\zeta_k.\end{cases}
 \end{aligned} \tag{2}
$$

这里 $\rho=\beta+i\gamma$ 遍历相应函数的非平凡零点并计重数。对任意固定 $v>0$，同一条 $\Re s=-1/2$ 移线都成立。式中余项的系数在 $v\ge1$ 时小于 30。以下将各误差除以 $E=e^{u/2+v/8}$。[移线、余项及零点总量的推导](cubic23-gaussian-bias.md)

对 $\zeta_k$，Dirichlet GRH 以及 $\max_t(4+t^2)e^{-vt^2/2}=4$（$v\ge1$）给出 $|S_{\zeta_k}^v-e^{u+v/2}|/E\le12+30e^{-u}$。对 $L$，证书表明低零点 $|\gamma|\le T$ 具有 $\beta=1/2$、$|\gamma|>5$。它们总共至多贡献 $5(4+25)e^{-25/2}$。其余零点用 $\beta\le1$，总共至多贡献

$$
 5(4+T^2)\exp\big((u-v(T^2-3/4))/2\big)
 \le 5(4+T^2)e^{-99/2-v/8}.
$$

再计入常数和左线余项，即得 $|S_L^v|/E<0.001$。

对素数幂，若 $A_r(x)=\sum_{p^j\le x,\,j\ge r}\log p\le c_r x^{1/r}$，则有对任意 $v>0$ 有效的估计

$$
 \sum_p\sum_{j\ge r}\log p\, k_v(u-j\log p)
 \le c_r\left(\frac1r+\sqrt{\frac2{\pi v}}\right)
             e^{u/r+v/(2r^2)}. \tag{3}
$$

对 $A_r(e^w)$ 分部积分，取 $|k_v'(u-w)|$ 后，得到系数 $v^{-1/2}\mathbb E(|Z|e^{-\sqrt v Z/r})$，其中 $Z$ 是标准正态变量。平移正态密度后用 $\mathbb E|Z-\sqrt v/r|\le\sqrt{2/\pi}+\sqrt v/r$ 即得 (3)。伴篇已证 $c_2=15$ 可取。

从 $\zeta_M=\zeta_k L^2$ 提取未分歧素数一次项，有 $6P_v=S_{\zeta_k}^v+2S_L^v$ 减去 23 的一次项和全部高次素数幂。23 的一次项小于 4。高次幂系数至多为 $6\log p$，故 (1) 左侧除以 $E$ 后至多为

$$
 \frac{12+30e^{-30}+2(0.001)+4e^{-15}}6
      +15\left(\frac12+\sqrt{2/\pi}\right)<22.
$$

## 最小斜率与剩余零点恰好等价

固定 $\lambda>0$，取 $v=\lambda u$，置 $\sigma_\lambda=1/2+\lambda/8$。令 $\mathcal Z$ 为 $L$ 的非平凡零点多重集，定义

$$
 \lambda_{23}=
 \sup\!\left(\{0\}\cup
 \left\{\frac{2\beta-1}{\gamma^2-\beta^2+1/4}:
        \beta+i\gamma\in\mathcal Z,\ \beta>1/2\right\}\right). \tag{4}
$$

高度证书保证式中的分母严格为正。任选正整数值函数 $\ell(N)\to\infty$，例如 $\ell(N)=\lceil\log\log(N+3)\rceil$。

**定理 2。** 对每个固定的 $\lambda>0$，以下三项等价。

1. $\lambda\ge\lambda_{23}$。
2. 对所有充分大的实数 $u$，有 $\left|P_{\lambda u}(u)-e^{u+\lambda u/2}/6\right|=O_\lambda(e^{\sigma_\lambda u})$。
3. 对每个 $\varepsilon>0$ 都有任意大的整数 $N$，使每个 $N\le n\le N+\ell(N)$ 的整数都满足
   $P_{\lambda n}(n)\ge e^{n+\lambda n/2}/6-e^{(\sigma_\lambda+\varepsilon)n}$。

所以 $\zeta_K$ 的 ERH 等价于 $\lambda_{23}=0$，也等价于条件 3 对一列趋于零的正数 $\lambda$ 成立。每个斜率可以选自己的窗口，常数与窗口起点都允许依赖斜率。已验证的完整零点计数无条件给出关于 $L$ 的谱界

$$
 0\le\lambda_{23}\le\frac1{1024^2-3/4}<9.537\cdot10^{-7}. \tag{5}
$$

若已有独立的全高度实部上界 $\beta\le b$，$1/2<b\le1$，则 (5) 的右端可改为 $(2b-1)/(1024^2-b^2+1/4)$。用 Cureton–Thorner 的 $b=31/32$ 得 $\lambda_{23}<8.941\cdot10^{-7}$。若再采用 Liu 预印本定理 1.1 的 $b=(777-\sqrt{1081})/852$，上界为 $7.122\cdot10^{-7}$。式 (5) 与 (1) 的常数略有不同，因为 (1) 还在每个尺度明确控制绝对误差常数。[Cureton–Thorner 定理及其形式化](https://github.com/AxiomMath/DedekindZetaZeroFree) · [Liu 对有限阶 Hecke 特征的预印本](https://hexagonmath.org/pdf/2611.00082v1)

### 有限组相位如何强迫素数短缺

取 $\Phi_\lambda(s)=s+\lambda s^2/2$。从 (2)、(3) 和六次闭包的分解，对充分大的 $u$ 得

$$
 \begin{aligned}
 6P_{\lambda u}(u)-e^{u+\lambda u/2}
   &=2S_L^{\lambda u}(u)+O(e^{\sigma_\lambda u}),\\
 S_L^{\lambda u}(u)
   &=-1-\sum_{\rho\in\mathcal Z} e^{u\Phi_\lambda(\rho)}
               +O(e^{(-1/2+\lambda/8)u}).
 \end{aligned} \tag{6}
$$

实部满足 $\Re\Phi_\lambda(\rho)=\beta+\lambda(\beta^2-\gamma^2)/2$，随 $|\gamma|\to\infty$ 一致趋于负无穷。条件 1 正好说每个零点都有 $\Re\Phi_\lambda(\rho)\le\sigma_\lambda$。对于 $\beta\le1/2$，这个不等式直接成立。对其余零点整理不等式即得 (4)。若条件 1 成立，在任一固定 $u_0>0$，(2) 保证 $\sum_\rho e^{u_0(\Re\Phi_\lambda(\rho)-\sigma_\lambda)}<\infty$。对所有 $u\ge u_0$，该常数控制归一化零点和，于是 (6) 给出条件 2。条件 2 蕴涵条件 3。

若条件 1 失败，存在最大实部 $B=\max_\rho\Re\Phi_\lambda(\rho)>\sigma_\lambda$。只有有限多个零点达到 $B$，记其多重集为 $J$。其余零点的实部至多为 $B-\delta$，其中 $\delta>0$。这些结论来自 $\Re\Phi_\lambda(\rho)\to-\infty$ 及原零点的离散性。由同一绝对收敛估计，在整数 $n\to\infty$ 时有

$$
 e^{-Bn}S_L^{\lambda n}(n)
   =-\sum_{\rho\in J}e^{in\Im\Phi_\lambda(\rho)}+o(1). \tag{7}
$$

存在一个固定整数 $m$，使每个长度为 $m$ 的连续整数区段都含有一个 $n$，对全部 $\rho\in J$ 同时有 $\Re(e^{in\Im\Phi_\lambda(\rho)})>3/4$。为见此，令 $x=(e^{i\Im\Phi_\lambda(\rho)})_{\rho\in J}$，取它在有限维圆环中生成的闭子群 $H$。正整数幂在 $H$ 中稠密，因为紧性给出趋于单位元的任意大正整数幂。令 $W$ 为各坐标实部都大于 $3/4$ 的单位元邻域，则 $x^{-j}W$（$j\ge0$）覆盖 $H$。取有限子覆盖，令 $m$ 大于这些 $j$ 的最大值，作用于 $x^N$ 即得所需整数。

这些 $n$ 上，(7) 的负号和零点重数使 $S_L^{\lambda n}(n)\le-|J|e^{Bn}/2$ 对充分大的 $n$ 成立。由 (6)，完全分裂素数的短缺至少为 $c e^{Bn}$，其中 $c>0$。取 $0<\varepsilon<B-\sigma_\lambda$，则每个充分靠后的长度为 $m$ 的区段都违反条件 3 的下界。一旦 $\ell(N)\ge m$ 就与条件 3 矛盾。最后，$\lambda_{23}=0$ 恰好排除全部 $\beta>1/2$ 的零点。函数方程和 $\zeta_K=\zeta L$ 给出所述 ERH 等价性。

原窗口引理通过 Laplace 变换中的极点，处理固定权和任意大的潜在零点高度。上面的变方差让高处零点的指数实部多出 $-\lambda\gamma^2/2$。于是在每个固定正斜率下，极大的那组零点有限，可以直接在整数尺度上让它们同相，窗口长度只须趋于无穷。将斜率压到任意小仍要求对这个固定域取得新的素数下界。有限高度计算给出了 (5) 的数值上界，式 (4) 则说明任何真正的偏线零点怎样阻止它降到零。[原窗口定理](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/DedekindLFunctions/LandauTheorems/landauWindow.lean) · [固定方差的离散窗口判据](cubic23-sparse-windows.md) · [不等式常数复算](gaussian23_bounds.py)
