---
title: "在 ζ 零点对的过渡尺度上分离固定差 2"
date: "2026-11-24"
abstract: "在 ζ-RH 下，构造三频带的带符号二次统计，使它等于固定差 2 的光滑素数幂相关，绝对误差为 O(X⁻²)。给出剔除高次幂后的单对精度判据及有限算例。"
---

# 在 ζ 零点对的过渡尺度上分离固定差 2

假设 Riemann ζ 函数满足 RH。令零点 $1/2+i\gamma_r$ 遍历正、负高度，按重数编号。以下给出一个不同编号零点的带符号二次统计。它在 $X=T$ 的尺度上，等于差恰好为 2 的加权素数幂对之和，绝对误差为 $O(X^{-2})$。这比单个平衡孪生素数对在窗口中心的贡献 $\asymp\log^2 X/X$ 小。本文由 von Mangoldt 显式公式构造三个频带，使整数差筛选、极点误差和同编号零点项可以分别计算。Montgomery 的对关联和固定移位的素数对渐近为此统计提供了两种比较标准。[Montgomery](https://public.websites.umich.edu/hlm/paircor1.pdf) · [Gou–Wang–Xu，第 4 节](https://arxiv.org/abs/2611.18825)

## 1. 要测试的零点和

固定 $b=1/32$。取非负实偶的 $f\in C_c^\infty((-b,b))$，$f(0)>0$，写 $\widehat f(t)=\int f(v)e^{itv}\,dv$。取实函数 $\phi$ 支撑于 $[7/8,9/8]$，在 $[15/16,17/16]$ 上恒等于 1，$0\le\phi\le1$，可先取光滑。定义 $k_a(s)=\phi(|s|/a)$ 和

$$
K(s)=k_2(s)-7k_{1/2}(s)+6k_{1/4}(s),\qquad
W(t)=\frac1{2\pi}\int_{\mathbb R}K(s)e^{-its}\,ds .
$$

对每个零点置 $z_r(t)=X^{i\gamma_r}\widehat f(\gamma_r-t)$。定义

$$
\mathcal Q_2^\circ(X)=\frac1X\int_{\mathbb R}W(t/X)
\sum_{r\ne s}z_r(t)\overline{z_s(t)}\,dt . \tag{1}
$$

不同编号的重零点全部保留，和为实数。单位区间上零点数 $\ll\log(2+|t|)$，加上 $\widehat f$ 的快速衰减，使双和与积分绝对可交换。对每个 $X$，(1) 是 ζ 零点的函数。高度权是平滑的，$X=T$ 对应标准参数 $\alpha=\log X/\log T=1$。

**固定差公式。** 记 $v_n=\log(n/X)$。当 $X\to\infty$ 时，有

$$
\boxed{\displaystyle
\mathcal Q_2^\circ(X)=2\sum_{n\ge1}
\frac{\Lambda(n)\Lambda(n+2)}{\sqrt{n(n+2)}}f(v_n)f(v_{n+2})
+O_{f,\phi}(X^{-2}).} \tag{2}
$$

例如第 5 节的具体权给出 $0.05325<\mathcal Q_2^\circ(10^6)<0.05338$。积中含真正的素数以及高次素数幂，二者在第 3 节分开。

**需要知道多高的零点。** 对光滑 $f,\phi$，任取固定 $\varepsilon,A>0$，只计 $|\gamma_r|,|\gamma_s|\le X^{1+\varepsilon}$，产生的额外误差为 $O_{f,\phi,\varepsilon,A}(X^{-A})$。证明可把 $t$ 分为 $|t|\le L/2$ 和 $|t|>L/2$，$L=X^{1+\varepsilon}$。前者中，高度超过 $L$ 的 $\widehat f(\gamma-t)$ 绝对质量由任意负幂 $O_j(L^{-j})$ 控制。后者中，全部零点的绝对双和至多为 $O_f(\log^2(2+|t|))$，而 $W(t/X)$ 在 $|t|/X\gg X^\varepsilon$ 具有任意次幂衰减。

选取具体光滑权还可以将高度收紧为 $L=C_A X\log X\,(\log\log X)^2$。用下引论文中以均匀分布无穷卷积构造的紧支概率密度，将其缩放至 $[-1/32,1/32]$，与 $[1-3/32,1+3/32]$ 的指示函数作卷积，所得 $\phi$ 正好有上述支撑和全高平台。该密度的 Fourier 变换满足 $\ll\exp(-c|t|/\log^2(2+|t|))$，缩放、卷积和三项线性组合后 $W$ 也满足这个界。代入上一段的第二个积分并选大 $C_A$ 即得。取 $A=3$ 便保留 (2) 的精度。第 5 节数表另选分段多项式权，并通过算术侧包围完整的零点和。[紧支密度及衰减证明，第 5 节](../Goldbach_TwinPrimes/goldbach-low-zeros.md)

## 2. 三个权分别消去了什么

**选出整数差。** 令 $a_n=\Lambda(n)f(v_n)/\sqrt n$，$D_X(t)=\sum_n a_ne^{-itv_n}$。Fourier 反演给出

$$
\frac1X\int W(t/X)|D_X(t)|^2dt
=\sum_{m,n} a_ma_n K(X\log(m/n)).
$$

在窗口内，由中值定理，$X|\log(m/n)|$ 落在 $|m-n|[e^{-b},e^b]$ 内。差为 2 时，$k_2=1$。其余整数差都不在 $k_2$ 的支撑内。最小非零差的下边界又大于 $9/16$，所以两个小频带没有任何整数对落入。于是这里右边**精确**等于 (2) 的素数幂和，因子 2 来自两个次序。

**控制极点和平凡零点。** 由于 $\int_0^\infty k_a(s)\,ds$ 与 $a$ 成正比，而 $\int_0^\infty k_a(s)\,ds/s$ 与 $a$ 无关，所选系数满足

$$
\int_0^\infty K(s)\,ds=0,\qquad \int_0^\infty K(s)\,\frac{ds}{s}=0 .
\tag{3}
$$

由第一式及偶性，$W(0)=W'(0)=0$。对经典显式公式在 $u=\log x>0$ 上作分布微分，再乘 $e^{-u/2}$，得到
$$
\sum_r e^{i\gamma_r u}=J_\zeta(u)-
 \sum_{n\ge2}\frac{\Lambda(n)}{\sqrt n}\,\delta_{\log n}(u).
$$

这里每个测试函数紧支于开半轴。极点贡献 $e^{u/2}$，平凡零点 $-2,-4,\ldots$ 给出负的几何级数。以 $f(u-\log X)e^{-it(u-\log X)}$ 测试即得

$$
Z_X(t):=\sum_r z_r(t)=B_X(t)-D_X(t),\qquad
B_X(t)=\int f(v)J_\zeta(\log X+v)e^{-itv}dv,\qquad
J_\zeta(y)=e^{y/2}-\frac{e^{-5y/2}}{1-e^{-2y}} .
$$

因为 $|W(t/X)|\ll(t/X)^2$，$|B_X(t)|\ll_{f,j}\sqrt X(1+|t|)^{-j}$，且 $|Z_X(t)|\ll_f\log(2+|t|)$，以 $|Z_X|^2$ 替换 $|D_X|^2$ 在 $X^{-1}\int W(t/X)\,dt$ 下产生的误差为 $O(X^{-2})$。这里保留了 ζ 在 1 处极点的贡献再估计其大小。

**扣除同一编号的零点。** 写 $g(u)=\int f(v)f(v-u)\,dv$，$I_f=g(0)$。相同编号在零点二次型中的贡献可单独计算

$$
\mathcal D_X=\frac1X\int W(t/X)\sum_r|\widehat f(\gamma_r-t)|^2dt
=\frac2X\int_0^\infty g(s/X)K(s)J_\zeta(s/X)\,ds . \tag{4}
$$

第一步用 $|\widehat f(y)|^2=\int g(u)e^{iyu}du$。再用上述分布显式公式测试 $g(u)K(Xu)$。这个测试远离 $u=0$，其正半支撑在 $(0,\log 2)$ 内（$X$ 充分大），其中没有 $\log n$ 处的素数幂。负半轴由 ζ 零点的共轭对称性给出相等的一份。这个步骤消除了高处零点密度近似所带来的误差。

当 $u\downarrow0$ 时，$J_\zeta(u)=-1/(2u)+7/4+O(u)$，而 $g(u)=I_f+O_f(u^2)$。代回 (4)，常数级项正比于 $\int K(s)\,ds/s$，$X^{-1}$ 级项正比于 $\int K(s)\,ds$。两者都由 (3) 消失，剩余为 $O(X^{-2})$。从完整二次型扣掉 $\mathcal D_X$，即证 (2)。

## 3. 这给孪生素数一个多强的判据

令 $C_{\rm pow}(X)$ 为 (2) 右边主和中至少有一项是高次素数幂的部分，含相同的系数 2。窗口内至多有 $O(\sqrt X)$ 个高次素数幂，每个对固定差至多有两个位置，所以 $0\le C_{\rm pow}(X)\ll_f X^{-1/2}\log^2 X$。

- 若 $\limsup_{X\to\infty}\mathcal Q_2^\circ(X)>0$，则在无穷多个固定比例的窗口内有 $\gg_f X/\log^2 X$ 对孪生素数。
- 在单个素数对的精度下，可以得到等价式

$$
\text{有无穷多对孪生素数}\quad\Longleftrightarrow\quad
\limsup_{X\to\infty}\frac{X}{\log^2 X}
\bigl(\mathcal Q_2^\circ(X)-C_{\rm pow}(X)\bigr)>0 . \tag{5}
$$

有无穷多对时，沿 $X=p$ 取 $p,p+2$，一对的归一化贡献趋于 $2f(0)^2>0$，其他真素数对贡献非负。若仅有有限对，则窗口最终离开它们，(2) 的余项归一化后趋零。这里写出 $C_{\rm pow}$ 很重要，它的粗上界大于单对素数的信号。(2) 和 (5) 尚未给出右边所需的无穷次正下界。

如希望所有输入都来自零点，可定义 $\mathcal L(B)=\sum_m\Lambda(m)B(m)$，由第 2 节的分布显式公式计算每个支撑于 $(1,\infty)$ 的光滑紧支测试 $B$（在 $u$ 变量取测试 $e^{u/2}B(e^u)$）。在本篇的窗口内有有限 Möbius 反演

$$
\sum_p(\log p)B(p)=
\sum_{k\le\log(Xe^b)/\log2}\mu(k)\,\mathcal L\bigl(B(\,\cdot\,^k)\bigr).
$$

这是因为 $p^a$ 的总系数为 $\log p\sum_{k\mid a}\mu(k)$。取 $B(n)=f(v_n)e^{-itv_n}/\sqrt n$ 后，可分别对纯素数多项式和 $D_X$ 作相同的 $W$ 积分。两者之差恰好给出 $C_{\rm pow}$（以素数幂版本减纯素数版本）。这样获得完全用显式公式表示的判据，仍需对随 $X$ 变化的值取得下界。

若另假定固定差 2 的 Hardy–Littlewood 渐近式 $\sum_{n\le x}\Lambda(n)\Lambda(n+2)\sim\mathfrak S(2)x$，分部求和给 $\mathcal Q_2^\circ(X)\to2\mathfrak S(2)I_f$，其中 $\mathfrak S(2)=2\prod_{p>2}(1-1/(p-1)^2)$。这里是从算术预测推出下阶零点信号。Chan 已从定量素数对假设研究过普通形状因子的较低阶项。[Chan 的原文](https://arxiv.org/abs/math/0206293)

## 4. 为什么 GUE 的主项读不出差 2

在具有恒定平均密度的平稳正弦核点过程中，把所有点同时平移 $t$，联合分布不变，而 $X^{i(\gamma_r-\gamma_s)}$ 只依赖差。因此 (1) 被积二次和的期望与 $t$ 无关。又 $\int W(t/X)\,dt/X=K(0)=0$，所以这个模型对 $\mathcal Q_2^\circ(X)$ 的期望**精确为零**。这对使局部密度取 $\log X/(2\pi)$ 的平稳模型同样成立。真实 ζ 测试还保留随高度变化的信息。它可以在常数级看到素数差，而普通 Montgomery 主项的归一化要除以 $\log X$。因此 (5) 要求的是主项极限以外的精度。[Montgomery 关于随机矩阵对应的原文](https://public.websites.umich.edu/hlm/paircor1.pdf)

同样的构造能分辨其他固定整数差 $d\ge1$。把 $\phi$ 及 $b$ 选得足够窄，使第 $d$ 个整数带与其邻带分离，再取

$$
K_{[d]}=k_d+(1-4d)k_{1/2}+(4d-2)k_{1/4}.
$$

这保留 (3) 和 (2) 的证明。对于固定奇差，其中一端必须是 2 的幂，故谱统计为 $O_{f,d}(\log X/X)$。对于偶差，利用 Polymath8b 的有界间隙定理，至少存在一个固定偶数 $2\le d\le246$ 使 (5) 的相应上极限为正。这里先以抽屉原理固定 $d$，再取以相应素数为中心的窗口。选定 $d=2$ 时正性仍是开放问题。[Polymath8b](https://arxiv.org/abs/1407.4897)

## 5. 有限高度尺度上的数值包围

数值例取 $f(v)=(1-(v/b)^2)_+^4$（即在 $|v|<b$ 取四次幂，在其余处取零），$\phi$ 在支撑边缘与平台间线性过渡。此 $f$ 为 $C^3$，所得 $W(t)=O((1+|t|)^{-2})$，零点积分仍然绝对收敛。(4) 可先平滑 $\phi$、再取极限，或直接用连续紧支的分段光滑核，所用支撑离 $0,\log2$ 有正距离。以下有限误差估计只使用 $f$ 的两阶导数。

对上述 $f,\phi$ 和 $d=1,2$，若 $X\ge1000$，(2) 的误差绝对值小于

$$
4\,000\,000\,\frac{\log(2X)}{X^2}. \tag{6}
$$

为复核这个常数，置 $h_X(v)=f(v)J_\zeta(\log X+v)$，以及有符号测度 $\nu_X=h_X(v)\,dv-\sum_n a_n\delta_{v_n}$。它的 Fourier 变换（取负号）为 $Z_X$。完整的零点二次型因此等于 $\iint K(X(v-w))\,d\nu_X(v)d\nu_X(w)$。记 $L_X(v)=\int K(X(v-w))h_X(w)dw$。由 (3) 与偶性，Taylor 公式给出

$$
\|L_X\|_\infty\le\frac{\|h_X''\|_\infty}{2X^3}\int s^2|K(s)|ds .
$$

具体有 $\|f'\|_\infty\le256,\ \|f''\|_\infty\le57344$。在这扇窗上 $J_\zeta$ 及其前两阶导数的绝对值均至多 $2\sqrt X$，可逐项求导平凡零点的几何级数（$Xe^v>90$）。从而 $\|h_X''\|_\infty<116000\sqrt X$，$\|h_X\|_1\le\sqrt X$，且直接数整数得 $\sum |a_n|\le\sqrt X\log(2X)$。两种 $K$ 均满足 $\int_{\mathbb R}|K|\le3.5$，$\int_{\mathbb R}s^2|K|<18$。完整二次型与素数幂主和之差因此至多为 $(\|h_X\|_1+2\sum|a_n|)\|L_X\|_\infty$。

对 (4) 的误差，$|g(u)-I_f|\le u^2\|f'\|_2^2/2\le2048u^2$，$I_f\le1/16$。在 $0<u\le0.1$ 有 $|J_\zeta(u)+1/(2u)-7/4|\le3u$。例如将第二项写成 $e^{-3u/2}/(2\sinh u)$，用 $0\le1-u/\sinh u\le u^2$ 和指数函数的二阶 Taylor 界即得。这也给出 $|J_\zeta(u)|\le1/u$。于是 (4) 扣掉前两项后，绝对值至多 $2(2049)X^{-2}\int_0^\infty s|K(s)|ds<17000/X^2$。合并两种误差即得 (6)。

| $X$ | 差 $d$ | 窗内的无序素数幂对 | 其中两端都为素数 | 真实 ζ 零点统计的严格包围 |
| ---: | ---: | ---: | ---: | --- |
| $10^5$ | 1 | 0 | 0 | $\lvert\mathcal Q_1^\circ\rvert<0.004883$ |
| $10^5$ | 2 | 67 | 66 | $0.03642<\mathcal Q_2^\circ<0.04620$ |
| $10^6$ | 1 | 0 | 0 | $\lvert\mathcal Q_1^\circ\rvert<0.000059$ |
| $10^6$ | 2 | 451 | 450 | $0.05325<\mathcal Q_2^\circ<0.05338$ |

这些区间由精确的素数幂枚举、96 比特 Arb 球算术和 (6) 导出。$X$ 是频率及平滑高度尺度，(1) 中所有纵坐标都计入。代码最后还在 $X=4096$ 的整个整数窗口上逐对检查两个 Fourier 权是否只选出指定的差。

```python
from math import isqrt
from flint import arb, ctx

ctx.prec = 96
def prime_powers(stop):
    flags = bytearray([1]) * (stop+1)
    flags[:2] = b"\x00\x00"
    for p in range(2, isqrt(stop)+1):
        if flags[p]:
            flags[p*p::p] = bytes(len(range(p*p, stop+1, p)))
    base = {}
    for p in range(2, stop+1):
        if flags[p]:
            n = p
            while n <= stop:
                base[n] = p
                n *= p
    return base

xs = (100_000, 1_000_000)
base = prime_powers(104*max(xs)//100+2)
b = arb(1)/32
for x in xs:
    weight = {}
    for n in range(96*x//100, 104*x//100+1):
        if n not in base:
            continue
        v = (arb(n)/x).log()
        if abs(v)>b:
            continue
        assert abs(v)<b
        weight[n] = arb(base[n]).log()*(1-(v/b)**2)**4/arb(n).sqrt()
    err = 4_000_000*arb(2*x).log()/arb(x)**2
    for d in (1,2):
        pairs = [n for n in weight if n+d in weight]
        prime_count = sum(base[n]==n and base[n+d]==n+d for n in pairs)
        mass = 2*sum((weight[n]*weight[n+d] for n in pairs), arb(0))
        lower, upper = {
            (100_000, 1): ("-0.004883", "0.004883"),
            (100_000, 2): ("0.03642", "0.04620"),
            (1_000_000, 1): ("-0.000059", "0.000059"),
            (1_000_000, 2): ("0.05325", "0.05338"),
        }[x, d]
        assert arb(lower) < mass-err and mass+err < arb(upper)
        print(x, d, len(pairs), prime_count, mass, err)

from fractions import Fraction

def band_value(y):
    if y < arb(7)/8 or y > arb(9)/8:
        return 0
    assert arb(15)/16 < y < arb(17)/16
    return 1

x = 4096
logs = {n: (arb(n)/x).log() for n in range(3950, 4240)}
logs = {n: v for n, v in logs.items() if abs(v) < b}
for d in (1, 2):
    coeffs = (1, 1-4*d, 4*d-2)
    scales = (d, Fraction(1, 2), Fraction(1, 4))
    assert sum(coeffs)==0 and sum(c*s for c, s in zip(coeffs, scales))==0
    for n, v in logs.items():
        for m, w in logs.items():
            arg = x*abs(v-w)
            actual = (band_value(arg/d) + coeffs[1]*band_value(2*arg)
                      + coeffs[2]*band_value(4*arg))
            assert actual == int(abs(n-m)==d)
print("integer-shift and moment checks passed")
```
