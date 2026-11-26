---
title: "主项为 x/48 的素数上界能在哪些区间证明？"
date: "2026-11-26"
abstract: "对判别式 229 的全实三次域与 Q(i) 的复合域，证明双惰性素数平均在任意远的、长度为常数乘 √Y log log log Y 的整个区间上严格低于 x/48。原幂次窗口可用超额正部的一个加权积分检验。给出保持原窗口的已知误差指数，并用两种独立算法核对一亿以内完全分裂素数的逐点平方根误差。"
---

# 主项为 x/48 的素数上界能在哪些区间证明？

设 $F=\mathbb Q(\alpha)$，其中 $f(\alpha)=\alpha^3-4\alpha+1=0$。记 $N$ 为 $f$ 的分裂域，$E=F(i)$，$M=N(i)$。对 $p\nmid458$，同时满足 $p\equiv1\bmod4$ 且 $f\bmod p$ 有三个根，恰好表示 $p$ 在 $M$ 完全分裂。记这些素数的集合为 $S$，写 $\Theta_S(x)=\sum_{p\le x,\,p\in S}\log p$。

令 $A_Y(x)$ 表示按前篇规则抽取辅助素数 $q,d$ 后，对 $S$ 中在两个虚二次域同时惰性的素数取对数权之和的期望。固定 $0<e<1/52$，置 $z=Y^{e/5}$。给定 $Y$ 后，整个窗口都使用同一个分布。前篇利用有理 Dirichlet GRH 得到，对每个 $B>0$ 一致地有

$$
 A_Y(x)=\frac{\Theta_S(x)-\Theta_S(z)}4+O_{e,B}(Y^{-B}),
 \qquad Y\le x\le Y^{3/2}. \tag{1}
$$

具体地，令 $P$ 为不超过 $z$ 的奇素数之积，$U=\exp(3z/4)$，$D=\exp(3z(\log Y)^2/4)$。先在 $(U,2U]$ 抽取 $q\equiv3\bmod8$，再在 $(D,2D]$ 抽取 $d\equiv7\bmod16,\ d\equiv-1\bmod qP$，两步都用素数的对数权归一化。原问题要求对**每个** $\delta>0$，有与 $Y$ 无关的 $C_\delta$ 及任意大的 $Y$，使 $A_Y(x)\le x/48+C_\delta x^{1/2+\delta}$ 在整个 $[Y,Y^{3/2}]$ 上成立。这个命题等价于 $\operatorname{ERH}(E)$。下文给出两个可证明的范围，以及检验原命题所需的一个较弱的算术条件。[原估计、抽样分布与等价性](winpair-inert-transfer.md) · [Dirichlet GRH 的独立核验](https://math.mit.edu/~drew/grh/)

## 任意远的短区间上可以去掉误差项

写 $\ell(x)=\log\log\log x$，只在充分大的 $x$ 使用。存在常数 $h,c>0$ 和任意大的 $Y$，使

$$
 A_Y(x)\ \le\ \frac{x}{48}-c\sqrt Y\,\ell(Y)
 \qquad\bigl(Y\le x\le Y+h\sqrt Y\,\ell(Y)\bigr). \tag{2}
$$

另有任意大的 $Y$，使同样长的整个区间都满足 $A_Y(x)\ge x/48+c\sqrt Y\,\ell(Y)$。两条结论均无需假设 $E$ 的 ERH。特别地，即使允许任意固定的平方根误差常数，也不能在**所有**充分大的 $Y$ 上断言 $A_Y(Y)\le Y/48+C\sqrt Y$。

**平方项的计算。** 令 $\Psi_M(x)=\sum_{N\mathfrak p^j\le x}\log N\mathfrak p$，其中 $j\ge1$。记 $g_p$ 为未分歧素数在 $G=S_3\times C_2$ 中的 Frobenius 共轭类，令 $\mathcal C=\{g\in G:g^2=1\}$。若 $g_p$ 的阶为 $r$，则 $p$ 的欧拉因子在每个 $p^{jr}$ 处给 $\Psi_M$ 贡献 $12\log p$。分开一次项和二次项即得

$$
 \Psi_M(x)=12\Theta_S(x)+12\Theta_{\mathcal C}(\sqrt x)
                +O_M(x^{1/3}\log(2x)).
$$

这里 $\Theta_{\mathcal C}$ 按 $g_p\in\mathcal C$ 计未分歧有理素数。三次及以上幂由 $\sum_{j\ge3,\,2^j\le x}\sum_{p\le x^{1/j}}\log p$ 控制，两个分歧素数只花 $O_M(\log(2x))$。$S_3$ 中单位元和三个换位的平方为一，$C_2$ 的两个元素均符合，故 $|\mathcal C|=8$。Chebotarev 素数定理给出 $\Theta_{\mathcal C}(t)=(2/3)t+o(t)$。由 (1) 且 $\Theta_S(z)\ll z=o(\sqrt Y)$，在原窗口上一致地有

$$
 \Delta_Y(x):=A_Y(x)-\frac{x}{48}
   =\frac{\Psi_M(x)-x}{48}-\frac{\sqrt x}{6}+o(\sqrt x). \tag{3}
$$

式 (3) 的 $-\sqrt x/6$ 是由删除素数平方产生的确定项。决定 $\Delta_Y$ 能在哪些区间保持同号的另一项是 $\Psi_M(x)-x$。[Chebotarev 密度定理与 Frobenius](https://www.jmilne.org/math/CourseNotes/CFT.pdf)

**从一个点到整个短区间。** Landau 关于素理想定理余项的振荡定理给出

$$
 \Psi_M(X)-X+\sum_{\substack{\rho\ \mathrm{real}\\0<\rho<1}}
                 \frac{X^\rho}{\rho}
     =\Omega_\pm\bigl(\sqrt X\,\ell(X)\bigr).
$$

实零点按重数计。若采用只计素理想的版本，补回素理想高次幂花 $O_M(\sqrt X)$，由有理素数的 Chebyshev 上界即得。这里 $\zeta_M=\zeta_{\mathbb Q(\sqrt{229},i)}(L_+L_-)^2$，阿贝尔部分的实非平凡零点只能在 $1/2$。两个非阿贝尔因子在包含整个实临界带的高度 64 矩形内已有完整零点计数，均无实零点。故实零点校正总共是 $O_M(\sqrt X)$。结合平方项公式，存在 $a>0$，使 $D_S(X)=\Theta_S(X)-X/12$ 沿两列无界的 $X$ 分别小于 $-a\sqrt X\ell(X)$ 和大于 $a\sqrt X\ell(X)$。[Landau 原文](https://geodesic.mathdoc.fr/item/MZ_1918__2_167475/) · [Kaczorowski–Staś 式 (1.1)、(1.2)](https://doi.org/10.4064/cm-56-1-185-197) · [本域的分解与实零点证书](cm229-heights-and-zeros.md)

对负偏差点取 $Y=X-h\sqrt X\,\ell(X)$，其中 $0<h<6a$。$\sqrt t\,\ell(t)$ 最终递增，故 (2) 的区间包含在 $[Y,X]$。对其中每个 $x$，单调性给出 $D_S(x)\le D_S(X)+(X-x)/12\le-(a/2)\sqrt X\ell(X)$。取 (1) 并缩小 $c$ 即得 (2)。对正偏差点取 $Y=X$，利用 $D_S(x)\ge D_S(X)-(x-X)/12$，再用 $\Theta_S(z)=o(\sqrt Y)$ 即可。这里 $h,c$ 是存在常数，证明不指定第一个起点。

若保留原问题的误差，短区间还可加长。对每个 $0<\delta<1/2$，都有任意大的 $Y$，使

$$
 A_Y(x)\le\frac{x}{48}+\frac{x^{1/2+\delta}}{48}
 \qquad\bigl(Y\le x\le Y+Y^{1/2+\delta}\bigr). \tag{4}
$$

对上面的负偏差点解 $X=Y+Y^{1/2+\delta}$。在 $[Y,X]$ 上将 $\Theta_S(x)\le\Theta_S(X)$ 代入 (1)，由 $(X-x)/48\le x^{1/2+\delta}/48$ 即得。式 (4) 的右端与原问题精度相同，窗口右端仍小于 $Y^{3/2}$。

## 保持整个幂次窗口时可以达到的指数

采用 Liu 预印本定理 7.33，其中每个数域上的酉 Hecke $L$ 函数在 $\Re s>\theta_G$ 无零，且 $\theta_G<0.8733326374459051924$。取向上舍入的 $b=0.8733326374459052$，则对每个充分大的 $Y$，有

$$
 A_Y(x)=\frac{x}{48}+O_{M,e}(x^b)
 \qquad(Y\le x\le Y^{3/2}). \tag{5}
$$

特别地，$\delta=3/8$ 的原估计已经成立。更一般地，每个 $\delta>\theta_G-1/2$ 都可用。式 (5) 引用 Liu 预印本的无零定理。采用已有 Lean 形式化的 Cureton–Thorner 定理，则相应门槛为 $\delta>15/32$。[Liu 定理 7.33](https://hexagonmath.org/pdf/2611.00741v1) · [Cureton–Thorner 定理及形式化](https://hexagonmath.org/2610.00200v2)

为看清指数的余量，将 Liu 定理用于 $M$ 上的平凡 Hecke 特征，选定 $\theta_G<\sigma_0<b$。在 $\sigma\ge\sigma_0$ 的固定条带上，扣除 $s=1$ 的极点后，$\zeta'_M/\zeta_M$ 由单位高度零点计数及至零点的水平距离界为 $O_{M,\sigma_0}(\log(2+|t|))$。对半整数 $u$ 用高度 $u^2$ 的截断 Perron 公式，系数逐项受 $12\Lambda(n)$ 控制。移线越过 1 给出 $\Psi_M(u)=u+O_{M,\sigma_0}(u^{\sigma_0}\log^2(2u))$。半整数之间用单调性，且 $x^{\sigma_0}\log^2(2x)\ll_b x^b$。代入 (3) 即得 (5)。常数依赖所选的指数间隔，本式不指定数值起点。[截断和移线的原窗口推导](winpair-inert-transfer.md)

## 原问题只须控制超额的正部积分

对 $1/2\le\theta<1$，写 $[v]_+=\max(v,0)$，定义

$$
 P_\theta^\pm(Y)=\int_Y^{Y^{3/2}}
        [\,\pm\Delta_Y(x)\,]_+\,\frac{dx}{x^{1+\theta}} .
$$

**积分判据。** 以下两个条件等价。

1. $L_+L_-$ 在 $\Re s>\theta$ 没有零点。
2. 对每个 $\eta>0$，存在任意大的 $Y$ 使 $P_\theta^+(Y)\le Y^\eta$。

把第二条换成负部 $P_\theta^-$ 也等价。取 $\theta=1/2$，仅正部的积分条件便足以推出 $E$ 的 ERH 和原来的逐点上界。若确有 $\Re\rho>\theta$ 的零点，则存在 $\eta>0$，使**每个**充分远的原窗口同时满足 $P_\theta^+(Y)>Y^\eta$、$P_\theta^-(Y)>Y^\eta$。

**证明。** 设 $L=L_+L_-$，写 $-L'/L=\sum b(n)n^{-s}$，$T(x)=\sum_{n\le x}b(n)$。记 $S_0$ 为在 $\mathbb Q(\sqrt{229},i)$ 完全分裂的未分歧素数。由标准表示在单位元、换位、三循环的迹 $2,0,-1$ 以及模 4 的扭转，$b(p)/\log p=6\mathbf1_S(p)-2\mathbf1_{S_0}(p)$。阿贝尔域的 ERH 给出 $\Theta_{S_0}(x)=x/4+O(\sqrt x\log^2(2x))$。加回高次幂并用 (1)，在原窗口上一致得到

$$
 T(x)=24\Delta_Y(x)+O_{M,e}(\sqrt x\log^2(2x)).
$$

正部运算是 1-Lipschitz 的。上述误差按 $dx/x^{1+\theta}$ 在整窗的积分为 $O(\log^3(2Y))$。故第二条转成 $T$ 的相应积分条件时，只要先将 $\eta$ 换成 $\eta/2$ 即可吸收固定因子和对数损失。

令 $g(u)=T(e^u)$。它可测且在 $u\ge0$ 有 $O(e^{2u})$ 的增长。其 Laplace 变换是 $-L'(s)/(sL(s))$，在 $\Re s>0$ 的极点是 $L$ 的孤立零点，留数为 $-m/\rho$。$L$ 是两个非平凡有限阶 Hecke 特征的整 $L$ 函数之积。前述矩形计数排除了全部实非平凡零点，$\Re s\ge1$ 由 Hecke 非消失性。故 $g$ 和 $-g$ 均满足原 Lean 窗口引理的极点条件。取 $V=\log Y$，在 $[V,3V/2]$ 上取非负上界 $R(u)=[g(u)]_+$，则原引理需要的 $\int R(u)e^{-\theta u}\,du$ 恰是 $[T(x)]_+$ 的加权积分。引理即排除 $\Re\rho>\theta$。换成 $-g$ 得负部结论。[原始极点条件和积分量词](https://github.com/leanprover-community/riemann-hypothesis/blob/9c9f665ed79ce0fa8da61bd803b7a4bf1ad62ed0/src/Development/AutomorphicArithmetic/DedekindLFunctions/LandauZeroFreeCriterion.lean#L57) · [原窗口引理](https://github.com/leanprover-community/riemann-hypothesis/blob/9c9f665ed79ce0fa8da61bd803b7a4bf1ad62ed0/src/Development/AutomorphicArithmetic/DedekindLFunctions/LandauTheorems/landauWindow.lean#L125)

反向若 $L$ 在 $\Re s>\theta$ 无零，阿贝尔因子和 $\zeta_M$ 的分解给出 $\zeta_M$ 同一无零区。依 (5) 的移线证明，对任意小的 $\varepsilon>0$ 有 $\Delta_Y(x)\ll_\varepsilon x^{\theta+\varepsilon}$。选 $3\varepsilon/2<\eta$，便在每个充分远的窗口上得到两个正部积分界。对两条等价命题分别取逆否命题并取两个 $\eta$ 的较小者，即得有越线零点时的结论。已有证书把这样的零点限制在 $|\Im\rho|>64$。证明原上界所缺的正是 $\theta=1/2$ 时在任意远的幂次窗口上控制 $P_{1/2}^+(Y)$。

## 有理同余上界会损失三倍

取任意固定模数 $m$。设 $w$ 是模 $m$ 可逆剩余类上的函数，且除有限多个素数外，有 $w(p\bmod m)\ge\mathbf1_S(p)$。则

$$
 \frac{1}{\varphi(m)}\sum_{a\bmod m}^{*} w(a)
 \ \ge\ \frac1{[M\cap\mathbb Q(\zeta_m):\mathbb Q]}
 \ \ge\ \frac14 .
$$

事实上，对每个在交域上作用平凡的 $a$，合成域中存在限制为 $(1,a)$ 的元素。Chebotarev 定理给出无穷多个 $p\in S$ 落在这个类，故 $w(a)\ge1$。其他可逆类上也由素数存在性得 $w(a)\ge0$。$S_3\times C_2$ 的交换子群为 $A_3\times\{1\}$，所以 $M$ 的最大阿贝尔子域为 $\mathbb Q(\sqrt{229},i)$。取 $m=916$，用这个双二次域的完全分裂指示函数便达到 $1/4$。于是先用固定有理同余权逐点控制 $\mathbf1_S$ 再用 (1)，最好的主项是 $x/16$。要达到 $x/48$，必须在这些阿贝尔分裂素数中辨认剩下的三阶理想类信息。此论证适用于任意有限组固定模数，取其最小公倍数即可。

## 一亿以内的逐点计算

计算针对 (1) 右端的固定域计数。对每个实数 $2\le x\le10^8$，置 $H(x)=(\Theta_S(x)-x/12)/\sqrt x$。两种独立的分裂判据对全部 $2\,880\,503$ 个 $p\equiv1\bmod4,\ p\ne229$ 的候选素数给出相同结论，其中 $479\,605$ 个属于 $S$。128 位 Arb 球算术给出

| 量 | 严格包围 | 位置 |
| --- | --- | --- |
| $\sup H$ | $(0.4619508,\ 0.4619509)$ | $x=873157$，取到 |
| $\inf H$ | $(-2.1541849,\ -2.1541848)$ | $x$ 从左趋于 $55511413$ |

因此全范围有 $\Theta_S(x)/4\le x/48+0.11549\sqrt x$。第一次高于主项在 $p=37$，最后一次跳跃后仍高于主项在 $p=7207729$。最后的过零点是 $12\Theta_S(7207729)=7207871.98766513\ldots$。整个 $[7207872,10^8]$ 上均有 $\Theta_S(x)/4<x/48$。将这些有限数字用于移动的 $A_Y$ 还需要 (1) 的显式起点及误差常数。

第一种算法在 $\mathbb F_p[X]/(f)$ 内检验 $X^p=X$。因 $p\nmid458$，这等价于 $f$ 有三个不同的一次因子。第二种算法枚举实二次域 $k=\mathbb Q(\sqrt{229})$ 的主理想范数。这个域的 Hilbert 类域正是 $N$，相应的二次型判据也见于剑桥 2004 年试题。[剑桥 Part III，第 5(ii) 题](https://www.maths.cam.ac.uk/postgrad/part-iii/files/pastpapers/2004/Paper25.pdf) · [LMFDB 的实二次域数据](https://www.lmfdb.org/NumberField/2.2.229.1)

简述该判据的证明。$f$ 不可约且判别式为素数 229，因此 $\operatorname{Gal}(N/\mathbb Q)=S_3$ 且 $N/k$ 为循环三次扩张。229 处为驯分歧，惯性群为循环群，在 $k$ 上限制非平凡，故由一个换位生成。别处没有有限分歧，且 $N$ 全实，故 $N/k$ 处处不分歧。Minkowski 界为 $\sqrt{229}/2<8$，类群由 3 和 5 上方的素理想生成，因为 2、7 惰性。$(17+\sqrt{229})/2$ 的范数为 15，给出 3 与 5 上方类的关系。$(11+\sqrt{229})/2$ 的范数为 $-27$，且该元素不被 3 整除，故 3 上方素理想的类的阶整除 3。类数至多为 3，不分歧三次扩张的存在给出等号。Artin 互反遂将完全分裂判据化为存在范数为 $\pm p$ 的元素。

单位 $u=(15+\sqrt{229})/2$ 的范数为 $-1$。借助 $u$ 和 $-1$ 可以选取范数为 $p$ 的全正生成元 $\pi$。再乘以 $u^2$ 的整数幂并视需要作共轭，可使第一嵌入满足 $1\le\pi/\sqrt p\le u$。写 $\pi=(m+b\sqrt{229})/2$，由 $u-u^{-1}=15$ 得到完整的有限搜索条件

$$
 m^2-229b^2=4p,\qquad m\equiv b\pmod2,\qquad
 m>0,\quad b\ge0,\quad 229b^2\le225p .
$$

反向由该式立即得到范数为 $p$ 的代数整数。程序先用 Eratosthenes 筛判素，再逐素数比较两个判据。两次跳跃之间 $H(x)=\Theta_S(x)x^{-1/2}-\sqrt x/12$ 严格递减，故只需以向外舍入核验所有跳点、左极限及两个端点。运行下列命令即可复算表中数字，分裂结果或球区间比较有异会报错。[完整核验程序](verify_cm229_split_primes.py)

~~~sh
python -m pip install python-flint==0.9.0
python verify_cm229_split_primes.py
~~~

后续工作给出了 $P(Y)$ 在完整窗口上的精确增长指数，将正部实现为 Hermitian 格的正次数，并用混合 CM 型的 $a$ 数识别双惰性素数。[完整窗口与整值多项式格](cm229-binomial-slopes.md) · [局部长度和两尺度斜坡](cm229-local-length-window.md)
