---
title: "完整窗口上的正部积分与整值多项式格的正次数"
date: "2026-11-27"
abstract: "在两个固定数域上，用整值多项式从 n 次到 2n 次的乘法映射构造可复算的有限余核。其对数指数之差给出全窗口上的离散正部和，与 x/48 积分有相同的精确增长指数。把它实现为一个多项式秩 Hermitian 格的正次数，构造一个全局整映射并计算现有几何界的损失。"
---

# 完整窗口上的正部积分与整值多项式格的正次数

令 $F=\mathbb Q[X]/(X^3-4X+1)$，$N$ 为其正规闭包，$M=N(i)$，$B=\mathbb Q(\sqrt{229},i)$。已知 $\zeta_M/\zeta_B=L^2$，其中 $L=L_+L_-$ 是六次 CM 域 $F(i)$ 的两个非阿贝尔整因子之积。写 $\beta_*=\sup_{L(\rho)=0,\ 0<\Re\rho<1}\Re\rho$。由函数方程及已隔离的零点，$\beta_*\ge1/2$。本文使用已核验的有理 Dirichlet GRH。[数域分解与零点证书](cm229-heights-and-zeros.md) · [上游定理的独立核验](https://math.mit.edu/~drew/grh/)

记 $A_Y$ 为前篇的双惰性素数对数权平均，$P(Y)=\int_Y^{Y^{3/2}}[A_Y(x)-x/48]_+\,dx/x^{3/2}$。下文构造一个由格映射给出的离散量 $\mathcal B(Y)$，并证明

$$
 \lim_{Y\to\infty}\frac{\log(1+P(Y))}{\log Y}
 =\lim_{Y\to\infty}\frac{\log(1+\mathcal B(Y))}{\log Y}
 =\frac32\bigl(\beta_*-\tfrac12\bigr). \tag{1}
$$

两个极限都沿全部实数 $Y$ 取。若有临界线右侧零点，每个充分远的完整窗口都有同一正的极限指数。若 $\beta_*=1/2$，则每个 $\eta>0$ 都能取到。将 Liu 预印本定理 7.33 的上界代入 (1)，原积分对每个 $\eta>0.5599989561688578$ 都在所有充分大的窗口小于 $Y^\eta$。[原始正部积分判据](cm229-one-sided-primes.md) · [Liu 定理 7.33](https://hexagonmath.org/pdf/2611.00741v1)

## 一个由乘法给出的全局整映射

对固定数域 $K$，令 $V_{K,n}$ 为次数至多 $n$ 且在 $\mathcal O_K$ 上取整值的 $K$ 系数多项式。系数 $[X^n]$ 将 $V_{K,n}/V_{K,n-1}$ 标识为可逆分式理想 $I_{K,n}$。整数值多项式相乘仍取整值，因而有自然单射

$$
 \mu_{K,n}: I_{K,n}\otimes_{\mathcal O_K} I_{K,n}
                     \longrightarrow I_{K,2n}.
$$

写 $\mathfrak f_K(n)=I_{K,n}^{-1}$。在 $\mathfrak p$ 的剩余域大小为 $q$ 时，阶乘理想满足 $v_\mathfrak p(\mathfrak f_K(n))=\sum_{j\ge1}\lfloor n/q^j\rfloor$。按均匀的 $q$ 进制剩余类顺序列点，从 $a_0$ 开始，在 $\prod_{i<n}(a_n-a_i)$ 中模 $\mathfrak p^j$ 的同余恰贡献 $\lfloor n/q^j\rfloor$ 个因子。局部整值多项式的逐次 Newton 基给出全部首项系数，局部化再拼合理想即得全局公式。[Bhargava 的素理想排序定理](https://collaborate.princeton.edu/en/publications/p-orderings-and-polynomial-functions-on-arbitrary-subsets-of-dede/) · [Byszewski–Frączyk–Szumowicz 的数域阶乘理想](https://arxiv.org/abs/1506.02696)

秩一映射 $\mu_{K,n}$ 的有限余核有 Fitting 理想

$$
 \mathfrak b_K(n)=\mathfrak f_K(2n)\mathfrak f_K(n)^{-2},\qquad
 \#\operatorname{coker}\mu_{K,n}=\mathrm N_K\mathfrak b_K(n).
$$

每个局部指数是 $\sum_{j\ge1}(\lfloor 2n/q^j\rfloor-2\lfloor n/q^j\rfloor)$，各项均为 0 或 1。在 $K=\mathbb Q$ 时，余核的阶恰为 $\binom{2n}{n}$。对本题的两个固定数域，定义

$$
 h(n)=\frac12\log\frac{\#\operatorname{coker}\mu_{M,n}}
                            {\#\operatorname{coker}\mu_{B,n}},\qquad
 \mathcal B(Y)=\sum_{Y\le n<Y^{3/2}}\frac{[h(n)]_+}{n^{3/2}} . \tag{2}
$$

这里大小关系由两个自然格的指数决定。每个域的阶乘渐近式中的线性常数都在本域 $2n$ 与两个 $n$ 的比较中消掉。[数域广义阶乘的 Dirichlet 级数与渐近式](https://arxiv.org/abs/2603.13720)

## 传给零点时损失多少

设 $-L'/L(s)=\sum b(m)m^{-s}$，$a(m)=\sum_{j\mid m}b(j)$，$D(x)=\sum_{m\le x}a(m)$。由局部阶乘公式及 $\zeta_M/\zeta_B=L^2$，

$$
 D(x)=\frac12\log
 \frac{\mathrm N_M\mathfrak f_M(\lfloor x\rfloor)}
      {\mathrm N_B\mathfrak f_B(\lfloor x\rfloor)},\qquad
 H(x)=D(2x)-2D(x),\qquad H(n)=h(n).
$$

$b(1)=0$。在 $\Re s>1$ 逐项积分，可得

$$
 \int_1^\infty H(x)x^{-s-1}\,dx
   =\frac{(2^s-2)\zeta(s)}s\left(-\frac{L'}L(s)\right). \tag{3}
$$

$2^s-2$ 消去 $\zeta$ 在 1 的极点，在开临界带内无零。由有理 RH，$\zeta$ 在 $\Re s>1/2$ 无零，因此 $L$ 在右半临界带的每个零点都在 (3) 留下极点。两个因子的高度 64 完整计数排除了全部实非平凡零点，(3) 在正实轴上无极点。它满足原 Lean 窗口引理的极点条件。$|b(m)|\le4\Lambda(m)$ 可由两个二维表示的局部迹得到，进而 $|a(m)|\le4\log m$。另一方面，地板差只有 0 和 1，由 $\sum_{m\le 2x}|b(m)|\ll x$ 得 $H(x)\ll x$。[原极点与增长条件](https://github.com/leanprover-community/riemann-hypothesis/blob/9c9f665ed79ce0fa8da61bd803b7a4bf1ad62ed0/src/Development/AutomorphicArithmetic/DedekindLFunctions/LandauZeroFreeCriterion.lean#L57)

在 $[n,n+1)$ 的前半段 $H(x)=h(n)$，后半段多出 $a(2n+1)$。因此以 $n^{-3/2}$ 取代单位段上的积分权并处理两个端点，总误差满足

$$
 \left|\mathcal B(Y)-\int_Y^{Y^{3/2}}[H(x)]_+\,\frac{dx}{x^{3/2}}\right|
                 \ll \frac{\log(2Y)}{\sqrt Y}. \tag{4}
$$

事实上半步的修正由 $\sum_{n\ge Y-1}\log(2n)n^{-3/2}$ 控制，权的变化每段为 $O(n^{-5/2})$，端点的不完整单位段各为 $O(Y^{-1/2})$。

**证明 (1) 的上界。** 对任意充分小的 $\varepsilon>0$，在 $\Re s\ge\beta_*+\varepsilon/3$、$\Re s\le2$，对数导数 $-L'/L$ 由单位高度零点数和离零点的距离界为 $O_\varepsilon(\log(2+|\Im s|))$。在该条带且 $|t|\ge2$，有理 RH 的 Lindelöf 推论给出 $|\zeta(\sigma+it)|\ll_{\varepsilon,\nu}(1+|t|)^\nu$，其中 $\nu>0$ 可任意小。对系数 $a(m)$ 在半整数截断点 $u$ 用高度 $u^2$ 的 Perron 公式，$|a(m)|\le4\log m$ 控制截断误差。移线时只跨过 1 处的极点，取 $\nu<\varepsilon/6$，并吸收对数，得到

$$
 D(x)=-\frac{L'}{L}(1)x+O_\varepsilon(x^{\beta_*+\varepsilon}).
$$

阶梯函数使半整数估计延至所有 $x$。从而 $H(x)\ll_\varepsilon x^{\beta_*+\varepsilon}$。对原 $A_Y$，同样的估计由 $\zeta_M=\zeta_B L^2$、阿贝尔域的 ERH 及原素数转移式得到。将两条界积分到 $Y^{3/2}$，再令 $\varepsilon\downarrow0$，即得两个上极限。[Titchmarsh 关于 RH 下的临界带增长，第 14 章](https://sites.math.rutgers.edu/~zeilberg/EM18/TitchmarshZeta.pdf) · [固定域素数计数与 $A_Y$ 的转移](winpair-inert-transfer.md)

**下界用窗口的末段。** 当 $\beta_*=1/2$，$\log(1+P)$ 及 $\log(1+\mathcal B)$ 非负已足够。否则固定 $1/2<\theta<\beta_*$ 与 $1<a<3/2$。每个充分大的 $Y$ 都有

$$
 \int_{Y^a}^{Y^{3/2}} [H(x)]_+\,\frac{dx}{x^{1+\theta}}>1.
$$

若有无界的反例序列，以 $V=a\log Y$ 为左端，窗口右端为 $3V/(2a)$，正部积分不超过 1 就满足原窗口引理的全部次指数预算。由 (3) 将排除 $\Re\rho>\theta$ 的零点，与 $\beta_*>\theta$ 矛盾。同一论证适用于 $T(x)=\sum_{m\le x}b(m)$，其变换为 $-L'/\bigl(sL\bigr)$。已有 $T(x)=24(A_Y(x)-x/48)+O(\sqrt x\log^2(2x))$，在此末段按 $dx/x^{1+\theta}$ 积分的误差趋于零。故用 $A_Y(x)-x/48$ 取代 $H(x)$ 时，同样得到固定正下界。再用 $x^{\theta-1/2}\ge Y^{a(\theta-1/2)}$，结合 (4)，并令 $a\uparrow3/2$、$\theta\uparrow\beta_*$，即得两个下极限。[Levent–Claude 的抽象窗口引理](https://github.com/leanprover-community/riemann-hypothesis/blob/9c9f665ed79ce0fa8da61bd803b7a4bf1ad62ed0/src/Development/AutomorphicArithmetic/DedekindLFunctions/LandauTheorems/landauWindow.lean#L125)

## 用一个 Hermitian 格记录整个窗口的正部

取 $\mathcal L_n=\operatorname{Nm}_{M/B}(\mathfrak b_M(n))\mathfrak b_B(n)^{-1}$，它是 $\mathcal O_B$ 上的可逆分式理想，$\log\mathrm N_B\mathcal L_n=2h(n)$。所有分式理想用其在 $B$ 的各无穷嵌入下的通常绝对值作度量。采用未除以 $[B:\mathbb Q]$ 的 Arakelov 次数，则 $\widehat{\deg}(I)=-\log\mathrm N_B I$。设 $X=Y^{3/2}$，并置

$$
 W=\lceil Y^3\rceil,\quad k_n=\lfloor W/n^{3/2}\rfloor,\qquad
 \mathcal W_Y=\bigoplus_{Y\le n<X}\mathcal L_n^{-k_n}.
$$

无穷处取正交直和。记 $\widehat{\deg}_+$ 为所有饱和子格的最大次数，允许零子格。对上述直和，它等于各条线的正次数之和。确实，任一秩 $r$ 子格的行列式至少有一个非零的 $r$ 阶坐标投影，投影在有限处是整映射、在无穷处是压缩。由线丛次数不等式，该子格的次数不超过最大 $r$ 条线的次数之和。取所有正次数直和即得等号。这给出

$$
 0\le\mathcal B(Y)-\frac{\widehat{\deg}_+(\mathcal W_Y)}{2W}
       \le \frac1W\sum_{Y\le n<X}[h(n)]_+=O(1). \tag{5}
$$

末步用 $h(n)=O(n)$、$X^2\le W$。因此在多项式秩 $O(Y^{3/2})$ 的这个显式格上，只需对每个 $\eta>0$ 在任意远窗口取得 $\widehat{\deg}_+\le W Y^\eta$。反向也成立，固定因子可用较小的 $\eta$ 吸收。仅总次数的估计会丢掉正部。例如在 $\mathbb Z$ 上，$m\mathbb Z\oplus m^{-1}\mathbb Z$ 的总次数恒为 0，正次数却等于 $\log m$。

**实际构造一个整映射。** 对所有素理想 $\mathfrak l$，局部地板差均非负。若 $\mathfrak l$ 在 $M/B$ 中的剩余次数为 $f$，上方各素处理想范数的幂次只从 $(\mathrm N_B\mathfrak l)^{fj}$ 取样，且 $\sum_{\mathfrak P\mid\mathfrak l} f_{\mathfrak P/\mathfrak l}\le3$。故 $\mathfrak b_B(n)^3\subseteq\operatorname{Nm}_{M/B}\mathfrak b_M(n)$。同理 $\binom{2n}{n}\mathcal O_B\subseteq\mathfrak b_B(n)$。于是 $\binom{2n}{n}^2\in\mathcal L_n$，逐行乘以 $\binom{2n}{n}^{2k_n}$ 就给出 $\mathcal W_Y$ 到同秩标准格的整单射。

一般地，对一个这样的整单射 $j_Y$，设 $s_l(j_{Y,v})$ 是无穷处的奇异值，$d_v=[B_v:\mathbb R]$。将外幂映射作用到一个取得正次数的子格的行列式，利用标准格的子格次数不大于零，得到

$$
 \widehat{\deg}_+(\mathcal W_Y)
 \le \sum_{v\mid\infty}d_v\sum_l[\log s_l(j_{Y,v})]_+. \tag{6}
$$

上面的二项式单射给出的右端大小为 $\asymp W Y^{3/4}$，因为 $\log\binom{2n}{n}\asymp n$ 且 $\sum_{Y\le n<X}k_n n\asymp W\sqrt X$。这一步量出了直接用有理整值多项式比较时的损失。Minkowski 定理也可对各个分式理想 $\mathcal L_n^{k_n}$ 分别取非零元素，使每个无穷位的绝对值至多 $C_B(\mathrm N_B\mathcal L_n^{k_n})^{1/[B:\mathbb Q]}$。用这些元素作对角线，(6) 的右端至多为 $\widehat{\deg}_+(\mathcal W_Y)+O_B(X)$，但此估计的右端仍含待控制的理想范数。要由 (6) 证明所求积分界，尚须用区分 $M/B$ 三阶分裂信息的几何估计，把奇异值代价压到 $W Y^{o(1)}$。[原几何证明中的行列式与局部增益](https://hexagonmath.org/pdf/2611.00812) · [用混合 CM 型实际提取所需的分裂信息](cm229-local-length-window.md)

后续工作从一个混合 CM 型求出了 $\mathcal L_n$ 的正有理生成元，使 (6) 在每个窗口严格取等。只读取平方根以上的素数，额外代价为 $O(W\log Y)$。在本例的 ERH 下，最小代价与 $W\log Y$ 的比值有正的极限。[最优整同构](cm229-optimal-splitting-map.md) · [临界线上的对数渐近式](cm229-critical-logarithm.md)

## 三个完整窗口的计算

下表的第三列是**固定域**的 $\int_Y^{Y^{3/2}}[\Theta_S(x)/4-x/48]_+x^{-3/2}\,dx$，$S$ 为 $M$ 的完全分裂素数。计算移动宿主平均 $P(Y)$ 另需相应的有限抽样误差。第二列是 (2) 的离散积分。

| 窗口 $[Y,Y^{3/2}]$ | $\mathcal B(Y)$ 的严格区间 | 固定域正部积分 |
| --- | --- | --- |
| $[1024,32768]$ | $(1.955656,1.955657)$ | 精确为 0 |
| $[4096,262144]$ | $(5.191168,5.191169)$ | 精确为 0 |
| $[10000,1000000]$ | $(8.864300,8.864301)$ | $(0.00855272,0.00855273)$ |

程序逐素数幂比较两种系数算法。一种使用 $M,B$ 的剩余次数及分歧指数，另一种使用三次多项式的独立因式分解和两个 Artin 表示的迹。分歧处 $p=2$ 在两域的 $(e,f)$ 均为 $(2,2)$，$p=229$ 均为 $(2,1)$。再用除数卷积累加 $D$，并在每个窗口端点直接用局部阶乘赋值复核。固定域积分在每两个相邻完全分裂素数之间对正部作精确分段，再用 160 位 Arb 球算术向外舍入。前两行的零通过逐段严格比较得出。将两份程序放在同一目录即可复算。[本篇核验程序](verify_cm229_binomial.py) · [共享素数筛与 Frobenius 判据](verify_cm229_split_primes.py)

~~~sh
python -m pip install python-flint==0.9.0
python -B verify_cm229_binomial.py
~~~
