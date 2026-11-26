---
title: "三阶 Hecke 特征在平方根以上进位素数中的抵消区间"
date: "2026-11-27"
abstract: "在判别式 229 的三次域所给的数域塔中，将相对分裂与惰性的权重写为一个 Hecke 特征的迹。用加法大筛证明，从平方根到 4/7 次幂的素数段在完整窗口上的绝对值累计量有幂次节省。靠近两倍次数的一小段素数仍能检出原和的零点指数。用次数二 L 函数的系数均方估计剔除小因子后，所需抵消归结为单个双线性和。"
---

# 三阶 Hecke 特征在平方根以上进位素数中的抵消区间

令 $F=\mathbb Q(\alpha)$，$\alpha^3-4\alpha+1=0$，$N=F^{\rm gal}$，$k=\mathbb Q(\sqrt{229})$，$B=k(i)$，$M=N(i)$，$E=F(i)$。$M/B$ 是处处无分歧的循环三次扩张。记 $S_0,S$ 为分别在 $B,M$ 完全分裂的未分歧有理素数。对全部有理素数置 $w(p)=3\mathbf1_S(p)-\mathbf1_{S_0}(p)$，并约定 $\mathbf1_S,\mathbf1_{S_0}$ 在 $2,229$ 处为零。在 $S_0$ 中，$w$ 对相对分裂取 2，对相对惰性取 $-1$。

令 $\epsilon_q(n)=\lfloor 2n/q\rfloor-2\lfloor n/q\rfloor$。取 $Y\ge2$，$X=Y^{3/2}$，$W=\lceil Y^3\rceil$，$k_n=\lfloor W/n^{3/2}\rfloor$。以下关于 $n$ 的窗口和均遍历 $Y\le n<X$ 的所有整数。定义

$$
 D_\lambda(n)=\sum_{(2n)^\lambda<p\le 2n}w(p)\epsilon_p(n)\log p,
 \qquad P_\lambda(Y)=\frac1W\sum_n k_n[D_\lambda(n)]_+ .
 \tag{1}
$$

所求的是每个 $\eta>0$ 在任意远的完整窗口上有 $P_{1/2}(Y)\le Y^\eta$。本篇证明在所有充分大的窗口上，

$$
 \frac1W\sum_n k_n\,|D_{1/2}(n)-D_{4/7}(n)|\ll Y^{-9/10000}. \tag{2}
$$

式 (2) 采用 Liu 预印本的无零半平面。只采用 Hu 的 $7/8$ 定理，右端也已对每个 $\eta>0$ 为 $O_\eta(Y^\eta)$。超过 $(2n)^{4/7}$ 的部分仍有同一个 ERH 障碍。下文进一步将它定位在靠近 $2n$ 的单个素数区间。[全局格映射与原始进位和](cm229-optimal-splitting-map.md) · [Hu 的八分之七定理](https://hexagonmath.org/2610.00216v2) · [Liu 定理 7.33](https://hexagonmath.org/pdf/2611.00741v1)

## 理想类特征如何给出这两个权重

取 $H=\operatorname{Gal}(M/B)\simeq C_3$ 的一个非平凡特征 $\psi$，亦记作 $B$ 的非分歧 Hecke 特征。可从 $N/k$ 的三阶理想类特征经 $N_{B/k}$ 拉回。对 $p\in S_0$ 及任一 $\mathfrak p\mid p$，有

$$
 w(p)=\psi(\mathfrak p)+\bar\psi(\mathfrak p).
$$

设 $\rho$ 为 $S_3$ 的二维标准表示，$\chi_{-4}$ 为高斯二次特征。由于 $\operatorname{Gal}(M/\mathbb Q)=S_3\times C_2$，从 $A_3\times\{1\}$ 诱导 $\psi$，得到 $\rho\oplus(\rho\otimes\chi_{-4})$。故

$$
 L_B(s,\psi)=L_B(s,\bar\psi)=L_+(s)L_-(s)=:L(s).
 \tag{3}
$$

换位将 $\psi$ 送到 $\bar\psi$，并置换同一有理素数上方的素理想。在相对惰性的 $p\in S_0$ 上，四个素理想中两个取 $\exp(2\pi i/3)$，另两个取其共轭。因此 $-L'/L=\sum_m b(m)m^{-s}$ 的一次素数系数为 $b(p)=2w(p)\log p$（$p\nmid458$）。对任意只依赖范数的有限权重，两个共轭特征的素理想和逐项按范数相等，故两者的平方模平均就等于其中一个和的平方模。[类域互反与理想范数](https://www.jmilne.org/math/CourseNotes/CFT.pdf) · [本例两个 Artin 因子的分解](cm229-heights-and-zeros.md)

写 $\beta_*=\sup_{L(\rho')=0,\,0<\Re\rho'<1}\Re\rho'$，以及 $\Theta(x)=\sum_{p\le x}w(p)\log p$。对 $\beta_*<b<1$ 有

$$
 \Theta(x)\ll_b x^b. \tag{4}
$$

在 $\beta_*$ 与 $b$ 之间取直线对 $-L'/L$ 使用截断 Perron 公式，单位高度零点数和至零点的距离给出该线上对数导数的对数增长。再删去高次素数幂，代价为 $O(\sqrt x)$，两个分歧素数贡献 $O(\log(2x))$。这也解释了 (3) 中两个共轭特征所共享的零点为何直接影响此素数和。[相同系数的移线估计](cm229-binomial-slopes.md)

## 移动进位的加法大筛

对奇素数 $p$，$\epsilon_p(n)=\mathbf1_{\{n\bmod p\ge(p+1)/2\}}$，其模 $p$ 平均值为 $\mu_p=(p-1)/(2p)$。令

$$
 R(n,V)=\sum_{2<p\le V}w(p)\log p\,(\epsilon_p(n)-\mu_p).
$$

**引理。** 对 $2\le P\le 4T$，有

$$
 \sum_{T\le n<2T}\max_{V\le P}|R(n,V)|^2
 \ll (T+P^2)P\log^3(2P). \tag{5}
$$

将中心化进位展开在频率 $h/p$，$1\le h<p$。由有限 Fourier 的 Parseval 等式，每个 $p$ 的归一化系数模方和为 $\mu_p(1-\mu_p)\le1/4$。这些既约分数在圆周上的间距至少 $P^{-2}$。对任一素数子集应用 Montgomery–Vaughan 的对偶加法大筛，所得平方和至多为常数乘 $(T+P^2)\sum|w(p)\log p|^2$，求和限于该子集。按大小排列素数，并按下标作二进分块。任一前缀每层至多用一块，Cauchy–Schwarz 花费 $O(\log(2P))$。将各块的平方和相加，同一层内互不交，再花 $O(\log(2P))$。最后用 $\sum_{p\le P}\log^2p\ll P\log(2P)$，得到 (5)。所以移动终点可以依赖 $n$。[Montgomery–Vaughan，加法大筛及其对偶](https://personal.science.psu.edu/rcv4/personal/Publications/large_sieve.pdf)

取 $1/2\le\lambda\le2/3$，在一个二进块上用 $P=(4T)^\lambda$。由 (4) 和 $\sum_{p\le V}\log p/p\ll\log(2V)$，均值项是 $\Theta(V)/2+O(\log(2V))$。用 (5) 和 Cauchy–Schwarz 处理两个端点，得到

$$
 \sum_{T\le n<2T}\frac{|D_{1/2}(n)-D_\lambda(n)|}{n^{3/2}}
 \ll_b T^{\lambda b-1/2}
       +\bigl(T^{(\lambda-1)/2}+T^{3\lambda/2-1}\bigr)\log^{3/2}(2T).
 \tag{6}
$$

**(2) 的证明。** 由 Liu 的界可取 $b=0.87334$。当 $\lambda=4/7$，$\tfrac12-4b/7>0.00094$，而 $3\lambda/2-1=-1/7$。从 $Y$ 起对二进块求和，并用 $k_n/W\le n^{-3/2}$，即得 (2)。若取 $b=7/8+\varepsilon$，可由同式得到开头引用 Hu 定理的结论。一般地，任何固定 $\lambda<\min(2/3,1/(2b))$ 都有归一化误差的幂次衰减。

到达 $\lambda=2/3$ 时，(5) 仍将**中心化进位**在整窗的归一化绝对值和控制在 $O(\log^{5/2}(2Y))$。进入 (6) 第一项的算术量只剩两个普通素数和 $\Theta((2n)^{2/3})$ 与 $\Theta(\sqrt{2n})$。如果进一步要求 $\lambda=2/3$ 的整段未经中心化就有幂次节省，则 (6) 中的均值需要额外控制。

稍微降低终点，就能用一个普通 Hecke 素数和替换几乎整段进位。具体地，取 $V(n)=(2n)^{2/3}/\log^2(2n)$，令 $D_{>V}(n)$ 表示 (1) 在 $V(n)<p\le2n$ 的和。则

$$
 \frac1W\sum_n k_n\left|D_{1/2}(n)-D_{>V}(n)-\tfrac12\Theta(V(n))\right|
                         \ll (\log Y)^{-1/2}. \tag{6a}
$$

事实上每块用 $P\ll T^{2/3}/\log^2 T$，(5) 给出的归一化误差为 $O(T^{-1/6}\log^{1/2}(2T)+\log^{-3/2}(2T))$。另一端的均值由 $\Theta(\sqrt{2n})\ll n^{4/9}$ 控制。整窗有 $O(\log Y)$ 块，即得 (6a)。

## 哪些大素数仍然检出全部零点指数

记 $\delta=\tfrac32(\beta_*-\tfrac12)$。对每个 $1/2\le\lambda\le2/3$，

$$
 \lim_{Y\to\infty}\frac{\log(1+P_\lambda(Y))}{\log Y}=\delta . \tag{7}
$$

当 $\lambda=1/2$ 时，前篇得到 $P_{1/2}(Y)=C_*(Y)/(4W)+O(\log(2Y))$，其中 $C_*$ 已有精确指数 $\delta$。对一般的 $\lambda$，由正部的 1-Lipschitz 性，误差受 (6) 的整窗和控制。以 $b\downarrow\beta_*$，该误差加 1 的上指数至多为 $\max(0,\tfrac32(\lambda\beta_*-1/2))$。若 $\beta_*>1/2$，它严格小于 $\delta$。若 $\beta_*=1/2$，两者都有任意小正指数的上界，而 $\log(1+P_\lambda)\ge0$。此即 (7)。[最优格映射的完整窗口指数](cm229-optimal-splitting-map.md)

还有一个只需检验**单个素数区间**的判据。每个窗口选定一个正整数 $m=m_Y$，满足 $\log m=o(\log Y)$，令 $r=2^{1/m}$，并置

$$
 I_Y(x)=\Theta(2x)-\Theta(2x/r),\qquad
 P_{\rm end}(Y)=W^{-1}\sum_n k_n[I_Y(n)]_+ .
$$

在 $(2n/r,2n]$ 内每个素数的进位为 1。即使 $m_Y$ 趋于无穷，例如取 $\lceil\exp(\sqrt{\log Y})\rceil$，仍然有

$$
 \lim_{Y\to\infty}\frac{\log(1+P_{\rm end}(Y))}{\log Y}=\delta. \tag{8}
$$

**证明。** $(4)$ 直接给出上界。设 $\beta_*>1/2$，固定 $1/2<\theta<\beta_*$ 及 $1<a<3/2$。令 $T_L(x)=\sum_{m'\le x}b(m')$。$T_L(2x)-T_L(x)$ 的 Mellin 变换为 $(2^s-1)(-L'/L(s))/s$。在 $\Re s>0$，因子 $2^s-1$ 不消去任何零点带来的极点，且低零点证书已排除实非平凡零点。由 Alpöge–Claude 窗口引理，每个充分大的 $Y$ 都有 $\int_{Y^a}^{X}[T_L(2x)-T_L(x)]_+x^{-1-\theta}\,dx>1$，否则反例序列就满足该引理的任意次幂预算。将被积差换为 $2(\Theta(2x)-\Theta(x))$ 的误差是 $O(\sqrt{x})$，在该权下费用趋于零。故对权 $x^{-3/2}dx$，这个末段上的正部积分至少为常数乘 $Y^{a(\theta-1/2)}$。[Levent–Claude 的完整窗口引理](https://github.com/leanprover-community/riemann-hypothesis/blob/9c9f665ed79ce0fa8da61bd803b7a4bf1ad62ed0/src/Development/AutomorphicArithmetic/DedekindLFunctions/LandauTheorems/landauWindow.lean#L125) · [本域的实零点证书](cm229-heights-and-zeros.md)

精确分割给出 $\Theta(2x)-\Theta(x)=\sum_{j=0}^{m-1} I_Y(x/r^j)$。当 $x\in[Y^a,X]$ 时，各 $x/r^j$ 均在原窗口内。取正部，用次可加性并换元，每项权系数为 $r^{-j/2}\le1$，故上述末段积分至多为 $m\int_Y^X[I_Y(t)]_+t^{-3/2}dt$。在每个单位区间内，$I_Y$ 的两个端点各自移动至多 2，取整误差为 $O(\log(2t))$，与 $m$ 无关。由 $|I_Y(t)|\ll t$，处理 $k_n$ 的取整及端点后，连续积分与 $P_{\rm end}(Y)$ 相差 $O(1)$。用 $m=Y^{o(1)}$，再令 $a\uparrow3/2$、$\theta\uparrow\beta_*$ 即得下界。

所以欲得原目标，只需对任一上述选择的窄区间，证明三阶特征的相对分裂盈余在任意远整窗上有任意小正指数的正部上界。式 (8) 表明，这一步保留 $L_B(s,\psi)$ 在临界线右侧的所有零点障碍。目前对完整素数范围可代入 (7) 的外部无零界给出 $P_{1/2}(Y)\ll Y^{0.56}$。下面将这个窄区间再分解一次。伴篇另给出原和的对数级下界及三个完整窗口的数据。[平方项与实际累计正部](cm229-prime-carry-drift.md)

## 将剩余估计写成单个双线性和

固定 $m_Y=\lceil\log^2(2Y)\rceil$，$r=2^{1/m_Y}$，在每个二进块 $T\le n<2T$ 上取 $z=\max(1,\lfloor T^{1/6}/\log^3(2T)\rfloor)$，其中 $T=2^jY$。令 $a(h)$、$u(h)$、$b_+(h)$ 分别为 $L_+$、$1/L_+$、$-L_+' / L_+$ 在 $\Re s>1$ 的 Dirichlet 系数。$u$ 是 $a$ 的卷积逆。定义

$$
 c_z(h)=\mathbf1_{\{h=1\}}-\sum_{\substack{d\mid h\\d\le z}}u(d)a(h/d),\qquad
 R_Y(n)=\frac12
 \sum_{\substack{h,\ell>z\\2n/r<h\ell\le2n}}
   (1+\chi_{-4}(h\ell))\,c_z(h)b_+(\ell).
 \tag{9}
$$

**定理。** 在所有充分远的完整窗口上有

$$
 \frac1W\sum_n k_n\,|I_Y(n)-R_Y(n)|\ll\frac1{\log(2Y)}. \tag{10}
$$

所以只须估计 (9) 在任意远完整窗口上的累计正部。两个因子的范数下界随 $T$ 增长，$b_+(\ell)$ 由三阶理想类特征给出，且仅支撑于素数幂。这是一条可独立研究的双线性估计。式 (10) 已控制分解中带小因子的项。

后续工作在保持 $c_z$ 不变的条件下将 $z<\ell\le\sqrt T/z$ 的贡献控制在多对数预算内。只保留 $\ell=p>T^{2/3}$ 的互素素数对仍检出原指数。在本域 ERH 下，全部 $R_Y$ 的归一化绝对值累计量为 $O(\log(2+\log Y))$。[双线性和内部的抵消区间](cm229-type-ii-ranges.md) · [互素大素因子判据](cm229-type-ii-large-prime.md) · [全部因子的有限高度公式](cm229-type-ii-zeros.md)

**证明。** 分别对 $L_\nu$（$\nu=+,-$）定义 $a_\nu,u_\nu,b_\nu$。两个 $L_\nu$ 都是次数为二、具实系数的整 Hecke $L$ 函数，满足标准函数方程，Gamma 位移为实数，且 $|a_\nu(h)|\le d(h)$。因此它们满足 Cao–Tanigawa–Zhai 的次数二均方定理（推论 1，取其参数 $\theta=0$）。设 $A_\nu(t)=\sum_{h\le t}a_\nu(h)$，有

$$
 \int_Z^{2Z}|A_\nu(t)|^2dt\ll Z^{3/2}\qquad(Z\ge1). \tag{11}
$$

该定理中 $L_\nu(s)t^s/s$ 的留数项至多是常数 $L_\nu(0)$，不影响 (11)。两因子均满足 $|u_\nu(h)|\le d(h)$ 及 $|b_\nu(h)|\le2\Lambda(h)$，这是由次数至多二且 Satake 参数绝对值至多 1 的局部因子逐项得到的。[次数二部分和的均方定理，推论 1](https://arxiv.org/pdf/1501.04269) · [本域实系数与 Gamma 因子](cm229-heights-and-zeros.md)

对每个因子暂略 $\nu$，用 $*$ 表示 Dirichlet 卷积，$a\log$ 表示 $a(h)\log h$。因 $a*b=a\log$，按 $z$ 截断后有精确的 Vaughan 恒等式

$$
 b=u_{\le z}*(a\log)-(u_{\le z}*b_{\le z})*a
        +b_{\le z}+c_z*b_{>z},\qquad c_z=u_{>z}*a . \tag{12}
$$

事实上，拆分 $b=u*(a\log)$ 的 $u$，然后将第二部分改为 $(u_{>z}*a)*b$，再拆分 $b$ 即得。$c_z$ 在 $h\le z$ 上为零。窄区间的左端点至少为 $n$，故第三项对大 $T$ 没有贡献。

对 $h\le z^2$，(11) 和 Cauchy–Schwarz 给出

$$
 \frac1T\sum_{T\le n<2T}
 \left|\sum_{2n/(rh)<j\le2n/h} a_\nu(j)(\log j)^e\right|
     \ll (T/h)^{1/4}\log^e(2T)\qquad(e=0,1). \tag{13}
$$

先对连续 $n$ 在两个端点分别用 (11)。$e=1$ 时分部求和，额外的积分由 $\int_1^{4Z}|A_\nu(t)|dt/t\ll Z^{1/4}$ 控制，取 $Z=T/h$。从连续 $n$ 改为整数，每个单位区间内每一端在 $j$ 坐标移动至多 2，所跨系数各为 $O(Z^{1/8}\log^e(2T))$。故 (13) 对所有 $1\le r\le2$ 一致。

第一项外因子支撑于 $h\le z$。第二项外因子支撑于 $h\le z^2$，在 $h\le t$ 上的绝对系数和由 $2\sum_{ij\le t}d(i)\Lambda(j)\ll t\log^2(2t)$ 控制。以 $h^{-1/4}$ 分部求和，结合 (13)，(12) 前两项在一块上的平均绝对值为

$$
 \ll T^{1/4}(z^{3/4}+z^{3/2})\log^2(2T).
$$

乘以 $T^{-1/2}$ 即得该块按 $n^{-3/2}$ 加权的费用。代入 $z$，再对整窗各块求和，总费用为 $O((\log Y)^{-3/2})$。若直接取 $z=\lfloor T^{1/6}\rfloor$，同一证明给出 $O(\log^3(2Y))$，也在所需预算内。

最后比较有理素数与素数幂。$b_-(h)=\chi_{-4}(h)b_+(h)$，同样的扭转关系对 $a,u,c_z$ 均成立，包括偶数处。因此 (12) 的最后两因子卷积对 $+,-$ 两支取和再除以二，正好是 (9)。从 $\tfrac12\sum_{2n/r<h\le2n}(b_+(h)+b_-(h))$ 删除素数幂时，平方项为 $\sqrt{2n}-\sqrt{2n/r}+O(n^{4/9})$。这是因为平方权为 $\operatorname{tr}\rho(g^2)=1-\operatorname{sgn}g+\operatorname{tr}\rho(g)$，其素数定理误差为 $O(t^{8/9})$。其余高次幂的费用为 $O(n^{1/3}\log(2n))$。归一化整窗费用于是为 $O(\log Y/m_Y+Y^{-1/18})$，与前两项的界合在一起即得 (10)。[平方权的推导与素数定理](cm229-prime-carry-drift.md)

如果保留平方项的显式主项，就能对窗口宽度给出一致的幂次误差。对每个窗口任选 $1\le r\le2$，在各块改取 $z=\lfloor T^{4/25}\rfloor$，仍按 (9) 定义 $R_Y$，则

$$
 \frac1W\sum_n k_n
 \left|I_Y(n)+\sqrt{2n}(1-r^{-1/2})-R_Y(n)\right|
          \ll Y^{-1/100}\log^2(2Y). \tag{14}
$$

这里 $I_Y(n)=\Theta(2n)-\Theta(2n/r)$。前两项卷积在每块的归一化费用至多为 $T^{-1/4+3(4/25)/2}\log^2(2T)=T^{-1/100}\log^2(2T)$，而扣除平方主项后的费用是 $O(T^{-1/18})$。沿二进块求和即得 (14)。尤其可直接取 $r=2$ 检验 $(n,2n]$，也可取前述随窗口收缩的区间。

核验程序在 $h\le4096$ 上对三个截断逐一比较 (12) 的每个 $\log p$ 系数，全部使用整数运算，并以两个 Artin 因子的独立局部迹检查合并后的系数。另在 $[4096,262144)$ 的每个整数上，取 $m=82$ 与 (14) 的截断，分别从素数分裂类型计算 $I_Y$、从乘积 $h\ell$ 计算 $R_Y$。归一化后，$I_Y$ 与 $R_Y$ 的累计正部依次为 $0.49468097\ldots$、$0.67649279\ldots$，式 (14) 左端为 $1.66678923\ldots$。每个真值均严格位于所写小数与该数加 $10^{-8}$ 之间。[双线性恒等式与窗口核验](verify_cm229_bilinear.py) · [Hecke 系数](verify_cm229.py) · [其公共依赖](verify_tetra3136.py) · [独立局部迹](verify_cm229_binomial.py) · [素数筛](verify_cm229_split_primes.py)
