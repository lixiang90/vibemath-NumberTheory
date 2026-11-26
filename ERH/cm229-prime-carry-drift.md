---
title: "素数平方的偏移与三阶 Hecke 进位和的累计正部"
date: "2026-11-27"
abstract: "将平方根以上素数的带符号进位和与整值多项式格的对数生成元作比较，算出被删去的平方根主项和完整窗口上幂次衰减的误差。在本域 ERH 下，原始素数正部有独立的对数渐近式。利用三个已隔离的 Hecke 零点及其整数线性组合，得到对所有充分远窗口成立的显式对数下界。"
---

# 素数平方的偏移与三阶 Hecke 进位和的累计正部

沿用 $F=\mathbb Q(\alpha)$、$\alpha^3-4\alpha+1=0$、$E=F(i)$、$M=F^{\rm gal}(i)$、$B=\mathbb Q(\sqrt{229},i)$。取 $M/B$ 的非平凡三阶 Hecke 特征 $\psi$，令 $L=L_B(s,\psi)=L_+L_-$。对在 $B$ 完全分裂的未分歧有理素数，$w(p)$ 按相对分裂或惰性取 $2,-1$，其余取 0。写 $\epsilon_q(n)=\lfloor2n/q\rfloor-2\lfloor n/q\rfloor$，

$$
 D(n)=\sum_{\sqrt{2n}<p\le2n}w(p)\epsilon_p(n)\log p,\qquad
 P(Y)=\frac1W\sum_{Y\le n<Y^{3/2}}k_n[D(n)]_+ ,
$$

其中 $W=\lceil Y^3\rceil$，$k_n=\lfloor W/n^{3/2}\rfloor$。此处 $P$ 指大素数进位和的正部。用 $\log r_n$ 表示伴篇最优整同构第 $n$ 列在取幂前的对数生成元。[有理生成元的局部公式](cm229-optimal-splitting-map.md) · [三阶特征与进位区间](cm229-cubic-hecke-carries.md)

**定理 1。** 令

$$
 \kappa=(\sqrt2-2)\zeta(1/2)=0.85545586538795644\ldots .
$$

对所有充分大的实数 $Y$ 都有

$$
 \sum_{Y\le n<Y^{3/2}}\frac{|D(n)-\log r_n+\kappa\sqrt n|}{n^{3/2}}
                      \ll Y^{-1/38}. \tag{1}
$$

前篇给出的逐点界为 $D(n)-\log r_n=O(\sqrt n)$。式 (1) 确定了这个差在整窗内的主项是 $-\kappa\sqrt n$，其余误差在取绝对值以后仍然可加。本式使用已公布的 $7/8$ 无零半平面以及有理 Dirichlet GRH。[Hu 的数域 Hecke 无零定理](https://hexagonmath.org/2610.00216v2) · [有理 Dirichlet GRH 的核验](https://math.mit.edu/~drew/grh/)

## 素数平方的平凡特征分量

对于 $p\nmid458$，写 $g\in S_3$ 为 $F^{\rm gal}$ 中的 Frobenius，$\rho$ 为二维标准表示。$-L'/L=\sum b(m)m^{-s}$ 的平方项满足

$$
 a_2(p):=\frac{b(p^2)}{2\log p}=\operatorname{tr}\rho(g^2)
      =1-\operatorname{sgn}(g)+\operatorname{tr}\rho(g). \tag{2}
$$

$L_-$ 的二次扭转在平方上取值为 1。右端在单位元、换位、三循环上依次取 $2,2,-1$，平凡表示的系数为 1。按有理范数收集欧拉因子时，$p^2$ 也包含 $B$ 中剩余次数为 2 的素理想的一次贡献。由有理素数定理、二次 Dirichlet 特征估计，以及 $L_+$ 的无零半平面，对 $b_0=8/9>7/8$ 有

$$
 \Theta_2(t)=\sum_{p\le t,\,p\nmid458}a_2(p)\log p
                  =t+O(t^{b_0}). \tag{3}
$$

具体地，对三个表示的对数导数在 $7/8$ 与 $8/9$ 间移线，平凡表示给出极点主项，删去素数高次幂只花 $O(\sqrt t)$。Lean 按六种置换检验了 (2)，区间程序另对每个用到的素数幂比较从 $B$ 的三阶特征和从两个有理 Artin 因子算出的迹。[Lean 群作用恒等式](CM229Frobenius.lean) · [本例的欧拉因子](cm229-heights-and-zeros.md)

平方造成进位恰当 $p$ 落在某个 $(\sqrt{2n/(2j+2)},\sqrt{2n/(2j+1)}]$。用 (3) 对前 $J$ 段求和，误差为 $O(n^{b_0/2}J^{1-b_0/2})$。其余素数都小于等于 $\sqrt{2n/(2J+1)}$，绝对权之和由 Chebyshev 界控制为 $O(\sqrt{n/J})$。主项中的交错级数尾也有同一上界。取 $J=\lfloor n^{1/19}\rfloor$，得到

$$
 \sum_{p\nmid458}a_2(p)\epsilon_{p^2}(n)\log p
 =\sqrt{2n}\sum_{j\ge0}\bigl((2j+1)^{-1/2}-(2j+2)^{-1/2}\bigr)
       +O(n^{9/19})
 =\kappa\sqrt n+O(n^{9/19}). \tag{4}
$$

在 $\log r_n=\tfrac12\sum_{p^j\le2n}b(p^j)\epsilon_{p^j}(n)$ 中，三次及以上幂的绝对贡献为 $O(n^{1/3}\log(2n))$，两个分歧素数总计 $O(\log(2n))$。还需扣去 $p\le\sqrt{2n}$ 的一次项。其模 $p$ 均值加总为 $O(n^{4/9})$，由 $\Theta(t)=\sum_{p\le t}w(p)\log p\ll t^{8/9}$ 即得。其中心化部分在 $T\le n<2T$ 按 $n^{-3/2}$ 计的绝对值和，由伴篇加法大筛 (5) 取 $P=2\sqrt T$ 得 $O(T^{-1/4}\log^{3/2}(2T))$。式 (4) 的逐点余项在同一块贡献 $O(T^{-1/38})$。合并并从 $Y$ 起求和，即证 (1)。[移动终点的加法大筛](cm229-cubic-hecke-carries.md)

## 原始大素数正部的临界极限

设 $A(x)=\sum_{m\le x}\sum_{d\mid m}b(d)$，$H(x)=A(2x)-2A(x)$，$G(u)=e^{-u/2}H(e^u)$。有 $H(n)=2\log r_n$。对 $U=\log Y$，(1) 推出

$$
 P(Y)=\int_U^{3U/2}[G(u)/2-\kappa]_+\,du+O(Y^{-1/38}). \tag{5}
$$

为检验误差，先以 $k_n/W\le n^{-3/2}$ 和正部的 1-Lipschitz 性使用 (1)。又 $H(n)\ll n^{8/9}$，因为前篇的 Perron 估计允许任意指数 $b>\beta_*$，而 $\beta_*\le7/8$。用 $n^{-3/2}$ 替换 $k_n/W$ 只花 $O(W^{-1}(Y^{3/2})^{1+8/9})=O(Y^{-1/6})$。在 $[n,n+1)$ 的后半段，$H$ 比前半段多 $\sum_{d\mid2n+1}b(d)=O(\log(2n))$，$\sqrt x$ 在单位段内只变动 $O(n^{-1/2})$。因此替换积分、权和两个端点的总误差为 $O(Y^{-1/2}\log(2Y))$。最后令 $x=e^u$ 即得 (5)。[阶乘格的 Mellin 变换与 Perron 估计](cm229-binomial-slopes.md)

**定理 2。** 若 $\operatorname{ERH}(E)$ 成立，则存在 $d_{229}>0.01606$ 使

$$
 \frac{P(Y)}{\log Y}\longrightarrow d_{229}\qquad(Y\to\infty). \tag{6}
$$

若 $E$ 有临界线外的非平凡零点，此比值趋于 $+\infty$。因此对每个充分远的完整窗口，无条件有 $\sum_n k_n[D(n)]_+>0.01606\,W\log Y$。

为证有限极限，假设本域 ERH，写 $\mathcal M(f)=\lim_{U\to\infty}(2/U)\int_U^{3U/2}f(u)\,du$，在极限存在时使用。前篇的均方逼近给出 $G$ 在这些区间上由有限和 $\sum_{|\gamma|<T}A_\gamma e^{i\gamma u}$ 逼近，系数是

$$
 A_\gamma=-m_\rho(2^\rho-2)\zeta(\rho)/\rho,\qquad
 \rho=1/2+i\gamma,\quad L(\rho)=0 .
$$

频率重合时使用总重数。由同一均方界及正部的 Lipschitz 性，(5) 给出 $d_{229}=\tfrac12\mathcal M([G/2-\kappa]_+)$。负偏移应放在取正部之前，故这里需要重新证明极限严格为正。若 ERH 不成立，伴篇的精确指数 $\tfrac32(\beta_*-1/2)>0$ 已给出第二种情形。[零点留数的均方逼近](cm229-critical-logarithm.md) · [带符号大素数正部的增长指数](cm229-cubic-hecke-carries.md)

## 三个零点足以越过负偏移

取高度最小的三个正零点，记其纵坐标为 $\gamma_1,\gamma_2,\gamma_3$。它们的隔离区间都为 $(a,a+1)/2^{24}$，分子 $a$ 依次为 $13347347,36732456,43857713$。三者在 $L$ 中均为单零点。设 $\varphi_j=\arg A_{\gamma_j}$，$v_j=\gamma_j u+\varphi_j$，用非负函数

$$
 K(u)=\tfrac23(1+\cos v_1)^2(1+\cos v_2)(1+\cos v_3),
 \qquad 0\le K\le32/3 .
$$

展开中的全部频率为 $\sum_j e_j\gamma_j$，其中 $|e_1|\le2$，$|e_2|,|e_3|\le1$。对 $\sum |e_j|\ge2$ 用整数区间运算取绝对值。所得区间距 0 至少 $6222089/2^{24}$，距全部 192 个已认证的正零点区间至少 $487800/2^{24}$，且上端点至多 $107284867/2^{24}<64$。所以这些频率在 $G$ 的均值傅里叶系数均为零。只用有限次整数比较即可检验这些关系。[192 个隔离区间](cm229-zeros64.csv) · [总零点数与单性证书](verify_cm229.py)

常数项给 $\mathcal M(K)=1$。三个一次余弦在 $K$ 中的系数为 $4/3,1,1$，而 $\mathcal M(G\cos v_j)=|A_{\gamma_j}|$。对 $J=G/2-\kappa$ 逐点有 $(32/3)[J]_+\ge JK$。因此

$$
 d_{229}\ge
 \frac3{64}\left(\frac43\frac{|A_{\gamma_1}|}{2}
        +\frac{|A_{\gamma_2}|}{2}+\frac{|A_{\gamma_3}|}{2}-\kappa\right)
 >0.01606 . \tag{7}
$$

按完整隔离区间在复球中计算留数，(7) 的中间量严格落在 $(0.0160647,0.0160652)$ 内。这是一个下界的证书，$d_{229}$ 本身由完整的零点谱决定。[区间算术与频率间隔核验](verify_cm229_hecke_carries.py)

## 三个完整窗口

各列都除以 $W$，每个表值都是严格开区间的下端点，上端点为表值加 $10^{-6}$。第二列是题设的实际累计正部，第三列先对 $\sqrt{2n}<p\le(2n)^{4/7}$ 求带符号和再取绝对值，第四列为 $p>(2n)^{2/3}$ 的累计正部。最后一列为 $\sum_n k_n|D(n)-\log r_n+\kappa\sqrt n|/W$。

| $[Y,Y^{3/2}]$ | 原始正部 | 中间素数段绝对值 | 三分之二次幂以上正部 | 扣除平方主项后的绝对误差 |
| --- | ---: | ---: | ---: | ---: |
| $[1024,32768]$ | 0.320451 | 0.304625 | 0.329555 | 0.692335 |
| $[4096,262144]$ | 1.074672 | 0.334594 | 1.203290 | 0.657061 |
| $[10000,1000000]$ | 2.209373 | 0.315753 | 2.493294 | 0.621041 |

程序以整数幂比较确定移动截断，遍历全部次数和进位跳点，用整数平方根计算 $k_n$。在每个窗口的首末整数上另按局部赋值和分裂类型直接计算。复算需要 python-flint 0.9.0，将程序、CSV 与下列三个被引用的程序放在同一目录即可。[本篇核验程序](verify_cm229_hecke_carries.py) · [进位跳点](verify_cm229_optimal_map.py) · [独立 Artin 迹](verify_cm229_binomial.py) · [素数筛](verify_cm229_split_primes.py) · [零点 CSV](cm229-zeros64.csv)

~~~sh
python -B verify_cm229_hecke_carries.py
~~~
