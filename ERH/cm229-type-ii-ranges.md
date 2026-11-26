---
title: "三阶 Hecke 双线性和中可以剔除的因子区间"
date: "2026-11-27"
abstract: "对判别式 229 的六次 CM 域，保持原双线性和的第一个系数不变，将第二个因子从六分之一次幂附近一直剔除到三分之一次幂以上。利用两个次数二函数的函数方程和 Dirichlet 多项式均方估计，在完整窗口上取得绝对值累计界。另给出两因子均超过四分之一次幂的等价估计，以及只含唯一大素因子的合数表述。"
---

# 三阶 Hecke 双线性和中可以剔除的因子区间

设 $F=\mathbb Q(\alpha)$，$\alpha^3-4\alpha+1=0$，$E=F(i)$，$k=\mathbb Q(\sqrt{229})$，$B=k(i)$，$M=F^{\rm gal}(i)$。循环三次扩张 $M/B$ 的非平凡 Hecke 特征 $\psi$ 满足 $L_B(s,\psi)=L_+(s)L_-(s)$。这里 $L_+$ 是 $S_3$ 标准表示的次数二函数，$L_-$ 是它的 $\chi_{-4}$ 扭转。对每个因子用 $a_\nu,u_\nu,b_\nu$ 依次表示 $L_\nu,1/L_\nu,-L_\nu'/L_\nu$ 的 Dirichlet 系数。无下标的 $a,u,b_+$ 取自 $L_+$。[本例的函数方程与分歧素数](cm229-heights-and-zeros.md) · [原双线性和式（9）](cm229-cubic-hecke-carries.md)

取 $X=Y^{3/2}$，$W=\lceil Y^3\rceil$，$k_n=\lfloor W/n^{3/2}\rfloor$，$m=\lceil\log^2(2Y)\rceil$，$r=2^{1/m}$。所有窗口包含每个 $Y\le n<X$ 的整数。令 $T=2^jY$，$T\le n<2T$，最后一块在 $X$ 截断。对这个块内的整数截断 $U,V\ge1$ 定义

$$
 \begin{aligned}
 c_U(h)&=\mathbf1_{h=1}-\sum_{\substack{d\mid h\\d\le U}}u(d)a(h/d),\\
 R_{U,V}(n)&=\frac12\sum_{\substack{h>U,\ \ell>V\\2n/r<h\ell\le 2n}}
 (1+\chi_{-4}(h\ell))c_U(h)b_+(\ell).
 \end{aligned} \tag{1}
$$

原式为 $R_Y=R_{z,z}$，其中 $z=\max(1,\lfloor T^{1/6}/\log^3(2T)\rfloor)$。以下写 $\|G\|_Y=W^{-1}\sum_n k_n|G(n)|$，以 $C$ 表示可随公式增大的绝对常数，数域固定。

**定理 1。** 保持 $h$ 侧系数 $c_z$ 不变，取 $V=\lfloor\sqrt T/z\rfloor$。原式中 $z<\ell\le V$ 的部分满足

$$
 \|R_{z,z}-R_{z,V}\|_Y\ll\log^C(2Y). \tag{2}
$$

若只删至 $V_0=\lfloor T^{3/10}\rfloor$，还有幂次节省

$$
 \|R_{z,z}-R_{z,V_0}\|_Y\ll Y^{-1/90}\log^C(2Y). \tag{3}
$$

在 $\ell>V$ 的余和中，所有 $\ell=p^j$、$j\ge3$ 的合计也满足 $O(Y^{-1/18}\log^C(2Y))$ 的归一化绝对值界。因此原式只需控制 $h>z$ 与 $\ell>\sqrt T/z$ 同时成立、且 $\ell$ 为素数或素数平方的部分。界 (2)、(3) 对所有充分大的实数 $Y$ 成立，证明不要求 $E$ 的 ERH。另一篇进一步用欧拉因子和截断倒数，证明只保留 $\ell=p>T^{2/3}$ 的素数仍检出全部正部幂指数。[原系数上的互素大素因子判据](cm229-type-ii-large-prime.md)

## 两个次数二函数的外层卷积估计

固定 $\nu\in\{+,-\}$，写 $a_e(n)=a_\nu(n)(\log n)^e$，其中 $e=0,1$，$a_0=a_\nu$。设 $f$ 支撑于 $d\le D$，$1\le D\le 4T$，$|f(d)|\le C_0 d_2(d)\log(2d)$，$d_j$ 是 $j$ 重除数函数。在任意 $1\le r\le2$ 下有

$$
 \frac1T\sum_{T\le n<2T}
 \left|\sum_{2n/r<q\le2n}(f*a_e)(q)\right|
 \ll_{C_0}\big((TD)^{1/3}+D\big)\log^C(2T). \tag{4}
$$

**证明。** 记 $g=f*a_e$，则 $|g(q)|\ll d_4(q)\log^2(2q)$，从而 $\sum_{q\le 8T}|g(q)|\ll T\log^5(2T)$。先在连续变量 $x\in[T,2T]$ 估计区间和 $S(x)=\sum_{2x/r<q\le 2x}g(q)$。取一个支撑于 $[-1,1]$、积分为 1 的非负光滑函数 $\phi$，用

$$
 S_K(x)=\int_{-1}^1\phi(v)S(xe^{v/K})\,dv,\qquad 2\le K\le T .
$$

对每个 $q$，当它跨过某端点时，相应 $x$ 的集合长度是 $O(T/K)$。因此 $T^{-1}\int_T^{2T}|S-S_K|\,dx\ll(T/K)\log^C(2T)$。同样，从连续平均换回整数时，每个 $q$ 至多影响常数个单位区间，归一化费用 $O(\log^C(2T))$。这些计数均包含上下两个端点。

设 $P(s)=\sum_{d\le D}f(d)d^{-s}$，$\Phi_K(s)=\int\phi(v)e^{sv/K}dv$。$S_K$ 的逆 Mellin 被积函数，去掉 $(2x)^s$，为

$$
 P(s)(-1)^e L_\nu^{(e)}(s)\,\frac{1-r^{-s}}s\,\Phi_K(s).
$$

区间乘子在 $s=0$ 可去，$L_\nu$ 整。将积分线移到 $\sigma=-1/\log(2T)$。固定次数二函数的函数方程、Stirling 公式及 $|a_\nu|\le d_2$ 给出

$$
 |L_\nu^{(e)}(\sigma+it)|
 \ll (1+|t|)^{1-2\sigma}\log^C(2T)\log(2+|t|).
$$

在反射线 $1-\sigma>1$ 上直接对 $\sum d_2(n)n^{\sigma-1}$ 及其对数加权和求界即可。小 $t$ 用整函数在紧集上的有界性。$\Phi_K$ 由分部积分有任意阶 $(1+|t|/K)^{-1}$ 衰减，而区间乘子一致为 $O((1+|t|)^{-1})$。

另一方面，Montgomery–Vaughan 均方不等式给出

$$
 \int_{-H}^{H}|P(\sigma+it)|^2dt
 \ll (H+D)\sum_{d\le D}|f(d)d^{-\sigma}|^2
 \ll (H+D)D\log^5(2T).
$$

这里使用 $d_2(d)^2\le d_4(d)$。在 $|t|\le K$ 及外侧二进带分别应用此界，$\Phi_K$ 的衰减吸收 $(1+|t|)^{-2\sigma}\log(2+|t|)$ 在外带的增长。对 $\log(2x)$ 作 Plancherel，再用 Cauchy–Schwarz，得到 $T^{-1}\int_T^{2T}|S_K(x)|dx\ll\sqrt{D(K+D)}\log^C(2T)$。取 $K=(T^2/D)^{1/3}$ 即得 (4)，有界 $T$ 可增大常数。[Montgomery–Vaughan 的均方不等式](https://personal.science.psu.edu/rcv4/personal/Publications/s2-8-1-73.pdf) · [本域两个 Gamma 因子](cm229-heights-and-zeros.md)

## 删去原和的一段与三次以上素数幂

逐支作卷积，写 $c_{U,\nu}=\delta_1-u_{\nu,\le U}*a_\nu$。因 $u_\nu*a_\nu=\delta_1$，在 $h\le U$ 上 $c_{U,\nu}(h)=0$。对于 $z<\ell\le V<T$，$\delta_1$ 项在所求区间内为零。故 $R_{z,z}-R_{z,V}$ 的每一支正好是 (4) 的负区间和，外层系数为 $f=u_{\nu,\le z}*b_{\nu,(z,V]}$。其长度至多 $D=zV$，并且

$$
 |f(d)|\le 2d_2(d)\sum_{\ell\mid d}\Lambda(\ell)=2d_2(d)\log d .
$$

这里 $|u_\nu|\le d_2$、$|b_\nu|\le2\Lambda$ 来自局部欧拉因子。二进块平均值乘 $T^{-1/2}$ 就控制该块的 $\|\cdot\|_Y$ 费用。取 $D\le\sqrt T$，每块费用为对数幂，整窗共有 $O(\log(2Y))$ 块，得到 (2)。取 $V_0=T^{3/10}$ 的整数部分时，$D\ll T^{7/15}$，而 $(1+7/15)/3-1/2=-1/90$。几何求和给出 (3)。

要处理 $\ell=p^j$、$j\ge3$，只须 $q=p^j\le4T/(z+1)$。对每个 $d\le z$，度数二部分和均方定理给出

$$
 \frac1T\sum_{T\le n<2T}
 \left|\sum_{2n/(rqd)<t\le2n/(qd)} a_\nu(t)\right|\ll (T/(qd))^{1/4} .
$$

此处 $T/(qd)\ge1/4$，小于 1 的情形有绝对常数界。其余情形由 $\int_Z^{2Z}|\sum_{t\le x}a_\nu(t)|^2dx\ll Z^{3/2}$ 得到。用 $d_2(t)\ll t^{1/8}$ 处理单位段上端点的移动。又 $\sum_{p^j\le L,j\ge3}\log p\ll L^{1/3}\log(2L)$，故权以 $(p^j)^{-1/4}$ 后为 $O(L^{1/12}\log(2L))$。对 $d$ 以 $d^{-1/4}|u_\nu(d)|$ 求和，得到该块平均绝对值

$$
 \ll T^{1/4}z^{3/4}(T/z)^{1/12}\log^C(2T)
 =T^{1/3}z^{2/3}\log^C(2T).
$$

$\delta_1$ 项的原始费用至多 $O(T^{1/3}\log(2T))$。代入 $z\le T^{1/6}$ 并加权即得定理。这里各项只对实际支撑 $q\le4T/(z+1)$ 求和，包含端点的最后短块可用全块的绝对值和上界。[Cao–Tanigawa–Zhai，推论 1](https://arxiv.org/abs/1501.04269)

## 改选两个截断与唯一大素因子

令 $J(n)=\tfrac12\sum_{2n/r<q\le2n}(b_+(q)+b_-(q))$。若 $UV\le 2\sqrt T$ 且 $V<T$，则

$$
 \|R_{U,V}-J\|_Y\ll\log^C(2Y). \tag{5}
$$

事实上，$b_\nu=u_{\nu,\le U}*(a_\nu\log)-(u_{\nu,\le U}*b_{\nu,\le V})*a_\nu+b_{\nu,\le V}+c_{U,\nu}*b_{\nu,>V}$。这是由 $a_\nu*b_\nu=a_\nu\log$ 得到的非对称 Vaughan 恒等式。前两项在 (4) 中的外层长度分别为 $U,UV$，第三项支撑在窄区间下方。取 $U=V=z$，同一计算还给出

$$
 \|R_Y-J\|_Y\ll Y^{-1/18}\log^C(2Y), \tag{6}
$$

因为 $z^2\le T^{1/3}$。所以可改取 $U=V=\lfloor T^{1/4}\rfloor$，此时两个因子都严格超过截断。在这一表述中 $c_U$ 也随之改变。式 (2) 给出的删除则逐项保留原来的 $c_z$。

还有一个只含有理合数的等价问题。取 $U=1$，$V=\lfloor2\sqrt T\rfloor$，于是对 $h>1$ 有 $c_1(h)=-a(h)$。从此和删除 $\ell=p^j$、$j\ge2$，归一化绝对值费用为 $O(\log^C(2Y))$。证明与上节相同，只需用 $\sum_{p^j\le L,j\ge2}(\log p)(p^j)^{-1/4}\ll L^{1/4}\log(2L)$，并对 $h=1$ 的端点修正用 $O(\sqrt T\log(2T))$。留下的 $\ell=p>2\sqrt T$ 严格大于 $h$，而且 $p^2>4T$。故 $p$ 是 $q=hp$ 的唯一大素因子，$p\nmid h$，有 $a(h)a(p)=a(q)$。用 $P^+(q)$ 表示最大素因子，得到

$$
 D_Y(n)=-\frac12
 \sum_{\substack{2n/r<q\le2n,\ q\ {\rm composite}\\P^+(q)>2\sqrt T}}
 (1+\chi_{-4}(q))a(q)\log P^+(q),\qquad
 \|R_Y-D_Y\|_Y\ll\log^C(2Y). \tag{7}
$$

特别地，(7) 的累计正部有 $Y^{o(1)}$ 的归一化上界当且仅当原式有这个上界。式 (7) 的负号来自 $c_1(h)$。对式 (1) 保持 $U=z$ 时，定理 1 留下的素数平方仍属于待估部分。

## 固定三阶特征提供的两个变量的关系

令 $\xi$ 是 $F^{\rm gal}/k$ 的三阶理想类特征，$\psi$ 由 $\xi$ 经 $N_{B/k}$ 拉回。按 $k$ 的理想范数收集，$a(h)=\sum_{N\mathfrak h=h}\xi(\mathfrak h)$，$u(h)$ 在同一和中再乘理想 Möbius 函数 $\mu_k$，$b_+(\ell)$ 再按素理想幂加权 $\Lambda_k$。因此 $c_U(h)$ 的理想权为

$$
 \xi(\mathfrak h)\left(\mathbf1_{\mathfrak h=(1)}
   -\sum_{\substack{\mathfrak d\mid\mathfrak h\\N\mathfrak d\le U}}\mu_k(\mathfrak d)\right).
$$

二次扭转给每项再乘 $\chi_{-4}(N\mathfrak h)$ 或 $\chi_{-4}(N\mathfrak l)$。对只依赖这两个范数的可分离矩形权，按三阶特征值 $1,\omega,\bar\omega$ 收集两变量各自除去 $\xi$ 的权，记作 $(A_0,A_1,A_2)$ 和 $(B_0,B_1,B_2)$。$k/\mathbb Q$ 的共轭保持范数且反转该三阶商，故 $A_1=A_2$、$B_1=B_2$。于是单支的双变量特征和恰为

$$
 (A_0-A_1)(B_0-B_1). \tag{8}
$$

这是 $\xi(\mathfrak h\mathfrak l)=\xi(\mathfrak h)\xi(\mathfrak l)$ 及 $\omega+\bar\omega=-1$ 的直接结果。矩形上的两份类差仍需算术估计。把 $\xi$ 和 $\bar\xi$ 再平均一遍会得到同一个乘积。Lean 对三类的配对恒等式作了整系数的符号核验。[类域互反与范数](https://www.jmilne.org/math/CourseNotes/CFT.pdf) · [Lean 三类配对](CM229Frobenius.lean)

在 Dirichlet 级数上还能看到剩余项为何保留每个零点。对固定有限截断置 $U_\nu(s)=\sum_{d\le U}u_\nu(d)d^{-s}$、$B_{\nu,V}(s)=\sum_{\ell\le V}b_\nu(\ell)\ell^{-s}$。式 (1) 单支卷积的级数是

$$
 (1-U_\nu(s)L_\nu(s))\big(-L_\nu'(s)/L_\nu(s)-B_{\nu,V}(s)\big). \tag{9}
$$

在 $L_\nu$ 的每个非平凡零点 $\rho$，第一因子精确地等于 1，整个乘积的留数等于零点重数的负值。两支相加时共同零点的重数也相加。对随块变化的截断，窗口结论由 (5)、(6) 与前篇的完整窗口定理传递。伴篇直接把这批零点在一个给定窗口上的贡献分离出来。[双线性和的有限高度公式](cm229-type-ii-zeros.md) · [原窗口的精确指数](cm229-cubic-hecke-carries.md)

## 完整窗口的逐项复算

在 $Y=4096$、$X=262144$、$m=82$ 时，原定义的 $z$ 在六个二进块上都等于 1。为核验两个截断实际增长的情形，另算一遍 $z=\lfloor T^{1/6}\rfloor$，六个值依次为 $4,4,5,5,6,7$。两行均取 $V=\lfloor\sqrt T/z\rfloor$。各列已除以 $W$，除第一列数值取正部外，其余各列都在每个 $n$ 先取绝对值。表值是真值的严格下界，加 $10^{-6}$ 即为严格上界。

| $z$ 的选择 | $R_{z,z}$ 正部 | 删去的 $z<\ell\le V$ | $\ell>V$ 的三次及更高幂 | 与四分之一次幂截断之差 | 与最大素因子和 (7) 之差 |
| --- | ---: | ---: | ---: | ---: | ---: |
| 原式 | 0.868852 | 0.745264 | 0.042464 | 1.792361 | 0.913569 |
| $\lfloor T^{1/6}\rfloor$ | 0.653600 | 0.276567 | 0.079605 | 0.911352 | 1.385341 |

程序用全部乘积的前缀和分别构造原和、被删去的两部分及留下的素数和素数平方部分，对每个前缀检验相加后的区间相交。式 (7) 另从每个整数的最大素因子构造，并在每个符合条件的合数上用整数检查 $a(q)=a(q/p)a(p)$。最后用整数平方根计算 $k_n$，用实球算术确定每个窄区间的左端并累加所有 $n$。从仓库根目录执行下列命令，依赖为 python-flint 0.9.0。[完整窗口核验程序](verify_cm229_type_ii.py) · [双线性卷积程序](verify_cm229_bilinear.py) · [系数与零点数据](cm229-heights-and-zeros.md)

~~~sh
python -B ERH/verify_cm229_type_ii.py
~~~
