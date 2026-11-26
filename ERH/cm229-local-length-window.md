---
title: "用 CM 约化的局部长度识别 x/48 素数，并消去几何斜坡的线性项"
date: "2026-11-27"
abstract: "由三个混合 CM 型的 a 数、一个 CM 椭圆曲线的 a 数及另一个的 Frobenius 平方，构造长度恰为目标素数指示函数的局部模。两个尺度的斜坡之差保留正权并消去线性项，其正部积分足以检验原完整窗口。高阶矩给出把稀少宿主上的几何界转给平均的定量代价。"
---

# 用 CM 约化的局部长度识别 x/48 素数，并消去几何斜坡的线性项

令 $f=X^3-4X+1$，$F=\mathbb Q[X]/(f)$，$N$ 为其分裂域，$E=F(i)$，$M=N(i)$。设 $S$ 是在 $M$ 完全分裂的未分歧有理素数。对两个辅助素数 $q,d$，记 $I_{q,d}(x)$ 为 $S$ 中同时在 $K_q=\mathbb Q(\sqrt{-q})$ 和 $K_d=\mathbb Q(\sqrt{-d})$ 惰性的 $p\le x$ 的 $\log p$ 之和。前篇的逐步对数权抽样给出 $A_Y(x)=\mathbb E_Y I_{q,d}(x)$。问题是在任意远的完整窗口上证明

$$
 P(Y)=\int_Y^{Y^{3/2}}[A_Y(x)-x/48]_+\,\frac{dx}{x^{3/2}}\le Y^\eta
 \qquad(\text{每个 }\eta>0),\qquad [v]_+=\max(v,0).
$$

抽样中固定 $0<e<1/52$，令 $z=Y^{e/5}$，$P_0$ 为不超过 $z$ 的奇素数之积。取 $U=\exp(3z/4)$、$D=\exp(3z(\log Y)^2/4)$。先在 $(U,2U]$ 中取 $q\equiv3\bmod8$，再在 $(D,2D]$ 中取 $d\equiv7\bmod16,\ d\equiv-1\bmod qP_0$，两步分别用素数的对数权归一化。下面把这笔超额关联到可以在约化上逐项算出的局部长度，并给出保持全部 $\eta$ 的两尺度积分判据。该判据的全局上界仍待证明。[原窗口定义与等价定理](cm229-one-sided-primes.md) · [辅助素数的转移证明](winpair-inert-transfer.md)

## 一个未分歧素数贡献的长度

固定 $p\nmid 2\cdot229qd$。在 $N$ 中选择一个 $p$ 上方的素处，记 $g\in S_3$ 为其 Frobenius。对 $FK_q$ 的三个混合 CM 型 $\Phi_j$，令第 $j$ 个实嵌入取正号，另两个取负号。各取一个带完整整数环作用、具有该型的 CM 三维阿贝尔簇 $A_j$。按同一标号及同一代数闭包上的素处约化，可在有限扩张后取得良约化并定义所有复乘。另取分别具有 $K_d$ 和 $\mathbb Q(i)$ 的完整复乘的椭圆曲线 $C_d,C_i$。

在 $\kappa=\overline{\mathbb F}_p$ 上记 $\mathcal D(A)$ 为 $p$-可除群的 Dieudonné 模，$F_p,V_p$ 为其 Frobenius、Verschiebung。置

$$
 Q(A)=\mathcal D(A)/(F_p\mathcal D(A)+V_p\mathcal D(A)),\qquad
 U(C)=\operatorname{im}(\overline F_p^{\,2}\mid\mathcal D(C)/p).
$$

因为 $F_pV_p=p$，$Q$ 被 $p$ 杀掉，它的 $\kappa$ 维数是 $a$ 数。用这两个模在 $\kappa$ 上构造

$$
 \mathcal T_{q,d,p}=
 \left(\bigotimes_{j=1}^3\bigwedge\nolimits^3 Q(A_j)\right)
       \otimes Q(C_d)\otimes U(C_i).
$$

**局部定理。** $\mathcal T_{q,d,p}$ 作为 Witt 环上的模长为

$$
 \ell_p=\mathbf1_{\{g=1,\ p\equiv1\bmod4,\ (-q/p)=(-d/p)=-1\}} .
 \tag{1}
$$

故 $\sum_{p\le x}\ell_p\log p=I_{q,d}(x)$，在分歧处约定 $\ell_p=0$。三种混合型一起取使公式与所选的上方素处无关。

**证明。** 未分歧处的完整 CM 作用把 Dieudonné 模分为六条秩一特征线。Frobenius 沿嵌入置换 $\sigma:(r,s)\mapsto(gr,ts)$，其中 $t=(-q/p)$，各边的乘数是单位或 $p$ 乘单位。指数为 1 的边由 CM 型给出，取对偶约定会互换两个符号。在每个圆周上，对 $Q$ 有贡献的位置恰是相邻边指数由 1 变 0 的位置。于是 $a_j=|\Phi_j\setminus\sigma\Phi_j|$。对符号向量 $s_r$，这就是满足 $s_{gr}\ne ts_r$ 的 $r$ 的个数。得到下表。[Smirnov–Zaytsev 关于 CM 约化的圆周字及 a 数，第 2–3 节](https://arxiv.org/abs/1809.10368)

| $g$ 对少数符号位置 $j$ 的作用 | $K_q$ 分裂 | $K_q$ 惰性 |
| --- | ---: | ---: |
| 固定 $j$ | $a_j=0$ | $a_j=3$ |
| 移动 $j$ | $a_j=2$ | $a_j=1$ |

$\bigwedge^3Q(A_j)$ 的维数等于 $\mathbf1_{\{a_j=3\}}$。三个位置都固定等价于 $g=1$。对 CM 椭圆曲线，惰性时 $a=1$、$\overline F_p^2=0$，分裂时 $a=0$、$\overline F_p^2$ 的秩为 1。它们也可用同样的二条边计算，遂得 (1)。还可任取其中两个不同的混合型，它们都在同一 Galois 轨道，可取互为共轭的簇。固定其中两个位置已迫使 $g=1$，故删去第三个顶外幂仍得 (1)。单独一个位置会让一个换位通过。若已在算术上限制 $(229/p)=+1$，$g$ 只能是单位元或三循环，此时一个混合型就足够。

取 $q=11,d=263$ 作局部例子。193 与 277 均为 $1\bmod4$，且两宿主均惰性。模 193 的三个根为 $6,36,151$，模 277 无根。在这两处八种 CM 型的 Newton 斜率都全是 $1/2$，三个混合型的 $a$ 数却分别是 $(3,3,3)$ 和 $(1,1,1)$。公式 (1) 给出长度 1 和 0。程序核验了 $6\times2^3=48$ 种 Frobenius 与二次特征组合，以及这两个素数和一个换位例子。Dieudonné 模与 CM 型的对应使用上引的约化定理。[Lean 有限计算](CM229Frobenius.lean) · [此前的 Newton 斜率表](cm229-frobenius-and-heights.md)

## 两个尺度保留正权

原 Lean 的算术请求逐项使用 $(t/p^j-w)_+$。取基本斜坡 $r(v)=(v-1)_+$，定义

$$
 K(v)=r(2v)-2r(v)=
 \begin{cases}0&v\le1/2,\\2v-1&1/2<v<1,\\1&v\ge1.\end{cases}
$$

记 $R_{q,d}(t)=\sum_p\ell_p(\log p)r(t/p)$，$J_{q,d}(t)=R_{q,d}(2t)-2R_{q,d}(t)$。每个 $p$ 的新权在 $[0,1]$ 内，且有精确恒等式

$$
 J_{q,d}(t)=\sum_p\ell_p(\log p)K(t/p)
            =2t\int_t^{2t} I_{q,d}(s)\frac{ds}{s^2}. \tag{2}
$$

对阶梯函数逐个素数积分即得。$Ct$ 型的未知高度项在差中消失。将 $s/48$ 代入右侧，留下的几何主项是 $(\log2)t/24$。小素数的线性中心化也不再需要。[原 ASK、MAIN 和各素数的斜坡定义](https://github.com/leanprover-community/riemann-hypothesis/blob/9c9f665ed79ce0fa8da61bd803b7a4bf1ad62ed0/src/Formalism/SizeWindowBrackets/FirstPrimeNeSecondPrime.lean#L118)

**窗口定理。** 记 $J_Y=\mathbb E_Y J_{q,d}$。原来的 $P(Y)\le Y^\eta$ 条件等价于

$$
 \forall\eta>0,\quad\text{存在任意大的 }Y,\qquad
 \int_Y^{Y^{3/2}}\left[J_Y(t)-\frac{\log2}{24}t\right]_+
                        \frac{dt}{t^{3/2}}\le Y^\eta . \tag{3}
$$

在 (2) 中对 $t,2t$ 使用同一批宿主。原抽样转移在 $[Y,Y^2]$ 上同样成立，所以覆盖这里用到的 $2Y^{3/2}$。设 $L=L_+L_-$ 是 $E$ 的两个非阿贝尔因子之积，$-L'/L(s)=\sum b(n)n^{-s}$，$T(s')=\sum_{n\le s'}b(n)$。已知 $T(s')=24(A_Y(s')-s'/48)+O(\sqrt{s'}\log^2(2s'))$。对两边用 (2) 的积分算子，得到的误差在 (3) 中总共花 $O(\log^3(2Y))$。而 $T$ 的变换 $H(t)=\sum b(n)K(t/n)$ 满足

$$
 \int_1^\infty H(t)t^{-s-1}\,dt
      =-\frac{L'}L(s)\frac{2^s-2}{s(s-1)}
       \qquad(\Re s>1). \tag{4}
$$

这里 $b(1)=0$，故积分下端没有截断修正。乘子在 $s=1$ 可去，在 $0<\Re s<1$ 无零点。$L$ 的高度 64 证书排除了整个实临界带上的实零点，故变换的所有极点均符合原窗口引理的留数条件。将 (3) 代入该引理得到 $L$ 在 $\Re s>1/2$ 无零，函数方程及前篇的判据给出原目标。反向，由原条件推出 $\operatorname{ERH}(E)$，再对 $T$ 的平方根级估计使用 (2) 即得 (3)。这条推理经过零点，两个正部积分可以在某个具体窗口上取不同值。[原窗口引理](https://github.com/leanprover-community/riemann-hypothesis/blob/9c9f665ed79ce0fa8da61bd803b7a4bf1ad62ed0/src/Development/AutomorphicArithmetic/DedekindLFunctions/LandauTheorems/landauWindow.lean#L125) · [实零点证书与分解](cm229-heights-and-zeros.md)

整数重数也能保留 (2)。给定整数 $m\ge1$，以 $\lfloor m r(2t/p)\rfloor-2\lfloor m r(t/p)\rfloor$ 为局部长度的重数。它非负，与 $mK(t/p)$ 的差的绝对值小于 2。累加后除以 $m$，代价为 $O(t/m)$。取 $m=\lceil\sqrt t\rceil$，整窗积分代价为 $O(\log Y)$。要用斜率不等式估计这个差，须有两尺度之间的格映射或分母的下界。单独给 $R(2t)$ 和 $R(t)$ 两个上界时，第二项的负系数妨碍直接相减。[一个在特征零中由乘法实现两尺度比较的格](cm229-binomial-slopes.md)

## 能从多稀少的宿主转回原积分

几何构造若只能在一部分宿主上实现，需要控制对这部分宿主作选择的代价。令 $\alpha$ 为所选宿主在原条件抽样中的总概率。存在常数 $c,C>0$，使对充分大的 $Y$ 且 $\log(2/\alpha)\le c z/\log Y$，任何概率至少为 $\alpha$ 的宿主集合中，至少有一对满足

$$
 \int_Y^{Y^{3/2}} |J_{q,d}(t)-J_Y(t)|\frac{dt}{t^{3/2}}
 \le C(\log Y)^{3/2}\sqrt{\log(2/\alpha)}. \tag{5}
$$

**证明。** 对 $z<p\le2Y^{3/2}$，置 $u_p=(-q/p)$、$v_p=(-d/p)$，$W_p=(-u_p-v_p+u_pv_p)/4$。它是双惰性指示函数减 $1/4$。在 $2k\le c_0 z/\log Y$ 时，至多 $2k$ 个不同素数的符号单项式 $\prod_{p\in S_1}u_p\prod_{p\in S_2}v_p$，只要 $S_1,S_2$ 至少一个非空，其期望为 $O(e^{-c_1 z})$。若 $S_2$ 非空，先对给定 $q$ 的 $d$ 取期望，否则对 $q$ 取期望。写 $m$ 为测试素数之积，$h=16qP_0$。Dirichlet GRH 特征估计除以抽样总权，在 $d$ 这一步给出 $O(\varphi(h)D^{-1/2}\log^2(2Dhm))$，在 $q$ 这一步给出 $O(U^{-1/2}\log^2(2Um))$。因 $\log h=O(z)$ 且 $\log m\le2k\log(2Y^{3/2})=O(z)$，两项都满足所需指数界。辅助素数远大于所有测试素数，$p\le z$ 在每对宿主中贡献均为零。[同一抽样的特征展开与分母估计](winpair-inert-transfer.md)

在任意 $2k$ 个 $W_p$ 的乘积中，用 $u_p^2=v_p^2=1$ 合并重复变量。展开系数绝对值之和至多 $(3/4)^{2k}$，故该乘积的期望与独立公平符号模型相差 $O(e^{-c_1z})$。取系数 $a_p=\mathbf1_S(p)(\log p)K(t/p)$，有 $\sum|a_p|\ll t$、$\sum a_p^2\ll t\log(2t)$。独立模型中的中心化 Bernoulli 和用指数矩界有 $2k$ 次矩至多 $(Ck\sum a_p^2)^k$。原模型的附加项至多 $O(e^{-c_1z}(Ct)^{2k})$，取小 $c_0$ 即可使它小于 1。用真实均值中心化至多将 $L^{2k}$ 范数乘 2，遂得 $\|J_{q,d}(t)-J_Y(t)\|_{2k}\ll\sqrt{k t\log(2t)}$。积分 Minkowski 不等式给出 (5) 左端的 $L^{2k}$ 范数 $\ll\sqrt{k}(\log Y)^{3/2}$。取 $k=\lceil\log(2/\alpha)\rceil$，由 Markov 不等式将超出 (5) 右端的概率压到 $\alpha/2$ 以下。

特别地，若 $\log(2/\alpha)=Y^{o(1)}$，(5) 只花 $Y^{o(1)}$。若对每个 $\eta>0$，在任意远的窗口上，能对这一总权的宿主把 (3) 中**单独宿主**的积分压到 $Y^{\eta/2}$，正部的 1-Lipschitz 性就给出所需的 $P(Y)$ 界。原 Lean 的宿主算术数据以有限交换群为参数。本题的局部模还记录 $S_3$ 的作用，式 (1) 给出了各素数所需要的整数长度。尚缺的是由全局格比较给出这个宿主集合及其高度预算。[原宿主数据的群假设](https://github.com/leanprover-community/riemann-hypothesis/blob/9c9f665ed79ce0fa8da61bd803b7a4bf1ad62ed0/src/Formalism/SizeWindowBrackets/FirstPrimeNeSecondPrime.lean#L83)
