---
title: "用混合 CM 型构造达到最小奇异值代价的整同构"
date: "2026-11-27"
abstract: "对判别式 229 的三次域，从单个混合 CM 型及两个二次特征读出相对三次分裂，显式求出窗口中每条分式理想的正有理生成元。相应对角映射在任意靶秩的所有整单射中达到最小代价。只读取大于平方根的素数也只多付 O(W log Y)。若要求矩阵系数本身为代数整数，最小代价渐近于 (8 log 2 / 3) W Y^(3/4)。"
---

# 用混合 CM 型构造达到最小奇异值代价的整同构

令 $F=\mathbb Q(\alpha)$，$\alpha^3-4\alpha+1=0$，$N$ 为其正规闭包，$E=F(i)$，$M=N(i)$，$B=\mathbb Q(\sqrt{229},i)$。在 $M/B$ 这个循环三次扩张中，每个未分歧有理素数的相对分裂行为，都可以由 $E$ 的一个固定混合 CM 型及两个二次符号读出。把这些数据写入前篇的整值多项式格，可逐项给出一个 $\mathcal O_B$-线性同构。它在每个完整窗口上都达到最小正奇异值代价。[窗口格与度量的原定义](cm229-binomial-slopes.md)

这个最小值的全高度次幂界仍由 $E$ 的两个非阿贝尔 $L$ 因子决定。用已有的无零半平面，显式映射的代价从 $\asymp WY^{3/4}$ 降至 $O(WY^{0.56})$。若 $E$ 满足 ERH，伴篇证明最小代价渐近于 $c_{229}W\log Y$，其中 $c_{229}>0.61070$。[无零半平面定理 7.33](https://hexagonmath.org/pdf/2611.00741v1) · [临界线情形的对数渐近式](cm229-critical-logarithm.md)

## 一个混合型读出相对 Frobenius

对 $p\nmid458$，选定一个上方素处，写有理 Frobenius 为 $(g,t)\in S_3\times\{\pm1\}$。其中 $t=(-4/p)$，$\operatorname{sgn}(g)=(229/p)$。在 $E$ 的三个实嵌入位置中固定一个 $j$，让 CM 型 $\Phi_j$ 只在该位置取正号。取带完整 $\mathcal O_E$ 作用的相应 CM 三维阿贝尔簇，在获得良约化的有限扩张上计算其几何 $a$ 数 $a_j(p)$。若各位置的符号为 $s_r$，约化的 Dieudonné 模给出 $a_j(p)=\#\{r:s_{gr}\ne t s_r\}$。当 $t=+1$ 时，固定、移动 $j$ 分别给出 $0,2$。当 $t=-1$ 时，这两个值分别为 $3,1$。[圆周字与 CM 约化的计算](https://arxiv.org/abs/1809.10368) · [该混合型的局部长度证明](cm229-local-length-window.md)

设 $\chi=(229/p)$。$B$ 中的剩余次数 $f_B$ 在 $\chi=t=1$ 时等于 1，其余情形等于 2。$B$ 中 $p$ 上方素理想的相对 Frobenius 为 $(g,t)^{f_B}$。若 $\chi=-1$，则 $g$ 是换位，其平方为一。若 $\chi=1$，则 $g$ 是恒等或三循环，固定一个 $j$ 就足以区分二者。因此相对完全分裂的判据为

$$
 \chi=-1\quad\text{或}\quad(t=1,\ a_j(p)=0)
              \quad\text{或}\quad(t=-1,\ a_j(p)=3). \tag{1}
$$

这个答案与所选上方素处和少数符号位置无关。相对分裂时 $f_M=f_B$，惰性时 $f_M=3f_B$。Lean 枚举了六种置换、两个高斯符号、三个少数符号位置，将 (1) 与相对 Frobenius 在六个嵌入上的作用逐一比较。[有限群作用的 Lean 证明](CM229Frobenius.lean)

## 分式理想的正有理生成元

令 $I_{K,n}$ 为次数至多 $n$ 的整值多项式的 $X^n$ 系数理想。以 $\mathfrak b_K(n)$ 记乘法 $I_{K,n}^{\otimes2}\longrightarrow I_{K,2n}$ 的余核的 Fitting 理想，置 $\mathcal L_n=\operatorname{Nm}_{M/B}\mathfrak b_M(n)\,\mathfrak b_B(n)^{-1}$。对整数 $q\ge2$ 定义

$$
 \epsilon_q(n)=\lfloor 2n/q\rfloor-2\lfloor n/q\rfloor,\qquad
 S_f(p,n)=\sum_{j\ge1}\epsilon_{p^{fj}}(n).
$$

局部阶乘公式给出 $v_{\mathfrak p}(\mathfrak b_K(n))=S_{f_K}(p,n)$。记 $e_p=e_B(p)=e_M(p)$。对 $\mathfrak l\mid p$，相对扩张无分歧，$\sum_{\mathfrak P\mid\mathfrak l} f_{\mathfrak P/\mathfrak l}=3$，故 $v_{\mathfrak l}(\mathcal L_n)=3S_{f_M}(p,n)-S_{f_B}(p,n)$。[素理想排序与阶乘公式](https://collaborate.princeton.edu/en/publications/p-orderings-and-polynomial-functions-on-arbitrary-subsets-of-dede/)

仅有 $2,229$ 需要检验除以 $e_p$ 后的整性。模 2 的三次多项式分解为一次和不可约二次因子，$\mathbb Q(i)$ 在 2 分歧，$\mathbb Q(\sqrt{229})$ 在 2 惰性，两域 $M,B$ 的 $(e,f)$ 均为 $(2,2)$。在 229，三次域的判别式只有一次因子 229，正规闭包中的惯性群由换位生成，其在 $S_3$ 中的正规化子就是它自身，故剩余次数为 1。又 $229\equiv1\bmod4$，两域均有 $(e,f)=(2,1)$。两处都在 $M/B$ 完全分裂，$\mathcal L_n$ 的局部赋值为偶数。其余素数有 $e_p=1$。

**定理 1。** 对每个 $n$，理想 $\mathcal L_n$ 都由如下正有理数生成。

$$
 d_p(n)=\frac{3S_{f_M}(p,n)-S_{f_B}(p,n)}{e_p},\qquad
 r_n=\prod_{p\le 2n}p^{d_p(n)},\qquad
 \mathcal L_n=r_n\mathcal O_B. \tag{2}
$$

事实上，每个 $\mathfrak l\mid p$ 都有 $v_{\mathfrak l}(r_n)=e_pd_p(n)$，这就证明了理想等式。例如 $r_5=1/3$，$r_{15}=4/255$，$r_{20}=37^2/5$。第一例的源格为 $3\mathcal O_B$，乘 $1/3$ 正好将它送到 $\mathcal O_B$。分母记录了惰性处的负赋值，可直接用于整映射。

## 奇异值的逐窗最小值

取实数 $Y\ge2$，$X=Y^{3/2}$，$W=\lceil Y^3\rceil$，$k_n=\lfloor W/n^{3/2}\rfloor$。下文的 $n$ 遍历 $Y\le n<X$ 的所有整数，与前篇对连续闭窗口的离散化相同。源格为 $\mathcal W_Y=\bigoplus_n\mathcal L_n^{-k_n}$，度量由分式理想在 $B$ 的标准绝对值及正交直和给定。靶取任意秩的标准 Hermitian 格 $\mathcal O_B^m$。对整单射 $J$ 记

$$
 C(J)=\sum_{v\mid\infty}[B_v:\mathbb R]\sum_l[\log s_l(J_v)]_+,\qquad
 [z]_+=\max(z,0).
$$

**定理 2。** 在每个这样的窗口取

$$
 (J_*x)_n=r_n^{k_n}x_n,\qquad
 C_*(Y)=\min_{m,J} C(J)=4\sum_n k_n[\log r_n]_+. \tag{3}
$$

式 (2) 说明 $J_*$ 是同秩整同构。两个复位上各坐标的奇异值相同，局部次数之和为 4，故它的代价就是 (3)。为证最小性，取所有 $r_n^{k_n}>1$ 的列。任意整单射的第 $n$ 列均等于 $r_n^{k_n}$ 乘一个整列。选定各列除去这些因子后仍然线性无关，所以存在一个非零方子式 $a\in\mathcal O_B$。在各无穷位，原子式的绝对值至多为相应个数的最大奇异值之积。取对数并按局部次数求和，又 $|N_{B/\mathbb Q}a|\ge1$，即得 (3) 的下界。空列集的情形由 $C(J)\ge0$ 给出。这个论证也适用于任意指定的正有理列因子。

若 $X$ 恰为整数且将上端点也列为一条格线，则 (2)、(3) 对该列仍然成立。因 $|\log r_X|=O(X)$，其额外代价为 $O(W/\sqrt X)$，以下全部渐近式对两种端点约定相同。

令 $L=L_+L_-$，其中 $\zeta_M=\zeta_B L^2$，$\beta_*$ 为 $L$ 的非平凡零点实部的上确界。因 $h(n)=\tfrac12\log N_B\mathcal L_n=2\log r_n$，将 (3) 代入前篇的窗定理，得到

$$
 \lim_{Y\to\infty}\frac{\log(1+C_*(Y)/W)}{\log Y}
      =\delta,\qquad \delta=\tfrac32(\beta_*-\tfrac12). \tag{4}
$$

该式对所有实数左端点成立。Hu 的 $7/8$ 无零半平面给出 $C_*(Y)\ll_\varepsilon WY^{9/16+\varepsilon}$，只需将其 Hecke 定理用于 $M$ 的平凡特征。采用 Liu 预印本定理 7.33，可得 $C_*(Y)=O(WY^{0.5599989561688578})$，这里将指数向上舍入。取得任意远窗口上的 $WY^{o(1)}$ 仍等价于 $\beta_*=1/2$，也就是本例的 $\operatorname{ERH}(E)$。上述两条数值上界均使用外部无零定理。[完整窗口的精确指数定理](cm229-binomial-slopes.md) · [Hu 定理及形式化说明](https://hexagonmath.org/2610.00216v2) · [Liu 定理 7.33](https://hexagonmath.org/pdf/2611.00741v1)

## 只读取平方根以上的素数

令 $b_n=\binom{2n}{n}$。定义正有理数 $s_n$，在 $p\le\sqrt{2n}$ 取 $v_p(s_n)=2v_p(b_n)$，在 $p>\sqrt{2n}$ 取 $v_p(s_n)=d_p(n)$。写 $S_0,S$ 分别为在 $B,M$ 完全分裂的未分歧有理素数。对于后一区间的未分歧素数，

$$
 d_p(n)=(3\mathbf1_S(p)-\mathbf1_{S_0}(p))\epsilon_p(n). \tag{5}
$$

因而 $s_n$ 的既约分母是平方自由的，仅由 $S_0\setminus S$ 中大于 $\sqrt{2n}$ 且产生进位的素数组成。此处只需对 $S_0$ 中的素数读取 (1) 的混合型数据。分歧素数若在此区间，按 (2) 取值即可。

**定理 3。** 逐列乘 $s_n^{k_n}$ 得到整单射 $J_{\rm large}$，且一致地有

$$
 0\le C(J_{\rm large})-C_*(Y)\ll W\log(2Y). \tag{6}
$$

事实上，(2) 中每一个进位的系数介于 $-1$ 与 2 之间，所以 $s_n/r_n$ 为正整数，且

$$
 0\le\log(s_n/r_n)\le3\sum_{p\le\sqrt{2n}}v_p(b_n)\log p
   \le 3\pi(\sqrt{2n})\log(2n)\ll\sqrt n .
$$

最后用了 Chebyshev 上界。利用正部的 1-Lipschitz 性，并在整窗上求和 $k_n\sqrt n\le W/n$ 即得 (6)。特别地，小素数和全部高次素数幂可以按此方法统一处理。对于 (5) 的带符号和

$$
 D_n=\sum_{\sqrt{2n}<p\le2n}
       (3\mathbf1_S(p)-\mathbf1_{S_0}(p))\epsilon_p(n)\log p
$$

有 $\log r_n=D_n+O(\sqrt n)$，分歧素数的贡献包含在误差内。于是 $C_*(Y)/(4W)$ 与 $W^{-1}\sum_n k_n[D_n]_+$ 之差为 $O(\log(2Y))$。原目标归结为对每个 $\eta>0$，在任意远的完整窗口上取得 $\sum_n k_n[D_n]_+\le WY^\eta$。

## 清除分母的代价

将 $r_n=u_n/v_n$ 写为互素正整数之比。由 (5) 和小素数估计，

$$
 \log u_n=2\sum_{p\in S}\epsilon_p(n)\log p+O(\sqrt n),\qquad
 \log v_n=\sum_{p\in S_0\setminus S}\epsilon_p(n)\log p+O(\sqrt n).
$$

对 Chebotarev 密度为 $d$ 的 Frobenius 共轭类并集 $A$，有 $\sum_{p\in A}\epsilon_p(n)\log p=2d(\log2)n+o(n)$。确实，进位素数位于不交区间 $(2n/(2j+2),\,2n/(2j+1)]$。先固定 $J$ 对前 $J$ 段用 Chebotarev 素数定理，余段由 Chebyshev 界贡献 $O(n/J)$，再令 $J\to\infty$ 即得交错调和级数。$S$ 与 $S_0\setminus S$ 在全部有理素数中的密度分别为 $1/12$ 和 $1/6$，从而

$$
 \log u_n\sim\log v_n\sim \frac{\log 2}{3}\,n . \tag{7}
$$

这里 $u_n$ 和 $v_n$ 来自同一个有理生成元。Chebotarev 的密度只给出线性主项抵消，余项的正部正是 (3) 所需的量。[Lagarias–Odlyzko 的 Chebotarev 素数定理](https://aareyanmanzoor.github.io/assets/articles/lagarias-odlyzko.pdf)

**定理 4。** 在原分式理想坐标中，若另要求矩阵每一项属于 $\mathcal O_B$，则任意靶秩下的最小奇异值代价为

$$
 C_{\rm int}(Y)=4\sum_n k_n\log u_n
       \sim \frac{8\log2}{3}W Y^{3/4}. \tag{8}
$$

第 $n$ 列的每一项必须属于 $(r_n^{k_n})\cap\mathcal O_B=u_n^{k_n}\mathcal O_B$。以这些列因子使用定理 2 的子式下界，对角映射取等。又 $\sum_n k_n n\sim 2W Y^{3/4}$，其中取整误差 $O(X^2)=O(W)$，得到 (8)。

更一般地，令 $0\le\ell_n\le k_n$ 为整数，只允许第 $n$ 列的分母整除 $v_n^{\ell_n}$。两个分式理想的交给出列因子 $u_n^{k_n}/v_n^{\ell_n}$，所以最小值恰为 $4\sum_n[k_n\log u_n-\ell_n\log v_n]_+$。若 $\ell_n=\lfloor\lambda k_n\rfloor$，其中 $0\le\lambda<1$ 固定，(7) 给出最小值渐近于 $(1-\lambda)(8\log2/3)WY^{3/4}$。这说明保留每个分母的任一固定比例仍有四分之三次幂的代价。定理 3 则指出，可将分母集中在平方根以上的特定素数处。

允许比例随 $Y$ 变化，还能算出接近最优值所需的精度。令 $0\le d_Y\le1$，取 $\ell_n=\lfloor(1-d_Y)k_n\rfloor$，用 $C_d$ 表示此时的最小值，$R=4\sum_n k_n\log v_n$，$C_-=4\sum_n k_n[-\log r_n]_+$。逐项用 $a-[-b]_+\le[a+b]_+\le a+[b]_+$（$a\ge0$），并以 $\sum_n\log v_n=O(W)$ 处理取整，有

$$
 \max(C_*,d_Y R-C_-)\ \le C_d\ \le C_*+d_Y R+O(W),\qquad
 R\sim (8\log2/3)WY^{3/4}. \tag{9}
$$

例如 $d_Y=Y^{-\gamma}$、$\gamma\ge0$ 时，$C_d/W$ 的增长指数，即 $\log(1+C_d/W)/\log Y$ 的极限，为 $\max(\delta,3/4-\gamma)$。前篇的点态估计 $|\log r_n|\ll_\varepsilon n^{\beta_*+\varepsilon}$ 给出 $C_-/W\ll_\varepsilon Y^{\delta+\varepsilon}$（重选 $\varepsilon$）。若 $3/4-\gamma>\delta$，(9) 的 $d_Y R$ 给出主项。其余情形用 $C_*$ 的 (4) 即得。若 $\operatorname{ERH}(E)$ 成立，则 $C_*,C_-=O(W\log Y)$，后一个界由伴篇对 $|G|$ 的均值及相同的取整估计得到。此时 (9) 还给出 $C_d=O(W\log Y)$ 当且仅当 $d_Y=O(Y^{-3/4}\log Y)$。[临界线上的均方估计与正部均值](cm229-critical-logarithm.md)

## 三个完整窗口的复算

各列都除以 $2W$。表中每个数都是严格区间的下端点，上端点由它加 $10^{-6}$ 得到。最后一列取 (8) 的最小值。

| $[Y,Y^{3/2}]$ | 最优整同构 | 只读取大素数 | 矩阵系数均为代数整数 |
| --- | ---: | ---: | ---: |
| $[1024,32768]$ | 1.953553 | 11.287457 | 131.589139 |
| $[4096,262144]$ | 5.190081 | 19.506001 | 412.129950 |
| $[10000,1000000]$ | 8.863359 | 27.461267 | 835.188995 |

程序按 $\epsilon_q(n)$ 的全部跳点累加，先对同一素数合并各次幂，再对局部赋值取正部以计算 $u_n$。对 $r_n$ 则先合并全部素数的对数，再取正部。$k_n$ 用整数平方根精确计算，最后用 160 位 Arb 球算术验证九个开区间。另在八个次数上按阶乘赋值直接重建有理数进行交叉检查。复算需将下列程序和 CSV 放在同一目录，使用 python-flint 0.9.0。[本篇程序](verify_cm229_optimal_map.py) · [剩余次数与 Artin 迹](verify_cm229_binomial.py) · [素数筛](verify_cm229_split_primes.py) · [零点区间](cm229-zeros64.csv)

~~~sh
python -B verify_cm229_optimal_map.py
lean +leanprover/lean4:v4.35.0 CM229Frobenius.lean
~~~
