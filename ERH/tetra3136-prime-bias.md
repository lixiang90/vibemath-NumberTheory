---
title: "一个四面体数域的 123 个零点与完全分裂素数的负偏差"
date: "2026-11-26"
abstract: "给出一维仿射群扩张的全部中间域中非阿贝尔因子的重数公式。对判别式 3136 的四次域，严格验证三维因子在正高度 64 内的全部 123 个零点。证明完全分裂素数相对模 7 素数的高斯加权负次项，刻画渐近式成立的全部正斜率，并从临界斜率的下极限读出极值零点的总重数。"
---

# 一个四面体数域的 123 个零点与完全分裂素数的负偏差

取 $f(X)=X^4-2X^3+2X^2+2$，$K=\mathbb Q[X]/(f)$，$M$ 为其分裂域。伽罗瓦群为 $A_4$。域 $K$ 的 Dedekind ζ 函数除去 $\zeta(s)$ 后，剩下一个三维因子 $L(s)$。在十二次闭包 $M$ 中，这个因子出现三次。以下两个结果分别涉及零点和有理素数。

**零点。** 在 $0<\Re s<1$、$0<\Im s\le64$ 中，$L(s)$ 恰有 **123 个零点**，全部是中心线上的单零点。实轴上没有非平凡零点。结合 Alpöge–Claude 的 Dirichlet GRH，$M$ 的每个中间域 $E$ 的 $\zeta_E$ 在 $|\Im s|\le64$ 的所有非平凡零点都在中心线上。123 只计算三维因子的正高度零点。首个和末个正纵坐标落在

$$
 \frac{29430853}{2^{24}}<\gamma_1<\frac{29430854}{2^{24}},\qquad
 \frac{1071339820}{2^{24}}<\gamma_{123}<\frac{1071339821}{2^{24}} .
$$

逐一变号和独立的全矩形计数共同证明了完备性。[123 对整数端点](tetra3136-zeros64.csv) · [区间算术复核程序](verify_tetra3136.py)

**素数。** 对 $u,v>0$ 记 $g_v(x)=e^{-x^2/(2v)}/\sqrt{2\pi v}$，每个 $p\nmid14$ 的权为 $(\log p)g_v(\log p-u)$。令 $P_v(u)$ 是在 $M$ 中完全分裂的素数的加权和，$Q_v(u)$ 是 $p\equiv\pm1\bmod 7$ 的加权和。前者在后者中的相对密度为 $1/4$。对每个固定的 $\lambda\ge4/16381$，有

$$
 \lim_{u\to\infty}
 \frac{4P_{\lambda u}(u)-Q_{\lambda u}(u)}
      {\exp((1/2+\lambda/8)u)}=-\frac12 . \tag{1}
$$

因此在这种平滑下，对所有充分大的实数 $u$ 都有 $P_{\lambda u}(u)<Q_{\lambda u}(u)/4$。这个具体范围只用有限高度证书和经典素数定理。加入 Cureton–Thorner 的全高度无零半平面，还可取 $\lambda=2^{-12}$。下文给出使 (1) 成立的斜率的精确下确界 $\lambda_*$。上述四次域、六次域、十二次域的全高度 ERH 都等价于 $\lambda_*=0$。[Cureton–Thorner 的定理及 Lean 形式化](https://github.com/AxiomMath/DedekindZetaZeroFree)

## 从 Dirichlet 特征到一个三维因子

Alpöge 的原始发布与社区维护版的最终定理遍历每个模数的 Dirichlet 特征。对任意有限阿贝尔域 $F/\mathbb Q$，由 Kronecker–Weber 和包含分歧处的 Artin 分解，$\zeta_F$ 是一组本原 Dirichlet $L$ 函数的乘积，故得到全高度 ERH。对任意数域 $E$ 和它的正规闭包 $N$，$N/E$ 为伽罗瓦扩张，Aramata–Brauer 定理给出 $\zeta_N/\zeta_E$ 整。因而从所有伽罗瓦域下降到所有数域不另需零点估计。固定闭包的难处是由它的非阿贝尔表示引入的新因子。[原公告](https://alpo.ge/remark.pdf) · [核读版本的 Lean 题面](https://github.com/leanprover-community/riemann-hypothesis/blob/9c9f665ed79ce0fa8da61bd803b7a4bf1ad62ed0/src/Statements.lean) · [Sutherland 对原发布的独立核验](https://math.mit.edu/~drew/grh/) · [已有的分裂素数判据与下降证明](splitting-prime-criterion.md)

本例中 $f$ 在 2 处 Eisenstein，多项式判别式为 $56^2$。模 7 有 $f=(X+2)(X+1)^3$，且 $f(-1)=7$，故由 Dedekind 指标判据，7 不整除幂基指标。2 处的 Eisenstein 性同样排除 2，得到 $\mathcal O_K=\mathbb Z[X]/(f)$。又 $f(x)=x^2((x-1)^2+1)+2>0$，故 $K$ 的签名为 $(0,2)$。若 $\alpha_i$ 是四个根，三种配对的 $z=\alpha_1\alpha_2+\alpha_3\alpha_4$ 满足 $z^3-2z^2-8z+8=0$。写 $z=2r$，则

$$
 r^3-r^2-2r+1=0,\qquad k=\mathbb Q(r)=\mathbb Q(\zeta_7)^+ .
$$

例如可取 $r=-(\zeta_7+\zeta_7^{-1})$。该三次多项式不可约，结合平方判别式可知 $\operatorname{Gal}(M/\mathbb Q)=A_4$，且 $k=M^{V_4}$，其中 $V_4$ 为 Klein 四元群。[LMFDB 的四次域和整基](https://www.lmfdb.org/NumberField/4.0.3136.1) · [循环三次子域](https://www.lmfdb.org/NumberField/3.3.49.1)

记 $\tau$ 为四个顶点的置换表示去掉常数线所得的三维表示，$\chi$ 为其迹。$V_4$ 的三个非平凡特征被 $A_4/V_4$ 循环置换。任取其中一个 $\eta$，则 $\tau=\operatorname{Ind}_{V_4}^{A_4}\eta$。因此 $L(s)=L_k(s,\eta)$ 是非平凡二次 Hecke 特征的整函数。记 $S=M^{C_2}$，其中 $C_2$ 为任一二阶子群。$A_4$ 在顶点、六条边、三组对边以及群自身的作用给出

$$
 \zeta_K=\zeta L,\qquad
 \zeta_S=\zeta_k L,\qquad
 \zeta_M=\zeta_k L^3 . \tag{2}
$$

例如三组对边的置换迹在 $V_4$ 上为 3，在三循环上为 0。六条边的置换迹等于该迹加 $\chi$，正则表示的迹则等于该迹加 $3\chi$。这覆盖了 $A_4$ 的全部子群共轭类型。式 (2) 是表示恒等式，因而也给出分歧素数的正确局部因子。在开带内，每个因子全纯。$\zeta$ 和 $\zeta_k$ 的零点已有 Dirichlet 定理处理，余下的全高度问题就是 $L$ 的 RH。附件在四点仿射平面 $\mathbb F_2^2$ 上逐个验证这些特征标等式。[有限群恒等式的 Lean 证明](Tetra3136.lean) · [Artin 形式论与 Hecke 特征](https://www.jmilne.org/math/CourseNotes/CFT.pdf)

## 三次子域分裂时的三个二次剩余符号

对 $p\nmid14$ 写 $g=\operatorname{Frob}_p$。两类三循环在下表有相同的迹。$g\in V_4$ 等价于 $p\equiv\pm1\bmod7$。

| $g$ 的类型 | 元素数 | $\chi(g)$ | $\chi(g^2)$ |
|---|---:|---:|---:|
| 恒等元 | 1 | 3 | 3 |
| 两个不相交的对换 | 3 | $-1$ | 3 |
| 三循环 | 8 | 0 | 0 |

因此一次素数项的迹为 $4\mathbf1_{\{g=1\}}-\mathbf1_{p\equiv\pm1\bmod7}$。素数平方的迹恰等于 $3\mathbf1_{p\equiv\pm1\bmod7}$。后一个等式使 (1) 中的常数可以只靠有理等差数列的素数定理计算。

对前一个迹也有显式判别。在 $p\equiv\pm1\bmod7$ 时，取上述三次方程在 $\mathbb F_p$ 的三个根 $r_1,r_2,r_3$。则

$$
 \chi(g)=\sum_{j=1}^3\left(\frac{2r_j-1}{p}\right). \tag{3}
$$

事实上 $(\alpha_1+\alpha_2-1)^2=z-1=2r-1$。三种配对给出的三个 $2r_j-1$ 的乘积为 1。知道相应的带符号平方根，就能线性恢复四次方程的全部根。于是 $M/k$ 是由这些平方根给出的双二次扩张。模 $p$ 的三个符号乘积为 1，全正时 Frobenius 为单位元，否则有两个负号。例如 $13$ 和 $167$ 同余于 $-1\bmod7$。前者的三次根为 $3,5,6$，给出 $-,+,-$。后者的根为 $45,142,148$，给出 $+,+,+$。

原 Lean 证明的 `QH.MeanSq` 对任意有界权成立，但和式仍含辅助 Jacobi 扭转，并对辅助素数求均方。`QH.Ind` 用的是有理整数上的本原 Dirichlet 特征。这里的 $\eta$ 定义在 $k$ 的理想上。若要沿该路线得到全高度结论，一个充分也必要的目标是对每个 $\varepsilon>0$ 控制固定域的无扭转和

$$
 \sum_{p\le x,\ p\nmid14}\chi(\operatorname{Frob}_p)\log p
       =O_\varepsilon(x^{1/2+\varepsilon}). \tag{4}
$$

正向由 $L$ 的显式公式得出。反向对一次项分部求和，高次素数幂在 $\Re s>1/2$ 绝对收敛，故 $-L'/L$ 在该半平面全纯。函数方程再处理左半带。此前对判别式 $-23$ 的工作使用三阶理想类特征。式 (3) 把同一问题落实为一个实三次基域上的二次特征。[Lean 中的均方、诱导和窗口条件](https://github.com/leanprover-community/riemann-hypothesis/blob/9c9f665ed79ce0fa8da61bd803b7a4bf1ad62ed0/src/Development/AutomorphicArithmetic/MeanSquareBounds/PrimeSetTwoBounds.lean) · [判别式 $-23$ 的二维因子及其计数](cubic23-zeros.md)

## 可推广的数域族

判别式 $-23$ 的群 $S_3$ 和这里的 $A_4$ 分别是 $\operatorname{AGL}(1,3)$ 和 $\operatorname{AGL}(1,4)$。对这一整族，可以给出每个中间域中新增因子的精确重数。设 $q>2$ 为素数幂，$N/\mathbb Q$ 的伽罗瓦群是
$$
 G=\{x\mapsto ax+b:\ a\in\mathbb F_q^\times,\ b\in\mathbb F_q\}
   =V\rtimes A,\qquad k_0=N^V ,
$$
其中 $V$ 是平移群，$A=\mathbb F_q^\times$，记 $\pi:G\to A$。$k_0/\mathbb Q$ 是 $q-1$ 次循环扩张。

**仿射群的归约。** 取 $V$ 的任一非平凡复特征 $\eta$，经类域论看作 $k_0$ 上的 Hecke 特征，记 $L_N=L_{k_0}(s,\eta)$。对每个 $H\le G$ 和 $E=N^H$，
$$
 \zeta_E(s)=\zeta_{E\cap k_0}(s)\,L_N(s)^{m_H},\qquad
 m_H=\frac{q/|H\cap V|-1}{|\pi(H)|}. \tag{5}
$$
特别地，$E\not\subseteq k_0$ 时 $m_H>0$。在 Dirichlet GRH 下，这些 $E$ 的 ERH 各自等价于同一个 $L_N$ 的 RH。$q$ 为 2 的幂时，这唯一需要新增控制的函数是循环基域上的**二次** Hecke $L$ 函数。

**证明。** $A$ 在 $V$ 的非平凡加性特征上传递且自由，故 $\operatorname{Ind}_V^G\eta$ 不可约，次数为 $q-1$。它在 $V$ 上的迹是非平凡加性特征之和，在 $G\setminus V$ 上为零，正是 $q$ 点置换表示去掉常数线的迹 $\chi_q=q\mathbf1_{\{1\}}-\mathbf1_V$。$A$ 的 $q-1$ 个线性特征与它已用尽全部不可约表示，因为次数平方和为 $(q-1)+(q-1)^2=|G|$。在 $G/H$ 的置换表示中，由 Frobenius 互反，平凡于 $\pi(H)$ 的线性特征各出现一次，其余线性特征出现零次。这一部分对应 $N^{VH}=E\cap k_0$。余下维数 $[G:H]-[A:\pi(H)]$ 除以 $q-1$ 就是 (5) 的 $m_H$。Artin 诱导给出所述 Hecke 因子及包含分歧素数的恒等式。

若 $g(x)=ax+b$，则 $g^2(x)=a^2x+(1+a)b$。由此得到
$$
 \chi_q(g^2)=
 \begin{cases}
 (q-1)\mathbf1_{a=1}, & q\text{ 为偶数},\\
 \chi_q(g)+(q-1)\mathbf1_{a=-1}, & q\text{ 为奇数}.
 \end{cases} \tag{6}
$$
这说明两个例子的素数平方项为何不同。偶特征时，平方迹全部来自循环商 $A$。特别地，对任何这样的 $N$，若已核准 $L_N$ 在 $|\Im s|\le T$ 的非平凡零点均在线上，$T>\sqrt{3}/2$，令 $m_0=\operatorname{ord}_{1/2}L_N$。以相同高斯权记 $P_{N,v}$、$P_{k_0,v}$ 为分别在 $N$、$k_0$ 完全分裂的素数和，两者均略去 $N$ 的分歧素数。那么对每个固定的 $\lambda\ge (T^2-3/4)^{-1}$ 有
$$
 \frac{qP_{N,\lambda u}(u)-P_{k_0,\lambda u}(u)}
      {e^{(1/2+\lambda/8)u}}\ \longrightarrow\ -\frac12-m_0
 \qquad(q=2^d).
$$
证明见下节的显式公式。式 (6) 把平方项换成 $(q-1)$ 倍的循环子域分裂素数和，素数定理给出的密度为 $1/(q-1)$。中心零点贡献 $-m_0$，其余中心线零点经平滑后趋零，高度 $T$ 之外的估计只用 $\beta<1$。这里的有限高度前提已在本文给定的 $A_4$ 域核验。

## 从有限高度到所有充分大尺度

设 $\rho=\beta+i\gamma$ 遍历 $L$ 的非平凡零点，按重数计算。定义

$$
 \lambda_*=\sup\left(\{0\}\cup
 \left\{\frac{2\beta-1}{\gamma^2+1/4-\beta^2}\ :\
                      \beta>1/2\right\}\right). \tag{7}
$$

**定理。** 对正斜率 $\lambda$，式 (1) 成立当且仅当 $\lambda>\lambda_*$。若 $0<\lambda<\lambda_*$，(1) 左侧所取函数的下极限为 $-\infty$。若 $\lambda_*>0$，设 $m_*$ 为 (7) 取到最大值的零点总重数，计入共轭，则在 $\lambda=\lambda_*$ 的下极限为 $-1/2-m_*$。特别地，偏差在一列趋零正斜率上具有 (1) 的渐近式，就推出 $L$ 的全高度 RH。

**证明。** 置 $v=\lambda u$，$\sigma_\lambda=1/2+\lambda/8$。对 $-L'/L$ 的素数幂系数用 $g_v(\log n-u)$ 加权，记和为 $D_v(u)$。在 $\Re s=2$ 逐项 Mellin 反演，把积分线移到 $-1/2$，得到

$$
 D_v(u)=-\sum_\rho e^{u\rho+v\rho^2/2}
          +o_\lambda(e^{\sigma_\lambda u}). \tag{8}
$$

平移跨过的平凡零点在 0，其贡献为常数。函数方程与右半平面的绝对收敛给出左线上 $L'/L=O(\log(2+|t|))$，高斯乘子的模为 $e^{(-1/2+\lambda/8)u-vt^2/2}$。可沿避开零点的水平边取极限，固定导子的零点计数 $N(T)=O(T\log(2+T))$ 保证零点和绝对收敛。这也是前一判别式使用的高斯显式公式。[已有的高斯平滑与斜率判据](cubic23-gaussian-threshold.md) · [Booker 对此类完备 $L$ 函数的显式公式与零点计数](https://arxiv.org/abs/math/0507502)

令 $A(x)=3\sum_{p\le x,\,p\nmid14,\, p\equiv\pm1\bmod7}\log p=x+o(x)$。由上表，未分歧素数平方对 $D_v$ 的贡献是 $\int g_v(2\log x-u)\,dA(x)$。换成 $dx$ 后积分恰为 $\tfrac12e^{u/2+v/8}$（积分从 0 到无穷）。换元 $y=\log x$ 后，连续主项归一化得到均值 $u/2+v/4$、方差 $v/4$ 的正态分布。分部求和时，权导数的归一化绝对积分是该分布下 $4|y-u/2|/v$ 的期望，至多 $1+2/\sqrt v$。任何固定有界的 $x$ 区间贡献趋于零，其余部分可用 $|A(x)-x|\le\varepsilon x$。所以替换误差是 $o(e^{\sigma_\lambda u})$。

三次及以上幂在 $p^j\le x$ 的总绝对权 $\sum_{j\ge3,\ p^j\le x}3\log p\ll x^{1/3}\log^2(2x)$。同样分部求和并配方，高斯加权后至多为 $e^{u/3+v/18}$ 乘 $u$ 的一个固定次幂（$\lambda$ 固定）。分歧素数的对数幂是等间距数列，高斯权之和为 $O(1+v^{-1/2})$。从 (8) 减去这些项，得到

$$
 \frac{4P_v(u)-Q_v(u)}{e^{\sigma_\lambda u}}
  =-\frac12-\sum_\rho e^{(b_\lambda(\rho)-\sigma_\lambda)u}
              e^{i\gamma(1+\lambda\beta)u}+o_\lambda(1),\qquad
 b_\lambda(\rho)=\beta+\tfrac{\lambda}{2}(\beta^2-\gamma^2). \tag{9}
$$

高度证书排除了实零点，且 $\beta>1/2$ 的零点必有 $|\gamma|>64$，所以 (7) 分母为正。对 $\beta\le1/2$，总有 $b_\lambda(\rho)<\sigma_\lambda$。对 $\beta>1/2$，这个严格不等式恰等价于 $\lambda$ 大于 (7) 中该零点的比值。因比值在 $|\gamma|\to\infty$ 时趋零，正的上确界一定取到。若 $\lambda>\lambda_*$，在 (9) 中由高斯衰减逐项取极限，便得 (1)。

在临界斜率，只剩 $m_*$ 项的相位和。其和的实部至多为 $m_*$，又可在任意远的整数 $u$ 上由有限个频率的同时 Dirichlet 逼近使所有相位任意接近 1。故 (9) 的下极限恰为 $-1/2-m_*$。在临界斜率以下，取所有使 $b_\lambda$ 最大的零点，同样使它们的相位同时回到 1。其他项除以最大指数的指数函数后趋零，因 $b_\lambda\to-\infty$ 且零点局部有限。此时最大指数严格超过 $\sigma_\lambda$，得下极限 $-\infty$。整个论证容许重零点及频率间的有理关系。[判别式 $-23$ 的素数平方偏差](cubic23-gaussian-bias.md) · [Ng 关于 Artin 零点与素数竞赛的研究](https://www.cs.uleth.ca/~nathanng/RESEARCH/phd.thesis.pdf)

若所有零点满足 $\beta\le b$，其中 $1/2<b\le1$，则由 (7) 及 $|\gamma|>64$ 有 $\lambda_*<(2b-1)/(64^2+1/4-b^2)$。单用 $b=1$ 得 $4/16381$。Cureton–Thorner 适用于 $k$ 上的有限阶 Hecke 特征，取 $b=31/32$ 得 $960/4193599<2^{-12}$。由此得到开头所述的两个范围。将已经证明的正斜率继续压到零，仍需对 (3) 中的二次剩余符号取得新的算术估计。[全高度无零半平面的定理 1.1 和推论 1.2](https://github.com/AxiomMath/DedekindZetaZeroFree)

## 高度 64 的证书怎样复算

写 $\Gamma_{\mathbb R}(s)=\pi^{-s/2}\Gamma(s/2)$。由 (2)、导子判别式公式和签名，
$$
 G(s)=56^s\Gamma_{\mathbb R}(s)\Gamma_{\mathbb R}(s+1)^2,\qquad
 \Lambda(s)=G(s)L(s)=\Lambda(1-s).
$$
根数 $+1$ 也可从 $\zeta_K/\zeta$ 的函数方程直接读出。$\Lambda$ 是整函数，$W(z)=\Lambda(1/2+iz)$ 为偶函数，且在实轴上取实值。计算将 $W$ 作周期化，沿绝对收敛线取傅里叶系数，再回到中心线取符号。下面列出截断和边界计数所需的误差。[Booker 的傅里叶算法](https://arxiv.org/abs/math/0507502) · [Arb 的区间算术](https://fredrikj.net/arb/)

**欧拉系数。** $L$ 的局部因子为 $(1-X)\prod_{\mathfrak p\mid p}(1-X^{f_{\mathfrak p}})^{-1}$，$X=p^{-s}$。程序在每个 $p\le32768$ 上分解整基多项式，逐个使用不同不可约因子的次数。在 2 和 7 处分别得到局部因子 $1$ 和 $(1-X)^{-1}$。写 $L(s)=\sum a(n)n^{-s}$，由三个局部特征根的模至多为 1，有 $|a(n)|\le d_3(n)$。在 $\Re s=c=48$，截到 $H=32768$ 的绝对尾界为 $\zeta(2)^3 H^{2-c}<8.7\times10^{-208}$。程序先将这项误差乘入各点各自的 $G(c+it)$。

**傅里叶系数。** 取周期 $B=256$，$\omega=2\pi/B$，令 $b_m=B^{-1}\int_{\mathbb R}W(t)e^{-im\omega t}\,dt$。$\sum_jW(z+jB)$ 的傅里叶系数为 $b_m=b_{-m}\in\mathbb R$。对 $c>1$、$0<\alpha<\pi/2$ 记
$$
 C_c(\alpha)=\zeta(c)^3G(c)(\cos\alpha)^{-(3c+2)/2},\qquad
 U_c=\frac{2C_c(1)}{(3/2)B}.
$$
把三个 Gamma 积分分别旋转 $\pm\alpha$，可得 $|\Lambda(c+it)|\le C_c(\alpha)e^{-3\alpha|t|/2}$。函数方程给出 $1-c$ 上的同一界。对 $e^{\pm(3\alpha/2)z}W(z)$ 用条带 Phragmén–Lindelöf，该界在中间各线仍成立（完备 Hecke 函数为有限阶）。傅里叶积分移到 $z=t-i(c-1/2)$，于是 $|b_m|\le U_c e^{-(c-1/2)|m|\omega}$。

程序取 $J=2048$，$h=B/J=1/8$，沿 $\Re s=48$ 在 $|jh|\le R=384$ 取样，以离散傅里叶变换求 $0\le m\le D=480$ 的系数。记 $a=c-1/2$，$\kappa=3/2$。在乘 $e^{-am\omega}$ 之前，漏掉的采样点和折入第 $m$ 项的其他频率分别有界于

$$
 \frac{2C_c(1)e^{-\kappa(R+h)}}{J(1-e^{-\kappa h})}
 <1.9\times10^{-117},\qquad
 \frac{2U_{c+22}e^{-22\omega(J-D)}}{1-e^{-22\omega J}}
 <2.0\times10^{-155}.
$$

第二式来自沿右线周期化后的系数 $b_\ell e^{a\ell\omega}$，将 $\ell=m+jJ,\ j\ne0$ 的项按等比级数求和。在 $|\Re z|\le64,\ |\Im z|\le3/2$，舍去 $|m|>D$ 的总误差与撤销周期化的误差分别为至多

$$
 \frac{2U_c e^{-(a-3/2)\omega(D+1)}}{1-e^{-(a-3/2)\omega}}
 <1.25\times10^{-102},\qquad
 \frac{2C_2(1)e^{-\kappa(B-64)}}{1-e^{-\kappa B}}
 <2.22\times10^{-122}.
$$

**全矩形计数。** 取顶半边 $s=\sigma+64i$，$1/2\le\sigma\le2$，以 $h_0=1/64$ 分段。每段用连续 21 个节点作 Newton 插值，从一端向前或向后选点。在节点的单位圆盘上可取模上界 $V=C_3(3/2)e^{-(9/4)(64-1)}$。Cauchy 不等式与插值余项给出每个实分量的误差 $\le V\,20!\,h_0^{21}$。因 $|\binom{x}{j}|\le1/j$（$0\le x\le1,\ j\ge1$），各阶有限差分除以 $j$ 即可控制该分量相对端点的变化。在每一段至少有一个坐标与零分离，最小余量大于 $2.68\times10^{-64}$。因此该段的辐角变化是两端值之比的主辐角。记这一顶半边从左到右的辐角增量为 $\Delta$。

在 $\Re s=2$ 从高度 0 起连续取辐角，因为 $|\Im\log L|\le3\log\zeta(2)<\pi/2$，右上角的 $L$ 可以取主辐角。$G$ 的连续辐角则从三个对数 Gamma 因子计算。由实对称和函数方程，矩形 $[-1,2]\times[-64,64]$ 内按重数计的零点数之半等于
$$
 \frac{\arg_{\rm cont}\Lambda(2+64i)-\Delta}{\pi}.
$$
区间计算给出此数与 123 的距离小于 $2\times10^{-26}$。同时，CSV 中 123 个互不相交的正高度区间都在中心线上有严格变号，负高度各有一个共轭区间。它们用尽全部 246 个名额，因而每个区间内恰有一个单零点，也排除了实零点与其他偏线零点。

在本文目录运行下列命令。第一条核验表内所有端点和矩形绕数，第二条从 $1/64$ 网格重新隔离零点并写出 CSV，第三条核验两条有限群定理。数值证明使用 Python 3 与 python-flint 0.9.0 的 896 位球算术。Lean 文件只需标准库。上述解析推导为纸面证明。[数值复核程序](verify_tetra3136.py) · [区间端点](tetra3136-zeros64.csv) · [Lean 源文件](Tetra3136.lean)

```sh
python verify_tetra3136.py
python verify_tetra3136.py --regenerate
lean +leanprover/lean4:v4.35.0 Tetra3136.lean
```
