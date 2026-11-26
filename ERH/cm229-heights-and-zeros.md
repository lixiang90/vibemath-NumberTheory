---
title: "判别式 229 的全实三次域及其 CM 扩张的两类高度与 192 个零点"
date: "2026-11-26"
abstract: "对非伽罗瓦全实三次域 F，分离 F(i) 的诱导型与混合型 Faltings 高度。对判别式 229 给出两类高度差的严格区间、三个无穷位的单位平衡常数，以及两个非阿贝尔因子在高度 64 内全部正高度零点的证书。"
---

# 判别式 229 的全实三次域及其 CM 扩张的两类高度与 192 个零点

设 $\alpha^3-4\alpha+1=0$，$F=\mathbb Q(\alpha)$，$E=F(i)$。从 $\mathbb Q$ 上的 Dirichlet GRH 推广到这对数域，要分别控制一个偶的和一个奇的二维 Artin 表示。我们证明了三个可分别复核的结论。

- $E$ 的八个 CM 型有两种稳定 Faltings 高度。混合型比诱导型高 $1.1323555935686836618893\ldots$，并可把差值写成两个 $L$ 函数在 1 处的对数导数。
- 只用 $F$ 中两个显式单位，就能把任意 $x\in E^\times$ 在三个无穷位的对数绝对值都调到距其平均值小于 $1.297$。同一对单位也适用于 $F^\times$。
- 在 $0<\Im s\le64$ 的开临界带，两个非阿贝尔因子分别有 82 和 110 个零点，全部是临界线上的单零点。两个因子在这段高度没有共同零点，也都没有实的非平凡零点。[192 个隔离区间](cm229-zeros64.csv) · [区间算术核验程序](verify_cm229.py)

## 一个偶因子和一个奇因子

多项式 $f(X)=X^3-4X+1$ 没有有理根，判别式 $229$ 是素数。故 $\mathcal O_F=\mathbb Z[\alpha]$。三个实根分别在 $(-3,-2),(0,1),(1,2)$，正规闭包 $N$ 的群为 $S_3$，其二次子域为 $k=\mathbb Q(\sqrt{229})$。$N$ 全实，故 $\operatorname{Gal}(N(i)/\mathbb Q)=S_3\times C_2$。$F$ 与 $\mathbb Q(i)$ 的判别式互素，乘积整基给出 $|D_E|=4^3\cdot229^2=3\,356\,224$。[LMFDB 的三次域记录，取其生成元的相反数](https://www.lmfdb.org/NumberField/3.3.229.1)

写 $\chi=\chi_{-4}$，$D(s)=L(s,\chi)$。令 $\rho$ 是 $S_3$ 的二维标准表示，$L_+=L(s,\rho)$，$L_-=L(s,\rho\otimes\chi)$。三点置换表示为 $\mathbf1\oplus\rho$，复共轭交换每一对延拓到 $E$ 的嵌入。Artin 形式性给出含所有分歧素数在内的恒等式

$$
 \zeta_F=\zeta L_+,\qquad \zeta_E=\zeta D L_+L_-,\qquad
 L_F(s,\eta_{E/F})=D L_- . \tag{1}
$$

$\rho$ 是由 $\operatorname{Gal}(N/k)=C_3$ 的非平凡三阶特征诱导而来。扭转后仍由 $k$ 上非平凡有限阶 Hecke 特征诱导，故 $L_+,L_-$ 都是整函数。$N$ 和 $N(i)$ 的正则表示又给出

$$
 \zeta_N=\zeta_k L_+^2,\qquad
 \zeta_{N(i)}=\zeta_{k(i)} L_+^2L_-^2. \tag{2}
$$

一般地，在一个正规 CM 域 $M/\mathbb Q$ 中，复共轭 $c$ 位于群的中心。若 $K=M^H$ 是 CM 子域，则 $K^+=M^{\langle H,c\rangle}$。在 $K$ 的嵌入置换表示中，把每一对共轭嵌入的基向量相加、相减，分别得到 $+1$ 和 $-1$ 特征空间。和空间就是 $K^+$ 的嵌入置换表示，差空间给出相对因子 $\zeta_K/\zeta_{K^+}$。因此在这个群的表示环内，任意多个相对 CM 因子的虚特征线性组合都只能落在奇部分。式 (1) 的 $L_+$ 位于偶部分。这说明即使已控制相对因子的全部零点，仍要为全实基域的因子提供估计。

## 八个 CM 型的两种高度

这一段允许 $F$ 是**任意非伽罗瓦全实三次域**，仍取 $E=F(i)$，以同样的置换表示定义 $L_-$。每个 CM 型在 $F$ 的三个实嵌入上各选 $i$ 或 $-i$。$S_3$ 排列三个选择，$C_2$ 同时换号。因此八个选择恰分为两个 Galois 轨道。两个同号型来自 $\mathbb Q(i)$，记其高度为 $h_I$。其余六个混合型的高度记为 $h_M$。这里用 Andreatta–Goren–Howard–Madapusi Pera 的稳定高度归一化，复乘环取最大整环。高度对共轭不变，故 $\bar h_E=(h_I+3h_M)/4$。[平均 Colmez 定理 A 与归一化](https://www.math.mcgill.ca/goren/PAPERSpublic/Colmez_v3.5.pdf) · [Yang–Yin 关于酉型 $(n-1,1)$ 的逐型 Colmez 定理](https://arxiv.org/abs/1711.02988)

设 $h_i$ 为 $\mathbb C/\mathbb Z[i]$ 的稳定高度。对诱导型，$i$ 在切空间上以同一标量作用。周期格作为 $\mathbb Z[i]$ 模无挠，因 $\mathbb Z[i]$ 是主理想整环而自由，秩为 3。它的三个基向量在 $\mathbb C$ 上线性无关，因为它们和其 $i$ 倍数实张成整个 $\mathbb C^3$。因此底层复阿贝尔簇同构于 $(\mathbb C/\mathbb Z[i])^3$。高度的乘积可加性给出 $h_I=3h_i$。

写 $\Delta_F=|D_E|/(4^3 D_F)$。对 $E/F$ 以及 $\mathbb Q(i)/\mathbb Q$ 分别应用平均 Colmez 公式的 $s=1$ 形式。将后者乘 3 后相减，次数所决定的常数恰好消去。式 (1) 给出

$$
 3(h_M-h_I)=2\,\frac{L_-'}{L_-}(1)-4\,\frac{D'}D(1)+\log\Delta_F . \tag{3}
$$

例如 $D_F=229$ 时 $\Delta_F=229$。利用 $D(s)=4^{-s}(\zeta(s,1/4)-\zeta(s,3/4))$、Hurwitz zeta 在 $s=0$ 处的导数公式和函数方程，有

$$
 \frac{D'}D(1)=\gamma+\log(2\pi)
          -2\log\frac{\Gamma(1/4)}{\Gamma(3/4)} .
$$

从下文的完备奇因子得到 $-0.5271084418694659866091<(L_-'/L_-)(1)<-0.5271084418694659866089$。代入 (3) 的严格区间为

$$
 1.1323555935686836618893<h_M-h_I<1.1323555935686836618894 . \tag{4}
$$

在上述归一化中，$h_I=-4.971319548574\ldots$，$h_M=-3.838963955005\ldots$。从而混合型高度高于全型平均 $(h_M-h_I)/4$，诱导型高度低于全型平均 $3(h_M-h_I)/4$。原 Lean 周期命题固定了曲线理想、模型和指标 $\tau$，平均的是同一几何对象的共轭度量。若要在高维有同样的周期等式，必须先确定该对象选取的 CM 型及其权重。式 (4) 给出了把两种平均混用所带来的误差。它通过解析公式计算高度，对个别嵌入的周期仍需几何估计。[Hurwitz zeta 导数公式](https://dlmf.nist.gov/25.11#E18) · [平均高度的 $s=1$ 公式](real-cm-prime-windows.md) · [原证明对固定对象的周期平均](https://github.com/leanprover-community/riemann-hypothesis/blob/9c9f665ed79ce0fa8da61bd803b7a4bf1ad62ed0/src/Development/AbelianVarietyTheory/Proofs/FoundationProofs/MeanLogPeriodConjugatesEqChowlaSelbergStatement.lean#L45)

这也给出对提升链选型的一个数值条件。给各步赋非负权，令 $W_I,W_M$ 分别为两类的总权，$W=W_I+W_M$。逐型高度之和与采用平均高度所得主项之间的差恰为

$$
 W_I h_I+W_M h_M-W\bar h_E
       =(W/4-W_I)(h_M-h_I).
$$

由 (4)，在这个六次域中，两个高度主项相等当且仅当诱导型占总权的四分之一。等权的四步中取一个诱导型和三个混合型，就能在高度之和上精确消去差值。要将这种选型带入实际提升，还须证明这些步可满足同一组有限处整除条件。

## 三个无穷位的单位平衡

回到判别式 229，记三个递增的实根为 $\alpha_j$，选定 $\sigma_j(\alpha)=\alpha_j$。在 $E$ 中任选每对共轭嵌入中的一个延拓，绝对值与此选择无关。对 $x\ne0$ 令 $r_j=\log|\sigma_j x|$ 及 $\bar r=\frac13\sum_j r_j$。这分别是 $F^\times$ 上的 $\frac13\log|N_{F/\mathbb Q}x|$ 和 $E^\times$ 上的 $\frac16\log|N_{E/\mathbb Q}x|$。

因为 $f(0)=f(2)=1$，$\alpha$ 和 $\alpha-2$ 都是范数 $-1$ 的单位。令 $u_j=\log|\alpha_j|$，$v_j=\log|\alpha_j-2|$。两列坐标和均为零。又 $u_1,v_1,v_2>0$，$u_2<0$，故 $u_1v_2-u_2v_1>0$。于是 $u,v$ 张成坐标和为零的平面。将 $\bar r\mathbf1-r$ 在这两列下的系数就近取整，得到一个单位 $\epsilon=\alpha^m(\alpha-2)^n$，满足

$$
 |\log|\sigma_j(\epsilon x)|-\bar r|\le B_j,\qquad
 B_j=\tfrac12\bigl(|u_j|+|v_j|\bigr),\qquad \max_j B_j<1.297 . \tag{5}
$$

三个 $B_j$ 约为 $1.081814,0.963645,1.296448$，程序用实代数根与区间对数核验严格上界。若要求单位全正，在两个单位的平方张成的格上取整，得到 $\max_j 2B_j<2.594$。这两种操作均保持所有有限处的赋值。

对任意实权 $w_j$，写 $\bar w=\frac13\sum_j w_j$。因三个偏差之和为零，(5) 推出

$$
 \left|\sum_j w_j\log|\sigma_j(\epsilon x)|-\bar r\sum_j w_j\right|
 \le \sum_j |w_j-\bar w|B_j . \tag{6}
$$

所有权相等时，误差精确为零。更一般地，若 $m$ 步中每步满足 $\max_j|w_j-\bar w|\le H$，则这 $m$ 个标量平衡误差的绝对值之和小于 $3.342\,H m$（$H m>0$），因为 $\sum B_j<3.342$。这一部分在权重偏差有界时只需线性的误差预算。一般单位格上的取整和中心化加权不等式各有一条 Lean 证明。提升链中随秩变化的 Hermitian Gram 矩阵还需要单独估计。[Lean 单位格引理](CMUnitBalance.lean)

## 高度 64 的完整零点计数

在本例中，写 $\Gamma_{\mathbb R}(s)=\pi^{-s/2}\Gamma(s/2)$。完备因子为 $G(s)=q^{s/2}\Gamma_{\mathbb R}(s+e)^2$，两种情形的参数和计算结果如下。$N^+(64)$ 按重数计开临界带内全部正高度不超过 64 的零点。

| 因子 | 导子 $q$ | 实位移位 $e$ | $N^+(64)$ | 顶边每段的最小坐标余量 |
| --- | ---: | ---: | ---: | ---: |
| $L_+$ | $229$ | $0$ | $82$ | $>7.51\cdot10^{-44}$ |
| $L_-$ | $3664$ | $1$ | $110$ | $>9.45\cdot10^{-43}$ |

实位因子由复共轭在两个表示上的标量作用给出。导子可由 (1) 的判别式读出，$\zeta_F/\zeta$ 给出 $q_+=D_F$，$\zeta_E/(\zeta_F D)$ 给出 $q_-=|D_E|/(4D_F)=16\cdot229=3664$。两个完备函数 $\Lambda=G L$ 均有根数 $+1$，由 Dedekind zeta 及 $D$ 的函数方程即得。

**局部输入。** 对 $p\nmid 2\cdot229$，按照 $f\bmod p$ 的分解型 $(1,1,1),(1,2),(3)$，$L_+$ 的局部因子分母依次为 $(1-T)^2,1-T^2,1+T+T^2$，$T=p^{-s}$。在 $229$ 处有两个不同的一次因子，故局部因子为 $(1-T)^{-1}$。在 2 处 $f$ 的分解型为 $(1,2)$，奇扭转在此处的局部因子为 1。因此对所有正整数，$a_-(n)=\chi(n)a_+(n)$，且 $|a_\pm(n)|\le d_2(n)$。

**截断与计数。** 令 $W(z)=\Lambda(1/2+iz)$。它在实轴上取实值，且为偶整函数。对 $c>1$ 和 $0<\theta<\pi/2$ 置 $C_c(\theta)=\zeta(c)^2G(c)(\cos\theta)^{-(c+e)}$。两个 Gamma 积分各旋转 $\theta$ 后，有 $|\Lambda(c+it)|\le C_c(\theta)e^{-\theta|t|}$。函数方程与条带 Phragmén–Lindelöf 原理将此界延至 $1-c\le\Re s\le c$。程序沿 $c=48$ 截取 $n\le32768$ 的 Dirichlet 级数，右线上的级数尾界为 $\zeta(2)^2\,32768^{-46}$。取周期 256、步长 $1/8$、采样半径 512，以 2048 点离散傅里叶变换求前 480 个正频率及常数项，用 $c+22$ 界定混叠。以共轭和偶性恢复其余频率。误差按前篇的等比级数逐项向外取整，其中 Gamma 移位取本篇表中数值。[傅里叶截断与插值误差的推导](tetra3136-prime-bias.md) · [Booker 的傅里叶算法](https://arxiv.org/abs/math/0507502) · [共用的球算术实现](verify_tetra3136.py)

在 $s=\sigma+64i,\ 1/2\le\sigma\le2$ 上取步长 $1/64$。对每段以 21 个连续节点作 Newton 插值，用 $C_3(3/2)e^{-(3/2)(64-1)}20!/64^{21}$ 控制每个实分量的余项。上表的正余量保证每段至少有一个坐标始终与零分离，故可累加端点比值的主辐角。右边界上 $|\Im\log L(2+it)|\le 2\log\zeta(2)<\pi/2$，对数 Gamma 给出 $G$ 的连续辐角。利用共轭与函数方程还原矩形 $[-1,2]\times[-64,64]$ 的绕数，其一半分别落在以 82 和 110 为中心、半径 $10^{-25}$ 的区间内。

CSV 各行的两个端点分母均为 $2^{24}$。每行在临界线上严格变号，再计入负高度的共轭区间，已经用尽两个矩形的总零点数。这同时证明单性、完备性和实高度零点的排除。两个不同因子的相邻区间端点之间最少相距 $22898/2^{24}$。由 (1)、(2) 以及已核验的 Dirichlet GRH，$F,E,N,N(i)$ 的 Dedekind zeta 在 $|\Im s|\le64$ 满足 ERH。式 (2) 的平方保留了正规闭包中相应零点的重数。[逐点隔离表](cm229-zeros64.csv) · [Dirichlet GRH 的独立内核核验](https://math.mit.edu/~drew/grh/)

同一个傅里叶多项式还可计算 (3) 中的导数。取 $h=2^{-64}$，用 $\Lambda(1\pm h)$ 的中心差商近似 $\Lambda'(1)$。以 1 为心的单位圆盘上有 $|\Lambda|\le C_2(1)$。Cauchy 系数界给出差商的绝对误差至多 $C_2(1)h^2/(1-h^2)$。除以严格与零分离的 $\Lambda(1)$，再减去 $\tfrac12\log3664-\log\pi+\psi(1)$ 即得 $L_-'/L_-(1)$ 的上述区间。

## 全高度问题的一个具体等价式

对 $p\nmid 2\cdot229$ 写 $a_p=\#\{t\bmod p:f(t)=0\}-1$。这是 $\rho$ 的 Frobenius 迹，奇扭转的迹为 $\chi(p)a_p$。已有 Dirichlet GRH 后，$E$ 的全高度 ERH 等价于对 $b=1,3$ 同时有

$$
 \sum_{\substack{p\le x,\ p\equiv b\bmod 4\\p\nmid 2\cdot229}} a_p\log p
       = O_\varepsilon(x^{1/2+\varepsilon})\qquad(\varepsilon>0). \tag{7}
$$

事实上两和的和与差就是两个因子的一次素数和。两因子若满足 RH，标准截断显式公式给出 $O(\sqrt{x}\log^2(2x))$，去掉素数高次幂只花 $O(\sqrt{x}\log^2(2x))$。反过来，(7) 经分部求和使两个对数导数的一次素数级数在 $\Re s>1/2$ 局部一致收敛，高次幂级数在那里绝对收敛。与 $\Re s>1$ 的欧拉乘积接合，便排除右半临界带的零点，再用函数方程。式 (7) 因而指出这对**固定数域**还缺哪两项抵消。伴篇的大筛定理控制随导子移动的辅助二次特征，在同一窗口内同时给出全实与 CM 两侧的估计。[二次 Hecke 素数窗口与平均高度](real-cm-prime-windows.md)

将两份 Python 文件和 CSV 放在同一目录，安装 python-flint 0.9.0 后执行前两条命令，分别复核已发布证书和从欧拉因子重建区间表。第三条在 Lean 4.35.0、mathlib 提交 56a005855ea17e0b0958f73ec801ae140b0a1679 的 Lake 环境中执行，参数指向本篇的 Lean 文件。数值证明使用 896 位 Arb 球算术，式 (1) 至 (7) 的推导见正文。[核验程序](verify_cm229.py) · [共享傅里叶程序](verify_tetra3136.py) · [区间 CSV](cm229-zeros64.csv) · [单位格的 Lean 证明](CMUnitBalance.lean)

~~~sh
python verify_cm229.py
python verify_cm229.py --regenerate
lake env lean CMUnitBalance.lean
~~~
