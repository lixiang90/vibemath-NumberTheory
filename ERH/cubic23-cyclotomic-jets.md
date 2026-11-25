---
title: "圆分 Hasse 导数与判别式 −23 的趋零斜率窗口"
date: "2026-11-25"
abstract: "用相邻两阶单位根上的 Hasse 导数精确恢复整数多项式的系数公因子，包括小剩余特征。在高斯窗口中按权的二进位分层，取样点数与导数阶数均可压至 Y 的一个趋近 1/2 的幂。给出未约分取样高度的双向估计，并将权 7 的首个非零导数与目标素数和之间的额外高度压到平方根量级。"
---

# 圆分 Hasse 导数与判别式 −23 的趋零斜率窗口

设 $k=\mathbb Q(\sqrt{-23})$，$F$ 为三次多项式 $X^3-X-1$ 完全分裂的素数集，$C$ 为其不可约的未分歧素数集。对偶数 $d\ge2$，前篇用两个 theta 权构造整数 $U_p=A_d(p)^3,V_p=A_{3d}(p)$。行 $(U_p,V_p)$ 的公因子在 $C$ 上是 $p^{3d}$，在 $F$ 上是 3，在换位处是 1。其最大绝对值对所有好素数都在 $p^{3d}/8$ 与 $27p^{3d}$ 之间。本篇考虑将这一局部结论送入 Stepanov 重数论证，需要多大的取样以及多精确的高度界。[局部整数与公因子定理](cubic23-stepanov-content.md)

## 两个圆分网格同时处理所有素数

设 $W(T)=\sum_{i=0}^D c_iT^i\in\mathbb Z[T]$ 非零，$D\ge1$，$s=\lceil\sqrt{D+1}\rceil$，$h=\lfloor D/s\rfloor$。第 $j$ 个 Hasse 导数 $W^{[j]}(T)=[Z^j]W(T+Z)$ 可在任意剩余特征使用。取以下整数多项式的全部系数，称其最大绝对值为 $J(W)$，

$$
 W^{[j]}(T)\bmod (T^q-1),\qquad
 q\in\{s,s+1\},\quad 0\le j\le h . \tag{1}
$$

取余使用首一多项式，在 $\mathbb Z[T]$ 内完成。每个网格最多 $s+1$ 个点，导数阶数最多 $h$，总共有至多 $(2s+1)(h+1)$ 个整数坐标。

**定理 1。** (1) 的全部整数坐标的正最大公因子恰为 $\operatorname{cont}W=\gcd_i|c_i|$，在每个素数处的赋值都相同。

**证明。** 先除掉 $\operatorname{cont}W$。对任一素数 $\ell$，约化得到 $\mathbb F_\ell[T]$ 上非零多项式。在 $s,s+1$ 中选与 $\ell$ 互素的一项 $q$，则 $T^q-1$ 在代数闭包有 $q$ 个不同的根。如果 (1) 的坐标全为零，每个根处的头 $h+1$ 个 Hasse 导数都为零，由 Taylor 展开其重数至少 $h+1$。但 $q(h+1)\ge s(\lfloor D/s\rfloor+1)>D$。故每个素数都有某个非零模 $\ell$ 的坐标，乘回内容即得。[有限点集上的 Hasse 重数界的 Lean 证明](Theta23Local.lean)

只取连续整数点会遇到小特征重合。例如 $W=(T^2-T)^2$ 本原，在模 2 下，整数点 $0,1,2$ 处的前两个 Hasse 导数（阶 $0,1$）全为零。相邻两阶的构造使这种损失在每个素数处都可修复。

还可在有理数上直接看见分母。固定一个网格阶 $q$ 和 $0\le r<q$、$r\le D$，令 $n=\lfloor(D-r)/q\rfloor$。在第 $j$ 个余式中，指数同余于 $r-j\bmod q$ 的坐标为
$$
 H_{r,j}=\sum_{k=0}^n {r+kq\choose j}c_{r+kq}. \tag{2}
$$
取 $0\le j\le n$，系数矩阵行列式的绝对值为 $q^{n(n+1)/2}$。这是在等差点列上取二项式基底的 Vandermonde 行列式，其首项系数分别为 $1/j!$。所以仅用一个网格时，可能损失的素因子都整除 $q$。

## 取样数减少后，阿基米德高度付多少

**定理 2。** 若 $W=\prod_p(U_p+TV_p)^{m_p}$、$D=\sum m_p\ge1$，则有绝对常数一致的双向估计
$$
 \log J(W)=3d\sum_p m_p\log p+O(D). \tag{3}
$$

为证下界，取 $q=s$ 并沿用 (2)。令 $L_k(x)$ 为节点 $r,r+q,\ldots,r+nq$ 上第 $k$ 个 Lagrange 基多项式。把它写成 $\sum_{j=0}^n\Delta^j L_k(0){x\choose j}$，其中 $\Delta$ 为向前差分，则
$$
 c_{r+kq}=\sum_{j=0}^n\Delta^jL_k(0)\,H_{r,j}.
$$
对 $0\le b\le n<q$，有 $|(b-r)/q|\le1$。由基多项式的乘积公式，
$$
 |L_k(b)|\le\frac{\prod_{0\le l\le n,\ l\ne k}(l+1)}{k!(n-k)!}
 ={n+1\choose k+1}\le2^{n+1}.
$$
故 $|\Delta^jL_k(0)|\le2^{j+n+1}$，$\max_i|c_i|\le2^{2h+2}J(W)$。Jensen 公式给 $\mathcal M(W)\le\sum_i|c_i|$，其中 $\mathcal M$ 是 Mahler 测度。反过来，由一次因子的乘积，$\sum_i|c_i|\le2^D\mathcal M(W)$。式 (2) 的系数至多 $2^D$，于是
$$
 \frac{\mathcal M(W)}{(D+1)2^{2h+2}}\le J(W)\le4^D\mathcal M(W).
$$
代入前篇的 $\log\mathcal M(W)=3d\sum m_p\log p+O(D)$ 即得 (3)。所以改在单位根处取样可以保护所有有限处的重数，直接估计这些导数的最大绝对值仍会付出每一个被取入的素数的高度。

## 对高斯权作二进分层

取 $u=\log Y$，$v=\lambda u$，$0<\lambda\le1/16$，$1\le v$，记
$$
 k_v(t)=\frac{e^{-t^2/(2v)}}{\sqrt{2\pi v}},\quad
 B_v(u)=\sum_{p\in C}(\log p)k_v(u-\log p),\quad
 P_v(u)=\sum_{p\in F}(\log p)k_v(u-\log p).
$$
$\chi(p)=1$ 的素数权和 $Q_v=B_v+P_v$ 由 Dirichlet GRH 满足 $Q_v=e^{u+v/2}/2+O_\lambda(e^{u/2+v/8})$。前稿给出了此等式所用的素数幂估计。[显式公式和阈值](cubic23-gaussian-threshold.md) · [类特征均方的对角项 $3B_v$](cubic23-class-character-mean-square.md)

置
$$
 a_\lambda=\lambda+\sqrt{\lambda+\lambda^2},\quad X=Y^{1+a_\lambda},
 \quad H=\lceil X/\sqrt Y\rceil,\quad
 m_p=\lfloor H k_v(u-\log p)\rfloor .
$$
只取 $p\le X,\ \chi(p)=1,\ p\ne3,23$。将 $m_p=\sum_b2^b e_{p,b}$ 作二进展开，各位 $e_{p,b}$ 为 0 或 1。对每个非空位集取
$$
 W_b(T)=\prod_{p:\,e_{p,b}=1}(U_p+TV_p),\qquad D_b=\deg W_b.
$$
令 $g_b$ 为从 (1) 得到的正公因子，$g_{b,3'}=g_b/3^{v_3(g_b)}$。则对上述 $Y,\lambda,v$ 一致成立
$$
 \boxed{\quad
 \frac1{3dH}\sum_b2^b\log g_{b,3'}=B_v(u)+O(\sqrt Y).\quad} \tag{4}
$$

**证明。** 由定理 1 与局部公式，左侧恰等于 $\sum_{p\in C,\ p\le X,\ p\ne3}(m_p/H)\log p$。取整误差至多 $\sum_{p\le X}\log p/H\ll\sqrt Y$。截尾取 $t=a_\lambda/\lambda>5$，用 $k_v(w)\le (2\pi v)^{-1/2}\exp(tw+t^2v/2)$ 和 $\sum_{p>X}\log(p)p^{-t}\ll X^{1-t}$。后式由 Chebyshev 界分部求和得到，常数对 $t\ge5$ 一致。指数为
$$
 t+t^2\lambda/2+(1+a_\lambda)(1-t)
 =1+a_\lambda-a_\lambda^2/(2\lambda)=1/2 .
$$
删去 3 至多影响 $O(1)$。全程只用初等尾界，(4) 本身不需要 Dirichlet GRH。

因为 $k_v\le1$，只有 $O(\log Y)$ 层。每层 $D_b\le\pi(X)\ll X/\log X$，所以在每个网格中，**取样点数和所需最高导数阶数**皆为
$$
 O\!\left(\frac{Y^{(1+a_\lambda)/2}}{\sqrt{\log Y}}\right).
$$
当 $\lambda\downarrow0$，这个幂指数趋于 $1/2$。若把全部整数权直接当成单个多项式的次数，次数会额外受 $H$ 放大。二进分层消去了这一损失。这里计量的是两个方向各自的规模，每层的坐标总数是 $O(D_b+1)$。Alpöge 原说明中约 $\sqrt X$ 的量是提升链长度；此处得到的是可严格实现的导数阶数和取样点数。[原说明](https://alpo.ge/remark.pdf)

若希望继续只用辅助虚二次域的惰性位置，取 $\mathcal Q=\{Y^5<q\le Y^{10}:q\equiv3\bmod8,\ q\text{ 为素数}\}$。在 $W_b$ 中只保留奇素数 $p$ 满足 $(p/q)=-1$，所得去 3 内容记作 $g_{b,q,3'}$。二次互反律将这个条件解释为 $p$ 在 $\mathbb Q(\sqrt{-q})$ 惰性。对每个奇 $p\le X$，模 8 特征乘 $q\mapsto(p/q)$ 仍在 $p$ 处有非平凡分歧。和前稿同一条 Dirichlet GRH 素数和估计给出一致概率 $1/2+O(Y^{-5}(1+u)^3)$。因此
$$
 \frac{2}{3dH}\frac1{|\mathcal Q|}\sum_{q\in\mathcal Q}\sum_b2^b\log g_{b,q,3'}
 =B_v(u)+O(\sqrt Y). \tag{5}
$$
这里舍去 2 只影响 $O(1)$，概率误差乘上的总权为 $O(X)$。式 (5) 只把这个内容抽样到原方法所用类型的惰性位置。它尚未把 $W_b$ 实现成原酉群提升链上的整截面。[辅助导子均方的证明](cubic23-class-character-mean-square.md)

## 要得到零斜率还缺哪条高度不等式

固定 $d=6$ 就可使用 (4)，对应的模形式权是 7 和 19。写 $\sigma_\lambda=1/2+\lambda/8$。在已核验的 Dirichlet GRH 下，以下是一个充分且必要的条件：存在 $\lambda_j\downarrow0$，使对每个 $\lambda=\lambda_j$ 和每个 $\varepsilon>0$，有任意大的整数 $N$，在 $N\le u\le N+\ell(N)$ 内每个整数 $u$ 都有
$$
 \frac1{3dH}\sum_b2^b\log g_{b,3'}
 \le \frac{e^{u+v/2}}3+O_{\lambda,\varepsilon}(e^{(\sigma_\lambda+\varepsilon)u}). \tag{6}
$$
其中 $\ell(N)$ 可预先选为任意趋于无穷的整值函数。由 (4) 和 $Q_v$ 的估计，(6) 等价于前稿的完全分裂素数单侧下界，常数可由略缩小 $\varepsilon$ 吸收。前稿的斜率定理再给 $\lambda_{23}=0$ 等价于 (6)，亦等价于 $\zeta_K$ 的 ERH。[区段的量词及两方向证明](cubic23-gaussian-threshold.md) · [已获证书的正斜率范围](cubic23-height-8192.md)

可给出一次直接的高度比较。在 $\chi(p)=1$ 的行上，把 (6) 的 $g_{b,3'}$ 换成未约分取样的 $J(W_b)$，由 (3) 得归一化结果
$$
 \frac1{3dH}\sum_b2^b\log J(W_b)
 =Q_v(u)+O\!\left(\sqrt Y+\frac{e^{u+v/2}}{du}\right). \tag{7}
$$
为界定最后一项，$\sum_b2^bD_b/H\le\sum_{p\le X,\,\chi(p)=1} k_v(u-\log p)$。在 $p>\sqrt Y$ 上插入 $2\log p/u$，对分裂素数用 $Q_v$ 的界；在 $p\le\sqrt Y$ 上用 $k_v\le1$，总量为 $O(\sqrt Y)$。式 (7) 的主项系数为 $1/2$，(6) 需要 $1/3$。这个差在任意远窗口都存在，增大权或缩小取样网格无法从 (7) 补回它。

等价地，也可寻求本原多项式的**高度下界**。令
$$
 R(u)=\frac1{3dH}\sum_b2^b\log\mathcal M(W_b/g_b).
$$
前篇的本原高度公式和同样的取整、截尾给出
$$
 R(u)=P_v(u)+O\!\left(\sqrt Y+\frac{e^{u+v/2}}{du}\right).
$$
若选 $d=6\lceil Y/6\rceil$，两项误差合计为 $O(\sqrt Y)$。于是只要在 (6) 的同样区段证明 $R(u)\ge e^{u+v/2}/6-O_{\lambda,\varepsilon}(e^{(\sigma_\lambda+\varepsilon)u})$ 即可。完全分裂行在本原高度中承担这项下界，三循环行的本原 Mahler 测度精确为 1。这给混合 theta 提升的度量估计指定了应控制的对象和方向。本文证明了取样规模、整数重数及未约分高度的上下界，(6) 的中心系数仍是待证的全局算术估计。

## 在一个非零导数上把额外高度降到平方根量级

还有一种方式可以避开 (7) 的主项损失，且只需权 7 的 $\Theta_6$。沿用相同的位集和素数截尾，在每一非空层内按有理素数的大小排序，每 $s_b=\lceil\sqrt{D_b}\rceil$ 个素数分成一组，末组可较短。这个分组只用 $\chi$ 和 $m_p$，无需先知道 $F,C$。取
$$
 n=\lceil e^{u/2+3v/8}\rceil,\quad b_p=B_6(p),\qquad
 Z_I(T)=\prod_{p\in I}\big(p^n(1+T)+b_p\big)
$$
并记 $J_I$ 为 $Z_I$ 在 $-1$ 处**首个非零 Hasse 导数的绝对值**。每组依次取 $0,\ldots,|I|$ 阶必能找到它。每层的组数、组内最高阶数均至多 $s_b$。

**定理 3。** 对一个整数尺度，置 $\mathcal J=\prod_b\prod_{I\ {\rm in\ layer}\ b}J_I^{\,2^b}$，则
$$
 \frac{\log\mathcal J}{nH}=B_v(u)+O(e^{u/2+v/8}). \tag{8}
$$
这个误差界可用绝对常数，在前述 $1\le v=\lambda u\le u/16$ 上一致，对 $L_k(s,\psi)$ 没有零点假设。指数 $n$ 的大小是 $Y^{1/2+3\lambda/8}$。它和每层组数、最高导数阶数的幂指数，在 $\lambda\downarrow0$ 时都趋近 $1/2$。

**证明。** 在 $p\in C$ 时 $b_p=0$。在 $p\in F$ 时 $b_p$ 是非零整数，$|b_p|\le 2p^3\le p^4$。因此首个非零阶恰为 $|I\cap C|$，Taylor 展开给出精确公式
$$
 J_I=\prod_{p\in I\cap C}p^n\prod_{p\in I\cap F}|b_p|. \tag{9}
$$
例如 $I=\{13,59\}$ 时，$b_{13}=0$，$59=(6+\sqrt{-23})(6-\sqrt{-23})$ 给 $b_{59}=-253942$。首个非零导数为一阶，$J_I=253942\cdot13^n$。再合并所有组和二进位，
$$
 0\le \frac{\log\mathcal J}{nH}
 -\frac1H\sum_{\substack{p\in C\\p\le X,\ p\ne3}}m_p\log p
 \le\frac{4Q_v(u)}n .
$$
Chebyshev 界已给 $Q_v\ll e^{u+v/2}$。具体地，对 $\sum_{p\le x}\log p\ll x$ 分部积分，权导数的积分至多为 $e^{u+v/2}(1+1/\sqrt v)$ 乘绝对常数；这是高斯倾斜后均值为 $v$、方差为 $v$ 的绝对一阶矩除以 $v$。代入 $n$，再用 (4) 的取整和截尾误差，得 (8)。

这次选取首个非零导数，在完全分裂处只付出 $\log|b_p|\le4\log p$，三循环素数仍各留下 $n\log p$。若取全部导数的最大值，则每个素数都可能留下约 $n\log p$。权固定而深度增大的来源是显式插入的有理整数 $p^n$；这个参数尚未被识别为原提升链的长度。

这些组多项式还可以用前篇的 $P(T_p+p^n(1+T)\operatorname{Id})P$ 在固定的三维整数格上依次实现。若逐个使用 Petersson 算子范数，每个非空组的投影数为 $|I|+1\le2|I|$。所以所有组纯由投影引入的额外对数上界，在归一化后至多为
$$
 \frac{2\log\|P\|_{\rm Pet}}{nH}\sum_b2^bD_b
 \ll \frac{e^{u/2+v/8}}{u}+\frac{\sqrt Y}{n}.
$$
这里用的是 (7) 证明中对 $D_b$ 之和的估计，以及 (8) 所选的 $n$。这个界只计算投影步骤，Hecke 因子组合后首个非零导数的整体上界仍须单独估计。[整数格和投影范数的证明](cubic23-regular-theta-packet.md)

由固定域的有效 Chebotarev 定理，$\sum_{p\le x,\ p\in C}\log p=x/3+O_M(xe^{-c\sqrt{\log x}})$，其中 $c>0$。该固定域的任何例外实零点项可吸收到此余项。对 $x\ge\sqrt Y$ 分部积分，取出 $e^{-c\sqrt{u/2}}$ 并用上一段的一阶矩；对 $x<\sqrt Y$ 用 $O(\sqrt Y)$ 界定素数和与连续主项。由 (8) 得某个 $c'>0$ 满足
$$
 \frac{\log\mathcal J}{nH}
 =\frac{e^{u+v/2}}3+O_M(e^{u+v/2-c'\sqrt u}). \tag{10}
$$
因此这个整数的正确主项已经可证。要排除余下的偏线零点，(6) 还要求在趋零斜率的一列远端区段，把 (10) 的**上侧**误差压到 $O_{\lambda,\varepsilon}(e^{(1/2+\lambda/8+\varepsilon)u})$。式 (8) 表明，局部杂因子和所选导数本身只贡献允许的误差。真正缺少的是不先知道 $|I\cap C|$ 就能给这些非零导数的整体大小定上界的算术论证。原链的同一步内积公式还未提供这项估计。[Lagarias–Odlyzko，有效 Chebotarev 定理 1.3](https://aareyanmanzoor.github.io/assets/articles/lagarias-odlyzko.pdf) · [已有零点证书给出的更细斜率依赖界](cubic23-class-character-mean-square.md)

结合仓库已有的零点证书，(8) 在 $\lambda=2^{-26}$ 已给出 (6) 所需精度的导数上界，并且对所有充分大尺度成立。对几乎每个 $t\in[1,2]$，沿 $\lambda_j=t2^{-j}$（$j\ge6$），每个 $\varepsilon>0$ 下，$\log\mathcal J/(nH)\le e^{u+v/2}/3+e^{(\sigma_{\lambda_j}+\varepsilon)u}$ 的整数尺度下密度至少为 $1/2$。证明是将前稿的完全分裂下界取容差 $\varepsilon/2$，再用 $B_v=Q_v-P_v$ 和 (8) 吸收常数。尚须将这些尺度连成长达 $\ell(N)$ 的任意远区段，或按前稿的密度判据取得严格超过一半的平均比例。此归约沿用原有斜率界。[高度 8192 的证书](cubic23-height-8192.md) · [密度与斜率的精确关系](cubic23-class-character-mean-square.md)

复算程序对次数 8、包含三种 Frobenius 类型的整数乘积逐一取 $1\le s\le9$，用两组圆分余式的实际系数核对精确内容，并用二进分层复原同一整数。另以 $n=11$ 在七个分裂于 $k$ 的素数上分组，逐组实际计算导数，重组后的 592 位整数与 (9) 完全一致。[复算程序](verify_theta23_stepanov.py) · [已核验的局部代数与重数引理](Theta23Local.lean)
