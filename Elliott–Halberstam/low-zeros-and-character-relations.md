---
title: "Elliott–Halberstam 猜想需要哪些 Dirichlet 零点？"
date: "2026-11-25"
abstract: "从 Alpöge–Claude 归档的 Dirichlet-GRH 出发，将最大剩余类、最大端点的误差归约到大导子、高阶特征在对数幂高度内的有限零点和。四阶矩中满足短乘法关系的部分可按绝对值估计。若在三分之二以下出现足够大的等差数列误差，其余四特征组合必在一个共同端点产生正的幂次超额。文末以七个小导子特征的完整零点计数核对显式公式。"
---

# Elliott–Halberstam 猜想需要哪些 Dirichlet 零点？

设 $\Lambda$ 是 von Mangoldt 函数，$\psi(y;q,a)=\sum_{n\le y,\,n\equiv a\bmod q}\Lambda(n)$。对每个模数分别挑选剩余类和端点，记
$$
E(x,Q)=\sum_{q\le Q}\max_{a\in(\mathbb Z/q\mathbb Z)^\times}
  \sup_{2\le y\le x}\left|\psi(y;q,a)-\frac{y}{\varphi(q)}\right|.
$$
Elliott–Halberstam 猜想（EH）要求对任意固定的 $0<\varepsilon<1$ 和 $A>0$，有 $E(x,x^{1-\varepsilon})\ll_{A,\varepsilon}x/(\log x)^A$。Bombieri–Vinogradov 定理已无条件达到 $Q=\sqrt{x}/(\log x)^{B(A)}$。逐个特征使用 GRH 给出的经典界为 $O(\sqrt{x}\log^2 x)$，对模数直接求和得到 $O(Q\sqrt{x}\log^2 x)$。固定剩余类并对模数的因子结构加以限制时，Maynard 已在超过平方根的范围证明等分布。[等差数列中的素数，第 20 章，式 (20.17)–(20.19) 及猜想 20.4](https://public.websites.umich.edu/~hlm/math775/ch20.pdf) · [Maynard，固定剩余类](https://arxiv.org/abs/2006.06572)

下文把 Dirichlet-GRH 作为解析输入。**定理 1** 给出达到误差 $x/(\log x)^A$ 所需保留的零点、特征和端点。**定理 2** 进一步删去一批可按绝对值求和的四特征组合。剩余四特征和只要节省任一固定幂，就能在这一判据中越过 $1/2$。能否取得这个节省仍是开放的估计。

## 1. 原证明的窗口如何依赖特征

Alpöge 于 2026 年 11 月 20 日发布了与 Claude 合作的 Lean 归档。其 `LFunction_eq_zero_re_of_isPrimitive` 把本原 Dirichlet 特征的每个非平凡零点写成 $\rho=1/2+i\gamma$，包括可能的 $\gamma=0$。Sutherland 给出了五条最终声明的核验记录。[作者声明](https://alpo.ge/remark.pdf) · [原始版本的最终定理](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Main.lean#L4-L12) · [核验记录](https://math.mit.edu/~drew/grh/)

原稿的 `winPairs` 先固定本原特征 $\psi$，随后选择参数 $J$，并在任意大的某些 $[Z,Z^{3/2}]$ 上为每个 $x$ 找到窗口与辅助整数。`plugWin` 对每个 $\psi$ 分别调用 `supplyChi`。因此这一层声明允许不同特征使用不同的窗口。下文给出了适合回到这一步检验的联合估计，其中四个特征在同一个端点 $y$ 取值。[窗口声明及量词](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/CharacterBoxes/WinPairs.lean#L5-L10) · [逐特征调用](https://github.com/leanprover-community/riemann-hypothesis/blob/40e56bab026552e089d7afcfeff08d2defb9c005/src/Development/AutomorphicArithmetic/CharacterBoxes/EveryCharacter.lean#L8-L18)

已有的两条相关路线给出了参照。Kandhil、Languasco 与 Moree 在 GRH 和其加权零点对关联猜想下推出了校正的 Montgomery 猜想及 EH。de la Bretèche 与 Fiorilli 研究了等差数列素数的平滑高阶矩及其特征组合。本文给出的结论是固定误差精度下的有限截取，以及第四矩中短乘法关系的总量估计。全文的 EH 分布指数推论均注明所需的附加估计。[Kandhil–Languasco–Moree，定理 4 与推论 5](https://doi.org/10.1007/s00208-026-03383-y) · [de la Bretèche–Fiorilli](https://arxiv.org/abs/2010.05944)

## 2. 有限零点和与最大剩余类

固定 $A>0$，取 $x$ 充分大，令 $L=\log x$，并设
$$
 h=L^{-A-4},\qquad T=4h^{-1}\sqrt L,\qquad R=\sqrt{x}\,L^{-A-4},
 \qquad J=\lceil L^{A+4}\rceil,\qquad y_j=jx/J\quad(1\le j\le J).
$$
记 $\omega(q)$ 为 $q$ 的不同素因子数，另设
$$
 K=\left\lceil (A+6)\frac{\log L}{\log 2}\right\rceil,\qquad
 V=x^{1/[4(K+2)]}.
$$
模 $q$ 的特征 $\chi$ 由本原特征 $\chi^*$ 诱导，导子记作 $r_\chi$，阶记作 $\operatorname{ord}(\chi)$。对 $\omega(q)\le K$、$r_\chi>R$、$\operatorname{ord}(\chi)>V$ 的特征置
$$
 z_\chi(y)=
 \sum_{\substack{\rho=1/2+i\gamma\ \text{of }L(s,\chi^*)\\|\gamma|\le T}}
     \frac{y^{i\gamma}\exp(h^2\rho^2/2)}{\rho}.
$$
其余特征取 $z_\chi(y)=0$。零点计重数，两种符号的高度分别保留。特别地，中心零点的每一重贡献 $2\exp(h^2/8)$。对复特征不能用本特征的正高度代替全部高度。定义
$$
 Z_q(a,y)=\frac{1}{\varphi(q)}\sum_{\chi\bmod q}\overline{\chi(a)}\,z_\chi(y).
$$

**定理 1（有限截取）**　在 Dirichlet-GRH 下，对 $1\le Q\le x^{1-\varepsilon}$ 一致有
$$
 \left|E(x,Q)-
 \sum_{\substack{q\le Q\\\omega(q)\le K}}
 \max_{\substack{a\in(\mathbb Z/q\mathbb Z)^\times\\1\le j\le J}}
       \sqrt{y_j}\,|Z_q(a,y_j)|\right|
 \ll_{A,\varepsilon} xL^{-A}. \tag{1}
$$
所以在给定精度下，余下的问题只涉及 $|\gamma|\le4(\log x)^{A+9/2}$。导子至少接近 $\sqrt{x}$，可保留的特征阶下界 $V$ 最终超过 $\log x$ 的任一固定次幂。模数仍覆盖原猜想的整个范围。只取 $y=x$ 时，同一证明给出固定端点版，无需端点网格。

### 在剩余类中平滑

取标准正态变量 $U$。把区间 $(1,y]$ 换成随机区间 $(e^{hU},ye^{hU}]$ 并取期望。一个长度为 $b$ 的区间内，单个模 $q$ 剩余类至多有 $b/q+1$ 个整数。由 $\Lambda(n)\le\log n$、$\mathbb E|e^{hU}-1|\ll h$ 以及正态分布的指数矩，平滑前后单个剩余类的差为 $O((hx/q+1)L)$，对 $2\le y\le x$ 一致。具体地，上端点处可逐样本用 $(y|e^{hU}-1|/q+1)(L+h|U|)$ 控制，下端点的期望为 $O(1)$。同一误差界适用于对既约剩余类取平均后的量。

先扣掉既约剩余类的实际平均值。这正好除去主特征。用 $\zeta$ 的 RH 估计 $\psi(y)-y\ll\sqrt{y}\log^2(2y)$，再扣除 $p\mid q$ 的素数幂，将实际平均值换回 $y/\varphi(q)$ 的代价是 $O((\sqrt{x}L^2+\omega(q)L)/\varphi(q))$。任意端点移动到最近的网格点，长度至多 $x/J$，代价是 $O((x/(Jq)+1)L+x/(J\varphi(q)))$。利用
$$
 \sum_{n\le u}\frac1{\varphi(n)}\ll\log(2u)
$$
求和，上述误差合计为 $O(x(h+J^{-1})L^2+QL+\sqrt{x}L^3)$。这里的倒数和界可由 $n/\varphi(n)=\sum_{d\mid n}\mu^2(d)/\varphi(d)$ 及 $\sum_d\mu^2(d)/(d\varphi(d))<\infty$ 得到。

### 移线与高度尾项

对每个非主特征，从 $\chi$ 换成 $\eta=\chi^*$ 只补入 $p\mid q,\ p\nmid r_\chi$ 的素数幂。每个素数的总权重在取期望后为 $O(L)$，故每个特征的误差为 $O(L^2)$。反演时有因子 $1/\varphi(q)$，所有模数的这一步代价为 $O(QL^2)$。

本原非主特征的平滑和是 $-L'(s,\eta)/L(s,\eta)$ 乘以下式在 $\Re s=2$ 上的 Mellin 积分，积分测度为 $ds/(2\pi i)$。
$$
 F_y(s)=\frac{y^s-1}{s}\exp(h^2s^2/2),\qquad F_y(0)=\log y.
$$
移到 $\Re s=-1/2$，每个非平凡零点贡献 $-F_y(\rho)$。偶本原非主特征在 $s=0$ 的平凡零点贡献 $-\log y$。在新直线上，由函数方程、$\Re(1-s)=3/2$ 上绝对收敛的对数导数以及 Gamma 因子，有 $|L'/L|\ll\log(2r_\chi(2+|t|))$。又有 $|F_y(-1/2+it)|\ll e^{-h^2t^2/2}/(1+|t|)$，所以新直线积分是 $O_A(L^2)$。移线可先在避开零点纵坐标的高度截断，再用高斯衰减令高度趋于无穷。具体地，单位高度零点数为 $O(\log(2r_\chi(2+|t|)))$，在每个充分高的单位区间可选一条距附近纵坐标 $\gg1/\log(2r_\chi(2+|t|))$ 的横线，对数导数的局部部分分式展开在其上给出对数平方界。[Bennett–Martin–O’Bryant–Rechnitzer 的零点计数估计](https://doi.org/10.1090/mcom/3599)

同一计数估计按单位高度求和，给出
$$
 \sum_\rho\frac{e^{-h^2\gamma^2/4}}{1+|\gamma|}\ll_A L^2.
$$
因 $|\rho|\asymp 1+|\gamma|$，删去 $F_y(\rho)$ 中来自下端点 $1$ 的部分，代价为 $O_A(L^2)$。在 $|\gamma|>T$ 的余项中另取出 $e^{-h^2T^2/4}=x^{-4}$，故上端点的尾项是 $O_A(\sqrt{x}L^2 x^{-4})$。这样每个本原非主特征的平滑和等于 $-\sqrt{y}$ 乘以未施加导子和阶筛选的 $z_\chi(y)$，误差为 $O_A(L^2)$。每个零点系数和的绝对值亦为 $O_A(L^2)$。

### 小导子、低阶与少数模数的贡献

每个导子 $r\mid q$ 至多对应 $\varphi(r)$ 个模 $q$ 特征。由于 $\varphi(rm)\ge\varphi(r)\varphi(m)$，所有 $r\le R$ 的贡献，连同模数 $q$ 的全部倍数重数，至多为
$$
 \ll_A \sqrt{x}L^2\sum_{r\le R}\varphi(r)
             \sum_{m\le Q/r}\frac1{\varphi(rm)}
 \ll_A\sqrt{x}L^3R=xL^{-A-1}. \tag{2}
$$

对 $\omega(q)>K$ 先在素数一侧估计整个模数。单模数误差 $\ll (x/q+1)L+x/\varphi(q)$。用 $q/\varphi(q)\le2^{\omega(q)}$、$2^{\omega(q)}\le d_2(q)$、$4^{\omega(q)}\le d_4(q)$，其中 $d_j$ 是 $j$ 重除数函数，有
$$
 \sum_{\substack{q\le Q\\\omega(q)>K}}
 \left((x/q+1)L+x/\varphi(q)\right)\ll 2^{-K}x(L^3+L^4)\ll_A xL^{-A-2}.
$$
这里用到了 $q\le x$ 和 $\sum_{q\le Q}d_j(q)/q\le(1+\log Q)^j$。

在其余模数上，特征群至多由 $\omega(q)+1\le K+1$ 个元素生成。这来自奇素数幂的单位群为循环群，$2$ 的幂至多多出一个生成元。阶整除 $d$ 的特征至多有 $d^{K+1}$ 个，因此阶不超过 $V$ 的特征至多有 $V^{K+2}=x^{1/4}$ 个。它们在所有这些模数上的总贡献是 $O_A(x^{3/4}L^3)$。把各项合并，由 $Q\le x^{1-\varepsilon}$ 得 (1)。所有界在单个模数内都对剩余类和端点一致，取最大值保留误差界。

## 3. 四个特征之间的短乘法关系

在模 $q$ 的特征群 $G$ 中令 $g=|G|=\varphi(q)$，在同一端点省略 $y$。正交性给出精确恒等式
$$
 \frac1g\sum_{a\in(\mathbb Z/q\mathbb Z)^\times}|Z_q(a,y)|^4
 =\frac1{g^4}\sum_{\chi,\eta,\nu\in G}
 z_\chi z_{\chi\eta\nu}\overline{z_{\chi\eta}z_{\chi\nu}}. \tag{3}
$$
任意满足 $\chi_1\chi_2=\chi_3\chi_4$ 的有序四元组，在右侧恰出现一次。$\eta,\nu$ 分别是 $\chi_3/\chi_1$ 与 $\chi_4/\chi_1$。

取 $S=\exp(\sqrt L)$。对每个三元组 $(\chi,\eta,\nu)$，检验是否存在整数 $m,n$ 满足
$$
 0<\max(|m|,|n|)\le S,\qquad
 \min\bigl(r_{\eta^m\nu^n},\operatorname{ord}(\eta^m\nu^n)\bigr)\le S. \tag{4}
$$
这包括有特征重复的项、两个特征相差一个小导子扭转的项以及 $\eta=\nu^2$ 之类的关系。记 $H_q(y)$ 为 (3) 右侧其余项的实部之和，包含因子 $g^{-4}$。

**定理 2（可按绝对值估计的四元组）**　对 $\omega(q)\le K$，(3) 中满足 (4) 的所有项，其绝对值之和至多为
$$
 \frac{\exp(O_A(\sqrt L\log L))}{\varphi(q)^2}.
 \tag{5}
$$
因而在 $D<q\le\min(2D,Q)$ 上求和，代价是 $D^{-1}x^{o_A(1)}$，对所有网格端点一致。

**证明。** 固定 $m,n$，令 $d=\gcd(|m|,|n|)>0$。群同态 $(\eta,\nu)\mapsto\eta^m\nu^n$ 的像是 $G^d$，因为 Bézout 等式可以产生任意 $d$ 次幂。每个非空纤维的大小恰为 $g|G[d]|$，其中 $G[d]=\{\tau:\tau^d=1\}$。因此每个纤维至多有 $gS^{K+1}$ 个点。目标特征中导子不超过 $S$ 的至多有 $\sum_{r\le S}\varphi(r)\ll S^2$ 个，阶不超过 $S$ 的至多有 $S^{K+2}$ 个。再乘以 $O(S^2)$ 组整数 $(m,n)$ 和 $g$ 种基准特征 $\chi$，合计至多 $g^2\exp(O_A(\sqrt L\log L))$ 个三元组。每项的未归一化绝对值为 $O_A(L^8)$，得到 (5)。最后 $q/\varphi(q)\le 2^K\ll_A L^{A+6}$，故模数块上的总量有上述界。

## 4. 四特征估计与分布指数

**条件推论。** 固定 $0\le\lambda\le1$ 及 $0<\theta<2/(3+\lambda)$。若对每个 $A>0$，在 $Q=x^\theta$ 的各个二倍模数区间、每个 $y_j$ 上，一致有单侧上界
$$
 \sum_{\substack{D<q\le\min(2D,Q)\\\omega(q)\le K}}
 H_q(y_j)\ \le\ D^{-1+\lambda}x^{o(1)}, \tag{6}
$$
则 $E(x,x^\theta)\ll_{A,\theta}xL^{-A}$。这里 $x^{o(1)}$ 表示对每个固定 $\delta>0$ 都可用 $O_{A,\theta,\delta}(x^\delta)$ 控制，且对 $D,j$ 一致。

**证明。** 由 (3)、(5)、(6)，同一块上各模数的归一化第四矩之和至多 $D^{-1+\lambda}x^{o(1)}$。用 Hölder 不等式对模数求和，再将每个模数的最大四次幂放入对 $a,j$ 的总和。去掉剩余类平均的代价至多为 $2D$，端点数 $J$ 是对数幂。该块对 (1) 中零点和的贡献至多为
$$
 \sqrt{x}\, D^{3/4}\big(2D\,J\,D^{-1+\lambda}x^{o(1)}\big)^{1/4}
 \ll \sqrt{x}\,D^{(3+\lambda)/4}x^{o(1)}.
$$
对至多 $O(\log x)$ 个块求和，严格不等式 $\theta<2/(3+\lambda)$ 留出的固定幂足以吸收任意对数幂。

| 剩余四特征和 (6) 的大小，略去 $x^{o(1)}$ | 得到的全部模数、最大剩余类范围 |
| --- | --- |
| $1$（直接对三个特征求和即可得到） | $\theta<1/2$ |
| $D^{-1/2}$（待证） | $\theta<4/7$ |
| $D^{-1}$（待证） | $\theta<2/3$ |

还有一条通向完整 EH 的条件推论。若将 (3) 左侧换成归一化的 $2k$ 阶矩，并能对每个固定整数 $k\ge2$ 把其模数块之和控制到 $D^{1-k}x^{o(1)}$，则同一 Hölder 步骤给出块界 $\sqrt{x}D^{1/2+1/(2k)}x^{o(1)}$，可取任意 $\theta<k/(k+1)$。让 $k$ 随目标指数选取，即覆盖任意固定 $\theta<1$。定理 2 对短关系的计数专用于四阶。更高阶矩所需的估计与 (6) 都有待证明。

## 5. 如果在分布指数三分之二以下出现大误差

还有一个反向结论，直接限制零点在共同端点上的叠加。固定 $B>0$ 和 $1/2<\theta<2/3$。若某个充分大的 $x$ 满足 $E(x,x^\theta)\ge x/L^B$，用 $A=B+2$ 构造前述零点和。则必有一个模数块和**一个共同的** $y_j$ 使
$$
 D\gg_{B,\theta}\frac{\sqrt{x}}{L^{B+3}},\qquad D\le x^\theta,\qquad
 \sum_{\substack{D<q\le\min(2D,x^\theta)\\\omega(q)\le K}} H_q(y_j)
 \gg_{B,\theta}\frac{x^2}{D^4 L^{5B+10}}. \tag{7}
$$
例如取 $\theta=3/5$，(7) 的下界至少是 $D^{-1}$ 的 $x^{1/5}/L^{5B+10}$ 倍（差一个固定常数）。这里 $H_q$ 只计入不满足 (4) 的四特征组合。因此若在此范围出现 EH 量级的大误差，零点必须在同一尺度上产生一个正的、超过 (5) 尺度固定幂的四阶和。

**证明。** (1) 的误差为 $o(x/L^B)$。至多 $O(L)$ 个模数块中，有一块使 $\sum_q b_q\gg\sqrt{x}/L^{B+1}$，其中 $b_q=\max_{a,j}|Z_q(a,y_j)|$，因为 $\sqrt{y_j}\le\sqrt{x}$。由 $b_q\ll_A L^2$ 即得 $D$ 的下界。同一块上 Hölder 不等式给出 $\sum_q b_q^4\gg x^2/(D^3L^{4B+4})$。将最大值放入对 $a,j$ 的总和，除以最多 $2D$ 个剩余类和 $J\ll L^{B+6}$ 个端点，至少一个 $y_j$ 处的归一化第四矩之和达到 (7) 右侧的量级。扣去满足 (4) 的项至多损失 $D^{-1}x^{o(1)}$。两者之比由 $x^2/(D^3L^{5B+10})\ge x^{2-3\theta}/L^{5B+10}$ 控制，后者比任意次幂的对数都大，且超过该 $x^{o(1)}$ 因子。故余下的实部之和仍有 (7) 的正下界。

## 6. 有限零点计数与复算

取 $h=0.30$、$y=1024$，使用 $(y/2,y]$ 的高斯平滑权。这等于两个单端点权之差，零点项多出因子 $1-2^{-\rho}$。复算程序从 Hardy $Z$ 函数的变号区间求零点，又用辐角原理在整个临界带独立计数。对这里七个本原特征，区间算术在两条水平线 $t=\pm35$ 的 $1/2\le\sigma\le2$ 上均验证 $\Re L(\sigma+it,\chi)>0$。在 $\sigma=2$ 上也有 $|L-1|\le\zeta(2)-1<1$。于是原始 $L$ 函数在右半边界上的辐角可连续取主值。用函数方程反射补全矩形的左半边界，零点总数为下式包住的唯一整数，其中 $\kappa$ 是特征的奇偶性，并取 $t=35$。
$$
 \frac{t\log(r/\pi)+2\Im\log\Gamma((1/2+\kappa+it)/2)
 +\arg L(1/2+it,\chi)-\arg L(1/2-it,\chi)}{\pi}.
$$
这里取右半平面上的解析 $\log\Gamma$。每一个计数单位都对应临界线上一个独立的变号区间，故该有限范围内的零点也都是单零点。

| 本原导子 | Conrey 特征编号 | 每个特征在 $\lvert\gamma\rvert\le35$ 内的零点数 |
| --- | --- | --- |
| $3$ | $2$ | $20$ |
| $5$ | $2,3,4$ | 各 $26$ |
| $15$ | $2,8,14$ | 各 $38$ |

共计 212 个，按特征分别计数。将这些纵坐标代入零点侧，补上移线至 $\Re s=-7/2$ 跨过的平凡零点及诱导特征缺少的 Euler 因子，再与直接枚举素数幂所得的各既约剩余类中心化计数比较。模 $5,10,15$ 的最大绝对差依次为 $4.98\times10^{-12}$、$5.00\times10^{-12}$、$7.12\times10^{-12}$。小样本取 $r_\chi>3$ 且 $\operatorname{ord}(\chi)>2$ 检验 (3)，直接剩余类平均与特征卷积的差小于 $4\times10^{-21}$。程序还在六个模数上以整数运算核对了 30 组幂映射的纤维。辐角计数使用球算术。两侧求和的差值是浮点数值核对。检验范围是上述有限样本，(6) 所需的是随 $x$ 增长的一致估计。

使用 Python 3 和 `python-flint==0.9.0`，在本目录运行 `python verify.py` 即可重算零点、剩余类和四阶矩。[完整复算程序](verify.py)
