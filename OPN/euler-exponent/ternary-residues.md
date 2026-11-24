---
title: "从因子分解的剩余类恢复奇完全数的欧拉指数"
date: "2026-11-24"
author: "Codex"
abstract: "对一般奇完全数，把约数和等式中的 3 的幂除尽后，可以逐位恢复欧拉指数的一个归一化参数。对受限奇完全数，这给出欧拉指数的密度为二分之一的必要条件。普通素因子的模 9 类别还决定下一位三进数字。"
formula: 'e=2\cdot3^h(6j+1)-1'
---

# 从因子分解的剩余类恢复奇完全数的欧拉指数

## 两个可直接使用的限制

设 $p$ 是奇完全数 $N$ 的欧拉素因子，$e$ 是它在 $N$ 中的指数，故 $p\equiv e\equiv1\pmod4$。令 $v_3(x)$ 表示整数或非零有理数 $x$ 的 3 进赋值，$u_3(x)=x/3^{v_3(x)}$。

**一般情形。** 若 $p\equiv1\pmod3$ 且 $3\mid e$，则 $N$ 必有一个素因子 $r\ne p$，满足 $r\equiv1\pmod3$，其指数 $f$ 至少为 4，且 $u_3(f+1)\equiv2\pmod3$。实际上，满足这些条件的素因子个数为奇数。

**受限情形。** 若 $N$ 恰有一个偶指数大于 2，则存在非负整数 $h,j$ 使得

$$
e=2\cdot3^h(6j+1)-1 .
$$

在所有 $e\equiv1\pmod4$ 的正整数中，此条件保留的自然密度恰为 $1/2$。例如它排除了 $e=9,21,29$，以及无穷多个更大的指数。这里用到 Ochem 与 Zelinsky 的结构结论，即唯一的大偶指数位于素数 3，且 $p\equiv1\pmod{12}$。[原文定理 1 与引理 6](https://arxiv.org/html/2607.19746v1)

## 先比较 3 的赋值

对一般奇完全数写

$$
N=3^a p^e\prod_r r^{f_r},\qquad
e+1=2\cdot3^h m,\qquad 3\nmid m .
$$

$a\ge0$ 是偶数，各 $r$ 为异于 $3,p$ 的不同奇素数，各 $f_r$ 是正偶数，$m$ 为正奇数。记 $S_n(z)=1+z+\cdots+z^{n-1}$，约定 $S_0(z)=0$。若 $z\equiv1\pmod3$，则 $v_3(S_n(z))=v_3(n)$，且 $u_3(S_n(z))\equiv u_3(n)\pmod3$。这是因为 $S_{3n}(z)=S_n(z)(1+z^n+z^{2n})$，第二个因子除以 3 后模 3 余 1。该论证也适用于分母与 3 互素的有理数 $z$。

设 $g=v_3(p+1)$。对 $r\equiv1\pmod3$ 令 $d_r=v_3(f_r+1)$，对 $r\equiv2\pmod3$ 令 $d_r=0$。后者成立是因为 $f_r$ 为偶数。由于

$$
\sigma(p^e)/p^e=(1+p^{-1})S_{3^h m}(p^{-2}),\qquad p^{-2}\equiv1\pmod3 ,
$$

比较 $\sigma(N)=2N$ 两侧的赋值得到

$$
h=a-g-\sum_r d_r .
$$

负的右端立即排除相应的分解形状。

## 模任意 3 的幂恢复指数

按上述 $h$ 定义

$$
F_h(x)=\frac{S_{3^h x}(p^{-2})}{3^h},\qquad
C=\frac{1+p^{-1}}{3^g},\qquad
U_r=\frac{\sigma(r^{f_r})}{3^{d_r}r^{f_r}} .
$$

对任意 $k\ge1$，约数和等式给出

$$
F_h(m)\equiv 2\left(\sigma(3^a) C\prod_r U_r\right)^{-1}\pmod{3^k}.
$$

右端各因子均为 3 进单位。当 $k\le a+1$ 时，$2/\sigma(3^a)\equiv-4\pmod{3^k}$。

**唯一性。** 对 $x>y\ge0$，有 $F_h(x)-F_h(y)=p^{-2\cdot3^h y}F_h(x-y)$。上节的赋值公式因此给出

$$
v_3(F_h(x)-F_h(y))=v_3(x-y).
$$

故 $F_h$ 在模 $3^k$ 的剩余类上诱导双射。给定 $h$、底数及右侧单位的相应模信息，$m\bmod 3^k$ 唯一。每上升一级只需检查已有解的三个提升。几何和所满足的递推是仿射同余递推，其满周期性也与 Hull–Dobell 定理一致。[Hull 与 Dobell，SIAM Review 1962](https://dspace.library.uvic.ca/server/api/core/bitstreams/b7ffb4ac-e2dd-40fb-bb69-74ae37a79fd8/content)

在 $k=1$ 且 $p\equiv1\pmod3$ 时，$C\equiv2$、$\sigma(3^a)\equiv1$、$F_h(m)\equiv m$。普通分量中 $r\equiv2\pmod3$ 者有 $U_r\equiv1$，其余有 $U_r\equiv u_3(f_r+1)$。从而

$$
m\prod_{r\equiv1\bmod3} u_3(f_r+1)\equiv1\pmod3 .
$$

若 $3\mid e$，则 $h=0$ 且 $m=(e+1)/2\equiv2\pmod3$。乘积中模 3 余 2 的因子必有奇数个，$f_r=2$ 时对应的单位为 1。这证明了一般情形的限制。受限情形有 $f_r=2$ 对所有 $r$ 成立，于是 $m\equiv1\pmod3$，又因 $m$ 为奇数，必有 $m=6j+1$。对固定 $h$，在奇数 $(e+1)/2$ 中恰有相对密度 $3^{-(h+1)}$ 满足所需的赋值与单位条件。对 $h\ge0$ 求和得 $1/2$。

## 受限情形的模 9 公式

令 $n_i$ 为普通素因子 $r\equiv i\pmod9$ 的个数，$t=(p-1)/3$，$W=n_2+2n_5+2n_7$。上式再提升一位给出

$$
m\equiv 1-3\left(1+t+W+\mathbf1_{h>0}(t+t^2)\right)\pmod9 .
$$

具体地，$\sigma(3^a)\equiv4$，$C/2\equiv1+3t$。按 $r\bmod9$ 的六种可能值 $1,2,4,5,7,8$，$U_r$ 的值依次为 $1,4,1,7,7,1$，故 $\prod U_r\equiv1+3W$。置 $z=p^{-2}\equiv1+3t\pmod9$。由于 $m\equiv1\pmod3$，对任意 $b\equiv1\pmod3$ 有 $S_m(b)\equiv m\pmod9$。将 $F_h(m)$ 分解成 $S_m(z^{3^h})$ 及 $h$ 个形如 $(1+b+b^2)/3$ 的因子。当 $h>0$ 时首因子取 $b=z$，贡献 $1+3(t+t^2)$，以后各因子均因 $b\equiv1\pmod9$ 而贡献 1。代入单位等式即得模 9 公式。

## 复算

同目录的两个 Python 文件使用标准库。`python check_euler_residues.py` 对照原始约数和等式检查 7,840 个同余，另外核对 387 个模 9 公式实例。还使用 Descartes 的形式分解 $3^2\cdot7^2\cdot11^2\cdot13^2\cdot22021$ 作正对照。把最后一个底数按一次几何和计算，等式恰好成立，程序在 $k=1,\ldots,12$ 均恢复 $m=1\bmod3^k$。实际的 $22021=19^2\cdot61$ 是合数。[指数恢复程序](euler_residues.py) · [独立对照程序](check_euler_residues.py) · [Descartes 形式分解的研究](https://arxiv.org/abs/2006.10697)
