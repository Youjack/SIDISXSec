# Design of `RC_conv_xsec` with SIDIS threshold

## The integration region

The integration region of $\xi$ and $\zeta$ is restricted by
$$
\xi_m(\zeta) = \frac{B+R\zeta}{A\zeta-1},\quad \zeta_m(\xi) = \frac{B+\xi}{A\xi-R},
$$
where
$$
A \equiv \frac{l\cdot(P-P_h)}{l\cdot l'}, \quad B \equiv \frac{l'\cdot(P-P_h)}{l\cdot l'}, \quad R \equiv \frac{z_h}{x_B} + \frac{M_\mathrm{th}^2-(M^2+M_h^2)}{Q^2}.
$$
Both $\xi_m(\zeta)$ and $\zeta_m(\xi)$ should be monotonically decreasing functions to meet physical conditions.

We define $\tilde{\xi}$ and $\tilde{\zeta}$ such that
$$\begin{aligned}
&\xi \equiv \frac{\tilde{\xi}-\xi_m(1)}{1-\xi_m(1)} + \xi_m(\zeta) \frac{1-\tilde{\xi}}{1-\xi_m(1)}, \\
&\zeta \equiv \frac{\tilde{\zeta}-\zeta_m(1)}{1-\zeta_m(1)} + \zeta_m(\xi) \frac{1-\tilde{\zeta}}{1-\zeta_m(1)}.
\end{aligned}$$
With these variables, we can map the integration region to two rectangular regions:
$$\begin{aligned}
&(\tilde{\xi},\zeta) \in [\xi_m(1),1]\times[\zeta_m(1),1], \\
&(\xi,\tilde{\zeta}) \in [\xi_m(1),1]\times[\zeta_m(1),1].
\end{aligned}$$
The Jacobians are
$$\begin{aligned}
&\left\lVert\frac{\partial(\xi,\zeta)}{\partial(\tilde{\xi},\zeta)}\right\rVert = \frac{1-\xi_m(\zeta)}{1-\xi_m(1)}, \\
&\left\lVert\frac{\partial(\xi,\zeta)}{\partial(\xi,\tilde{\zeta})}\right\rVert = \frac{1-\zeta_m(\xi)}{1-\zeta_m(1)}.
\end{aligned}$$

## The subtraction trick

LDF and LFF are singular at momentum fraction goes to $1$.
To evaluate their convolution with the partonic cross section, we use the subtraction trick:
$$
\int_{x_\mathrm{min}}^1 f(x) H(x) dx = \int_{x_\mathrm{min}}^1 f(x) \big(H(x)-H(1)\big) dx + H(1) \int_{x_\mathrm{min}}^1 f(x) dx,
$$
where the integrand in the first term is suppressed at $x\rightarrow1$ since $H(x\rightarrow1)-H(1)\rightarrow0$.

It can be generalized to two-dimensional case as:
$$\begin{aligned}
\iint f(\xi)D(\zeta) H(\xi,\zeta) d\xi d\zeta
&= H(1,1) \left(\int f(\xi) d\xi\right) \left(\int D(\zeta) d\zeta\right) \\
&+ \left(\int f(\xi) \big(H(\xi,1)-H(1,1)\big) d\xi\right) \left(\int D(\zeta) d\zeta\right) \\
&+ \left(\int D(\zeta) \big(H(1,\zeta)-H(1,1)\big) d\zeta\right) \left(\int f(\xi) d\xi\right) \\
&+ \iint f(\xi) D(\zeta) \big(H(\xi,\zeta)-H(\xi,1)-H(1,\zeta)+H(1,1)\big) d\xi d\zeta.
\end{aligned}$$

### An improvement

Although we hope that the subtraction trick can suppress the integrand at $x\rightarrow1$, it actually can still be divergent.
It turns out that in most cases (not proved), we can suppress the integrand to $0$ at $x\rightarrow1$ by the following change of variable, in addition to the subtraction trick:
$$
\int_{x_\mathrm{min}}^1 F(x) dx = \int_0^1 F\left(x(2-x)+x_\mathrm{min}(1-x)^2\right) \cdot 2(1-x)(1-x_\mathrm{min}) dx.
$$

# Design of `RC_conv_xsec` with loose constraints

`RC_conv_xsec` is based on the 2D subtraction trick implemented in `qed_conv` in `QEDFactorization.jl`, but use a different method to calculate the regular 2D integral.
The idea is to transform the integration region to a retangular one so that `hcubature` can be applied to obtain high precision.

## The $x$- and $z$- constraints

For cross-section calculations, the integration region is restricted by the conditions
$$
\begin{aligned}
&\xi \geq \frac{1-y}{\zeta-xy} \quad\text{or equiv.}\quad \zeta \geq xy+\frac{1-y}{\xi}, \\
&\xi \geq zy+\frac{1-y}{\zeta} \quad\text{or equiv.}\quad \zeta \geq \frac{1-y}{\xi-zy}.
\end{aligned}
$$
We call the dirst line $x$-constraints and the second line $z$-constraints.
The minimun values of $\xi$ and $\zeta$ are given by
$$
\begin{aligned}
&\xi_{mx} = \frac{1-y}{1-xy},\quad \zeta_{mx} = 1-y+xy, \\
&\xi_{mz} = 1-y+zy,\quad \zeta_{mz} = \frac{1-y}{1-zy},
\end{aligned}
$$
where the subscripts $x$ and $z$ indicate the corresponding constraints.
We always have the following nice properties
$$
\xi_{mx} \leq \zeta_{mx} \quad\text{and}\quad \xi_{mz} \geq \zeta_{mz}.
$$

## Mapping to a rectangular region

The above formulas can be reversed to get
$$
xy = \frac{\zeta_{mx}-\xi_{mx}}{1-\xi_{mx}},\quad 1-y = \frac{1-\zeta_{mx}}{1-\xi_{mx}}\xi_{mx},
$$
from the $x$-constraints, and
$$
xy = \frac{\xi_{mz}-\zeta_{mz}}{1-\zeta_{mz}},\quad 1-y = \frac{1-\xi_{mz}}{1-\zeta_{mz}}\zeta_{mz}.
$$
from the $z$-constraints.
The $x$-constraint can now be written as
$$
\zeta \geq \frac{\zeta_{mx}-\xi_{mx}}{1-\xi_{mx}} + \frac{1-\zeta_{mx}}{1-\xi_{mx}}\frac{\xi_{mx}}{\xi},
$$
and the $z$-constraint can be written as
$$
\xi \geq \frac{\xi_{mz}-\zeta_{mz}}{1-\zeta_{mz}} + \frac{1-\xi_{mz}}{1-\zeta_{mz}}\frac{\zeta_{mz}}{\zeta}.
$$
We can now replace the $\zeta_{mx}$ in the $x$-constraint with a free variable $\tilde{\zeta}_m$, and replace the $\xi_{mz}$ in the $z$-constraint with a free variable $\tilde{\xi}_m$ to obtain a coordinate system on a rectangular region, that is,
$$
\left\{\begin{aligned}
&\xi \equiv \frac{\tilde{\xi}_m-\zeta_{mz}}{1-\zeta_{mz}} + \frac{1-\tilde{\xi}_m}{1-\zeta_{mz}}\frac{\zeta_{mz}}{\zeta} \\
&\zeta \equiv \frac{\tilde{\zeta}_m-\xi_{mx}}{1-\xi_{mx}} + \frac{1-\tilde{\zeta}_m}{1-\xi_{mx}}\frac{\xi_{mx}}{\xi}
\end{aligned}\right.,\quad
\left\{\begin{aligned}
&\tilde{\xi}_m \in [\max\{\xi_{mx},\xi_{mz}\},1] \\
&\tilde{\zeta}_m \in [\max\{\zeta_{mx},\zeta_{mz}\},1]
\end{aligned}\right..
\tag{$*$}
$$

## Transformation between the coordinate systems

We write $(*)$ in terms of auxiliary variables
$$
\left\{\begin{aligned}
&\xi \equiv A+\frac{C}{\zeta} \\
&\zeta \equiv B+\frac{D}{\xi}
\end{aligned}\right.,
$$
so that
$$
\left\{\begin{aligned}
&\xi = \frac{AB+C-D+\Delta}{2B} \\
&\zeta = \frac{AB-C+D+\Delta}{2A}
\end{aligned}\right.,\quad
\Delta \equiv \sqrt{A^2B^2+(C-D)^2+2AB(C+D)}.
$$

To calculate the Jacobian, we first differentiate $(*)$:
$$
\left\{\begin{aligned}
&d\xi = \frac{1}{\zeta}\left( \frac{\zeta-\zeta_{mz}}{1-\zeta_{mz}}d\tilde{\xi}_m - \frac{1-\tilde{\xi}_m}{1-\zeta_{mz}}\frac{\zeta_{mz}}{\zeta}d\zeta \right) \equiv \frac{1}{\zeta}\left( F~d\tilde{\xi}_m - C \frac{d\zeta}{\zeta} \right) \\
&d\zeta = \frac{1}{\xi}\left( \frac{\xi-\xi_{mx}}{1-\xi_{mx}}d\tilde{\zeta}_m - \frac{1-\tilde{\zeta}_m}{1-\xi_{mx}}\frac{\xi_{mx}}{\xi}d\xi \right) \equiv \frac{1}{\xi}\left( E~d\tilde{\zeta}_m - D \frac{d\xi}{\xi} \right)
\end{aligned}\right..
$$
Solving the equations gives
$$
\left\{\begin{aligned}
&d\xi = \frac{\xi^2\zeta F ~d\tilde{\xi}_m - \xi CE ~d\tilde{\zeta}_m}{\xi^2\zeta^2-CD} \\
&d\zeta = \frac{\xi\zeta^2 E ~d\tilde{\zeta}_m - \zeta DF ~d\tilde{\xi}_m}{\xi^2\zeta^2-CD}
\end{aligned}\right..
$$
The Jacobian is thus
$$
\left\lVert\frac{\partial(\xi,\zeta)}{\partial(\tilde{\xi}_m,\tilde{\zeta}_m)}\right\rVert
= \frac{\xi\zeta EF}{\xi^2\zeta^2-CD}.
$$

