# Design of `collgrid`

The integral in collinear factorized structure function depends only on $x_B$, $z_h$, and $q_T^2/Q^2$, so we can store it as a function of these variables. The ranges of these variables are restricted by
$$
\frac{q_T^2}{Q^2} \leq \frac{1-x_B}{x_B} \frac{1-z_h}{z_h} \tag{$*$}.
$$
We will call this the master relation. Due to this relation, it is not suitable to store structure functions on a Cartesian $\{x_B\}\times\{z_h\}\times\{q_T^2/Q^2\}$ grid. Here we use a polar-like coordinate for the $x_B$-$z_h$ plane.

## Bounds of variables

Let's determine some bounds of the variables. We make use of the relations of QED distorted variables:
$$
x_B \leq \hat{x}_B, \quad z_h \leq \hat{z}_h,
$$
which implies that we can set lower bounds for $x_B$ and $z_h$:
$$
x_{B,\text{low}} \leq x_B, \quad z_{h,\text{low}} \leq z_h.
$$
With these lower bounds, the master relation $(*)$ gives an upper bound of $q_T^2/Q^2$:
$$
\frac{q_T^2}{Q^2} \leq \frac{1-x_{B,\text{low}}}{x_{B,\text{low}}} \frac{1-z_{h,\text{low}}}{z_{h,\text{low}}} \equiv \left(\frac{q_T^2}{Q^2}\right)_\text{up}.
$$
We also know that collinear factorization is not applicable at small $q_T$, and we can set a lower bound for $q_T^2/Q^2$:
$$
\left(\frac{q_T^2}{Q^2}\right)_\text{low} \leq \frac{q_T^2}{Q^2}.
$$
This lower bound we choose should of course satisfies
$$
\left(\frac{q_T^2}{Q^2}\right)_\text{low} < \left(\frac{q_T^2}{Q^2}\right)_\text{up}. \tag{$1$}
$$
With the lower bound of $q_T^2/Q^2$, the master relation $(*)$ also gives upper bounds of $x_B$ and $z_h$:
$$
x_B \leq \left( 1 + \frac{z_{h,\text{low}}}{1-z_{h,\text{low}}} \left(\frac{q_T^2}{Q^2}\right)_\text{low} \right)^{-1} \equiv x_{B,\text{up}},
$$
$$
z_h \leq \left( 1 + \frac{x_{B,\text{low}}}{1-x_{B,\text{low}}} \left(\frac{q_T^2}{Q^2}\right)_\text{low} \right)^{-1} \equiv z_{h,\text{up}}.
$$
It follows from $(1)$ that
$$
x_{B,\text{low}} < x_{B,\text{up}}, \quad z_{h,\text{low}} < z_{h,\text{up}}.
$$

## $q_T^2/Q^2$ grid

We first fix the grid of $q_T^2/Q^2$.

$q_T^2/Q^2$ usually won't get as large as $(q_T^2/Q^2)_\text{up}$. So we can choose another bound $(q_T^2/Q^2)_{\text{up}'}$ that satisfies
$$
\left(\frac{q_T^2}{Q^2}\right)_\text{low} < \left(\frac{q_T^2}{Q^2}\right)_{\text{up}'} < \left(\frac{q_T^2}{Q^2}\right)_\text{up}.
$$
The interval $[(q_T^2/Q^2)_\text{low},(q_T^2/Q^2)_{\text{up}'}]$ should cover the range of interest. We then fix the grid of $q_T^2/Q^2$ on this interval and use `logrange` for this grid.

## Radial grid on $x_B$-$z_h$ plane

We define
$$
\left(\frac{q_T^2}{Q^2}\right)_\text{max} \equiv \frac{1-x_B}{x_B} \frac{1-z_h}{z_h}.
$$
Given a value of $q_T^2/Q^2$, $(q_T^2/Q^2)_\text{max}$ should satisfy
$$
\frac{q_T^2}{Q^2} \leq \left(\frac{q_T^2}{Q^2}\right)_\text{max} \leq \left(\frac{q_T^2}{Q^2}\right)_\text{up}.
$$
We define a new variable $R$ so that
$$
\left(\frac{q_T^2}{Q^2}\right)_\text{max} \equiv \left(\frac{q_T^2}{Q^2}\right)^R \left(\frac{q_T^2}{Q^2}\right)_\text{up}^{1-R}.
$$
We fix the grid of $R$ on the interval $[0,1]$ and use a linear `range` for this grid.

Note that at $R=1$, or equivalently $q_T^2/Q^2=(q_T^2/Q^2)_\text{max}$, structure functions should be exactly $0$. At $R=0$, or equivalently $(q_T^2/Q^2)_\text{max}=(q_T^2/Q^2)_\text{up}$, $x_B$ is constantly $x_{B,\text{low}}$, and $z_h$ is constantly $z_{h,\text{low}}$.

## Angular grid on $x_B$-$z_h$ plane

We define a new variable
$$
D \equiv \left. \frac{1-x_B}{x_B} \middle/ \frac{1-z_h}{z_h} \right..
$$
Given a value of $(q_T^2/Q^2)_\text{max}$, thee bounds of $D$ are given by
$$
D \geq \left. \left(\frac{q_T^2}{Q^2}\right)_\text{max} \middle/ \left(\frac{1-z_{h,\text{low}}}{z_{h,\text{low}}}\right)^2 \right. \equiv D_\text{min}
$$
$$
D \leq \left. \left(\frac{1-x_{B,\text{low}}}{x_{B,\text{low}}}\right)^2 \middle/ \left(\frac{q_T^2}{Q^2}\right)_\text{max} \right. \equiv D_\text{max}
$$

We define another variable $A$ so that
$$
D \equiv \left(D_\text{min}\right)^{1-A} \left(D_\text{max}\right)^A.
$$
We fix the grid of $A$ on the interval $[0,1]$ and use a linear `range` for this grid.

## Summary

In summary, we first manually choose four bounds:
$$
x_{B,\text{low}}, \quad z_{h,\text{low}}, \quad \left(\frac{q_T^2}{Q^2}\right)_\text{low}, \quad \left(\frac{q_T^2}{Q^2}\right)_{\text{up}'}.
$$
We then calculate the bound
$$
\left(\frac{q_T^2}{Q^2}\right)_{\text{up}} = \frac{1-x_{B,\text{low}}}{x_{B,\text{low}}} \frac{1-z_{h,\text{low}}}{z_{h,\text{low}}}.
$$
We next
- fix a `logrange` grid for $q_T^2/Q^2$ on the interval $[(q_T^2/Q^2)_\text{low},(q_T^2/Q^2)_{\text{up}'}]$,
- fix a linear `range` grid for $R$ on the interval $[0,1]$,
- and fix a linear `range` grid for $A$ on the interval $[0,1]$.

$x_B$ and $z_h$ are related to $R$ and $A$ as
$$
x_B = \left( 1 + \frac{z_{h,\text{low}}}{1-z_{h,\text{low}}} \left(\frac{q_T^2}{Q^2}\right)_\text{max}^{1-A} \left(\frac{q_T^2}{Q^2}\right)_\text{up}^A \right)^{-1},
$$
$$
z_h = \left( 1 + \frac{x_{B,\text{low}}}{1-x_{B,\text{low}}} \left(\frac{q_T^2}{Q^2}\right)_\text{max}^A \left(\frac{q_T^2}{Q^2}\right)_\text{up}^{1-A} \right)^{-1},
$$
where
$$
\left(\frac{q_T^2}{Q^2}\right)_\text{max} \equiv \left(\frac{q_T^2}{Q^2}\right)^R \left(\frac{q_T^2}{Q^2}\right)_\text{up}^{1-R}.
$$
$R$ and $A$ are related to $x_B$ and $z_h$ as
$$
R = \frac{\log\left[ \left(\frac{q_T^2}{Q^2}\right)_\text{up} \middle/ \left(\frac{q_T^2}{Q^2}\right)_\text{max} \right]}{\log\left[ \left(\frac{q_T^2}{Q^2}\right)_\text{up} \middle/ \left(\frac{q_T^2}{Q^2}\right) \right]},
$$
$$
A = \frac{\log\left[ \left(\frac{q_T^2}{Q^2}\right)_\text{up} \middle/ \left( \frac{1-x_{B,\text{low}}}{x_{B,\text{low}}} \frac{1-z_h}{z_h} \right) \right]}{\log\left[ \left(\frac{q_T^2}{Q^2}\right)_\text{up} \middle/ \left(\frac{q_T^2}{Q^2}\right)_\text{max} \right]}.
$$
where
$$
\left(\frac{q_T^2}{Q^2}\right)_\text{max} = \frac{1-x_B}{x_B} \frac{1-z_h}{z_h}.
$$

## Example bounds

For example, we can choose
$$
x_{B,\text{low}} = \texttt{1.0000E-03}, \quad z_{h,\text{low}} = \texttt{1.0000E-01},
$$
$$
\left(\frac{q_T^2}{Q^2}\right)_\text{low} = \texttt{1.0000E-04}, \quad \left(\frac{q_T^2}{Q^2}\right)_{\text{up}'} = \texttt{1.0000E+02},
$$
so that
$$
\left(\frac{q_T^2}{Q^2}\right)_{\text{up}} \approx \texttt{8.9910E+03} \quad\text{(round down)}.
$$
