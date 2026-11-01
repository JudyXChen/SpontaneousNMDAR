# NMDAR magnesium-block kinetic calibration

## Scope

Only the reversible open-channel magnesium transition was recalibrated:

$$
O \mathop{\rightleftharpoons}^{k_B(V,M)}_{k_U(V)} B.
$$

Here, $O$ is conducting, $B$ is nonconducting, $V$ is in mV, and
$M=[\mathrm{Mg}^{2+}]_o$ is in mM. All other transition rates in the
ten-state receptor model were held fixed. The calibration was performed for
the stationary model at 100 µM glutamate.

This was an analytic constrained calibration, not a numerical least-squares
optimization. The Jahr--Stevens equation fixed the required ratio $k_B/k_U$.
Preserving the original product $k_Bk_U$ then selected one rate pair from the
infinitely many pairs with that ratio.

## 1. Original rates

The original block and unblock rates were

$$
k_B^{(0)}(V,M)=1200M\exp(-V/17)\quad\mathrm{s}^{-1},
$$

$$
k_U^{(0)}(V)=10800\exp(V/47)\quad\mathrm{s}^{-1}.
$$

Their product was therefore

$$
p_0(V,M)=k_B^{(0)}k_U^{(0)}
=12{,}960{,}000M
\exp\left[\left(-\frac{1}{17}+\frac{1}{47}\right)V\right].
$$

## 2. Stationary open probability of the ten-state model

At fixed glutamate, all states other than the $O\rightleftharpoons B$
transition can be represented by stationary weights. Taking the unblocked
state $C_2$ to have weight one, the unblocked weights at
$G=100\ \mu\mathrm{M}=10^{-4}\ \mathrm{M}$ are

$$
\begin{aligned}
w_{C_2}&=1,\\
w_{C_1}&=\frac{k_{C_2C_1}}{k_{C_1C_2}G}=0.022,\\
w_{C_0}&=w_{C_1}\frac{k_{C_1C_0}}{k_{C_0C_1}G}=0.000121,\\
w_D&=\frac{k_{C_2D}}{k_{DC_2}}=4.6666667,\\
w_O&=\frac{k_{C_2O}}{k_{OC_2}}=0.5076419.
\end{aligned}
$$

Thus,

$$
Z_U=w_{C_0}+w_{C_1}+w_{C_2}+w_D+w_O=6.1964296
$$

and the zero-magnesium stationary open probability is

$$
P_O(0)=\frac{w_O}{Z_U}=0.08192491.
$$

For every unit of stationary weight in $B$, the total weight in the blocked
branch $(B,C_{2B},D_B,C_{1B},C_{0B})$ is

$$
\begin{aligned}
F_B={}&1+\frac{k_{BC_{2B}}}{k_{C_{2B}B}}
\left[
1+\frac{k_{C_{2B}D_B}}{k_{D_BC_{2B}}}
+\frac{k_{C_{2B}C_{1B}}}{k_{C_{1B}C_{2B}}G}
\left(1+\frac{k_{C_{1B}C_{0B}}}{k_{C_{0B}C_{1B}}G}\right)
\right]\\
={}&34.66354764.
\end{aligned}
$$

In the model's parameter names, this calculation uses

$$
\frac{k_{BC_{2B}}}{k_{C_{2B}B}}=\frac{k_{0C2b}}{k_{C20b}}
=\frac{574}{97}.
$$

Because stationary balance across $O\rightleftharpoons B$ gives

$$
\frac{P_B}{P_O}=\frac{k_B}{k_U},
$$

the normalized stationary open probability is

$$
\frac{P_O(V,M)}{P_O(V,0)}
=\frac{1}{1+K\,k_B(V,M)/k_U(V)},
$$

where

$$
K=P_O(0)F_B=2.8398079368402547.
$$

The factor $K$ is necessary because blocking $O$ gives access to the
entire blocked branch, not only to state $B$. It depends on the fixed
glutamate concentration and the other receptor rates.

## 3. Ratio required by the Jahr--Stevens curve

The comparison target was the Jahr--Stevens Eq. 5 approximation

$$
g_{JS}(V,M)=
\frac{1}{1+(M/3.57)\exp(-0.062V)}.
$$

Equating the ten-state normalized open probability to this target gives

$$
K\frac{k_B}{k_U}=\frac{M}{3.57}\exp(-0.062V),
$$

so the required rate ratio is

$$
r_*(V,M)=\frac{k_B}{k_U}
=\frac{M}{3.57K}\exp(-0.062V).
$$

This equation fixes the ratio, but it does not fix the two rates separately.

## 4. Product-preserving constraint

The additional constraint was

$$
k_B(V,M)k_U(V)=p_0(V,M).
$$

Combining this product with the target ratio gives the rates directly:

$$
k_B(V,M)=\sqrt{p_0(V,M)r_*(V,M)},
$$

$$
k_U(V)=\sqrt{\frac{p_0(V,M)}{r_*(V,M)}}.
$$

The magnesium factors cancel from $k_U$, whereas $k_B$ remains linear in
$M$, as required for the pseudo-first-order association rate.

Substitution and simplification give

$$
\boxed{
k_B(V,M)=1130.6388608604M\exp(-V/20.0910256733)\ \mathrm{s}^{-1}
}
$$

and

$$
\boxed{
k_U(V)=11462.5460424539\exp(V/81.7893336063)\ \mathrm{s}^{-1}.
}
$$

The exponential scales follow from

$$
\begin{aligned}
\frac{1}{20.0910256733}
&=\frac{1}{2}\left[
0.062+\left(\frac{1}{17}-\frac{1}{47}\right)
\right],\\
\frac{1}{81.7893336063}
&=\frac{1}{2}\left[
0.062-\left(\frac{1}{17}-\frac{1}{47}\right)
\right].
\end{aligned}
$$

The prefactors follow from

$$
1130.6388608604
=\sqrt{\frac{1200\times10800}{3.57K}},
$$

$$
11462.5460424539
=\sqrt{(1200\times10800)(3.57K)}.
$$

## 5. Reproducible calculation

The following self-contained Python calculation reproduces $K$ and all four
calibrated constants. No optimizer is required.

```python
from math import sqrt

# Fixed glutamate concentration, M
G = 100e-6

# Unchanged ten-state rates from archive/NMDAR_old/model/evoked/homogenous/model.bngl
k_C0C1, k_C1C0 = 2e7, 11.0
k_C1C2, k_C2C1 = 1e7, 22.0
k_C20, k_0C2 = 93.0, 183.2
k_C20b, k_0C2b = 97.0, 574.0
k_C2D, k_DC2 = 16.8, 3.6

# Stationary unblocked weights relative to C2
w_C2 = 1.0
w_C1 = w_C2 * k_C2C1 / (k_C1C2 * G)
w_C0 = w_C1 * k_C1C0 / (k_C0C1 * G)
w_D = w_C2 * k_C2D / k_DC2
w_O = w_C2 * k_C20 / k_0C2
Z_U = w_C0 + w_C1 + w_C2 + w_D + w_O
p_O_zero_Mg = w_O / Z_U

# Total blocked-branch weight per unit weight in B
w_C2B_per_B = k_0C2b / k_C20b
w_DB_per_B = w_C2B_per_B * k_C2D / k_DC2
w_C1B_per_B = w_C2B_per_B * k_C2C1 / (k_C1C2 * G)
w_C0B_per_B = w_C1B_per_B * k_C1C0 / (k_C0C1 * G)
F_B = (1.0 + w_C2B_per_B + w_DB_per_B
       + w_C1B_per_B + w_C0B_per_B)
K = p_O_zero_Mg * F_B

# Original prefactors and voltage coefficients
b0_prefactor = 1200.0       # per mM per s
u0_prefactor = 10800.0      # per s
product_prefactor = b0_prefactor * u0_prefactor
product_voltage_slope = -1.0 / 17.0 + 1.0 / 47.0

# Jahr--Stevens target: exp(-target_slope * V)
target_slope = 0.062
target_ratio_prefactor = 1.0 / (3.57 * K)

# b = sqrt(product * ratio), u = sqrt(product / ratio)
b_prefactor = sqrt(product_prefactor * target_ratio_prefactor)
u_prefactor = sqrt(product_prefactor / target_ratio_prefactor)
b_voltage_coefficient = 0.5 * (target_slope - product_voltage_slope)
u_voltage_coefficient = 0.5 * (target_slope + product_voltage_slope)
b_voltage_scale = 1.0 / b_voltage_coefficient
u_voltage_scale = 1.0 / u_voltage_coefficient

print(f"K = {K:.16g}")
print(f"block prefactor = {b_prefactor:.16g} per mM per s")
print(f"block voltage scale = {b_voltage_scale:.16g} mV")
print(f"unblock prefactor = {u_prefactor:.16g} per s")
print(f"unblock voltage scale = {u_voltage_scale:.16g} mV")
```

Expected output:

```text
K = 2.839807936840254
block prefactor = 1130.638860860404 per mM per s
block voltage scale = 20.0910256732631 mV
unblock prefactor = 11462.54604245389 per s
unblock voltage scale = 81.78933360630565 mV
```

## 6. How the result was checked in the repository

The stored implementation is
`scripts/nmdar/compare_conductance.py`. It:

1. defines the calibrated rates;
2. constructs the full ten-state continuous-time Markov generator;
3. solves $Qp=0$, with $\sum_i p_i=1$, for the stationary occupancy;
4. calculates $P_O(V,M)/P_O(V,0)$;
5. compares that result with the Jahr--Stevens curve; and
6. checks that the MCell variable-rate tables implement the calibrated laws.

Across -100 to +60 mV, the numerical stationary solution agrees with the
target equation to machine precision (maximum absolute error below
$1.3\times10^{-15}$ for the checked 0.1 and 1.25 mM curves).

Because both sides have the same linear dependence on $M$ in their blocking
odds, 0.1 and 1.25 mM were verification curves rather than two independent
pieces of information in a regression. Within the Jahr--Stevens Eq. 5
approximation, the analytic equality holds for any positive $M$ at the
calibration glutamate concentration. The 0.1 mM curve is nevertheless an
extrapolation of that approximation outside its preferred high-magnesium
range.

## 7. Interpretation boundary

The procedure calibrates stationary normalized conductance. It does not fit
absolute open- or blocked-state dwell times, voltage-step relaxation, or
single-channel block/unblock events. Those transient measurements would be
needed to test or replace the product-preserving constraint. The subsequent
one-cluster MCell run tests the consequence of the calibrated rates; it was
not used to estimate them.
