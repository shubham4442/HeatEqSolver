# Numerical Methods: Stability, Convergence, and Constraints

## Overview

This document explains the numerical methods used in the heat equation solver, with particular focus on stability conditions, convergence properties, and practical constraints.

## Finite Difference Method

### Discretization

We discretize the continuous domain into a grid:

```
Space: xᵢ = x₀ + i·Δx,  i = 0, 1, ..., Nₓ-1
Time:  tⁿ = n·Δt,       n = 0, 1, 2, ...
```

The solution is approximated at grid points: `uᵢⁿ ≈ u(xᵢ, tⁿ)`

### Second Derivative Approximation

The spatial second derivative is approximated using the **central difference**:

```
∂²u/∂x² ≈ (uᵢ₊₁ - 2uᵢ + uᵢ₋₁) / Δx²
```

This approximation has **second-order accuracy**: O(Δx²)

### Time Derivative Approximation

**Forward difference** (explicit method):
```
∂u/∂t ≈ (uᵢⁿ⁺¹ - uᵢⁿ) / Δt
```

This has **first-order accuracy**: O(Δt)

## Explicit FTCS Scheme

### The Scheme

FTCS = Forward Time, Central Space

Combining the approximations:

```
(uᵢⁿ⁺¹ - uᵢⁿ)/Δt = α(uᵢ₊₁ⁿ - 2uᵢⁿ + uᵢ₋₁ⁿ)/Δx²
```

Solving for the new time level:

```
uᵢⁿ⁺¹ = uᵢⁿ + r(uᵢ₊₁ⁿ - 2uᵢⁿ + uᵢ₋₁ⁿ)
```

where `r = αΔt/Δx²` is the **mesh ratio** or **Fourier number**.

### Multi-dimensional Extensions

**2D:**
```
uᵢⱼⁿ⁺¹ = uᵢⱼⁿ + rₓ(uᵢ₊₁,ⱼ - 2uᵢⱼ + uᵢ₋₁,ⱼ) + rᵧ(uᵢ,ⱼ₊₁ - 2uᵢⱼ + uᵢ,ⱼ₋₁)
```

**3D:**
```
uᵢⱼₖⁿ⁺¹ = uᵢⱼₖⁿ + rₓ(Δ²ₓu) + rᵧ(Δ²ᵧu) + r_z(Δ²_zu)
```

where `Δ²` denotes the second difference operator in each direction.

## Stability Analysis

### Von Neumann Stability Analysis

Assume a Fourier mode solution:
```
uᵢⁿ = ξⁿ exp(ikxᵢ)
```

Substituting into the FTCS scheme:
```
ξ = 1 + r(exp(ikΔx) - 2 + exp(-ikΔx))
  = 1 + 2r(cos(kΔx) - 1)
  = 1 - 4r·sin²(kΔx/2)
```

For stability, we need `|ξ| ≤ 1` for all wave numbers k.

### Stability Conditions

The **amplification factor** `ξ` satisfies:
- Maximum value: ξ = 1 (when k = 0)
- Minimum value: ξ = 1 - 4r (when kΔx = π)

For stability: `-1 ≤ 1 - 4r ≤ 1`

This gives the **stability criterion**:

| Dimension | Stability Condition | Physical Meaning |
|-----------|---------------------|------------------|
| 1D | r ≤ 1/2 | αΔt/Δx² ≤ 0.5 |
| 2D | r ≤ 1/4 | αΔt(1/Δx² + 1/Δy²) ≤ 0.5 |
| 3D | r ≤ 1/6 | αΔt(1/Δx² + 1/Δy² + 1/Δz²) ≤ 0.5 |

### What Happens When Unstable?

If r > 1/2 (in 1D), the amplification factor |ξ| > 1 for high-frequency modes:
- Errors grow exponentially in time
- Solution exhibits oscillations
- Eventually blows up to infinity

Example of instability (r = 0.6):
```
Step 0: [0, 0, 1, 0, 0]
Step 1: [0, 0.6, -0.2, 0.6, 0]      # Oscillations appear
Step 2: [0.36, -0.56, 1.08, -0.56, 0.36]  # Growing!
Step 3: ...                          # Explosion
```

## Convergence

### Lax Equivalence Theorem

For a well-posed linear initial value problem, a consistent finite difference scheme is convergent **if and only if** it is stable.

### Consistency

A scheme is consistent if the truncation error → 0 as Δt, Δx → 0.

For FTCS, the truncation error is:
```
τ = O(Δt) + O(Δx²)
```

The scheme is **first-order in time, second-order in space**.

### Convergence Rate

The global error satisfies:
```
|u(xᵢ, tⁿ) - uᵢⁿ| ≤ C(Δt + Δx²)
```

where C depends on the smoothness of the exact solution.

### Richardson Extrapolation

To verify convergence, run with grid spacings h and h/2:
```
Order p ≈ log₂(E_h / E_{h/2})
```

Expected: p ≈ 2 (dominated by spatial error when Δt is small)

## Practical Constraints

### Time Step Selection

Given spatial resolution Δx and thermal diffusivity α, the maximum stable time step is:

**1D:**
```
Δt_max = 0.5 × Δx² / α
```

**2D:**
```
Δt_max = 0.25 × min(Δx², Δy²) / α
```

**3D:**
```
Δt_max = (1/6) × min(Δx², Δy², Δz²) / α
```

### Safety Factor

In practice, use a **safety factor** of 0.9 or less:
```
Δt = safety × Δt_max
```

This accounts for:
- Numerical roundoff
- Non-uniform grids
- Variable coefficients

### Grid Resolution

Choose Δx based on:

1. **Feature resolution**: Δx << characteristic length of initial condition
   ```
   Δx ≤ σ/5  (for Gaussian with width σ)
   ```

2. **Accuracy requirements**: Smaller Δx → more accuracy
   ```
   Error ∝ Δx²
   ```

3. **Computational cost**: More points → more computation
   ```
   Cost ∝ N × (T/Δt)
   ```

### Memory Requirements

| Dimension | Memory (doubles) |
|-----------|------------------|
| 1D | 2 × Nₓ |
| 2D | 2 × Nₓ × Nᵧ |
| 3D | 2 × Nₓ × Nᵧ × N_z |

Example: 3D with N = 100 in each direction:
```
Memory = 2 × 100³ × 8 bytes = 16 MB
```

### Computational Cost per Time Step

| Dimension | Operations |
|-----------|------------|
| 1D | O(Nₓ) |
| 2D | O(Nₓ × Nᵧ) |
| 3D | O(Nₓ × Nᵧ × N_z) |

Total cost for simulation to time T:
```
Cost = O(N^d × T/Δt) = O(N^d × T × N² / α) = O(T × N^{d+2} / α)
```

## Error Sources

### Truncation Error
- Inherent in the discretization
- Controlled by grid refinement
- Order: O(Δt) + O(Δx²)

### Roundoff Error
- Due to finite precision arithmetic
- Accumulates over many time steps
- Typically ~10⁻¹⁴ per operation (double precision)

### Initial Condition Error
- Discretization of continuous initial data
- Aliasing of high-frequency components

### Boundary Condition Error
- Numerical implementation of BCs
- Ghost point interpolation

## Verification Strategies

### Method of Manufactured Solutions (MMS)

1. Choose an exact solution: `u_exact(x,t) = sin(πx)exp(-απ²t)`
2. Verify it satisfies the PDE
3. Run solver and compare to exact solution
4. Check convergence rates

### Conservation Tests

For periodic BC, total heat should be conserved:
```
∫ u(x,t) dx = constant
```

Numerical test:
```
|Σᵢ uᵢⁿ - Σᵢ uᵢ⁰| < tolerance
```

### Symmetry Tests

If initial condition is symmetric, solution should remain symmetric:
```
u(x,t) = u(L-x,t)  for all t
```

## Common Pitfalls

### 1. Ignoring Stability
**Symptom**: Solution explodes
**Fix**: Check r ≤ r_max before running

### 2. Insufficient Resolution
**Symptom**: Solution is too smooth, features missing
**Fix**: Increase N or decrease Δx

### 3. Wrong Boundary Conditions
**Symptom**: Solution behaves unexpectedly near boundaries
**Fix**: Verify BC implementation, check for off-by-one errors

### 4. Accumulating Roundoff
**Symptom**: Solution drifts over very long simulations
**Fix**: Use higher precision, periodic restarts

### 5. CFL Condition Confusion
**Symptom**: Stability works for some α but not others
**Fix**: Remember r depends on α, recalculate Δt for each case

## Summary Table

| Parameter | Symbol | Constraint | Effect of Violation |
|-----------|--------|------------|---------------------|
| Mesh ratio | r | r ≤ 1/(2d) | Instability, explosion |
| Time step | Δt | Δt ≤ Δx²/(2αd) | Instability |
| Grid spacing | Δx | Δx << L_feature | Poor resolution |
| Total time | T | T > 0 | Invalid simulation |
| Diffusivity | α | α > 0 | Ill-posed problem |

where d = dimension (1, 2, or 3).

## References

1. Strikwerda, J.C. (2004). *Finite Difference Schemes and PDEs*
2. LeVeque, R.J. (2007). *Finite Difference Methods for Ordinary and Partial Differential Equations*
3. Thomas, J.W. (1995). *Numerical Partial Differential Equations: Finite Difference Methods*
4. Morton, K.W. & Mayers, D.F. (2005). *Numerical Solution of Partial Differential Equations*
