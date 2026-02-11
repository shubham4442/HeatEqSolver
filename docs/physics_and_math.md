# Heat Equation: Physics and Mathematics

## Introduction

The heat equation is one of the most fundamental partial differential equations (PDEs) in mathematical physics. It describes how heat (or more generally, any diffusive quantity) spreads through a medium over time.

## Physical Background

### What is Heat Conduction?

Heat conduction (or thermal diffusion) is the transfer of thermal energy through a material without bulk motion of the material itself. At the microscopic level, heat is transferred through:

1. **In solids**: Lattice vibrations (phonons) and free electron movement
2. **In fluids**: Molecular collisions and energy exchange

### Fourier's Law

The foundation of heat conduction theory is **Fourier's Law** (1822), which states that the heat flux is proportional to the negative temperature gradient:

```
q = -k ∇T
```

where:
- `q` is the heat flux vector (W/m²)
- `k` is the thermal conductivity (W/(m·K))
- `T` is the temperature (K)
- `∇T` is the temperature gradient

The negative sign indicates that heat flows from hot to cold regions.

### Conservation of Energy

Combining Fourier's Law with the principle of energy conservation, we get:

```
ρcₚ ∂T/∂t = -∇·q + Q
```

where:
- `ρ` is the density (kg/m³)
- `cₚ` is the specific heat capacity (J/(kg·K))
- `Q` is the volumetric heat source (W/m³)

## The Heat Equation

### Derivation

Substituting Fourier's Law into the energy conservation equation:

```
ρcₚ ∂T/∂t = ∇·(k∇T) + Q
```

For a homogeneous material with constant thermal properties and no heat sources:

```
∂T/∂t = α ∇²T
```

where `α = k/(ρcₚ)` is the **thermal diffusivity** (m²/s).

### Mathematical Forms

#### 1D Heat Equation
```
∂u/∂t = α ∂²u/∂x²
```

#### 2D Heat Equation
```
∂u/∂t = α (∂²u/∂x² + ∂²u/∂y²)
```

#### 3D Heat Equation
```
∂u/∂t = α (∂²u/∂x² + ∂²u/∂y² + ∂²u/∂z²)
```

Or using the Laplacian operator:
```
∂u/∂t = α ∇²u
```

## Physical Interpretation

### Thermal Diffusivity (α)

The thermal diffusivity `α = k/(ρcₚ)` represents how quickly heat spreads through a material:

| Material       | α (m²/s)        |
|----------------|-----------------|
| Silver         | 1.66 × 10⁻⁴     |
| Copper         | 1.11 × 10⁻⁴     |
| Aluminum       | 9.7 × 10⁻⁵      |
| Steel          | 1.2 × 10⁻⁵      |
| Concrete       | 7.5 × 10⁻⁷      |
| Water          | 1.43 × 10⁻⁷     |
| Wood           | 8.2 × 10⁻⁸      |
| Air            | 2.2 × 10⁻⁵      |

- **High α**: Heat spreads quickly (metals)
- **Low α**: Heat spreads slowly (insulators)

### Time Scales

The characteristic diffusion time for heat to travel a distance `L` is:

```
τ = L² / α
```

Example: For a 1 cm aluminum rod:
```
τ = (0.01)² / (9.7 × 10⁻⁵) ≈ 1 second
```

## Boundary Conditions

### Dirichlet (First Kind)
Temperature is specified at the boundary:
```
u(x_boundary, t) = g(t)
```
**Physical meaning**: The boundary is in contact with a heat reservoir at a fixed temperature.

### Neumann (Second Kind)
Heat flux is specified at the boundary:
```
-k ∂u/∂n = q(t)
```
**Physical meaning**: A fixed rate of heat flow through the boundary.

Special case - **Insulated boundary** (q = 0):
```
∂u/∂n = 0
```

### Robin (Third Kind)
Convective heat transfer at the boundary:
```
-k ∂u/∂n = h(u - u_∞)
```
where `h` is the heat transfer coefficient and `u_∞` is the ambient temperature.

### Periodic
The solution wraps around:
```
u(0, t) = u(L, t)
```
**Physical meaning**: Ring-like or infinite periodic structures.

## Initial Conditions

The temperature distribution at t = 0 must be specified:
```
u(x, 0) = f(x)
```

Common initial conditions:
- **Uniform**: `f(x) = T₀`
- **Gaussian pulse**: `f(x) = A exp(-(x-x₀)²/σ²)`
- **Step function**: `f(x) = T₁ for x < L/2, T₂ otherwise`
- **Sinusoidal**: `f(x) = A sin(πx/L)`

## Analytical Solutions

### 1D with Dirichlet Boundaries

For `u(0,t) = u(L,t) = 0` and initial condition `u(x,0) = f(x)`:

```
u(x,t) = Σₙ Bₙ sin(nπx/L) exp(-α(nπ/L)²t)
```

where:
```
Bₙ = (2/L) ∫₀ᴸ f(x) sin(nπx/L) dx
```

### Gaussian Pulse Diffusion

For an initial Gaussian pulse in an infinite domain:
```
u(x,t) = (σ/√(σ² + 4αt)) exp(-(x-x₀)²/(σ² + 4αt))
```

Key observations:
- Peak amplitude decreases as `1/√t`
- Width increases as `√t`
- Total heat (integral) is conserved

## Well-Posedness

The heat equation is **well-posed** in the sense of Hadamard:

1. **Existence**: A solution exists for reasonable initial/boundary conditions
2. **Uniqueness**: The solution is unique (maximum principle)
3. **Stability**: Small changes in initial conditions produce small changes in the solution

### Maximum Principle

The maximum (or minimum) temperature in the domain occurs either:
- At the initial time (t = 0), or
- On the boundary

This means heat cannot spontaneously concentrate—it only diffuses.

## Related Equations

| Equation | Form | Application |
|----------|------|-------------|
| Heat equation | ∂u/∂t = α∇²u | Heat conduction |
| Diffusion equation | ∂c/∂t = D∇²c | Mass diffusion |
| Schrödinger equation | iℏ∂ψ/∂t = -ℏ²∇²ψ/2m | Quantum mechanics |
| Black-Scholes | ∂V/∂t = ½σ²S²∂²V/∂S² | Option pricing |

## References

1. Fourier, J. (1822). *Théorie analytique de la chaleur*
2. Carslaw, H.S. & Jaeger, J.C. (1959). *Conduction of Heat in Solids*
3. Incropera, F.P. et al. (2007). *Fundamentals of Heat and Mass Transfer*
4. Evans, L.C. (2010). *Partial Differential Equations*
