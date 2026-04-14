# 2D Wave Packet on a Rectangle

This script solves the **2D scalar wave equation** on a rectangular domain with **fixed-wall (Dirichlet) boundary conditions** using a **double sine-series expansion**.

The initial condition is a **Gaussian-modulated carrier wave packet** with an initial velocity chosen so that the packet begins traveling in a prescribed direction.

The code is written to run in both **MATLAB** and **GNU Octave**.

---

## Equation solved

The script solves

\[
u_{tt} = c^2 (u_{xx} + u_{yy})
\]

on the rectangle

\[
0 \le x \le L_x, \qquad 0 \le y \le L_y
\]

with boundary conditions

\[
u = 0
\]

on all four walls.

These are fixed-end or hard-wall boundary conditions.

---

## Method

The solution is expanded in the sine basis

\[
u(x,y,t) = \sum_{m=1}^{M}\sum_{n=1}^{N} C_{mn}(t)
\sin\left(\frac{m\pi x}{L_x}\right)
\sin\left(\frac{n\pi y}{L_y}\right)
\]

which automatically satisfies the homogeneous Dirichlet boundary conditions.

The script:

1. Defines the initial displacement `u0` as a Gaussian-envelope wave packet.
2. Computes the initial velocity `v0` so the packet initially propagates in a chosen direction.
3. Projects `u0` and `v0` onto the sine basis using numerical quadrature.
4. Evolves each mode analytically in time.
5. Reconstructs the field on a display grid and animates it.

This is a standard **spectral separation-of-variables** solution with numerical projection of the initial data.

---

## File

- `wave2d_rectangle_packet.m`

---

## Features

- 2D wave equation on a rectangle
- Fixed-wall boundaries
- Traveling Gaussian wave packet initial condition
- Adjustable propagation angle and carrier wavelength
- Spectral solution using sine modes
- Animation of the time evolution
- Compatible with MATLAB and Octave

---

## Main parameters

### Domain and wave speed

```matlab
Lx = 1.0;   % domain length in x
Ly = 0.8;   % domain length in y
c  = 1.0;   % wave speed
```

### Mode truncation

```matlab
M = 40;     % number of sine modes in x
N = 40;     % number of sine modes in y
```

Higher values improve representation of narrow packets or short wavelengths, but increase cost.

### Quadrature grid for projection

```matlab
Nxq = 220;
Nyq = 180;
```

This grid is used only to project the initial conditions onto the sine basis.

### Packet parameters

```matlab
A0 = 1.0;
x0 = 0.28 * Lx;
y0 = 0.35 * Ly;

sigx = 0.08 * Lx;
sigy = 0.08 * Ly;

theta   = pi/6;
lambda0 = 0.12;
k0      = 2*pi/lambda0;
```

These control the amplitude, initial location, spatial width, propagation direction, and carrier wavelength of the packet.

### Time and animation

```matlab
Tfinal  = 1.6;
nFrames = 220;
```

---

## Initial conditions

The initial displacement is a Gaussian-envelope carrier wave:

\[
u(x,y,0) = \text{env}(x,y)\cos(\phi(x,y))
\]

where the Gaussian envelope is centered at `(x0, y0)` and the carrier wavevector is determined by `theta` and `lambda0`.

The initial velocity is chosen as

\[
u_t(x,y,0) = -c\, \hat d \cdot \nabla u_0
\]

so the packet initially moves in the direction

\[
\hat d = (\cos\theta,\sin\theta)
\]

This produces a clean traveling-packet launch, subject to the constraints of the bounded domain and the reflected-wave dynamics.

---

## How to run

### MATLAB

Save the script as:

```text
wave2d_rectangle_packet.m
```

and run:

```matlab
wave2d_rectangle_packet
```

### GNU Octave

From the Octave prompt:

```octave
wave2d_rectangle_packet
```

or from the shell:

```bash
octave --gui wave2d_rectangle_packet.m
```

---

## Output

The script opens a figure window and animates the wave field over time using `imagesc`.

The color scale is chosen from the initial amplitude estimate, and the title displays the current simulation time.

---

## What to expect physically

Because the boundaries are fixed, the wave packet reflects from the walls with sign change consistent with Dirichlet conditions. Over time, the packet disperses into the modal content supported by the rectangle, and the evolution becomes a superposition of reflected and interfering standing-wave components.

This is not a free-space propagation model. It is a bounded-domain modal solution.

---

## Notes on accuracy

There are three main numerical resolution controls:

1. **Mode truncation**: `M`, `N`
2. **Projection grid**: `Nxq`, `Nyq`
3. **Display grid**: `Nx`, `Ny`

If the packet is too narrow or the carrier wavelength too short, increase:

- `M`, `N`
- `Nxq`, `Nyq`

Otherwise the initial condition may be under-resolved in the modal projection.

---

## Typical modifications

You can easily adapt the script to:

- change the domain size
- change the wave speed
- make the packet narrower or wider
- rotate the propagation direction
- change the carrier wavelength
- increase the number of modes
- save frames or export a movie
- replace the initial condition with another function

---

## Limitations

- Boundary conditions are homogeneous Dirichlet only.
- The geometry is a rectangle only.
- The solution is based on a truncated modal expansion.
- Very narrow or high-frequency packets require more modes and finer projection grids.
- No damping, forcing, or nonlinear terms are included.

---

## Possible extensions

Natural next steps include:

- Neumann or mixed boundary conditions
- external forcing
- damping
- multiple packets
- obstacle masks or irregular geometries using other numerical methods
- movie export to GIF or MP4
- energy diagnostics
- comparison with finite-difference or finite-element solutions

---

## License / reuse

Adapt freely for teaching, experimentation, and numerical PDE demonstrations.
