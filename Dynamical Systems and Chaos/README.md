# Dynamical Systems and Chaos coursework

##  Stability of a Fractional Lotka-Volterra Model
A vector field is given as
```math
\begin{equation}
    \left\{
\begin{array}{l}
    {}^{C}_{0}D_t^{\beta} u(t) = \rho u(t) \left( 1 - \dfrac{u(t)}{k} \right) - \dfrac{\theta u(t) v(t)}{\delta + u(t)}, \\
    {}^{C}_{0}D_t^{\beta} v(t) = \dfrac{\theta_1 u(t) v(t)}{\delta + u(t)} - a v(t) - \dfrac{h v(t)}{h + v(t)},
\end{array}
\right.
\end{equation}
```
where $u(t), v(t) \geq $ are prey and predator population densities, $0 < \beta \leq 1$ is the order of the Caputo fractional differential operator 
${}^{C}_{0}D_t^{\beta}(\cdot)$, and $\rho, k, \theta, \theta_1, \delta, a, h$ are positive parameters that determine the dynamics of the system.<br>





## Analysis of a Chaotic Continuous-time Dynamical System
The dynamical system needed to be analyzed is
```math
\begin{equation}
    \begin{cases}
        \dot{x} = a(y-\frac{1}{16}x^3-\frac{1}{6}x),\\
        \dot{y} = x-y+z,\\
        \dot{z} = -by,
    \end{cases}
\end{equation}
```
where $a$ and $b$ are parameters.<br>
