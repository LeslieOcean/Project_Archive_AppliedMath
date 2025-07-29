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

### Bifurcation Anlaysis by varying $h$ and find the third and fourth fixed points
We already have the prey isocline $v=\frac{\rho(k-u)(\delta+u)}{k\theta}$, which is a convex function and passes through $(0,\frac{\rho\delta}{\theta})$ and $(k,0)$. We can also derive the predator isocline $v = \frac{h(\delta+u)}{-a(\delta+u)+\theta_1u}-h$, which approaches $\infty$ as $u\rightarrow\frac{a\delta}{\theta_1-a}$ and passes through point $(k, \frac{h(\delta+k)}{-a(\delta+k)+\theta_1k}-h)$. Let $h = h_{bif}$ when $\Delta=0$. Varying $h$, we have phase portraits shown in following figures and we can observe a saddle-node bifurcation. 

<img width="300" height="230" alt="Q1_0fp" src="https://github.com/user-attachments/assets/0c927ead-add2-4ad1-a5a0-ff0af682ef04" />
<img width="300" height="230" alt="Q1_1fp" src="https://github.com/user-attachments/assets/8caee070-af8b-44a1-91ab-a507b1ae26ba" />
<img width="300" height="230" alt="Q1_2fp" src="https://github.com/user-attachments/assets/48893e10-aa8e-40a6-9e56-882bbed63379" />

According to the physical meaning of the model, at least we know the conversion efficiency $\theta_1$ should be within $[0, 1]$ normally and the equilibrium points $(u^\*, v^\*)$ that are biologically meaningful should be positive and satisfying $u^*\leq k$. Hence we assume that $\rho=1.76$, $k=5$, $\theta=1.33$, $\delta=1.5$, $a=0.02$ and $\theta_1=0.21$. With these parameters, we can determine that $h_{bif}\approx0.3324$.

### Stability Analysis
And two fixed points can be obtained when $h=0.24<h_{bif}$. The first one is $(u_1^\*, v_1^\*)=(3.8891, 1.5844)$ with corresponding eigenvalues $\lambda_1=-1.0729$, $\lambda_2=0.1003$. It is evident that this is a saddle point. The other one is $(u_2^\*, v_2^\*)= (1.3661, 2.7565)$, corresponding eigenvalues are complex $\lambda_{1.2}=0.1012\pm0.2574 i$ and $|arg(\lambda_{1,2})|=1.1961$. If the equilibrium point is asymptotically stable, it must satisfy that $|arg(\lambda_{1,2})|=1.1961>\frac{\beta\pi}{2}$. Thus, in this case, if $0<\beta<0.7615$, this fixed point is asymptotically stable, otherwise it is unstable.<br>

<img width="700" height="350" alt="Q1vectorfiled" src="https://github.com/user-attachments/assets/6ca2fbf0-193f-45f6-9cca-eaeb00212482" />


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
Here we use the Routh table and the Routh-Hurwitz Test mentioned in Wiggins' book to analyze the roots of this third-order characteristic polynomial. The associated Routh table is given by:
```math 
\begin{equation}
    \begin{array}{cc}
 1 & -\frac{5a}{6} + b \\
 \frac{a}{6} + 1 & \frac{ab}{6} \\
 \displaystyle\frac{\left( \frac{a}{6} + 1 \right)\left( -\frac{5a}{6} + b \right) - \frac{ab}{6}}{ \frac{a}{6} + 1 } &  \\
 \frac{ab}{6} & \\
\end{array}
\end{equation}
```
Based on inequalities from Routh-Hurwitz test, we visualize the stable region in $(a,b)$ space in the figure:
<img width="350" height="265" alt="Q4stableregion" src="https://github.com/user-attachments/assets/3f80c760-c3a9-459e-8ffc-e2cebb13dc29" />

If $a$ and $b$ use the values falling out of the green region, the equilibrium point is an unstable point or saddle point, which is necessary for chaos. To better understand the type of the fixed point, we can also visualize the region that there are one real eigenvalue and a complex conjugate pair in $(a,b)$ space, where discriminant $\Delta<0$:

<img width="350" height="265" alt="Q4complexroot" src="https://github.com/user-attachments/assets/5e1e15a6-bf91-4fdb-8336-45a26f749a31" />

### Occurance of Bifurcation
Root locus analysis - Hopf Bifurcation when $a\approx7.478$

<img width="350" height="265" alt="Q4rootlocus" src="https://github.com/user-attachments/assets/397c49e8-89d5-48ed-a081-9835ee2b1d82" />

Appearance of a stable limit cycle and an unstable limit cycle when $a=6.63$ - Subcritical Hopf Bifurcation but stablized by the negative higher order term like $-x^3$

<img width="350" height="265" alt="Q46_63" src="https://github.com/user-attachments/assets/9f8c93a4-bc86-41cf-bd2f-2c251c6e9364" />
