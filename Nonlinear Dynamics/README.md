# Synchronization of Coupled Biological Oscillators

## Kuramoto Model
```math
\begin{equation}
    \theta_i'=\omega_i+\frac{k}{n} \sum\sin{(\theta_j-\theta_i)}, \ i= 1, ..., n
\end{equation}
```
where $\theta_i$ is the phase of the oscillator, $\omega_i$ is its natural frequency, $k$ is the coupling strength of the oscillator and $n$ is the total number of oscillators.<br>
Phase difference and frequency mismatch can be defined as:
```math
\begin{equation}
  \phi=\theta_1-\theta_2
\end{equation}
```
```math
\begin{equation}
  \Delta\omega=\omega_1-\omega_2  
\end{equation} 
```
 The phase difference for the system is: 
```math
\begin{equation}
  \phi'=\Delta\omega-k \ sin(\phi)    
\end{equation}
```
When $\Delta\omega > k$, there are no fixed points and the phase difference keeps increasing. The slowest rate of phase difference occurs at $\phi = \frac{\pi}{2}$. As increasing $k$ till it is equal to $\Delta \omega$,  $\phi = \frac{\pi}{2}$ becomes the semi-stable fixed point. If $\Delta\omega < k$, the two fixed points appear and $\phi = arcsin(\frac{\Delta\omega}{k})$ is the stable one and $\phi = \pi - arcsin(\frac{\Delta\omega}{k})$ is the unstable one. Results for negative $\Delta\omega$ are similar.<br>
<img width="250" height="200" alt="omegasmallerk" src="https://github.com/user-attachments/assets/fc3a42e0-dbdf-4324-8c41-87dab9bf57f1" />
<img width="250" height="200"  alt="omegaequalk" src="https://github.com/user-attachments/assets/f2b7ceea-2855-4fe0-b894-84f88db87f55" />
<img width="250" height="200"  alt="omegalargerk" src="https://github.com/user-attachments/assets/c721ce76-574c-4b74-95ed-81ba4ecacce0" /><br>
It's evident that a saddle-node bifurcation occurs at $\phi = \frac{\pi}{2}$ (or  $\phi = \pm\frac{\pi}{2}n$ for $n\in\mathbb{N}$).

## Two Oscillators
<!---### Phase locked analysis $\Delta \omega$ = 0.5 rad/s-->


<!---### Frequency locked analysis $k=0.25 \ rad / s^{-1}$-->
Arnold Tougue plot shows the synchronization region ($\phi'=0$) <br>
<img width="400" height="250" alt="arnold" src="https://github.com/user-attachments/assets/d898c310-8dab-4f98-9e36-c9faf9f5ea42" />

## Three Oscillators
Two possible configurations for three coupled oscillators <br>
<img width="465" height="150" alt="Configurations" src="https://github.com/user-attachments/assets/b306c1d4-7431-4a97-9f59-9e22ddc81379" /><br>
In the first case ($A$), each oscillator is connected and influenced by the other two. Assuming the coupling strength $k$ is the same for each couple, the equations of motion are:
```math
\begin{equation}
    \theta_1'=\omega_1+\frac{k}{3}(\sin{(\theta_2-\theta_1)}+\sin{(\theta_3-\theta_1)}) 
\end{equation}
```
```math
\begin{equation}
    \theta_2'=\omega_2+\frac{k}{3}(\sin{(\theta_3-\theta_2)}+\sin{(\theta_1-\theta_2)}) 
\end{equation}
```
```math
\begin{equation}
    \theta_3'=\omega_3+\frac{k}{3}(\sin{(\theta_2-\theta_3)}+\sin{(\theta_1-\theta_3)}) 
\end{equation}
```
Finding the equations for phase differences:
```math
\begin{equation}
 \phi'_{12} =\Delta\omega_{12}+\frac{k}{3}(-2\cdot\sin(\phi_{12})+\sin{(\phi_{23})}-\sin{(\phi_{12}+\phi_{23})}) 
 \end{equation}
```
```math
\begin{equation}
    \phi'_{23}=\Delta\omega_{23}+\frac{k}{3}(-2\cdot\sin(\phi_{23})+\sin{(\phi_{12})}-\sin{(\phi_{12}+\phi_{23})}) 
\end{equation}
```
For identical oscillators $\Delta\omega_{12}=\Delta\omega_{23}=0$, five fixed points can be identified in the space $[0,\pi]\times[0,\pi]$.
<img width="350" height="250" alt="configA" src="https://github.com/user-attachments/assets/b7783852-9ce5-482a-b5bd-3a030bad9c9a" /><br>

In the second configuration ($B$), the oscillators are connected in series; there is no direct influence of oscillator 1 on oscillator 3. 
```math
\begin{equation}
    \theta_1'=\omega_1+\frac{k}{2}(\sin{(\theta_2-\theta_1)})
\end{equation}
```
```math
\begin{equation}
    \theta_2'=\omega_2+\frac{k}{3}(\sin{(\theta_3-\theta_2)}+\sin{(\theta_1-\theta_2)}) 
\end{equation}
```
```math
\begin{equation}
    \theta_3'=\omega_3+\frac{k}{2}(\sin{(\theta_2-\theta_3)})
\end{equation}
```
Finding the equations for phase differences for configuration $B$::
```math
\begin{equation}
  \phi'_{12}= \Delta\omega_{12}+ \frac{k}{6}(-5\cdot \sin{\phi_{12} + 2\cdot\sin{\phi_{23}}})
 \end{equation}
```
```math
\begin{equation}
  \phi'_{23}= \Delta\omega_{23}+ \frac{k}{6}(-5\cdot \sin{\phi_{23} + 2\cdot\sin{\phi_{12}}})
 \end{equation}
```
From the phase plane $\phi_{23} \ vs \ \phi_{12}$, four fixed points are identified.<br>
<img width="350" height="250" alt="configB" src="https://github.com/user-attachments/assets/901936ba-2c2a-4f9c-aaff-2324fb4b76e5" /><br>

Arnold Tongue plots for both configurations:


<img width="350" height="250" alt="configB_arnold" src="https://github.com/user-attachments/assets/130e8548-77fb-4db1-9260-05f6b122635b" />
<img width="350" height="250" alt="configA_arnold" src="https://github.com/user-attachments/assets/5d64e95f-3f1c-4193-878c-41b97e01ca7b" />


## Multiple Oscillators
