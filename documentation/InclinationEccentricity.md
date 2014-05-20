Inclined and Eccentric Orbital Evolution (NOT IMPLEMENTED YET!!!) {#InclinationEccentricity}
=================================================================

The implementation of the evolution of inclined/eccentric orbits is based on
the formalism of Lai 2012, hence we will use the same notation.

The first thing is to derive the spherical harmonic expansion of the tidal
potential. 
\f[
	U(\mathbf{r}, t) = \frac{GM'}{r(t)}\left(
		1-\frac{r(t)}{\left|\mathbf{r}_{M'}\right|}\right)
\f]

where \f$r(t)\f$ is the distance between the centers of \f$M\f$ and \f$M'\f$
and \f$\left|\mathbf{r}_{M'}\right|\f$ is the distance between \f$M'\f$ and
the point \f$\mathbf{r}\f$ where the potential is being evaluated. Obviously
for a circular orbit \f$r(t)\f$ is a constant.

We start with the expansion in a coordinate system with
\f$\mathbf{\hat{z}}=\mathbf{\hat{L}}\f$, and
\f$\mathbf{\hat{y}}=\mathbf{\hat{S}}\times\mathbf{\hat{L}}\f$.

Let us define this expansion as:
\f[
	U(\mathbf{r}, t) = \frac{GM'}{r(t)}\sum_{l,m} c_{l,m}(\rho) 
		Y_{l,m}(\theta, \phi)
\f]
With \f$ \mathbf{r}=(\rho, \theta, \phi)\f$ in spherical coordinates.

Clearly:
\f[
	c_{l,m}(\rho, t) = \int_{0}^{2\pi} d\phi \int_{0}^{\pi} d\theta
		\tilde{U}(\mathbf{r}, t) Y_{l,m}^*(\theta, \phi)
\f]
where \f$\tilde{U}(\mathbf{r}, t)\equiv U(\mathbf{r}, t)r(t)/(GM')\f$.

The integrals are calulated the easiest by writing both
\f$U(\mathbf{r}, t)\f$ and \f$Y_{l,m}(\theta, \phi)\f$ in a coondinate
system where \f$\hat{\tilde{z}}\f$ points from \f$M\f$ to \f$M'\f$. Tilde
will be used to identify coordinates in that system.

\f[
	\tilde{U}(\mathbf{r}, t) = 1 - \frac{1}{\sqrt{(1-\tilde{z})^2 +
		\tilde{x}^2 + \tilde{y}^2}}
\f]
where \f$\tilde{x}\f$, \f$\tilde{y}\f$ and \f$\tilde{z}\f$ are the
coordinates of the point the potential is being evaluated scaled by
\f$r(t)\f$:

\f{eqnarray*}{
	\tilde{x}&=&\tilde{\rho}\sin\tilde{\theta}\cos\tilde{\phi}\\
	\tilde{y}&=&\tilde{\rho}\sin\tilde{\theta}\sin\tilde{\phi}\\
	\tilde{z}&=&\tilde{\rho}\cos\tilde{\theta}
\f}
with \f$\tilde{\rho}\f$ also scaled by \f$r(t)\f$.

With these we have:
\f[
	\tilde{U}(\tilde{\rho}, \tilde{\theta}, t) = 1 - \frac{1}{
		\sqrt{1-2\tilde{\rho}\cos\tilde{\theta} + \tilde{\rho}^2}}
\f]

As expected there is no \f$\tilde{\phi}\f$ dependence.

Now we need to write \f$Y_{l,m}\f$ in the tilde coordinate system: 
\f{eqnarray*}{
	x_{rot} &=& \tilde{z}\\
	y_{rot} &=& \tilde{y}\\
	z_{rot} &=& -\tilde{x}\\
	\rho=\rho_{rot} &=& \tilde{\rho}\\
	\theta = \theta_{rot} &\Rightarrow& \sin\theta = \sqrt{
		\sin^2\tilde{\theta}\sin^2\tilde{\phi} + \cos^2\tilde{\theta}},\quad
		\cos\theta = -\sin\tilde{\theta}\cos\tilde{\phi}\\
	\phi&=&\phi_{rot}-\Delta\phi(t)\\
	\sin\phi_{rot}&=& \frac{\sin\tilde{\theta}\sin\tilde{\phi}}
		{
			\sqrt{\sin^2\tilde{\theta}\sin^2\tilde{\phi}
			+ 
			\cos^2\tilde{\theta}}
		}\\
	\cos\phi_{rot}&=& \frac{\cos\tilde{\theta}}
		{
			\sqrt{\sin^2\tilde{\theta}\sin^2\tilde{\phi}
			+ 
			\cos^2\tilde{\theta}}
		}	
\f}
Where coordinates with "rot" subscript are in a coordinate system where
\f$\mathbf{\hat{z}_{rot}}=\mathbf{\hat{L}}\f$ and \f$\hat{x}_{rot}\f$ points
from the center of \f$M\f$ to the center of \f$M'\f$, and \f$\Delta\phi(t)\f$
is the angle between \f$\hat{x}_{rot}\f$ and \f$\hat{x}\f$
(\f$\Delta\phi(t)=\Omega t\f$ for a circular orbit with angular velocity
\f$\Omega\f$).

Note that:
\f{eqnarray*}{
	c_{l,m}(\rho, t) &=& \int_{\Delta\phi(t)}^{2\pi+\Delta\phi(t)}
		d\phi_{rot} \int_{0}^{\pi} d\theta \tilde{U}(\mathbf{r}, t)
		Y_{l,m}^*(\theta,\phi_{rot}-\Delta\phi(t))\\
		&=& exp(-im\Delta\phi(t))\int_{\Delta\phi(t)}^{2\pi+\Delta\phi(t)}
			d\phi_{rot} \int_{0}^{\pi} d\theta \tilde{U}(\mathbf{r}, t)
				Y_{l,m}^*(\theta, \phi_{rot})\\
		&=& \exp(-im\Delta\phi(t))\left\{
			\int_{\Delta\phi(t)}^{2\pi}  d\phi_{rot}
				\int_{0}^{\pi} d\theta \tilde{U}(\mathbf{r}, t)
				Y_{l,m}^*(\theta, \phi_{rot})
			+
			\int_{2\pi}^{2\pi+\Delta\phi(t)}  d\phi_{rot}
				\int_{0}^{\pi} d\theta \tilde{U}(\mathbf{r}, t)
				Y_{l,m}^*(\theta, \phi_{rot})\right\}\\
		&=& \exp(-im\Delta\phi(t))\left\{
			\int_{\Delta\phi(t)}^{2\pi}  d\phi_{rot}
				\int_{0}^{\pi} d\theta \tilde{U}(\mathbf{r}, t)
				Y_{l,m}^*(\theta, \phi_{rot})
			+
			\int_{0}^{\Delta\phi(t)}  d\phi_{rot}
				\int_{0}^{\pi} d\theta \tilde{U}(\mathbf{r}, t)
				Y_{l,m}^*(\theta, \phi_{rot}+2\pi)\right\}\\
		&=& \exp(-im\Delta\phi(t)) \int_{0}^{2\pi}
			d\phi_{rot} \int_{0}^{\pi} d\theta \tilde{U}(\mathbf{r}, t)
			Y_{l,m}^*(\theta, \phi_{rot})
\f}
We can then use Mathematica to evaluate:
\f[
	c_{0,0}=c_{1,0}=c_{2,\pm1}=0\quad,
	\quad c_{1,\pm 1}=\pm\sqrt{\frac{2\pi}{3}}\tilde{\rho}
		\exp(-im\Delta\phi(t))\quad,
	\quad c_{2,0}=\sqrt{\frac{\pi}{5}}\rho^2\exp(-im\Delta\phi(t))\quad,
	\quad c_{2,\pm 2}=-\sqrt{\frac{3\pi}{10}}\rho^2\exp(-im\Delta\phi(t))
\f]
The \f$c_{1,\pm 1}\f$ coefficients represent the gravitational acceleration
of the center of mass of \f$M\f$ due to \f$m\f$, and \f$c_{2,*}\f$ are the
lowest order tidal potential terms, so we have reproduced eq. (4) of Lai
2012.

Now we need to transform the \f$Y_{2,m}\f$ functions to a coordinate system where
the \f$z\f$ axis is along \f$\hat{S}\f$ and the \f$y\f$ axis is along
\f$\hat{S}\times\hat{L}\f$. We will use primes for the coordinates in the new
system. We have the following relations:
\f{eqnarray*}{
	y&=&y'=\rho'\sin\theta'\sin\phi'\\
	z&=&x'\sin\Theta + z'\cos\Theta
		=\rho'\sin\theta'\cos\phi'\sin\Theta + \rho'\cos\theta'\cos\Theta\\
	x&=&x'\cos\Theta - z'\sin\Theta
		=\rho'\sin\theta'\cos\phi'\cos\Theta - \rho'\cos\theta'\sin\Theta\\
	\rho&=&\rho'\\
	\cos\theta&=&z/\rho=\sin\theta'\cos\phi'\sin\Theta+\cos\theta'\cos\Theta\\
	\exp(i\phi)&=&\frac{x+iy}{\rho\sin\theta}
		=\frac{\sin\theta'\sin\phi'+i\left(\sin\theta'\cos\phi'\cos\Theta -
									\cos\theta'\sin\Theta\right)}
			{\sqrt{1-\cos^2\theta}}
\f}
Then using mathematica we show that the transformation between the prime and
non-prime coordinate system is indeed given by the Wigner D matrices quoted
in Lai 2012 (Equations 6-11).

The equivalent of Lai 2014 eq. 12 is then:
\f[
	U(\mathbf{r}, t)=-\sum_{m,m'} U_{m,m'}\rho^2 Y_{2,m}(\theta',\phi')
		\frac{\exp\left(-im'\Delta\phi(t)\right)a^3}{r^3(t)}
\f]
with
\f[
	U_{m,m'} \equiv \frac{GM'}{a^3}W_{2,m'}D_{m,m'}(\Theta)
\f]

This is where we diverge from Lai 2012 because we wish to consider elliptical
orbits. For general elliptical orbits it is convenient to define the origin
of time so that the planet is at periastron at \f$t=0\f$. Further we will take
\f$\Delta\phi(t)=\phi(t)+\phi_0\f$ where \f$\phi(0)=0\f$, which implies that
\f$\phi_0\f$ is the angle between periastron and \f$\hat{x}\f$ or
\f$90^\circ\f$ less than the angle between periastron and
\f$\hat{S}\times\hat{L}\f$. With these definitions:
\f{eqnarray*}{
	r(t)&=&a(1-e\cos u)\\
	\tan\left(\frac{\phi}{2}\right)&=&\sqrt{\frac{1+e}{1-e}}
		\tan\left(\frac{u}{2}\right)\\
	\Rightarrow
	1+\tan^2\left(\frac{\phi}{2}\right)&=&\frac{1-e\cos u}
		{(1-e)\cos^2\left(\frac{u}{2}\right)}\\
	1-\tan^2\left(\frac{\phi}{2}\right)&=&\frac{\cos u - e}
		{(1-e)\cos^2\left(\frac{u}{2}\right)}\\
	\Rightarrow
	\cos\phi &=& \frac{\cos u - e}{1 - e\cos u}\\
	\sin\phi &=& \frac{\sqrt{1-e^2}\sin u}{1 - e\cos u}
\f}
Where \f$u\f$ is the eccentric anomaly:
\f[
	u-e\sin u = \Omega t,\quad\Omega=\sqrt{\frac{G(M+M')}{a^3}}
\f]
Differentiating:
\f[
	dt=\frac{1-e\cos u}{\Omega}du
\f]

Since the orbital solution \f$r(t)\f$ and \f$\Delta \phi(t)\f$ is periodic
with a period of \f$2\pi/\Omega\f$, we can expand:
\f[
	\frac{a^3\exp(-im'\Delta \phi(t))}{r^3(t)}=\sum_s p_{m',s}
		\exp\left(-i s \Omega t\right)
\f]

Expressions for the \f$p_{m',s}(e)\f$ coefficients are derived 
(\ref InclinationEccentricity_pms1 "here") or 
(\ref InclinationEccentricity_pms2 "here").

Hence, our tidal potential can be written exactly as in Lai (2012), eq. 12,
except with \f$m'\f$ not limited to only 0 and 2:
\f[
	U(\mathbf{r}, t)=-\sum_{m,m'} U_{m,m'}\rho^2 Y_{2,m}(\theta',\phi')
		\exp\left(-im'\Omega t\right)
\f]
with
\f[
	U_{m,m'} \equiv \frac{GM'}{a^3} \mathcal{U}_{m,m'}\equiv
		\frac{GM'}{a^3} \sum_s W_{2,s}D_{m,s}(\Theta) p_{s,m'}
\f]
where we have switched the \f$m'\f$ and \f$s\f$ coefficients. Here is a (\ref
InclinationEccentricity_Ummtable "table") of \f$W_{2,m'}D_{m,m'}\f$.

From here we proceed following Lai (2012) again, but we have more than 6
independent timelags if the orbit is eccentric (for circular orbits,
\f$p_{m',s}=\delta_{m',s}\f$). 

The ansatz:
\f{eqnarray*}{
	\mathbf{\xi}_{m,s}(\mathbf{r},t)&=&
		\frac{U_{m,s}}{\omega_0^2}\mathbf{\bar{\xi}}_{m,s}(\mathbf{r})
		\exp(-is\Omega t + i\Delta_{m,s})\\
	\delta\rho(\mathbf{r},t)&=&\frac{U_{m,s}}{\omega_0^2}
		\delta\bar{\rho}_{m,s}(\mathbf{r}) 
		\exp(-is\Omega t + i\Delta_{m,s})\\
	\Delta_{m,s}&=&\tilde{\omega}_{m,s}t_{m,s}
\f}
with \f$\delta\bar{\rho}_{m,s}=-\nabla\cdot(\rho\mathbf{\bar{\xi}}_{m,s})\f$,
\f$\tilde{\omega}{m,s}=s\Omega-m\Omega_s\f$, where \f$\Omega_s\f$ is the
spin angular velocity of \f$M\f$, and \f$\omega_0\equiv\sqrt{GM/R^3}\f$ is the
dynamical frequency of \f$M\f$.

Here are the detailed devirations of the tidal (\ref
InclinationEccentricity_torque "torque") and (\ref
InclinationEccentricity_power "power").

As noted before, for general eccentric orbits, the number of timelags is not
only six, like in Lai (2012), but could be arbitrarily large, depending on
the precision required of the expansion and the value of the eccentricity.
Luckily, the fully general solution can be preserved, except we will need to
specify the dissipation as 3 timelags that are functions of frequency: one
for each \f$m=0,1,2\f$. The negative \f$m\f$ values can be handled by
inverting the sign of \f$m'\f$ and the forcing frequency.

Finally, we use the orbital energy and angular momentum:
\f{eqnarray*}{
	E_{orb}&=&-\frac{GMM'}{2a}\\
	L_{orb}&=&\frac{MM'}{M+M'}a^2\Omega\sqrt{1-e^2}=GMM'\sqrt{\frac{(1-e^2)MM'}{2E(M+M')}}
\f}
To derive the rate of change of the orbit and spin of \f$M\f$:
\f{eqnarray*}{
	\dot{S}&=&T_z\\
	\dot{a}&=&\frac{GMM'}{2E^2}\dot{E}\\
	\dot{\Theta}&=&-\frac{T_x}{S} - \frac{T_x\cos\Theta}{L} +
					\frac{T_z\sin\Theta}{L}\\
	\dot{e}&=&\frac{2(\dot{E}L+2E\dot{L})L(M+M')}{G(MM')^3}
\f}
