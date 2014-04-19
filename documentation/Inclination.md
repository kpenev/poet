Inclined and Eccentric Orbital Evolution (NOT IMPLEMENTED YET!!!)
=================================================================

The implementation of the evolution of inclined/eccentric orbits is based on
the formalism of Lai 2012, hence we will use the same notation.

The first thing is to derive the spherical harmonic expansion of the tidal
potential. 
\f[
	U(\mathbf{r}, t) = \frac{GM'}{a(t)}\left(
		1-\frac{a(t)}{\left|\mathbf{r}_{M'}\right|}\right)
\f]

where \f$a(t)\f$ is the distance between the centers of \f$M\f$ and \f$M'\f$
and \f$\left|\mathbf{r}_{M'}\right|\f$ is the distance between \f$M'\f$ and
the point \f$\mathbf{r}\f$ where the potential is being evaluated. Obviously
for a circular orbit \f$a(t)\f$ is a constant.

We start with the expansion in a coordinate system with
\f$\mathbf{\hat{z}}=\mathbf{\hat{L}}\f$, and
\f$\mathbf{\hat{y}}=\mathbf{\hat{S}}\times\mathbf{\hat{L}}\f$.

Let us define this expansion as:
\f[
	U(\mathbf{r}, t) = \frac{GM'}{a(t)}\sum_{l,m} c_{l,m}(\rho) 
		Y_{l,m}(\theta, \phi)
\f]
With \f$ \mathbf{r}=(\rho, \theta, \phi)\f$ in spherical coordinates.

Clearly:
\f[
	c_{l,m}(\rho, t) = \int_{0}^{2\pi} d\phi \int_{0}^{\pi} d\theta
		\tilde{U}(\mathbf{r}, t) Y_{l,m}^*(\theta, \phi)
\f]
where \f$\tilde{U}(\mathbf{r}, t)\equiv U(\mathbf{r}, t)a(t)/(GM')\f$.

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
\f$a(t)\f$:

\f{eqnarray*}{
	\tilde{x}&=&\tilde{\rho}\sin\tilde{\theta}\cos\tilde{\phi}\\
	\tilde{y}&=&\tilde{\rho}\sin\tilde{\theta}\sin\tilde{\phi}\\
	\tilde{z}&=&\tilde{\rho}\cos\tilde{\theta}
\f}
with \f$\tilde{\rho}\f$ also scaled by \f$a(t)\f$.

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
from the center f \f$M\f$ to the center of \f$M'\f$, and \f$\Delta\phi(t)\f$
is the angle between \f$\hat{x}_{rot}\f$ and \f$\hat{x}\f$
(\f$\Delta\phi(t)=\Omega t\f$ for a circular orbit with angular velocity
\f$\Omega\f$).

Note that:
\f{eqnarray*}{
	c_{l,m}(\rho, t) &=& \int_{\Delta\phi(t)}^{2\pi+\Delta\phi(t)}
		d\phi_{rot} \int_{0}^{\pi} d\theta \tilde{U}(\mathbf{r}, t)
		Y_{l,m}^*(\theta,\phi_{rot}-\Delta\phi(t))\\
		&=& e^{-im\Delta\phi(t)}\int_{\Delta\phi(t)}^{2\pi+\Delta\phi(t)}
			d\phi_{rot} \int_{0}^{\pi} d\theta \tilde{U}(\mathbf{r}, t)
				Y_{l,m}^*(\theta, \phi_{rot})\\
		&=& e^{-im\Delta\phi(t)}\left\{
			\int_{\Delta\phi(t)}^{2\pi}  d\phi_{rot}
				\int_{0}^{\pi} d\theta \tilde{U}(\mathbf{r}, t)
				Y_{l,m}^*(\theta, \phi_{rot})
			+
			\int_{2\pi}^{2\pi+\Delta\phi(t)}  d\phi_{rot}
				\int_{0}^{\pi} d\theta \tilde{U}(\mathbf{r}, t)
				Y_{l,m}^*(\theta, \phi_{rot})\right\}\\
		&=& e^{-im\Delta\phi(t)}\left\{
			\int_{\Delta\phi(t)}^{2\pi}  d\phi_{rot}
				\int_{0}^{\pi} d\theta \tilde{U}(\mathbf{r}, t)
				Y_{l,m}^*(\theta, \phi_{rot})
			+
			\int_{0}^{\Delta\phi(t)}  d\phi_{rot}
				\int_{0}^{\pi} d\theta \tilde{U}(\mathbf{r}, t)
				Y_{l,m}^*(\theta, \phi_{rot}+2\pi)\right\}\\
		&=& e^{im\Delta\phi(t)} \int_{0}^{2\pi}
			d\phi_{rot} \int_{0}^{\pi} d\theta \tilde{U}(\mathbf{r}, t)
			Y_{l,m}^*(\theta, \phi_{rot})
\f}
We can then use Mathematica to evaluate:
\f[
	c_{0,0}=c_{1,0}=c_{2,\pm1}=0\quad,
	\quad c_{1,\pm 1}=\pm\sqrt{\frac{2\pi}{3}}\tilde{\rho}\quad,
	\quad c_{2,0}=\sqrt{\frac{\pi}{5}}\rho^2\quad,
	\quad c_{2,\pm 2}=-\sqrt{\frac{3\pi}{10}}\rho^2
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
	e^{i\phi}&=&\frac{x+iy}{\rho\sin\theta}
		=\frac{\sin\theta'\sin\phi'+i\left(\sin\theta'\cos\phi'\cos\Theta -
									\cos\theta'\sin\Theta\right)}
			{\sqrt{1-\cos^2\theta}}
\f}
Then using mathematica we show that the transformation between the prime and
non-prime coordinate system is indeed given by the Wigner D matrices quoted
in Lai 2012 (Equations 6-10).
