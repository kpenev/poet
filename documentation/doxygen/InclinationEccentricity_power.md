Tidal Power for Inclined and Eccentric Orbits {#InclinationEccentricity_power}
=============================================
Starting from eq. 21 of Lai (2012)
\f{eqnarray*}{
	\dot{E}&=&-\int d^3 x \rho(\mathbf{r}) 
		\frac{\partial \mathbf{\xi}(\mathbf{r}, t)}{\partial t}
		\cdot\nabla U*(\mathbf{r}, t)\\
	&=& \int d^3 x \nabla\cdot
		\left(\rho(\mathbf{r})
			\frac{\partial \mathbf{\xi}(\mathbf{r}, t)}{\partial t}
		\right)U*(\mathbf{r},t)
\f}

Where the surface term is zero because the density is zero at the surface.
\f[
	\dot{E}=-\left(\frac{GM'}{\omega_0 a^3}\right)^2\Omega
		\sum_{m,m',\mu,\mu'} \mathcal{U}_{m,m'}\mathcal{U}_{\mu,\mu'}im'
		\exp(i(\mu'-m')\Omega t + i\Delta_{m,m'}) \int d^3 x 
			\delta\bar{\rho}_{m,m'}(\mathbf{r}) r^2 Y_{2,\mu}^*(\theta, \phi)
\f]
Since \f$\delta\bar{\rho}_{m,m'}(\mathbf{r})\propto\exp{im\phi}\f$ we must
have \f$m=\mu\f$, further, averaging over an orbit we must have the time
dependence term in the exponent vanish, \f$\Rightarrow \mu'=m'\f$:
\f{eqnarray*}{
	\dot{E}&=&-\left(\frac{GM'}{\omega_0 a^3}\right)^2\Omega \sum_{m,m'}
		\mathcal{U}_{m,m'}^2 im'
		\exp(i\Delta_{m,m'}) \int d^3 x 
			\delta\bar{\rho}_{m,m'}(\mathbf{r}) r^2 Y_{2,m}^*(\theta, \phi)\\
		&=&-T_0\Omega \sum_{m,m'}
		\mathcal{U}_{m,m'}^2im'\exp(i\Delta_{m,m'}) \kappa_{m,m'}\\
\f}
So taking the real part:
\f[
	\dot{E}=T_0\Omega \sum_{m,m'}
		\mathcal{U}_{m,m'}^2m'\sin(\Delta_{m,m'}) \kappa_{m,m'}\\
\f]
