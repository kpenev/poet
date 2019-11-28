Tidal Power for Inclined and Eccentric Orbits
*********************************************

Starting from eq. 21 of Lai (2012)

.. math::

	\dot{E}&=&-\int d^3 x \rho(\mathbf{r}) 
		\frac{\partial \mathbf{\xi}(\mathbf{r}, t)}{\partial t}
		\cdot\nabla U*(\mathbf{r}, t)\\
	&=& \int d^3 x \nabla\cdot
		\left(\rho(\mathbf{r})
			\frac{\partial \mathbf{\xi}(\mathbf{r}, t)}{\partial t}
		\right)U*(\mathbf{r},t)

Where the surface term is zero because the density is zero at the surface.

.. math::

	\dot{E}=-\left(\frac{GM'}{\omega_0 a^3}\right)^2\Omega
		\sum_{m,m',\mu,\mu'} \mathcal{U}_{m,m'}\mathcal{U}_{\mu,\mu'}im'
		\exp(i(\mu'-m')\Omega t + i\Delta_{m,m'}) \int d^3 x 
			\delta\bar{\rho}_{m,m'}(\mathbf{r}) r^2 Y_{2,\mu}^*(\theta, \phi)

Since :math:`\delta\bar{\rho}_{m,m'}(\mathbf{r})\propto\exp{im\phi}` we must
have :math:`m=\mu`, further, averaging over an orbit we must have the time
dependence term in the exponent vanish, :math:`\Rightarrow \mu'=m'`:

.. math::

	\dot{E}&=&-\left(\frac{GM'}{\omega_0 a^3}\right)^2\Omega \sum_{m,m'}
		\mathcal{U}_{m,m'}^2 im'
		\exp(i\Delta_{m,m'}) \int d^3 x 
			\delta\bar{\rho}_{m,m'}(\mathbf{r}) r^2 Y_{2,m}^*(\theta, \phi)\\
		&=&-T_0\Omega \sum_{m,m'}
		\mathcal{U}_{m,m'}^2im'\exp(i\Delta_{m,m'}) \kappa_{m,m'}\\

So taking the real part:

.. math::

	\dot{E}=T_0\Omega \sum_{m,m'}
		\mathcal{U}_{m,m'}^2m'\sin(\Delta_{m,m'}) \kappa_{m,m'}\\
