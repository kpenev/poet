**********************************************
Tidal Torque for Inclined and Eccentric Orbits
**********************************************

Starting from eq. 20 of Lai (2012):

.. math::

	\mathbf{T}&=&-\int d^3 x \delta\rho(\mathbf{r}, t) \mathbf{r}\times
		\nabla U*(\mathbf{r}, t)\\
		&=&-\left(\frac{G M'}{a^3 \omega_0}\right)^2
			\int d^3 x 
			\left\{
				\sum_{m,m'}\mathcal{U}_{m,m'}
					\delta\bar{\rho}_{m,m'}(\mathbf{r})
					\exp(-im'\Omega t + i\Delta_{m,m'})
			\right\}
			\mathbf{r}\times
			\left\{
				\sum_{\mu,\mu'}\mathcal{U}_{\mu,\mu'}
				\exp(i \mu'\Omega t) 
				\nabla \left[r^2 Y_{2,\mu}^*(\theta, \phi)\right]
			\right\}

Since :math:`\mathbf{r}\times\nabla r^2=0`:

.. math::

	\mathbf{T}&=&-\left(\frac{G M'}{a^3 \omega_0}\right)^2
			\int d^3 x 
			\left\{
				\sum_{m,m'}\mathcal{U}_{m,m'}
					\delta\bar{\rho}_{m,m'}(\mathbf{r})
					\exp(-im'\Omega t + i\Delta_{m,m'})
			\right\}
			\mathbf{r}\times
			\left\{
				\sum_{\mu,\mu'}\mathcal{U}_{\mu,\mu'}
				\exp(i \mu'\Omega t) r^2 
				\nabla Y_{2,\mu}^*(\theta, \phi)
			\right\}\\
		&=& -\left(\frac{G M'}{a^3 \omega_0}\right)^2
			\sum_{m,m',\mu,\mu'} \mathcal{U}_{m,m'} \mathcal{U}_{\mu,\mu'}
				\exp(i(\mu'-m')\Omega t + i\Delta_{m,m'})
				\int d^3 x 
					\delta\bar{\rho}_{m,m'}(\mathbf{r})
					r^2 \mathbf{r}\times
					\nabla Y_{2,\mu}^*(\theta, \phi)

The torque in the z direction
=============================

Taking the dot product of the torque with :math:`\hat{z}=\cos\theta \hat{r} -
\sin\theta \hat{\theta}`, and noting that the :math:`\hat{r}` part of
:math:`\hat{z}` does not contribute, since at any point it is orthogonal to
:math:`\mathbf{r}\times \nabla Y_{2,\mu}^*(\theta, \phi)`, we get:

.. math::

	T_z&=& -\left(\frac{G M'}{a^3 \omega_0}\right)^2
			\sum_{m,m',\mu,\mu'} \mathcal{U}_{m,m'} \mathcal{U}_{\mu,\mu'}
				\exp(i(\mu'-m')\Omega t + i\Delta_{m,m'})
				\int d^3 x 
					\delta\bar{\rho}_{m,m'}(\mathbf{r})
					r^2 \sin\theta \hat{\theta}\cdot \left[
						\mathbf{r}\times
						\nabla Y_{2,\mu}^*(\theta, \phi)
					\right]\\
	&=& \left(\frac{G M'}{a^3 \omega_0}\right)^2
			\sum_{m,m',\mu,\mu'} \mathcal{U}_{m,m'} \mathcal{U}_{\mu,\mu'}
				\exp(i(\mu'-m')\Omega t + i\Delta_{m,m'})
				\int d^3 x 
					\delta\bar{\rho}_{m,m'}(\mathbf{r})
					r^2 \sin \theta (\mathbf{r}\times\hat{\theta})\cdot
					\nabla Y_{2,\mu}^*(\theta, \phi)\\
	&=& -\left(\frac{G M'}{a^3 \omega_0}\right)^2
			\sum_{m,m',\mu,\mu'} \mathcal{U}_{m,m'} \mathcal{U}_{\mu,\mu'}
				\exp(i(\mu'-m')\Omega t + i\Delta_{m,m'})
				\int d^3 x 
					\delta\bar{\rho}_{m,m'}(\mathbf{r})
					r^3 \sin\theta
					\hat{\phi}\cdot\nabla Y_{2,\mu}^*(\theta, \phi)\\
	&=& -\left(\frac{G M'}{a^3 \omega_0}\right)^2
			\sum_{m,m',\mu,\mu'} \mathcal{U}_{m,m'} \mathcal{U}_{\mu,\mu'}
				\exp(i(\mu'-m')\Omega t + i\Delta_{m,m'})
				\int d^3 x 
					\delta\bar{\rho}_{m,m'}(\mathbf{r})
					r^2
					\frac{\partial Y_{2,\mu}^*(\theta, \phi)}{\partial\phi}\\
	&=& -\left(\frac{G M'}{a^3 \omega_0}\right)^2
			\sum_{m,m',\mu,\mu'} i\mu\mathcal{U}_{m,m'} \mathcal{U}_{\mu,\mu'}
				\exp(i(\mu'-m')\Omega t + i\Delta_{m,m'})
				\int d^3 x 
					\delta\bar{\rho}_{m,m'}(\mathbf{r})
					r^2 Y_{2,\mu}^*(\theta, \phi)

Using the fact that
:math:`\delta\bar{\rho}_{m,m'}(\mathbf{r})\propto\exp(im\phi)`, we must have
:math:`m=\mu`, further if we average over over an orbit we must have
:math:`m'=\mu'`:

.. math::

	T_z= -\left(\frac{G M'}{a^3 \omega_0}\right)^2
			\sum_{m,m'} i m \mathcal{U}_{m,m'}^2
				\exp(i\Delta_{m,m'})
				\int d^3 x \delta\bar{\rho}_{m,m'}(\mathbf{r}) r^2 
						Y_{2,m}^*(\theta, \phi)

The real part of which is:

.. math::

	T_z&=&\frac{G R^3}{M}\left(\frac{M'}{a^3}\right)^2
			\sum_{m,m'} m \mathcal{U}_{m,m'}^2
				\sin(\Delta_{m,m'})
				\int d^3 x \delta\bar{\rho}_{m,m'}(\mathbf{r}) r^2 
						Y_{2,m}^*(\theta, \phi)\\
	T_z&=&T_0 \sum_{m,m'} \mathcal{U}_{m,m'}^2 m \kappa_{m,m'}
			\sin(\Delta_{m,m'})\\
	T_0&\equiv& G R^5\left(\frac{M'}{a^3}\right)^2\\
	\kappa&\equiv&\frac{1}{MR^2}
			\int d^3 x \delta\bar{\rho}_{m,m'}(\mathbf{r}) r^2 
						Y_{2,m}^*(\theta, \phi)

The torque in the x direction
=============================

Now we dot with :math:`\hat{x}=\sin\theta\cos\phi\hat{r} +
\cos\theta\cos\phi\hat{\theta} - \sin\phi \hat{\phi}`. Again, the
:math:`\hat{r}` term does not contribute:

.. math::

	T_x &=& -\left(\frac{G M'}{a^3 \omega_0}\right)^2
			\sum_{m,m',\mu,\mu'} \mathcal{U}_{m,m'} \mathcal{U}_{\mu,\mu'}
				\exp(i(\mu'-m')\Omega t + i\Delta_{m,m'})
				\int d^3 x 
					\delta\bar{\rho}_{m,m'}(\mathbf{r})
					r^2 (\cos\theta\cos\phi\hat{\theta} -
						 \sin\phi\hat{\phi}) \left[
						\mathbf{r}\times
						\nabla Y_{2,\mu}^*(\theta, \phi)
					\right]\\
	&=& -\left(\frac{G M'}{a^3 \omega_0}\right)^2
			\sum_{m,m',\mu,\mu'} \mathcal{U}_{m,m'} \mathcal{U}_{\mu,\mu'}
				\exp(i(\mu'-m')\Omega t + i\Delta_{m,m'})
				\int d^3 x 
					\delta\bar{\rho}_{m,m'}(\mathbf{r})
					r^3 (\cos\theta\cos\phi\hat{\phi} +
						 \sin\phi\hat{\theta})\cdot
						\nabla Y_{2,\mu}^*(\theta, \phi)\\
	&=& -\left(\frac{G M'}{a^3 \omega_0}\right)^2
			\sum_{m,m',\mu,\mu'} \mathcal{U}_{m,m'} \mathcal{U}_{\mu,\mu'}
				\exp(i(\mu'-m')\Omega t + i\Delta_{m,m'})
				\int d^3 x 
					\delta\bar{\rho}_{m,m'}(\mathbf{r})
					r^2 \left(
						\cot\theta\cos\phi\frac{\partial}{\partial \phi} +
						 \sin\phi\frac{\partial}{\partial \theta}
					\right) Y_{2,\mu}^*(\theta, \phi)

Averaging over an orbit:

.. math::

	T_x &=& \left(\frac{G M'}{a^3 \omega_0}\right)^2
			\sum_{m,m',\mu} \mathcal{U}_{m,m'} \mathcal{U}_{\mu,m'}
				\exp(i\Delta_{m,m'})
				\int d^3 x 
					\delta\bar{\rho}_{m,m'}(\mathbf{r})
					r^2 \left(
						i\mu\cot\theta\cos\phi Y_{2,\mu}^*(\theta, \phi) -
						 \sin\phi\frac{\partial Y_{2,\mu}^*(\theta, \phi)}
									  {\partial \theta}
					\right)\\

Now we use:

.. math::

	\frac{\partial Y_{2,\mu}^*}{\partial \theta}=\mu\cot\theta Y_{2,\mu}^* +
	\sqrt{(2-\mu)(3+\mu)}\exp(i\phi) Y_{2,\mu+1}^*

To get:

.. math::

	T_x &=& \frac{T_0}{MR^2}
			\sum_{m,m',\mu} \mathcal{U}_{m,m'} \mathcal{U}_{\mu,m'}
				i\exp(i\Delta_{m,m'})
				\int d^3 x 
					\delta\bar{\rho}_{m,m'}(\mathbf{r})
					r^2 \left(
						\mu\cot\theta\exp(i\phi) Y_{2,\mu}^*(\theta, \phi) +
						 i\sin\phi\sqrt{(2-\mu)(3+\mu)}\exp(i\phi) Y_{2,\mu+1}^*
					\right)\\
		&=& \frac{T_0}{MR^2}
			\sum_{m,m',\mu} \mathcal{U}_{m,m'} \mathcal{U}_{\mu,m'}
				i\exp(i\Delta_{m,m'})
				\int d^3 x 
					\delta\bar{\rho}_{m,m'}(\mathbf{r})
					r^2 \left(
						\mu\cot\theta\exp(i\phi) Y_{2,\mu}^*(\theta, \phi) +
						 \sqrt{(2-\mu)(3+\mu)}[\exp(2i\phi)-1] Y_{2,\mu+1}^*
					\right)\\

Since :math:`\delta\bar{\rho}_{m,m'}(\mathbf{r})\propto\exp(im\phi)` the real
part of the above expression is:

.. math::

	T_x &=& -\frac{T_0}{MR^2}
			\sum_{m,m',\mu} \mathcal{U}_{m,m'} \mathcal{U}_{\mu,m'}
				\sin(\Delta_{m,m'})
				\int d^3 x 
					\delta\bar{\rho}_{m,m'}(\mathbf{r})
					r^2 \left(
						\mu\cot\theta\exp(i\phi) Y_{2,\mu}^*(\theta, \phi) +
						 \frac{\sqrt{(2-\mu)(3+\mu)}}{2}[\exp(2i\phi)-1] Y_{2,\mu+1}^*
					\right)\\
		&=& T_0
			\sum_{m,m'} \mathcal{U}_{m,m'} \sin(\Delta_{m,m'})(
				\kappa_{m,m'}^-\mathcal{U}_{m-1,m'}+
				\kappa_{m,m'}^+\mathcal{U}_{m+1,m'})\\
	\kappa_{m,m'}^-&\equiv& \frac{\sqrt{(3-m)(2+m)}}{2MR^2}
		\int d^3 x \delta\bar{\rho}_{m,m'}(\mathbf{r})
		r^2 Y_{2,m}^*=\frac{\sqrt{(3-m)(2+m)}}{2}\kappa_{m,m'}\\
	\kappa_{m,m'}^+&\equiv&-\frac{1}{MR^2} \int d^3 x 
		\delta\bar{\rho}_{m,m'}(\mathbf{r}) r^2\left[
			(m+1)\cot\theta\exp(i\phi) Y_{2,m+1}^*(\theta,\phi)+
			\frac{\sqrt{(1-m)(4+m)}}{2}\exp(2i\phi) Y_{2,m+2}^*
		\right]

We have already expressed :math:`\kappa_{m,m'}^-` in terms of
:math:`\kappa_{m,m'}`, and just as in Lai (2012), we only need
:math:`\kappa_{m,m'}^+` for :math:`m=0, \pm 1, 2`. 

.. math::

	\kappa_{0,m'}^+&=&-\frac{1}{MR^2} \int d^3 x 
		\delta\bar{\rho}_{0,m'}(\mathbf{r}) r^2\left[
			\cot\theta\exp(i\phi) Y_{2,1}^*(\theta,\phi)+
			\exp(2i\phi) Y_{2,2}^*(\theta,\phi)
		\right]\\
	&=&-\frac{1}{MR^2} \int d^3 x 
		\delta\bar{\rho}_{0,m'}(\mathbf{r}) r^2\left[
			-\frac{1}{2}\sqrt{\frac{15}{2\pi}}\cos^2\theta+
			\frac{1}{4}\sqrt{\frac{15}{2\pi}}\sin^2\theta
		\right]\\
	&=&\frac{1}{MR^2} \int d^3 x 
		\delta\bar{\rho}_{0,m'}(\mathbf{r}) r^2
		\frac{1}{4}\sqrt{\frac{15}{2\pi}}(3\cos^2\theta-1)\\
	&=&\frac{\sqrt{3/2}}{MR^2}\int d^3 x 
		\delta\bar{\rho}_{0,m'}(\mathbf{r}) r^2 Y_{2,0}^*(\theta, \phi)\\
	\Rightarrow \kappa_{0,m'}^+&=&\sqrt{3/2}\kappa_{0,m'}\\
	\kappa_{-1,m'}^+&=&-\frac{1}{MR^2} \int d^3 x 
		\delta\bar{\rho}_{-1,m'}(\mathbf{r}) r^2
		\frac{\sqrt{6}}{2}\exp(2i\phi) Y_{2,1}^*\\
	&=&\frac{1}{MR^2} \int d^3 x 
		\delta\bar{\rho}_{-1,m'}(\mathbf{r}) r^2
		\frac{\sqrt{6}}{2}\exp(2i\phi) Y_{2,-1}^*\\
	\Rightarrow\kappa_{-1,m'}^+&=&\sqrt{3/2}\kappa_{-1,m'}\\
	\kappa_{1,m'}^+&=&-\frac{1}{MR^2} \int d^3 x 
		\delta\bar{\rho}_{1,m'}(\mathbf{r}) r^2
		2\cot\theta\exp(i\phi) Y_{2,2}^*(\theta,\phi)\\
	&=&-\frac{1}{MR^2} \int d^3 x 
		\delta\bar{\rho}_{1,m'}(\mathbf{r}) r^2
		\frac{1}{2}\sqrt{\frac{15}{2\pi}}\sin\theta\cos\theta\exp(-i\phi)\\
	&=&\frac{1}{MR^2} \int d^3 x 
		\delta\bar{\rho}_{1,m'}(\mathbf{r}) r^2 Y_{2,1}^*(\theta,\phi)\\
	\Rightarrow\kappa_{1,m'}^+&=&\kappa_{1,m'}\\
	\kappa_{-2,m'}^+&=&-\frac{1}{MR^2} \int d^3 x 
		\delta\bar{\rho}_{-2,m'}(\mathbf{r}) r^2\left[
			-\cot\theta\exp(i\phi) Y_{2,-1}^*(\theta,\phi)+
			\sqrt{\frac{3}{2}}\exp(2i\phi) Y_{2,0}^*
		\right]\\
	&=&-\frac{1}{MR^2} \int d^3 x 
		\delta\bar{\rho}_{-2,m'}(\mathbf{r}) r^2\left[
			-\frac{1}{2}\sqrt{\frac{15}{2\pi}}\cos^2\theta\exp(2i\phi)+
			\frac{1}{4}\sqrt{\frac{15}{2\pi}}\exp(2i\phi)(3\cos^2\theta-1)
		\right]\\
	&=&-\frac{1}{MR^2} \int d^3 x 
		\delta\bar{\rho}_{-2,m'}(\mathbf{r}) r^2
		\frac{1}{4}\sqrt{\frac{15}{2\pi}}\exp(2i\phi)\left[
			-\cos^2\theta+3\cos^2\theta-1
		\right]\\
	&=&\frac{1}{MR^2} \int d^3 x 
		\delta\bar{\rho}_{-2,m'}(\mathbf{r}) r^2 Y_{-2,m}^*(\theta, phi)\\
	\Rightarrow \kappa_{-2,m'}^+&=&\kappa_{-2,m'}\\

The torque in the y direction
=============================

Now we dot with :math:`\hat{y}=\sin\theta\sin\phi\hat{r} +
\cos\theta\sin\phi\hat{\theta} - \cos\phi \hat{\phi}`. Again, the
:math:`\hat{r}` term does not contribute:

.. math::

	T_y &=& -\left(\frac{G M'}{a^3 \omega_0}\right)^2
			\sum_{m,m',\mu,\mu'} \mathcal{U}_{m,m'} \mathcal{U}_{\mu,\mu'}
				\exp(i(\mu'-m')\Omega t + i\Delta_{m,m'})
				\int d^3 x 
					\delta\bar{\rho}_{m,m'}(\mathbf{r})
					r^2 (\cos\theta\sin\phi\hat{\theta} -
						 \cos\phi\hat{\phi}) \left[
						\mathbf{r}\times
						\nabla Y_{2,\mu}^*(\theta, \phi)
					\right]\\
	&=& -\left(\frac{G M'}{a^3 \omega_0}\right)^2
			\sum_{m,m',\mu,\mu'} \mathcal{U}_{m,m'} \mathcal{U}_{\mu,\mu'}
				\exp(i(\mu'-m')\Omega t + i\Delta_{m,m'})
				\int d^3 x 
					\delta\bar{\rho}_{m,m'}(\mathbf{r})
					r^2 \left(
						\cot\theta\sin\phi\frac{\partial}{\partial \phi} +
						 \cos\phi\frac{\partial}{\partial \theta}
					\right) Y_{2,\mu}^*(\theta, \phi)

Averaging over an orbit:

.. math::

	T_y &=& -\left(\frac{G M'}{a^3 \omega_0}\right)^2
		\sum_{m,m',\mu} \mathcal{U}_{m,m'} \mathcal{U}_{\mu,m'}
			\exp(i\Delta_{m,m'})
			\int d^3 x 
				\delta\bar{\rho}_{m,m'}(\mathbf{r})
				r^2 \left(
					i\mu\cot\theta\sin\phi Y_{2,\mu}^*(\theta, \phi) +
					 \cos\phi\frac{\partial Y_{2,\mu}^*(\theta, \phi)}
								  {\partial \theta}
				\right)\\
	&=& -\left(\frac{G M'}{a^3 \omega_0}\right)^2
		\sum_{m,m',\mu} \mathcal{U}_{m,m'} \mathcal{U}_{\mu,m'}
			\exp(i\Delta_{m,m'})
			\int d^3 x 
				\delta\bar{\rho}_{m,m'}(\mathbf{r})
				r^2 \left\{
					i\mu\cot\theta\sin\phi Y_{2,\mu}^*(\theta, \phi) +
					\cos\phi\left[\mu\cot\theta Y_{2,\mu}^* +
						\sqrt{(2-\mu)(3+\mu)}\exp(i\phi) Y_{2,\mu+1}^*
					\right]
				\right\}\\
	&=& -\left(\frac{G M'}{a^3 \omega_0}\right)^2
		\sum_{m,m',\mu} \mathcal{U}_{m,m'} \mathcal{U}_{\mu,m'}
			\exp(i\Delta_{m,m'})
			\int d^3 x 
				\delta\bar{\rho}_{m,m'}(\mathbf{r})
				r^2 \left[
					\mu\cot\theta\exp(i\phi) Y_{2,\mu}^*(\theta, \phi) +
					\cos\phi\sqrt{(2-\mu)(3+\mu)}\exp(i\phi) Y_{2,\mu+1}^*
				\right]\\
	&=& -\left(\frac{G M'}{a^3 \omega_0}\right)^2
		\sum_{m,m',\mu} \mathcal{U}_{m,m'} \mathcal{U}_{\mu,m'}
			\exp(i\Delta_{m,m'})
			\int d^3 x 
				\delta\bar{\rho}_{m,m'}(\mathbf{r})
				r^2 \left[
					\mu\cot\theta\exp(i\phi) Y_{2,\mu}^*(\theta, \phi) +
					\frac{\sqrt{(2-\mu)(3+\mu)}}{2}(\exp(2i\phi)+1)
					Y_{2,\mu+1}^*
				\right]\\
	&=& -\frac{T_0}{MR^2}
		\sum_{m,m'} \mathcal{U}_{m,m'}
						\left(\kappa^-_{m,m'}\mathcal{U}_{m-1,m'}+
					     \kappa^+_{m,m'}\mathcal{U}_{m+1,m'}\right)
						\exp(i\Delta_{m,m'})

With the same :math:`\kappa^+_{m,m'}` and :math:`\kappa^-_{m,m'}` as for
:math:`T_x`. So taking the real part gives:

.. math::

	T_y = -T_0
		\sum_{m,m'} \mathcal{U}_{m,m'}
						\left(\kappa^-_{m,m'}\mathcal{U}_{m-1,m'}+
					     \kappa^+_{m,m'}\mathcal{U}_{m+1,m'}\right)
						\cos(\Delta_{m,m'})

So the real part is proportional to :math:`\cos\Delta_{m,m'}` and so for small
tidal dissipation it is independent of the dissipation, as expected since
this term is responsible for the precession.
