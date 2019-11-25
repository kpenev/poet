****************************************
Inclined and Eccentric Orbital Evolution
****************************************

The implementation of the evolution of inclined/eccentric orbits is based on
the formalism of Lai 2012, hence we will use the same notation.

The first thing is to derive the spherical harmonic expansion of the tidal
potential. 

.. math::

	U(\mathbf{r}, t) = \frac{GM'}{|\mathbf{r}_{M'}|}\left(
		\frac{\mathbf{r}\cdot\mathbf{r}_{M'}}{\left|\mathbf{r}_{M'}\right|^2}
        -
        \frac{\left|\mathbf{r_{M'}}\right|}{\left|\mathbf{r} - \mathbf{r}_{M'}\right|}
    \right)

where :math:`\mathbf{r}_{M'}(t)` is the vector from the centers of :math:`M` to
:math:`M'` and :math:`\mathbf{r}` is the point where the potential is being
evaluated. Obviously for a circular orbit
:math:`\left|\mathbf{r}_{M'}(t)\right|` is a constant.

Letting :math:`r(t) \equiv \left|\mathbf{r}_{M'}(t)\right|` We start with the
expansion in a coordinate system with :math:`\mathbf{\hat{z}}=\mathbf{\hat{L}}`,
and :math:`\mathbf{\hat{y}}=\mathbf{\hat{S}}\times\mathbf{\hat{L}}`.

Let us define this expansion as:

.. math::

	U(\mathbf{r}, t) = \frac{GM'}{r(t)}\sum_{l,m} c_{l,m}(\rho) 
		Y_{l,m}(\theta, \phi)

With :math:`\mathbf{r}=(\rho, \theta, \phi)` in spherical coordinates.

Clearly:

.. math::

	c_{l,m}(\rho, t) = \int_{0}^{2\pi} d\phi \int_{0}^{\pi} d\theta
		\tilde{U}(\mathbf{r}, t) Y_{l,m}^*(\theta, \phi)

where :math:`\tilde{U}(\mathbf{r}, t)\equiv U(\mathbf{r}, t)r(t)/(GM')`\ .

The integrals are calulated the easiest by writing both :math:`U(\mathbf{r}, t)`
and :math:`Y_{l,m}(\theta, \phi)` in a coondinate system where
:math:`\hat{\tilde{z}}` points from :math:`M` to :math:`M'`. Tilde will be used
to identify coordinates in that system.

.. math::

	\tilde{U}(\mathbf{r}, t) = \tilde{z} - \frac{1}{\sqrt{(1-\tilde{z})^2 +
		\tilde{x}^2 + \tilde{y}^2}}

where :math:`\tilde{x}`\ , :math:`\tilde{y}` and :math:`\tilde{z}` are the
coordinates of the point the potential is being evaluated scaled by
:math:`r(t)`\ :

.. math::

	\tilde{x} & = & \tilde{\rho}\sin\tilde{\theta}\cos\tilde{\phi}\\
	\tilde{y} & = &\tilde{\rho}\sin\tilde{\theta}\sin\tilde{\phi}\\
	\tilde{z} & = &\tilde{\rho}\cos\tilde{\theta}

with :math:`\tilde{\rho}` also scaled by :math:`r(t)`.

With these we have:

.. math::

	\tilde{U}(\tilde{\rho}, \tilde{\theta}, t)
    =
    \tilde{\rho}\cos\tilde{\theta}
    -
    \frac{1}{\sqrt{1-2\tilde{\rho}\cos\tilde{\theta} + \tilde{\rho}^2}}

As expected there is no :math:`\tilde{\phi}` dependence.

Now we need to write :math:`Y_{l,m}` in the tilde coordinate system: 

.. math::

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

Where coordinates with "rot" subscript are in a coordinate system where
:math:`\mathbf{\hat{z}_{rot}}=\mathbf{\hat{L}}` and :math:`\hat{x}_{rot}` points
from the center of :math:`M` to the center of :math:`M'`, and
:math:`\Delta\phi(t)` is the angle between :math:`\hat{x}_{rot}` and
:math:`\hat{x}` (:math:`\Delta\phi(t)=\Omega t` for a circular orbit with
angular velocity :math:`\Omega`).

Note that:

.. math::

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

We can then use Mathematica to evaluate:

.. math::

	c_{0,0} = c_{1,0} = c_{2,\pm1} & = & 0\\
	\quad c_{1,\pm 1} & = & \pm\sqrt{\frac{2\pi}{3}}\tilde{\rho}
		\exp(-im\Delta\phi(t))\\
	\quad c_{2,0} & = & \sqrt{\frac{\pi}{5}}\rho^2\exp(-im\Delta\phi(t))\\
	\quad c_{2,\pm 2} & = & -\sqrt{\frac{3\pi}{10}}\rho^2\exp(-im\Delta\phi(t))

The :math:`c_{1,\pm 1}` coefficients represent the gravitational acceleration of
the center of mass of :math:`M` due to :math:`m`, and :math:`c_{2,*}` are the
lowest order tidal potential terms, so we have reproduced eq. (4) of Lai 2012.

Now we need to transform the :math:`Y_{2,m}` functions to a coordinate system where
the :math:`z` axis is along :math:`\hat{S}` and the :math:`y` axis is along
:math:`\hat{S}\times\hat{L}`\ . We will use primes for the coordinates in the new
system. We have the following relations:

.. math::

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

Then using mathematica we show that the transformation between the prime and
non-prime coordinate system is indeed given by the Wigner D matrices quoted
in Lai 2012 (Equations 6-11).

The equivalent of Lai 2014 eq. 12 is then:

.. math::

	U(\mathbf{r}, t)=-\sum_{m,m'} U_{m,m'}\rho^2 Y_{2,m}(\theta',\phi')
		\frac{\exp\left(-im'\Delta\phi(t)\right)a^3}{r^3(t)}

with

.. math::
	U_{m,m'} \equiv \frac{GM'}{a^3}W_{2,m'}D_{m,m'}(\Theta)

This is where we diverge from Lai 2012 because we wish to consider elliptical
orbits. For general elliptical orbits it is convenient to define the origin of
time so that the planet is at periastron at :math:`t=0`\ . Further we will take
:math:`\Delta\phi(t)=\phi(t)+\phi_0` where :math:`\phi(0)=0`\ , which implies that
:math:`\phi_0` is the angle between periastron and :math:`\hat{x}` or
:math:`90^\circ` less than the angle between periastron and
:math:`\hat{S}\times\hat{L}`\ . With these definitions:

.. math::

	r(t) & = & a(1-e\cos u)\\
	\tan\left(\frac{\phi}{2}\right)&=&\sqrt{\frac{1+e}{1-e}}
		\tan\left(\frac{u}{2}\right)\\
	\Rightarrow
	1+\tan^2\left(\frac{\phi}{2}\right)&=&\frac{1-e\cos u}
		{(1-e)\cos^2\left(\frac{u}{2}\right)}\\
	1-\tan^2\left(\frac{\phi}{2}\right)&=&\frac{\cos u - e}
		{(1-e)\cos^2\left(\frac{u}{2}\right)}\\
	\Rightarrow
	\cos\phi & = & \frac{\cos u - e}{1 - e\cos u}\\
	\sin\phi & = & \frac{\sqrt{1-e^2}\sin u}{1 - e\cos u}

Where :math:`u` is the eccentric anomaly:

.. math::

	u-e\sin u = \Omega t,\quad\Omega=\sqrt{\frac{G(M+M')}{a^3}}

Differentiating:

.. math::

	dt=\frac{1-e\cos u}{\Omega}du

Since the orbital solution :math:`r(t)` and :math:`\Delta \phi(t)` is periodic
with a period of :math:`2\pi/\Omega`\ , we can expand:

.. math::

	\frac{a^3\exp(-im'\Delta \phi(t))}{r^3(t)}=\sum_s p_{m',s}
		\exp\left(-i s \Omega t\right)

Expressions for the :math:`p_{m',s}(e)` coefficients are derived :doc:`here
<inclination_eccentricity_pms1>` or :doc:`here <inclination_eccentricity_pms2>`\
.

Hence, our tidal potential can be written exactly as in Lai (2012), eq. 12,
except with :math:`m'` not limited to only 0 and 2:

.. math::

	U(\mathbf{r}, t)=-\sum_{m,m'} U_{m,m'}\rho^2 Y_{2,m}(\theta',\phi')
		\exp\left(-im'\Omega t\right)

with

.. math::

	U_{m,m'} \equiv \frac{GM'}{a^3} \mathcal{U}_{m,m'}\equiv
		\frac{GM'}{a^3} \sum_s W_{2,s}D_{m,s}(\Theta) p_{s,m'}

where we have switched the :math:`m'` and :math:`s` coefficients. Here is a
:doc:`table <inclination_eccentricity_Ummtable>` of :math:`W_{2,m'}D_{m,m'}`\ .

From here we proceed following Lai (2012) again, but we have more than 6
independent timelags if the orbit is eccentric (for circular orbits,
:math:`p_{m',s}=\delta_{m',s}`\ ). 

The ansatz:

.. math::

	\mathbf{\xi}_{m,s}(\mathbf{r},t)&=&
		\frac{U_{m,s}}{\omega_0^2}\mathbf{\bar{\xi}}_{m,s}(\mathbf{r})
		\exp(-is\Omega t + i\Delta_{m,s})\\
	\delta\rho(\mathbf{r},t)&=&\frac{U_{m,s}}{\omega_0^2}
		\delta\bar{\rho}_{m,s}(\mathbf{r}) 
		\exp(-is\Omega t + i\Delta_{m,s})\\
	\Delta_{m,s}&=&\tilde{\omega}_{m,s}t_{m,s}

with :math:`\delta\bar{\rho}_{m,s}=-\nabla\cdot(\rho\mathbf{\bar{\xi}}_{m,s})`,
:math:`\tilde{\omega}_{m,s}=s\Omega-m\Omega_s`\ , where :math:`\Omega_s` is the
spin angular velocity of :math:`M`\ , and :math:`\omega_0\equiv\sqrt{GM/R^3}` is
the dynamical frequency of :math:`M`\ .

Here are the detailed devirations of :doc:`the tidal torque
<inclination_eccentricity_torque>` and :doc:`the tidal power
<inclination_eccentricity_power>`\ .

As noted before, for general eccentric orbits, the number of timelags is not
only six, like in Lai (2012), but could be arbitrarily large, depending on
the precision required of the expansion and the value of the eccentricity. In
order to preserve full generality, we allow the user to specify each tidal
lag :math:`\Delta'_{m,m'}\equiv\kappa_{m,m'}\sin(\Delta_{m,m'})` and each love
coefficient :math:`\kappa'_{m,m'}\equiv\kappa_{m,m'}\cos(\Delta_{m,m'})`
separately.

The variables evolved will be the usual orbital
elements (the semimajor axis - :math:`a` and eccentricity :math:`e`\ ), and for
each zone we will use the inclination relative to the orbit - :math:`i`\ .
Finally, one zone will be designated as a reference and for all other zones,
we will follow the evolution of the difference between their argument of
periapsis and that of the referenc zone - :math:`\Delta\omega`\ . Thus, if the
two bodies are split into n zones, the evolution of 1+2n variables will be
followed. The equations for the evolution of these variables are derived 
[here](@ref EccentricEvolutionEquations).

The collected equations are:

.. math::

	\dot{a} & = & a\frac{-\dot{E}}{E}\\
	\dot{e} & = & \frac{2(\dot{E}L+2E\dot{L})L(M+M')}{G(MM')^3}\\
	\dot{\theta} & = & \frac{(T_z+\tilde{T}_z)\sin\theta}{L} 
					 - \frac{(T_x+\tilde{T}_x)\cos\theta}{L}
					 - \frac{T_x+\mathscr{T}_x}{S}\\
	\dot{\omega} & = & \frac{(T_y+\tilde{T}_y)\cos\theta}{L\sin\theta}
				     + \frac{T_y+\mathscr{T}_y}{S\sin\theta}\\
	\mathbf{\hat{\tilde{x}}} & = & \left(\sin\theta\sin\tilde{\theta}
							  + \cos\theta\cos\tilde{\theta}\cos\Delta\omega
						\right)\mathbf{\hat{x}}
						+ \cos\tilde{\theta}\sin\Delta\omega\mathbf{\hat{y}}
						+ \left(\cos\theta\sin\tilde{\theta}
								-
								\sin\theta\cos\tilde{\theta}\cos\Delta\omega
						\right)\mathbf{\hat{z}}\\
	\mathbf{\hat{\tilde{y}}} & = & -\cos\theta\sin\Delta\omega\mathbf{\hat{x}}
					    + \cos\Delta\omega\mathbf{\hat{y}}
					    + \sin\theta\sin\Delta\omega\mathbf{\hat{z}}\\
	\mathbf{\hat{\tilde{z}}} & = & \left(\sin\theta\cos\tilde{\theta}
							  - \cos\theta\sin\tilde{\theta}\cos\Delta\omega
						\right)\mathbf{\hat{x}}
						- \sin\tilde{\theta}\sin\Delta\omega\mathbf{\hat{y}}
						+ \left(\cos\theta\cos\tilde{\theta}
								+
								\sin\theta\sin\tilde{\theta}\cos\Delta\omega
						\right)\mathbf{\hat{z}}

If the system ever gets in a state where the forcing frequency
:math:`\tilde{\omega}\equiv m'\Omega-mS^{conv}/I_{conv}` for some (m,m')
combination from the expansion of the potential (see above), and the
corresponding :math:`\sin(\Delta_{m,m'}(\tilde{\omega})\kappa_{m,m'}` is
discontinuous at zero, it is possible that locks between the spin of some zones
of the bodies and the orbit will be established. In that case, for each locked
zone, let us  split the tidal dissipation torque and power into components
:math:`T_x^\pm`\ , :math:`T_z^\pm` and :math:`\dot{E}^\pm` with the (+) terms
assuming that the spin frequency is just above the lock and the (-) terms
assuming it is just below. Further, let :math:`T_x^0`\ , :math:`T_z^0` and
:math:`\dot{E}^0` be the tidal torques and power due to all other zones. We can
imagine taking an infinitesimally small timestep, during which the not-locked
components will contribute as usual, but over a fraction (:math:`\lambda`\ ) the
timestep the (+) locked components contribute and over the remaining fraction
(:math:`1-\lambda`\ ) the (-) locked components contribute.

In order to maintain the lock, we must have for each locked zone (denoted by
index i):

.. math::

	&&m'\frac{\partial}{\partial t}\sqrt{\frac{G(M+M')}{a^3}}=
		m\frac{\partial}{\partial t}\frac{S_i}{I_i}\\
	\Rightarrow && -\frac{3m'}{2}\sqrt{\frac{G(M+M')}{a^5}}\dot{a}=
		m\frac{\dot{S}_i}{I_i}
		-
		m\frac{S_i\dot{I}_i}{I_i^2}\\
	\Rightarrow && -\frac{3 m S_i}{2 I_i}\frac{\dot{a}}{a} =
		m\frac{\dot{S}_i}{I_i}
		-
		m\frac{S_i\dot{I}_i}{I_i^2}\\
	\Rightarrow && -\frac{3}{2}\frac{\dot{a}}{a}=
		\frac{\dot{S}_i}{S_i} - \frac{\dot{I}_i}{I}\\
	\Rightarrow && -\frac{3}{2}\frac{\dot{E}^0 + 
									 \sum_k\left[\lambda_k\dot{E}_k^+ +
									 (1-\lambda_k)\dot{E}_k^-\right]}{E}
				   =
				   \frac{\dot{I_i}}{I_i}
				   -
				   \frac{\dot{S}_i^0 + \lambda_i\dot{S}_i^+ 
						 + (1-\lambda_i)\dot{S}_i^-}{S_i}\\
	\Rightarrow && \lambda_i\left(\frac{\dot{S}_i^+ 
								  -
							      \dot{S}_i^-}{S_i}\right)
								-\sum_k \lambda_k
								\left(1.5\frac{\dot{E}_k^+ - \dot{E}_k^-}{E}
								\right)
					=\frac{\dot{I_i}}{I_i}
				   -
				   \frac{\dot{S}_i^0 + \dot{S}_i^-}{S_i}
				   +\frac{3}{2}\frac{\dot{E}^0+\sum_k\dot{E}_k^-}{E}

And the lock is maintained as long as :math:`0<\lambda<1`\ .

The resulting evolution equations are then the same as before, but with:

.. math::

	T_x & = & T_x^0 + \lambda T_x^+ + (1-\lambda)T_x^-\\
	T_z & = & T_z^0 + \lambda T_z^+ + (1-\lambda)T_z^-\\
	\dot{E} & = & \dot{E}^0 + \lambda\dot{E}^+ + (1-\lambda)\dot{E}^-
