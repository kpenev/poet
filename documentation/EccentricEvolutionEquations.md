Evolution Equations From Tidal Torque+Power {#EccentricEvolutionEquations}
===========================================

The orbital energy and angular momentum are:
\f{eqnarray*}{
	E&=&-\frac{GMM'}{2a}\\
	L&=&\frac{MM'}{M+M'}a^2\Omega\sqrt{1-e^2}
			=GMM'\sqrt{\frac{(1-e^2)MM'}{2E(M+M')}}
\f}
Hence:
\f{eqnarray*}{
	\dot{a}&=&-\frac{GMM'}{2E^2}\dot{E}\\
	\dot{e}&=&\frac{2(\dot{E}L+2E\dot{L})L(M+M')}{G(MM')^3}\\
\f}

Now consider a single zone subject to tidal torques. We will use the
following variables:
 - \f$\mathbf{S}\f$ and \f$S\f$: the spin angular momentum vector of the zone
   and its absolute value 
 - \f$\mathbf{L}\f$ and \f$L\f$: the orbital angular momentum vector and its
   absolute value 
 - \f$\theta\f$: angle between the angular momentum vector of the zone and
   the angular momentum of the orbit (inclination).
 - \f$\omega\f$: the argument of periapsis of the orbit with a plane of
   reference perpendicular to \f$\mathbf{S}\f$.
 - \f$\mathbf{T}\f$: the tidal torque vector acting on this zone (and of
   course with a negative sign on the orbit).
 - \f$\mathbf{\tilde{T}}\f$: the negative of the tidal torques on the orbit
   due to other zones.
 - \f$\bhat{p}\f$: a unit vector along the periapsis of the orbit
 - \f$\bhat{z}\f$: a unit vector along \f$\mathbf{S}\f$.
 - \f$\bhat{y}\f$: a unit vector along the ascending node of the orbit.
 - \f$\bhat{x}\f$: \f$\bhat{y}\times\bhat{z}\f$.

We will use primed quantities to denote the updated value of a quantity after
an infinitesimal times step, and will use x, y and z indices to indicate
projections of quantities alonge \f$\bhat{x}\f$, \f$\bhat{y}\f$ and
\f$\bhat{z}\f$ respectively.

\f{eqnarray*}{
	\mathbf{S} &=& S\bhat{z}\\
	\mathbf{L} &=& L\sin\theta\bhat{x}+L\cos\theta\bhat{z}\\
	\bhat{y} &=& \frac{\mathbf{S}\times\mathbf{L}}{LS\sin\theta}\\
	\bhat{p} &=& - \sin\omega\cos\theta\bhat{x} + \cos\omega\bhat{y}
				 + \sin\omega\sin\theta\bhat{z}\\
	\mathbf{S}' &=& T_xdt\bhat{x} + T_ydt\bhat{y} + (S+T_zdt)\bhat{z}\\
	\mathbf{L}' &=& (L\sin\theta-T_xdt-\tilde{T}_xdt)\bhat{x}
					- (T_y+\tilde{T}_y)dt\bhat{y}
					+ (L\cos\theta-T_zdt-\tilde{T}_zdt)\bhat{z}\\
	\frac{1}{S'} &=& \frac{1}{S} - \frac{T_z}{S^2}dt\\
	\frac{1}{L'} &=& \frac{1}{L} 
					 + \frac{(T_x+\tilde{T}_x)\sin\theta
							 + (T_z+\tilde{T}_z)\cos\theta}{L^2}dt\\
	\frac{1}{L'S'} &=& \frac{1}{LS}\left(1 - \frac{T_z}{S}dt
					 					+ \frac{(T_x+\tilde{T}_x)\sin\theta
										+ (T_z+\tilde{T}_z)\cos\theta}{L}dt
									\right)\\
	\Rightarrow \frac{d}{dt}\left(\frac{1}{LS}\right) &=&
		\frac{1}{LS}\left(\frac{(T_x+\tilde{T}_x)\sin\theta
						  + (T_z+\tilde{T}_z)\cos\theta}{L} - \frac{T_z}{S}
					\right)
\f}

In order to derive the evolution rate for the inclination:
\f{eqnarray*}{
	\cos\theta &=& \frac{\mathbf{S}\cdot\mathbf{L}}{LS}\\
	\cos\theta' &=& \frac{\mathbf{S}'\cdot\mathbf{L}'}{L'S'}\\
				&=& \frac{L\sin\theta T_xdt + LS\cos\theta
						  + L\cos\theta T_zdt - S(T_z+\tilde{T}_z)dt}{LS}
					\left(1 - \frac{T_z}{S}dt
						  + \frac{(T_x+\tilde{T}_x)\sin\theta
						  + (T_z+\tilde{T}_z)\cos\theta}{L}dt
					\right)\\
				&=& \cos\theta -\frac{\cos\theta T_z}{S}dt
					+ \frac{(T_x+\tilde{T}_x)\sin\theta\cos\theta
							+ (T_z+\tilde{T}_z)\cos^2\theta}{L}dt
					+\frac{\sin\theta T_x}{S}dt
					+ \frac{\cos\theta T_z}{S}dt
					- \frac{T_z+\tilde{T}_z}{L}dt\\
				&=& \cos\theta + \frac{(T_x+\tilde{T}_x)\sin\theta\cos\theta
							- (T_z+\tilde{T}_z)\sin^2\theta}{L}dt
					+\frac{\sin\theta T_x}{S}dt\\
	\Rightarrow \dot{\theta} &=& \frac{(T_z+\tilde{T}_z)\sin\theta}{L} 
								 - \frac{(T_x+\tilde{T}_x)\cos\theta}{L}
								 - \frac{T_x}{S}
\f}

Now to derive the evolution rate for the argument of periapsis:
\f[
	\dot{\omega} = -\frac{1}{\sin\omega}\frac{d(\bhat{p}\cdot\bhat{y})}{dt}
				 = -\frac{1}{\sin\omega}\left(
						\bhat{y}\cdot\frac{d\bhat{p}}{dt}
						+ \bhat{p}\cdot\frac{d\bhat{y}}{dt}
				   \right)
\f]

\f{eqnarray*}{
	\frac{d}{dt}\left(\mathbf{S}\times\mathbf{L}\right)
		&=& - S(T_x+\tilde{T}_x)\bhat{y} + S(T_y+\tilde{T}_y)\bhat{x}
			- L\cos\theta T_x\bhat{y} + L\cos\theta T_y\bhat{x}
			- L\sin\theta T_y\bhat{z} + L\sin\theta T_z\bhat{y}\\

		&=& \left[S(T_y+\tilde{T}_y) + L\cos\theta T_y\right]\bhat{x}
			+ \left[L\sin\theta T_z - L\cos\theta T_x
					-S(T_x+\tilde{T}_x)
			  \right]\bhat{y}
			- L\sin\theta T_y\bhat{z}
\f}

From \f$\dot{\theta}\f$:
\f{eqnarray*}{
	\frac{d}{dt}\left(\frac{1}{\sin\theta}\right) &=&
		\frac{(T_x+\tilde{T}_x)\cos^2\theta}{L\sin^2\theta}
		+ \frac{T_x\cos\theta}{S\sin^2\theta}
		- \frac{(T_z+\tilde{T}_z)\cos\theta}{L\sin\theta} 
\f}

Combining \f$\frac{d}{dt}\left(\frac{1}{LS}\right)\f$,
\f$\frac{d}{dt}\left(\mathbf{S}\times\mathbf{L}\right)\f$ and 
\f$\frac{d}{dt}\left(\frac{1}{\sin\theta}\right)\f$:
\f{eqnarray*}{
	\frac{d}{dt}\bhat{y}
		&=& \left\{\frac{(T_x+\tilde{T}_x)\sin\theta
				  + (T_z+\tilde{T}_z)\cos\theta}{L} - \frac{T_z}{S}
			\right\}\bhat{y}\\
		&& +
			\frac{1}{LS\sin\theta}\left\{
				\left[S(T_y+\tilde{T}_y) + L\cos\theta T_y\right]\bhat{x}
				+ \left[L\sin\theta T_z - L\cos\theta T_x
						- S(T_x+\tilde{T}_x)
				  \right]\bhat{y}
				- L\sin\theta T_y\bhat{z}
			\right\}\\
		&& + 
			\left\{\frac{(T_x+\tilde{T}_x)\cos^2\theta}{L\sin\theta}
				+ \frac{T_x\cos\theta}{S\sin\theta}
				- \frac{(T_z+\tilde{T}_z)\cos\theta}{L}
			\right\}\bhat{y}\\
		&=& \left[\frac{T_y+\tilde{T}_y}{L\sin\theta} + 
				  \frac{\cos\theta T_y}{S\sin\theta}\right]\bhat{x}\\
		&& +
			\left[\frac{(T_x+\tilde{T}_x)\sin\theta}{L}
				  + \frac{(T_z+\tilde{T}_z)\cos\theta}{L} 
				  - \frac{T_z}{S} - \frac{\cos\theta T_x}{S\sin\theta}
				  + \frac{T_z}{S} - \frac{T_x+\tilde{T}_x}{L\sin\theta}
				  + \frac{(T_x+\tilde{T}_x)\cos^2\theta}{L\sin\theta}
				  + \frac{T_x\cos\theta}{S\sin\theta}
				  - \frac{(T_z+\tilde{T}_z)\cos\theta}{L}
			\right]\bhat{y}\\
		&& - \frac{T_y}{S}\bhat{z}\\
		&=& \left[\frac{T_y+\tilde{T}_y}{L\sin\theta} + 
				  \frac{T_y\cos\theta}{S\sin\theta}\right]\bhat{x}
			- \frac{T_y}{S}\bhat{z}\\
	\Rightarrow -\frac{\bhat{p}}{\sin\omega}\cdot\frac{d}{dt}\bhat{y}
		&=& \frac{(T_y+\tilde{T}_y)\cos\theta}{L\sin\theta}
			+ \frac{T_y\cos^2\theta}{S\sin\theta}
			+ \frac{T_y\sin\theta}{S}\\
		&=& \frac{(T_y+\tilde{T}_y)\cos\theta}{L\sin\theta}
			+ \frac{T_y}{S\sin\theta}
\f}

The evolution of the direction of periapsis:
\f{eqnarray*}{
	\frac{d}{dt}\bhat{p} &=& -\mathbf{T}\cdot\bhat{p}\frac{\mathbf{L}}{L^2}\\
		&=& \left(\frac{T_x\sin\omega\cos\theta}{L}
				  - \frac{T_y\cos\omega}{L}
				  - \frac{T_z\sin\omega\sin\theta}{L}
			\right)
			\left(\sin\theta\bhat{x}+\cos\theta\bhat{z}\right)\\
	\Rightarrow -\frac{\bhat{y}}{\sin\omega}\cdot\frac{d}{dt}\bhat{p} &=& 0
\f}

So we get:
\f[
	\dot{\omega} = \frac{(T_y+\tilde{T}_y)\cos\theta}{L\sin\theta}
				   + \frac{T_y}{S\sin\theta}
\f]

Finally:
\f{eqnarray*}{
	\dot{S} &=& T_z\\
	\dot{L} &=& -T_x\sin\theta - T_z\cos\theta
\f}

What remanains is to find \f$\tilde{T}_x\f$, \f$\tilde{T}_y\f$ and
\f$\tilde{T}_z\f$. All that is necessary is to express the coornidate system
unit vectors of all other zones in terms of the ones for this zone. We will
use \f$\bhat{\tilde{x}}\f$, \f$\bhat{\tilde{y}}\f$ and \f$\bhat{\tilde{z}}\f$
to refer to the unit vectors of another zone, and we will denote the
difference between this zone's argument of periapsis and the second zone by
\f$\Delta\omega\f$.
Clearly:
\f[
	\bhat{\tilde{y}}=-\cos\theta\sin\Delta\omega\bhat{x}
					+ \cos\Delta\omega\bhat{y}
					+ \sin\theta\sin\Delta\omega\bhat{z}
\f]

Next:
\f{eqnarray*}{
	\bhat{\tilde{z}}&=&\cos\tilde{\theta}\bhat{L}+
					  \sin\tilde{\theta}\bhat{L}\times\bhat{\tilde{y}}\\
				    &=&\sin\theta\cos\tilde{\theta}\bhat{x}
						+ \cos\theta\cos\tilde{\theta}\bhat{z}
						+ \sin\tilde{\theta}\left(
							\sin\theta\cos\Delta\omega\bhat{z}
							- \sin^2\theta\sin\Delta\omega\bhat{y}
							- \cos^2\theta\sin\Delta\omega\bhat{y}
							- \cos\theta\cos\Delta\omega\bhat{x}
						\right)\\
					&=& \left(\sin\theta\cos\tilde{\theta}
							  - \cos\theta\sin\tilde{\theta}\cos\Delta\omega
						\right)\bhat{x}
						- \sin\tilde{\theta}\sin\Delta\omega\bhat{y}
						+ \left(\cos\theta\cos\tilde{\theta}
								+
								\sin\theta\sin\tilde{\theta}\cos\Delta\omega
						\right)\bhat{z}
\f}

Finally:
\f{eqnarray*}{
	\bhat{\tilde{x}}&=& \sin\tilde{\theta}\bhat{L} 
						 - \cos\tilde{\theta}\bhat{L}\times\bhat{\tilde{y}}\\
					&=& \left(\sin\theta\sin\tilde{\theta}
							  + \cos\theta\cos\tilde{\theta}\cos\Delta\omega
						\right)\bhat{x}
						+ \cos\tilde{\theta}\sin\Delta\omega\bhat{y}
						+ \left(\cos\theta\sin\tilde{\theta}
								-
								\sin\theta\cos\tilde{\theta}\cos\Delta\omega
						\right)\bhat{z}
\f}

Crosscheck that
\f$\bhat{\tilde{x}}\times\bhat{\tilde{y}}=\bhat{\tilde{z}}\f$:

\f{eqnarray*}{
	\bhat{\tilde{x}}\times\bhat{\tilde{y}}
		&=& \left(\sin\theta\sin\tilde{\theta}
							  + \cos\theta\cos\tilde{\theta}\cos\Delta\omega
			\right)\cos\Delta\omega\bhat{z}
			-
			\left(\sin\theta\sin\tilde{\theta}
				  + \cos\theta\cos\tilde{\theta}\cos\Delta\omega
			\right)\sin\theta\sin\Delta\omega\bhat{y}\\
			&& +
			\cos\tilde{\theta}\sin\Delta\omega \cos\theta\sin\Delta\omega
			\bhat{z}
			+
			\cos\tilde{\theta}\sin\Delta\omega \sin\theta\sin\Delta\omega
			\bhat{x}\\
			&& -
			\left(\cos\theta\sin\tilde{\theta}
				  - \sin\theta\cos\tilde{\theta}\cos\Delta\omega
			\right)\cos\theta\sin\Delta\omega\bhat{y}
			-
			\left(\cos\theta\sin\tilde{\theta}
				  - \sin\theta\cos\tilde{\theta}\cos\Delta\omega
			\right)\cos\Delta\omega\bhat{x}\\
		&=& \left(
			\sin\theta\cos\tilde{\theta}\sin^2\Delta\omega
			- \cos\theta\sin\tilde{\theta}\cos\Delta\omega
			+ \sin\theta\cos\tilde{\theta}\cos^2\Delta\omega\right)\bhat{x}
			-\left(\sin^2\theta\sin\tilde{\theta}\sin\Delta\omega
				  + \sin\theta\cos\theta\cos\tilde{\theta}
					\sin\Delta\omega\cos\Delta\omega
				  + \cos^2\theta\sin\tilde{\theta}\sin\Delta\omega
				  - \sin\theta\cos\theta\cos\tilde{\theta}
					\sin\Delta\omega\cos\Delta\omega
			\right)\bhat{y}
			+\left(\sin\theta\sin\tilde{\theta}\cos\Delta\omega
				   + \cos\theta\cos\tilde{\theta}\cos^2\Delta\omega
				   + \cos\theta\cos\tilde{\theta}\sin^2\Delta\omega
			\right)\bhat{z}\\
		&=& \left(
			\sin\theta\cos\tilde{\theta}
			- \cos\theta\sin\tilde{\theta}\cos\Delta\omega\right)\bhat{x}
			-\sin\tilde{\theta}\sin\Delta\omega\bhat{y}
			+\left(\cos\theta\cos\tilde{\theta}
				   + \sin\theta\sin\tilde{\theta}\cos\Delta\omega
			\right)\bhat{z}
\f}
Which is exactly \f$\bhat{\tilde{z}}\f$.

Crosscheck that the direction of the orbital angular momentum matches, i.e.
that
\f$
\sin\tilde{\theta}\bhat{\tilde{x}}+\cos\tilde{\theta}\bhat{\tilde{z}}=
\sin\theta\bhat{x}+\cos\theta\bhat{z}
\f$:

\f{eqnarray*}{
	\sin\tilde{\theta}\bhat{\tilde{x}}+\cos\tilde{\theta}\bhat{\tilde{z}}
		&=& \left(\sin\theta\sin^2\tilde{\theta}
				  + \cos\theta\sin\tilde{\theta}\cos\tilde{\theta}
					\cos\Delta\omega
			\right)\bhat{x}
			+ \sin\tilde{\theta}\cos\tilde{\theta}\sin\Delta\omega\bhat{y}
			+ \left(\cos\theta\sin^2\tilde{\theta}
					- \sin\theta\sin\tilde{\theta}\cos\tilde{\theta}
					  \cos\Delta\omega
			\right)\bhat{z}\\
			&& +
			\left(\sin\theta\cos^2\tilde{\theta}
				  - \cos\theta\sin\tilde{\theta}\cos\tilde{\theta}
					\cos\Delta\omega
			\right)\bhat{x}
			- \sin\tilde{\theta}\cos\tilde{\theta}\sin\Delta\omega\bhat{y}
			+ \left(\cos\theta\cos^2\tilde{\theta}
					+ \sin\theta\sin\tilde{\theta}\cos\tilde{\theta}
					  \cos\Delta\omega
			\right)\bhat{z}\\

		&=& \sin\theta\bhat{x} + \cos\theta\bhat{z}
\f}

Finally, crosscheck that the direction of periapsis is consistent, i.e. that
\f[ 
	- \sin(\omega-\Delta\omega)\cos\tilde{\theta}\bhat{\tilde{x}}
	+ \cos(\omega-\Delta\omega)\bhat{\tilde{y}}
	+ \sin(\omega-\Delta\omega)\sin\tilde{\theta}\bhat{\tilde{z}}=
	- \sin\omega\cos\theta\bhat{x} + \cos\omega\bhat{y}
	+ \sin\omega\sin\theta\bhat{z}
\f]

We will go component by component:
\f{eqnarray*}{
	\bhat{\tilde{p}}\cdot\bhat{x}
		&=& -\sin(\omega-\Delta\omega)\cos\tilde{\theta}
			\left(\sin\theta\sin\tilde{\theta}
				  + \cos\theta\cos\tilde{\theta}\cos\Delta\omega
			\right)\\
			&& -
			\cos(\omega-\Delta\omega)\cos\theta\sin\Delta\omega\\
			&& +
			\sin(\omega-\Delta\omega)\sin\tilde{\theta}
			\left(\sin\theta\cos\tilde{\theta}
				  - \cos\theta\sin\tilde{\theta}\cos\Delta\omega
			\right)\\
		&=& -\cos\theta\cos^2\tilde{\theta}
			\sin(\omega-\Delta\omega)\cos\Delta\omega
			-
			\cos\theta\cos(\omega-\Delta\omega)\sin\Delta\omega
			-
			\cos\theta\sin^2\tilde{\theta}
			\sin(\omega-\Delta\omega)\cos\Delta\omega\\
		&=& - \cos\theta\sin(\omega-\Delta\omega)\cos\Delta\omega
			- \cos\theta\cos(\omega-\Delta\omega)\sin\Delta\omega\\
		&=& - \sin\omega\cos\theta\\
		&=& \bhat{p}\cdot\bhat{x}
\f}

Next:
\f{eqnarray*}{
	\bhat{\tilde{p}}\cdot\bhat{y}
		&=& - \sin(\omega-\Delta\omega)\cos\tilde{\theta}
			  \cos\tilde{\theta}\sin\Delta\omega
			+ \cos(\omega-\Delta\omega)\cos\Delta\omega
			- \sin(\omega-\Delta\omega)\sin\tilde{\theta}
			  \sin\tilde{\theta}\sin\Delta\omega\\
		&=& - \sin(\omega-\Delta\omega)\sin\Delta\omega
			+ \cos(\omega-\Delta\omega)\cos\Delta\omega\\
		&=& \cos\omega\\
		&=& \bhat{p}\cdot\bhat{y}
\f}

Finally:
\f{eqnarray*}{
	\bhat{\tilde{p}}\cdot\bhat{z}
		&=& - \sin(\omega-\Delta\omega)\cos\tilde{\theta}
			  \left(\cos\theta\sin\tilde{\theta}
					- \sin\theta\cos\tilde{\theta}\cos\Delta\omega
			  \right)\\
			&&
			+ \cos(\omega-\Delta\omega)\sin\theta\sin\Delta\omega\\
			&&
			+ \sin(\omega-\Delta\omega)\sin\tilde{\theta}
			  \left(\cos\theta\cos\tilde{\theta}
					+ \sin\theta\sin\tilde{\theta}\cos\Delta\omega
			  \right)\\
		&=& + \sin\theta\cos^2\tilde{\theta}
			  \sin(\omega-\Delta\omega)\cos\Delta\omega
			+ \cos(\omega-\Delta\omega)\sin\theta\sin\Delta\omega
			+ \sin\theta\sin^2\tilde{\theta}
			  \sin(\omega-\Delta\omega)\cos\Delta\omega\\
		&=& \sin\theta\sin(\omega-\Delta\omega)\cos\Delta\omega
			+ \sin\theta\cos(\omega-\Delta\omega)\sin\Delta\omega\\
		&=& \sin\theta\sin\omega\\
		&=& \bhat{p}\cdot\bhat{z}
\f}
