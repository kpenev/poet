*******************************************
Evolution Equations From Tidal Torque+Power
*******************************************

The orbital energy and angular momentum are:

.. math::

	E&=&-\frac{GMM'}{2a}\\
	L&=&\frac{MM'}{M+M'}a^2\Omega\sqrt{1-e^2}
			=GMM'\sqrt{\frac{(1-e^2)MM'}{(-2E)(M+M')}}

Hence:

.. math::

	\dot{a}&=&a\frac{-\dot{E}}{E}\\
	\dot{e}&=&\frac{(M+M')}{G^2(MM')^3}\frac{(\dot{E}L+2E\dot{L})L}{e}\\

Now consider a single zone subject to tidal torques. We will use the
following variables:

  - :math:`\mathbf{S}` and :math:`S`: the spin angular momentum vector of the
    zone and its absolute value 

  - :math:`\mathbf{L}` and :math:`L`: the orbital angular momentum vector and
    its absolute value 

  - :math:`\theta`: angle between the angular momentum vector of the zone and
    the angular momentum of the orbit (inclination).

  - :math:`\omega`: the argument of periapsis of the orbit with a plane of
    reference perpendicular to :math:`\mathbf{S}`.

  - :math:`\mathbf{T}`: the tidal torque vector acting on this zone (and of
    course with a negative sign on the orbit).

  - :math:`\mathbf{\tilde{T}}`: the negative of the tidal torques on the orbit
    due to other zones.

  - :math:`\mathbf{\mathscr{T}}`: the torque on this zone due to coupling to
    other zones (e.g. due to differential rotation coupling or mass exchange).

  - :math:`\mathbf{\hat{p}}`: a unit vector along the periapsis of the orbit

  - :math:`\mathbf{\hat{z}}`: a unit vector along :math:`\mathbf{S}`.

  - :math:`\mathbf{\hat{y}}`: a unit vector along the ascending node of the
    orbit.

  - :math:`\mathbf{\hat{x}}`: :math:`\mathbf{\hat{y}}\times\mathbf{\hat{z}}`.

We will use primed quantities to denote the updated value of a quantity after an
infinitesimal times step, and will use x, y and z indices to indicate
projections of quantities alonge :math:`\mathbf{\hat{x}}`,
:math:`\mathbf{\hat{y}}` and :math:`\mathbf{\hat{z}}` respectively.

.. math::

	\mathbf{S} &=& S\mathbf{\hat{z}}\\
	\mathbf{L} &=& L\sin\theta\mathbf{\hat{x}}+L\cos\theta\mathbf{\hat{z}}\\
	\mathbf{\hat{y}} &=& \frac{\mathbf{S}\times\mathbf{L}}{LS\sin\theta}\\
	\mathbf{\hat{p}} &=& - \sin\omega\cos\theta\mathbf{\hat{x}}
                         + \cos\omega\mathbf{\hat{y}}
                         + \sin\omega\sin\theta\mathbf{\hat{z}}\\
	\mathbf{S}' &=& (T_x+\mathscr{T}_x)dt\mathbf{\hat{x}}
					+ (T_y+\mathscr{T}_y)dt\mathbf{\hat{y}}
					+ (S+T_zdt+\mathscr{T}_zdt)\mathbf{\hat{z}}\\
	\mathbf{L}' &=& (L\sin\theta-T_xdt-\tilde{T}_xdt)\mathbf{\hat{x}}
					- (T_y+\tilde{T}_y)dt\mathbf{\hat{y}}
					+ (L\cos\theta-T_zdt-\tilde{T}_zdt)\mathbf{\hat{z}}\\
	\frac{1}{S'} &=& \frac{1}{S} - \frac{T_z+\mathscr{T}_z}{S^2}dt\\
	\frac{1}{L'} &=& \frac{1}{L} 
					 + \frac{(T_x+\tilde{T}_x)\sin\theta
							 + (T_z+\tilde{T}_z)\cos\theta}{L^2}dt\\
	\frac{1}{L'S'} &=& \frac{1}{LS}\left(1 - \frac{T_z+\mathscr{T}_z}{S}dt
					 					+ \frac{(T_x+\tilde{T}_x)\sin\theta
										+ (T_z+\tilde{T}_z)\cos\theta}{L}dt
									\right)\\
	\Rightarrow \frac{d}{dt}\left(\frac{1}{LS}\right) &=&
		\frac{1}{LS}\left(\frac{(T_x+\tilde{T}_x)\sin\theta
						  + (T_z+\tilde{T}_z)\cos\theta}{L}
						  - \frac{T_z+\mathscr{T}_z}{S}
					\right)

In order to derive the evolution rate for the inclination:

.. math::

	\cos\theta &=& \frac{\mathbf{S}\cdot\mathbf{L}}{LS}\\
	\cos\theta' &=& \frac{\mathbf{S}'\cdot\mathbf{L}'}{L'S'}\\
				&=& \frac{L\sin\theta(T_x+\mathscr{T}_x)dt + LS\cos\theta
						  + L\cos\theta(T_z+\mathscr{T}_z)dt
						  - S(T_z+\tilde{T}_z)dt}{LS}
					\left(1 - \frac{T_z+\mathscr{T}_z}{S}dt
						  + \frac{(T_x+\tilde{T}_x)\sin\theta
						  + (T_z+\tilde{T}_z)\cos\theta}{L}dt
					\right)\\
				&=& \cos\theta -\frac{\cos\theta(T_z+\mathscr{T}_z)}{S}dt
					+ \frac{(T_x+\tilde{T}_x)\sin\theta\cos\theta
							+ (T_z+\tilde{T}_z)\cos^2\theta}{L}dt
					+\frac{\sin\theta(T_x+\mathscr{T}_x)}{S}dt
					+ \frac{\cos\theta(T_z+\mathscr{T}_z)}{S}dt
					- \frac{T_z+\tilde{T}_z}{L}dt\\
				&=& \cos\theta + \frac{(T_x+\tilde{T}_x)\sin\theta\cos\theta
							- (T_z+\tilde{T}_z)\sin^2\theta}{L}dt
					+\frac{\sin\theta(T_x+\mathscr{T}_x)}{S}dt\\
	\Rightarrow \dot{\theta} &=& \frac{(T_z+\tilde{T}_z)\sin\theta}{L} 
								 - \frac{(T_x+\tilde{T}_x)\cos\theta}{L}
								 - \frac{T_x+\mathscr{T}_x}{S}

Now to derive the evolution rate for the argument of periapsis:

.. math::

	\dot{\omega} = -\frac{1}{\sin\omega}\frac{d(\mathbf{\hat{p}}\cdot\mathbf{\hat{y}})}{dt}
				 = -\frac{1}{\sin\omega}\left(
						\mathbf{\hat{y}}\cdot\frac{d\mathbf{\hat{p}}}{dt}
						+ \mathbf{\hat{p}}\cdot\frac{d\mathbf{\hat{y}}}{dt}
				   \right)

.. math::

	\frac{d}{dt}\left(\mathbf{S}\times\mathbf{L}\right)
		&=& - S(T_x+\tilde{T}_x)\mathbf{\hat{y}}
            + S(T_y+\tilde{T}_y)\mathbf{\hat{x}}
			- L\cos\theta(T_x+\mathscr{T}_x)\mathbf{\hat{y}}
			+ L\cos\theta(T_y+\mathscr{T}_y)\mathbf{\hat{x}}
			- L\sin\theta(T_y+\mathscr{T}_y)\mathbf{\hat{z}}
			+ L\sin\theta(T_z+\mathscr{T}_z)\mathbf{\hat{y}}\\
		&=& \left[S(T_y+\tilde{T}_y)
				  + L\cos\theta(T_y+\mathscr{T}_y)\right]\mathbf{\hat{x}}
			+ \left[L\sin\theta(T_z+\mathscr{T}_z)
					- L\cos\theta(T_x+\mathscr{T}_x) - S(T_x+\tilde{T}_x)
			  \right]\mathbf{\hat{y}}
			- L\sin\theta(T_y+\mathscr{T}_y)\mathbf{\hat{z}}

From :math:`\dot{\theta}`:

.. math::

    \frac{d}{dt}\left(\frac{1}{\sin\theta}\right)
    =
    \frac{(T_x+\tilde{T}_x)\cos^2\theta}{L\sin^2\theta} +
    \frac{(T_x+\mathscr{T}_x)\cos\theta}{S\sin^2\theta} -
    \frac{(T_z+\tilde{T}_z)\cos\theta}{L\sin\theta}

Combining :math:`\frac{d}{dt}\left(\frac{1}{LS}\right)`,
:math:`\frac{d}{dt}\left(\mathbf{S}\times\mathbf{L}\right)` and
:math:`\frac{d}{dt}\left(\frac{1}{\sin\theta}\right)`:

.. math::

	\frac{d}{dt}\mathbf{\hat{y}}
		&=& \left\{\frac{(T_x+\tilde{T}_x)\sin\theta
					+ (T_z+\tilde{T}_z)\cos\theta}{L}
					- \frac{T_z+\mathscr{T}_z}{S}
			\right\}\mathbf{\hat{y}}\\
		&& +
			\frac{1}{LS\sin\theta}\left\{
				\left[S(T_y+\tilde{T}_y)
					  + L\cos\theta(T_y+\mathscr{T}_y)\right]\mathbf{\hat{x}}
				+ \left[L\sin\theta(T_z+\mathscr{T}_z)
						- L\cos\theta(T_x+\mathscr{T}_x)
						- S(T_x+\tilde{T}_x)
				  \right]\mathbf{\hat{y}}
				- L\sin\theta(T_y+\mathscr{T}_y)\mathbf{\hat{z}}
			\right\}\\
		&& + 
			\left\{\frac{(T_x+\tilde{T}_x)\cos^2\theta}{L\sin\theta}
				+ \frac{(T_x+\mathscr{T}_x)\cos\theta}{S\sin\theta}
				- \frac{(T_z+\tilde{T}_z)\cos\theta}{L}
			\right\}\mathbf{\hat{y}}\\
		&=& \left[\frac{T_y+\tilde{T}_y}{L\sin\theta} + 
				  \frac{\cos\theta(T_y+\mathscr{T}_y)}{S\sin\theta}
			\right]\mathbf{\hat{x}}\\
		&& +
			\left[\frac{(T_x+\tilde{T}_x)\sin\theta}{L}
				  + \frac{(T_z+\tilde{T}_z)\cos\theta}{L} 
				  - \frac{T_z+\mathscr{T}_z}{S}
				  - \frac{\cos\theta(T_x+\mathscr{T}_x)}{S\sin\theta}
				  + \frac{T_z+\mathscr{T}_z}{S}
				  - \frac{T_x+\tilde{T}_x}{L\sin\theta}
				  + \frac{(T_x+\tilde{T}_x)\cos^2\theta}{L\sin\theta}
				  + \frac{(T_x+\mathscr{T}_x)\cos\theta}{S\sin\theta}
				  - \frac{(T_z+\tilde{T}_z)\cos\theta}{L}
			\right]\mathbf{\hat{y}}\\
		&& - \frac{T_y+\mathscr{T}_y}{S}\mathbf{\hat{z}}\\
		&=& \left[
                \frac{T_y+\tilde{T}_y}{L\sin\theta}
                + 
                \frac{(T_y+\mathscr{T}_y)\cos\theta}{S\sin\theta}
			\right]\mathbf{\hat{x}}
            -
            \frac{T_y+\mathscr{T}_y}{S}\mathbf{\hat{z}}\\
	\Rightarrow 
        -\frac{\mathbf{\hat{p}}}{\sin\omega}\cdot\frac{d}{dt}\mathbf{\hat{y}}
		&=& \frac{(T_y+\tilde{T}_y)\cos\theta}{L\sin\theta}
			+ \frac{(T_y+\mathscr{T}_y)\cos^2\theta}{S\sin\theta}
			+ \frac{(T_y+\mathscr{T}_y)\sin\theta}{S}\\
		&=& \frac{(T_y+\tilde{T}_y)\cos\theta}{L\sin\theta}
			+ \frac{T_y+\mathscr{T}_y}{S\sin\theta}

The evolution of the direction of periapsis:

.. math::

	\frac{d}{dt}\mathbf{\hat{p}}
		&=& -(\mathbf{T}+\mathbf{\tilde{T}})
			\cdot\mathbf{\hat{p}}\frac{\mathbf{L}}{L^2}\\
		&=& \left(\frac{(T_x+\tilde{T}_x)\sin\omega\cos\theta}{L}
				  - \frac{(T_y+\tilde{T}_y)\cos\omega}{L}
				  - \frac{(T_z+\tilde{T}_z)\sin\omega\sin\theta}{L}
			\right)
			\left(\sin\theta\mathbf{\hat{x}}+\cos\theta\mathbf{\hat{z}}\right)\\
	\Rightarrow
    -
    \frac{\mathbf{\hat{y}}}{\sin\omega}\cdot\frac{d}{dt}\mathbf{\hat{p}}
    &=&
    0

So we get:

.. math::

	\dot{\omega} = \frac{(T_y+\tilde{T}_y)\cos\theta}{L\sin\theta}
				   + \frac{T_y+\mathscr{T}_y}{S\sin\theta}

Finally:

.. math::

	\dot{S} &=& T_z+\mathscr{T}_z\\
	\dot{L} &=& -T_x\sin\theta - T_z\cos\theta

What remanains is to find :math:`\tilde{T}_x`, :math:`\tilde{T}_y` and
:math:`\tilde{T}_z`. All that is necessary is to express the coornidate system
unit vectors of all other zones in terms of the ones for this zone. We will use
:math:`\mathbf{\hat{\tilde{x}}}`, :math:`\mathbf{\hat{\tilde{y}}}` and
:math:`\mathbf{\hat{\tilde{z}}}` to refer to the unit vectors of another zone,
and we will denote the difference between this zone's argument of periapsis and
the second zone by :math:`\Delta\omega`.

Clearly:

.. math::

	\mathbf{\hat{\tilde{y}}}=-\cos\theta\sin\Delta\omega\mathbf{\hat{x}}
					+ \cos\Delta\omega\mathbf{\hat{y}}
					+ \sin\theta\sin\Delta\omega\mathbf{\hat{z}}

Next:

.. math::

	\mathbf{\hat{\tilde{z}}}
    &=&
    \cos\tilde{\theta}\mathbf{\hat{L}}
    +
    \sin\tilde{\theta}\mathbf{\hat{L}}\times\mathbf{\hat{\tilde{y}}}\\
    &=&
    \sin\theta\cos\tilde{\theta}\mathbf{\hat{x}}
    +
    \cos\theta\cos\tilde{\theta}\mathbf{\hat{z}}
    +
    \sin\tilde{\theta}\left(
        \sin\theta\cos\Delta\omega\mathbf{\hat{z}}
        - \sin^2\theta\sin\Delta\omega\mathbf{\hat{y}}
        - \cos^2\theta\sin\Delta\omega\mathbf{\hat{y}}
        - \cos\theta\cos\Delta\omega\mathbf{\hat{x}}
    \right)\\
    &=&
    \left(
        \sin\theta\cos\tilde{\theta}
        -
        \cos\theta\sin\tilde{\theta}\cos\Delta\omega
    \right)\mathbf{\hat{x}}
    -
    \sin\tilde{\theta}\sin\Delta\omega\mathbf{\hat{y}}
    +
    \left(
        \cos\theta\cos\tilde{\theta}
        +
        \sin\theta\sin\tilde{\theta}\cos\Delta\omega
    \right)\mathbf{\hat{z}}

Finally:

.. math::

        \mathbf{\hat{\tilde{x}}}
    &=&
        \sin\tilde{\theta}\mathbf{\hat{L}}
        -
        \cos\tilde{\theta}\mathbf{\hat{L}}\times\mathbf{\hat{\tilde{y}}}\\
    &=&
        \left(
            \sin\theta\sin\tilde{\theta}
            +
            \cos\theta\cos\tilde{\theta}\cos\Delta\omega
        \right)\mathbf{\hat{x}}
        +
        \cos\tilde{\theta}\sin\Delta\omega\mathbf{\hat{y}}
        +
        \left(
            \cos\theta\sin\tilde{\theta}
            -
            \sin\theta\cos\tilde{\theta}\cos\Delta\omega
        \right)\mathbf{\hat{z}}

Crosscheck that :math:`\mathbf{\hat{\tilde{x}}}\times\mathbf{\hat{\tilde{y}}} =
\mathbf{\hat{\tilde{z}}}`:

.. math::

	\mathbf{\hat{\tilde{x}}}\times\mathbf{\hat{\tilde{y}}}
		&=& \left(\sin\theta\sin\tilde{\theta}
							  + \cos\theta\cos\tilde{\theta}\cos\Delta\omega
			\right)\cos\Delta\omega\mathbf{\hat{z}}
			-
			\left(\sin\theta\sin\tilde{\theta}
				  + \cos\theta\cos\tilde{\theta}\cos\Delta\omega
			\right)\sin\theta\sin\Delta\omega\mathbf{\hat{y}}\\
			&& +
			\cos\tilde{\theta}\sin\Delta\omega \cos\theta\sin\Delta\omega
			\mathbf{\hat{z}}
			+
			\cos\tilde{\theta}\sin\Delta\omega \sin\theta\sin\Delta\omega
			\mathbf{\hat{x}}\\
			&& -
			\left(\cos\theta\sin\tilde{\theta}
				  - \sin\theta\cos\tilde{\theta}\cos\Delta\omega
			\right)\cos\theta\sin\Delta\omega\mathbf{\hat{y}}
			-
			\left(\cos\theta\sin\tilde{\theta}
				  - \sin\theta\cos\tilde{\theta}\cos\Delta\omega
			\right)\cos\Delta\omega\mathbf{\hat{x}}\\
		&=& \left(
                \sin\theta\cos\tilde{\theta}\sin^2\Delta\omega
                -
                \cos\theta\sin\tilde{\theta}\cos\Delta\omega
                +
                \sin\theta\cos\tilde{\theta}\cos^2\Delta\omega
            \right)\mathbf{\hat{x}}
			-\left(
                \sin^2\theta\sin\tilde{\theta}\sin\Delta\omega
                +
                \sin\theta\cos\theta\cos\tilde{\theta}
                \sin\Delta\omega\cos\Delta\omega
                +
                \cos^2\theta\sin\tilde{\theta}\sin\Delta\omega
                -
                \sin\theta\cos\theta\cos\tilde{\theta}
                \sin\Delta\omega\cos\Delta\omega
			\right)\mathbf{\hat{y}}
			+\left(
                \sin\theta\sin\tilde{\theta}\cos\Delta\omega
                +
                \cos\theta\cos\tilde{\theta}\cos^2\Delta\omega
                +
                \cos\theta\cos\tilde{\theta}\sin^2\Delta\omega
			\right)\mathbf{\hat{z}}\\
		&=& \left(
                \sin\theta\cos\tilde{\theta}
                -
                \cos\theta\sin\tilde{\theta}\cos\Delta\omega
            \right)\mathbf{\hat{x}}
            -
            \sin\tilde{\theta}\sin\Delta\omega\mathbf{\hat{y}}
            +
            \left(
                \cos\theta\cos\tilde{\theta}
                +
                \sin\theta\sin\tilde{\theta}\cos\Delta\omega
			\right)\mathbf{\hat{z}}

Which is exactly :math:`\mathbf{\hat{\tilde{z}}}`.

Crosscheck that the direction of the orbital angular momentum matches, i.e.
that

:math:`\sin\tilde{\theta}\mathbf{\hat{\tilde{x}}}
+
\cos\tilde{\theta}\mathbf{\hat{\tilde{z}}}
=
\sin\theta\mathbf{\hat{x}}
+
\cos\theta\mathbf{\hat{z}}`:

.. math::

	\sin\tilde{\theta}\mathbf{\hat{\tilde{x}}}
    +
    \cos\tilde{\theta}\mathbf{\hat{\tilde{z}}}
    &=&
    \left(
        \sin\theta\sin^2\tilde{\theta}
        +
        \cos\theta\sin\tilde{\theta}\cos\tilde{\theta}\cos\Delta\omega
    \right)\mathbf{\hat{x}}
    +
    \sin\tilde{\theta}\cos\tilde{\theta}\sin\Delta\omega\mathbf{\hat{y}}
    +
    \left(
        \cos\theta\sin^2\tilde{\theta}
        -
        \sin\theta\sin\tilde{\theta}\cos\tilde{\theta}
        \cos\Delta\omega
    \right)\mathbf{\hat{z}}\\
    &&
    +
    \left(
        \sin\theta\cos^2\tilde{\theta}
        -
        \cos\theta\sin\tilde{\theta}\cos\tilde{\theta}\cos\Delta\omega
    \right)\mathbf{\hat{x}}
    -
    \sin\tilde{\theta}\cos\tilde{\theta}\sin\Delta\omega\mathbf{\hat{y}}
    +
    \left(
        \cos\theta\cos^2\tilde{\theta}
        +
        \sin\theta\sin\tilde{\theta}\cos\tilde{\theta}
        \cos\Delta\omega
    \right)\mathbf{\hat{z}}\\
    &=&
    \sin\theta\mathbf{\hat{x}} + \cos\theta\mathbf{\hat{z}}

Finally, crosscheck that the direction of periapsis is consistent, i.e. that

.. math::

	- \sin(\omega-\Delta\omega)\cos\tilde{\theta}\mathbf{\hat{\tilde{x}}}
	+ \cos(\omega-\Delta\omega)\mathbf{\hat{\tilde{y}}}
	+ \sin(\omega-\Delta\omega)\sin\tilde{\theta}\mathbf{\hat{\tilde{z}}}
    =
	- \sin\omega\cos\theta\mathbf{\hat{x}}
    + \cos\omega\mathbf{\hat{y}}
	+ \sin\omega\sin\theta\mathbf{\hat{z}}

We will go component by component:

.. math::

	\mathbf{\hat{\tilde{p}}}\cdot\mathbf{\hat{x}}
    &=&
    -\sin(\omega-\Delta\omega)\cos\tilde{\theta}
    \left(
        \sin\theta\sin\tilde{\theta}
        +
        \cos\theta\cos\tilde{\theta}\cos\Delta\omega
    \right)\\
    &&
    -\cos(\omega-\Delta\omega)\cos\theta\sin\Delta\omega\\
    &&
    +\sin(\omega-\Delta\omega)\sin\tilde{\theta}
    \left(
        \sin\theta\cos\tilde{\theta}
        -
        \cos\theta\sin\tilde{\theta}\cos\Delta\omega
    \right)\\
    &=&
    -
    \cos\theta\cos^2\tilde{\theta}\sin(\omega-\Delta\omega)\cos\Delta\omega
    -
    \cos\theta\cos(\omega-\Delta\omega)\sin\Delta\omega
    -
    \cos\theta\sin^2\tilde{\theta} \sin(\omega-\Delta\omega)\cos\Delta\omega\\
    &=&
    -
    \cos\theta\sin(\omega-\Delta\omega)\cos\Delta\omega
    -
    \cos\theta\cos(\omega-\Delta\omega)\sin\Delta\omega\\
    &=&
    - \sin\omega\cos\theta\\
    &=&
    \mathbf{\hat{p}}\cdot\mathbf{\hat{x}}

Next:

.. math::

    \mathbf{\hat{\tilde{p}}}\cdot\mathbf{\hat{y}}
    &=&
    -
    \sin(\omega-\Delta\omega)\cos\tilde{\theta}
    \cos\tilde{\theta}\sin\Delta\omega
    +
    \cos(\omega-\Delta\omega)\cos\Delta\omega
    -
    \sin(\omega-\Delta\omega)\sin\tilde{\theta}
    \sin\tilde{\theta}\sin\Delta\omega\\
    &=&
    -
    \sin(\omega-\Delta\omega)\sin\Delta\omega
    +
    \cos(\omega-\Delta\omega)\cos\Delta\omega\\
    &=&
    \cos\omega\\
    &=&
    \mathbf{\hat{p}}\cdot\mathbf{\hat{y}}

Finally:

.. math::

	\mathbf{\hat{\tilde{p}}}\cdot\mathbf{\hat{z}}
    &=&
    - \sin(\omega-\Delta\omega)\cos\tilde{\theta}
    \left(
        \cos\theta\sin\tilde{\theta}
        -
        \sin\theta\cos\tilde{\theta}\cos\Delta\omega
    \right)\\
    &&
    + \cos(\omega-\Delta\omega)\sin\theta\sin\Delta\omega\\
    &&
    + \sin(\omega-\Delta\omega)\sin\tilde{\theta}
    \left(
        \cos\theta\cos\tilde{\theta}
        +
        \sin\theta\sin\tilde{\theta}\cos\Delta\omega
    \right)\\
    &=&
    +
    \sin\theta\cos^2\tilde{\theta}\sin(\omega-\Delta\omega)\cos\Delta\omega
    +
    \cos(\omega-\Delta\omega)\sin\theta\sin\Delta\omega
    +
    \sin\theta\sin^2\tilde{\theta}
    \sin(\omega-\Delta\omega)\cos\Delta\omega\\
    &=&
    \sin\theta\sin(\omega-\Delta\omega)\cos\Delta\omega
    +
    \sin\theta\cos(\omega-\Delta\omega)\sin\Delta\omega\\
    &=&
    \sin\theta\sin\omega\\
    &=&
    \mathbf{\hat{p}}\cdot\mathbf{\hat{z}}
