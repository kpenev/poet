************************************
Calculation of the Pm,s Coefficients
************************************

We need only :math:`m=0` and :math:`m=\pm2`.
Clearly:

.. math::

	p_{m,s}&=&\int_0^{2\pi/\omega} 
		\frac{a^3\exp(-im\Delta \phi(t))}{r^3(t)}e^{i s \omega t}dt\\
	&=& \int_0^{2\pi/\omega} 
		a^3\exp(-im\phi_0)\frac{\cos(m\phi(t))-i\sin(m\phi(t))}{r^3(t)} 
		e^{i s \omega t}dt

For :math:`m=0`:

.. math::

	p_{0,s}&=& \int_{0}^{2\pi} 
		\frac{\exp(i s (u-e\sin u))}{\omega (1-e\cos u)^2} du\\
	&=& \frac{1}{\omega}
		\int_{0}^{2\pi} \frac{\exp(i s (u-e\sin u))} {(1-e\cos u)^2} du

From :math:`1/(1-x)^2=\sum_{k=0}^\infty (k+1)x^k`:

.. math::

	p_{0,s}&=&\frac{1}{\omega}\int_0^{2\pi} 
		\left[\cos u + i \sin u\right]^s 
		\left[\sum_{l=0}^\infty \frac{(-ise)^l\sin^l u}{l!}\right]
		\left[\sum_{p=0}^\infty (p+1) e^p\cos^p u\right]\\
	&=&\frac{1}{\omega}\int_0^{2\pi} 
		\left[\sum_{k=0}^s {s \choose k} i^k\cos^{s-k} u \sin^k u\right]
		\left[\sum_{l=0}^\infty \frac{(-ise)^l\sin^l u}{l!}\right]
		\left[\sum_{p=0}^\infty (p+1) e^p\cos^p u\right]\\
	&=&\frac{1}{\omega}\sum_{k=0}^s \sum_{l=0}^\infty \sum_{p=0}^\infty 
			(-1)^l i^{k+l} (p+1) {s \choose k} \frac{s^l}{l!} e^{l+p}
			I_{k+l,s-k+p}\\
	&=&\frac{1}{\omega}\sum_{n=0}^\infty \left[\sum_{k=0}^s {s \choose k}
			\sum_{l=0}^n (-1)^l i^{k+l} (n-l+1) \frac{s^l}{l!}
			I_{k+l,s-k+n-l}\right]e^n
			

with

.. math::

	I_{m,n} \equiv \int_0^{2\pi} \sin^m u \cos^n u

To find an expression for :math:`m=\pm2` consider:

.. math::

	Q_{p,q}&\equiv&
		\int_{0}^{2\pi} \frac{\exp[i s (u-e\sin u)]} {(1-e\cos u)^4} \sin^p u
					    \cos^q udu\\
		&=& \sum_{\lambda=0}^\infty \frac{(-ise)^\lambda}{\lambda!}
			\int_{0}^{2\pi} 
				\frac{\sin^{\lambda+p} u \cos^q u 
					\left(\cos u + i\sin u\right)^s}{(1-e\cos u)^4}
								
Using: :math:`1/(1-x)^4=\sum_{k=0}^\infty {{k+3} \choose 3} x^k`

.. math::
	Q_{p,q}&=&\sum_{\lambda=0}^\infty \frac{(-ise)^\lambda}{\lambda!}
			\sum_{k=0}^\infty {k+3 \choose 3} e^k
			\sum_{\sigma=0}^{s} {s \choose \sigma} i^\sigma 
				I_{\lambda+\sigma+p, s-\sigma+k+q}\\
	Q_{p,q}&=&\sum_{n=0}^\infty \left[
			\sum_{\lambda=0}^n \frac{(-is)^\lambda}{\lambda!}
			{n-\lambda+3 \choose 3} 
			\sum_{\sigma=0}^{s} {s \choose \sigma} i^\sigma 
				I_{\lambda+\sigma+p, s+n+q-\sigma-\lambda}\right] e^n

:doc:`here <inclination_eccentricity_pms1>` we show:

.. math::
	p_{\pm2,s}&=&\exp\left(\mp 2i\phi_0\right) p_{0,s} -\\
		&&{}-\frac{\exp(\mp 2i\phi_0)(1-e^2)}{\omega}
		\int_{0}^{2\pi} \frac{\exp[i s (u-e\sin u)]} {(1-e\cos u)^4} du+\\
		&&{}+\frac{\exp(\mp 2i\phi_0)(1-e^2)}{\omega}
		\int_{0}^{2\pi} \frac{\exp[i s (u-e\sin u)]\cos 2u} {(1-e\cos u)^4}
		du\mp\\
		&&{}\mp i\frac{\exp(\mp 2i\phi_0)\sqrt{1-e^2}}{\omega}
		\int_{0}^{2\pi} \frac{\exp[i s (u-e\sin u)]\sin 2u} {(1-e\cos u)^4}
		du\pm\\
		&&{}\pm i\frac{2e\exp(\mp 2i\phi_0)\sqrt{1-e^2}}{\omega}
		\int_{0}^{2\pi} \frac{\exp[i s (u-e\sin u)]\sin u} {(1-e\cos u)^4}
		du\\
	&=&\exp\left(\mp 2i\phi_0\right) p_{0,s} -\\
		&&{}-\frac{\exp(\mp 2i\phi_0)(1-e^2)}{\omega}Q_{0,0}+\\
		&&{}+\frac{\exp(\mp 2i\phi_0)(1-e^2)}{\omega}(Q_{0,0}-2Q_{2,0})\mp\\
		&&{}\mp i\frac{\exp(\mp 2i\phi_0)\sqrt{1-e^2}}{\omega}2Q_{1,1}\pm\\
		&&{}\pm i\frac{2e\exp(\mp 2i\phi_0)\sqrt{1-e^2}}{\omega}Q_{1,0}\\
	&=&\exp\left(\mp 2i\phi_0\right) \left[p_{0,s} -
		2\frac{(1-e^2)Q_{2,0} \pm i\sqrt{1-e^2}(Q_{1,1} - e Q_{1,0})}{\omega}
		\right]

According to `this page <http://planetmath.org/taylorexpansionofsqrt1x>`_

.. math::

	\sqrt{1-e^2}=1-\frac{e^2}{2}-
		\sum_{n=2}^\infty\frac{(2n-3)!}{2^{2n-2} n! (n-2)!} e^{2n}

Also note that:

.. math::

	I_{m,n}&=&I_{m,n-2}-I_{m+2,n-2}=I_{m,n-2}-\frac{m+1}{n-1}I_{m,n}\\
	\Rightarrow I_{m,n}&=&\frac{n-1}{n+m}I_{m,n-2}

Similarly:

.. math::

	I_{m,n}&=&I_{m-2,n}-I_{m-2,n+2}=I_{m-2,n}-\frac{n+1}{m-1}I_{m,n}\\
	\Rightarrow I_{m,n}&=&\frac{m-1}{n+m}I_{m-2,n}

Further :math:`I_{m,1}=I_{1,n}=0`, so :math:`I_{m,n}` is non-zero only if both m
and n are even.

Finally, :math:`I_{m,n}=I_{n,m}`.
