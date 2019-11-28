************************************
Calculation of the Pm,s Coefficients
************************************

We need only :math:`m=0` and :math:`m=\pm2`.
Clearly:

.. math::

	p_{m,s}&=&\frac{a^3}{2\pi}\int_0^{2\pi/\omega} 
		\frac{e^{-im\Delta \phi(t)}}{r^3(t)}e^{i s \omega t}dt\\
	&=& a^3\int_0^{2\pi/\omega} 
		e^{-im\phi_0}\frac{\cos(m\phi(t))-i\sin(m\phi(t))}{r^3(t)} 
		e^{i s \omega t}dt

For :math:`m=0`:

.. math::

	p_{0,s}&=& \frac{1}{2\pi}\int_{0}^{2\pi} 
		\frac{e^{i s (u-e\sin u)}}{\omega (1-e\cos u)^2} du\\
	&=& \frac{1}{2\pi\omega}
		\int_{0}^{2\pi} \frac{e^{i s (u-e\sin u)}} {(1-e\cos u)^2} du

From :math:`1/(1-x)^2=\sum_{k=0}^\infty (k+1)x^k`:

.. math::

	p_{0,s}= \sum_{k=0}^\infty \frac{(k+1)e^k}{2\pi\omega}
		\int_{0}^{2\pi} e^{i s (u-e\sin u)} \cos^k u du

Which for :math:`s=0`, using :math:`\int_0^{2\pi} \cos^{2k} u du = \frac{2\pi
(2k)!}{2^{2k}(k!)^2}` gives:

.. math::

	p_{0,0}=\frac{1}{\omega}
		\sum_{k=0}^\infty \frac{(2k+1)!}{2^{2k}(k!)^2}e^k

And for :math:`s \neq 0`\ :

.. math::

	p_{0,s} & = & \sum_{k=0}^\infty \frac{(k+1)e^k}{2\pi\omega}
		\int_{0}^{2\pi} e^{i s (u-e\sin u)} \cos^k u du\\
	& = & \sum_{k=0}^\infty \frac{(k+1)e^k}{2^{k+1}\pi\omega}
		\int_{0}^{2\pi} e^{i s (u-e\sin u)}
			\left(e^{iu}+e^{-iu}\right)^k du\\
	& = & \sum_{k=0}^\infty \frac{(k+1)e^k}{2^{k+1}\pi\omega}\sum_{c=0}^k
		{k \choose c} \int_{0}^{2\pi} e^{i s (u-e\sin u)} 
									e^{icu}e^{-i(k-c)u} du\\
	& = & \sum_{k=0}^\infty \frac{(k+1)e^k}{2^{k+1}\pi\omega}\sum_{c=0}^k
		{k \choose c} \int_{0}^{2\pi} e^{i (s+2c-k) u} e^{-ies\sin u} du

If we change variable :math:`u=u'-\pi/2\Rightarrow \sin u = -\cos u'`\ :

.. math::

	2\pi p_{0,s}&=&\sum_{k=0}^\infty \frac{(k+1)e^k}{2^k\omega}
		\sum_{c=0}^k {k \choose c} e^{-i (s+2c-k)\pi/2}
			\int_{\pi/2}^{5\pi/2} e^{i (s+2c-k) u'} e^{ies\cos u'} du'\\
	&=&\sum_{k=0}^\infty \frac{(k+1)e^k}{2^k\omega}
		\sum_{c=0}^k {k \choose c} (-i)^{s+2c-k} \left\{
			\int_{\pi/2}^{2\pi} e^{i (s+2c-k) u'} e^{ies\cos u'} du'
			+
			\int_{2\pi}^{5\pi/2} e^{i (s+2c-k) u'} e^{ies\cos u'} du' 
		\right\}\\
	&=&\sum_{k=0}^\infty \frac{(k+1)e^k}{2^k\omega}
		\sum_{c=0}^k {k \choose c} (-i)^{s+2c-k} \left\{
			\int_{\pi/2}^{2\pi} e^{i (s+2c-k) u'} e^{ies\cos u'} du'
			+
			\int_{0}^{\pi/2} e^{i (s+2c-k) u'} e^{ies\cos u'} du' \right\}\\
	&=&\sum_{k=0}^\infty \frac{(k+1)e^k}{2^k\omega}
		\sum_{c=0}^k {k \choose c} (-i)^{s+2c-k}
			\int_{0}^{2\pi} e^{i (s+2c-k) u'} e^{ies\cos u'} du'\\
	&=&\sum_{k=0}^\infty \frac{2\pi(k+1)e^k}{2^k\omega}
		\sum_{c=0}^k {k \choose c} (-i)^{s+2c-k} i^{s+2c-k}
			J_{s+2c-k}(es)\\
	&=&\sum_{k=0}^\infty \frac{2\pi(k+1)e^k}{2^k\omega}
		\sum_{c=0}^k {k \choose c} J_{s+2c-k}(es)\\
	&=&\sum_{k=0}^\infty \frac{2\pi(k+1)(se)^k}{2^k\omega s^k}
		\sum_{c=0}^k {k \choose c} \sum_{\lambda=max(0,k-s-2c)}^{\infty}
			\frac{(-1)^\lambda (se)^{2\lambda+s+2c-k}}
				{2^{2\lambda+s+2c-k}\lambda!(\lambda+s+2c-k)!}\\
	&=&\sum_{k=0}^\infty \frac{2\pi(k+1)(se)^k}{2^k\omega s^k}
		\sum_{c=0}^k {k \choose c} \sum_{\lambda=max(0,k-s-2c)}^{\infty}
			\frac{(-1)^\lambda (se)^{2\lambda+s+2c-k}}
				{2^{2\lambda+s+2c-k}\lambda!(\lambda+s+2c-k)!}\\
	&=&\frac{2\pi}{\omega} \left(\frac{es}{2}\right)^s
		\sum_{k=0}^\infty \frac{(k+1)}{s^k}
		\sum_{c=0}^k {k \choose c} \sum_{\lambda=max(0,k-s-2c)}^{\infty}
			\frac{(-1)^\lambda (s^2e^2/4)^{\lambda+c}}
				{\lambda!(\lambda+s+2c-k)!}\\

If we now change indices to :math:`n=\lambda+c\Rightarrow
\lambda=n-c,\quad\lambda+s+2c-k=n+s+c-k`\ . The lower limit on :math:`\lambda`
gives: :math:`k-s-2c\leq n-c\Rightarrow c\geq k-n-s`\ . Finally, in order for the
range of :math:`c` to not be empty: :math:`k-n-s\leq n\Rightarrow k\leq 2n+s`. With
all these we can write:

.. math::

	p_{0,s}=\sum_{n=0}^\infty \alpha_{s,n}\left\{
		\begin{array}{l@{,\quad}l}
			e^{2n} & s=0\\
			(se/2)^{s+2n} & s \neq 0
		\end{array}\right.

with

.. math::

	\alpha_{s,n}\equiv\frac{1}{\omega}\left\{\begin{array}{l@{,\quad}l}
		\frac{(2n+1)!}{2^{2n}(n!)^2} & s=0\\
		(-1)^n \sum_{k=0}^{2n+s} \frac{k+1}{s^k}
		\sum_{c=max(0,k-n-s)}^{min(n,k)}
		{k \choose c} \frac{(-1)^c}{(n-c)!(n+s+c-k)!} & s \neq 0
	\end{array} \right.

Verified using Mathematica.

For :math:`m=\pm2` we need:

.. math::

	\cos2\phi &=& 1-2\sin^2\phi = 1-2\frac{(1-e^2)\sin^2u}{(1-e\cos u)^2} 
				= 1-\frac{(1-e^2)(1-\cos2u)}{(1-e\cos u)^2}\\
	\sin2\phi &=& 2\sin\phi\cos\phi = 2\sqrt{1-e^2}
				\frac{\sin u(\cos u - e)}{(1-e\cos u)^2}
				= 2\sqrt{1-e^2}\frac{\sin 2u - e\sin u}{(1-e\cos u)^2}

Plugging into the expression for :math:`p_{\pm2,s}`\ :

.. math::

	p_{\pm2,s}&=& \frac{1}{2\pi}\int_0^{2\pi/\omega} 
		a^3\exp(\mp 2i\phi_0)\frac{\cos(2\phi(t))\mp i\sin(2\phi(t))}{r^3(t)}
		\exp[i s (u-e\sin u)]dt\\
	&=& \frac{e^{\mp 2i\phi_0}}{2\pi\omega}\int_0^{2\pi} 
		\left[1-\frac{(1-e^2)(1-\cos2u)}{(1-e\cos u)^2}
				\mp
				i\sqrt{1-e^2}\frac{\sin 2u - 2e\sin u}{(1-e\cos u)^2}\right]
		\frac{\exp[i s (u-e\sin u)]}{(1-e\cos u)^2} du

Thus we need to evaluate 5 different integrals, the first of which was
already done while calculating :math:`p_{0,s}`\ :

.. math::

	p_{\pm2,s}&=&\exp\left(\mp 2i\phi_0\right) p_{0,s} -\\
		&&{}-\frac{\exp(\mp 2i\phi_0)(1-e^2)}{2\pi\omega}
		\int_{0}^{2\pi} \frac{\exp[i s (u-e\sin u)]} {(1-e\cos u)^4} du+\\
		&&{}+\frac{\exp(\mp 2i\phi_0)(1-e^2)}{2\pi\omega}
		\int_{0}^{2\pi} \frac{\exp[i s (u-e\sin u)]\cos 2u} {(1-e\cos u)^4}
		du\mp\\
		&&{}\mp i\frac{\exp(\mp 2i\phi_0)\sqrt{1-e^2}}{2\pi\omega}
		\int_{0}^{2\pi} \frac{\exp[i s (u-e\sin u)]\sin 2u} {(1-e\cos u)^4}
		du\pm\\
		&&{}\pm i\frac{2e\exp(\mp 2i\phi_0)\sqrt{1-e^2}}{2\pi\omega}
		\int_{0}^{2\pi} \frac{\exp[i s (u-e\sin u)]\sin u} {(1-e\cos u)^4}
		du\\

To solve them we will use :math:`1/(1-x)^4=\sum_{k=0}^\infty {{k+3} \choose 3}
x^k` and we will directly calculate the following general integral:


.. math::

	2\pi\omega I_{\lambda,s}&\equiv&
	\int_{0}^{2\pi} \frac{\exp[i s (u-e\sin u)]\exp(i\lambda u)}
						{(1-e\cos u)^4} du\\
	&=&\sum_{k=0}^\infty {{k+3} \choose 3} e^k
	\int_{0}^{2\pi} \exp[i (s+\lambda) u]\exp[-ise\sin u)]\cos^k u du\\
	&=&\sum_{k=0}^\infty {{k+3} \choose 3} \left(\frac{e}{2}\right)^k
	\sum_{c=0}^k {k \choose c} \int_{0}^{2\pi} 
						\exp[i (s+\lambda+2c-k) u]\exp(-ise\sin u) du\\
	&=&\sum_{k=0}^\infty {{k+3} \choose 3} 2\pi \left(\frac{e}{2}\right)^k
		\sum_{c=0}^k {k \choose c} J_{s+\lambda+2c-k}(es)\\
	&=&\sum_{k=0}^\infty {{k+3} \choose 3} 2\pi \left(\frac{e}{2}\right)^k
		\sum_{c=0}^k {k \choose c} \sum_{\nu=max(0,k-s-\lambda-2c)}^{\infty}
			\frac{(-1)^\nu (se)^{2\nu+s+\lambda+2c-k}}
				{2^{2\nu+s+\lambda+2c-k}\nu!(\nu+s+\lambda+2c-k)!}\\
	&=&2\pi\left(\frac{se}{2}\right)^{s+\lambda}\sum_{k=0}^\infty
		{{k+3} \choose 3} s^{-k}\sum_{c=0}^k {k \choose c}
		\sum_{\nu=max(0,k-s-\lambda-2c)}^{\infty}
			\frac{(-1)^\nu (s^2e^2/4)^{\nu+c}}{\nu!(\nu+s+\lambda+2c-k)!}

Similarly to before we would like to group by powers of the eccentricity:
:math:`n=\nu+c`\ . This leads to the following constraints:

.. math::

	\nu>=0 & \Rightarrow & c \le n\\
	\nu>=k-s-\lambda-2c & \Rightarrow & c \ge k-s-\lambda-n\\
	k-s-\lambda-n \le n & \Rightarrow & k \le 2n+\lambda+s\\
	k-s-\lambda-n \le k & \Rightarrow & n \ge -s - \lambda

Plugging into the expression above:

.. math::

	I_{\lambda,s} = \sum_{n=\max(0,-s-\lambda)}^\infty \beta_{\lambda,s,n}
		(se/2)^{s+\lambda+2n}

with

.. math::

	\beta_{\lambda,s,n}\equiv \frac{(-1)^n}{\omega}
	\sum_{k=0}^{2n+\lambda+s}
		{{k+3} \choose 3} s^{-k}\sum_{c=\max(0,k-\lambda-s-n)}^{\min(n,k)}
			{k \choose c} \frac{(-1)^c}{(n-c)!(n+\lambda+s+c-k)!}

In terms of :math:`I_{\lambda,s}`\ :

.. math::

	p_{\pm2,s}=\exp(\mp 2i\phi_0)\left\{p_{0,s}
		+(1-e^2)\left[(I_{2,s}+I_{-2,s})/2-I_{0,s}\right]
		\mp \sqrt{1-e^2}(I_{2,s}-I_{-2,s})/2
		\pm e\sqrt{1-e^2}(I_{1,s}-I_{-1,s})
	\right\}

Verified using Methematica for :math:`s\neq0`\ .

Using:

.. math::

	\sqrt{1-e^2}=\sum_{n=0}^\infty \frac{(2n)!}{4^n (n!)^2(1-2n)} e^{2n}

we can rewrite:

.. math::

	p_{\pm2,s}=\exp(\mp 2i\phi_0)\sum_{n=-1}^\infty 
		\gamma^\pm_{s,n}\left(\frac{se}{2}\right)^{2n+s}

with:

.. math::

	\gamma^\pm_{s,n} &\equiv &\alpha_{s,n}
	+
	\frac{\beta_{2,s,n-1}+\beta_{-2,s,n+1}}{2}
	-
	\beta_{0,s,n}+\frac{4}{s^2}\beta_{0,s,n-1}
	-
	\frac{2}{s^2}\left(\beta_{2,s,n-2}+\beta_{-2,s,n}\right)\\
	&&{}\pm
	\sum_{k=0}^{n+1} \frac{(2k)!}{s^{2k}(k!)^2(2k-1)}
		\left[\frac{1}{2}\left(\beta_{2,s,n-k-1}-\beta_{-2,s,n-k+1}\right)+
		\frac{2}{s}\left(\beta_{-1,s,n-k}-\beta_{1,s,n-k-1}\right)\right]

Verified by Mathematica.

Plugging in the bessel function expressions:

.. math::

	p_{\pm2,s}&=&\frac{\exp(\mp 2i\phi_0)}{\omega}\sum_{k=0}^\infty 
		\left(\frac{e}{2}\right)^k \sum_{c=0}^k {k \choose c}
		\Bigg\{
			(k+1)J_{s+2c-k}(es) + \\
		&&{}+{{k+3} \choose 3} \Bigg[
			-(1-e^2)J_{s+2c-k}(es)
			+\frac{1-e^2}{2}\big[J_{s+2+2c-k}(es)+J_{s-2+2c-k}(es)\big]
			\mp\\
		&&\quad\quad\quad\quad\quad{}
			\mp\frac{\sqrt{1-e^2}}{2}\big[J_{s+2+2c-k}(es)-
											J_{s-2+2c-k}(es)\big]
			\pm e\sqrt{1-e^2}\big[J_{s+1+2c-k}(es)-J_{s-1+2c-k}(es)\big]
		\Bigg]
		\Bigg\}

For s=0 we need to go back to:

.. math::

	2\pi p_{\pm2,0}&=&\exp\left(\mp 2i\phi_0\right) \left\{2\pi p_{0,0} +
		\frac{1}{\omega}\left[
			(1-e^2)\int_{0}^{2\pi} \frac{\cos 2u -1} {(1-e\cos u)^4} du
			\mp
			i\sqrt{1-e^2}\int_{0}^{2\pi} \frac{\sin 2u} {(1-e\cos u)^4} du
			\pm
			i2e\sqrt{1-e^2} \int_{0}^{2\pi} \frac{\sin u} {(1-e\cos u)^4} du
		\right]\right\}\\
	&=&\exp\left(\mp 2i\phi_0\right) \left\{2\pi p_{0,0} +
		\frac{2}{\omega}\left[
			(1-e^2)\int_{0}^{2\pi} \frac{\cos^2 u -1} {(1-e\cos u)^4} du
			\pm
			i\sqrt{1-e^2}\int_{0}^{2\pi} \frac{\cos u}{(1-e\cos u)^4}d\cos u
			\mp
			ie\sqrt{1-e^2} \int_{0}^{2\pi} \frac{1}{(1-e\cos u)^4} d\cos u
		\right]\right\} 

.. math::

	\int_{0}^{2\pi}\frac{\cos^{2n} u} {(1-e\cos u)^4} du
	&=&\sum_{k=0}^\infty {2k+3 \choose 3} e^{2k}
			\int_{0}^{2\pi}\frac{\cos^{2n+2k} u} du\\
	&=&2\pi\sum_{k=0}^\infty {2k+3 \choose 3}
			\frac{(2k+2n)!}{2^{2k+2n}[(k+n)!]^2} e^{2k}

.. math::

	\int_{0}^{2\pi} \frac{1}{(1-e\cos u)^4} d\cos u
	&=&-\frac{1}{e}\int_{0}^{2\pi} \frac{1}{(1-e\cos u)^4} d(1-e\cos u)\\
	&=&\left.\frac{1}{3e(1-e\cos u)^3}\right|_{0}^{2\pi}\\
	&=&0

.. math::

	\int_{0}^{2\pi} \frac{\cos u}{(1-e\cos u)^4}d\cos u
	&=&\frac{1}{e^2}\int_{0}^{2\pi} \frac{1-e\cos u-1}{(1-e\cos u)^4}
		d(1-e\cos u)\\
	&=&\frac{1}{e^2}\int_{0}^{2\pi} \frac{1}{(1-e\cos u)^3} d(1-e\cos u)\\
	&=&-\frac{1}{2e^2(1-e\cos u)^2}\\
	&=&0

So we are left with:

.. math::

	p_{\pm2,0}&=&\exp\left(\mp 2i\phi_0\right) \left\{p_{0,0} +
		\frac{2(1-e^2)}{\omega}\left\{\sum_{k=0}^\infty {2k+3 \choose 3}
			\frac{(2k+2)!}{2^{2k+2}[(k+1)!]^2}-
			\frac{2k!}{2^{2k}(k!)^2}\right\}e^{2k}\right\}\\
	&=&\exp\left(\mp 2i\phi_0\right) \left\{p_{0,0} -
		\frac{2(1-e^2)}{\omega}\sum_{k=0}^\infty {2k+3 \choose 3}
			\frac{2k!}{2^{2k}(k!)^2(2k+2)}e^{2k}\right\}\\
	&=&\exp\left(\mp 2i\phi_0\right) \left\{p_{0,0} -
		\frac{2}{\omega}\sum_{k=0}^\infty \left[
			{2k+3 \choose 3} \frac{2k!}{2^{2k}(k!)^2(2k+2)}
			-
			{2k+1 \choose 3} \frac{2(k-1)!}{2^{2k-2}[(k-1)!]^2 2k}
		\right]e^{2k}\right\}\\
	&=&\exp\left(\mp 2i\phi_0\right) \left\{p_{0,0} -
		\frac{2}{\omega}\sum_{k=0}^\infty \left[
			\frac{(2k+3)(2k+1)!}{6\,2^{2k}(k!)^2}
			-
			\frac{4k^2(2k+1)(2k-1)!}{6\,2^2k(k!)^2}
		\right]e^{2k}\right\}\\
	&=&\exp\left(\mp 2i\phi_0\right) \left\{p_{0,0} -
		\frac{2}{\omega}\sum_{k=0}^\infty \left[
			\frac{(2k+3)(2k+1)!}{6\,2^{2k}(k!)^2}
			-
			\frac{2k(2k+1)!}{6\,2^2k(k!)^2}
		\right]e^{2k}\right\}\\
	&=&\exp\left(\mp 2i\phi_0\right) \left\{p_{0,0} -
		\frac{2}{\omega}\sum_{k=0}^\infty \frac{(2k+1)!}{2^{2k+1}(k!)^2}
		e^{2k}\right\}\\
	&=&0

Confirmed by Mathematica.
