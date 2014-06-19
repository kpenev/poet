Calculation of the Pm,s Coefficients {#InclinationEccentricity_pms1}
====================================

We need only \f$m=0\f$ and \f$m=\pm2\f$.
Clearly:
\f{eqnarray*}{
	p_{m,s}&=&a^3\int_0^{2\pi/\omega} 
		\frac{e^{-im\Delta \phi(t)}}{r^3(t)}e^{i s \omega t}dt\\
	&=& a^3\int_0^{2\pi/\omega} 
		e^{-im\phi_0}\frac{\cos(m\phi(t))-i\sin(m\phi(t))}{r^3(t)} 
		e^{i s \omega t}dt
\f}
For \f$m=0\f$:
\f{eqnarray*}{
	p_{0,s}&=& \int_{0}^{2\pi} 
		\frac{e^{i s (u-e\sin u)}}{\omega (1-e\cos u)^2} du\\
	&=& \frac{1}{\omega}
		\int_{0}^{2\pi} \frac{e^{i s (u-e\sin u)}} {(1-e\cos u)^2} du
\f}
From \f$1/(1-x)^2=\sum_{k=0}^\infty (k+1)x^k\f$:
\f{eqnarray*}{
	p_{0,s}&=& \sum_{k=0}^\infty \frac{(k+1)e^k}{\omega}
		\int_{0}^{2\pi} e^{i s (u-e\sin u)} \cos^k u du\\
	&=&\sum_{k=0}^\infty \frac{(k+1)e^k}{2^k\omega}
		\int_{0}^{2\pi} e^{i s (u-e\sin u)}
			\left(e^{iu}+e^{-iu}\right)^k du\\
	&=&\sum_{k=0}^\infty \frac{(k+1)e^k}{2^k\omega}\sum_{c=0}^k
		{k \choose c} \int_{0}^{2\pi} e^{i s (u-e\sin u)} 
									e^{icu}e^{-i(k-c)u} du\\
	&=&\sum_{k=0}^\infty \frac{(k+1)e^k}{2^k\omega}\sum_{c=0}^k
		{k \choose c} \int_{0}^{2\pi} e^{i (s+2c-k) u} e^{-ies\sin u} du
\f}
If we change variable \f$u=u'-\pi/2\Rightarrow \sin u = -\cos u'\f$:
\f{eqnarray*}{
	p_{0,s}&=&\sum_{k=0}^\infty \frac{(k+1)e^k}{2^k\omega}
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
\f}
If we now change indices to
 \f$n=\lambda+c\Rightarrow \lambda=n-c,\quad\lambda+s+2c-k=n+s+c-k\f$. The
lower limit on \f$\lambda\f$ gives:
 \f$k-s-2c\leq n-c\Rightarrow c\geq k-n-s\f$. Finally, in order for the range
of \f$c\f$ to not be empty: \f$k-n-s\leq n\Rightarrow k\leq 2n+s\f$. With all
these we can write:
\f[
	p_{0,s}=\sum_{n=0}^\infty \alpha_{s,n}(se/2)^{s+2n}
\f]
with
\f[
	\alpha_{s,n}\equiv\frac{2\pi}{\omega}
		(-1)^n
		\sum_{k=0}^{2n+s} \frac{k+1}{s^k} \sum_{c=max(0,k-n-s)}^{min(n,k)}
		{k \choose c} \frac{(-1)^c}{(n-c)!(n+s+c-k)!}
\f]

For \f$m=\pm2\f$ we need:
\f{eqnarray*}{
	\cos2\phi &=& 1-2\sin^2\phi = 1-2\frac{(1-e^2)\sin^2u}{(1-e\cos u)^2} 
				= 1-\frac{(1-e^2)(1-\cos2u)}{(1-e\cos u)^2}\\
	\sin2\phi &=& 2\sin\phi\cos\phi = 2\sqrt{1-e^2}
				\frac{\sin u(\cos u - e)}{(1-e\cos u)^2}
				= 2\sqrt{1-e^2}\frac{\sin 2u - e\sin u}{(1-e\cos u)^2}
\f}
Plugging into the expression for \f$p_{\pm2,s}\f$:
\f{eqnarray*}{
	p_{\pm2,s}&=& \int_0^{2\pi/\omega} 
		a^3\exp(\mp 2i\phi_0)\frac{\cos(2\phi(t))\mp i\sin(2\phi(t))}{r^3(t)}
		\exp[i s (u-e\sin u)]dt\\
	&=& \frac{e^{\mp 2i\phi_0}}{\omega}\int_0^{2\pi} 
		\left[1-\frac{(1-e^2)(1-\cos2u)}{(1-e\cos u)^2}
				\mp
				i\sqrt{1-e^2}\frac{\sin 2u - 2e\sin u}{(1-e\cos u)^2}\right]
		\frac{\exp[i s (u-e\sin u)]}{(1-e\cos u)^2} du
\f}
Thus we need to evaluate 5 different integrals, the first of which was
already done while calculating \f$p_{0,s}\f$:
\f{eqnarray*}{
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
\f}
To solve them we will use 
 \f$1/(1-x)^4=\sum_{k=0}^\infty {{k+3} \choose 3} x^k\f$ and we will directly
calculate the following general integral:
\f{eqnarray*}{
	\omega I_{\lambda,s}&\equiv&
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
		{{k+3} \choose 3} \sum_{c=0}^k {k \choose c}
		\sum_{\nu=max(0,k-s-\lambda-2c)}^{\infty}
			\frac{(-1)^\nu (s^2e^2/4)^{\nu+c}}{\nu!(\nu+s+\lambda+2c-k)!}
\f}
Similarly to before we would like to group by powers of the eccentricity:
\f$n=\nu+c\f$. This leads to the following constraints:
\f{eqnarray*}{
	\nu>=0 & \Rightarrow & c \le n\\
	\nu>=k-s-\lambda-2c & \Rightarrow & c \ge k-s-\lambda-n\\
	k-s-\lambda-n \le n & \Rightarrow & k \le 2n+\lambda+s\\
	k-s-\lambda-n \le k & \Rightarrow & n \ge -s - \lambda
\f}
Plugging into the expression above:
\f[
	I_{\lambda,s} = \sum_{n=\max{0,-s-\lambda}}^\infty \beta_{s+\lambda,n}
		(se/2)^{s+\lambda+2n}
\f]

with
\f[
	\beta_{\lambda,n}\equiv \frac{2\pi(-1)^n}{\omega} \sum_{k=0}^{2n+\lambda}
		{{k+3} \choose 3} \sum_{c=\max(0,k-\lambda-n)}^{\min(n,k)}
			{k \choose c} \frac{(-1)^c}{(n-c)!(n+\lambda+c-k)!}
\f]

In terms of \f$I_{\lambda,s}\f$:
\f[
	p_{\pm2,s}=\exp(\mp 2i\phi_0)\left\{p_{0,s}
		+(1-e^2)\left[(I_{2,s}+I_{-2,s})/2-I_{0,s}\right]
		\mp \sqrt{1-e^2}(I_{2,s}-I_{-2,s})/2
		\pm e\sqrt{1-e^2}(I_{1,s}-I_{-1,s})
	\right\}
\f]
Plugging in the bessel function expressions:
\f{eqnarray*}{
	p_{\pm2,s}&=&\frac{\exp(\mp 2i\phi_0)}{\omega}\sum_{k=0}^\infty 
		2\pi \left(\frac{e}{2}\right)^k \sum_{c=0}^k {k \choose c}
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
\f}
