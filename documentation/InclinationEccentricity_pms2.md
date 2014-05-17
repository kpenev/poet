Calculation of the Pm,s Coefficients {#InclinationEccentricity_pms2}
=============================================

We need only \f$m=0\f$ and \f$m=\pm2\f$.
Clearly:
\f{eqnarray*}{
	p_{m,s}&=&\int_0^{2\pi/\omega} 
		\frac{e^{-im\Delta \phi(t)}}{r^3(t)}e^{i s \omega t}dt\\
	&=& \int_0^{2\pi/\omega} 
		e^{-im\phi_0}\frac{\cos(m\phi(t))-i\sin(m\phi(t))}{r^3(t)} 
		e^{i s \omega t}dt
\f}
For \f$m=0\f$:
\f{eqnarray*}{
	p_{0,s}&=& \int_{0}^{2\pi} 
		\frac{e^{i s (u-e\sin u)}}{\omega a^3(1-e\cos u)^2} du\\
	&=& \frac{1}{\omega a^3}
		\int_{0}^{2\pi} \frac{e^{i s (u-e\sin u)}} {(1-e\cos u)^2} du
\f}
From \f$1/(1-x)^2=\sum_{k=0}^\infty (k+1)x^k\f$:
\f{eqnarray*}{
	p_{0,s}&=&\int_0^{2\pi} 
		\left[\cos u + i \sin u\right]^s 
		\left[\sum_{l=0}^\infty \frac{(-ise)^l\sin^l u}{l!}\right]
		\left[\sum_{p=0}^\infty (p+1) e^p\cos^p u\right]\\
	&=&\int_0^{2\pi} 
		\left[\sum_{k=0}^s {s \choose k} i^s\cos^{s-k} u \sin^k u\right]
		\left[\sum_{l=0}^\infty \frac{(-ise)^l\sin^l u}{l!}\right]
		\left[\sum_{p=0}^\infty (p+1) e^p\cos^p u\right]\\
	&=& \sum_{k=0}^s \sum_{l=0}^\infty \sum_{p=0}^\infty 
			(-1)^l i^{s+l} (p+1) {s \choose k} \frac{s^l}{l!} e^{l+p}
			I_{k+l,s-k+p}
\f}
with
\f[
	I_{m,n} \equiv \int_0^{2\pi} \sin^m u \cos^n u
\f]
