#!/usr/bin/python3 -u

from fractions import Fraction
from math import factorial

def choose(a, b) :
    """ Compute a choose b as an exact rational number (Fraction). """

    return Fraction(factorial(a), factorial(a-b)*factorial(b))

class ExpansionCoefficients :
    """ The eccentricity expansion coefficients of the p_{m,s} (see doc). """

    def __compute(self, s, n, l=None) :
        """ Computes alpha_{s,n} and beta_{s,n} (without the 2pi/omega). 
        
        See the documentation of __call__ for more details. """

        if s==0 :
            return Fraction(factorial(2*n+1),
                            factorial(n)**2*2**(2*n))
        else :
            l_use=(0 if l is None else l)
            assert(n+s+l_use>=0)
            result=Fraction(0, 1)
            for k in range(2*n+s+l_use+1) :
                cmin=max(0,k-n-s-l_use)
                sign=(-1 if cmin%2 else 1)
                c_sum=Fraction(0, 1)
                for c in range(cmin, min(n,k)+1) :
                    c_sum+=(choose(k,c)*sign/factorial(n-c)/
                            factorial(n+s+l_use+c-k))
                    sign=-sign
                result+=Fraction((k+1 if l is None else
                                  (k+1)*(k+2)*(k+3)//6), s**k)*c_sum
            return (-1 if n%2 else 1)*result

    def __init__(self, max_power) :
        """ Create an instance for evaluating alpha values. 
        Arguments:
            - max_power : The maximum power of the eccentricity that will
            ever be requested.
        Returns: None.
        """

        self.__alpha=[]
        self.__beta=[[], [], [], [], []]
        for epower in range(max_power+1) :
            self.__alpha.append([])
            for l in range(-2, 3) : self.__beta[l].append([])
            for n in range(epower+1) :
                s=epower-2*n
                self.__alpha[-1].append(self.__compute(s, n))
                for l in range(-2, 3) :
                    s=epower-2*n-l
                    self.__beta[l][-1].append(self.__compute(s, n, l))
    
    def alpha_or_beta(self, s, n, l=None) :
        """ Returns the value of alpha or beta for the given indices.

        Arguments:
            - s, n, l: The indices of the alpha/beta (if l is/is not None)
                       value being evaluated. The total power of the
                       eccentricity is s+2n+l, where l=None is taken to be
                       zero (see documentation).
        Returns:
            omega/(2pi) (alpha or beta)_{s,n} as an exact rational number.
        """

        l_use=(0 if l is None else l)
        if(n+s+l_use<0) : return Fraction(0, 1)
        epower=s+2*n+l_use
        assert(epower<len(self.__alpha))
        assert(n<len(self.__alpha[epower]))
        if l is None : return self.__alpha[epower][n]
        else : return self.__beta[l][epower][n]

    def __call__(self, m, s, epower) :
        """ The coefficient in front of 2*pi*e^epower/omega in p_{m,s}. 
        
        Arguments:
            - m : The first index of the p_{m,s} being evaluated. Can only be
                  (0 or +-2).
            - s : The second index of the p_{m,s} being evaluated.
            - epower : The power of eccentricity of the term being evaluated.
        Returns: The specified coefficient as an exact rational number. """

        if(m==0) : 
            if epower%2!=s%2 or abs(s)>epower : return 0
            n=(epower-s)//2
            return self.alpha_or_beta(s, n)*(1 if s==0 else 
                                             Fraction(s,2)**epower)
        else :
            raise Exception("Not implemented yet.")

if __name__=='__main__' :
    coef=ExpansionCoefficients(10)
    for m in [-2, 0, 2] :
        if m!=0 : continue
        for s in range(-10, 11) :
            for epower in range(11) :
                print(m, s, epower, coef(m,s,epower),
                      "%25.16e"%coef(m,s,epower))

