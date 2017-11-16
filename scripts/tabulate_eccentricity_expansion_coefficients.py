#!/usr/bin/python3 -u

from fractions import Fraction
from math import factorial
from os.path import exists
from os import rename
try:
    from multiprocessing.pool import Pool, cpu_count
except:
    from multiprocessing.pool import Pool
    from multiprocessing import cpu_count

from argparse import ArgumentParser
from math import pi

def choose(a, b) :
    """ Compute a choose b as an exact rational number (Fraction). """

    return Fraction(factorial(a), factorial(a-b)*factorial(b))

def next_line(f) :
    """ Returns the next non-comment line from the file. """

    line='#'
    while line[0]=='#' : line=f.readline()
    return line

def parse_fraction(frac_str) :
    """ Parses the given fraction string into a Fraction. """

    return Fraction(*map(eval, frac_str.split('/')))

def compute_alpha_beta(task) :
    """ Computes the alpha} and beta coefficients without the 2pi/omega.

    Arguments:
        - task: a tuple of (s, n[, l]) values. If l is omitted computes
                alpha_{s,n}, else computes beta_{s,n,l}
    Returns: the computed value.
    """

    if len(task)==2 :
        s, n=task
        l=None
    else : s, n, l=task
    if s==0 :
        return Fraction(factorial(2*n+1),
                        factorial(n)**2*2**(2*n))
    else :
        l_use=(0 if l is None else l)
        if n+s+l_use<0 : return Fraction(0, 1)
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

class ExpansionCoefficients :
    """ The eccentricity expansion coefficients of the p_{m,s} (see doc). """

    def __read_alpha_beta_file(self, alpha_beta_fname) :
        """ Reads alpha and beta values from the given filename. """

        self.__alpha=[]
        self.__beta=[[], [], [], [], []]
        if not exists(alpha_beta_fname) : return
        alpha_beta_f=open(alpha_beta_fname, 'r')
        line=next_line(alpha_beta_f)
        assert(line.startswith('max e power='))
        self.max_e_power=eval(line[len('max e power='):].strip())
        for epower in range(self.max_e_power+1) :
            self.__alpha.append([])
            for l in range(-2, 3) : self.__beta[l].append([])
            for n in range(epower+1) :
                s=epower-2*n
                line=next_line(alpha_beta_f)
                values=map(parse_fraction, line.split())
                self.__alpha[-1].append(next(values))
                for l in range(-2, 3) :
                    s=epower-2*n-l
                    self.__beta[l][-1].append(next(values))

    def __write_alpha_beta_file(self, alpha_beta_fname) :
        """ Writes the current alpha and beta values to the given file. """

        f=open(alpha_beta_fname+'.temporary', 'w')
        f.write('max e power='+str(self.max_e_power)+'\n')
        for epower in range(self.max_e_power+1) :
            for n in range(epower+1) :
                s=epower-2*n
                line=str(self.__alpha[epower][n])
                for l in range(-2, 3) :
                    line+=' '+str(self.__beta[l][epower][n])
                f.write(line+'\n')
        f.close()
        rename(alpha_beta_fname+'.temporary', alpha_beta_fname)

    def __init__(self, max_power, alpha_beta_fname,
                 num_processes, chunksize=100) :
        """ Create an instance for evaluating alpha values. 
        Arguments:
            - max_power : The maximum power of the eccentricity that will
                          ever be requested.
            - alpha_beta_fname : The filename to read previously computed
                                 alpha/beta values from and store any newly 
                                 calculated ones into.
        Returns: None.
        """

        self.max_e_power=-1
        self.__read_alpha_beta_file(alpha_beta_fname)
        while True :
            power_step=min(self.max_e_power+11, max_power+1)
            to_compute=[]
            for epower in range(self.max_e_power+1, power_step) :
                for n in range(epower+1) :
                    s=epower-2*n
                    to_compute.append((s,n))
                    for l in range(-2, 3) :
                        to_compute.append((s-l, n, l))
            if not to_compute : return
            compute_pool=Pool(processes=num_processes)
            computed=iter(compute_pool.map(compute_alpha_beta, to_compute,
                                           chunksize))
            for epower in range(self.max_e_power+1, power_step) :
                self.__alpha.append([])
                for l in range(-2, 3) : self.__beta[l].append([])
                for n in range(epower+1) :
                    s=epower-2*n
                    self.__alpha[-1].append(next(computed))
                    for l in range(-2, 3) :
                        s=epower-2*n-l
                        self.__beta[l][-1].append(next(computed))
            self.max_e_power=power_step-1
            self.__write_alpha_beta_file(alpha_beta_fname)
            print('max_e_power=', self.max_e_power)
            if power_step==max_power+1 : return
    
    def alpha_or_beta(self, s, n, l=None) :
        """ Returns the value of alpha or beta for the given indices.

        Arguments:
            - s, n, l: The indices of the alpha/beta (if l is/is not None)
                       value being evaluated. The total power of the
                       eccentricity is s+2n+l, where l=None is taken to be
                       zero (see documentation).
        Returns:
            omega (alpha or beta)_{s,n} as an exact rational number.
        """

        l_use=(0 if l is None else l)
        if(n+s+l_use<0) : return Fraction(0, 1)
        epower=s+2*n+l_use
        assert(epower<len(self.__alpha))
        assert(n<len(self.__alpha[epower]))
        if l is None : return self.__alpha[epower][n]
        else : return self.__beta[l][epower][n]

    def alpha(self, s, n) : return self.alpha_or_beta(s, n)

    def beta(self, l, s, n) : return self.alpha_or_beta(s, n, l)

    def __call__(self, m, s, epower) :
        """ The coefficient in front of e^epower/omega in p_{m,s}. 
        
        Arguments:
            - m : The first index of the p_{m,s} being evaluated. Can only be
                  (0 or +-2).
            - s : The second index of the p_{m,s} being evaluated.
            - epower : The power of eccentricity of the term being evaluated.
        Returns: The specified coefficient as an exact rational number. """

        if epower%2!=s%2 : return 0
        n=(epower-s)//2
        if(m==0) : 
            if abs(s)>epower : return 0
            return self.alpha(s, n)*(1 if s==0 else Fraction(s,2)**epower)
        else :
            if s==0 or n<-1: return 0
            if m==-2 : msign=-1
            elif m==2 : msign=1
            else : raise ValueError('Asking for an expansion coefficient for'
                                    ' p_{m=%d, s=%d, p=%d}. The value of m '
                                    'must be one of -2, 0, 2!'%
                                    (m, s, epower))
            result=self.beta(-2, s, n+1)/2
            s2=s**2
            if n>=0 :
                result+=(self.alpha(s, n) - self.beta(0, s, n) - 
                         Fraction(2, s2)*self.beta(-2, s, n))
                if n>=1 : 
                    result+=(self.beta(2, s, n-1)/2 + 
                             Fraction(4, s2)*self.beta(0, s, n-1))
                    if n>=2 : result-=Fraction(2, s2)*self.beta(2, s, n-2)
            for k in range(0, n+2) :
                betas=-self.beta(-2, s, n-k+1)/2
                if n-k>=0 : 
                    betas+=Fraction(2, s)*self.beta(-1, s, n-k)
                    if n-k>=1 : 
                        betas+=(self.beta(2, s, n-k-1)/2
                                -
                                Fraction(2, s)*self.beta(1, s, n-k-1))
                result+=msign*Fraction(factorial(2*k),
                                       s2**k*factorial(k)**2*(2*k-1))*betas
        return result*Fraction(s, 2)**epower

def parse_command_line() :
    """ Parse the command line. """

    parser=ArgumentParser(description='Computes and stores the coefficients '
                          'defining the eccentricity expansion of the '
                          'various tidal terms. The coefficients are '
                          'computed exactly with no numerical roundoff as '
                          'rational numbers, but converted and stored as '
                          '16 significant figures floating point numbers. '
                          'Uses multiple CPUs to carry out the computation '
                          'and stores enough information to re-use a '
                          'previous computation to lower order.')
    parser.add_argument('-o', '--output', type=str,
                        default='eccentricity_expansion_coefficients',
                        help='The filename to store the computed expansion '
                        'coefficients to. The file is overwritten if it '
                        'exists. Default: "%(default)s".')
    parser.add_argument('-a', '--ab-file', type=str,
                        default='alpha_beta_values', help='The filename '
                        'where to store exact alpha and beta values '
                        'calculated while calculating the expansion '
                        'coefficients. Default: "%(default)s"')
    parser.add_argument('-c', '--cpus', type=int, default=cpu_count(),
                         help='The number of processes to use for the '
                         'calculation. Default: %(default)d.')
    parser.add_argument('-p', '--max-power', type=int, default=50,
                        help='Expansion is calculated up to this power of '
                        'the eccentricity. Default: %(default)d.')
    return parser.parse_args()

if __name__=='__main__' :
    '''
    usage example:

    >>> python tabulate_eccentricity_expansion_coefficients.py -c 16 > eccentricity_expansion_coef.txt
    '''
    options=parse_command_line()
    coef=ExpansionCoefficients(options.max_power, options.ab_file,
                               options.cpus)
    print(options.max_power)
    for epower in range(options.max_power+1) :
        for m in [-2, 0, 2] :
            for s in range(-epower+m, epower+m+1, 2) :
                if m==0 or s!=0 :
                    print("%25.16e"%(coef(m,s,epower)), end='')
        print()

