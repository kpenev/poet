import subprocess
import numpy
import psutil
import time

class Structure :
    """An empty class used only to hold user defined attributes."""

    def __init__(self, **initial_attributes) :
        """Create a class with (optionally) initial attributes."""

        for attribute_name, attribute_value in initial_attributes.items() :
            setattr(self, attribute_name, attribute_value)

    def format(self, prefix='') :
        """Generate a tree-like representation of self."""

        result=''
        for attr_name in dir(self) :
            if attr_name[0]!='_' :
                attribute=getattr(self, attr_name)
                if isinstance(attribute, Structure) :
                    result+=(prefix
                             +
                             '|-'
                             +
                             attr_name
                             +
                             '\n'
                             +
                             attribute.format(prefix + '| '))
                else : result+=(prefix
                                +
                                '|-'
                                +
                                attr_name
                                +
                                ': '
                                +
                                str(attribute)
                                +
                                '\n')
        return result

class PoetInterp :
    """Interface with POET to perform stellar evolution interpolation."""

    def start_poet(self) :
        """(Re-)start the poet process used for interpolation."""

        self.poet = subprocess.Popen(
            [
                self.poet_executable,
                '--input-columns', ','.join(self.input_columns),
                '--output-columns', ','.join(self.output_columns)
            ],
            stdin = subprocess.PIPE,
            bufsize = 0
        )
        self.poet_process_files = psutil.Process(self.poet.pid).open_files

    def __init__(self, poet_executable, output_columns) :
        """Start poet and wait for interpolation requests."""

        self.input_columns = ['M', 't0', 'tmax', 'maxdt', 'outf']
        self.input_line_format = ' '.join(
            ['%(' + c + ')s' for c in self.input_columns]
        ) + '\n'
        self.output_columns = output_columns
        self.poet_executable = poet_executable
        self.start_poet()

        self.output_fname = 'poet.evol'

    def __call__(self, star_mass, start_age, end_age, max_age_step) :
        """Return the interpolated stellar evolution per the arguments."""

        assert(not os.path.exists(self.output_fname))
        end_age = min(end_age, star_lifetime(star_mass))
        input_line = (self.input_line_format
                      %
                      dict(M = repr(star_mass),
                           t0 = repr(start_age),
                           tmax = repr(end_age),
                           maxdt = repr(max_age_step),
                           outf = self.output_fname)).encode('ascii')

        print(input_line)
        self.poet.stdin.write(input_line)
        while(not os.path.exists(self.output_fname)) : time.sleep(0.01)

        last_age = numpy.nan
        while not (end_age - last_age < 0.5 * max_age_step) :
            if self.poet.poll() is not None :
                print('Restarting poet!')
                self.start_poet()
                break
            try :
                result = numpy.genfromtxt(self.output_fname,
                                          names = True)
                last_age = result['t'][-1]
            except : 
                last_age = numpy.nan
                if self.poet.poll() is not None :
                    print('Restarting poet!')
                    self.start_poet()
                    self.poet.stdin.write(input_line)

        os.remove(self.output_fname)
        return result


