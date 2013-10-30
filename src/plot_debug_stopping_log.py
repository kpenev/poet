#!/usr/bin/python -u

from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import pyplot as plt
from optparse import OptionParser
from numpy import nan, inf, isfinite, array, arange

def parse_command_line() :
    """ Returns the command line options as members of a structure. """

    parser=OptionParser(usage="%prog <log file>", description="Generates a "
                        "PDF file with each page showing all the information"
                        " about an attempted step.")
    parser.add_option('-o', '--output', default='debug_stopping.pdf',
                      help='Filename to use for the plots. '
                      'Default: %default')
    options, args=parser.parse_args()
    return options, args[0]

def read_stored_stop_conditions(f) :
    """ Reads a single  'Stored stop condition information' section from the
    given file. """

    line=f.readline()
    entries=line.split()
    result=dict(age=dict(accepted=[], rejected=[]), conditions=[],
                derivatives=[], zero_skip=[], extremum_skip=[])
    assert entries[0].startswith('Age:')
    if entries[0][-1]=='|' : dest=result['age']['rejected']
    else : dest=result['age']['accepted']
    for e in entries[1:] :
        dest.append(eval(e.strip('|')))
        if e[-1]=='|' :
            dest=result['age']['rejected']
            e=e[:-1]
    line=f.readline()
    assert(line[0]=='_' and line.strip('_')=='\n')
    line=f.readline()
    while(line.startswith('   Condition[') or
          line.startswith('  Derivative[')) :
        if line.startswith('   Condition[') :
            what='conditions'
            if(len(result['zero_skip'])==len(result[what])-1) :
                result['zero_skip'].append(len(result[what][-1]['accepted']))
            if(len(result['extremum_skip'])==len(result[what])-1) :
                result['extremum_skip'].append(
                    len(result[what][-1]['accepted']))
            result['conditions'].append(dict(accepted=[], rejected=[]))
        elif line.startswith('  Derivative[') :
            result['derivatives'].append(dict(accepted=[], rejected=[]))
            what='derivatives'
        dest=result[what][-1]['accepted']
        assert(eval(line[line.find('[')+1:line.find(']')])==
               len(result[what])-1)
        line=line[line.find(']:')+2:]
        while line :
            if line[0]=='|' : 
                dest=result[what][-1]['rejected']
                line=line[1:]
            if line[0]=='z' :
                result['zero_skip'].append(len(result[what][-1]['accepted']))
                assert len(result['zero_skip'])==len(result[what])
            if line[1]=='e' :
                result['extremum_skip'].append(
                    len(result[what][-1]['accepted']))
                assert len(result['extremum_skip'])==len(result[what])
            if line[2]=='>' :
                line=line[3:]
            line=line.strip()
            cut1, cut2, cut3=line.find(' '), line.find('|'), line.find('z')
            if cut1==-1 : cut1=max(cut2, cut3)
            if cut2==-1 : cut2=max(cut1, cut3)
            if cut3==-1 : cut3=max(cut1, cut2)
            cut=min(cut1, cut2, cut3)
            if cut==-1 : cut=len(line)
            dest.append(eval(line[:cut]))
            line=line[cut:]
        line=f.readline()
    what='conditions'
    if(len(result['zero_skip'])==len(result[what])-1) :
        result['zero_skip'].append(len(result[what][-1]['accepted']))
    if(len(result['extremum_skip'])==len(result[what])-1) :
        result['extremum_skip'].append(
            len(result[what][-1]['accepted']))

    return line, result

def read_condition_stop(line, f) :
    """ Parses the information for where to stop because of the current
    condition. """

    assert(line.startswith('Condition['))
    result=dict(name=line[line.index('(')+1:line.index(')')],
                value=eval(line[line.index('=')+1:]))
    line=f.readline()
    assert(line.startswith('crossing: age='))
    result['crossing']=dict(age=eval(line[line.index('=')+1:]))
    line=f.readline()
    assert(line.startswith('extremum: age='))
    result['extremum']=dict(
        age=eval(line[line.index('=')+1:line.index(',')]))
    line=line[line.index(',')+1:]
    result['extremum']['value']=eval(line[line.index('=')+1:])
    return result

def read_final_stop(f) :
    """ Reads the final stop information for the step. """

    result=dict()
    line=f.readline()
    assert(line.startswith('Stop at t='))
    result['age']=eval(line[line.index('=')+1:line.index(',')])
    line=line[line.index(',')+1:]
    result['type'], skip, result['name'], line=line.split(None, 3)
    result['name']=result['name'][:-1]
    assert(line.startswith('precision='))
    result['precision']=eval(line[line.index('=')+1:])
    return result

def read_interval(f) :
    """ Reads an interval table from the given file. """

    line=f.readline()
    if(line.strip()=='') : line=f.readline()
    entries=line.replace('|', ' ').split()
    assert(entries[0]=='Age:')
    result=dict(ages=map(eval, entries[1:]), conditions=[], derivatives=[])
    line=f.readline()
    while(line.strip()) :
        if line.strip().startswith('Condition[') :
            dest=result['conditions']
        else :
            assert line.strip().startswith('Derivative[')
            dest=result['derivatives']
        index=eval(line[line.index('[')+1:line.index(']')])
        assert(index==len(dest))
        dest.append(
            map(eval, line[line.index(':')+1:].replace('|', ' ').split()))
        assert(len(dest[-1])==len(result['ages']))
        line=f.readline()
    return result

def parse_comma_separated_list(csl) :
    """ Parses a list of values of the form '<value1>, <value2>[, ...]. """

    print 'csl=', csl
    return map(lambda s: eval(s.strip()), csl.split(','))

def parse_zerocrossing(line) :
    """ Parses a line containing the information about how a zerocrossing is
    found. """

    line=line.strip()
    if(line.startswith('Linear')) : order=1;
    elif(line.startswith('Quadratic')) : order=2;
    elif(line.startswith('Cubic')) : order=3;
    else : assert(False)
    result=dict(order=order)
    between=parse_comma_separated_list(
        line[line.index('(')+1:line.index(')')])
    if(len(between)==3) :
        assert(order==3)
        has_deriv=True
    else :
        assert(len(between)==2)
        has_deriv=False
    result=dict(ages=[between[0]], values=[between[1]],
                derivs=([between[2]] if has_deriv else []))
    line=line[line.find(')')+1:]
    between=parse_comma_separated_list(
        line[line.index('(')+1:line.index(')')])
    result['ages'].append(between[0])
    result['values'].append(between[1])
    if has_deriv :
        assert(len(between)==3)
        result['derivs'].append(between[2])
    else : assert(len(between)==2)
    line=line[line.find(')')+1:].strip()
    if(order==2 or (order==3 and not has_deriv)) :
        between=parse_comma_separated_list(
            line[line.index('(')+1:line.index(')')])
        assert(len(between)==2)
        result['ages'].append(between[0])
        result['values'].append(between[1])
        line=line[line.find(')')+1:].strip()
    if(order==3 and not has_deriv) :
        between=parse_comma_separated_list(
            line[line.index('(')+1:line.index(')')])
        assert(len(between)==2)
        result['ages'].append(between[0])
        result['values'].append(between[1])
        line=line[line.find(')')+1:].strip()
    if(order==1) :
        result['coef']=None
        assert(line.startswith('='))
        result['selected']=eval(line[1:])
        result['solutions']=[result['selected']]
        result['fallback']=None
    else :
        if(order==2 or (order==3 and not has_deriv)) :
            assert(line.startswith('in range ('))
            result['range']=tuple(parse_comma_separated_list(
                line[line.index('(')+1:line.index(')')]))
            assert(len(result['range'])==2)
            line=line[line.find(')')+1:].strip()
        assert(line.startswith(', coef=('))
        result['coef']=parse_comma_separated_list(
                line[line.index('(')+1:line.index(')')])
        line=line[line.find(')')+1:].strip()
        assert(line.startswith(', solutions=('))
        result['solutions']=parse_comma_separated_list(
                line[line.index('(')+1:line.index(')')])
        line=line[line.find(')')+1:].strip()
        if(line.startswith(', selected')) :
            result['fallback']=None
            result['selected']=eval(line[line.index(':')+1:])
        else :
            assert(line.startswith(', fallback to linear:'))
            result['selected']=None
            if(order==3 and has_deriv) :
                result['fallback']=eval(line[line.index(':')+1:])
            else :
                line=line[line.index(':')+1:]
                result['fallback']=parse_zerocrossing(line)

    if(order==1 or (order==3 and has_deriv)) :
        result['range']=tuple(result['ages'])
    return result

def parse_extremum(line) :
    """ Parses a line containing the information about how a zerocrossing is
    found. """

    line=line.strip()
    if(line.startswith('Quadratic')) : order=2;
    elif(line.startswith('Cubic')) : order=3;
    else : assert(False)
    result=dict(order=order)
    between=parse_comma_separated_list(
        line[line.index('(')+1:line.index(')')])
    if(len(between)==3) :
        assert(order==3)
        has_deriv=True
    else :
        assert(len(between)==2)
        has_deriv=False
    result=dict(ages=[between[0]], values=[between[1]],
                derivs=([between[2]] if has_deriv else []))
    line=line[line.find(')')+1:]
    between=parse_comma_separated_list(
        line[line.index('(')+1:line.index(')')])
    result['ages'].append(between[0])
    result['values'].append(between[1])
    if has_deriv :
        assert(len(between)==3)
        result['derivs'].append(between[2])
    else : assert(len(between)==2)
    line=line[line.find(')')+1:].strip()
    if(not has_deriv) :
        between=parse_comma_separated_list(
            line[line.index('(')+1:line.index(')')])
        assert(len(between)==2)
        result['ages'].append(between[0])
        result['values'].append(between[1])
        line=line[line.find(')')+1:].strip()
    if(order==3 and not has_deriv) :
        between=parse_comma_separated_list(
            line[line.index('(')+1:line.index(')')])
        assert(len(between)==2)
        result['ages'].append(between[0])
        result['values'].append(between[1])
        line=line[line.find(')')+1:].strip()

    if(order==3 and not has_deriv) :
        assert(line.startswith('in range ('))
        result['range']=tuple(parse_comma_separated_list(
            line[line.index('(')+1:line.index(')')]))
        assert(len(result['range'])==2)
        line=line[line.find(')')+1:].strip()
    else :
        result['range']=(result['ages'][0], result['ages'][-1])
    assert(line.startswith(', coef=('))
    result['coef']=parse_comma_separated_list(
            line[line.index('(')+1:line.index(')')])
    line=line[line.find(')')+1:].strip()
    if order==2 :
        assert(line.startswith(': ('))
        result['selected']=tuple(parse_comma_separated_list(
            line[line.index('(')+1:line.index(')')]))
        result['solutions']=[result['selected']]
        result['fallback']=None
    else :
        assert(line.startswith(', solutions: ('))
        result['solutions']=[tuple(parse_comma_separated_list(
            line[line.index('(')+1:line.index(')')]))]
        line=line[line.find(')')+1:].strip()
        result['solutions'].append(tuple(parse_comma_separated_list(
            line[line.index('(')+1:line.index(')')])))
        line=line[line.find(')')+1:].strip()
        if(line.startswith(', selected')) :
            result['fallback']=None
            line=line[line.index(':')+1:]
            result['selected']=tuple(parse_comma_separated_list(
                line[line.index('(')+1:line.index(')')]))
        else :
            assert(line.startswith(', fallback to '))
            result['selected']=None
            if(has_deriv) :
                result['fallback']=tuple(parse_comma_separated_list(
                    line[line.index('(')+1:line.index(')')]))
            else :
                line=line[line.index(':')+1:]
                result['fallback']=parse_extremum(line)
    return result

def plot_step(step_info) :
    """ Plots all the information contained in the given step_info on a page
    in the given PdfPages file. """

    def plot_polynomial(coefficients, axis, style, resolution=200) :
        """ Plots the polynomial defined by the given coefficients on the
        given axis. """

        xlim, ylim=axis.get_xlim(), axis.get_ylim()
        x=arange(xlim[0], xlim[1], (xlim[1]-xlim[0])/resolution)
        xpow=x**0
        y=None
        for c in reversed(coefficients) :
            if y is None : y=c*xpow
            else : y+=c*xpow
            xpow*=x
        axis.plot(x, y, style)
        axis.set_xlim(xlim)
        axis.set_ylim(ylim)

    def make_full_plot(ax_full, cond_ind) :
        """ Makes the full plot. """

        if(not step_info['stored stop conditions']['age']['accepted']) :
            switch=step_info['stored stop conditions']['age']['rejected'][0]
            ax_full.axvspan(
                step_info['stored stop conditions']['age']['rejected'][0],
                step_info['stored stop conditions']['age']['rejected'][-1],
                facecolor='#FFAAAA', linewidth=0, edgecolor='#FFAAAA')
        else :
            switch=0.5*(
                step_info['stored stop conditions']['age']['accepted'][-1]+
                step_info['stored stop conditions']['age']['rejected'][0])
            ax_full.axvspan(
                step_info['stored stop conditions']['age']['accepted'][0],
                switch, facecolor='#AAFFAA', linewidth=0,
                edgecolor='#AAFFAA')
            ax_full.axvspan(
                switch,
                step_info['stored stop conditions']['age']['rejected'][-1],
                facecolor='#FFAAAA', linewidth=0, edgecolor='#FFAAAA')
        zero_first=(
            step_info['stored stop conditions']['zero_skip'][cond_ind]<=
            step_info['stored stop conditions']['extremum_skip'][cond_ind])
        all_skip=min(
            step_info['stored stop conditions']['zero_skip'][cond_ind],
            step_info['stored stop conditions']['extremum_skip'][cond_ind])
        one_skip=max(
            step_info['stored stop conditions']['zero_skip'][cond_ind],
            step_info['stored stop conditions']['extremum_skip'][cond_ind])
        ax_full.plot(step_info['stored stop conditions']['age']['accepted']\
                     [:all_skip],
                     step_info['stored stop conditions']['conditions']\
                     [cond_ind]['accepted'][:all_skip], 'o', 
                     markeredgecolor='black', markerfacecolor='black')
        ax_full.plot(step_info['stored stop conditions']['age']['accepted']\
                     [all_skip:one_skip],
                     step_info['stored stop conditions']['conditions']\
                     [cond_ind]['accepted'][all_skip:one_skip], 'o', 
                     markeredgecolor=('black' if zero_first else 'blue'),
                     markerfacecolor=('none' if zero_first else 'blue'))
        ax_full.plot(step_info['stored stop conditions']['age']['accepted']\
                     [one_skip:]+
                     step_info['stored stop conditions']['age']['rejected'],
                     step_info['stored stop conditions']['conditions']\
                     [cond_ind]['accepted'][one_skip:]+
                     step_info['stored stop conditions']['conditions']\
                     [cond_ind]['rejected'], 'o',
                     markeredgecolor='blue', markerfacecolor='none')
        if(step_info['condition_details'][cond_ind]) :
            details=step_info['condition_details'][cond_ind]
            if 'zerocrossing' in details :
                coef=details['zerocrossing']['coef']
                if coef is not None :
                    plot_polynomial(coef, ax_full, '-r')
            if 'extremum' in details :
                coef=details['extremum']['coef']
                if coef is not None :
                    plot_polynomial(coef, ax_full, '-b')
            if 'crossing_interval' in details :
                ax_full.plot(
                    details['crossing_interval']['ages'],
                    details['crossing_interval']['conditions'][cond_ind],
                    'xr')
            if 'extremum_interval' in details :
                ax_full.plot(
                    details['extremum_interval']['ages'],
                    details['extremum_interval']['conditions'][cond_ind],
                    'xg')
        condition_stop=step_info['condition_stops'][cond_ind]
        if(isfinite(condition_stop['crossing']['age'])) :
            ax_full.axvline(condition_stop['crossing']['age'], color='black')
        if(isfinite(condition_stop['extremum']['age'])) :
            ax_full.axvline(condition_stop['extremum']['age'], color='red')
            ax_full.axhline(condition_stop['extremum']['value'], color='red')
        return switch

    def make_zoom_plot(ax_zoom, cond_ind, reject_boundary) :
        """ Makes a plot that zooms only on the search interval. """

        details=step_info['condition_details'][cond_ind]
        if 'zerocrossing' in details :
            xmin=details['zerocrossing']['ages'][0]
            xmax=details['zerocrossing']['ages'][-1]
            ymin=min(details['zerocrossing']['values'])
            ymax=max(details['zerocrossing']['values'])
        else : xmin, xmax, ymin, ymax=inf, -inf, inf, -inf
        if 'extremum' in details :
            xmin=min(details['extremum']['ages'][0], xmin)
            xmax=max(details['extremum']['ages'][-1], xmax)
            ymin=min(min(details['extremum']['values']), ymin)
            ymax=max(max(details['extremum']['values']), ymax)
        if reject_boundary>xmin and reject_boundary<xmax :
            ax_zoom.axvspan(xmin, reject_boundary, facecolor='#AAFFAA',
                            linewidth=0, edgecolor='#AAFFAA')
            ax_zoom.axvspan(reject_boundary, xmax, facecolor='#FFAAAA',
                            linewidth=0, edgecolor='#FFAAAA')
        elif reject_boundary<xmin :
            ax_zoom.axvspan(xmin, xmax, facecolor='#FFAAAA',
                            linewidth=0, edgecolor='#FFAAAA')
        else : 
            ax_zoom.axvspan(xmin, xmax, facecolor='#AAFFAA',
                            linewidth=0, edgecolor='#AAFFAA')
        ax_zoom.set_xlim(xmin-0.05*(xmax-xmin), xmax+0.05*(xmax-xmin))
        ax_zoom.set_ylim(ymin-0.05*(ymax-ymin), ymax+0.05*(ymax-ymin))
        if 'extremum' in details :
            plot_polynomial(details['extremum']['coef'], ax_zoom, '-b')
        if 'zerocrossing' in details :
            if details['zerocrossing']['coef'] :
                plot_polynomial(details['zerocrossing']['coef'], ax_zoom,
                                '-r')
            else :
                assert len(details['zerocrossing']['ages'])==2
                ax_zoom.plot(details['zerocrossing']['ages'], 
                             details['zerocrossing']['values'], '-r')
        if 'extremum' in details :
            if details['extremum']['selected'] is not None :
                extremum=details['extremum']['selected']
            else :
                if type(details['extremum']['fallback']) is tuple :
                    extremum=details['extremum']['fallback']
                else : extremum=details['extremum']['fallback']['selected']
            ax_zoom.plot([extremum[0]], [extremum[1]], 'o',
                         markeredgecolor='none', markerfacecolor='blue')
            ax_zoom.plot(details['extremum']['ages'], 
                         details['extremum']['values'], '+b')
        if 'zerocrossing' in details :
            if details['zerocrossing']['selected'] is not None :
                crossing=details['zerocrossing']['selected']
            elif type(details['zerocrossing']['fallback']) is dict :
                crossing=details['zerocrossing']['fallback']['selected']
            else : crossing=details['zerocrossing']['fallback']
            ax_zoom.plot([crossing], [0], 'o',  markeredgecolor='none',
                         markerfacecolor='red')
            ax_zoom.plot(details['zerocrossing']['ages'], 
                         details['zerocrossing']['values'], 'xr')

    plt.clf()
    xpad=0.05
    ypad=0.05
    full_plot_fraction=0.66

    num_cond=len(step_info['stored stop conditions']['conditions'])
    width=full_plot_fraction*(1.0-3.0*xpad)
    for cond_ind in range(num_cond) :
        height=(1.0-ypad*(num_cond+1))/num_cond
        ax_full=plt.axes([xpad, ypad*(cond_ind+1)+height*cond_ind, width,
                          height])
        reject_boundary=make_full_plot(ax_full, cond_ind)
        if(step_info['condition_details'][cond_ind]) :
            ax_zoom=plt.axes([xpad*2+width,
                              ypad*(cond_ind+1)+height*cond_ind,
                              1.0-xpad*3-width, height])
            make_zoom_plot(ax_zoom, cond_ind, reject_boundary)

if __name__=='__main__' :
    options, log_file=parse_command_line()
    f=open(log_file, 'r')
    line=f.readline()
    step_info=dict()
    pdf=PdfPages(options.output)
    while line :
        if line[0]=='@' and line.strip('@')=='\n' :
            if(step_info and
               len(step_info['stored stop conditions']['conditions'])) :
                print 'plotting'
                plot_step(step_info)
            else : plt.plot([], [])
            print 'saving'
            pdf.savefig()
            print 'done'
            step_info=dict(condition_stops=[], condition_details=[])
        elif line.strip()=='Stored stop condition information:' :
            line, step_info['stored stop conditions']=\
                    read_stored_stop_conditions(f)
            condition_details=dict()
            continue
        elif line.strip()==('Crossing search interval without derivative '
                            'information:') :
            condition_details['crossing_interval']=read_interval(f)
        elif line.strip()==('Extremum search interval without derivative '
                            'information:') :
            condition_details['extremum_interval']=read_interval(f)
        elif line.startswith('Condition[') :
            step_info['condition_stops'].append(read_condition_stop(line, f))
            step_info['condition_details'].append(dict(condition_details))
            condition_details=dict()
        elif line.strip()=='Final stop information:' :
            step_info['final_stop']=read_final_stop(f)
            line=f.readline()
            if(line.strip()=='Accepted') : step_info['accepted']=True
            else :
                assert(line.strip()=='Rejected')
                step_info['accepted']=False
        elif (line.strip().startswith('Linear zerocrossing between') or 
              line.strip().startswith('Quadratic zerocrossing between') or 
              line.strip().startswith('Cubic zerocrossing between')) :
            condition_details['zerocrossing']=parse_zerocrossing(line)
        elif (line.strip().startswith('Quadratic extremum between') or 
              line.strip().startswith('Cubic extrema between')) :
            condition_details['extremum']=parse_extremum(line)
        line=f.readline()
    pdf.close()