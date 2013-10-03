#!/usr/bin/python -u

from optparse import OptionParser
from numpy import nan, inf

def parse_command_line() :
    """ Returns the command line options as members of a structure. """

    parser=OptionParser(usage="%prog <log file>", description="Generates a "
                        "PDF file with each page showing all the information"
                        " about an attempted step.")
    options, args=parser.parse_args()
    return options, args[0]

def read_stored_stop_conditions(f) :
    """ Reads a single  'Stored stop condition information' section from the
    given file. """

    line=f.readline()
    entries=line.split()
    result=dict(age=dict(accepted=[], rejected=[]), conditions=[],
                derivatives=[])
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
            result['conditions'].append(dict(accepted=[], rejected=[]))
            what='conditions'
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
                result['zero skip']=len(result[what][-1]['accepted'])
            if line[1]=='e' :
                result['extremum skip']=len(
                    result[what][-1]['accepted'])
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
    if(line.strip()=='') line=f.readline()
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

    return map(lambda s: eval(s.strip()), cls.split(','))

def parse_zerocrossing(line) :
    """ Parses a line containing the information about how a zerocrossing is
    found. """

    line=line.strip()
    if(line.startswith('Linear')) order=1;
    elif(line.startswith('Quadratic')) order=2;
    elif(line.startswith('Cubic')) order=3;
    else : assert(False)
    result['order']=order
    between=parse_comma_separated_list(line[line.index('(')+1:line.index(')')])
    if(len(between)==3) :
        assert(order==3)
        has_deriv=True
    else :
        assert(len(between)==2)
        has_deriv=False
    result=dict(ages=[between[0]], values=[between[1]],
                derivs=([] if has_deriv else [between[2]]))
    line=line[line.find(')'):]
    between=parse_comma_separated_list(line[line.index('(')+1:line.index(')')])
    result['ages'].append(between[0])
    result['values'].appned(between[1])
    if has_deriv :
        assert(len(between)==3)
        result['derivs'].append(between[2])
    else : assert(len(between)==2)
    line=line[line.find(')'):].strip()
    if(order==2 or (order==3 and not has_deriv)) :
        between=parse_comma_separated_list(line[line.index('(')+1:line.index(')')])
        assert(len(between)==2)
        result['ages'].append(between[0])
        result['values'].appned(between[1])
        line=line[line.find(')'):].strip()
    if(order==3 and not has_deriv) :
        between=parse_comma_separated_list(line[line.index('(')+1:line.index(')')])
        assert(len(between)==2)
        result['ages'].append(between[0])
        result['values'].appned(between[1])
        line=line[line.find(')'):].strip()
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
            line=line[line.find(')'):].strip()
        assert(line.startswith(', coef=('))
        result['coef']=parse_comma_separated_list(
                line[line.index('(')+1:line.index(')')])
        line=line[line.find(')'):].strip()
        assert(line.startswith(', solutions=('))
        result['solutions']=parse_comma_separated_list(
                line[line.index('(')+1:line.index(')')])
        line=line[line.find(')'):].strip()
        if(line.startswith(', selected')) :
            result['fallback']=None
            result['selected']=eval(line[line.index(':'):])
        else :
            assert(line.startswith(', fallback to linear:'))
            if(order==3 and has_deriv) :
                result['fallback']=eval(line[line.index(':'):])
            else :
                line=line[line.index(':'):]
                result['fallback']=parse_zerocrossing(line)

    if(order==1 or (order==3 and has_deriv))
        result['range']=tuple(result['ages'])
    return result

if __name__=='__main__' :
    options, log_file=parse_command_line()
    f=open(log_file, 'r')
    line=f.readline()
    step_info=dict()
    while line :
        if line[0]=='@' and line.strip('@')=='\n' :
            print 'step_info=', step_info
            step_info=dict(condition_stops=[], condition_details=[])
        elif line.strip()=='Stored stop condition information:' :
            line, step_info['stored stop conditions']=\
                    read_stored_stop_conditions(f)
            condition_details=dict()
            continue
        elif line.strip()==('Crossing search interval without derivative '
                            'information:') :
            condition_details['crossing_interval']=read_interval(f))
        elif line.strip()==('Extremum search interval without derivative '
                            'information:') :
            condition_details['extremum_interval']=read_interval(f)))
        elif line.startswith('Condition[') :
            step_info['condition_stops'].append(read_condition_stop(line, f))
            step_info['condition_details'].append(dict(condition_details))
            condition_details=dict()
        elif line.strip()=='Final stop information:' :
            step_info['final_stop']=read_final_stop(f)
            line=f.readline()
            if(line.strip()=='Accepted') step_info['accepted']=True
            else :
                assert(line.strip()=='Rejected')
                step_info['accepted']=False
        elif line.strip().startswith('Linear zerocrossing between') :
            condition_details['linear_zerocrossing']=parse_zerocrossing(line)
        line=f.readline()
