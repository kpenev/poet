#!/usr/bin/python -u

from optparse import OptionParser

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
    dest=result['age']['accepted']
    assert entries[0].startswith('Age:')
    for e in entries :
        if e[-1]=='|' :
            dest=result['age']['rejected']
            e=e[:-1]
        else : dest.append(eval(e))
    line=f.readline()
    assert(line[0]=='_' and line.strip('_')=='\n')
    while(line.strip() not in ['Accepted', 'Rejected']) :
        if line.startswith('Condition[') :
            result['conditions'].append(dict(accepted=[], rejected=[]))
            dest=result['conditions'][-1]['accepted']
            assert(eval(line[line.find('[')+1:line.find(']')])==
                   len(result['conditions']-1))
            line=line[line.find(']:')+2:]
            if line[0]=='|' : 
                dest=result['conditions'][-1]['rejected']
                line=line[1:]
            if line[0]=='z' :
                result['zero skip']=len(result['conditions'][-1]['accepted'])
            if line[1]=='e' :
                result['extremum skip']=len(
                    result['conditions'][-1]['accepted'])
            if line[2]=='>' :
                line=line[3:]
            line=line.strip()
        line=f.readline()

if __name__=='__main__' :
    options, log_file=parse_command_line()
    f=open(log_file, 'r')
    line=f.readline()
    while line :
        if line[0]=='@' and line.strip('@')=='\n' :
            step_info=dict()
        elif line.strip()=='Stored stop condition information:' :
            step_info['stored stop conditions']=\
                    read_stored_stop_conditions(f)
        line=f.readline()
