#!/usr/bin/env python
"""
To allow command-line updating of MESA inlist files

DLK 2010-12-17
"""

import getopt,sys,os,logging,shutil,datetime,re,subprocess,math,tempfile


######################################################################
def main():

    filelist=[]
    params=[]
    commands=[]
    values=[]
    sections=[]

    i=1
    section='any'
    while (i<len(sys.argv)):
        arg=sys.argv[i]
        isarg=0
        if (arg in ("-h", "--help")):
            usage()
            sys.exit(1)
        elif (arg in ("-m", "-a")):
            # modify or add
            params.append(sys.argv[i+1])
            if (arg in ("-m")):
                commands.append('modify')
            elif (arg in ("-a")):
                commands.append('add')            
            i+=1
            isarg=1
            if (i+1 < len(sys.argv)):
                values.append(sys.argv[i+1])
                i+=1
            else:
                values.append(None)
            sections.append(section)
        elif (arg in ("-d","-u","-c")):
            # delete, uncomment, comment
            params.append(sys.argv[i+1])
            if (arg in ("-d")):
                commands.append('delete')
            if (arg in ("-c")):
                commands.append('comment')
            if (arg in ("-u")):
                commands.append('uncomment')
            values.append(None)
            sections.append(section)            
            i+=1
            isarg=1
        elif (arg in ('-s')):
            # section
            section=sys.argv[i+1]
            i+=1
            isarg=1
        elif (arg[0] == '-'):
            print 'Argument: \"%s\" not recognized' % arg
            usage()
            sys.exit(1)
        if (not isarg):
            filelist.append(arg)
        i+=1


    if (len(filelist) == 0):
        print 'Must supply filename'
        usage()
        sys.exit(1)

    for filename in filelist:
        try:
            param_file=open(filename,'r')
        except IOError, err:
            print 'Unable to open inlist file: %s\n%s' % (filename,err)
            sys.exit(1)
        lines=param_file.readlines()
        param_file.close()

        for param,command,value,section in zip(params,commands,values,sections):
            lines=update_param(lines, param, value, command, section)
            if (not lines):
                print 'Unable to parse results of update_param()'
                sys.exit(1)

        fd,tmpfile=tempfile.mkstemp(suffix='.param', prefix='addtoparam', dir=os.path.abspath(os.getcwd()), text=True)
        os.close(fd)
        try:
            fout=open(tmpfile,'w')
        except IOError, err:
            print 'Unable to open temporary file: %s\n%s'% (tmpfile,err)
            sys.exit(1)
        fout.write('\n'.join(lines))
        fout.write('\n')
        fout.close()
        os.rename(tmpfile,filename)

        print "File %s updated" % filename
    sys.exit(0)

    



######################################################################
def update_param(filecontents, param=None, value=None, command=None, section=None):
    (xdir,xname)=os.path.split(sys.argv[0])
    now=datetime.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%S")
    blamestring='Modified by %s, %s' % (xname, now)

    lines=filecontents
    currentsection=None
    for i in xrange(len(lines)):
        line=lines[i]
        line=line.rstrip()
        lines[i]=line
        if (line.startswith('&')):
            currentsection=line[1:]
            if (command == 'add' and section == currentsection):
                newline='\n      %s = %s ! %s' % (param,value,blamestring)
                lines[i]+=newline
        if (section == currentsection or section == 'any'):
            if (param in line and command != 'add'):            
                if (command == 'modify'):
                    newline='      %s = %s ! %s' % (param,value,blamestring)
                elif (command == 'comment'):
                    newline=re.sub('^(\s+)','\\1!',line)
                    newline=newline[:-1] + (' ! %s' % blamestring)
                elif (command == 'uncomment'):
                    newline=line
                    #while ("!" in newline):
                    newline=re.sub('^(\s+)!','\\1',newline)
                    newline=newline[:-1] + (' ! %s' % blamestring)            
                elif (command == 'delete'):
                    newline='      ! %s' % blamestring
                else:
                    print "Unknown command: %s" % command
                    return None
                lines[i]=newline
    return lines

######################################################################
def usage():
    (xdir,xname)=os.path.split(sys.argv[0])
    print "Usage:"
    print "\t%s\t [-h] [-u <param>] [-m <param> <value>] [-a <param> <value>] [-d <param>] [-c <param>] [-s <section>] <filename> [<filename2> ...]" % xname
    print "\tUpdates MESA inlist file(s) specified"
    print "\tPossible commands:"
    print "\t\t-m: modify (supply parameter and value)"
    print "\t\t-a: add (supply parameter and value)"
    print "\t\t-c: comment (supply parameter)"
    print "\t\t-u: uncomment (supply parameter)"
    print "\t\t-d: delete (supply parameter)"
    print "\t\t-s: specify section"
     


######################################################################

if __name__=="__main__":
    main()
