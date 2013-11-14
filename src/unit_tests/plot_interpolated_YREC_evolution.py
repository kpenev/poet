#!/usr/bin/python -u

from matplotlib import pyplot
from matplotlib.backends.backend_pdf import PdfPages 
from scipy import array, arange

YREC_files={'0.5': '../YREC/m0p5_alpha1p9_Y0p28_GN93.track',
            '0.6': '../YREC/m0p6_alpha1p9_Y0p28_GN93.track',
            '0.7': '../YREC/m0p7_alpha1p9_Y0p28_GN93.track',
            '0.8': '../YREC/m0p8_alpha1p9_Y0p28_GN93.track',
            '0.9': '../YREC/m0p9_alpha1p9_Y0p28_GN93.track',
            '1.0': '../YREC/m1p0_alpha1p9_Y0p28_GN93.track',
            '1.05': '../YREC/m1p05_alpha1p9_Y0p28_GN93.track',
            '1.1': '../YREC/m1p1_alpha1p9_Y0p28_GN93.track',
            '1.15': '../YREC/m1p15_alpha1p9_Y0p28_GN93.track',
            '1.2': '../YREC/m1p2_alpha1p9_Y0p28_GN93.track'}

plot_colors=['#FF0000', '#03B402', '#0005F7', '#FF6D00', '#FF00FF',
             '#8F7000', '#0000FF']

def evol_fname(mass) :
    """ Returns the filename with the interpolated track corresponding to the
    given mass. """

    return "test_YREC_interp/test_YREC_interp_M%.3f.evol"%mass

def read_YREC_file(fname, mass) :
    """ Reads a YREC track. """

    Msun=1.98892e30
    Rsun=6.955e8
    Inorm=Msun*Rsun*Rsun*1e7

    f=open(fname)
    result={'t':[], 'R':[], 'Iconv':[], 'Irad':[], 'Mrad':[], 'Rcore':[]}
    for l in f :
        if l[0]=='#' : continue
        entries=map(eval, l.split())
        result['t'].append(entries[2])
        result['R'].append(10.0**entries[7])
        result['Mrad'].append(entries[10 if mass<1.01 else 12]*mass)
        result['Rcore'].append(entries[11 if mass<1.01 else 13]*10.0**entries[7])
        result['Irad'].append(entries[15 if mass<=1.01 else 17]/Inorm)
        result['Iconv'].append(entries[16 if mass<=1.01 else 18]/Inorm)
    return result

def read_evol_file(fname) :
    """ Reads an output interpolated evolution file. """

    f=open(fname)
    result={'t':[], 'R':[[], [], []], 'L':[[], [], []], 'Iconv':[[], [], []],
            'Irad':[[], [], []], 'Mrad':[[], [], []], 'Rcore':[[], [], []]}
    for l in f :
        if l[0]=='#' : continue
        entries=map(eval, l.split())
        result['t'].append(entries[0])
        for i in range(3) :
            result['R'][i].append(entries[1+i])
            result['L'][i].append(entries[4+i])
            result['Iconv'][i].append(entries[7+i])
            result['Mrad'][i].append(entries[13+i])
            result['Rcore'][i].append(entries[16+i])
            result['Irad'][i].append(entries[19+i])
    return result

def plot_quantity(quantity, low_mass, high_mass, YREC, evol, pdf, logx=True,
              logy=True) :
    """ Plots the tracks of the given quantity (a strintg key) for masses
    in the given range as a separate page in the given PDF. """

    if logx and logy : plot_cmd=pyplot.loglog
    elif logx==True : plot_cmd=pyplot.semilogx
    elif logy==True : plot_cmd=pyplot.semilogy
    else : plot_cmd=pyplot.plot
    plot_cmd(YREC[str(low_mass)]['t'], YREC[str(low_mass)][quantity], '-r')
    plot_cmd(YREC[str(high_mass)]['t'], YREC[str(high_mass)][quantity], '-b')
    color_ind=0
    line_style='--'
    for m in arange(low_mass, high_mass, 0.001) :
        plot_cmd(evol[str(m)]['t'], evol[str(m)][quantity][0],
                 color=plot_colors[color_ind], linestyle=line_style)
        first_t=0
        while evol[str(m)]['t'][first_t]<0.005 : first_t+=1
        min_ind=evol[str(m)][quantity][0][first_t:].index(min(evol[str(m)][quantity][0][first_t:]))
        print q, m, evol[str(m)]['t'][min_ind], evol[str(m)][quantity][0][min_ind]
        color_ind+=1
        if color_ind==len(plot_colors) :
           color_ind=0
           if line_style=='--' : line_style=':'
           else : line_style='--'
    pyplot.suptitle("M="+str(low_mass)+"$M_\odot$ to "+str(high_mass)+
                    "$M_\odot$ "+quantity)
    pdf.savefig()
    pyplot.close()

if __name__=='__main__' :
    YREC=dict()
    for mass, fname in YREC_files.iteritems() :
        YREC[mass]=read_YREC_file(fname, eval(mass))
    evol=dict()
    for mass in arange(0.5, 1.205, 0.001) :
        evol[str(mass)]=read_evol_file(evol_fname(mass))
    pdf = PdfPages('test_YREC_interp.pdf')
    YREC_masses=sorted(map(eval, YREC_files.keys()))
    for q in ['R', 'Iconv', 'Irad', 'Mrad', 'Rcore'] :
        for i in range(len(YREC_masses)-1) :
            plot_quantity(q, YREC_masses[i], YREC_masses[i+1], YREC,
                          evol, pdf, True, q in ['R', 'Iconv'])
    pdf.close()
    exit(1)
