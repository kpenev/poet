import sys, getopt
import numpy as np
import scipy.integrate as integrate
#import matplotlib.pyplot as plt
from multiprocessing import Pool
#import multiprocessing as mp
from functools import partial

def p_0s(u,s,e):
    return (1-e**2)**(3/2)*(np.exp(1j*s*(u-e*np.sin(u)))/(1-e*np.cos(u))**2).real
def t1(u,s,e):
    return ((1-e**2)**(3/2)*(np.exp(1j*s*(u-e*np.sin(u)))/(1-e*np.cos(u))**4)).real
def t2(u,s,e):
    return (1-e**2)**(3/2)*(np.exp(1j*s*(u-e*np.sin(u)))*np.cos(2*u)/(1-e*np.cos(u))**4).real
def t3(u,s,e):
    return -((1-e**2)**(3/2)*(np.exp(1j*s*(u-e*np.sin(u)))*np.sin(2*u)/(1-e*np.cos(u))**4)).imag
def t4(u,s,e):
    return -((1-e**2)**(3/2)*(np.exp(1j*s*(u-e*np.sin(u)))*np.sin(u)/(1-e*np.cos(u))**4)).imag

def integrize(kind,func,s,e,lB=50):
    y=p_0s
    if func == 1:
        y=t1
    if func == 2:
        y=t2
    if func == 3:
        y=t3
    if func == 4:
        y=t4
    
    if kind == 0:
        return integrate.quad(y,0,2*np.pi,args=(s,e),limit=lB,points=np.linspace(0,1,4*s))[0]
    if kind == 1:
        return integrate.romberg(y,0,2*np.pi,args=(s,e),divmax=40,tol=1e-10,rtol=1e-10)
    if kind == 2:
        return integrate.quadrature(y,0,2*np.pi,args=(s,e),maxiter=10000,tol=1e-10,rtol=1e-10)[0]
    if kind == 3:
        return integrate.fixed_quad(y,0,2*np.pi,args=(s,e),n=10000)[0]
    #return

def p_MS(e,m,s,kind,shouldPrint=0,objetFile=None,itera=0,toler=1e-9):
    m=int(m)
    s=int(s)
    #kind=whatKind
    limitBreak = 10*s
    if limitBreak < 50:
        limitBreak = 50
    
    #i=0
    #value = 0
    
    value = 0
    if e == 1: # Special case of we-know-what-this-should-be
        #print("hello", e)
        if m == 0:
            value=1
    
    else:
        p0s = integrize(kind,0,s,e,limitBreak)
        #print("also hello", e)
        if abs(m) == 2:
        
            # Define a couple of constants that show up across terms
            sC1 = 1-e**2
            sC2 = np.sqrt(sC1)
            
            term1 = -sC1*integrize(kind,1,s,e,limitBreak)
            term2 = sC1*integrize(kind,2,s,e,limitBreak)
            term3 = sC2*integrize(kind,3,s,e,limitBreak)
            term4 = 2*e*sC2*integrize(kind,4,s,e,limitBreak)
            
            if m > 0:
                term3 = -term3
            elif m < 0:
                term4 = -term4
            
        else:
            term1 = term2 = term3 = term4 = 0
        
        result = (1/(2*np.pi))*(p0s+term1+term2+term3+term4)
        
        if np.abs(result)>=toler:
            value = result
        
    
    return value

# Gets the specified coefficient in the specified range
def getCoefficient(em,es,eList,GRID,ekind,fileName,eitera,etoler):
    
    procesNum = 8
    
    doPrint=1
    
    part_pms = partial(p_MS,m=em,s=es,kind=ekind,shouldPrint=0,objetFile=fileName,itera=eitera,toler=etoler)
    
    yList=np.array(())
    
    with Pool(procesNum) as p:
        yList=p.map(part_pms,eList)
    
    return yList