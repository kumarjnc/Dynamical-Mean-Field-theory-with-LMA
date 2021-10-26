#!/usr/bin/python3
import numpy as np
from numpy import tile
from scipy import linalg as LA
import math,cmath
pi = math.pi
cos = math.cos
sin = math.sin
exp=np.exp   
sqrt = math.sqrt
log=np.log

class commonfn:
    def __init__(self, x, y):
        self.x=x
        self.y=y
        
    def mycomplex(self):
        z=complex(self.x, self.y)
        return z

    def sem(self):
        global s
        z=self.mycomplex()
        z1=self.mycomplex()
        z2=z**2-z1**2
        sz2=cmath.sqrt(z2)
        sem1=2.0/(z+sz2)
        sem2=2.0/(z-sz2)
        
        if(sem1.imag < 0.0):
            s=sem1
        elif(sem2.imag < 0):
            s=sem2
        else:
            print('no casual root ')
            s=0.0
        return s
    
    def theta(self):
        if(self.x < 0.0):
            theta=0.0
        elif ( self.x == 0.0):
            theta=0.5
        else:
            theta=1.0
        return theta