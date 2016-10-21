
from sympy import *
from sympy.core.sympify import sympify as parser
from sympy.solvers.ode import dsolve
import matplotlib.pyplot as plt
from numpy import linspace
class U:    #Potential nergy
    def __init__(self,forces,m,k,xx,hh):
        x=symbols('x')
        self.forces=forces
        self.m=m
        self.k=k
        self.xx=parser(xx)
        self.hh=parser(hh)
        self.X=x
        self.forces=self.forces.split(',')
        
    
        
    def expre(self,m,k):#Expression of U
        m,g,h,k,x=symbols('m g h k x')  #The symbols needed in expressions
        pre_made_forces={       #Common forces, only weight and spring force for now.
            'p':(1/2)*m*g*h,
            'el':(1/2)*k*x**2}
        print(self.forces)
        for force in self.forces:
            print('checking for ' + force)
            if force in pre_made_forces:
                print('DETECTED')
                print(pre_made_forces[force])
                self.X=self.X+pre_made_forces[force]
            else:
                pass #custom forces TODO
        self.X=self.X.subs({x:self.xx,h:self.hh,m:self.m,k:self.k,g:9.8})
        X=self.X-self.xx
        
        return X
        
class T:    #Kinetic energy
    def __init__(self,x,y,z,m):
        t=symbols('t')
        f=Function('f')
        self.x=parser(x)
        self.y=parser(y)
        self.z=parser(z)
        self.m=m
        
        
    
    def calc(self):
        t=symbols('t')
        f=Function('f')
        
        
        xd=diff(self.x,f(t)) * Derivative(f(t),t)
        yd=diff(self.y,f(t)) * Derivative(f(t),t)
        zd=diff(self.z,f(t)) * Derivative(f(t),t)
        print(self.m)
        print(type(self.m))
        print(xd)
        print(1*(xd**2))
          
        T=(1/2)*self.m*(xd**2+yd**2+zd**2)
        
        T=simplify(T)
        return T
        
    
class ode:
    def __init__(self,forces,m,k,x,y,z,xx,hh):
        t=symbols('t')
        u=U(forces,m,k,xx,hh)
        print(x)
        
        t=T(x,y,z,m)
        print(z)
        self.L=t.calc()-u.expre(m,k)
        print(self.L)
    def calc(self):
        f=Function('f')
        expr = 1.0*cos(f(t))**2*Derivative(f(t), t)**2 + 29.4*cos(f(t))
        print((self.L)==expr)
        f=Function('f')        
        print(type(f))        
        de = diff(self.L,Derivative(f(t),t))
        de=diff(de,t)
        dw=diff(self.L,f(t))
        eq=de-dw
        eq=eq.subs({cos(f(t)):1-(((f(t))**2)/2),sin(f(t)):f(t)})
        print(eq)
        eq=simplify(eq)
        print(eq)
        solution=dsolve(eq,f(t))
        print(solution)


        
        
            