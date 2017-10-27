
from sympy import symbols,Function,diff,Derivative,simplify,sin,cos,Eq,init_printing,sqrt
from sympy.core.sympify import sympify as parser
import Utilities.ODE_solver as dsolv
from sympy.solvers.solveset import linsolve
from IPython.display import display
import Utilities.coefficient_finder as coeff
from sympy.solvers.ode import dsolve


init_printing()


import matplotlib.pyplot as plt
from numpy import linspace

class ode:
    def __init__(self,forces,m,k,x,y,z,xx,hh,equil=0,x0=0,xd0=1,alpha=0,v=''):
        '''
        This is the initialization function, The arguments are as follows :
            forces : The conservative forces applied on the system in the form of a stringed tuple
            m      : Mass(es), as float
            k      : Spring constant(s)
            x,y,z  : The expression(s) of the coordinates of the mass(s)
            xx     : Expression of the distance of the object from the rest position of the spring
            hh     : Expression of the height
            equil  : Value of f(t) in the equilibrium point
            x0     : Value of the generalized coordinate at t=0
            xd0    : Value of the derivative of the generalized coordinate at t=0
            alpha  : Value of the constant(s) of the Rayleigh function
            v      : Speed(s) of the free side(s) of the damper(s)
        The generalized coordinate is written f(t)
        In case of multiple masses,precise each mass's coordinates and weight(if the mass has no weight it's height is 0)
        Same with Alpha,v k,xx !!!Order matters!!!
        '''
        t=symbols('t')
        
        self.forces=forces.split(',')
        self.k=k
        self.y=parser(y)
        self.z=parser(z)
        self.xx=parser(xx)
        display(xx)
        self.hh=parser(hh)
                
        self.x=parser(x)
        self.m=m
        #print(type(self.x))
        arggs=[y,z,x,m]
        if 'p' in self.forces:
            arggs.append(hh)
            
           
        if all(isinstance(arg,(list,tuple))for arg in arggs):
            it=iter(arggs)
            the_len=len(next(it))
            if  all(len(1)== the_len for l in it):
                print('Checking complete for coordinates')
        elif all(isinstance(arg,(str,float,int))for arg in arggs):
            print('Checking complete for coordinates')
        else:
            raise ValueError('The arguments aren\'t of the same length, please check again')
        
            
            
        if 'el' in self.forces:
            arggs=(k,xx)
            
            if all(isinstance(arg,(list,tuple))for arg in arggs):
                
                if  len(xx)==len(k):
                    print('Checking complete for elastic arguments')
                else:
                    raise ValueError('The arguments of the coordinates aren\'t of the same length, please check again')
            elif all(isinstance(arg,(str,float,int))for arg in arggs):
                print('Checking completefor elastic arguments')
            else:
                raise ValueError('The arguments k and xx aren\'t of the same length, please check again')
                
                
        if alpha is list or alpha is tuple:
            self.alpha=list(alpha)
        else:
            self.alpha=parser(alpha)
        
        if v=='':
            self.v=0
        else:
            if v is list or v is tuple:
                self.v=list(v)
            else:
                self.v=parser(v)
        if alpha is list and v is list:
            if len(alpha) != len(v):
                raise ValueError('The arguments alpha and v aren\'t of the same length')
            else:
                print('Cheching complete for Rayleigh Dissipation')
        elif alpha is list or v is list:
            raise ValueError('The arguments alpha and v aren\'t of the same type')
        
            
        if self.v !=0:
            self.v = self.v.replace(sin,lambda *args: args[0])
        
            self.v = self.v.replace(cos,lambda *args: args[0])
        
        
        
        self.x0=x0
        self.xd0=xd0
        self.L=0
        
        
        
        
    def kinetic(self):      #Kinetic Energy of the system
        f=Function('f')
        t=symbols('t')
        ###Differentials of each coordinate by the time###
        if not self.x is list:
            xd = diff(self.x,t)
            yd = diff(self.y,t)
            zd = diff(self.z,t)
        else:
            xd=[] 
            yd=[]
            zd=[]
            for xc in self.x:
                xd.append(diff(xc,t))
            for yc in self.y:
                yd.append(diff(yc,t))
            for zc in self.z:
                zd.append(diff(zc,t))
        #####--------############
        
        if xd is list: 
            T=[]
            for l in len(self.x):
                T.append(simplify((xd[l]**2+yd[l]**2+zd[l]**2)*(1/2)*self.m[l]))  # Expression of the Kinetic energy
        else:
            T = (xd**2+yd**2+zd**2)*(1/2)*self.m
            T=simplify(T)
        
        ##Replace cos and sin by their approximations for low amplitudes##
        if T is list:
            for tt in T:
                tt = tt.replace(sin,lambda *args: args[0])
        
                tt = tt.replace(cos,lambda *args: args[0])
        else:
            T = T.replace(sin,lambda *args: args[0])
        
            T = T.replace(cos,lambda *args: args[0])
        #######-------------#########
        
        print('Kinetic Energy expression :')
        display(T)
        return T
        
        
    def potential(self):
        f=Function('f')
        t=symbols('t')
        m,g,h,k,x,l0=symbols('m g h k x l0')  #The symbols needed in expressions
        pre_made_forces={       #Common potential forces.
            'p':(1/2)*m*g*h,
            'el':(1/2)*k*x**2}
        xel=[] #Placebo values
        xp=[]
        X=0
        if 'el' in self.forces:
            if self.k is list:
                for i in range(len(self.k)):
                    xel.append(pre_made_forces['el'])
            else:
                xel.append(pre_made_forces['el'])
            if self.k is list:
                for i in range(len(self.k)):
                    temp=xel[i]
                    temp=temp.subs({k:self.k[i],x:self.xx[i]})
                    X+=temp
            else:
                temp=xel[0]
                temp=temp.subs({k:self.k,x:self.xx})
                X+=temp
        if 'p' in self.forces:
            if self.m is list:
                for i in range(len(self.m)):
                    xp.append(pre_made_forces['p'])
            else:
                xp.append(pre_made_forces['p'])
            if self.m is list:
                for i in range(len(self.m)):
                    temp=xp[i]
                    temp=temp.subs({m:self.m[i],h:self.hh[i],g:9.8})
                    X+=temp
                else:
                    temp=xp[0]
                    temp=temp.subs({m:self.m,h:self.hh,g:9.8})
                    X+=temp
        
        
        ##Replace cos and sin by their approximations for low amplitudes##
        X=X.subs({sin(f(t)):f(t)})
        
        X=X.subs({cos(f(t)):1-(f(t))**2/2})
        
        #######-------------#########
        
        print('Potential energy:')
        display(X)
        
        return X
        
    def calc(self): # The function that calculates the expression of f(t)
        
        f=Function('f')
        t,delta,omega_0,phi,A,B=symbols('t delta omega_0 phi A B')
        
        #The lagrangian
        self.L=ode.kinetic(self)-ode.potential(self)
        print('the lagrangian is : ')
        display(self.L)
        
        
        ##The different placebos for the Lagrange equation##
        de = diff(self.L,Derivative(f(t),t))
        de=diff(de,t)
        dw=diff(self.L,f(t))
        
        
        da=0.5*self.alpha*(self.v**2)#Dissipation function
        da=da.subs({sin(f(t)):f(t)})
        da=da.subs({cos(f(t)):1-(f(t))**2/2})
        
        da=diff(da,Derivative(f(t),t))
        eq=de-dw + da       #LAgrange Equation
        ##############----------##############
        
        
        
        
        #######----------#######
        
        delta_v,omega_v=coeff.get_coeff(eq,f)
        display(sp.Eq(delta,delta_v))
        display(sp.Eq(omega_0,omega_v))
        
        solution,case =  dsolv.odesolver((delta_v,omega_v))
        print
        #Initial conditions:
        cnd0=Eq(solution.subs({t:0}), self.x0)
        cnd1=Eq(solution.diff(t).subs({t:0}),self.xd0)
        
        if case==True:   #Solve for C1 and C2:
            
                C1,C2=symbols('C1 C2')  #Temporary symbols
                C1C2_sl=linsolve([cnd0,cnd1], (C1,C2))
        
        
                #Substitute into the Solution of the equation
                solution3=simplify(solution.subs(C1C2_sl))
                
                #Sympy's linsolve might return a constant as a function of the other,
                #In such case the previous line will only sub the one constant as a function of the other
                #Thus the necessity or replacing the variable constant by a real number 
                # (1 because it's what is used usually in exercises)
                
                if 'C1' in str(solution3):
                    solution3=solution3.subs({C1:1})
                elif 'C2' in str(solution3):
                    solution3=solution3.subs({C2:1})
        else:
                f_ret=dsolv.ps_periodic_csts((delta_v,omega_v),(self.x0,self.xd0))
                solution3=f_ret[1]
                    
        print('The solution is :')
        display(solution3)  #The Solution
        
        if case==False:
            print('The pulsation of this oscillation is :')
            display(f_ret[0][0])
            print(f_ret[0][1])
        return solution
        
    def graph(self):        #A graph function
        lspace=linspace(0,100,200)  #Temporary array
        t=symbols('t')
        f_t=[]
        solution=od.calc()
        solution=str(solution)  #The only way I could find to seperate the two parts of the equation
        solution=solution[9:-1]
        
        solution=parser(solution)

        
        for number in lspace:
            f_t.append(solution.subs(t,number).evalf())
        
        plt.plot(lspace,f_t)
        plt.show()
            
        
if __name__ == '__main__':
    m,k,l=symbols('m k l')
    
    od=ode(',el',10,0.5,'1*sin(f(t))','1*cos(f(t))','0','1*sin(f(t))','0')

    od.calc()
    


        
        
            