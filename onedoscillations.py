
from sympy import symbols,Function,diff,Derivative,simplify,sin,cos,Eq,init_printing
from sympy.core.sympify import sympify as parser
from sympy.solvers.ode import dsolve
from sympy.solvers.solveset import linsolve
from IPython.display import display
init_printing()


import matplotlib.pyplot as plt
from numpy import linspace

class ode:
    def __init__(self,forces,m,k,x,y,z,xx,hh,x0=0,xd0=0,alpha=0,v=''):
        '''
        This is the initialization function, The arguments are as follows :
            forces : The conservative forces applied on the system in the form of a stringed tuple
            m     : Mass, as float
            k     : Spring constant
            x,y,z : The expression of the coordinates of the mass
            xx    : Expression of the distance of the object from the rest position of the spring
            hh    : Expression of the height
            x0    : Value of the generalized coordinate at t=0
            xd0   : Value of the derivative of the generalized coordinate at t=0
            alpha : Value of the Rayleigh constant
            v     : Speed of the free side of the damper
        The generalized coordinate is written f(t)
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
        
        
        self.alpha = parser(alpha)
        print(self.alpha)
        if v=='':
            self.v=0
        else:
            self.v=parser(v)
        
        
        
        self.x0=x0
        self.xd0=xd0
        self.L=0
        
    def kinetic(self):      #Kinetic Energy of the system
        f=Function('f')
        t=symbols('t')
        ###Differentials of each coordinate by the time###
        xd = diff(self.x,t)
        yd = diff(self.y,t)
        zd = diff(self.z,t)
        #####--------############
          
        T = (xd**2+yd**2+zd**2)*(.5)*self.m # Expression of the Kinetic energy
        T = simplify(T)
        display(T)
        ##Replace cos and sin by their approximations for low amplitudes##
        T = T.replace(sin,lambda *args: args[0])
        
        T = T.replace(cos,lambda *args: args[0])
        #######-------------#########
        
        print('Kinetic Energy expression :')
        display(T)
        return T
        
    def potential(self):
        f=Function('f')
        t=symbols('t')
        m,g,h,k,x=symbols('m g h k x')  #The symbols needed in expressions
        pre_made_forces={       #Common potential forces.
            'p':(.5)*m*g*h,
            'el':(.5)*k*x**2}
        X=0 #Placebo value
        ###Checking each of the forces in the dictionary if their in the forces expression
        for force in self.forces:
            print('checking for ' + force)
            if force in pre_made_forces:
                
                X=X+pre_made_forces[force]
            else:
                pass #custom forces TODO
        
        X=X.subs({x:self.xx,h:self.hh,m:self.m,k:self.k,g:9.8})
        X=X
        
        ##Replace cos and sin by their approximations for low amplitudes##
        X=X.subs({sin(f(t)):f(t)})
        
        X=X.subs({cos(f(t)):1-(f(t))**2/2})
        #######-------------#########
        
        
        display(X)
        
        return X
        
    def calc(self): # The function that calculates the expression of f(t)
        
        f=Function('f')
        t=symbols('t')
        
        #The lagrangian
        self.L=ode.kinetic(self)-ode.potential(self)
        print('the lagrangian is : ')
        display(self.L)
        
        
        ##The different placebos for the Lagrange equation##
        de = diff(self.L,Derivative(f(t),t))
        de=diff(de,t)
        dw=diff(self.L,f(t))
        
        
        da=0.5*self.alpha*(self.v**2)   #Dissipation function
        da=da.subs({sin(f(t)):f(t)})
        
        da=da.subs({cos(f(t)):1-(f(t))**2/2})
        
        da=diff(da,Derivative(f(t),t))
        eq=de-dw + da       #LAgrange Equation
        ##############----------##############
        
        
        ##Solving the equation##
        print('The equation is : ')
        eq=simplify(eq)
        display(eq)
        solution=dsolve(eq,f(t))
        print('The solution is : ')
        display(solution)
        #######----------#######
        
        #Initial conditions:
        cnd0=Eq(solution.subs({t:0}), self.x0)
        cnd1=Eq(solution.diff(t).subs({t:0}),self.xd0)
        display(cnd0)
        
        #Solve for C1 and C2:
        C1,C2=symbols('C1 C2')  #Temporary symbols
        C1C2_sl=linsolve([cnd0,cnd1], (C1,C2))
        
        
        
        #Substitute into the Solution of the equation
        solution2=simplify(solution.subs(C1C2_sl))
        if 'C1' in str(solution2) :
            solution3=solution2.subs({C1:1})
        if 'C2' in str(solution2) :
            solution3=solution2.subs({C2:2})
        print('The solution is :')
        display(solution3)  #The Solution
        return solution3
        
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
            
        

od=ode(',el',2,200,'1*sin(f(t))','1*cos(f(t))','0','1*sin(f(t))','1*sin(f(t))')


od.graph()





        
        
            
