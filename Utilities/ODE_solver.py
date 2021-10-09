import sympy as sp
from sympy.solvers.solveset import linsolve
from IPython.display import display
sp.init_printing()


t,C1,C2,phi,A,u,v,dephasing=sp.symbols('t C1 C2 phi A u v dephasing')
f=sp.Function('f')

def odesolver(coef):
    '''
    The solver of the f"+2*deltaf'+omega^2f =0
    '''
    delta=coef[0]
    omega=coef[1]
    
    if delta > omega:
        case=True
        s1=-delta+sp.sqrt(delta**2-omega**2)
        s2=-delta-sp.sqrt(delta**2-omega**2)       
        sol=sp.Eq(f(t),C1*sp.exp(s1*t)+C2*sp.exp(s2*t))
    elif delta==omega:
        case=True
        sol=sp.Eq(f(t),(C1+C2*t)*sp.exp(-delta*t))
    elif delta < omega:
        case=False
        omega_d=sp.sqrt(omega**2-delta**2)
        sol=sp.Eq(f(t),C1*sp.exp(-delta*t)*sp.cos((omega_d*t)+C2))
    return (sol,case)
    
def ps_periodic_csts(coef,iniconditions=(0,1)):
    
    
    
    delta=coef[0]
    omega=coef[1]
    omega_d=sp.sqrt(omega**2-delta**2)
    #print('OMEGA HERE',omega_d)
    expr1=sp.Eq(A*u,iniconditions[0])
    expr2=sp.Eq(-omega_d*A*v-delta*A*u,iniconditions[1])
    expr3=sp.Eq(u**2+v**2,1)
    solu=sp.solve([expr1,expr2,expr3],(A,u,v))
    if len(solu)>1:
        taken_solution=solu[1]
        uu=taken_solution[1]
        vv=taken_solution[2]
        if uu !=0:
            dephasing=sp.solve(sp.Eq(sp.cos(phi),uu),phi)[0]
        elif vv!=0:
            dephasing=sp.solve(sp.Eq(sp.sin(phi),vv),phi)[0]
            #print(sp.simplify(dephasing))
            #dephasing=sp.pi*(dephasing/sp.pi)
    elif len(solu)==1:
        taken_solution=solu[0]
        uu=taken_solution[1]
        vv=taken_solution[2]
        if uu !=0:
            dephasing=sp.solveset(sp.Eq(sp.cos(phi),uu),phi)
        elif vv!=0:
            dephasing=sp.solveset(sp.Eq(sp.sin(phi),vv),phi)
    
    dephasing=sp.pi/round(1/float(dephasing/sp.pi))
    finalized_equation=A*sp.exp(-delta*t)*sp.cos((omega_d*t)+phi)
    finalized_equation=finalized_equation.subs({A:taken_solution[0],phi:dephasing})
    return ((omega_d,dephasing),finalized_equation)
    
if __name__=='__main__':
    print(ps_periodic_csts((0,0.2234))[1])
