import sympy as sp
f=sp.Function('f')
x=sp.Symbol('x')
t=sp.Symbol('t')
'''
f''+2*delta*f'+(omega**2)*f=0
'''
def get_coeff(expr,function):
    '''
    This function finds the delta and omega of the differential equation
    '''
    expr=expr.subs({sp.Derivative(function(t),(t,2)):x**2})
    expr=expr.subs({sp.Derivative(function(t),(t,2)):x**2,sp.Derivative(function(t),(t,1)):x,function(t):1})
    #return expr
    expr=sp.Poly(expr,x)        #Turns the differential equation into a polynomial equation of order 2
    coef= expr.all_coeffs()
    
    
    '''
    The differential equation comes in the form ((d/dt)**2+2*delta*(d/dt)+omega**2)f(t)=0
    This line below cuts off the potential coefficient of (d/dt)**2
    by dividing all other coefficients by it. af''+bf'+c=0 ==> f''+b'f'+c'f=0
    Thus 'normalizing' it as in QM
    '''
    
    
    coef_normalized=[ x/coef[0] for x in coef ] 
    delta=coef_normalized[1]/2             #Delta=b/2
    omega=sp.sqrt(coef_normalized[2])      #Omega=sqrt(c)
    
    return (delta,omega)
    
if __name__ == '__main__':
    exp=40*f(t) + 3*sp.Derivative(f(t), t) + 10.0*sp.Derivative(f(t), t, t)
    print(get_coeff(exp,f))
