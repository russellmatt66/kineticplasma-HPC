"""
Piecewise functions don't seem to be straightforward to convolve in Sympy
"""
from sympy import oo, Symbol, symbols, integrate, Piecewise
from math import fabs
# https://codereview.stackexchange.com/questions/174538/implementing-convolution-using-sympy
def convolve(f, g, t, lower_limit=-oo, upper_limit=oo):
    tau = Symbol('__very_unlikely_name__', real=True)
    return integrate(f.subs(t, tau) * g.subs(t, t - tau), 
                     (tau, lower_limit, upper_limit))

# 0th-order shape function
def sqwave(x, h, lb, rb):
    if (x < rb and x > rb):
        return h
    else:
        return 0.0

a_0 = 1.0
h_0 = 1.0 / a_0

x, xp = symbols('x xp')

S_0 = Piecewise((h_0, x < fabs(a_0 / 2.0)))

result = convolve(S_0, S_0, x)

print(result)