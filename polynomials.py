import sympy
from sympy import poly, div, degree, plot
from sympy.abc import x,a,n,k
import json
import matplotlib.pyplot as plt

def is_integer_and_natural(num):
    if num>0:
        try:
            int(num)
        except:
            return False
        else:
            return True
    else:
        return False

def is_integer(num):
	try:
		int(num)
	except:
		return False
	else:
		return True

def GegenbauerPoly(n,i,t):
	j = i-1	
	if i==0:
		return 1
	elif i==1:
		return t
	else:
		if type(n) == sympy.core.symbol.Symbol or type(t) == sympy.core.symbol.Symbol:
			return ((2*j + n - 2)*t*GegenbauerPoly(n,j,t) - i*GegenbauerPoly(n, j-1, t))/(j + n - 2)
		else:
			if is_integer_and_natural(n) and is_integer_and_natural(i):
				return ((2*j + n - 2)*t*GegenbauerPoly(n,j,t) - i*GegenbauerPoly(n, j-1, t))/(j + n - 2)
			else:
				raise Exception

def GegenbauersCoeffs(g, n):
	k = degree(g, gen = x)
	coeffs = {}
	current_g = g
	while k >= 0:
		if is_integer(current_g) and is_integer(GegenbauerPoly(n,k,x)):
			quotient = current_g//GegenbauerPoly(n,k,x)
			remainder = current_g%GegenbauerPoly(n,k,x)
		else:
			quotient, remainder = div(current_g, GegenbauerPoly(n, k, x))
		current_g = remainder
		coeffs['f_{}'.format(str(k))] = str(quotient)
		k-=1
	return coeffs



if __name__ == "__main__":
	g = 6*a*x**5 + (7*a-89)*x**4 - (a**2-1)*x**3 + (a*8)*x**2 - a*x + 10
	coefficients = GegenbauersCoeffs(g, n=2)

	for i in range(20):
		print(GegenbauerPoly(n, i, x))








