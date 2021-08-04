import sympy
from sympy import poly, div, degree, plot, Poly, simplify
from sympy.abc import x,a,n,k
import json
import matplotlib.pyplot as plt
import pprint
import time
from functools import wraps
from time import process_time
import json


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

def DivideByGegenbauer(f_x, n_):
	poly_f_x = Poly(f_x, x)
	f_x_coeffs  = poly_f_x.coeffs()
	f_x_degree = degree(poly_f_x)
	leading_coeff_f_x = f_x_coeffs[0]
	gegenbauer = Poly(GegenbauerPoly(n_, f_x_degree, x), x)
	gegenbauer_coeffs = gegenbauer.coeffs()
	leading_coeff_gegenbauer = gegenbauer_coeffs[0]

	quotient =  leading_coeff_f_x/ leading_coeff_gegenbauer
	remainder = poly_f_x - quotient* GegenbauerPoly(n_,f_x_degree,x)
	return (quotient, remainder)


def GenerateGegenbauerLinear(polynomials_arr ,n, i ,t):
	j = i-1
	if i == 0:
		polynomials_arr.append(1)
	elif i == 1:
		polynomials_arr.append(t)
	else:
		polynomials_arr.append(sympy.poly((2*j + n - 2)*t*polynomials_arr[j] - i* polynomials_arr[j-1]/(j + n - 2), x))
	return polynomials_arr




def GegenbauersCoeffs(g, n_):
	k = degree(g, gen = x)
	coeffs = {}
	current_g = g
	while k >= 0:
		if is_integer(current_g) and is_integer(GegenbauerPoly(n_,k,x)):
			quotient = current_g//GegenbauerPoly(n_,k,x)
			remainder = current_g%GegenbauerPoly(n_,k,x)
		else:
			quotient, remainder = DivideByGegenbauer(current_g, n_)
		current_g = remainder
		coeffs['f_{}'.format(str(k))] = quotient
		k-=1
	return coeffs

def CalculateGegenbauer():
	gegenbauer_polynomials = np.array([])
	for i in range(50):
		gegenbauer_polynomials =  np.append(gegenbauer_polynomials, sympy.polys.polytools.poly(GegenbauerPoly(n, i, x), x))
		if i%5 == 0:
			with open('gegenbauer{}.npy'.format(str(i)), 'wb') as f:
				np.save(f, gegenbauer_polynomials, allow_pickle=True)


if __name__ == "__main__":
	polynomials_arr = []
	for i in range(100):
		polynomials_arr = GenerateGegenbauerLinear(polynomials_arr, n ,i ,x)
		print(i)

	with open('gegenbauers_polynomials.txt', 'w') as f:
	    json.dump(polynomials_arr, f)










