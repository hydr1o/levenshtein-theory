import sympy
from sympy import poly, div, degree, plot, Poly, simplify, factorial
from sympy.abc import x,n,k,y,i
import json
import matplotlib.pyplot as plt
import pprint
import time
from functools import wraps
from time import process_time
import json
import numpy as np
import pickle
import json
import scipy.special
from sympy import S, imageset, Interval



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
			return ((2*j + n - 2)*t*GegenbauerPoly(n,j,t) - j*GegenbauerPoly(n, j-1, t))/(j + n - 2)
		else:
			if is_integer_and_natural(n) and is_integer_and_natural(i):
				return ((2*j + n - 2)*t*GegenbauerPoly(n,j,t) - j*GegenbauerPoly(n, j-1, t))/(j + n - 2)
			else:
				raise Exception


def GegenbauerPolyWikipediaFormula(a,n,x):
	if n==0:
		return 1
	elif n ==1:
		return x
	else:
		return 1/n*(2*x*(n+a-1)*GegenbauerPolyWikipediaFormula(a, n-1, x) - (n+2*a-2)*GegenbauerPolyWikipediaFormula(a, n-2, x))



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





def GegenbauersCoeffs(g, n_):
	k = degree(g, gen = x)
	coeffs = []
	current_g = g
	while k >= 0:
		if is_integer(current_g) and is_integer(GegenbauerPoly(n_,k,x)):
			quotient = current_g//GegenbauerPoly(n_,k,x)
			remainder = current_g%GegenbauerPoly(n_,k,x)
		else:
			quotient, remainder = DivideByGegenbauer(current_g, n_)
		current_g = remainder
		coeffs.append(quotient)
		k-=1
	return coeffs

def BinomialCoeff(a, b):
	return factorial(a)/(factorial(b)*factorial(a-b))


def ComputeR_i_sympy(n, i):
	return ((n + 2*i - 2)/i)*sympy.binomial(n+i-3, i-1)

def ComputeR_i(n, i):
	return ((n + 2*i - 2)/i)*BinomialCoeff(n+i-3, i-1)
	

def CalculateT_k(x, y, k, n):
	output = 0
	for i in range(1,k+1):
		output+=ComputeR_i(n, i)*GegenbauerPoly(n,i,x)*GegenbauerPoly(n,i,y)
	return output

def LimitSphericalCode(g, n ,s):
	if not(n>=3 and s>=-1 and s<1):
		raise Exception('n must be >= 3 and s ∈ [-1; 1)')
	else:
		g_interval = imageset(x, g, Interval(-1, s))
		first_condition = g_interval.end<=0	
		gegenbauer_coeffs = GegenbauersCoeffs(g, n)
		second_condition = all(elem >= 0 for elem in gegenbauer_coeffs[1:]) and (gegenbauer_coeffs[0]>0)
		if first_condition and second_condition:
			return (sympy.poly(g).subs(1, x))/ gegenbauer_coeffs[0]
		else:
			raise Exception('conditions not completed')


def LimitSphericalDesign(g, n , τ):
	if not(n>=3 and τ>=1 and τ < sympy.poly(g).degree()):
		raise Exception('n must be >= 3 and τ >= 1')
	else:
		g_interval = imageset(x, g, Interval(-1, 1))
		first_condition = g_interval.start>=0
		gegenbauer_coeffs = GegenbauersCoeffs(g, n)
		second_condition = all(elem <= 0 for elem in gegenbauer_coeffs[τ+1:])
		print(gegenbauer_coeffs[τ+1:])
		if first_condition and second_condition:
			return (sympy.poly(g).subs(1, x))/ gegenbauer_coeffs[0]
		else:
			raise Exception('conditions not completed')

if __name__ == "__main__":
	LimitSphericalDesign(x**2 + x - 5, 3, 1)










