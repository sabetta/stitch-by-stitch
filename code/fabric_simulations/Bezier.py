'''
* Some Bezier curve helper functions
'''

import numpy as np
from scipy.special import binom
from numpy.linalg import norm
from scipy.integrate import quad

def Brn(t,nu,n):
	# Bernstein basis polynomials
	if nu < 0 or nu > n:
		return 0
	else:
		return binom(n,nu)*((1-t)**(n-nu))*(t**nu)

def d_Brn(t,nu,n):
	return n*(Brn(t,nu-1,n-1) - Brn(t,nu,n-1))

def d2_Brn(t,nu,n):
	return n*(d_Brn(t,nu-1,n-1) - d_Brn(t,nu,n-1))

def d3_Brn(t,nu,n):
	return n*(d2_Brn(t,nu-1,n-1) - d2_Brn(t,nu,n-1))


class Bezier:
	def __init__(self, knots):
		self.k = knots
		self.ntot = len(knots)-1

	def at(self, t):
		tot = 0
		n = 0
		ntot = self.ntot
		for p in self.k:
			tot += Brn(t,n,ntot)*p
			n += 1
		return(tot)

	def d_at(self, t):
		tot = 0
		n = 0
		ntot = len(self.k)-1
		for p in self.k:
			tot += d_Brn(t,n,ntot)*p
			n += 1
		return(tot)

	def d2_at(self, t):
		tot = 0
		n = 0
		ntot = len(self.k)-1
		for p in self.k:
			tot += d2_Brn(t,n,ntot)*p
			n += 1
		return(tot)

	def d3_at(self, t):
		tot = 0
		n = 0
		ntot = len(self.k)-1
		for p in self.k:
			tot += d3_Brn(t,n,ntot)*p
			n += 1
		return(tot)


	def tan_at(self,t):
		dr = self.d_at(t)
		return dr/norm(dr,2)


	def curv_at(self, t):
		dr 	= self.d_at(t)
		d2r = self.d2_at(t)
		return norm(np.cross(dr,d2r))/(norm(dr,2)**3.)


	def tors_at(self, t):
		d3r = self.d3_at(t)
		cr  = np.cross(self.d_at(t),self.d2_at(t))
		if norm(cr,2) != 0:
			return np.dot(cr,d3r)/(norm(cr,2)**2.)
		else:
			return 0

	def curv_vec_at(self, t):
		dr 	 = self.d_at(t)
		d2r  = self.d2_at(t)
		tan  = self.tan_at(t)
		return ( (d2r - np.dot(d2r,tan)*tan)/(norm(dr,2)**2.) )


	def curv_deriv_at(self, t):
		dr 		= self.d_at(t)
		d2r 	= self.d2_at(t)
		d3r 	= self.d3_at(t)
		tan 	= self.tan_at(t)
		cv 		= self.curv_vec_at(t)

		ret = (d3r - np.dot(tan,d3r)*tan)/(norm(dr,2)**3.)
		ret += -3*cv*np.dot(tan,d2r)/(norm(dr,2)**2.) - tan*np.dot(cv,cv)
		return ret

	def d_len_at(self,t,n):
		ntot = self.ntot
		return d_Brn(t,n,ntot)*self.tan_at(t)


	def d_curv_at(self, t, n):
		ntot = self.ntot
		dr   = self.d_at(t)
		d2r  = self.d2_at(t)
		cr   = np.cross(dr,d2r)

		ret = (1./((norm(dr,2)**3.)*norm(cr,2)))*(
			dr*(d_Brn(t,n,ntot)*(norm(d2r,2)**2.) - d2_Brn(t,n,ntot)*np.dot(dr,d2r))
			-d2r*(d_Brn(t,n,ntot)*np.dot(dr,d2r) - d2_Brn(t,n,ntot)*(norm(dr,2)**2.)))
		ret += -3*self.curv_at(t)*self.d_len_at(t,n)/norm(dr,2)
		return ret


	def d_tors_at(self, t, n):
		ntot =  self.ntot
		dr   = self.d_at(t)
		d2r  = self.d2_at(t)
		d3r  = self.d3_at(t)
		cr12 = np.cross(dr,d2r)
		cr23 = np.cross(d2r,d3r)
		cr31 = np.cross(d3r,dr)
		tors = self.tors_at(t)

		ret = (1./(norm(cr12,2)**2.))*(d_Brn(t,n,ntot)*cr23 + d2_Brn(t,n,ntot)*cr31 + d3_Brn(t,n,ntot)*cr12
			-2*tors*(dr*(d_Brn(t,n,ntot)*(norm(d2r,2)**2.) - d2_Brn(t,n,ntot)*np.dot(dr,d2r))
				-d2r*(d_Brn(t,n,ntot)*np.dot(dr,d2r) - d2_Brn(t,n,ntot)*(norm(dr,2)**2.))))
		return ret


	def length(self):
		return quad(lambda x: norm(self.d_at(x),2), 0, 1)[0]

	def print(self):
		ret = ""
		it =  0
		for p in self.k:
			for comp in p:
				ret += str(comp)
				it +=  1
				if it != self.k.size:
					ret += ", "
		ret += "\n"
		return ret



# returns a string formatted for printing to file
# reuquires an array of knots
def format_print(knots):
	ret = ""
	it =  0
	for k in knots:
		for comp in k:
			ret += str(comp)
			it +=  1
			if it != knots.size:
				ret += ", "
	ret += "\n"
	return ret