'''
* Class object for connecting curves
'''

import numpy as np
from numpy.linalg import norm
import Bezier as bez
from mathHelper import *
from scipy.integrate import quad

class ConnectingCurve3:
	# uses a cubic Bezier curve in order to allow for continuity of first derivative
	# while allowing for some degrees of freedom; 
	# cubic is completely determined by boundary conditions
	# quartic gives us 3 additional degrees of freedom
	# we have access to, e.g., length, curvature, and torsion (note: torsion is discontinuous at boundary)
	def __init__(self, 
		end0, tan_end0, tan_mag0, 
		end1, tan_end1, tan_mag1, 
		args):

		self.rad		= args['YarnRadius']
		self.numBeads	= args['NumberBeads']
		self.bendMod	= args['BendingModulus']
		self.twistMod	= args['TwistingModulus']
		self.twist 		= args['PreferredTwist']
		self.reg 		= args['Regularization']

		self.update(end0, tan_end0, tan_mag0,
					end1, tan_end1, tan_mag1)


	def update(self, 
		end0, tan_end0, tan_mag0, 
		end1, tan_end1, tan_mag1):

		self.end0   	= end0
		self.tan_end0 	= tan_end0
		self.tan_mag0	= tan_mag0
		self.end1  		= end1
		self.tan_end1 	= tan_end1
		self.tan_mag1	= tan_mag1

		k0 = end0
		k1 = k0 + (1./3)*tan_mag0*tan_end0
		k3 = end1
		k2 = k3 - (1./3)*tan_mag1*tan_end1

		self.curve = bez.Bezier(np.array([k0,k1,k2,k3]))

	def reg_energy(self):
		ret = 0
		# functional form of the energy
		f = lambda x: 1./x

		# access curve knots
		c_knots = self.curve.k
		n_knots = len(c_knots)

		for i in range(0,n_knots-1):
			ret += f(norm([a_i - b_i for a_i, b_i in zip(c_knots[i], c_knots[i+1])],2))

		return self.reg*ret

	def energy_density(self, t):
		
		cd = self.curve.curv_deriv_at(t)
		osc = self.reg*norm(self.curve.d_at(t),2)*(norm(cd,2)**2.)
		
		return norm(self.curve.d_at(t),2)*(self.bendMod*(self.curve.curv_at(t))**2.
			+ self.twistMod*(self.curve.tors_at(t) -  self.twist)**2.) + osc

	def energy(self):
		ed = lambda x: self.energy_density(x)
		tot = quad(ed, 0, 1, limit=40,epsrel=1e-6)[0] #+ pen + curv_constraint
		return tot


class ConnectingCurve5:
	# uses a quintic Bezier curve in order to allow for continuity of curvature
	def __init__(self, 
		end0, tan_end0, tan_mag0, curv_vec0, tan_b0,
		end1, tan_end1, tan_mag1, curv_vec1, tan_b1,
		args):
		
		self.rad		= args['YarnRadius']
		self.numBeads	= args['NumberBeads']
		self.bendMod	= args['BendingModulus']
		self.twistMod	= args['TwistingModulus']
		self.twist 		= args['PreferredTwist']
		self.reg 		= args['Regularization']

		self.update(end0, tan_end0, tan_mag0, curv_vec0, tan_b0,
					end1, tan_end1, tan_mag1, curv_vec1, tan_b1)

	def update(self, 
		end0, tan_end0, tan_mag0, curv_vec0, tan_b0,
		end1, tan_end1, tan_mag1, curv_vec1, tan_b1):


		self.end0   	= end0
		self.tan_end0 	= tan_end0
		self.tan_mag0	= tan_mag0
		self.curv_vec0  = curv_vec0
		self.tan_b0		= tan_b0
		self.end1  		= end1
		self.tan_end1 	= tan_end1
		self.tan_mag1	= tan_mag1
		self.curv_vec1  = curv_vec1
		self.tan_b1		= tan_b1

		t0 = tan_end0
		t1 = tan_end1

		alpha0 = tan_mag0
		alpha1 = tan_mag1

		B0 = tan_b0
		B1 = tan_b1

		A0 = (alpha0**2)/20.
		A1 = (alpha1**2)/20.

		cv0 = curv_vec0
		cv1 = curv_vec1

		k0 = end0
		k1 = k0 + (1./5)*alpha0*t0
		k2 = k1 + A0*cv0 + B0*t0
		
		k5 = end1
		k4 = k5 - (1./5)*alpha1*t1
		k3 = k4 + A1*cv1 - B1*t1

		self.curve = bez.Bezier(np.array([k0,k1,k2,k3,k4,k5]))

	def reg_energy(self):
		ret = 0
		# functional form of the energy
		f = lambda x: 1./x

		# access curve knots
		c_knots = self.curve.k
		n_knots = len(c_knots)

		for i in range(0,n_knots-1):
			ret += f(norm([a_i - b_i for a_i, b_i in zip(c_knots[i], c_knots[i+1])],2))

		return self.reg*ret

	def energy_density(self, t):
		
		cd = self.curve.curv_deriv_at(t)
		osc = self.reg*norm(self.curve.d_at(t),2)*(norm(cd,2)**2.)
		
		#d2 = self.curve.d2_at(t)
		#d3 = self.curve.d3_at(t)
		#osc = self.reg*norm(self.curve.d_at(t),2)*(np.dot(d2,d3)**2.)
		return norm(self.curve.d_at(t),2)*(self.bendMod*(self.curve.curv_at(t))**2.
			+ self.twistMod*(self.curve.tors_at(t) -  self.twist)**2.) + osc

	def energy(self):
		ed = lambda x: self.energy_density(x)
		tot = quad(ed, 0, 1, limit=40,epsrel=1e-6)[0] #+ pen + curv_constraint
		return tot



class ConnectingCurve7:
	# uses a septic Bezier curve in order to allow for continuity of curvature
	def __init__(self, 
		end0, tan_end0, tan_mag0, curv_vec0, tan_b0, curv_deriv0, tan_c0,
		end1, tan_end1, tan_mag1, curv_vec1, tan_b1, curv_deriv1, tan_c1,
		args):
		
		self.rad		= args['YarnRadius']
		self.numBeads	= args['NumberBeads']
		self.bendMod	= args['BendingModulus']
		self.twistMod	= args['TwistingModulus']
		self.twist 		= args['PreferredTwist']
		self.reg 		= args['Regularization']

		self.update(end0, tan_end0, tan_mag0, curv_vec0, tan_b0, curv_deriv0, tan_c0,
					end1, tan_end1, tan_mag1, curv_vec1, tan_b1, curv_deriv1, tan_c1)

	def update(self, 
		end0, tan_end0, tan_mag0, curv_vec0, tan_b0, curv_deriv0, tan_c0,
		end1, tan_end1, tan_mag1, curv_vec1, tan_b1, curv_deriv1, tan_c1):


		self.end0   	 = end0
		self.tan_end0 	 = tan_end0
		self.tan_mag0	 = tan_mag0
		self.curv_vec0   = curv_vec0
		self.tan_b0		 = tan_b0
		self.curv_deriv0 = curv_deriv0
		self.tan_c0		 = tan_c0
		self.end1  		 = end1
		self.tan_end1 	 = tan_end1
		self.tan_mag1	 = tan_mag1
		self.curv_vec1   = curv_vec1
		self.tan_b1		 = tan_b1
		self.curv_deriv1 = curv_deriv1
		self.tan_c1 	 = tan_c1

		t0 = tan_end0
		t1 = tan_end1

		alpha0 = tan_mag0
		alpha1 = tan_mag1

		B0 = tan_b0
		B1 = tan_b1

		A0 = (alpha0**2)/42.
		A1 = (alpha1**2)/42.

		C0 = tan_c0
		C1 = tan_c1

		cv0 = curv_vec0
		cv1 = curv_vec1

		dc0 = curv_deriv0
		dc1 = curv_deriv1

		k0 = end0
		k1 = k0 + (1./7)*alpha0*t0
		k2 = k1 + A0*cv0 + B0*t0
		k3 = k0 + ((alpha0**3.)/210.)*(dc0 + t0*np.dot(cv0,cv0)) - (alpha0/70. - (3./5)*B0)*alpha0*cv0 + C0*t0
		
		k7 = end1
		k6 = k7 - (1./7)*alpha1*t1
		k5 = k6 + A1*cv1 - B1*t1
		k4 = k7 - ((alpha1**3.)/210.)*(dc1 + t1*np.dot(cv1,cv1)) - (alpha1/70. - (3./5)*B1)*alpha1*cv1 - C1*t1

		self.curve = bez.Bezier(np.array([k0,k1,k2,k3,k4,k5,k6,k7]))

	def reg_energy(self):
		ret = 0
		# functional form of the energy
		f = lambda x: 1./x

		# access curve knots
		c_knots = self.curve.k
		n_knots = len(c_knots)

		for i in range(0,n_knots-1):
			ret += f(norm([a_i - b_i for a_i, b_i in zip(c_knots[i], c_knots[i+1])],2))

		return self.reg*ret

	def energy_density(self, t):
		
		cd = self.curve.curv_deriv_at(t)
		osc = self.reg*norm(self.curve.d_at(t),2)*(norm(cd,2)**2.)

		return norm(self.curve.d_at(t),2)*(self.bendMod*(self.curve.curv_at(t))**2.
			+ self.twistMod*(self.curve.tors_at(t) -  self.twist)**2.) + osc

	def energy(self):
		ed = lambda x: self.energy_density(x)
		tot = quad(ed, 0, 1, limit=40,epsrel=1e-6)[0] #+ pen + curv_constraint
		return tot


