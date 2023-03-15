'''
* Class object for cross-over curves
'''

import numpy as np
from numpy.linalg import norm
import Bezier as bez
from mathHelper import *
from scipy.integrate import quad




class CrossOver3:
	def __init__(self, 
		pos, frame, z, phi,
		tan_frameA0, tan_magA0, tan_frameA1, tan_magA1,
		tan_frameB0, tan_magB0, tan_frameB1, tan_magB1,
		tw,
		parity, args):

		self.rad		= args['YarnRadius']
		self.numBeads	= args['NumberBeads']
		self.bendMod	= args['BendingModulus']
		self.twistMod	= args['TwistingModulus']
		self.twist 		= args['PreferredTwist']
		self.parity		= parity
		self.args       = args
		self.reg 		= args['Regularization']
		
		self.update(pos, frame, z, phi,
					tan_frameA0, tan_magA0, tan_frameA1, tan_magA1,
					tan_frameB0, tan_magB0, tan_frameB1, tan_magB1,
					tw)


	def update(self,
		pos, frame, z, phi,
		tan_frameA0, tan_magA0, tan_frameA1, tan_magA1,
		tan_frameB0, tan_magB0, tan_frameB1, tan_magB1,
		tw):

		self.pos 		 = pos
		self.frame 		 = frame
		self.z 			 = z
		self.phi		 = phi
		self.tan_frameA0 = tan_frameA0
		self.tan_magA0 	 = tan_magA0
		self.tan_frameA1 = tan_frameA1
		self.tan_magA1 	 = tan_magA1
		self.tan_frameB0 = tan_frameB0
		self.tan_magB0 	 = tan_magB0
		self.tan_frameB1 = tan_frameB1
		self.tan_magB1 	 = tan_magB1
		self.tw 		 = tw

		parity = self.parity
		par = (-1)**parity

		rot = np.dot(rotx(frame[0]),np.dot(roty(frame[1]),rotz(frame[2])))

		e3 = np.dot(rot,np.array([1.,0,0]))
		e1 = np.dot(rot,np.array([0,1.,0]))
		e2 = np.dot(rot,np.array([0,0,1.]))

		self.e3 = e3
		self.e1 = e1
		self.e2 = e2

		if par==1:
			self.E10 =  e1*np.cos(tw) + e2*np.sin(tw)
			self.E20 = -e1*np.sin(tw) + e2*np.cos(tw)
			self.E11 = e1
			self.E21 = e2
		else:
			self.E10 = e1
			self.E20 = e2
			self.E11 =  e1*np.cos(tw) + e2*np.sin(tw)
			self.E21 = -e1*np.sin(tw) + e2*np.cos(tw)

		E10 = self.E10
		E20 = self.E20
		E11 = self.E11
		E21 = self.E21

		psA0_1 = tan_frameA0[0]
		psA0_2 = tan_frameA0[1]
		psA1_1 = tan_frameA1[0]
		psA1_2 = tan_frameA1[1]

		psB0_1 = tan_frameB0[0]
		psB0_2 = tan_frameB0[1]
		psB1_1 = tan_frameB1[0]
		psB1_2 = tan_frameB1[1]

		self.tanA0 =   np.abs(np.cos(psA0_1))*E10 + par*np.abs(np.sin(psA0_1)*np.sin(psA0_2))*E20 + np.abs(np.sin(psA0_1)*np.cos(psA0_2))*e3
		self.tanA1 =  -np.abs(np.cos(psA1_1))*E11 + par*np.abs(np.sin(psA1_1)*np.sin(psA1_2))*E21 + np.abs(np.sin(psA1_1)*np.cos(psA1_2))*e3

		self.tanB0 =  -np.abs(np.cos(psB0_1))*E10 - par*np.abs(np.sin(psB0_1)*np.sin(psB0_2))*E20 + np.abs(np.sin(psB0_1)*np.cos(psB0_2))*e3
		self.tanB1 =   np.abs(np.cos(psB1_1))*E11 - par*np.abs(np.sin(psB1_1)*np.sin(psB1_2))*E21 + np.abs(np.sin(psB1_1)*np.cos(psB1_2))*e3

		kA0 = pos - z*np.cos(phi)*e3 - par*z*np.sin(phi)*E20
		kA1 = kA0 + (1./3)*tan_magA0*self.tanA0
		kA3 = pos + z*np.cos(phi)*e3 + par*z*np.sin(phi)*E21
		kA2 = kA3 - (1./3)*tan_magA1*self.tanA1

		self.kA0 = kA0
		self.kA1 = kA1
		self.kA2 = kA2
		self.kA3 = kA3

		kB0 = pos - z*np.cos(phi)*e3 + par*z*np.sin(phi)*E20
		kB1 = kB0 + (1./3)*tan_magB0*self.tanB0
		kB3 = pos + z*np.cos(phi)*e3 - par*z*np.sin(phi)*E21
		kB2 = kB3 - (1./3)*tan_magB1*self.tanB1

		self.kB0 = kB0
		self.kB1 = kB1
		self.kB2 = kB2
		self.kB3 = kB3

		self.curve = np.array([bez.Bezier(np.array([kA0,kA1,kA2,kA3])),
							   bez.Bezier(np.array([kB0,kB1,kB2,kB3]))])

		#self.find_contacts()

	'''
	def print_beads(self):
		# prints position of bead and radius of bead on each line
		# alter to include a third output indicating state
		ret = ""
		for j in range(0,2):
			for i in range(0,self.numBeads+1):
				ret += vec_print(self.curve[j].at(i/self.numBeads))
				ret += ", " + str(self.rad) +  "\n"
		return ret


	def find_contacts(self):
		self.overlap_loc = np.empty(shape=[0,2])
		self.overlap_val = np.array([])
		for i in range(0,self.numBeads+1):
			for j in range(0,self.numBeads+1):
				t1 = i/self.numBeads
				t2 = j/self.numBeads
				dist = norm(self.curve[0].at(t1) - self.curve[1].at(t2),2)
				if dist < 2*self.rad:
					self.overlap_loc = np.append(self.overlap_loc,[[t1,t2]], axis=0)
					self.overlap_val = np.append(self.overlap_val,[2*self.rad-dist], axis=0)
		return self.overlap_loc


	def print_beads_contact(self):
		# prints position of beads that are in contact and radius of bead on each line
		# alter to include a third output indicating state
		ret = ""
		it =  0
		maxval = 0
		for [t1,t2] in self.overlap_loc:
			if it==0:
				maxval = np.abs(self.overlap_val[0])
			elif it==len(self.overlap_loc)-1:
				ret += vec_print(self.curve[0].at(self.overlap_loc[it,0])) + ", " + str(self.rad) + ", " + str(maxval) + "\n"
				ret += vec_print(self.curve[1].at(1-self.overlap_loc[it,0])) + ", " + str(self.rad) + ", " + str(maxval) + "\n"
			elif it>0:
				if self.overlap_loc[it,0] == self.overlap_loc[it-1,0] and np.abs(self.overlap_val[it]) > maxval:
					maxval = np.abs(self.overlap_val[it])
				elif self.overlap_loc[it,0] != self.overlap_loc[it-1,0]:
					ret += vec_print(self.curve[0].at(self.overlap_loc[it-1,0])) + ", " + str(self.rad) + ", " + str(maxval) + "\n"
					ret += vec_print(self.curve[1].at(1-self.overlap_loc[it-1,0])) + ", " + str(self.rad) + ", " + str(maxval) + "\n"
					maxval = np.abs(self.overlap_val[it])
			it+=1

		return ret


	def contact_energy(self, floof_fct):
		ret = 0
		for i in range(0,self.numBeads+1):
			for j in range(0,self.numBeads+1):
				t1 = i/self.numBeads
				t2 = j/self.numBeads
				dist = norm(self.curve[0].at(t1) - self.curve[1].at(t2),2) - 2*self.rad
				if dist < 0:
					ret += floof_fct(np.abs(dist))*norm(self.curve[0].d_at(t1))*norm(self.curve[1].d_at(t2))/(self.numBeads**2.)
		return ret

	def d_contact_energy(self, n):
		k = self.hertzMod
		ret = 0
		it = 0
		ntot = 4
		e1 = self.e1
		e2 = self.e2
		e3 = self.e3
		for t1,t2 in self.overlap_loc:
			vec_diff = self.curve[0].at(t1) -  self.curve[1].at(t2)
			uvec = vec_diff/norm(vec_diff,2)
			ret += -k*(self.overlap_val[it]**(3./2))*(bez.Brn(t1,n,4)*uvec 
				- bez.Brn(t2,n,4)*(2*e3*np.dot(e3,uvec) - uvec))
			it+=1
		return ret


	def reg_energy(self):
		ret = 0
		# functional form of the energy
		f = lambda x: 1./x

		# access curve knots
		c0_knots = self.curve[0].k
		c1_knots = self.curve[1].k
		n_knots = len(c0_knots)

		for i in range(0,n_knots-1):
			ret += f(norm([a_i - b_i for a_i, b_i in zip(c0_knots[i], c0_knots[i+1])],2))
			ret += f(norm([a_i - b_i for a_i, b_i in zip(c1_knots[i], c1_knots[i+1])],2))

		return self.reg*ret
	'''


	def energy(self):
		energy_density = lambda x: (
			norm(self.curve[0].d_at(x),2)*(self.bendMod*(self.curve[0].curv_at(x))**2. 
			+ self.twistMod*(self.curve[0].tors_at(x) -  self.twist)**2.
			+ self.reg*(norm(self.curve[0].curv_deriv_at(x))**2.) )
			+ norm(self.curve[1].d_at(x),2)*(self.bendMod*(self.curve[1].curv_at(x))**2. 
			+ self.twistMod*(self.curve[1].tors_at(x) -  self.twist)**2.
			+ self.reg*(norm(self.curve[1].curv_deriv_at(x))**2.) )
			)
		return quad(energy_density, 0, 1, limit=40,epsrel=1e-6)[0] #+ self.reg_energy()






class CrossOver5:
	def __init__(self, 
		pos, frame, z, phi,
		tan_frameA0, tan_magA0, curv_A0, normal_A0, tan_mag2A0,
		tan_frameA1, tan_magA1, curv_A1, normal_A1, tan_mag2A1,
		tan_frameB0, tan_magB0, curv_B0, normal_B0, tan_mag2B0,
		tan_frameB1, tan_magB1, curv_B1, normal_B1, tan_mag2B1,
		tw,
		parity, args):

		self.rad		= args['YarnRadius']
		self.numBeads	= args['NumberBeads']
		self.bendMod	= args['BendingModulus']
		self.twistMod	= args['TwistingModulus']
		self.twist 		= args['PreferredTwist']
		self.parity		= parity
		self.args       = args
		self.reg 		= args['Regularization']
		
		self.update(pos, frame, z, phi,
					tan_frameA0, tan_magA0, curv_A0, normal_A0, tan_mag2A0,
					tan_frameA1, tan_magA1, curv_A1, normal_A1, tan_mag2A1,
					tan_frameB0, tan_magB0, curv_B0, normal_B0, tan_mag2B0,
					tan_frameB1, tan_magB1, curv_B1, normal_B1, tan_mag2B1,
					tw)


	def update(self,
		pos, frame, z, phi,
		tan_frameA0, tan_magA0, curv_A0, normal_A0, tan_mag2A0,
		tan_frameA1, tan_magA1, curv_A1, normal_A1, tan_mag2A1,
		tan_frameB0, tan_magB0, curv_B0, normal_B0, tan_mag2B0,
		tan_frameB1, tan_magB1, curv_B1, normal_B1, tan_mag2B1,
		tw):

		self.pos 		 = pos
		self.frame 		 = frame
		self.z 			 = z
		self.phi		 = phi
		self.tan_frameA0 = tan_frameA0
		self.tan_magA0 	 = tan_magA0
		self.curv_A0	 = curv_A0
		self.normal_A0	 = normal_A0
		self.tan_mag2A0	 = tan_mag2A0
		self.tan_frameA1 = tan_frameA1
		self.tan_magA1 	 = tan_magA1
		self.curv_A1	 = curv_A1
		self.normal_A1	 = normal_A1
		self.tan_mag2A1	 = tan_mag2A1
		self.tan_frameB0 = tan_frameB0
		self.tan_magB0 	 = tan_magB0
		self.curv_B0	 = curv_B0
		self.normal_B0	 = normal_B0
		self.tan_mag2B0	 = tan_mag2B0
		self.tan_frameB1 = tan_frameB1
		self.tan_magB1 	 = tan_magB1
		self.curv_B1	 = curv_B1
		self.normal_B1	 = normal_B1
		self.tan_mag2B1	 = tan_mag2B1
		self.tw 		 = tw

		parity = self.parity
		par = (-1)**parity

		rot = np.dot(rotx(frame[0]),np.dot(roty(frame[1]),rotz(frame[2])))

		e3 = np.dot(rot,np.array([1.,0,0]))
		e1 = np.dot(rot,np.array([0,1.,0]))
		e2 = np.dot(rot,np.array([0,0,1.]))

		self.e3 = e3
		self.e1 = e1
		self.e2 = e2

		if par==1:
			self.E10 =  e1*np.cos(tw) + e2*np.sin(tw)
			self.E20 = -e1*np.sin(tw) + e2*np.cos(tw)
			self.E11 = e1
			self.E21 = e2
		else:
			self.E10 = e1
			self.E20 = e2
			self.E11 =  e1*np.cos(tw) + e2*np.sin(tw)
			self.E21 = -e1*np.sin(tw) + e2*np.cos(tw)

		E10 = self.E10
		E20 = self.E20
		E11 = self.E11
		E21 = self.E21

		psA0_1 = tan_frameA0[0]
		psA0_2 = tan_frameA0[1]
		psA1_1 = tan_frameA1[0]
		psA1_2 = tan_frameA1[1]

		psB0_1 = tan_frameB0[0]
		psB0_2 = tan_frameB0[1]
		psB1_1 = tan_frameB1[0]
		psB1_2 = tan_frameB1[1]

		self.tanA0 =   np.abs(np.cos(psA0_1))*E10 + par*np.abs(np.sin(psA0_1)*np.sin(psA0_2))*E20 + np.abs(np.sin(psA0_1)*np.cos(psA0_2))*e3
		self.tanA1 =  -np.abs(np.cos(psA1_1))*E11 + par*np.abs(np.sin(psA1_1)*np.sin(psA1_2))*E21 + np.abs(np.sin(psA1_1)*np.cos(psA1_2))*e3

		self.tanB0 =  -np.abs(np.cos(psB0_1))*E10 - par*np.abs(np.sin(psB0_1)*np.sin(psB0_2))*E20 + np.abs(np.sin(psB0_1)*np.cos(psB0_2))*e3
		self.tanB1 =   np.abs(np.cos(psB1_1))*E11 - par*np.abs(np.sin(psB1_1)*np.sin(psB1_2))*E21 + np.abs(np.sin(psB1_1)*np.cos(psB1_2))*e3

		self.cvA0 = curv_A0*(np.cross(E10,self.tanA0)*np.cos(normal_A0)/np.abs(np.sin(psA0_1)) + np.cross(self.tanA0,np.cross(E10,self.tanA0))*np.sin(normal_A0)/np.abs(np.sin(psA0_1)))
		self.cvA1 = curv_A1*(np.cross(E11,self.tanA1)*np.cos(normal_A1)/np.abs(np.sin(psA1_1)) + np.cross(self.tanA1,np.cross(E11,self.tanA1))*np.sin(normal_A1)/np.abs(np.sin(psA1_1)))

		self.cvB0 = curv_B0*(np.cross(E10,self.tanB0)*np.cos(normal_B0)/np.abs(np.sin(psB0_1)) + np.cross(self.tanB0,np.cross(E10,self.tanB0))*np.sin(normal_B0)/np.abs(np.sin(psB0_1)))
		self.cvB1 = curv_B1*(np.cross(E11,self.tanB1)*np.cos(normal_B1)/np.abs(np.sin(psB1_1)) + np.cross(self.tanB1,np.cross(E11,self.tanB1))*np.sin(normal_B1)/np.abs(np.sin(psB1_1)))

		kA0 = pos - z*np.cos(phi)*e3 - par*z*np.sin(phi)*E20
		kA1 = kA0 + (1./5)*tan_magA0*self.tanA0
		kA2 = kA1 + ( (tan_magA0**2)/20. )*self.cvA0 + tan_mag2A0*self.tanA0

		kA5 = pos + z*np.cos(phi)*e3 + par*z*np.sin(phi)*E21
		kA4 = kA5 - (1./5)*tan_magA1*self.tanA1
		kA3 = kA4 + ( (tan_magA1**2)/20. )*self.cvA1 - tan_mag2A1*self.tanA1

		self.kA0 = kA0
		self.kA1 = kA1
		self.kA2 = kA2
		self.kA3 = kA3
		self.kA4 = kA4
		self.kA5 = kA5

		kB0 = pos - z*np.cos(phi)*e3 + par*z*np.sin(phi)*E20
		kB1 = kB0 + (1./5)*tan_magB0*self.tanB0
		kB2 = kB1 + ( (tan_magB0**2)/20. )*self.cvB0 + tan_mag2B0*self.tanB0

		kB5 = pos + z*np.cos(phi)*e3 - par*z*np.sin(phi)*E21
		kB4 = kB5 - (1./5)*tan_magB1*self.tanB1
		kB3 = kB4 + ( (tan_magB1**2)/20. )*self.cvB1 - tan_mag2B1*self.tanB1

		self.kB0 = kB0
		self.kB1 = kB1
		self.kB2 = kB2
		self.kB3 = kB3
		self.kB4 = kB4
		self.kB5 = kB5

		self.curve = np.array([bez.Bezier(np.array([kA0,kA1,kA2,kA3,kA4,kA5])),
							   bez.Bezier(np.array([kB0,kB1,kB2,kB3,kB4,kB5]))])

		#self.find_contacts()

	'''
	def print_beads(self):
		# prints position of bead and radius of bead on each line
		# alter to include a third output indicating state
		ret = ""
		for j in range(0,2):
			for i in range(0,self.numBeads+1):
				ret += vec_print(self.curve[j].at(i/self.numBeads))
				ret += ", " + str(self.rad) +  "\n"
		return ret


	def find_contacts(self):
		self.overlap_loc = np.empty(shape=[0,2])
		self.overlap_val = np.array([])
		for i in range(0,self.numBeads+1):
			for j in range(0,self.numBeads+1):
				t1 = i/self.numBeads
				t2 = j/self.numBeads
				dist = norm(self.curve[0].at(t1) - self.curve[1].at(t2),2)
				if dist < 2*self.rad:
					self.overlap_loc = np.append(self.overlap_loc,[[t1,t2]], axis=0)
					self.overlap_val = np.append(self.overlap_val,[2*self.rad-dist], axis=0)
		return self.overlap_loc


	def print_beads_contact(self):
		# prints position of beads that are in contact and radius of bead on each line
		# alter to include a third output indicating state
		ret = ""
		it =  0
		maxval = 0
		for [t1,t2] in self.overlap_loc:
			if it==0:
				maxval = np.abs(self.overlap_val[0])
			elif it==len(self.overlap_loc)-1:
				ret += vec_print(self.curve[0].at(self.overlap_loc[it,0])) + ", " + str(self.rad) + ", " + str(maxval) + "\n"
				ret += vec_print(self.curve[1].at(1-self.overlap_loc[it,0])) + ", " + str(self.rad) + ", " + str(maxval) + "\n"
			elif it>0:
				if self.overlap_loc[it,0] == self.overlap_loc[it-1,0] and np.abs(self.overlap_val[it]) > maxval:
					maxval = np.abs(self.overlap_val[it])
				elif self.overlap_loc[it,0] != self.overlap_loc[it-1,0]:
					ret += vec_print(self.curve[0].at(self.overlap_loc[it-1,0])) + ", " + str(self.rad) + ", " + str(maxval) + "\n"
					ret += vec_print(self.curve[1].at(1-self.overlap_loc[it-1,0])) + ", " + str(self.rad) + ", " + str(maxval) + "\n"
					maxval = np.abs(self.overlap_val[it])
			it+=1

		return ret


	def contact_energy(self, floof_fct):
		ret = 0
		for i in range(0,self.numBeads+1):
			for j in range(0,self.numBeads+1):
				t1 = i/self.numBeads
				t2 = j/self.numBeads
				dist = norm(self.curve[0].at(t1) - self.curve[1].at(t2),2) - 2*self.rad
				if dist < 0:
					ret += floof_fct(np.abs(dist))*norm(self.curve[0].d_at(t1))*norm(self.curve[1].d_at(t2))/(self.numBeads**2.)
		return ret

	def d_contact_energy(self, n):
		k = self.hertzMod
		ret = 0
		it = 0
		ntot = 4
		e1 = self.e1
		e2 = self.e2
		e3 = self.e3
		for t1,t2 in self.overlap_loc:
			vec_diff = self.curve[0].at(t1) -  self.curve[1].at(t2)
			uvec = vec_diff/norm(vec_diff,2)
			ret += -k*(self.overlap_val[it]**(3./2))*(bez.Brn(t1,n,4)*uvec 
				- bez.Brn(t2,n,4)*(2*e3*np.dot(e3,uvec) - uvec))
			it+=1
		return ret


	def reg_energy(self):
		ret = 0
		# functional form of the energy
		f = lambda x: 1./x

		# access curve knots
		c0_knots = self.curve[0].k
		c1_knots = self.curve[1].k
		n_knots = len(c0_knots)

		for i in range(0,n_knots-1):
			ret += f(norm([a_i - b_i for a_i, b_i in zip(c0_knots[i], c0_knots[i+1])],2))
			ret += f(norm([a_i - b_i for a_i, b_i in zip(c1_knots[i], c1_knots[i+1])],2))

		return self.reg*ret
	'''


	def energy(self):
		energy_density = lambda x: (
			norm(self.curve[0].d_at(x),2)*(self.bendMod*(self.curve[0].curv_at(x))**2. 
			+ self.twistMod*(self.curve[0].tors_at(x) -  self.twist)**2.
			+ self.reg*(norm(self.curve[0].curv_deriv_at(x))**2.) )
			+ norm(self.curve[1].d_at(x),2)*(self.bendMod*(self.curve[1].curv_at(x))**2. 
			+ self.twistMod*(self.curve[1].tors_at(x) -  self.twist)**2.
			+ self.reg*(norm(self.curve[1].curv_deriv_at(x))**2.) )
			)
		return quad(energy_density, 0, 1, limit=40,epsrel=1e-6)[0] #+ self.reg_energy()

