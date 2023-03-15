'''
* Some mathematical helper functions
'''

import numpy as np

def rotx(theta):
	return np.array([
		[1., 0. ,0.],
		[0.,  np.cos(theta), -np.sin(theta)],
		[0.,  np.sin(theta),  np.cos(theta)]])

def roty(theta):
	return np.array([
		[ np.cos(theta), 0., -np.sin(theta)],
		[0., 1., 0.],
		[ np.sin(theta), 0.,  np.cos(theta)]])

def rotz(theta):
	return np.array([
		[ np.cos(theta), -np.sin(theta), 0.],
		[ np.sin(theta),  np.cos(theta), 0.],
		[0., 0., 1.]])

def vec_print(vec):
	ret = ""
	it  = 0
	for comp in vec:
		ret += str(comp)
		it += 1
		if it != vec.size:
			ret += ", "
	return ret