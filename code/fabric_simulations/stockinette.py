'''
* Main code for the knitting simulation
* M. Dimitriyev, S. Gonzalez
'''
# command to run in terminal:
# python3 stockinette.py chosen_config_file.dat

# needed to run: config file, init file (both are dat file types)
# in same folder: Bezier.py, ConnectingCurve.py, Crossover.py, mathhelper.py
# Bezier, ConnectingCurve, Crossover all construct and update the bezier curve splines
# this code utilizes fifth order bezier curves


# Sweep over cell dimension values in either x or y direction
# analogous to stetching the fabric in respective direction
# step_size is set for both x,y directions

import sys
import subprocess
import numpy as np
import Bezier as bez
import CrossOver as co
import ConnectingCurve as cc
from mathHelper import *
from scipy.optimize import fmin_cg, minimize
from numpy.linalg import norm
from scipy.integrate import quad
from numpy.linalg import matrix_power
from copy import deepcopy

PI = np.pi


# initialize parameters from config file
init_file = ""
init_file_flag = False

output_prefix = "stockinette_out_" # default
output_prefix_flag = True

stitch_length = 0.0
stitch_length_flag = False

bending_modulus = 0.0
bending_modulus_flag = False

regularization = 0.0
regularization_flag = False

yarn_radius = 0.0
yarn_radius_flag = False

core_radius = 0.0
core_radius_flag = False

interactions_per_segment = 0
interactions_per_segment_flag = False

max_twist_constraint = 0.9 # default
max_twist_constraint_flag = True

contact_rigidity = 0.0
contact_rigidity_flag = False

contact_exponent = 0.0
contact_exponent_flag = False 

ftol = 1e-5 # default
ftol_flag = True

sweep_direction = ""
sweep_direction_flag = False

init_dimension = 0.0
init_dimension_flag = False

final_dimension = 0.0
final_dimension_flag = False



# read in parameter file
config_file_name = sys.argv[1]

with open(config_file_name) as config_file:
	for line in config_file:
		words = line.split()
		if len(words) == 0: continue
		
		if words[0] == "*init_file":
			init_file = words[1]
			init_file_flag = True
		elif words[0] == "*output_prefix":
			output_prefix = words[1]
			output_prefix_flag = True
		elif words[0] == "*stitch_length":
			stitch_length = float(words[1])
			stitch_length_flag = True
		elif words[0] == "*bending_modulus":
			bending_modulus = float(words[1])
			bending_modulus_flag = True
		elif words[0] == "*regularization":
			regularization = float(words[1])
			regularization_flag = True
		elif words[0] == "*yarn_radius":
			yarn_radius = float(words[1])
			yarn_radius_flag = True
		elif words[0] == "*core_radius":
			core_radius = float(words[1])
			core_radius_flag = True
		elif words[0] == "*interactions_per_segment":
			interactions_per_segment = int(words[1])
			interactions_per_segment_flag = True
		elif words[0] == "*max_twist_constraint":
			max_twist_constraint = float(words[1])
			max_twist_constraint_flag = True
		elif words[0] == "*contact_rigidity":
			contact_rigidity = float(words[1])
			contact_rigidity_flag = True
		elif words[0] == "*contact_exponent":
			contact_exponent = float(words[1])
			contact_exponent_flag = True
		elif words[0] == "*ftol":
			ftol = float(words[1])
			ftol_flag = True
		elif words[0] == "*sweep_direction":
			sweep_direction = words[1]
			if sweep_direction == 'y' or sweep_direction == 'Y':
				sweep_direction = 'y'
			else:
				sweep_direction = 'x'
			sweep_direction_flag = True
		elif words[0] == "*init_dimension":
			init_dimension = float(words[1])
			init_dimension_flag = True
		elif words[0] == "*final_dimension":
			final_dimension = float(words[1])
			final_dimension_flag = True


# if not all params are specified, report error and exit
if not init_file_flag:
	print("ERROR: missing init_file")
	exit()
if not output_prefix_flag:
	print("ERROR: missing output_prefix")
	exit()
if not stitch_length_flag:
	print("ERROR: missing stitch_length")
	exit()
if not bending_modulus_flag:
	print("ERROR: missing bending_modulus")
	exit()
if not regularization_flag:
	print("ERROR: missing regularization")
	exit()
if not yarn_radius_flag:
	print("ERROR: missing yarn_radius")
	exit()
if not core_radius_flag:
	print("ERROR: missing core_radius")
	exit()
if not interactions_per_segment_flag:
	print("ERROR: missing interactions_per_segment")
	exit()
if not max_twist_constraint_flag:
	print("ERROR: missing max_twist_constraint")
	exit()
if not contact_rigidity_flag:
	print("ERROR: missing contact_rigidity")
	exit()
if not contact_exponent_flag:
	print("ERROR: missing contact_exponent")
	exit()
if not ftol_flag:
	print("ERROR: missing ftol")
	exit()
if not sweep_direction_flag:
	print("ERROR: missing sweep_direction")
	exit()
if not init_dimension_flag:
	print("ERROR: missing init_dimension")
	exit()
if not final_dimension_flag:
	print("ERROR: missing final_dimension")
	exit()

# create output directory
output_directory_name = output_prefix + "_" + init_file +"_" + sweep_direction
subprocess.run(["mkdir",output_directory_name])

# write current configuration to directory and echo to terminal
ConfigOut = open(output_directory_name + "/config","w")
print("---- simulation's configuration ----")
ConfigOut.write("---- simulation's configuration ---- \n")
print("init_file " + init_file)
ConfigOut.write("*init_file " + init_file + "\n")
print("output_prefix " + output_prefix)
ConfigOut.write("*output_prefix " + output_prefix + "\n")
print("stitch_length " + str(stitch_length))
ConfigOut.write("*stitch_length " + str(stitch_length) + "\n")
print("bending_modulus " + str(bending_modulus))
ConfigOut.write("*bending_modulus " + str(bending_modulus) + "\n")
print("regularization " + str(regularization))
ConfigOut.write("*regularization " + str(regularization) + "\n")
print("yarn_radius " + str(yarn_radius))
ConfigOut.write("*yarn_radius " + str(yarn_radius) + "\n")
print("core_radius " + str(core_radius))
ConfigOut.write("*core_radius " + str(core_radius) + "\n")
print("interactions_per_segment " + str(interactions_per_segment))
ConfigOut.write("*interactions_per_segment " + str(interactions_per_segment) + "\n")
print("max_twist_constraint " + str(max_twist_constraint))
ConfigOut.write("*max_twist_constraint " + str(max_twist_constraint) + "\n")
print("contact_rigidity " + str(contact_rigidity))
ConfigOut.write("*contact_rigidity " + str(contact_rigidity) + "\n")
print("contact_exponent " + str(contact_exponent))
ConfigOut.write("*contact_exponent " + str(contact_exponent) + "\n")
print("ftol " + str(ftol))
ConfigOut.write("*ftol " + str(ftol) + "\n")
print("sweep_direction " + sweep_direction)
ConfigOut.write("*sweep_direction " + sweep_direction + "\n")
print("init_dimension " + str(init_dimension))
ConfigOut.write("*init_dimension " + str(init_dimension) + "\n")
print("final_dimension " + str(final_dimension))
ConfigOut.write("*final_dimension " + str(final_dimension) + "\n")
ConfigOut.close()


Consts = {
	'BendingModulus': 	bending_modulus,
	'TwistingModulus': 	0,
	'Regularization':	regularization,
	'PreferredTwist':	0.,		
	'CoreRadius': 		core_radius,
	'TotalLength':		stitch_length,
	'YarnRadius':		yarn_radius,
	'NumberBeads':		interactions_per_segment,
	'TwMax':            max_twist_constraint,
	'CompSpringConst':	contact_rigidity,
	'Power':			contact_exponent
}

# provides rotation for coursewise neighboring cells; x
def rotate_cell_c(cell_angle, join_angle):
	ca = cell_angle
	ja = join_angle

	return np.array([
		[np.cos(ja), 0, -np.sin(ja)],
		[0, 1, 0],
		[np.sin(ja), 0 , np.cos(ja)]
		])


# provides translation for coursewise neighboring cells; x
def translate_cell_c(c, w, cell_angle, join_angle):
	return c/2*(np.array([1.,0,0]) + np.dot(rotate_cell_c(cell_angle, join_angle),np.array([1,0,0])))


# provides rotation for walewise neighboring cells; y
def rotate_cell_w(cell_angle, join_angle):
	ca = cell_angle
	ja = join_angle

	return np.array([
		[1, 0, 0],
		[0, np.cos(ja), -np.sin(ja)],
		[0, np.sin(ja), np.cos(ja)]
		])


# provides translation for walewise neighboring cells; y
def translate_cell_w(c, w, cell_angle, join_angle):
	return w/2*(np.array([0,1.,0]) + np.dot(rotate_cell_w(cell_angle, join_angle),np.array([0,1,0])))



# define energy cost for deforming yarn; an energy density
# all of this is model dependent
def floof(x):
	# scale of the energy density
	k = contact_rigidity
	# power that shows up in energy density model
	p 	= contact_exponent
	dia_outer = 2*yarn_radius
	dia_core  = 2*core_radius
	dia_diff  = dia_outer - dia_core
	EPS = 1e-5

	dist = x + dia_outer
	z = (dist - dia_core)/dia_diff

	#energy density
	v = lambda y: k*( ((dia_diff)**2)/(p*(p - 1)) )*( y**(1-p) - 1 + (1-p)*(1-y) )

	if z >= EPS:
		return v(z)
	else:
		# ensure that the core is extremely stiff
		#return v(EPS) + 10.*(dist**(-1.) - (dia_core + EPS*dia_diff)**(-1))
		return v(EPS) + 100.*((dist - (dia_core + EPS*dia_diff))/(dia_core + EPS*dia_diff))**2


# numerical integration using Simpson's rule
def int_simpson(func_array, dx, limit_start=0, limit_end=-1):
	if limit_end == -1:
		limit_end = len(func_array)-1

	tot = func_array[limit_start] + func_array[limit_end]

	for n in np.arange(1,limit_end-limit_start):
		coeff=4
		if n%2==0:
			coeff = 2

		tot += coeff*func_array[limit_start+n]
	
	return tot*(dx/3.)



# read init file
x_init_l = []
with open("Init_" + str(init_file) + ".dat") as f:
	for line in f:
		x_init_l.append(float(line))

x_init = np.asarray(x_init_l)

# cell dimensions
ax = 0
ay = 0

# last two entries of init are cell dimensions
x0 = deepcopy(x_init[:-2])

if sweep_direction=='y':
	ax = x_init[-2]
	ay = init_dimension
	x0 = np.append(x0,ax)
	step_size = 0.03
else:
	ax = init_dimension
	ay = x_init[-1]
	x0 = np.append(x0,ay)
	step_size = 0.09


cell_angle = 0
join_angle_c = 0.
join_angle_w = 0.

xres_extend1 = []

# start sweep
sweep_flag = True
step_num = 0
while sweep_flag:

	print("\n-------------------")
	print("ax = " + str(ax))
	print("ay = " + str(ay))

	print(x0)

	# create the simulations output files
	if sweep_direction=='y':
		out_subdir = output_directory_name + "/yresults_" + str(ay)
	else:
		out_subdir = output_directory_name + "/xresults_" + str(ax)
	subprocess.run(["mkdir",out_subdir])

	BezierOut = open(out_subdir + "/BezierOut.dat","w")
	BezierOutCourse = open(out_subdir + "/BezierOutCourse.dat","w")
	BezierOutWale = open(out_subdir + "/BezierOutWale.dat","w")
	CellPropsOut = open(out_subdir + "/CellPropsOut.dat","w")
	ResultsOut = open(out_subdir + "/ResultsOut.dat","w")
	StretchOut =  open(out_subdir + "/StretchOut.dat","w")
	EnergyOut = open(out_subdir + "/EnergyOut.dat","w")
	ContactDensityOut = open(out_subdir + "/ContactDensityOut.dat","w")
	BendingDensityOut = open(out_subdir + "/BendingDensityOut.dat","w")
	ContactMapOut = open(out_subdir + "/ContactMapOut.dat", "w")


	tr_c  = translate_cell_c(ax, ay, cell_angle, join_angle_c)
	tr_w  = translate_cell_w(ax, ay, cell_angle, join_angle_w)
	
	# make the splines
	co1 = co.CrossOver5(np.array([-ax/4,0,0]), 
						np.array([x0[0],x0[1],x0[2]]), x0[3], x0[4],
					   	np.array([x0[5],x0[6]]),   x0[7],  x0[8],  x0[9],  x0[10],
					   	np.array([x0[11],x0[12]]), x0[13], x0[14], x0[15], x0[16], 
					   	np.array([x0[17],x0[18]]), x0[19], x0[20], x0[21], x0[22],
					   	np.array([x0[5],x0[6]]),   x0[7],  x0[8],  x0[9],  x0[10],
					   	x0[23],
					   	0,Consts)

	co2 = co.CrossOver5(np.array([ax/4,0,0]), 
						np.array([x0[0],-x0[1],-x0[2]]), x0[3], x0[4],
					   	np.array([x0[11],x0[12]]), x0[13], x0[14], x0[15], x0[16],
					   	np.array([x0[5],x0[6]]),   x0[7],  x0[8],  x0[9],  x0[10],
					   	np.array([x0[5],x0[6]]),   x0[7],  x0[8],  x0[9],  x0[10],
					   	np.array([x0[17],x0[18]]), x0[19], x0[20], x0[21], x0[22],
					   	x0[23],
					   	1,Consts)

	cc12 = cc.ConnectingCurve5(co1.curve[0].at(1), 
							   co1.curve[0].tan_at(1), 		x0[24],
							   co1.curve[0].curv_vec_at(1), x0[25],
					   		   co2.curve[0].at(0), 
					   		   co2.curve[0].tan_at(0), 		x0[24],
					   		   co2.curve[0].curv_vec_at(0), x0[25],
							   Consts)
	
	cc21 = cc.ConnectingCurve5(co2.curve[1].at(1),     	
							   co2.curve[1].tan_at(1), 		x0[24],
							   co2.curve[1].curv_vec_at(1), x0[25],
							   co1.curve[1].at(0)+tr_c, 
							   co1.curve[1].tan_at(0), 		x0[24],
							   co1.curve[1].curv_vec_at(0), x0[25],
							   Consts)
	
	cc11 = cc.ConnectingCurve5(co1.curve[1].at(1), 		
							   co1.curve[1].tan_at(1), 		x0[26],
							   co1.curve[1].curv_vec_at(1), x0[27],
							   co1.curve[0].at(0)+tr_w, 
							   co1.curve[0].tan_at(0), 		x0[26],
							   co1.curve[0].curv_vec_at(0), x0[27],
							   Consts)

	cc22 = cc.ConnectingCurve5(co2.curve[0].at(1)+tr_w, 
							   co2.curve[0].tan_at(1), 		x0[26],
							   co2.curve[0].curv_vec_at(1), x0[27],
							   co2.curve[1].at(0),      
							   co2.curve[1].tan_at(0), 		x0[26],
							   co2.curve[1].curv_vec_at(0), x0[27],
							   Consts)		

	#update the splines
	def update(x, *args):
		ax1 = ax
		ay1 = ay
		if sweep_direction=='y':
			ax1 = x[-1]
		else:
			ay1 = x[-1]

		tr_c  = translate_cell_c(ax1,ay1, cell_angle, join_angle_c)
		tr_w  = translate_cell_w(ax1,ay1, cell_angle, join_angle_w)

		co1.update(np.array([-ax1/4,0,0]), 
				   np.array([x[0],x[1],x[2]]), x[3], x[4],
				   np.array([x[5],x[6]]),   x[7],  x[8],  x[9],  x[10],
				   np.array([x[11],x[12]]), x[13], x[14], x[15], x[16], 
				   np.array([x[17],x[18]]), x[19], x[20], x[21], x[22],
				   np.array([x[5],x[6]]),   x[7],  x[8],  x[9],  x[10],
				   x[23])

		co2.update(np.array([ax1/4,0,0]), 
				   np.array([x[0],-x[1],-x[2]]), x[3], x[4],
				   np.array([x[11],x[12]]), x[13], x[14], x[15], x[16], 
				   np.array([x[5],x[6]]),   x[7],  x[8],  x[9],  x[10],
				   np.array([x[5],x[6]]),   x[7],  x[8],  x[9],  x[10],
				   np.array([x[17],x[18]]), x[19], x[20], x[21], x[22],
				   x[23])

		cc12.update(co1.curve[0].at(1), 
					co1.curve[0].tan_at(1), 	 x[24],
					co1.curve[0].curv_vec_at(1), x[25],
					co2.curve[0].at(0), 
					co2.curve[0].tan_at(0), 	 x[24],
					co2.curve[0].curv_vec_at(0), x[25])

		cc21.update(co2.curve[1].at(1),		 
					co2.curve[1].tan_at(1), 	 x[24],
					co2.curve[1].curv_vec_at(1), x[25],
					co1.curve[1].at(0)+tr_c, 
					co1.curve[1].tan_at(0), 	 x[24],
					co1.curve[1].curv_vec_at(0), x[25])

		cc11.update(co1.curve[1].at(1),		 
					co1.curve[1].tan_at(1), 	 x[26],
					co1.curve[1].curv_vec_at(1), x[27],
					co1.curve[0].at(0)+tr_w, 
					co1.curve[0].tan_at(0), 	 x[26],
					co1.curve[0].curv_vec_at(0), x[27])

		cc22.update(co2.curve[0].at(1)+tr_w, 
					co2.curve[0].tan_at(1), 	 x[26],
					co2.curve[0].curv_vec_at(1), x[27],
					co2.curve[1].at(0),		 
					co2.curve[1].tan_at(0), 	 x[26],
					co2.curve[1].curv_vec_at(0), x[27])

	# used in calculating compression
	def reparam_curve_t(x, t):
		ax1 = ax
		ay1 = ay
		if sweep_direction=='y':
			ax1 = x[-1]
		else:
			ay1 = x[-1]
		reps = np.floor(t/8)
		t_rem = t%8

		tr_c  = translate_cell_c(ax1,ay1, cell_angle, join_angle_c)
		tr_w  = translate_cell_w(ax1,ay1, cell_angle, join_angle_w)

		# 0 <= t < 1: co1[1]
		# 1 <= t < 2: cc11
		# 2 <= t < 3: co1[0]+tr_w
		# 3 <= t < 4: cc12+tr_w
		# 4 <= t < 5: co2[0]+tr_w
		# 5 <= t < 6: cc22
		# 6 <= t < 7: co2[1]
		# 7 <= t < 8: cc21

		if t_rem >= 0 and t_rem < 1:
			return (co1.curve[1].at(t_rem%1) + reps*tr_c)

		if t_rem >= 1 and t_rem < 2:
			return (cc11.curve.at(t_rem%1) + reps*tr_c)

		if t_rem >= 2 and t_rem < 3:
			return (co1.curve[0].at(t_rem%1) +tr_w + reps*tr_c)

		if t_rem >= 3 and t_rem < 4:
			return (cc12.curve.at(t_rem%1) +tr_w + reps*tr_c)

		if t_rem >= 4 and t_rem < 5:
			return (co2.curve[0].at(t_rem%1) +tr_w + reps*tr_c)

		if t_rem >= 5 and t_rem < 6:
			return (cc22.curve.at(t_rem%1) + reps*tr_c)

		if t_rem >= 6 and t_rem < 7:
			return (co2.curve[1].at(t_rem%1) + reps*tr_c)

		if t_rem >= 7:
			return (cc21.curve.at(t_rem%1) + reps*tr_c)

	# used for compression calculation
	def reparam_curve_dt(x, t):
		ax1 = ax
		ay1 = ay
		if sweep_direction=='y':
			ax1 = x[-1]
		else:
			ay1 = x[-1]
		reps = np.floor(t/8)
		t_rem = t%8

		# 0 <= t < 1: co1[1]
		# 1 <= t < 2: cc11
		# 2 <= t < 3: co1[0]+tr_w
		# 3 <= t < 4: cc12+tr_w
		# 4 <= t < 5: co2[0]+tr_w
		# 5 <= t < 6: cc22
		# 6 <= t < 7: co2[1]
		# 7 <= t < 8: cc21

		if t_rem >= 0 and t_rem < 1:
			return norm(co1.curve[1].d_at(t%1))

		if t_rem >= 1 and t_rem < 2:
			return norm(cc11.curve.d_at(t%1))

		if t_rem >= 2 and t_rem < 3:
			return norm(co1.curve[0].d_at(t%1))

		if t_rem >= 3 and t_rem < 4:
			return norm(cc12.curve.d_at(t%1))

		if t_rem >= 4 and t_rem < 5:
			return norm(co2.curve[0].d_at(t%1))

		if t_rem >= 5 and t_rem < 6:
			return norm(cc22.curve.d_at(t%1))

		if t_rem >= 6 and t_rem < 7:
			return norm(co2.curve[1].d_at(t%1))

		if t_rem >= 7:
			return norm(cc21.curve.d_at(t%1))

	# compute contact energy from compression
	def contact_energy(x, print_to_file=False):
		ax1 = ax
		ay1 = ay
		if sweep_direction=='y':
			ax1 = x[-1]
		else:
			ay1 = x[-1]
		ret = 0
		rad = yarn_radius

		tr_c  = translate_cell_c(ax1,ay1, cell_angle, join_angle_c)
		tr_w  = translate_cell_w(ax1,ay1, cell_angle, join_angle_w)

		num_contact_1 = interactions_per_segment
		num_contact_2 = interactions_per_segment

		dt1 = 1./num_contact_1
		dt2 = 1./num_contact_2

		# evaluatee coordinates of points bisecting the egdes between mesh vertices
		curve_points_mid = [reparam_curve_t(x,t) for t in np.arange(dt1/2, 1.0, dt1)]
		curve_points_mid =  np.append(curve_points_mid,[reparam_curve_t(x,t) for t in np.arange(1.0+(dt2/2), 2.0, dt2)], axis=0)
		curve_points_mid =  np.append(curve_points_mid,[reparam_curve_t(x,t) for t in np.arange(2.0+(dt1/2), 3.0, dt1)], axis=0)
		curve_points_mid =  np.append(curve_points_mid,[reparam_curve_t(x,t) for t in np.arange(3.0+(dt2/2), 4.0, dt2)], axis=0)
		curve_points_mid =  np.append(curve_points_mid,[reparam_curve_t(x,t) for t in np.arange(4.0+(dt1/2), 5.0, dt1)], axis=0)
		curve_points_mid =  np.append(curve_points_mid,[reparam_curve_t(x,t) for t in np.arange(5.0+(dt2/2), 6.0, dt2)], axis=0)
		curve_points_mid =  np.append(curve_points_mid,[reparam_curve_t(x,t) for t in np.arange(6.0+(dt1/2), 7.0, dt1)], axis=0)
		curve_points_mid =  np.append(curve_points_mid,[reparam_curve_t(x,t) for t in np.arange(7.0+(dt2/2), 8.0, dt2)], axis=0)

		# evaluate "speed" of points bisecting the egdes between mesh vertices
		curve_points_mid_dt = [reparam_curve_dt(x,t) for t in np.arange(dt1/2, 1.0, dt1)]
		curve_points_mid_dt =  np.append(curve_points_mid_dt,[reparam_curve_dt(x,t) for t in np.arange(1.0+(dt2/2), 2.0, dt2)], axis=0)
		curve_points_mid_dt =  np.append(curve_points_mid_dt,[reparam_curve_dt(x,t) for t in np.arange(2.0+(dt1/2), 3.0, dt1)], axis=0)
		curve_points_mid_dt =  np.append(curve_points_mid_dt,[reparam_curve_dt(x,t) for t in np.arange(3.0+(dt2/2), 4.0, dt2)], axis=0)
		curve_points_mid_dt =  np.append(curve_points_mid_dt,[reparam_curve_dt(x,t) for t in np.arange(4.0+(dt1/2), 5.0, dt1)], axis=0)
		curve_points_mid_dt =  np.append(curve_points_mid_dt,[reparam_curve_dt(x,t) for t in np.arange(5.0+(dt2/2), 6.0, dt2)], axis=0)
		curve_points_mid_dt =  np.append(curve_points_mid_dt,[reparam_curve_dt(x,t) for t in np.arange(6.0+(dt1/2), 7.0, dt1)], axis=0)
		curve_points_mid_dt =  np.append(curve_points_mid_dt,[reparam_curve_dt(x,t) for t in np.arange(7.0+(dt2/2), 8.0, dt2)], axis=0)

		# assign differentials to mesh edges
		dt_array = [dt1 for t in np.arange(dt1/2, 1.0, dt1)]
		dt_array =  np.append(dt_array,[dt2 for t in np.arange(1.0+(dt2/2), 2.0, dt2)])
		dt_array =  np.append(dt_array,[dt1 for t in np.arange(2.0+(dt1/2), 3.0, dt1)])
		dt_array =  np.append(dt_array,[dt2 for t in np.arange(3.0+(dt2/2), 4.0, dt2)])
		dt_array =  np.append(dt_array,[dt1 for t in np.arange(4.0+(dt1/2), 5.0, dt1)])
		dt_array =  np.append(dt_array,[dt2 for t in np.arange(5.0+(dt2/2), 6.0, dt2)])
		dt_array =  np.append(dt_array,[dt1 for t in np.arange(6.0+(dt1/2), 7.0, dt1)])
		dt_array =  np.append(dt_array,[dt2 for t in np.arange(7.0+(dt2/2), 8.0, dt2)])

		# use the midpoint rule to integrate arclength parameter onto mesh vertices
		curve_arclength = [sum(np.multiply(dt_array[:idx],curve_points_mid_dt[:idx])) for idx in range(1,len(curve_points_mid)+1)]

		# approximate arclength evaluated at midpoints as the average of the arclength at two bounding vertices
		curve_arclength_mid = [0.5*(curve_arclength[idx] + (0 if idx==0  else curve_arclength[idx-1]) ) for idx in range(0,len(curve_points_mid))]


		sTot = curve_arclength[-1]
		s_excl = 2.1*rad


		for idx1 in range(0,len(curve_points_mid)):
			pt1 = curve_points_mid[idx1]
			s1 = curve_arclength_mid[idx1]

			en_dens = 0.0

			# sweep around surrounding points and calculate compression energy density
			for n in range(-3,4):
				for m in range(-1,2):
					for idx2 in range(0,len(curve_points_mid)):
						pt2 = curve_points_mid[idx2] + m*tr_c + n*tr_w
						s2 = curve_arclength_mid[idx2] + m*sTot

						if ( (n==0 and abs(s1-s2) > s_excl) or n != 0):
							dist = norm(pt1 - pt2) - 2*rad
							
							if dist < 0:
								en_dens_contribution = floof(dist)*curve_points_mid_dt[idx1]*curve_points_mid_dt[idx2]*dt_array[idx1]*dt_array[idx2]
								en_dens += en_dens_contribution

								if print_to_file:
									ContactMapOut.write("%6.5f, %8.6f, %8.6f, %8.6f, %8.6f, %8.6f, %8.6f, %10.8f\n" % (s1, pt1[0], pt1[1], pt1[2], pt2[0], pt2[1], pt2[2], en_dens_contribution) )

			ret += en_dens

			if print_to_file:
				ContactDensityOut.write("%6.5f, %8.6f, %8.6f, %8.6f, %10.8f\n" % (s1, pt1[0], pt1[1], pt1[2], en_dens) )

		return ret

	# the bending energy from the splines
	def curve_energy(x):
		ret =  2*co1.energy()+2*cc11.energy()+cc12.energy()+cc21.energy()
		return ret

	# the energy functional being minimized
	def f(x, *args):		
		update(x,*args)

		cont_e  = contact_energy(x)
		curve_e = curve_energy(x) 

		return(cont_e + curve_e)


	# CONSTRAINT FUNCTIONS
	def lNew(x):
		update(x,0)
		return(2*co1.curve[0].length()+2*co1.curve[1].length()+2*cc11.curve.length()
			+cc12.curve.length()+cc21.curve.length())

	def tw_fct(x):
		update(x,0)
		return np.cos(x[23])


	# CONSTRAINED MINIMIZATION OCCURS HERE
	opts = {'maxiter' : None,    # default value None
			'disp' : False,    # non-default value.
			'ftol' : ftol,    
			'norm' : np.inf,  # default value.
			'eps' : 1.4901161193847656e-09}  # default value  1.4901161193847656e-08

	cons = ({'type': 'eq', 'fun': lambda x: lNew(x)-stitch_length},
			{'type': 'ineq', 'fun': lambda x: tw_fct(x) - max_twist_constraint})
	res = minimize(f, x0, constraints=cons,method='SLSQP')
	xres = res['x']
	print(xres)


	#writes to contact_density_out, see fnc
	contact_energy(xres, print_to_file=True)

	# writes to bezier out
	for i in range(0,2):
		BezierOut.write(co1.curve[i].print())
		BezierOut.write(co2.curve[i].print())
	BezierOut.write(cc12.curve.print())

	BezierOutWale.write(cc11.curve.print())
	BezierOutWale.write(cc22.curve.print())
	BezierOutCourse.write(cc21.curve.print())

	if sweep_direction=='y':
		ax = xres[-1]
	else:
		ay = xres[-1]

	# writes to cell props
	CellPropsOut.write(str(xres[0]) + ", " + str(cell_angle) + ", "  + str(ax) 
		+ ", " + str(ay) + ", " + str(join_angle_c) + ", " + str(join_angle_w))

	# write to resultsout, (can be used as new init file if you include the missing cell dimension)
	xres_extend = deepcopy(xres)
	if step_num == 0:
		xres_extend1 = xres_extend
	for el in xres_extend:
		ResultsOut.write(str(el) + "\n")

	#writes to stretchout
	StretchOut.write( str(Consts["TotalLength"]) + ", " + str(ax) + ", " + str(ay) )

	#writes to energy out
	EnergyOut.write(str(contact_energy(xres)) + "\n")
	EnergyOut.write(str(curve_energy(xres)) + "\n")
	EnergyOut.write(str(f(xres,0)))


	#writes to bending density out
	curve_arr = np.array([co1.curve[1],
					  cc11.curve,
					  co1.curve[0],
					  cc12.curve,
					  co2.curve[0],
					  cc22.curve,
					  co2.curve[1],
					  cc21.curve])

	for it1 in [m%8 for m in range(-1,4)]:
		l_offset = 0
		for it2 in [(n-1)%8 for n in range(0,(it1+1)%8)]:
			jac = lambda x: norm(curve_arr[it2].d_at(x))
			l_offset += quad(jac, 0, 1, limit=40,epsrel=1e-6)[0]
		for i in range(0,interactions_per_segment):
			t1 = i/interactions_per_segment
			jac = lambda x: norm(curve_arr[it1].d_at(x))
			BendingDensityOut.write(str(it1) + ", " + str(l_offset+quad(jac, 0, t1, limit=40,epsrel=1e-6)[0]) + ", " + str(.5*bending_modulus*(curve_arr[it1].curv_at(t1)**2.)) + "\n")

	#closes all the files
	BezierOut.close()
	BezierOutWale.close()
	BezierOutCourse.close()
	CellPropsOut.close()
	ResultsOut.close()
	StretchOut.close()
	EnergyOut.close()
	ContactDensityOut.close()
	BendingDensityOut.close()
	ContactMapOut.close()

	# iterate sweep
	current_dimension = 0
	if sweep_direction=='y':
		current_dimension = ay
	else:
		current_dimension = ax


	if current_dimension <= final_dimension:
		current_dimension += step_size
	else:
		sweep_flag = False
		break


	x0 = deepcopy(xres_extend)
	if sweep_direction=='y':
		ay = current_dimension
	else:
		ax = current_dimension

	step_num += 1

# end loop
print("END")