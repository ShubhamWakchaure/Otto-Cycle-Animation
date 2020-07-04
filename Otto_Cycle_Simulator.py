# Otto cycle Simulator

import math
import matplotlib.pyplot as plt
import numpy as np


def engine_kinematics(bore,stroke,con_rod,cr,start_crank,end_crank):
	# geometric parametrs
	a = stroke/2
	R = con_rod/a

	#Volume paramaters
	V_s  = (math.pi)*(1/4)*pow(bore,2)*stroke
	V_c = V_s/(cr-1)

	sc = math.radians(start_crank) # Start of crank angle
	ec = math.radians(end_crank) # End of crank angle

	num_values = 50

	dtheta = (ec-sc)/(num_values-1)
	V =[]
	for i in range(0,num_values):
		theta = sc + +i*dtheta
		term1 = 0.5*(cr-1)
		term2 = R+1-math.cos(theta)
		term3 = pow(R,2)-pow(math.sin(theta),2)
		term3 = pow(term3,0.5)
		V.append((1+term1*(term2-term3))*V_c)

	return V


#input

t3 = 2300 #in kelvin
gamma = 1.4

#geometric parameters
bore = 0.1 # in meters
stroke = 0.1 # in meters 
con_rod = 0.15 # in meters
cr = np.arange(1.5,15,0.25)

#volume computation
v_s  = (math.pi/4)*pow(bore,2)*stroke
ct = 1 # counter variable for saving each plot
for m in cr:
	#state point 1 start of compression stroke
	v_c = v_s/(m-1) 
	v1 = v_c+v_s
	p1 = 101325
	t1 = 500

	# # state point 2 start of heat addition
	v2 = v_c
	p2 = p1*pow(v1,gamma)/pow(v2,gamma) #p2v2^gamma  = p1v1^gamme
	rhs = p1*v1/t1 #p2v2/t2 = p1v1/t1
	t2 = p2*v2/rhs
	V_compression = engine_kinematics(bore,stroke,con_rod,m,0,180)
	constant = p1*pow(v1,gamma)
	P_compression = []

	for v in V_compression:
		P_compression.append(constant/pow(v,gamma))

	# comp_x = (min(V_compression)+max(V_compression))/2
	# comp_y = constant/pow(comp_x,gamma)

	#state point 3
	v3 = v2
	rhs = p2*v2/t2 #p3v3/t3 = p2v2/t2
	p3 = rhs*t3/v3
	V_expansion = engine_kinematics(bore,stroke,con_rod,m,0,180)
	constant = p3*pow(v3,gamma)
	P_expansion = []

	for v in V_expansion:
		P_expansion.append(constant/pow(v,gamma))


	# state point 4
	v4 = v1
	p4 = p3*pow(v3,gamma)/pow(v4,gamma)#p4v4^gamma = p3v3^gamma
	t4 = p4*v4/rhs #p4v4/t4 = p3v3/t3
	rp = p3/p2
	efficiency  = 1-1/pow(m,gamma)
	filename = 'Otto_Cycle_%05d.png' %ct
	ct = ct+1
	# Plotting the otto cycle
	plt.figure()
	plt.plot([v2 ,v3],[p2,p3]) # Contant pressure heat addition
	plt.plot(V_compression,P_compression) # Conpression
	plt.plot(V_expansion,P_expansion) # Expansion
	plt.plot([v1 ,v4],[p1,p4])#Constant pressure heat rejection
	plt.plot(v1,p1,'o') # Plotting first point
	plt.plot(v2,p2,'o') # Plotting second point
	plt.plot(v3,p3,'o') # Plotting third point
	plt.plot(v4,p4,'o') # Plotting fourth point
	plt.title('Otto Cycle for CR =%f ' %m)
	plt.xlabel('Volume in $m^3$')
	plt.ylabel('Pressure in MPa')
	plt.legend(('Heat addition','Compression','Expansion','Heat rejection'))
	#plt.show()
	plt.text((v1+0.5*v2)/2,(0.7*p2+p3)/2,'Efficiency is %f'%efficiency)
	plt.savefig(filename)


# After this stitch all the images using Imagemagick software 




