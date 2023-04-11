#visualization of particle motion for viscek flocking. rectangular box with reflecting or periodic boundaries

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

np.random.seed(00)
L = 1; #size of box

N = 5; #number of particles
v = 1; #velocity
dt = 0.001;
R = 0.01; #alignment distance
eta = 20; #noise accumulated in t=1s, i.e. after 1/dt timesteps is -eta pi to  eta pi


pos0 = L*np.random.rand(N,2);

#pos0 = np.concatenate((L*np.ones((N,1)), L*np.random.rand(N,1)), axis=1);
pos1 = np.zeros((N,2));

theta0 = 2*np.pi*np.random.rand(N,1);
theta1 = np.zeros((N,1));

fig, ax = plt.subplots()
t=0;

counter=0;
counterprint = 0.01/dt;

def update(frame):
	global pos0,pos1,theta0,theta1, L,N,v,dt,R,eta,t, counter, counterprint

	if counter%counterprint==0:
		ax.clear()
		ax.scatter(pos0[:,0], pos0[:,1],color='blue',marker='o')
		ax.set_xlim(0,L)
		ax.set_ylim(0,L)
		ax.set_title(t)
		#plt.pause(0.001)

	pos1 = pos0+v*np.concatenate((np.cos(theta0), np.sin(theta0)), axis=1 )*dt;

	for i in range(N):
		#find indices of particles within R
		dist = np.sqrt(np.square(pos0[:,0] - pos0[i,0]) + np.square(pos0[:,1]-pos0[i,1]));
		ind = dist<R;
		theta1[i] = np.mean(theta0[ind]) + eta*2*np.pi*(0.5-np.random.rand())*dt; #random noise (-eta*pi, eta*pi)
		# ind[i] = False; #don't include the particle itself
		# if np.any(ind):
		# 	theta1[i] = np.mean(theta0[ind]) + eta*2*np.pi*(0.5-np.random.rand())*dt; #random noise (-eta*pi, eta*pi)
		# else:
		# 	theta1[i] = theta0[i] + eta*2*np.pi*(0.5-np.random.rand())*dt;

	# #periodic boundary conditions, elicits no change in velocity direction
	# bnd1 = pos1<0;
	# pos1[bnd1] = L+pos1[bnd1];
	# bnd2 = pos1>L;
	# pos1[bnd2] = pos1[bnd2]-L;

	#reflective boundary conditions, changes velocity direction too
	
	bnd = pos1[:,0]<0; #left boundary
	pos1[bnd,0] = -pos1[bnd,0];
	theta1[bnd] = 3*np.pi - theta1[bnd];

	bnd = pos1[:,0]>L; #right boundary
	pos1[bnd,0] = 2*L-pos1[bnd,0];
	theta1[bnd] = 3*np.pi - theta1[bnd];

	bnd = pos1[:,1]<0; #bottom boundary
	pos1[bnd,1] = -pos1[bnd,1];
	theta1[bnd] =  2*np.pi- theta1[bnd];

	bnd = pos1[:,1]>L; #top boundary
	pos1[bnd,1] = 2*L-pos1[bnd,1];
	theta1[bnd] =  2*np.pi- theta1[bnd];

	pos0 = pos1;
	theta0 = np.mod(theta1,2*np.pi);
	t=t+dt;
	counter +=1;

ani = FuncAnimation(fig, update, interval=1)
plt.show()