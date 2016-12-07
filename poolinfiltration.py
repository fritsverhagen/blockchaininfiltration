import numpy as np
from scipy.optimize import minimize
import time
import copy

def kd(i,j):
	if (i==j):
		return 1
	else:
		return 0

def r(x,m,poolamount,p):
	A=np.array([[kd(i,j)-x[i][j]/(m[i]+sum([x[k][i] for k in range(poolamount)])) for j in range(poolamount)] for i in range(poolamount)])
	Ar=np.array([[(-kd(j,p)+1)*A[i][j]+kd(j,p)*((m[i]-sum(x[i]))/((1-sum([sum(x[k]) for k in range(poolamount)]))*(m[i]+sum([x[k][i] for k in range(poolamount)])))) for j in range(poolamount)] for i in range(poolamount)])
	return np.linalg.det(Ar)/np.linalg.det(A)

def simulate(poolamount,m,x):
	for n in range(int(1.5*poolamount**2)):
		p=n%poolamount
		def func(y):
			return -r([x[i] for i in range(p)]+[[y[i] for i in range(p)]+[0.0]+[y[p+i] for i in range(poolamount-1-p)]]+[x[p+1+i] for i in range(poolamount-p-1)],m,poolamount,p)
		cons=[]
		for j in range(poolamount-1):
			def f(y,a=j): return y[a]
			def g(y,a=j): return np.array([0.0 for i in range(a)]+[1.0]+[0.0 for i in range(poolamount-2-a)])
			cons.append({'type':'ineq','fun':f,'jac':g})
		cons.append({'type':'ineq','fun':lambda y:m[p]-sum(y),'jac':lambda y:np.array([-1.0 for i in range(poolamount-1)])})
		
		res=minimize(func,[0 for i in range(poolamount-1)],jac=False,constraints=cons,method='SLSQP')
		for i in range(poolamount):
			if i<p:
				x[p][i]=res.x[i]
			if i>p:
				x[p][i]=res.x[i-1]
						
#	out=""
#	for p in range(poolamount):
#		for pp in range(poolamount):
#			out+=" %f"%x[p][pp]
#	return out
	
#	print m
#	print ""
#	for i in range(poolamount):
#		out=""
#		for j in range(poolamount):
#			out+="%f "%x[i][j]
#		print out
#	print ""
#	#print data
#	print [r(x,m,poolamount,i) for i in range(poolamount)]f=open("data.dat","w")

	
	return x[0][1],x[0][poolamount-1],r(x,m,poolamount,0),r(x,m,poolamount,poolamount-1)
	
poolamount=7 #amount of pools
#m=[0.1,0.1,0.1,0.1,0.1,0.1,0.1] #list of power per pool
#x=[[0.0 for j in range(poolamount)] for i in range(poolamount)] #initial infiltration power
#
#print sum(m)

#print time.strftime("%H:%M:%S")
#
f=open("data.dat","w")
#
#for poolamount in [5]:
#	ix=[[0 for j in range(poolamount)] for i in range(poolamount)]
#	ran=[0.015*(i+1) for i in range(66)]
#	m=[[rani] for rani in ran]
#	for i in range(poolamount-1):
#		mtemp=[]
#		for x in m:
#			for rani in ran:
#				if rani<=x[-1] and sum(x)+rani<=1.0:
#					mtemp.append(x+[rani])
#		m=mtemp
#		
#	for x in m:
#		f.write("%i"%poolamount)
#		for xi in x:
#			f.write(" %f"%xi)
#		f.write(simulate(poolamount,x,ix))
#		f.write("\n")
#		
#
#f.close()
#
#print time.strftime("%H:%M:%S")
	
#	m=[ran for n in range(poolamount)]
#	mt=[m[i][j] for ]
#	for a in [ran for n in range(poolamount)]:

ran=[0.01*(i+1) for i in range(100)]
for poolamount in [3,4,5,6,7,8,9,10,11,12]:
	x=[[0.0 for j in range(poolamount)] for i in range(poolamount)] #initial infiltration power
	
	for b in ran:
		m=[]
		for i in range(poolamount-1):
			m+=[0.03]
		d1,d2,d3,d4=simulate(poolamount,m+[b],x)
		f.write("%f %f %f %f %f %f %i\n"%(0.03,b,d1,d2,d3,d4,poolamount))

#for a in ran:
#	for b in ran:
#		if 6.0*a+b<=1.0:
#			d1,d2,d3,d4=simulate(7,[a,a,a,a,a,a,b],x)
#			f.write("%f %f %f %f %f %f\n"%(a,b,d1,d2,d3,d4))
			
f.close()

#f.write(simulate(poolamount,m,x))

#f.close()