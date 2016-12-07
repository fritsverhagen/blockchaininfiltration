import numpy as np
from scipy.optimize import minimize
import time

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
	
	print m
	print ""
	for i in range(poolamount):
		out=""
		for j in range(poolamount):
			out+="%f "%x[i][j]
		print out
	print ""
	
print time.strftime("%H:%M:%S")

poolamount=5 #amount of pools
m=[] #list of power per pool
x=[[0.0 for j in range(poolamount)] for i in range(poolamount)] #initial infiltration power

print sum(m)

simulate(poolamount,m,x)

print time.strftime("%H:%M:%S")